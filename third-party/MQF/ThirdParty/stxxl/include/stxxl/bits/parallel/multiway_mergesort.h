/***************************************************************************
 *  include/stxxl/bits/parallel/multiway_mergesort.h
 *
 *  Parallel multiway mergesort.
 *  Extracted from MCSTL - http://algo2.iti.uni-karlsruhe.de/singler/mcstl/
 *
 *  Part of the STXXL. See http://stxxl.sourceforge.net
 *
 *  Copyright (C) 2007 Johannes Singler <singler@ira.uka.de>
 *  Copyright (C) 2014 Timo Bingmann <tb@panthema.net>
 *
 *  Distributed under the Boost Software License, Version 1.0.
 *  (See accompanying file LICENSE_1_0.txt or copy at
 *  http://www.boost.org/LICENSE_1_0.txt)
 **************************************************************************/

#ifndef STXXL_PARALLEL_MULTIWAY_MERGESORT_HEADER
#define STXXL_PARALLEL_MULTIWAY_MERGESORT_HEADER

#include <vector>
#include <iterator>
#include <algorithm>

#include <stxxl/bits/config.h>
#include <stxxl/bits/parallel/compiletime_settings.h>
#include <stxxl/bits/parallel/equally_split.h>
#include <stxxl/bits/parallel/multiway_merge.h>
#include <stxxl/bits/parallel/multiseq_selection.h>
#include <stxxl/bits/parallel/settings.h>
#include <stxxl/bits/parallel/timing.h>

STXXL_BEGIN_NAMESPACE

#if STXXL_PARALLEL

namespace parallel {

//! Subsequence description.
template <typename DiffType>
struct PMWMSPiece
{
    //! Begin of subsequence.
    DiffType begin;
    //! End of subsequence.
    DiffType end;
};

/*!
 * Data accessed by all threads.
 *
 * PMWMS = parallel multiway mergesort
 */
template <typename RandomAccessIterator>
struct PMWMSSortingData
{
    typedef typename std::iterator_traits<RandomAccessIterator>::value_type
        ValueType;
    typedef typename std::iterator_traits<RandomAccessIterator>::difference_type
        DiffType;

    //! Input begin.
    RandomAccessIterator source;
    //! Start indices, per thread.
    DiffType* starts;

    /*!
     *  Temporary arrays for each thread.
     *
     * Indirection Allows using the temporary storage in different ways,
     * without code duplication.  \see STXXL_MULTIWAY_MERGESORT_COPY_LAST
     */
    ValueType** temporaries;
#if STXXL_MULTIWAY_MERGESORT_COPY_LAST
    /** Storage in which to sort. */
    RandomAccessIterator* sorting_places;
    /** Storage into which to merge. */
    ValueType** merging_places;
#else
    /** Storage in which to sort. */
    ValueType** sorting_places;
    /** Storage into which to merge. */
    RandomAccessIterator* merging_places;
#endif
    /** Samples. */
    ValueType* samples;
    /** Offsets to add to the found positions. */
    DiffType* offsets;
    /** PMWMSPieces of data to merge \c [thread][sequence] */
    std::vector<PMWMSPiece<DiffType> >* pieces;
};

//! Thread local data for PMWMS.
template <typename RandomAccessIterator>
struct PMWMSSorterPU
{
    /** Total number of thread involved. */
    thread_index_t num_threads;
    /** Number of owning thread. */
    thread_index_t iam;
    /** Pointer to global data. */
    PMWMSSortingData<RandomAccessIterator>* sd;
};

/*!
 * Select samples from a sequence.
 * \param d Pointer to thread-local data. Result will be placed in \c d->ds->samples.
 * \param num_samples Number of samples to select.
 */
template <typename RandomAccessIterator, typename DiffType>
inline void determine_samples(PMWMSSorterPU<RandomAccessIterator>* d,
                              DiffType& num_samples)
{
    PMWMSSortingData<RandomAccessIterator>* sd = d->sd;

    num_samples = SETTINGS::sort_mwms_oversampling * d->num_threads - 1;

    std::vector<DiffType> es(num_samples + 2);
    equally_split(sd->starts[d->iam + 1] - sd->starts[d->iam],
                  (thread_index_t)(num_samples + 1), es.begin());

    for (DiffType i = 0; i < num_samples; i++)
        sd->samples[d->iam * num_samples + i] = sd->source[sd->starts[d->iam] + es[i + 1]];
}

/*!
 * PMWMS code executed by each thread.
 * \param d Pointer to thread-local data.
 * \param comp Comparator.
 */
template <bool Stable, typename RandomAccessIterator, typename Comparator>
inline void parallel_sort_mwms_pu(PMWMSSorterPU<RandomAccessIterator>* d,
                                  Comparator& comp)
{
    typedef typename std::iterator_traits<RandomAccessIterator>::value_type
        ValueType;
    typedef typename std::iterator_traits<RandomAccessIterator>::difference_type
        DiffType;

    Timing<inactive_tag> t;

    t.tic();

    PMWMSSortingData<RandomAccessIterator>* sd = d->sd;
    thread_index_t iam = d->iam;

    // length of this thread's chunk, before merging
    DiffType length_local = sd->starts[iam + 1] - sd->starts[iam];

#if STXXL_MULTIWAY_MERGESORT_COPY_LAST
    typedef RandomAccessIterator SortingPlacesIterator;
    // sort in input storage
    sd->sorting_places[iam] = sd->source + sd->starts[iam];
#else
    typedef ValueType* SortingPlacesIterator;
    // sort in temporary storage, leave space for sentinel
    sd->sorting_places[iam] = sd->temporaries[iam] = static_cast<ValueType*>(::operator new (sizeof(ValueType) * (length_local + 1)));
    // copy there
    std::uninitialized_copy(sd->source + sd->starts[iam], sd->source + sd->starts[iam] + length_local, sd->sorting_places[iam]);
#endif

    // sort locally
    if (Stable)
        std::stable_sort(sd->sorting_places[iam], sd->sorting_places[iam] + length_local, comp);
    else
        std::sort(sd->sorting_places[iam], sd->sorting_places[iam] + length_local, comp);

    STXXL_DEBUG_ASSERT(stxxl::is_sorted(sd->sorting_places[iam], sd->sorting_places[iam] + length_local, comp));

    // invariant: locally sorted subsequence in sd->sorting_places[iam], sd->sorting_places[iam] + length_local

    t.tic("local sort");

    if (SETTINGS::sort_splitting == SETTINGS::SAMPLING)
    {
        DiffType num_samples;
        determine_samples(d, num_samples);

#pragma omp barrier

        t.tic("sample/wait");

#pragma omp single
        std::sort(sd->samples, sd->samples + (num_samples * d->num_threads), comp);

#pragma omp barrier

        for (int s = 0; s < d->num_threads; s++)
        {
            // for each sequence
            if (num_samples * iam > 0)
                sd->pieces[iam][s].begin =
                    std::lower_bound(sd->sorting_places[s],
                                     sd->sorting_places[s] + sd->starts[s + 1] - sd->starts[s],
                                     sd->samples[num_samples * iam],
                                     comp)
                    - sd->sorting_places[s];
            else
                // absolute beginning
                sd->pieces[iam][s].begin = 0;

            if ((num_samples * (iam + 1)) < (num_samples * d->num_threads))
                sd->pieces[iam][s].end =
                    std::lower_bound(sd->sorting_places[s],
                                     sd->sorting_places[s] + sd->starts[s + 1] - sd->starts[s],
                                     sd->samples[num_samples * (iam + 1)],
                                     comp)
                    - sd->sorting_places[s];
            else
                // absolute end
                sd->pieces[iam][s].end = sd->starts[s + 1] - sd->starts[s];
        }
    }
    else if (SETTINGS::sort_splitting == SETTINGS::EXACT)
    {
#pragma omp barrier

        t.tic("wait");

        std::vector<std::pair<SortingPlacesIterator, SortingPlacesIterator> > seqs(d->num_threads);
        for (int s = 0; s < d->num_threads; s++)
            seqs[s] = std::make_pair(sd->sorting_places[s], sd->sorting_places[s] + sd->starts[s + 1] - sd->starts[s]);

        std::vector<SortingPlacesIterator> offsets(d->num_threads);

        // if not last thread
        if (iam < d->num_threads - 1)
            multiseq_partition(seqs.begin(), seqs.end(), sd->starts[iam + 1], offsets.begin(), comp);

        for (int seq = 0; seq < d->num_threads; seq++)
        {
            // for each sequence
            if (iam < (d->num_threads - 1))
                sd->pieces[iam][seq].end = offsets[seq] - seqs[seq].first;
            else
                // absolute end of this sequence
                sd->pieces[iam][seq].end = sd->starts[seq + 1] - sd->starts[seq];
        }

#pragma omp barrier

        for (int seq = 0; seq < d->num_threads; seq++)
        {
            // for each sequence
            if (iam > 0)
                sd->pieces[iam][seq].begin = sd->pieces[iam - 1][seq].end;
            else
                // absolute beginning
                sd->pieces[iam][seq].begin = 0;
        }
    }

    t.tic("split");

    // offset from target begin, length after merging
    DiffType offset = 0, length_am = 0;
    for (int s = 0; s < d->num_threads; s++)
    {
        length_am += sd->pieces[iam][s].end - sd->pieces[iam][s].begin;
        offset += sd->pieces[iam][s].begin;
    }

#if STXXL_MULTIWAY_MERGESORT_COPY_LAST
    // merge to temporary storage, uninitialized creation not possible since
    // there is no multiway_merge calling the placement new instead of the
    // assignment operator
    sd->merging_places[iam] = sd->temporaries[iam] = new ValueType[length_am];
#else
    // merge directly to target
    sd->merging_places[iam] = sd->source + offset;
#endif
    std::vector<std::pair<SortingPlacesIterator, SortingPlacesIterator> > seqs(d->num_threads);

    for (int s = 0; s < d->num_threads; s++)
    {
        seqs[s] = std::make_pair(sd->sorting_places[s] + sd->pieces[iam][s].begin, sd->sorting_places[s] + sd->pieces[iam][s].end);

        STXXL_DEBUG_ASSERT(stxxl::is_sorted(seqs[s].first, seqs[s].second, comp));
    }

    sequential_multiway_merge<Stable, false>(seqs.begin(), seqs.end(), sd->merging_places[iam], length_am, comp);

    t.tic("merge");

    STXXL_DEBUG_ASSERT(stxxl::is_sorted(sd->merging_places[iam], sd->merging_places[iam] + length_am, comp));

#pragma omp barrier

#if STXXL_MULTIWAY_MERGESORT_COPY_LAST
    // write back
    std::copy(sd->merging_places[iam], sd->merging_places[iam] + length_am, sd->source + offset);
#endif

    delete sd->temporaries[iam];

    t.tic("copy back");

    t.print();
}

/*!
 * PMWMS main call.
 * \param begin Begin iterator of sequence.
 * \param end End iterator of sequence.
 * \param comp Comparator.
 * \param num_threads Number of threads to use.
 * \tparam Stable Stable sorting.
 */
template <bool Stable,
          typename RandomAccessIterator, typename Comparator>
inline void
parallel_sort_mwms(RandomAccessIterator begin,
                   RandomAccessIterator end,
                   Comparator comp,
                   int num_threads)
{
    STXXL_PARALLEL_PCALL(end - begin)

    typedef typename std::iterator_traits<RandomAccessIterator>::value_type
        ValueType;
    typedef typename std::iterator_traits<RandomAccessIterator>::difference_type
        DiffType;

    DiffType n = end - begin;

    if (n <= 1)
        return;

    if (num_threads > n)           // at least one element per thread
        num_threads = static_cast<thread_index_t>(n);

    PMWMSSortingData<RandomAccessIterator> sd;

    sd.source = begin;
    sd.temporaries = new ValueType*[num_threads];
#if STXXL_MULTIWAY_MERGESORT_COPY_LAST
    sd.sorting_places = new RandomAccessIterator[num_threads];
    sd.merging_places = new ValueType*[num_threads];
#else
    sd.sorting_places = new ValueType*[num_threads];
    sd.merging_places = new RandomAccessIterator[num_threads];
#endif
    if (SETTINGS::sort_splitting == SETTINGS::SAMPLING)
        sd.samples = new ValueType[num_threads * (SETTINGS::sort_mwms_oversampling * num_threads - 1)];
    else
        sd.samples = NULL;
    sd.offsets = new DiffType[num_threads - 1];
    sd.pieces = new std::vector<PMWMSPiece<DiffType> >[num_threads];
    for (int s = 0; s < num_threads; s++)
        sd.pieces[s].resize(num_threads);
    PMWMSSorterPU<RandomAccessIterator>* pus = new PMWMSSorterPU<RandomAccessIterator>[num_threads];
    DiffType* starts = sd.starts = new DiffType[num_threads + 1];

    DiffType chunk_length = n / num_threads, split = n % num_threads, start = 0;
    for (int i = 0; i < num_threads; i++)
    {
        starts[i] = start;
        start += (i < split) ? (chunk_length + 1) : chunk_length;
        pus[i].num_threads = num_threads;
        pus[i].iam = i;
        pus[i].sd = &sd;
    }
    starts[num_threads] = start;

    //now sort in parallel
#pragma omp parallel num_threads(num_threads)
    parallel_sort_mwms_pu<Stable>(&(pus[omp_get_thread_num()]), comp);

    delete[] starts;
    delete[] sd.temporaries;
    delete[] sd.sorting_places;
    delete[] sd.merging_places;

    if (SETTINGS::sort_splitting == SETTINGS::SAMPLING)
        delete[] sd.samples;

    delete[] sd.offsets;
    delete[] sd.pieces;

    delete[] pus;
}

} // namespace parallel

#endif // STXXL_PARALLEL

STXXL_END_NAMESPACE

#endif // !STXXL_PARALLEL_MULTIWAY_MERGESORT_HEADER
