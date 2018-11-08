/***************************************************************************
 *  include/stxxl/bits/parallel.h
 *
 *  Part of the STXXL. See http://stxxl.sourceforge.net
 *
 *  Copyright (C) 2008-2010 Andreas Beckmann <beckmann@cs.uni-frankfurt.de>
 *  Copyright (C) 2011 Johannes Singler <singler@kit.edu>
 *  Copyright (C) 2015 Timo Bingmann <tb@panthema.net>
 *
 *  Distributed under the Boost Software License, Version 1.0.
 *  (See accompanying file LICENSE_1_0.txt or copy at
 *  http://www.boost.org/LICENSE_1_0.txt)
 **************************************************************************/

#ifndef STXXL_PARALLEL_HEADER
#define STXXL_PARALLEL_HEADER

#include <stxxl/bits/config.h>

#include <cassert>

#if STXXL_PARALLEL
 #include <omp.h>
#endif

#include <stxxl/bits/namespace.h>
#include <stxxl/bits/common/settings.h>
#include <stxxl/bits/verbose.h>

#if defined(_GLIBCXX_PARALLEL)
//use _STXXL_FORCE_SEQUENTIAL to tag calls which are not worthwhile parallelizing
#define _STXXL_FORCE_SEQUENTIAL , __gnu_parallel::sequential_tag()
#else
#define _STXXL_FORCE_SEQUENTIAL
#endif

#if 0
// sorting triggers is done sequentially
#define _STXXL_SORT_TRIGGER_FORCE_SEQUENTIAL _STXXL_FORCE_SEQUENTIAL
#else
// sorting triggers may be parallelized
#define _STXXL_SORT_TRIGGER_FORCE_SEQUENTIAL
#endif

#if !STXXL_PARALLEL
#undef STXXL_PARALLEL_MULTIWAY_MERGE
#define STXXL_PARALLEL_MULTIWAY_MERGE 0
#endif

#if defined(STXXL_PARALLEL_MODE) && ((__GNUC__ * 10000 + __GNUC_MINOR__ * 100) < 40400)
#undef STXXL_PARALLEL_MULTIWAY_MERGE
#define STXXL_PARALLEL_MULTIWAY_MERGE 0
#endif

#if !defined(STXXL_PARALLEL_MULTIWAY_MERGE)
#define STXXL_PARALLEL_MULTIWAY_MERGE 1
#endif

#if !defined(STXXL_NOT_CONSIDER_SORT_MEMORY_OVERHEAD)
#define STXXL_NOT_CONSIDER_SORT_MEMORY_OVERHEAD 0
#endif

#if STXXL_WITH_GNU_PARALLEL
#include <parallel/algorithm>
#else
#include <algorithm>
#endif

#include <stxxl/bits/parallel/multiway_merge.h>

STXXL_BEGIN_NAMESPACE

inline unsigned sort_memory_usage_factor()
{
#if STXXL_PARALLEL && !STXXL_NOT_CONSIDER_SORT_MEMORY_OVERHEAD && defined(STXXL_PARALLEL_MODE)
    return (__gnu_parallel::_Settings::get().sort_algorithm == __gnu_parallel::MWMS && omp_get_max_threads() > 1) ? 2 : 1;   //memory overhead for multiway mergesort
#else
    return 1;                                                                                                                //no overhead
#endif
}

inline void check_sort_settings()
{
#if STXXL_PARALLEL && defined(STXXL_PARALLEL_MODE) && !defined(STXXL_NO_WARN_OMP_NESTED)
    static bool did_warn = false;
    if (!did_warn) {
        if (__gnu_parallel::_Settings::get().sort_algorithm != __gnu_parallel::MWMS) {
            if (omp_get_max_threads() <= 2) {
                did_warn = true;  // no problem with at most 2 threads, no need to check again
            }
            else if (!omp_get_nested()) {
                STXXL_ERRMSG("Inefficient settings detected. To get full potential from your CPU it is recommended to set OMP_NESTED=TRUE in the environment.");
                did_warn = true;
            }
        }
    }
#else
    // nothing to check
#endif
}

inline bool do_parallel_merge()
{
#if STXXL_PARALLEL_MULTIWAY_MERGE && defined(STXXL_PARALLEL_MODE)
    return !stxxl::SETTINGS::native_merge && omp_get_max_threads() >= 1;
#else
    return false;
#endif
}

//! this namespace provides parallel or sequential algorithms depending on the
//! compilation settings. it should be used by all components, where
//! parallelism is optional.
namespace potentially_parallel {

#if STXXL_WITH_GNU_PARALLEL

using __gnu_parallel::sort;
using __gnu_parallel::random_shuffle;

#elif STXXL_PARALLEL

using std::sort;
using std::random_shuffle;

#else

using std::sort;
using std::random_shuffle;

#endif

/*! Multi-way merging dispatcher.
 * \param seqs_begin Begin iterator of iterator pair input sequence.
 * \param seqs_end End iterator of iterator pair input sequence.
 * \param target Begin iterator out output sequence.
 * \param comp Comparator.
 * \param length Maximum length to merge.
 * \return End iterator of output sequence.
 */
template <typename RandomAccessIteratorPairIterator,
          typename RandomAccessIterator3, typename DiffType, typename Comparator>
RandomAccessIterator3
multiway_merge(RandomAccessIteratorPairIterator seqs_begin,
               RandomAccessIteratorPairIterator seqs_end,
               RandomAccessIterator3 target, DiffType length,
               Comparator comp)
{
#if STXXL_PARALLEL
    return stxxl::parallel::multiway_merge(
        seqs_begin, seqs_end, target, length, comp);
#else
    return stxxl::parallel::sequential_multiway_merge<false, false>(
        seqs_begin, seqs_end, target, length, comp);
#endif
}

/*! Multi-way merging dispatcher.
 * \param seqs_begin Begin iterator of iterator pair input sequence.
 * \param seqs_end End iterator of iterator pair input sequence.
 * \param target Begin iterator out output sequence.
 * \param comp Comparator.
 * \param length Maximum length to merge.
 * \return End iterator of output sequence.
 */
template <typename RandomAccessIteratorPairIterator,
          typename RandomAccessIterator3, typename DiffType, typename Comparator>
RandomAccessIterator3
multiway_merge_stable(RandomAccessIteratorPairIterator seqs_begin,
                      RandomAccessIteratorPairIterator seqs_end,
                      RandomAccessIterator3 target, DiffType length,
                      Comparator comp)
{
#if STXXL_PARALLEL
    return stxxl::parallel::multiway_merge_stable(
        seqs_begin, seqs_end, target, length, comp);
#else
    return stxxl::parallel::sequential_multiway_merge<true, false>(
        seqs_begin, seqs_end, target, length, comp);
#endif
}

/*! Multi-way merging front-end.
 * \param seqs_begin Begin iterator of iterator pair input sequence.
 * \param seqs_end End iterator of iterator pair input sequence.
 * \param target Begin iterator out output sequence.
 * \param comp Comparator.
 * \param length Maximum length to merge.
 * \return End iterator of output sequence.
 * \pre For each \c i, \c seqs_begin[i].second must be the end marker of the sequence, but also reference the one more sentinel element.
 */
template <typename RandomAccessIteratorPairIterator,
          typename RandomAccessIterator3, typename DiffType, typename Comparator>
RandomAccessIterator3
multiway_merge_sentinels(RandomAccessIteratorPairIterator seqs_begin,
                         RandomAccessIteratorPairIterator seqs_end,
                         RandomAccessIterator3 target, DiffType length,
                         Comparator comp)
{
#if STXXL_PARALLEL
    return stxxl::parallel::multiway_merge_sentinels(
        seqs_begin, seqs_end, target, length, comp);
#else
    return stxxl::parallel::sequential_multiway_merge<false, true>(
        seqs_begin, seqs_end, target, length, comp);
#endif
}

/*! Multi-way merging front-end.
 * \param seqs_begin Begin iterator of iterator pair input sequence.
 * \param seqs_end End iterator of iterator pair input sequence.
 * \param target Begin iterator out output sequence.
 * \param comp Comparator.
 * \param length Maximum length to merge.
 * \return End iterator of output sequence.
 * \pre For each \c i, \c seqs_begin[i].second must be the end marker of the sequence, but also reference the one more sentinel element.
 */
template <typename RandomAccessIteratorPairIterator,
          typename RandomAccessIterator3, typename DiffType, typename Comparator>
RandomAccessIterator3
multiway_merge_stable_sentinels(RandomAccessIteratorPairIterator seqs_begin,
                                RandomAccessIteratorPairIterator seqs_end,
                                RandomAccessIterator3 target, DiffType length,
                                Comparator comp)
{
#if STXXL_PARALLEL
    return stxxl::parallel::multiway_merge_stable_sentinels(
        seqs_begin, seqs_end, target, length, comp);
#else
    return stxxl::parallel::sequential_multiway_merge<true, true>(
        seqs_begin, seqs_end, target, length, comp);
#endif
}

} // namespace potentially_parallel

STXXL_END_NAMESPACE

#endif // !STXXL_PARALLEL_HEADER
// vim: et:ts=4:sw=4
