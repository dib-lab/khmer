/***************************************************************************
 *  include/stxxl/bits/parallel/merge.h
 *
 *  Parallel implementation of std::merge().
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

#ifndef STXXL_PARALLEL_MERGE_HEADER
#define STXXL_PARALLEL_MERGE_HEADER

#include <stxxl/bits/namespace.h>
#include <stxxl/bits/parallel/compiletime_settings.h>
#include <iterator>
#include <cassert>

STXXL_BEGIN_NAMESPACE

namespace parallel {

/*!
 * Merge routine being able to merge only the \c max_length smallest elements.
 *
 * The \c begin iterators are advanced accordingly, they might not reach \c
 * end, in contrast to the usual variant.
 *
 * \param begin1 Begin iterator of first sequence.
 * \param end1 End iterator of first sequence.
 * \param begin2 Begin iterator of second sequence.
 * \param end2 End iterator of second sequence.
 * \param target Target begin iterator.
 * \param max_length Maximum number of elements to merge.
 * \param comp Comparator.
 * \return Output end iterator.
 */
template <typename RandomAccessIterator1, typename RandomAccessIterator2,
          typename OutputIterator,
          typename DiffType, typename Comparator>
OutputIterator
merge_advance_usual(RandomAccessIterator1& begin1, RandomAccessIterator1 end1,
                    RandomAccessIterator2& begin2, RandomAccessIterator2 end2,
                    OutputIterator target, DiffType max_length,
                    Comparator comp)
{
    while (begin1 != end1 && begin2 != end2 && max_length > 0)
    {
        // array1[i1] < array0[i0]
        if (comp(*begin2, *begin1))
            *target++ = *begin2++;
        else
            *target++ = *begin1++;
        --max_length;
    }

    if (begin1 != end1)
    {
        target = std::copy(begin1, begin1 + max_length, target);
        begin1 += max_length;
    }
    else
    {
        target = std::copy(begin2, begin2 + max_length, target);
        begin2 += max_length;
    }
    return target;
}

/*!
 * Merge routine being able to merge only the \c max_length smallest elements.
 *
 * The \c begin iterators are advanced accordingly, they might not reach \c
 * end, in contrast to the usual variant.  Specially designed code should allow
 * the compiler to generate conditional moves instead of branches.
 *
 * \param begin1 Begin iterator of first sequence.
 * \param end1 End iterator of first sequence.
 * \param begin2 Begin iterator of second sequence.
 * \param end2 End iterator of second sequence.
 * \param target Target begin iterator.
 * \param max_length Maximum number of elements to merge.
 * \param comp Comparator.
 * \return Output end iterator.
 */
template <typename RandomAccessIterator1, typename RandomAccessIterator2,
          typename OutputIterator,
          typename DiffType, typename Comparator>
OutputIterator
merge_advance_movc(RandomAccessIterator1& begin1, RandomAccessIterator1 end1,
                   RandomAccessIterator2& begin2, RandomAccessIterator2 end2,
                   OutputIterator target,
                   DiffType max_length, Comparator comp)
{
    typedef typename std::iterator_traits<RandomAccessIterator1>::value_type ValueType1;
    typedef typename std::iterator_traits<RandomAccessIterator2>::value_type ValueType2;

    while (begin1 != end1 && begin2 != end2 && max_length > 0)
    {
        RandomAccessIterator1 next1 = begin1 + 1;
        RandomAccessIterator2 next2 = begin2 + 1;
        ValueType1 element1 = *begin1;
        ValueType2 element2 = *begin2;

        if (comp(element2, element1))
        {
            element1 = element2;
            begin2 = next2;
        }
        else
        {
            begin1 = next1;
        }

        *target = element1;

        ++target;
        --max_length;
    }

    if (begin1 != end1)
    {
        target = std::copy(begin1, begin1 + max_length, target);
        begin1 += max_length;
    }
    else
    {
        target = std::copy(begin2, begin2 + max_length, target);
        begin2 += max_length;
    }

    return target;
}

/*!
 * Merge routine being able to merge only the \c max_length smallest elements.
 *
 * The \c begin iterators are advanced accordingly, they might not reach \c
 * end, in contrast to the usual variant.  Static switch on whether to use the
 * conditional-move variant.
 *
 * \param begin1 Begin iterator of first sequence.
 * \param end1 End iterator of first sequence.
 * \param begin2 Begin iterator of second sequence.
 * \param end2 End iterator of second sequence.
 * \param target Target begin iterator.
 * \param max_length Maximum number of elements to merge.
 * \param comp Comparator.
 * \return Output end iterator.
 */
template <typename RandomAccessIterator1, typename RandomAccessIterator2,
          typename OutputIterator,
          typename DiffType, typename Comparator>
OutputIterator
merge_advance(RandomAccessIterator1& begin1, RandomAccessIterator1 end1,
              RandomAccessIterator2& begin2, RandomAccessIterator2 end2,
              OutputIterator target,
              DiffType max_length, Comparator comp)
{
    STXXL_PARALLEL_PCALL(max_length);

    return merge_advance_movc(begin1, end1, begin2, end2, target, max_length, comp);
}

/*!
 * Merge routine fallback to sequential in case the iterators of the two input
 * sequences are of different type.
 *
 * \param begin1 Begin iterator of first sequence.
 * \param end1 End iterator of first sequence.
 * \param begin2 Begin iterator of second sequence.
 * \param end2 End iterator of second sequence.
 * \param target Target begin iterator.
 * \param max_length Maximum number of elements to merge.
 * \param comp Comparator.
 * \return Output end iterator.
 */
template <typename RandomAccessIterator1, typename RandomAccessIterator2,
          typename RandomAccessIterator3, typename Comparator>
RandomAccessIterator3
parallel_merge_advance(
    RandomAccessIterator1& begin1, RandomAccessIterator1 end1,
    RandomAccessIterator2& begin2, RandomAccessIterator2 end2,
    // different iterators, parallel implementation not available
    RandomAccessIterator3 target,
    typename std::iterator_traits<RandomAccessIterator1>::difference_type max_length,
    Comparator comp)
{
    return merge_advance(begin1, end1, begin2, end2, target, max_length, comp);
}

/*!
 * Parallel merge routine being able to merge only the \c max_length smallest
 * elements.
 *
 * The \c begin iterators are advanced accordingly, they might not reach \c
 * end, in contrast to the usual variant.  The functionality is projected onto
 * parallel_multiway_merge.
 *
 * \param begin1 Begin iterator of first sequence.
 * \param end1 End iterator of first sequence.
 * \param begin2 Begin iterator of second sequence.
 * \param end2 End iterator of second sequence.
 * \param target Target begin iterator.
 * \param max_length Maximum number of elements to merge.
 * \param comp Comparator.
 * \return Output end iterator.
 */
template <typename RandomAccessIterator1, typename RandomAccessIterator3,
          typename Comparator>
RandomAccessIterator3
parallel_merge_advance(
    RandomAccessIterator1& begin1, RandomAccessIterator1 end1,
    RandomAccessIterator1& begin2, RandomAccessIterator1 end2,
    RandomAccessIterator3 target,
    typename std::iterator_traits<RandomAccessIterator1>::difference_type max_length,
    Comparator comp)
{
    std::pair<RandomAccessIterator1, RandomAccessIterator1> seqs[2] = {
        std::make_pair(begin1, end1), std::make_pair(begin2, end2)
    };
    RandomAccessIterator3 target_end = parallel_multiway_merge(
        seqs, seqs + 2, target, comp, max_length, true, false
        );

    return target_end;
}

} // namespace parallel

STXXL_END_NAMESPACE

#endif // !STXXL_PARALLEL_MERGE_HEADER
