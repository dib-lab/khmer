/***************************************************************************
 *  include/stxxl/bits/parallel/multiway_merge.h
 *
 *  Implementation of sequential and parallel multiway merge.
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

#ifndef STXXL_PARALLEL_MULTIWAY_MERGE_HEADER
#define STXXL_PARALLEL_MULTIWAY_MERGE_HEADER

#include <vector>
#include <iterator>
#include <algorithm>

#include <stxxl/bits/verbose.h>
#include <stxxl/bits/common/is_sorted.h>
#include <stxxl/bits/common/utils.h>
#include <stxxl/bits/parallel/merge.h>
#include <stxxl/bits/parallel/losertree.h>
#include <stxxl/bits/parallel/settings.h>
#include <stxxl/bits/parallel/equally_split.h>
#include <stxxl/bits/parallel/multiseq_selection.h>
#include <stxxl/bits/parallel/timing.h>
#include <stxxl/bits/parallel/tags.h>

#if defined(_MSC_VER) && STXXL_DEBUG_ASSERTIONS
#include <BaseTsd.h>
typedef SSIZE_T ssize_t;
#endif

STXXL_BEGIN_NAMESPACE

namespace parallel {

//! Length of a sequence described by a pair of iterators.
template <typename RandomAccessIteratorPair>
typename std::iterator_traits<
    typename RandomAccessIteratorPair::first_type
    >::difference_type
iterpair_size(const RandomAccessIteratorPair& p)
{
    return p.second - p.first;
}

/*!
 * Iterator wrapper supporting an implicit supremum at the end of the sequence,
 * dominating all comparisons.  Deriving from RandomAccessIterator is not
 * possible since RandomAccessIterator need not be a class.
*/
template <typename RandomAccessIterator, typename Comparator>
class guarded_iterator
{
public:
    //! Our own type
    typedef guarded_iterator<RandomAccessIterator, Comparator> self_type;

    //! Value type of the iterator
    typedef typename std::iterator_traits<RandomAccessIterator>::value_type value_type;

protected:
    //! Current iterator position.
    RandomAccessIterator current;
    //! End iterator of the sequence.
    RandomAccessIterator end;
    //! Comparator.
    Comparator& comp;

public:
    /*!
     * Constructor. Sets iterator to beginning of sequence.
     * \param begin Begin iterator of sequence.
     * \param end End iterator of sequence.
     * \param comp Comparator provided for associated overloaded compare
     * operators.
     */
    guarded_iterator(RandomAccessIterator begin, RandomAccessIterator end,
                     Comparator& comp)
        : current(begin), end(end), comp(comp)
    { }

    /*!
     * Pre-increment operator.
     * \return This.
     */
    self_type& operator ++ ()
    {
        ++current;
        return *this;
    }

    /*!
     * Dereference operator.
     * \return Referenced element.
     */
    value_type& operator * ()
    {
        return *current;
    }

    /*!
     * Convert to wrapped iterator.
     * \return Wrapped iterator.
     */
    RandomAccessIterator & iterator()
    {
        return current;
    }

    /*!
     * Compare two elements referenced by guarded iterators.
     * \param bi1 First iterator.
     * \param bi2 Second iterator.
     * \return \c True if less.
     */
    friend bool operator < (self_type& bi1, self_type& bi2)
    {
        if (bi1.current == bi1.end)             // bi1 is sup
            return bi2.current == bi2.end;      // bi2 is not sup
        if (bi2.current == bi2.end)             // bi2 is sup
            return true;
        return bi1.comp(*bi1, *bi2);            // normal compare
    }

    /*!
     * Compare two elements referenced by guarded iterators.
     * \param bi1 First iterator.
     * \param bi2 Second iterator.
     * \return \c True if less equal.
     */
    friend bool operator <= (self_type& bi1, self_type& bi2)
    {
        if (bi2.current == bi2.end)         //bi1 is sup
            return bi1.current != bi1.end;  //bi2 is not sup
        if (bi1.current == bi1.end)         //bi2 is sup
            return false;
        return !bi1.comp(*bi2, *bi1);       //normal compare
    }
};

template <typename RandomAccessIterator, typename Comparator>
class unguarded_iterator
{
public:
    //! Our own type
    typedef unguarded_iterator<RandomAccessIterator, Comparator> self_type;

    //! Value type of the iterator
    typedef typename std::iterator_traits<RandomAccessIterator>::value_type value_type;

protected:
    //! Current iterator position.
    RandomAccessIterator current;
    //! Comparator.
    Comparator& comp;

public:
    /*!
     * Constructor. Sets iterator to beginning of sequence.
     * \param begin Begin iterator of sequence.
     * param end Unused, only for compatibility.
     * \param comp Unused, only for compatibility.
     */
    unguarded_iterator(RandomAccessIterator begin,
                       RandomAccessIterator /* end */,
                       Comparator& comp)
        : current(begin), comp(comp)
    { }

    /*!
     * Pre-increment operator.
     * \return This.
     */
    self_type& operator ++ ()
    {
        ++current;
        return *this;
    }

    /*!
     * Dereference operator.
     * \return Referenced element.
     */
    value_type& operator * ()
    {
        return *current;
    }

    /*!
     * Convert to wrapped iterator.
     * \return Wrapped iterator.
     */
    RandomAccessIterator & iterator()
    {
        return current;
    }

    /*!
     * Compare two elements referenced by unguarded iterators.
     * \param bi1 First iterator.
     * \param bi2 Second iterator.
     * \return \c True if less.
     */
    friend bool operator < (self_type& bi1, self_type& bi2)
    {
        return bi1.comp(*bi1, *bi2);    // normal compare, unguarded
    }

    /*!
     * Compare two elements referenced by unguarded iterators.
     * \param bi1 First iterator.
     * \param bi2 Second iterator.
     * \return \c True if less equal.
     */
    friend bool operator <= (self_type& bi1, self_type& bi2)
    {
        return !bi1.comp(*bi2, *bi1);   // normal compare, unguarded
    }
};

/*!
 * Prepare a set of sequences to be merged without a (end) guard
 *
 * \param seqs_begin
 * \param seqs_end
 * \param comp
 * \param min_sequence
 * \tparam Stable
 * \pre (seqs_end - seqs_begin > 0)
 */
template <bool Stable, typename RandomAccessIteratorIterator, typename Comparator>
typename std::iterator_traits<
    typename std::iterator_traits<RandomAccessIteratorIterator>::value_type::first_type
    >::difference_type
prepare_unguarded(RandomAccessIteratorIterator seqs_begin,
                  RandomAccessIteratorIterator seqs_end,
                  Comparator comp,
                  int& min_sequence)
{
    STXXL_PARALLEL_PCALL(seqs_end - seqs_begin);

    typedef typename std::iterator_traits<RandomAccessIteratorIterator>
        ::value_type::first_type RandomAccessIterator;
    typedef typename std::iterator_traits<RandomAccessIterator>
        ::value_type value_type;
    typedef typename std::iterator_traits<RandomAccessIterator>
        ::difference_type diff_type;

    if ((*seqs_begin).first == (*seqs_begin).second)
    {
        // empty sequence found, it's the first one
        min_sequence = 0;
        return -1;
    }

    // last element in sequence
    value_type min = *((*seqs_begin).second - 1);
    min_sequence = 0;
    for (RandomAccessIteratorIterator s = seqs_begin + 1; s != seqs_end; ++s)
    {
        if ((*s).first == (*s).second)
        {
            // empty sequence found
            min_sequence = static_cast<int>(s - seqs_begin);
            return -1;
        }
        const value_type& v = *((*s).second - 1);
        if (comp(v, min))
        {
            // last element in sequence is strictly smaller
            min = v;
            min_sequence = static_cast<int>(s - seqs_begin);
        }
    }

    diff_type overhang_size = 0;

    int s = 0;
    for (s = 0; s <= min_sequence; ++s)
    {
        RandomAccessIterator split;
        if (Stable)
            split = std::upper_bound(seqs_begin[s].first, seqs_begin[s].second,
                                     min, comp);
        else
            split = std::lower_bound(seqs_begin[s].first, seqs_begin[s].second,
                                     min, comp);

        overhang_size += seqs_begin[s].second - split;
    }

    for ( ; s < (seqs_end - seqs_begin); ++s)
    {
        RandomAccessIterator split =
            std::lower_bound(seqs_begin[s].first, seqs_begin[s].second,
                             min, comp);
        overhang_size += seqs_begin[s].second - split;
    }

    return overhang_size;       // so many elements will be left over afterwards
}

/*!
 * Prepare a set of sequences to be merged with a (end) guard (sentinel)
 * \param seqs_begin
 * \param seqs_end
 * \param comp
 */
template <typename RandomAccessIteratorIterator, typename Comparator>
typename std::iterator_traits<
    typename std::iterator_traits<RandomAccessIteratorIterator>::value_type::first_type
    >::difference_type
prepare_unguarded_sentinel(RandomAccessIteratorIterator seqs_begin,
                           RandomAccessIteratorIterator seqs_end,
                           Comparator comp)
{
    STXXL_PARALLEL_PCALL(seqs_end - seqs_begin);

    typedef typename std::iterator_traits<RandomAccessIteratorIterator>
        ::value_type::first_type RandomAccessIterator;
    typedef typename std::iterator_traits<RandomAccessIterator>
        ::value_type value_type;
    typedef typename std::iterator_traits<RandomAccessIterator>
        ::difference_type diff_type;

    value_type* max_value = NULL;   // last element in sequence
    for (RandomAccessIteratorIterator s = seqs_begin; s != seqs_end; ++s)
    {
        if ((*s).first == (*s).second)
            continue;
        value_type& v = *((*s).second - 1);     //last element in sequence
        if (!max_value || comp(*max_value, v))  //strictly greater
            max_value = &v;
    }

    diff_type overhang_size = 0;

    for (RandomAccessIteratorIterator s = seqs_begin; s != seqs_end; ++s)
    {
        RandomAccessIterator split = std::lower_bound((*s).first, (*s).second, *max_value, comp);
        overhang_size += (*s).second - split;
        *((*s).second) = *max_value; //set sentinel
    }

    return overhang_size;            // so many elements will be left over afterwards
}

/*!
 * Highly efficient 3-way merging procedure.
 *
 * Merging is done with the algorithm implementation described by Peter
 * Sanders.  Basically, the idea is to minimize the number of necessary
 * comparison after merging an element.  The implementation trick that makes
 * this fast is that the order of the sequences is stored in the instruction
 * pointer (translated into labels in C++).
 *
 * This works well for merging up to 4 sequences.
 *
 * Note that making the merging stable does \a not come at a performance hit.
 *
 * Whether the merging is done guarded or unguarded is selected by the used
 * iterator class.
 *
 * \param seqs_begin Begin iterator of iterator pair input sequence.
 * \param seqs_end End iterator of iterator pair input sequence.
 * \param target Begin iterator out output sequence.
 * \param length Maximum length to merge.
 * \param comp Comparator.
 * \return End iterator of output sequence.
 */
template <template <typename RAI, typename C> class Iterator,
          typename RandomAccessIteratorIterator,
          typename RandomAccessIterator3,
          typename DiffType, typename Comparator>
RandomAccessIterator3
multiway_merge_3_variant(RandomAccessIteratorIterator seqs_begin,
                         RandomAccessIteratorIterator seqs_end,
                         RandomAccessIterator3 target, DiffType length,
                         Comparator comp)
{
    STXXL_PARALLEL_PCALL(length);
    STXXL_ASSERT(seqs_end - seqs_begin == 3);

    typedef typename std::iterator_traits<RandomAccessIteratorIterator>
        ::value_type::first_type RandomAccessIterator;

    if (length == 0)
        return target;

#if STXXL_DEBUG_ASSERTIONS
    ssize_t orig_length = length;
#endif

    Iterator<RandomAccessIterator, Comparator>
    seq0(seqs_begin[0].first, seqs_begin[0].second, comp),
    seq1(seqs_begin[1].first, seqs_begin[1].second, comp),
    seq2(seqs_begin[2].first, seqs_begin[2].second, comp);

    if (seq0 <= seq1)
    {
        if (seq1 <= seq2)
            goto s012;
        else if (seq2 < seq0)
            goto s201;
        else
            goto s021;
    }
    else
    {
        if (seq1 <= seq2)
        {
            if (seq0 <= seq2)
                goto s102;
            else
                goto s120;
        }
        else
            goto s210;
    }

#define STXXL_MERGE3CASE(a, b, c, c0, c1)            \
    s ## a ## b ## c:                                \
    *target = *seq ## a;                             \
    ++target;                                        \
    --length;                                        \
    ++seq ## a;                                      \
    if (length == 0) goto finish;                    \
    if (seq ## a c0 seq ## b) goto s ## a ## b ## c; \
    if (seq ## a c1 seq ## c) goto s ## b ## a ## c; \
    goto s ## b ## c ## a;

    STXXL_MERGE3CASE(0, 1, 2, <=, <=);
    STXXL_MERGE3CASE(1, 2, 0, <=, <);
    STXXL_MERGE3CASE(2, 0, 1, <, <);
    STXXL_MERGE3CASE(1, 0, 2, <, <=);
    STXXL_MERGE3CASE(0, 2, 1, <=, <=);
    STXXL_MERGE3CASE(2, 1, 0, <, <);

#undef STXXL_MERGE3CASE

finish:
    ;

#if STXXL_DEBUG_ASSERTIONS
    STXXL_CHECK_EQUAL((seq0.iterator() - seqs_begin[0].first) +
                      (seq1.iterator() - seqs_begin[1].first) +
                      (seq2.iterator() - seqs_begin[2].first),
                      orig_length);
#endif

    seqs_begin[0].first = seq0.iterator();
    seqs_begin[1].first = seq1.iterator();
    seqs_begin[2].first = seq2.iterator();

    return target;
}

template <typename RandomAccessIteratorIterator,
          typename RandomAccessIterator3,
          typename DiffType, typename Comparator>
RandomAccessIterator3
multiway_merge_3_combined(RandomAccessIteratorIterator seqs_begin,
                          RandomAccessIteratorIterator seqs_end,
                          RandomAccessIterator3 target, DiffType length,
                          Comparator comp)
{
    STXXL_PARALLEL_PCALL(length);
    STXXL_ASSERT(seqs_end - seqs_begin == 3);

    int min_seq;
    RandomAccessIterator3 target_end;
    DiffType overhang = prepare_unguarded<true>(seqs_begin, seqs_end, comp, min_seq);

    DiffType total_length = 0;
    for (RandomAccessIteratorIterator s = seqs_begin; s != seqs_end; ++s)
        total_length += iterpair_size(*s);

    if (overhang != (DiffType)(-1))
    {
        DiffType unguarded_length = std::min(length, total_length - overhang);
        target_end = multiway_merge_3_variant<unguarded_iterator>
                         (seqs_begin, seqs_end, target, unguarded_length, comp);
        overhang = length - unguarded_length;
    }
    else
    {
        // empty sequence found
        overhang = length;
        target_end = target;
    }

    STXXL_DEBUG_ASSERT(target_end == target + length - overhang);
    STXXL_DEBUG_ASSERT(stxxl::is_sorted(target, target_end, comp));

    switch (min_seq)
    {
    case 0:
        // iterators will be advanced accordingly
        target_end = merge_advance(
            seqs_begin[1].first, seqs_begin[1].second,
            seqs_begin[2].first, seqs_begin[2].second,
            target_end, overhang, comp);
        break;
    case 1:
        target_end = merge_advance(
            seqs_begin[0].first, seqs_begin[0].second,
            seqs_begin[2].first, seqs_begin[2].second,
            target_end, overhang, comp);
        break;
    case 2:
        target_end = merge_advance(
            seqs_begin[0].first, seqs_begin[0].second,
            seqs_begin[1].first, seqs_begin[1].second,
            target_end, overhang, comp);
        break;
    default:
        assert(false);
    }

    STXXL_DEBUG_ASSERT(target_end == target + length);
    STXXL_DEBUG_ASSERT(stxxl::is_sorted(target, target_end, comp));

    return target_end;
}

/*!
 * Highly efficient 4-way merging procedure.
 *
 * Merging is done with the algorithm implementation described by Peter
 * Sanders. Basically, the idea is to minimize the number of necessary
 * comparison after merging an element.  The implementation trick that makes
 * this fast is that the order of the sequences is stored in the instruction
 * pointer (translated into goto labels in C++).
 *
 * This works well for merging up to 4 sequences.
 *
 * Note that making the merging stable does \a not come at a performance hit.
 *
 * Whether the merging is done guarded or unguarded is selected by the used
 * iterator class.
 *
 * \param seqs_begin Begin iterator of iterator pair input sequence.
 * \param seqs_end End iterator of iterator pair input sequence.
 * \param target Begin iterator out output sequence.
 * \param length Maximum length to merge.
 * \param comp Comparator.
 * \return End iterator of output sequence.
 */
template <template <typename RAI, typename C> class iterator,
          typename RandomAccessIteratorIterator,
          typename RandomAccessIterator3,
          typename DiffType, typename Comparator>
RandomAccessIterator3
multiway_merge_4_variant(RandomAccessIteratorIterator seqs_begin,
                         RandomAccessIteratorIterator seqs_end,
                         RandomAccessIterator3 target, DiffType length,
                         Comparator comp)
{
    STXXL_PARALLEL_PCALL(length);
    STXXL_ASSERT(seqs_end - seqs_begin == 4);

    typedef typename std::iterator_traits<RandomAccessIteratorIterator>
        ::value_type::first_type RandomAccessIterator;

    if (length == 0)
        return target;

#if STXXL_DEBUG_ASSERTIONS
    ssize_t orig_length = length;
#endif

    iterator<RandomAccessIterator, Comparator>
    seq0(seqs_begin[0].first, seqs_begin[0].second, comp),
    seq1(seqs_begin[1].first, seqs_begin[1].second, comp),
    seq2(seqs_begin[2].first, seqs_begin[2].second, comp),
    seq3(seqs_begin[3].first, seqs_begin[3].second, comp);

#define STXXL_DECISION(a, b, c, d) do {                      \
        if (seq ## d < seq ## a) goto s ## d ## a ## b ## c; \
        if (seq ## d < seq ## b) goto s ## a ## d ## b ## c; \
        if (seq ## d < seq ## c) goto s ## a ## b ## d ## c; \
        goto s ## a ## b ## c ## d;                          \
}                                                            \
    while (0)

    if (seq0 <= seq1)
    {
        if (seq1 <= seq2)
            STXXL_DECISION(0, 1, 2, 3);
        else if (seq2 < seq0)
            STXXL_DECISION(2, 0, 1, 3);
        else
            STXXL_DECISION(0, 2, 1, 3);
    }
    else
    {
        if (seq1 <= seq2)
        {
            if (seq0 <= seq2)
                STXXL_DECISION(1, 0, 2, 3);
            else
                STXXL_DECISION(1, 2, 0, 3);
        }
        else
            STXXL_DECISION(2, 1, 0, 3);
    }

#define STXXL_MERGE4CASE(a, b, c, d, c0, c1, c2)          \
    s ## a ## b ## c ## d:                                \
    if (length == 0) goto finish;                         \
    *target = *seq ## a;                                  \
    ++target;                                             \
    --length;                                             \
    ++seq ## a;                                           \
    if (seq ## a c0 seq ## b) goto s ## a ## b ## c ## d; \
    if (seq ## a c1 seq ## c) goto s ## b ## a ## c ## d; \
    if (seq ## a c2 seq ## d) goto s ## b ## c ## a ## d; \
    goto s ## b ## c ## d ## a;

    STXXL_MERGE4CASE(0, 1, 2, 3, <=, <=, <=);
    STXXL_MERGE4CASE(0, 1, 3, 2, <=, <=, <=);
    STXXL_MERGE4CASE(0, 2, 1, 3, <=, <=, <=);
    STXXL_MERGE4CASE(0, 2, 3, 1, <=, <=, <=);
    STXXL_MERGE4CASE(0, 3, 1, 2, <=, <=, <=);
    STXXL_MERGE4CASE(0, 3, 2, 1, <=, <=, <=);
    STXXL_MERGE4CASE(1, 0, 2, 3, <, <=, <=);
    STXXL_MERGE4CASE(1, 0, 3, 2, <, <=, <=);
    STXXL_MERGE4CASE(1, 2, 0, 3, <=, <, <=);
    STXXL_MERGE4CASE(1, 2, 3, 0, <=, <=, <);
    STXXL_MERGE4CASE(1, 3, 0, 2, <=, <, <=);
    STXXL_MERGE4CASE(1, 3, 2, 0, <=, <=, <);
    STXXL_MERGE4CASE(2, 0, 1, 3, <, <, <=);
    STXXL_MERGE4CASE(2, 0, 3, 1, <, <=, <);
    STXXL_MERGE4CASE(2, 1, 0, 3, <, <, <=);
    STXXL_MERGE4CASE(2, 1, 3, 0, <, <=, <);
    STXXL_MERGE4CASE(2, 3, 0, 1, <=, <, <);
    STXXL_MERGE4CASE(2, 3, 1, 0, <=, <, <);
    STXXL_MERGE4CASE(3, 0, 1, 2, <, <, <);
    STXXL_MERGE4CASE(3, 0, 2, 1, <, <, <);
    STXXL_MERGE4CASE(3, 1, 0, 2, <, <, <);
    STXXL_MERGE4CASE(3, 1, 2, 0, <, <, <);
    STXXL_MERGE4CASE(3, 2, 0, 1, <, <, <);
    STXXL_MERGE4CASE(3, 2, 1, 0, <, <, <);

#undef STXXL_MERGE4CASE
#undef STXXL_DECISION

finish:
    ;

#if STXXL_DEBUG_ASSERTIONS
    STXXL_CHECK_EQUAL((seq0.iterator() - seqs_begin[0].first) +
                      (seq1.iterator() - seqs_begin[1].first) +
                      (seq2.iterator() - seqs_begin[2].first) +
                      (seq3.iterator() - seqs_begin[3].first),
                      orig_length);
#endif

    seqs_begin[0].first = seq0.iterator();
    seqs_begin[1].first = seq1.iterator();
    seqs_begin[2].first = seq2.iterator();
    seqs_begin[3].first = seq3.iterator();

    return target;
}

template <typename RandomAccessIteratorIterator,
          typename RandomAccessIterator3,
          typename DiffType, typename Comparator>
RandomAccessIterator3
multiway_merge_4_combined(RandomAccessIteratorIterator seqs_begin,
                          RandomAccessIteratorIterator seqs_end,
                          RandomAccessIterator3 target, DiffType length,
                          Comparator comp)
{
    STXXL_PARALLEL_PCALL(length);
    STXXL_ASSERT(seqs_end - seqs_begin == 4);

    typedef typename std::iterator_traits<RandomAccessIteratorIterator>
        ::value_type RandomAccessIteratorPair;

    int min_seq;
    RandomAccessIterator3 target_end;
    DiffType overhang = prepare_unguarded<true>(seqs_begin, seqs_end, comp, min_seq);

    DiffType total_length = 0;
    for (RandomAccessIteratorIterator s = seqs_begin; s != seqs_end; ++s)
        total_length += iterpair_size(*s);

    if (overhang != (DiffType) - 1)
    {
        DiffType unguarded_length = std::min(length, total_length - overhang);
        target_end = multiway_merge_4_variant<unguarded_iterator>
                         (seqs_begin, seqs_end, target, unguarded_length, comp);
        overhang = length - unguarded_length;
    }
    else
    {
        // empty sequence found
        overhang = length;
        target_end = target;
    }

    STXXL_DEBUG_ASSERT(target_end == target + length - overhang);
    STXXL_DEBUG_ASSERT(stxxl::is_sorted(target, target_end, comp));

    std::vector<RandomAccessIteratorPair> one_missing(seqs_begin, seqs_end);
    one_missing.erase(one_missing.begin() + min_seq);                                               //remove

    target_end = multiway_merge_3_variant<guarded_iterator>(one_missing.begin(), one_missing.end(), target_end, overhang, comp);

    one_missing.insert(one_missing.begin() + min_seq, seqs_begin[min_seq]);                         //insert back again
    std::copy(one_missing.begin(), one_missing.end(), seqs_begin);                                  //write back modified iterators

    STXXL_DEBUG_ASSERT(target_end == target + length);
    STXXL_DEBUG_ASSERT(stxxl::is_sorted(target, target_end, comp));

    return target_end;
}

/*!
 * Basic multi-way merging procedure.
 *
 * The head elements are kept in a sorted array, new heads are inserted
 * linearly.
 *
 * \param seqs_begin Begin iterator of iterator pair input sequence.
 * \param seqs_end End iterator of iterator pair input sequence.
 * \param target Begin iterator out output sequence.
 * \param length Maximum length to merge.
 * \param comp Comparator.
 * \tparam Stable Stable merging incurs a performance penalty.
 * \return End iterator of output sequence.
 */
template <bool Stable,
          typename RandomAccessIteratorIterator,
          typename RandomAccessIterator3,
          typename DiffType, typename Comparator>
RandomAccessIterator3
multiway_merge_bubble(RandomAccessIteratorIterator seqs_begin,
                      RandomAccessIteratorIterator seqs_end,
                      RandomAccessIterator3 target, DiffType length,
                      Comparator comp)
{
    STXXL_PARALLEL_PCALL(length);

    typedef typename std::iterator_traits<RandomAccessIteratorIterator>
        ::value_type::first_type RandomAccessIterator;
    typedef typename std::iterator_traits<RandomAccessIterator>
        ::value_type value_type;

    // num remaining pieces
    int k = static_cast<int>(seqs_end - seqs_begin), nrp;

    value_type* pl = new value_type[k];
    int* source = new int[k];
    DiffType total_length = 0;

#define POS(i) seqs_begin[(i)].first
#define STOPS(i) seqs_begin[(i)].second

    //write entries into queue
    nrp = 0;
    for (int pi = 0; pi < k; ++pi)
    {
        if (STOPS(pi) != POS(pi))
        {
            pl[nrp] = *(POS(pi));
            source[nrp] = pi;
            ++nrp;
            total_length += iterpair_size(seqs_begin[pi]);
        }
    }

    if (Stable)
    {
        for (int k = 0; k < nrp - 1; ++k)
            for (int pi = nrp - 1; pi > k; --pi)
                if (comp(pl[pi], pl[pi - 1]) ||
                    (!comp(pl[pi - 1], pl[pi]) && source[pi] < source[pi - 1]))
                {
                    std::swap(pl[pi - 1], pl[pi]);
                    std::swap(source[pi - 1], source[pi]);
                }
    }
    else
    {
        for (int k = 0; k < nrp - 1; ++k)
            for (int pi = nrp - 1; pi > k; --pi)
                if (comp(pl[pi], pl[pi - 1]))
                {
                    std::swap(pl[pi - 1], pl[pi]);
                    std::swap(source[pi - 1], source[pi]);
                }
    }

    // iterate
    if (Stable)
    {
        int j;
        while (nrp > 0 && length > 0)
        {
            if (source[0] < source[1])
            {
                // pl[0] <= pl[1] ?
                while ((nrp == 1 || !(comp(pl[1], pl[0]))) && length > 0)
                {
                    *target = pl[0];
                    ++target;
                    ++POS(source[0]);
                    --length;
                    if (POS(source[0]) == STOPS(source[0]))
                    {
                        // move everything to the left
                        for (int s = 0; s < nrp - 1; ++s)
                        {
                            pl[s] = pl[s + 1];
                            source[s] = source[s + 1];
                        }
                        --nrp;
                        break;
                    }
                    else
                        pl[0] = *(POS(source[0]));
                }
            }
            else
            {
                // pl[0] < pl[1] ?
                while ((nrp == 1 || comp(pl[0], pl[1])) && length > 0)
                {
                    *target = pl[0];
                    ++target;
                    ++POS(source[0]);
                    --length;
                    if (POS(source[0]) == STOPS(source[0]))
                    {
                        for (int s = 0; s < nrp - 1; ++s)
                        {
                            pl[s] = pl[s + 1];
                            source[s] = source[s + 1];
                        }
                        --nrp;
                        break;
                    }
                    else
                        pl[0] = *(POS(source[0]));
                }
            }

            //sink down
            j = 1;
            while ((j < nrp) && (comp(pl[j], pl[j - 1]) ||
                                 (!comp(pl[j - 1], pl[j]) && (source[j] < source[j - 1]))))
            {
                std::swap(pl[j - 1], pl[j]);
                std::swap(source[j - 1], source[j]);
                ++j;
            }
        }
    }
    else
    {
        int j;
        while (nrp > 0 && length > 0)
        {
            // pl[0] <= pl[1] ?
            while ((nrp == 1 || !comp(pl[1], pl[0])) && length > 0)
            {
                *target = pl[0];
                ++target;
                ++POS(source[0]);
                --length;
                if (POS(source[0]) == STOPS(source[0]))
                {
                    for (int s = 0; s < (nrp - 1); ++s)
                    {
                        pl[s] = pl[s + 1];
                        source[s] = source[s + 1];
                    }
                    --nrp;
                    break;
                }
                else
                    pl[0] = *(POS(source[0]));
            }

            //sink down
            j = 1;
            while ((j < nrp) && comp(pl[j], pl[j - 1]))
            {
                std::swap(pl[j - 1], pl[j]);
                std::swap(source[j - 1], source[j]);
                ++j;
            }
        }
    }

    delete[] pl;
    delete[] source;

    return target;
}

/*!
 * Multi-way merging procedure for a high branching factor, guarded case.
 *
 * The head elements are kept in a loser tree.
 * \param seqs_begin Begin iterator of iterator pair input sequence.
 * \param seqs_end End iterator of iterator pair input sequence.
 * \param target Begin iterator out output sequence.
 * \param length Maximum length to merge.
 * \param comp Comparator.
 * \tparam Stable Stable merging incurs a performance penalty.
 * \return End iterator of output sequence.
 */
template <typename LoserTreeType,
          typename RandomAccessIteratorIterator,
          typename RandomAccessIterator3,
          typename DiffType, typename Comparator>
RandomAccessIterator3
multiway_merge_loser_tree(RandomAccessIteratorIterator seqs_begin,
                          RandomAccessIteratorIterator seqs_end,
                          RandomAccessIterator3 target, DiffType length,
                          Comparator comp)
{
    STXXL_PARALLEL_PCALL(length);

    typedef typename LoserTreeType::source_type source_type;
    typedef typename std::iterator_traits<RandomAccessIteratorIterator>
        ::value_type::first_type RandomAccessIterator;
    typedef typename std::iterator_traits<RandomAccessIterator>
        ::value_type value_type;

    source_type k = static_cast<source_type>(seqs_end - seqs_begin);

    LoserTreeType lt(k, comp);

    DiffType total_length = 0;

    const value_type* arbitrary_element = NULL;

    // find an arbitrary element to avoid default construction
    for (source_type t = 0; t < k; ++t)
    {
        if (!arbitrary_element && iterpair_size(seqs_begin[t]) > 0)
            arbitrary_element = &(*seqs_begin[t].first);

        total_length += iterpair_size(seqs_begin[t]);
    }

    for (source_type t = 0; t < k; ++t)
    {
        if (UNLIKELY(seqs_begin[t].first == seqs_begin[t].second))
            lt.insert_start(*arbitrary_element, t, true);
        else
            lt.insert_start(*seqs_begin[t].first, t, false);
    }

    lt.init();

    total_length = std::min(total_length, length);

    for (DiffType i = 0; i < total_length; ++i)
    {
        // take out
        source_type source = lt.get_min_source();

        *target = *seqs_begin[source].first;
        ++target;
        ++seqs_begin[source].first;

        // feed
        if (seqs_begin[source].first == seqs_begin[source].second)
            lt.delete_min_insert(*arbitrary_element, true);
        else
            // replace from same source
            lt.delete_min_insert(*seqs_begin[source].first, false);
    }

    return target;
}

/*!
 * Multi-way merging procedure for a high branching factor, unguarded case.
 * The head elements are kept in a loser tree.
 *
 * \param seqs_begin Begin iterator of iterator pair input sequence.
 * \param seqs_end End iterator of iterator pair input sequence.
 * \param target Begin iterator out output sequence.
 * \param length Maximum length to merge.
 * \param comp Comparator.
 * \tparam Stable Stable merging incurs a performance penalty.
 * \return End iterator of output sequence.
 * \pre No input will run out of elements during the merge.
 */
template <typename LoserTreeType,
          typename RandomAccessIteratorIterator,
          typename RandomAccessIterator3,
          typename DiffType,
          typename Comparator>
RandomAccessIterator3
multiway_merge_loser_tree_unguarded(
    RandomAccessIteratorIterator seqs_begin,
    RandomAccessIteratorIterator seqs_end,
    RandomAccessIterator3 target, DiffType length,
    Comparator comp)
{
    STXXL_PARALLEL_PCALL(length);

    int k = (int)(seqs_end - seqs_begin);

    // sentinel is item at end of first sequence.
    LoserTreeType lt(k, *(seqs_begin->second - 1), comp);

    DiffType total_length = 0;

    for (int t = 0; t < k; ++t)
    {
        assert(seqs_begin[t].first != seqs_begin[t].second);

        lt.insert_start(*seqs_begin[t].first, t);

        total_length += iterpair_size(seqs_begin[t]);
    }

    lt.init();

    // do not go past end
    length = std::min(total_length, length);

    int source;

#if STXXL_DEBUG_ASSERTIONS
    DiffType i = 0;
#endif

    RandomAccessIterator3 target_end = target + length;
    while (target < target_end)
    {
        // take out
        source = lt.get_min_source();

#if STXXL_DEBUG_ASSERTIONS
        assert(i == 0 || !comp(*(seqs_begin[source].first), *(target - 1)));
#endif

        *target = *seqs_begin[source].first;
        ++seqs_begin[source].first;
        ++target;

#if STXXL_DEBUG_ASSERTIONS
        assert((seqs_begin[source].first != seqs_begin[source].second) || (i == length - 1));
        ++i;
#endif
        // feed
        // replace from same source
        lt.delete_min_insert(*seqs_begin[source].first);
    }

    return target;
}

template <bool Stable, class ValueType, class Comparator>
struct loser_tree_traits
{
public:
    typedef LoserTreePointer<Stable, ValueType, Comparator> LT;
};

#define STXXL_NO_POINTER(T)                              \
    template <bool Stable, class Comparator>             \
    struct loser_tree_traits<Stable, T, Comparator>      \
    {                                                    \
        typedef LoserTreeCopy<Stable, T, Comparator> LT; \
    };

STXXL_NO_POINTER(unsigned char)
STXXL_NO_POINTER(char)
STXXL_NO_POINTER(unsigned short)
STXXL_NO_POINTER(short)
STXXL_NO_POINTER(unsigned int)
STXXL_NO_POINTER(int)
STXXL_NO_POINTER(unsigned long)
STXXL_NO_POINTER(long)
STXXL_NO_POINTER(unsigned long long)
STXXL_NO_POINTER(long long)

#undef STXXL_NO_POINTER

template <bool Stable, class ValueType, class Comparator>
class loser_tree_traits_unguarded
{
public:
    typedef LoserTreePointerUnguarded<Stable, ValueType, Comparator> LT;
};

#define STXXL_NO_POINTER_UNGUARDED(T)                             \
    template <bool Stable, class Comparator>                      \
    struct loser_tree_traits_unguarded<Stable, T, Comparator>     \
    {                                                             \
        typedef LoserTreeCopyUnguarded<Stable, T, Comparator> LT; \
    };

STXXL_NO_POINTER_UNGUARDED(unsigned char)
STXXL_NO_POINTER_UNGUARDED(char)
STXXL_NO_POINTER_UNGUARDED(unsigned short)
STXXL_NO_POINTER_UNGUARDED(short)
STXXL_NO_POINTER_UNGUARDED(unsigned int)
STXXL_NO_POINTER_UNGUARDED(int)
STXXL_NO_POINTER_UNGUARDED(unsigned long)
STXXL_NO_POINTER_UNGUARDED(long)
STXXL_NO_POINTER_UNGUARDED(unsigned long long)
STXXL_NO_POINTER_UNGUARDED(long long)

#undef STXXL_NO_POINTER_UNGUARDED

template <bool Stable,
          typename RandomAccessIteratorIterator,
          typename RandomAccessIterator3,
          typename DiffType, typename Comparator>
RandomAccessIterator3
multiway_merge_loser_tree_combined(
    RandomAccessIteratorIterator seqs_begin,
    RandomAccessIteratorIterator seqs_end,
    RandomAccessIterator3 target, DiffType length,
    Comparator comp)
{
    STXXL_PARALLEL_PCALL(length);

    typedef typename std::iterator_traits<RandomAccessIteratorIterator>
        ::value_type::first_type RandomAccessIterator;
    typedef typename std::iterator_traits<RandomAccessIterator>
        ::value_type value_type;

    int min_seq;
    RandomAccessIterator3 target_end;
    DiffType overhang = prepare_unguarded<Stable>(seqs_begin, seqs_end, comp, min_seq);

    DiffType total_length = 0;
    for (RandomAccessIteratorIterator s = seqs_begin; s != seqs_end; ++s)
        total_length += iterpair_size(*s);

    if (overhang != (DiffType)(-1))
    {
        DiffType unguarded_length = std::min(length, total_length - overhang);
        target_end = multiway_merge_loser_tree_unguarded
                     <typename loser_tree_traits_unguarded<Stable, value_type, Comparator>::LT>
                         (seqs_begin, seqs_end, target, unguarded_length, comp);
        overhang = length - unguarded_length;
    }
    else
    {
        // empty sequence found
        overhang = length;
        target_end = target;
    }

    STXXL_DEBUG_ASSERT(target_end == target + length - overhang);
    STXXL_DEBUG_ASSERT(stxxl::is_sorted(target, target_end, comp));

    target_end = multiway_merge_loser_tree
                 <typename loser_tree_traits<Stable, value_type, Comparator>::LT>
                     (seqs_begin, seqs_end, target_end, overhang, comp);

    STXXL_DEBUG_ASSERT(target_end == target + length);
    STXXL_DEBUG_ASSERT(stxxl::is_sorted(target, target_end, comp));

    return target_end;
}

template <bool Stable,
          typename RandomAccessIteratorIterator,
          typename RandomAccessIterator3,
          typename DiffType, typename Comparator>
RandomAccessIterator3
multiway_merge_loser_tree_sentinel(
    RandomAccessIteratorIterator seqs_begin,
    RandomAccessIteratorIterator seqs_end,
    RandomAccessIterator3 target, DiffType length,
    Comparator comp)
{
    STXXL_PARALLEL_PCALL(length);

    typedef typename std::iterator_traits<RandomAccessIteratorIterator>
        ::value_type::first_type RandomAccessIterator;
    typedef typename std::iterator_traits<RandomAccessIterator>
        ::value_type value_type;

    // move end of sequences to include the sentinel for merging
    for (RandomAccessIteratorIterator s = seqs_begin; s != seqs_end; ++s)
        ++(*s).second;

    RandomAccessIterator3 target_end
        = multiway_merge_loser_tree_unguarded
          <typename loser_tree_traits_unguarded<Stable, value_type, Comparator>::LT>
              (seqs_begin, seqs_end, target, length, comp);

    STXXL_DEBUG_ASSERT(target_end == target + length);
    STXXL_DEBUG_ASSERT(stxxl::is_sorted(target, target_end, comp));

    // restore end of sequences
    for (RandomAccessIteratorIterator s = seqs_begin; s != seqs_end; ++s)
        --(*s).second;

    return target_end;
}

/*!
 * Sequential multi-way merging switch.
 *
 * The decision if based on the branching factor and runtime settings.
 *
 * \param seqs_begin Begin iterator of iterator pair input sequence.
 * \param seqs_end End iterator of iterator pair input sequence.
 * \param target Begin iterator out output sequence.
 * \param length Maximum length to merge.
 * \param comp Comparator.
 * \tparam Stable Stable merging incurs a performance penalty.
 * \tparam Sentinels The sequences have a sentinel element.
 * \return End iterator of output sequence.
 */
template <bool Stable, bool Sentinels,
          typename RandomAccessIteratorIterator,
          typename RandomAccessIterator3,
          typename DiffType, typename Comparator>
RandomAccessIterator3
sequential_multiway_merge(RandomAccessIteratorIterator seqs_begin,
                          RandomAccessIteratorIterator seqs_end,
                          RandomAccessIterator3 target, DiffType length,
                          Comparator comp)
{
    STXXL_PARALLEL_PCALL(length);

    typedef typename std::iterator_traits<RandomAccessIteratorIterator>
        ::value_type::first_type RandomAccessIterator;
    typedef typename std::iterator_traits<RandomAccessIterator>
        ::value_type value_type;

    for (RandomAccessIteratorIterator s = seqs_begin; s != seqs_end; ++s)
        STXXL_DEBUG_ASSERT(stxxl::is_sorted((*s).first, (*s).second, comp));

    RandomAccessIterator3 return_target = target;
    int k = static_cast<int>(seqs_end - seqs_begin);

    SETTINGS::MultiwayMergeAlgorithm mwma = SETTINGS::multiway_merge_algorithm;

    if (!Sentinels && mwma == SETTINGS::LOSER_TREE_SENTINEL)
        mwma = SETTINGS::LOSER_TREE_COMBINED;

    switch (k)
    {
    case 0:
        break;
    case 1:
        return_target = std::copy(seqs_begin[0].first,
                                  seqs_begin[0].first + length,
                                  target);
        seqs_begin[0].first += length;
        break;
    case 2:
        return_target = merge_advance(
            seqs_begin[0].first, seqs_begin[0].second,
            seqs_begin[1].first, seqs_begin[1].second,
            target, length, comp);
        break;
    case 3:
        switch (mwma)
        {
        case SETTINGS::LOSER_TREE_COMBINED:
            return_target = multiway_merge_3_combined(
                seqs_begin, seqs_end, target, length, comp);
            break;
        case SETTINGS::LOSER_TREE_SENTINEL:
            return_target = multiway_merge_3_variant<unguarded_iterator>(
                seqs_begin, seqs_end, target, length, comp);
            break;
        default:
            return_target = multiway_merge_3_variant<guarded_iterator>(
                seqs_begin, seqs_end, target, length, comp);
            break;
        }
        break;
    case 4:
        switch (mwma)
        {
        case SETTINGS::LOSER_TREE_COMBINED:
            return_target = multiway_merge_4_combined(
                seqs_begin, seqs_end, target, length, comp);
            break;
        case SETTINGS::LOSER_TREE_SENTINEL:
            return_target = multiway_merge_4_variant<unguarded_iterator>(
                seqs_begin, seqs_end, target, length, comp);
            break;
        default:
            return_target = multiway_merge_4_variant<guarded_iterator>(
                seqs_begin, seqs_end, target, length, comp);
            break;
        }
        break;
    default:
    {
        switch (mwma)
        {
        case SETTINGS::BUBBLE:
            return_target = multiway_merge_bubble<Stable>(
                seqs_begin, seqs_end, target, length, comp);
            break;
        case SETTINGS::LOSER_TREE:
            return_target = multiway_merge_loser_tree<
                typename loser_tree_traits<Stable, value_type, Comparator>::LT>(
                seqs_begin, seqs_end, target, length, comp);
            break;
        case SETTINGS::LOSER_TREE_COMBINED:
            return_target = multiway_merge_loser_tree_combined<Stable>(
                seqs_begin, seqs_end, target, length, comp);
            break;
        case SETTINGS::LOSER_TREE_SENTINEL:
            return_target = multiway_merge_loser_tree_sentinel<Stable>(
                seqs_begin, seqs_end, target, length, comp);
            break;
        default:
            assert(0 && "multiway_merge algorithm not implemented");
            break;
        }
    }
    }

    STXXL_DEBUG_ASSERT(stxxl::is_sorted(target, target + length, comp));

    return return_target;
}

/*!
 * Splitting method for parallel multi-way merge routine: use sampling and
 * binary search for in-exact splitting.
 *
 * \param seqs_begin Begin iterator of iterator pair input sequence.
 * \param seqs_end End iterator of iterator pair input sequence.
 * \param length Maximum length to merge.
 * \param total_length Total length of all sequences combined.
 * \param comp Comparator.
 * \param chunks Output subsequences for num_threads.
 * \param num_threads Split the sequences into for num_threads.
 * \tparam Stable Stable merging incurs a performance penalty.
 * \return End iterator of output sequence.
 */
template <bool Stable,
          typename RandomAccessIteratorIterator,
          typename DiffType,
          typename Comparator>
void
parallel_multiway_merge_sampling_splitting(
    const RandomAccessIteratorIterator& seqs_begin,
    const RandomAccessIteratorIterator& seqs_end,
    DiffType length, DiffType total_length, Comparator comp,
    std::vector<typename std::iterator_traits<RandomAccessIteratorIterator>::value_type>* chunks,
    const thread_index_t num_threads)
{
    typedef typename std::iterator_traits<RandomAccessIteratorIterator>
        ::value_type::first_type RandomAccessIterator;
    typedef typename std::iterator_traits<RandomAccessIterator>
        ::value_type value_type;

    const DiffType num_seqs = seqs_end - seqs_begin;
    const DiffType num_samples = num_threads * SETTINGS::merge_oversampling;

    // pick samples
    value_type* samples = new value_type[num_seqs * num_samples];

    for (DiffType s = 0; s < num_seqs; ++s)
    {
        for (DiffType i = 0; i < num_samples; ++i)
        {
            DiffType sample_index = static_cast<DiffType>(
                double(iterpair_size(seqs_begin[s]))
                * (double(i + 1) / double(num_samples + 1))
                * (double(length) / double(total_length))
                );
            samples[s * num_samples + i] = seqs_begin[s].first[sample_index];
        }
    }

    if (Stable)
        std::stable_sort(samples, samples + (num_samples * num_seqs), comp);
    else
        std::sort(samples, samples + (num_samples * num_seqs), comp);

    // for each processor
    for (thread_index_t slab = 0; slab < num_threads; ++slab)
    {
        // for each sequence
        for (DiffType seq = 0; seq < num_seqs; ++seq)
        {
            if (slab > 0) {
                chunks[slab][seq].first =
                    std::upper_bound(
                        seqs_begin[seq].first, seqs_begin[seq].second,
                        samples[num_samples * num_seqs * slab / num_threads],
                        comp);
            }
            else    // absolute beginning
                chunks[slab][seq].first = seqs_begin[seq].first;

            if ((slab + 1) < num_threads) {
                chunks[slab][seq].second =
                    std::upper_bound(
                        seqs_begin[seq].first, seqs_begin[seq].second,
                        samples[num_samples * num_seqs * (slab + 1) / num_threads],
                        comp);
            }
            else    // absolute ending
                chunks[slab][seq].second = seqs_begin[seq].second;
        }
    }

    delete[] samples;
}

/*!
 * Splitting method for parallel multi-way merge routine: use multisequence
 * selection for exact splitting.
 *
 * \param seqs_begin Begin iterator of iterator pair input sequence.
 * \param seqs_end End iterator of iterator pair input sequence.
 * \param length Maximum length to merge.
 * \param total_length Total length of all sequences combined.
 * \param comp Comparator.
 * \param chunks Output subsequences for num_threads.
 * \param num_threads Split the sequences into for num_threads.
 * \tparam Stable Stable merging incurs a performance penalty.
 * \return End iterator of output sequence.
 */
template <bool Stable,
          typename RandomAccessIteratorIterator,
          typename DiffType,
          typename Comparator>
void
parallel_multiway_merge_exact_splitting(
    const RandomAccessIteratorIterator& seqs_begin,
    const RandomAccessIteratorIterator& seqs_end,
    DiffType length, DiffType total_length, Comparator comp,
    std::vector<typename std::iterator_traits<RandomAccessIteratorIterator>::value_type>* chunks,
    const thread_index_t num_threads)
{
    typedef typename std::iterator_traits<RandomAccessIteratorIterator>
        ::value_type RandomAccessIteratorPair;
    typedef typename RandomAccessIteratorPair
        ::first_type RandomAccessIterator;

    const size_t num_seqs = seqs_end - seqs_begin;
    const bool tight = (total_length == length);

    std::vector<RandomAccessIterator>* offsets
        = new std::vector<RandomAccessIterator>[num_threads];

    std::vector<DiffType> ranks(num_threads + 1);
    equally_split(length, num_threads, ranks.begin());

    for (thread_index_t s = 0; s < (num_threads - 1); ++s)
    {
        offsets[s].resize(num_seqs);
        multiseq_partition(seqs_begin, seqs_end,
                           ranks[s + 1], offsets[s].begin(), comp);

        if (!tight) // last one also needed and available
        {
            offsets[num_threads - 1].resize(num_seqs);
            multiseq_partition(seqs_begin, seqs_end,
                               length, offsets[num_threads - 1].begin(), comp);
        }
    }

    // for each processor
    for (thread_index_t slab = 0; slab < num_threads; ++slab)
    {
        // for each sequence
        for (size_t s = 0; s < num_seqs; ++s)
        {
            if (slab == 0) // absolute beginning
                chunks[slab][s].first = seqs_begin[s].first;
            else
                chunks[slab][s].first = offsets[slab - 1][s];

            if (!tight || slab < (num_threads - 1))
                chunks[slab][s].second = offsets[slab][s];
            else    // slab == num_threads - 1
                chunks[slab][s].second = seqs_begin[s].second;
        }
    }

    delete[] offsets;
}

#if STXXL_PARALLEL

/*!
 * Parallel multi-way merge routine.
 *
 * The decision if based on the branching factor and runtime settings.
 *
 * \param seqs_begin Begin iterator of iterator pair input sequence.
 * \param seqs_end End iterator of iterator pair input sequence.
 * \param target Begin iterator out output sequence.
 * \param length Maximum length to merge.
 * \param comp Comparator.
 * \tparam Stable Stable merging incurs a performance penalty.
 * \return End iterator of output sequence.
 */
template <bool Stable,
          typename RandomAccessIteratorIterator,
          typename RandomAccessIterator3,
          typename DiffType,
          typename Comparator>
RandomAccessIterator3
parallel_multiway_merge(RandomAccessIteratorIterator seqs_begin,
                        RandomAccessIteratorIterator seqs_end,
                        RandomAccessIterator3 target, const DiffType length,
                        Comparator comp)
{
    STXXL_PARALLEL_PCALL(length);

    typedef typename std::iterator_traits<RandomAccessIteratorIterator>
        ::value_type RandomAccessIteratorPair;

    for (RandomAccessIteratorIterator rii = seqs_begin; rii != seqs_end; ++rii)
        STXXL_DEBUG_ASSERT(stxxl::is_sorted((*rii).first, (*rii).second, comp));

    // leave only non-empty sequences
    std::vector<RandomAccessIteratorPair> seqs_ne;
    seqs_ne.reserve(seqs_end - seqs_begin);
    DiffType total_length = 0;

    for (RandomAccessIteratorIterator raii = seqs_begin; raii != seqs_end; ++raii)
    {
        DiffType length = iterpair_size(*raii);
        if (length > 0) {
            total_length += length;
            seqs_ne.push_back(*raii);
        }
    }

    size_t num_seqs = seqs_ne.size();

    STXXL_PARALLEL_PCALL(total_length);

    if (total_length == 0 || num_seqs == 0)
        return target;

    thread_index_t num_threads = static_cast<thread_index_t>(
        std::min(static_cast<DiffType>(SETTINGS::num_threads), total_length));

    Timing<inactive_tag>* t = new Timing<inactive_tag>[num_threads];

    for (int pr = 0; pr < num_threads; ++pr)
        t[pr].tic();

    // thread t will have to merge chunks[iam][0..k - 1]

    std::vector<RandomAccessIteratorPair>* chunks
        = new std::vector<RandomAccessIteratorPair>[num_threads];

    for (int s = 0; s < num_threads; ++s)
        chunks[s].resize(num_seqs);

#pragma omp parallel num_threads(num_threads)
    {
#pragma omp single
        {
            if (SETTINGS::multiway_merge_splitting == SETTINGS::SAMPLING)
            {
                parallel_multiway_merge_sampling_splitting<Stable>(
                    seqs_ne.begin(), seqs_ne.end(),
                    length, total_length, comp,
                    chunks, num_threads);
            }
            else // (SETTINGS::multiway_merge_splitting == SETTINGS::EXACT)
            {
                parallel_multiway_merge_exact_splitting<Stable>(
                    seqs_ne.begin(), seqs_ne.end(),
                    length, total_length, comp,
                    chunks, num_threads);
            }
        }

        thread_index_t iam = omp_get_thread_num();
        t[iam].tic();

        DiffType target_position = 0, local_length = 0;

        for (size_t s = 0; s < num_seqs; ++s)
        {
            target_position += chunks[iam][s].first - seqs_ne[s].first;
            local_length += iterpair_size(chunks[iam][s]);
        }

        sequential_multiway_merge<Stable, false>(
            chunks[iam].begin(), chunks[iam].end(),
            target + target_position,
            std::min(local_length, length - target_position),
            comp);

        t[iam].tic();
    }

    for (int pr = 0; pr < num_threads; ++pr)
        t[pr].tic();

    STXXL_DEBUG_ASSERT(stxxl::is_sorted(target, target + length, comp));

    //update ends of sequences
    size_t count_seqs = 0;
    for (RandomAccessIteratorIterator raii = seqs_begin; raii != seqs_end; ++raii)
    {
        DiffType length = iterpair_size(*raii);
        if (length > 0)
            raii->first = chunks[num_threads - 1][count_seqs++].second;
    }
    STXXL_DEBUG_ASSERT(count_seqs == num_seqs);

    delete[] chunks;

    for (int pr = 0; pr < num_threads; ++pr)
        t[pr].tic();
    for (int pr = 0; pr < num_threads; ++pr)
        t[pr].print();
    delete[] t;

    return target + length;
}

/*!
 * Multi-way merging front-end with unstable mode and without sentinels.
 *
 * \param seqs_begin Begin iterator of iterator pair input sequence.
 * \param seqs_end End iterator of iterator pair input sequence.
 * \param target Begin iterator out output sequence.
 * \param comp Comparator.
 * \param length Maximum length to merge.
 * \return End iterator of output sequence.
 */
template <typename RandomAccessIteratorPairIterator,
          typename RandomAccessIterator3,
          typename DiffType, typename Comparator>
RandomAccessIterator3
multiway_merge(RandomAccessIteratorPairIterator seqs_begin,
               RandomAccessIteratorPairIterator seqs_end,
               RandomAccessIterator3 target, DiffType length,
               Comparator comp)
{
    STXXL_PARALLEL_PCALL(seqs_end - seqs_begin);

    if (seqs_begin == seqs_end)
        return target;

    RandomAccessIterator3 target_end;
    if (STXXL_PARALLEL_CONDITION(
            ((seqs_end - seqs_begin) >= SETTINGS::multiway_merge_minimal_k) &&
            ((sequence_index_t)length >= SETTINGS::multiway_merge_minimal_n)
            ))
        target_end = parallel_multiway_merge<false>(
            seqs_begin, seqs_end, target, length, comp);
    else
        target_end = sequential_multiway_merge<false, false>(
            seqs_begin, seqs_end, target, length, comp);

    return target_end;
}

/*!
 * Multi-way merging front-end with unstable mode and without sentinels.
 *
 * \param seqs_begin Begin iterator of iterator pair input sequence.
 * \param seqs_end End iterator of iterator pair input sequence.
 * \param target Begin iterator out output sequence.
 * \param comp Comparator.
 * \param length Maximum length to merge.
 * \return End iterator of output sequence.
 */
template <typename RandomAccessIteratorPairIterator,
          typename RandomAccessIterator3,
          typename DiffType, typename Comparator>
RandomAccessIterator3
multiway_merge_stable(RandomAccessIteratorPairIterator seqs_begin,
                      RandomAccessIteratorPairIterator seqs_end,
                      RandomAccessIterator3 target, DiffType length,
                      Comparator comp)
{
    STXXL_PARALLEL_PCALL(seqs_end - seqs_begin);

    if (seqs_begin == seqs_end)
        return target;

    RandomAccessIterator3 target_end;
    if (STXXL_PARALLEL_CONDITION(
            ((seqs_end - seqs_begin) >= SETTINGS::multiway_merge_minimal_k) &&
            ((sequence_index_t)length >= SETTINGS::multiway_merge_minimal_n)
            ))
        target_end = parallel_multiway_merge<true>(
            seqs_begin, seqs_end, target, length, comp);
    else
        target_end = sequential_multiway_merge<true, false>(
            seqs_begin, seqs_end, target, length, comp);

    return target_end;
}

/*!
 * Multi-way merging front-end with unstable mode and sentinels.
 *
 * Each sequence must be suffixed with a sentinel as *end(), one item beyond
 * the end of each sequence.
 *
 * \param seqs_begin Begin iterator of iterator pair input sequence.
 * \param seqs_end End iterator of iterator pair input sequence.
 * \param target Begin iterator out output sequence.
 * \param comp Comparator.
 * \param length Maximum length to merge.
 * \return End iterator of output sequence.
 * \pre For each \c i, \c seqs_begin[i].second must be the end marker of the
 * sequence, but also reference the one more sentinel element.
 */
template <typename RandomAccessIteratorPairIterator,
          typename RandomAccessIterator3,
          typename DiffType, typename Comparator>
RandomAccessIterator3
multiway_merge_sentinels(RandomAccessIteratorPairIterator seqs_begin,
                         RandomAccessIteratorPairIterator seqs_end,
                         RandomAccessIterator3 target, DiffType length,
                         Comparator comp)
{
    if (seqs_begin == seqs_end)
        return target;

    STXXL_PARALLEL_PCALL(seqs_end - seqs_begin);

    if (STXXL_PARALLEL_CONDITION(
            ((seqs_end - seqs_begin) >= SETTINGS::multiway_merge_minimal_k) &&
            ((sequence_index_t)length >= SETTINGS::multiway_merge_minimal_n)
            ))
        return parallel_multiway_merge<false>(
            seqs_begin, seqs_end, target, length, comp);
    else
        return sequential_multiway_merge<false, true>(
            seqs_begin, seqs_end, target, length, comp);
}

/*!
 * Multi-way merging front-end with unstable mode and sentinels.
 *
 * Each sequence must be suffixed with a sentinel as *end(), one item beyond
 * the end of each sequence.
 *
 * \param seqs_begin Begin iterator of iterator pair input sequence.
 * \param seqs_end End iterator of iterator pair input sequence.
 * \param target Begin iterator out output sequence.
 * \param comp Comparator.
 * \param length Maximum length to merge.
 * \return End iterator of output sequence.
 * \pre For each \c i, \c seqs_begin[i].second must be the end marker of the
 * sequence, but also reference the one more sentinel element.
 */
template <typename RandomAccessIteratorPairIterator,
          typename RandomAccessIterator3,
          typename DiffType, typename Comparator>
RandomAccessIterator3
multiway_merge_stable_sentinels(RandomAccessIteratorPairIterator seqs_begin,
                                RandomAccessIteratorPairIterator seqs_end,
                                RandomAccessIterator3 target, DiffType length,
                                Comparator comp)
{
    if (seqs_begin == seqs_end)
        return target;

    STXXL_PARALLEL_PCALL(seqs_end - seqs_begin);

    if (STXXL_PARALLEL_CONDITION(
            ((seqs_end - seqs_begin) >= SETTINGS::multiway_merge_minimal_k) &&
            ((sequence_index_t)length >= SETTINGS::multiway_merge_minimal_n)
            ))
        return parallel_multiway_merge<true>(
            seqs_begin, seqs_end, target, length, comp);
    else
        return sequential_multiway_merge<true, true>(
            seqs_begin, seqs_end, target, length, comp);
}

#endif // STXXL_PARALLEL

} // namespace parallel

STXXL_END_NAMESPACE

#endif // !STXXL_PARALLEL_MULTIWAY_MERGE_HEADER
