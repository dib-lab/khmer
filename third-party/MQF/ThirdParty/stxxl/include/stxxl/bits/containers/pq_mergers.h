/***************************************************************************
 *  include/stxxl/bits/containers/pq_mergers.h
 *
 *  Part of the STXXL. See http://stxxl.sourceforge.net
 *
 *  Copyright (C) 1999 Peter Sanders <sanders@mpi-sb.mpg.de>
 *  Copyright (C) 2003, 2004, 2007 Roman Dementiev <dementiev@mpi-sb.mpg.de>
 *  Copyright (C) 2007-2009 Johannes Singler <singler@ira.uka.de>
 *  Copyright (C) 2007, 2008 Andreas Beckmann <beckmann@cs.uni-frankfurt.de>
 *
 *  Distributed under the Boost Software License, Version 1.0.
 *  (See accompanying file LICENSE_1_0.txt or copy at
 *  http://www.boost.org/LICENSE_1_0.txt)
 **************************************************************************/

#ifndef STXXL_CONTAINERS_PQ_MERGERS_HEADER
#define STXXL_CONTAINERS_PQ_MERGERS_HEADER

#include <stxxl/bits/containers/pq_helpers.h>

STXXL_BEGIN_NAMESPACE

//! \addtogroup stlcontinternals
//!
//! \{

/*! \internal
 */
namespace priority_queue_local {

////////////////////////////////////////////////////////////////////////////////
// auxiliary functions

// merge length elements from the two sentinel terminated input
// sequences source0 and source1 to target
// advance source0 and source1 accordingly
// require: at least length nonsentinel elements available in source0, source1
// require: target may overwrite one of the sources as long as
//   *(sourcex + length) is before the end of sourcex
template <class InputIterator, class OutputIterator, class CompareType>
void merge2_iterator(
    InputIterator& source0,
    InputIterator& source1,
    OutputIterator target, OutputIterator end,
    CompareType& cmp)
{
    while (target != end)
    {
        if (cmp(*source0, *source1))
        {
            *target = *source1;
            ++source1;
        }
        else
        {
            *target = *source0;
            ++source0;
        }
        ++target;
    }
}

// merge length elements from the three sentinel terminated input
// sequences source0, source1 and source2 to target
// advance source0, source1 and source2 accordingly
// require: at least length nonsentinel elements available in source0, source1 and source2
// require: target may overwrite one of the sources as long as
//   *(sourcex + length) is before the end of sourcex
template <class InputIterator, class OutputIterator,
          class CompareType>
void merge3_iterator(
    InputIterator& source0,
    InputIterator& source1,
    InputIterator& source2,
    OutputIterator target, OutputIterator end,
    CompareType& cmp)
{
    if (cmp(*source1, *source0)) {
        if (cmp(*source2, *source1)) {
            goto s012;
        }
        else {
            if (cmp(*source0, *source2)) {
                goto s201;
            }
            else {
                goto s021;
            }
        }
    }
    else {
        if (cmp(*source2, *source1)) {
            if (cmp(*source2, *source0)) {
                goto s102;
            }
            else {
                goto s120;
            }
        }
        else {
            goto s210;
        }
    }

#define Merge3Case(a, b, c)              \
    s ## a ## b ## c:                    \
    if (target == end)                   \
        return;                          \
    *target = *source ## a;              \
    ++target;                            \
    ++source ## a;                       \
    if (cmp(*source ## b, *source ## a)) \
        goto s ## a ## b ## c;           \
    if (cmp(*source ## c, *source ## a)) \
        goto s ## b ## a ## c;           \
    goto s ## b ## c ## a;

    // the order is chosen in such a way that
    // four of the trailing gotos can be eliminated by the optimizer
    Merge3Case(0, 1, 2);
    Merge3Case(1, 2, 0);
    Merge3Case(2, 0, 1);
    Merge3Case(1, 0, 2);
    Merge3Case(0, 2, 1);
    Merge3Case(2, 1, 0);

#undef Merge3Case
}

// merge length elements from the four sentinel terminated input
// sequences source0, source1, source2 and source3 to target
// advance source0, source1, source2 and source3 accordingly
// require: at least length nonsentinel elements available in source0, source1, source2 and source3
// require: target may overwrite one of the sources as long as
//   *(sourcex + length) is before the end of sourcex
template <class InputIterator, class OutputIterator,
          class CompareType>
void merge4_iterator(
    InputIterator& source0,
    InputIterator& source1,
    InputIterator& source2,
    InputIterator& source3,
    OutputIterator target, OutputIterator end,
    CompareType& cmp)
{
#define StartMerge4(a, b, c, d)               \
    if ((!cmp(*source ## a, *source ## b)) && \
        (!cmp(*source ## b, *source ## c)) && \
        (!cmp(*source ## c, *source ## d)))   \
        goto s ## a ## b ## c ## d;

    // b>a c>b d>c
    // a<b b<c c<d
    // a<=b b<=c c<=d
    // !(a>b) !(b>c) !(c>d)

    StartMerge4(0, 1, 2, 3);
    StartMerge4(1, 2, 3, 0);
    StartMerge4(2, 3, 0, 1);
    StartMerge4(3, 0, 1, 2);

    StartMerge4(0, 3, 1, 2);
    StartMerge4(3, 1, 2, 0);
    StartMerge4(1, 2, 0, 3);
    StartMerge4(2, 0, 3, 1);

    StartMerge4(0, 2, 3, 1);
    StartMerge4(2, 3, 1, 0);
    StartMerge4(3, 1, 0, 2);
    StartMerge4(1, 0, 2, 3);

    StartMerge4(2, 0, 1, 3);
    StartMerge4(0, 1, 3, 2);
    StartMerge4(1, 3, 2, 0);
    StartMerge4(3, 2, 0, 1);

    StartMerge4(3, 0, 2, 1);
    StartMerge4(0, 2, 1, 3);
    StartMerge4(2, 1, 3, 0);
    StartMerge4(1, 3, 0, 2);

    StartMerge4(1, 0, 3, 2);
    StartMerge4(0, 3, 2, 1);
    StartMerge4(3, 2, 1, 0);
    StartMerge4(2, 1, 0, 3);

#define Merge4Case(a, b, c, d)               \
    s ## a ## b ## c ## d:                   \
    if (target == end)                       \
        return;                              \
    *target = *source ## a;                  \
    ++target;                                \
    ++source ## a;                           \
    if (cmp(*source ## c, *source ## a))     \
    {                                        \
        if (cmp(*source ## b, *source ## a)) \
            goto s ## a ## b ## c ## d;      \
        else                                 \
            goto s ## b ## a ## c ## d;      \
    }                                        \
    else                                     \
    {                                        \
        if (cmp(*source ## d, *source ## a)) \
            goto s ## b ## c ## a ## d;      \
        else                                 \
            goto s ## b ## c ## d ## a;      \
    }

    Merge4Case(0, 1, 2, 3);
    Merge4Case(1, 2, 3, 0);
    Merge4Case(2, 3, 0, 1);
    Merge4Case(3, 0, 1, 2);

    Merge4Case(0, 3, 1, 2);
    Merge4Case(3, 1, 2, 0);
    Merge4Case(1, 2, 0, 3);
    Merge4Case(2, 0, 3, 1);

    Merge4Case(0, 2, 3, 1);
    Merge4Case(2, 3, 1, 0);
    Merge4Case(3, 1, 0, 2);
    Merge4Case(1, 0, 2, 3);

    Merge4Case(2, 0, 1, 3);
    Merge4Case(0, 1, 3, 2);
    Merge4Case(1, 3, 2, 0);
    Merge4Case(3, 2, 0, 1);

    Merge4Case(3, 0, 2, 1);
    Merge4Case(0, 2, 1, 3);
    Merge4Case(2, 1, 3, 0);
    Merge4Case(1, 3, 0, 2);

    Merge4Case(1, 0, 3, 2);
    Merge4Case(0, 3, 2, 1);
    Merge4Case(3, 2, 1, 0);
    Merge4Case(2, 1, 0, 3);

#undef StartMerge4
#undef Merge4Case
}

////////////////////////////////////////////////////////////////////////////////
// Loser tree data structure from Knuth, "Sorting and Searching", Section 5.4.1

/*!
 * Loser tree from Knuth, "Sorting and Searching", Section 5.4.1
 * \tparam Arity  maximum arity of merger, does not need to be a power of 2
 */
template <class ArraysType, class CompareType, unsigned Arity>
class loser_tree : private noncopyable
{
public:
    //! type of arrays container linked with loser tree
    typedef ArraysType arrays_type;
    //! comparator object type
    typedef CompareType compare_type;

    // arity_bound / 2  <  arity  <=  arity_bound
    enum { arity = Arity, max_arity = 1UL << (LOG2<Arity>::ceil) };

    //! type of values stored in the arrays container
    typedef typename arrays_type::value_type value_type;
    //! type of the ordered sequences in the arrays container
    typedef typename arrays_type::sequence_type sequence_type;

public:
    //! the comparator object
    compare_type cmp;

    //! current tree size, invariant (k == 1 << logK), always a power of two
    unsigned_type k;
    //! log of current tree size
    unsigned_type logK;

    // only entries 0 .. arity-1 may hold actual sequences, the other
    // entries arity .. max_arity-1 are sentinels to make the size of the tree
    // a power of 2 always

protected:
    //! reference to the linked arrays
    arrays_type& arrays;

    //! type of nodes in the loser tree
    struct Entry
    {
        value_type key;          //!< Key of Loser element (winner for 0)
        unsigned_type index;     //!< number of losing segment
    };

    //! levels of loser tree: entry[0] contains the winner info
    struct Entry entry[max_arity];

    //! stack of free player indices
    internal_bounded_stack<unsigned_type, arity> free_slots;

public:
    loser_tree(const compare_type& c, arrays_type& a)
        : cmp(c), k(1), logK(0), arrays(a)
    {
        // verify strict weak ordering
        assert(!cmp(cmp.min_value(), cmp.min_value()));
    }

    void initialize()
    {
        // initial state: one empty player slot
        free_slots.push(0);

        rebuild_loser_tree();

        assert(arrays.is_array_empty(0) && !arrays.is_array_allocated(0));
    }

    //! True if a is the sentinel value
    bool is_sentinel(const value_type& a) const
    {
        return !(cmp(cmp.min_value(), a)); // a <= cmp.min_value()
    }

    //! Allocate a free slot for a new player.
    unsigned_type new_player()
    {
        // get a free slot
        if (free_slots.empty()) {
            // tree is too small, attempt to enlarge
            double_k();
        }

        assert(!free_slots.empty());
        unsigned_type index = free_slots.top();
        free_slots.pop();

        return index;
    }

    //! Free a finished player's slot
    void free_player(unsigned_type slot)
    {
        // push on the stack of free segment indices
        free_slots.push(slot);
    }

    //! Whether there is still space for new array
    bool is_space_available() const
    {
        return (k < arity) || !free_slots.empty();
    }

    //! rebuild loser tree information from the values in current
    void rebuild_loser_tree()
    {
        unsigned_type winner = init_winner(1);
        entry[0].index = winner;
        entry[0].key = *arrays.get_array(winner);
    }

    // given any values in the leaves this
    // routing recomputes upper levels of the tree
    // from scratch in linear time
    // initialize entry[root].index and the subtree rooted there
    // return winner index
    unsigned_type init_winner(unsigned_type root)
    {
        if (root >= k || root >= max_arity)
        {       // leaf reached
            return root - k;
        }
        else
        {
            unsigned_type left = init_winner(2 * root);
            unsigned_type right = init_winner(2 * root + 1);
            value_type lk = *arrays.get_array(left);
            value_type rk = *arrays.get_array(right);
            assert(root < max_arity);

            if (!(cmp(lk, rk)))
            {
                // right subtree looses
                entry[root].index = right;
                entry[root].key = rk;
                return left;
            }
            else
            {
                entry[root].index = left;
                entry[root].key = lk;
                return right;
            }
        }
    }

    /*!
     * Update loser tree on insert or decrement of player index first go up the
     * tree all the way to the root hand down old winner for the respective
     * subtree based on new value, and old winner and loser update each node on
     * the path to the root top down.  This is implemented recursively
     */
    void update_on_insert(unsigned_type node,
                          const value_type& newKey, unsigned_type newIndex,
                          value_type* winner_key,
                          unsigned_type* winner_index,   // old winner
                          unsigned_type* mask)           // 1 << (ceil(log KNK) - dist-from-root)
    {
        if (node == 0)
        {
            // winner part of root
            *mask = (unsigned_type)(1) << (logK - 1);
            *winner_key = entry[0].key;
            *winner_index = entry[0].index;
            if (cmp(entry[node].key, newKey))
            {
                entry[node].key = newKey;
                entry[node].index = newIndex;
            }
        }
        else
        {
            update_on_insert(node >> 1, newKey, newIndex, winner_key, winner_index, mask);
            value_type loserKey = entry[node].key;
            unsigned_type loserIndex = entry[node].index;
            if ((*winner_index & *mask) != (newIndex & *mask)) {     // different subtrees
                // newKey will have influence here
                if (cmp(loserKey, newKey)) {
                    if (cmp(*winner_key, newKey)) {
                        // old winner loses here
                        entry[node].key = *winner_key;
                        entry[node].index = *winner_index;
                    }
                    else {
                        // new entry loses here
                        entry[node].key = newKey;
                        entry[node].index = newIndex;
                    }
                }
                *winner_key = loserKey;
                *winner_index = loserIndex;
            }
            // note that nothing needs to be done if the winner came from the
            // same subtree
            // a) newKey <= winner_key => even more reason for the other tree to lose
            // b) newKey >  winner_key => the old winner will beat the new
            //                           entry further down the tree
            // also the same old winner is handed down the tree

            *mask >>= 1;     // next level
        }
    }

    //! Initial call to recursive update_on_insert
    void update_on_insert(unsigned_type node,
                          const value_type& newKey, unsigned_type newIndex)
    {
        value_type dummyKey;
        unsigned_type dummyIndex, dummyMask;
        update_on_insert(node, newKey, newIndex,
                         &dummyKey, &dummyIndex, &dummyMask);
    }

    //! make the tree twice as wide
    void double_k()
    {
        STXXL_VERBOSE1("double_k (before) k=" << k << " logK=" << logK << " arity=" << arity << " max_arity=" << max_arity << " #free=" << free_slots.size());
        assert(k > 0);
        assert(k < arity);
        assert(free_slots.empty());                    // stack was free (probably not needed)

        // make all new entries free and push them on the free stack
        for (unsigned_type i = 2 * k - 1; i >= k; i--) //backwards
        {
            arrays.make_array_sentinel(i);
            if (i < arity)
                free_slots.push(i);
        }

        // double the size
        k *= 2;
        logK++;

        STXXL_VERBOSE1("double_k (after)  k=" << k << " logK=" << logK << " arity=" << arity << " max_arity=" << max_arity << " #free=" << free_slots.size());
        assert(!free_slots.empty());
        assert(k <= max_arity);

        // recompute loser tree information
        rebuild_loser_tree();
    }

    //! compact nonempty segments in the left half of the tree
    void compact_tree()
    {
        STXXL_VERBOSE3("compact_tree (before) k=" << k << " logK=" << logK << " #free=" << free_slots.size());
        assert(logK > 0);

        // compact all nonempty segments to the left
        unsigned_type last_empty = 0;
        for (unsigned_type pos = 0; pos < k; pos++)
        {
            if (!arrays.is_array_empty(pos))
            {
                assert(arrays.is_array_allocated(pos));
                if (pos != last_empty)
                {
                    assert(!arrays.is_array_allocated(last_empty));
                    arrays.swap_arrays(last_empty, pos);
                }
                ++last_empty;
            }
            /*
              else
              {
              if(segment[pos])
              {
              STXXL_VERBOSE2("int_arrays::compact_tree() deleting segment "<<pos<<
              " address: "<<segment[pos]<<" size: "<<segment_size[pos]);
              delete [] segment[pos];
              segment[pos] = 0;
              mem_cons_ -= segment_size[pos];
              }
              }*/
        }

        // half degree as often as possible
        while ((k > 1) && last_empty <= (k / 2))
        {
            k /= 2;
            logK--;
        }

        // overwrite garbage and compact the stack of free segment indices
        free_slots.clear(); // none free
        for ( ; last_empty < k; last_empty++)
        {
            assert(!arrays.is_array_allocated(last_empty));
            arrays.make_array_sentinel(last_empty);
            if (last_empty < arity)
                free_slots.push(last_empty);
        }

        STXXL_VERBOSE3("compact_tree (after)  k=" << k << " logK=" << logK << " #free=" << free_slots.size());

        // recompute loser tree information
        rebuild_loser_tree();
    }

    //! compact tree if it got considerably smaller
    void maybe_compact()
    {
        const unsigned_type num_segments_used = k - free_slots.size();
        const unsigned_type num_segments_trigger = k - (3 * k / 5);
        // using k/2 would be worst case inefficient (for large k)
        // for k \in {2, 4, 8} the trigger is k/2 which is good
        // because we have special mergers for k \in {1, 2, 4}
        // there is also a special 3-way-merger, that will be
        // triggered if k == 4 && is_array_atsentinel(3)
        STXXL_VERBOSE3("int_merger  compact? k=" << k << " #used=" << num_segments_used
                                                 << " <= #trigger=" << num_segments_trigger << " ==> "
                                                 << ((k > 1 && num_segments_used <= num_segments_trigger) ? "yes" : "no ")
                                                 << " || "
                                                 << ((k == 4 && !free_slots.empty() && !arrays.is_array_empty(3)) ? "yes" : "no ")
                                                 << " #free=" << free_slots.size());
        if (k > 1 &&
            ((num_segments_used <= num_segments_trigger) ||
             (k == 4 && !free_slots.empty() && !arrays.is_array_empty(3))))
        {
            compact_tree();
        }
    }

    //! multi-merge for arbitrary K
    template <class OutputIterator>
    void multi_merge_k(OutputIterator begin, OutputIterator end)
    {
        Entry* current_pos;
        value_type current_key;
        unsigned_type current_index; // leaf pointed to by current entry
        unsigned_type winner_index = entry[0].index;
        value_type winner_key = entry[0].key;

        while (begin != end)
        {
            // write result
            *begin++ = *arrays.get_array(winner_index);

            // advance winner segment
            ++(arrays.get_array(winner_index));

            winner_key = *arrays.get_array(winner_index);

            // remove winner segment if empty now
            if (is_sentinel(winner_key))
                arrays.free_array(winner_index);

            // go up the entry-tree
            for (unsigned_type i = (winner_index + k) >> 1; i > 0; i >>= 1)
            {
                current_pos = entry + i;
                current_key = current_pos->key;
                if (cmp(winner_key, current_key))
                {
                    current_index = current_pos->index;
                    current_pos->key = winner_key;
                    current_pos->index = winner_index;
                    winner_key = current_key;
                    winner_index = current_index;
                }
            }
        }
        entry[0].index = winner_index;
        entry[0].key = winner_key;
    }

    template <int LogK, typename OutputIterator>
    void multi_merge_f(OutputIterator begin, OutputIterator end)
    {
        unsigned_type winner_index = entry[0].index;
        value_type winner_key = entry[0].key;

        // TODO: reinsert assert(log_k >= LogK);
        while (begin != end)
        {
            // write result
            *begin++ = *arrays.get_array(winner_index);

            // advance winner segment
            ++(arrays.get_array(winner_index));

            winner_key = *arrays.get_array(winner_index);

            // remove winner segment if empty now
            if (is_sentinel(winner_key))
                arrays.free_array(winner_index);

            // update loser tree
#define TreeStep(L)                                                        \
    if (1 << LogK >= 1 << L) {                                             \
        int pos_shift = ((int(LogK - L) + 1) >= 0) ? ((LogK - L) + 1) : 0; \
        Entry* pos = entry + ((winner_index + (1 << LogK)) >> pos_shift);  \
        value_type key = pos->key;                                         \
        if (cmp(winner_key, key)) {                                        \
            unsigned_type index = pos->index;                              \
            pos->key = winner_key;                                         \
            pos->index = winner_index;                                     \
            winner_key = key;                                              \
            winner_index = index;                                          \
        }                                                                  \
    }
            TreeStep(10);
            TreeStep(9);
            TreeStep(8);
            TreeStep(7);
            TreeStep(6);
            TreeStep(5);
            TreeStep(4);
            TreeStep(3);
            TreeStep(2);
            TreeStep(1);
#undef TreeStep
        }
        entry[0].index = winner_index;
        entry[0].key = winner_key;
    }

    //! extract the (length = end - begin) smallest elements and write them to
    //! [begin..end) empty segments are deallocated. Requires: there are at
    //! least length elements and segments are ended by sentinels.
    template <class OutputIterator>
    void multi_merge(OutputIterator begin, OutputIterator end)
    {
        int_type length = end - begin;

        STXXL_VERBOSE3("multi_merge(length=" << length << ") from sequences k=" << k);

        if (begin == end)
            return;

        assert(k > 0);

        // This is the place to make statistics about external multi_merge calls.

        arrays.prefetch_arrays();

        switch (logK) {
        case 0: {
            assert(k == 1);
            assert(entry[0].index == 0);
            assert(free_slots.empty());

            // in int_merger:
            // memcpy(target, states[0], length * sizeof(value_type));

            sequence_type& seq = arrays.get_array(0);
            for (int_type i = 0; i < length; ++i, ++seq, ++begin)
                *begin = *seq;
            entry[0].key = *seq;

            if (arrays.is_array_empty(0) && arrays.is_array_allocated(0))
                arrays.free_array(0);

            break;
        }
        case 1:
            assert(k == 2);
            merge2_iterator(arrays.get_array(0), arrays.get_array(1),
                            begin, end, cmp);
            rebuild_loser_tree();

            if (arrays.is_array_empty(0) && arrays.is_array_allocated(0))
                arrays.free_array(0);

            if (arrays.is_array_empty(1) && arrays.is_array_allocated(1))
                arrays.free_array(1);

            break;
        case 2:
            assert(k == 4);
            if (arrays.is_array_empty(3))
                merge3_iterator(arrays.get_array(0), arrays.get_array(1),
                                arrays.get_array(2),
                                begin, end, cmp);
            else
                merge4_iterator(arrays.get_array(0), arrays.get_array(1),
                                arrays.get_array(2), arrays.get_array(3),
                                begin, end, cmp);

            rebuild_loser_tree();

            if (arrays.is_array_empty(0) && arrays.is_array_allocated(0))
                arrays.free_array(0);

            if (arrays.is_array_empty(1) && arrays.is_array_allocated(1))
                arrays.free_array(1);

            if (arrays.is_array_empty(2) && arrays.is_array_allocated(2))
                arrays.free_array(2);

            if (arrays.is_array_empty(3) && arrays.is_array_allocated(3))
                arrays.free_array(3);

            break;
        case  3: multi_merge_f<3>(begin, end);
            break;
        case  4: multi_merge_f<4>(begin, end);
            break;
        case  5: multi_merge_f<5>(begin, end);
            break;
        case  6: multi_merge_f<6>(begin, end);
            break;
        case  7: multi_merge_f<7>(begin, end);
            break;
        case  8: multi_merge_f<8>(begin, end);
            break;
        case  9: multi_merge_f<9>(begin, end);
            break;
        case 10: multi_merge_f<10>(begin, end);
            break;
        default: multi_merge_k(begin, end);
            break;
        }

        maybe_compact();

        //std::copy(target,target + length,std::ostream_iterator<ValueType>(std::cout, "\n"));
    }

    void swap(loser_tree& obj)
    {
        std::swap(free_slots, obj.free_slots);
        swap_1D_arrays(entry, obj.entry, max_arity);
    }
};

#if STXXL_PARALLEL && (STXXL_PARALLEL_PQ_MULTIWAY_MERGE_EXTERNAL || STXXL_PARALLEL_PQ_MULTIWAY_MERGE_INTERNAL)
/*!
 * Adapter from loser_tree to parallel merger
 * __gnu_parallel::multiway_merge. This class holds most attributes of
 * loser_tree, except for the tree itself: it thus basically only manages an
 * array of sequence.
 * \tparam Arity  maximum arity of merger, does not need to be a power of 2
 */
template <class ArraysType, class CompareType, unsigned Arity>
class parallel_merger_adapter : private noncopyable
{
public:
    //! type of arrays container linked with loser tree
    typedef ArraysType arrays_type;
    //! comparator object type
    typedef CompareType compare_type;

    // arity_bound / 2  <  arity  <=  arity_bound
    enum { arity = Arity, max_arity = 1UL << (LOG2<Arity>::ceil) };

    //! type of values stored in the arrays container
    typedef typename arrays_type::value_type value_type;
    //! type of the ordered sequences in the arrays container
    typedef typename arrays_type::sequence_type sequence_type;

public:
    //! the comparator object
    compare_type cmp;

    //! current tree size, invariant (k == 1 << logK), always a power of two
    unsigned_type k;
    //! log of current tree size
    unsigned_type logK;

protected:
    //! reference to the linked arrays
    arrays_type& arrays;

    //! stack of free player indices
    internal_bounded_stack<unsigned_type, arity> free_slots;

public:
    parallel_merger_adapter(const compare_type& c, arrays_type& a)
        : cmp(c), k(1), logK(0), arrays(a)
    {
        // verify strict weak ordering
        assert(!cmp(cmp.min_value(), cmp.min_value()));
    }

    void initialize()
    {
        // initial state: one empty player slot
        free_slots.push(0);
    }

    //! True if a is the sentinel value
    bool is_sentinel(const value_type& a) const
    {
        return !(cmp(cmp.min_value(), a)); // a <= cmp.min_value()
    }

    //! Allocate a free slot for a new player.
    unsigned_type new_player()
    {
        // get a free slot
        if (free_slots.empty()) {
            // tree is too small, attempt to enlarge
            double_k();
        }

        assert(!free_slots.empty());
        unsigned_type index = free_slots.top();
        free_slots.pop();

        return index;
    }

    //! Free a finished player's slot
    void free_player(unsigned_type slot)
    {
        free_slots.push(slot);
    }

    //! Whether there is still space for new array
    bool is_space_available() const
    {
        return (k < arity) || !free_slots.empty();
    }

    //! Initial call to recursive update_on_insert
    void update_on_insert(unsigned_type /* node */,
                          const value_type& /* newKey */,
                          unsigned_type /* newIndex */)
    { }

    //! make the tree twice as wide
    void double_k()
    {
        STXXL_VERBOSE1("double_k (before) k=" << k << " logK=" << logK << " arity=" << arity << " max_arity=" << max_arity << " #free=" << free_slots.size());
        assert(k > 0);
        assert(k < arity);
        assert(free_slots.empty());                    // stack was free (probably not needed)

        // make all new entries free and push them on the free stack
        for (unsigned_type i = 2 * k - 1; i >= k; i--) //backwards
        {
            arrays.make_array_sentinel(i);
            if (i < arity)
                free_slots.push(i);
        }

        // double the size
        k *= 2;
        logK++;

        STXXL_VERBOSE1("double_k (after)  k=" << k << " logK=" << logK << " arity=" << arity << " max_arity=" << max_arity << " #free=" << free_slots.size());
        assert(!free_slots.empty());
        assert(k <= max_arity);
    }

    //! compact nonempty segments in the left half of the tree
    void compact_tree()
    {
        STXXL_VERBOSE3("compact_tree (before) k=" << k << " logK=" << logK << " #free=" << free_slots.size());
        assert(logK > 0);

        // compact all nonempty segments to the left
        unsigned_type last_empty = 0;
        for (unsigned_type pos = 0; pos < k; pos++)
        {
            if (!arrays.is_array_empty(pos))
            {
                assert(arrays.is_array_allocated(pos));
                if (pos != last_empty)
                {
                    assert(!arrays.is_array_allocated(last_empty));
                    arrays.swap_arrays(last_empty, pos);
                }
                ++last_empty;
            }
            /*
              else
              {
              if(segment[pos])
              {
              STXXL_VERBOSE2("int_arrays::compact_tree() deleting segment "<<pos<<
              " address: "<<segment[pos]<<" size: "<<segment_size[pos]);
              delete [] segment[pos];
              segment[pos] = 0;
              mem_cons_ -= segment_size[pos];
              }
              }*/
        }

        // half degree as often as possible
        while ((k > 1) && last_empty <= (k / 2))
        {
            k /= 2;
            logK--;
        }

        // overwrite garbage and compact the stack of free segment indices
        free_slots.clear(); // none free
        for ( ; last_empty < k; last_empty++)
        {
            assert(!arrays.is_array_allocated(last_empty));
            arrays.make_array_sentinel(last_empty);
            if (last_empty < arity)
                free_slots.push(last_empty);
        }

        STXXL_VERBOSE3("compact_tree (after)  k=" << k << " logK=" << logK << " #free=" << free_slots.size());
    }

    //! compact tree if it got considerably smaller
    void maybe_compact()
    {
        const unsigned_type num_segments_used = k - free_slots.size();
        const unsigned_type num_segments_trigger = k - (3 * k / 5);
        // using k/2 would be worst case inefficient (for large k)
        // for k \in {2, 4, 8} the trigger is k/2 which is good
        // because we have special mergers for k \in {1, 2, 4}
        // there is also a special 3-way-merger, that will be
        // triggered if k == 4 && is_array_atsentinel(3)
        STXXL_VERBOSE3("int_merger  compact? k=" << k << " #used=" << num_segments_used
                                                 << " <= #trigger=" << num_segments_trigger << " ==> "
                                                 << ((k > 1 && num_segments_used <= num_segments_trigger) ? "yes" : "no ")
                                                 << " || "
                                                 << ((k == 4 && !free_slots.empty() && !arrays.is_array_empty(3)) ? "yes" : "no ")
                                                 << " #free=" << free_slots.size());
        if (k > 1 &&
            ((num_segments_used <= num_segments_trigger) ||
             (k == 4 && !free_slots.empty() && !arrays.is_array_empty(3))))
        {
            compact_tree();
        }
    }

    void swap(parallel_merger_adapter& obj)
    {
        std::swap(free_slots, obj.free_slots);
    }
};
#endif // STXXL_PARALLEL && (STXXL_PARALLEL_PQ_MULTIWAY_MERGE_EXTERNAL || STXXL_PARALLEL_PQ_MULTIWAY_MERGE_INTERNAL)

} // namespace priority_queue_local

//! \}

STXXL_END_NAMESPACE

#endif // !STXXL_CONTAINERS_PQ_MERGERS_HEADER
// vim: et:ts=4:sw=4
