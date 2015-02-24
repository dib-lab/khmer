// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2015, Knut Reinert, FU Berlin
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//     * Neither the name of Knut Reinert or the FU Berlin nor the names of
//       its contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL KNUT REINERT OR THE FU BERLIN BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
// LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
// OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
// DAMAGE.
//
// ==========================================================================

#ifndef SEQAN_HEADER_PRIORITY_TYPE_TREE_H
#define SEQAN_HEADER_PRIORITY_TYPE_TREE_H

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////

/*!
 * @class PriorityHeap
 * @headerfile <seqan/misc/priority_type_heap.h>
 * @extends PriorityType
 * @brief Stores the priority data on a heap.
 *
 * @signature template <[typename TValue[, typename TLess]]>
 *            class PriorityType<TValue, TLess, PriorityHeap>;
 *
 * @tparam TValue The value type.  Default: <tt>int</tt>.
 * @tparam TLess  The less-than comparator.  Default: <tt>std::less&lt;TValue&gt;</tt>.
 */

template < typename TValue, typename TLess>
class PriorityType<TValue,TLess,PriorityHeap>
{
public:
    typedef String<TValue> THeap;
//____________________________________________________________________________

    TLess less;
    THeap heap;
//____________________________________________________________________________

    inline PriorityType()
    {
SEQAN_CHECKPOINT
    }

    inline explicit PriorityType(TLess const & _less):
        less(_less)
    {
SEQAN_CHECKPOINT
    }

    inline PriorityType(PriorityType const & _other):
        less(_other.less),
        heap(_other.heap)
    {
SEQAN_CHECKPOINT
    }

//    inline PriorityType const &
//    operator = (PriorityType const & other_)
//    {
//        less = other_.less;
//        heap = other_.heap;
//        return *this;
//    }

}; // class PriorityType



/*!
 * @fn PriorityType#clear
 * @brief Remove all elements from the PriorityType.
 *
 * @signature void clear(pq);
 *
 * @param[in,out] pq PriorityType to clear.
 */

// Empty the priority queue
template <typename TValue, typename TLess>
inline void
clear (PriorityType<TValue,TLess, PriorityHeap> & me)
{
    clear(me.heap);
}

/*!
 * @fn PriorityType#empty
 * @headerfile <seqan/misc/priority_type_heap.h>
 * @brief Query priority queue for beging empty.
 *
 * @signature bool empty(pq);
 *
 * @param[in] pq The @link PriorityType @endlink to query.
 *
 * @return bool <tt>true</tt> if <tt>pq</tt> queue is empty.
 */

template <typename TValue, typename TLess>
inline bool
empty(PriorityType<TValue, TLess, PriorityHeap> const & me)
{
SEQAN_CHECKPOINT
    return empty(me.heap);
}

/*!
 * @fn PriorityType#length
 * @headerfile <seqan/misc/priority_type_heap.h>
 * @brief Return number of elements in priority queue.
 *
 * @signature TSize length(pq);
 *
 * @param[in] pq The PriorityType to query.
 * @return TSize Number of elements in priority queue.
 */

// Number of elements in the priority queue
template <typename TValue, typename TLess>
inline typename Size<PriorityType<TValue, TLess, PriorityHeap> >::Type
length( PriorityType<TValue, TLess, PriorityHeap> const & me)
{
SEQAN_CHECKPOINT
    return length(me.heap);
}






// Return the `best' element

/*!
 * @fn PriorityType#top
 * @brief Reference to the item with the highest priority.
 *
 * @signature TReference top(pq);
 *
 * @param[in] pq The PriorityType to query.
 *
 * @return TReference The result, reference to Value type.
 */

template <typename TValue, typename TLess>
inline TValue &
top(PriorityType<TValue, TLess, PriorityHeap> & me)
{
SEQAN_CHECKPOINT
    return value(me.heap, beginPosition(me.heap));
}

template <typename TValue, typename TLess>
inline TValue const &
top(PriorityType<TValue, TLess, PriorityHeap> const & me)
{
SEQAN_CHECKPOINT
    return value(me.heap, beginPosition(me.heap));
}

// Copy heap position i to heap position h.
template <typename TValue, typename TLess, typename TSize>
inline void
_copyHeapElement (PriorityType<TValue, TLess, PriorityHeap> & me, TSize i, TSize & h)
{
SEQAN_CHECKPOINT
    me.heap[h] = me.heap[i];
    h = i;
}

// Copy element to heap position h.
template <typename TValue, typename TLess, typename TSize>
inline void
_copyHeapElement (PriorityType<TValue, TLess, PriorityHeap> & me, TValue const & element, TSize h)
{
SEQAN_CHECKPOINT
    me.heap[h] = element;
}

/////////////////////////////////////////////////////////////////////////////////
//  lower priority of first element in queue

/*!
 * @fn PriorityType#adjustTop
 * @brief Adjusts the priority of the first item.
 *
 * @signature void adjustTop(pq);
 *
 * @param[in,out] pq The PriorityType to adjust.
 */

template <typename TValue, typename TLess>
inline void
adjustTop (PriorityType<TValue, TLess, PriorityHeap> & me)    // so könnte man es dann auch nennen
{
SEQAN_CHECKPOINT
    if (!empty(me.heap))
        _adjustHeapTowardLeaves (me, me.heap[0], 0, 2);
}

////////////
//?? adjustHeapElement...
/////////////

/////////////////////////////////////////////////////////////////////////////////
/// Push a new element

/*!
 * @fn PriorityType#push
 * @brief Inserts a new item and adjusts the priority queue if necessary.
 *
 * @signature void push(pq, element);
 *
 * @param[in,out] pq      The PriorityType to push to.
 * @param[in]     element The element to push.
 */

template <typename TValue, typename TLess>
inline void
push (PriorityType<TValue, TLess, PriorityHeap> & me, TValue const & element)
{
SEQAN_CHECKPOINT
    // root index is zero
    if (empty(me.heap)) {
        resize(me.heap, 1, Generous());
        _copyHeapElement (me, element, 0);
        return;
    }
    typedef typename Size<PriorityType<TValue, TLess, PriorityHeap> >::Type TSize;
    TSize h = length(me.heap);
    resize(me.heap, h + 1, Generous());
    _adjustHeapTowardRoot(me, element, h);
}

//////////////////////////////////////////////////////////////////////////////////
/// Priority got better.  Perform a cyclic shift along the tree edges toward root.
template <typename TValue, typename TLess, typename TSize>
inline void
_adjustHeapTowardRoot(
    PriorityType<TValue, TLess, PriorityHeap> & me,
    TValue const & element,
    TSize h )
{
SEQAN_CHECKPOINT
    // root index is zero
    while ( h > 0) {
        const TSize i = (h-1)/2;
        if ( me.less ( me.heap[i], element ) )
            _copyHeapElement ( me, i, h );
        else
            break;
    }
    _copyHeapElement ( me, element, h );
}

//////////////////////////////////////////////////////////////////////////////
/// Pop 'best' element

/*!
 * @fn PriorityType#pop
 * @brief Deletes item with the highest priority and adjusts the priority queue.
 *
 * @signature void push(pq);
 *
 * @param[in,out] pq      The PriorityType to pop from.
 */

template <typename TValue, typename TLess>
inline void
pop (PriorityType<TValue, TLess, PriorityHeap> & me)
{
SEQAN_CHECKPOINT
    // root index is zero
    TValue element = getValue(me.heap,endPosition(me.heap)-1);
    typedef typename Size<PriorityType<TValue, TLess, PriorityHeap> >::Type TSize;
    TSize heapsize = length(me.heap) - 1 ;
    resize(me.heap, heapsize, Generous());
    if ( heapsize > 0 )
        _adjustHeapTowardLeaves(me, element, 0, 2 );

}


//////////////////////////////////////////////////////////////////////////////////
/// Priority got worse. Perform a cyclic shift along the tree edges toward leaves.
template <typename TValue, typename TLess, typename TSize>
inline void
_adjustHeapTowardLeaves(
    PriorityType<TValue, TLess, PriorityHeap> & me,
    TValue element,
    TSize h,
    TSize i ) //für mich: h=0, i=1
{
SEQAN_CHECKPOINT
    // root index is zero
    const TSize heapsize = length(me.heap);
    TLess less = me.less;
    while ( i < heapsize )
    {
        if ( less ( element, me.heap[i] ) )
            if ( less ( me.heap[i-1], me.heap[i] ) )
                _copyHeapElement ( me, i, h );
            else
                _copyHeapElement ( me, i-1, h );
        else
            if ( less ( element, me.heap[i-1] ) )
                _copyHeapElement ( me, i-1, h );
            else
                break;
        i = 2*(h+1);
    }
    if ( i == heapsize && less ( element, me.heap[i-1] ) )
        _copyHeapElement ( me, i-1, h );
    _copyHeapElement ( me, element, h );
}


    //MetaFunctions

template < typename TValue, typename TLess>
struct Size<PriorityType<TValue, TLess, PriorityHeap> >
{
    typedef typename Size<typename PriorityType<TValue, TLess, PriorityHeap>::THeap>::Type Type;
};

template < typename TValue, typename TLess>
struct Value<PriorityType<TValue, TLess, PriorityHeap> >
{
    typedef TValue Type;
};




////////////////////////
// debug
//template <typename TValue, typename THeap, typename TLess>
//void check(PriorityType<TValue, THeap, TLess> & me) { // debug
//    typedef typename Size<PriorityType<TValue, THeap, TLess> >::Type TSize;
//    bool okay = true;
//    for ( TSize i = 1; i < length(me.heap)-1; ++i )
//        if ( me.less ( me.heap[(i-1)/2], me.heap[i] ) ) {
//            cout << '\n' << (i-1)/2 << " < " << i << " : "<< (me.heap[(i-1)/2]).value_ << " !< " << (me.heap[i]).value_;
//            okay = false;
//        }
//    if ( okay )
//        cout << " ... seems okay\n";
//    else
//        cout << "\n... there were errors!\n";
//}

    //////////////////////////////////////////////////



//////////////////////////////////////////////////////////////////////////////

} //namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
