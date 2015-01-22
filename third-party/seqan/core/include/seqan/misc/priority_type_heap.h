// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2013, Knut Reinert, FU Berlin
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
/**
.Spec.PriorityHeap
..cat:Miscellaneous
..general:Class.PriorityType
..summary:Stores the priority data on a heap.
..signature:PriorityType<TValue, TLess, PriorityHeap> >
..include:seqan/misc.h
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
		
//	inline PriorityType const &
//	operator = (PriorityType const & other_)
//	{
//		less = other_.less;
//		heap = other_.heap;
//		return *this;
//	}
	
}; // class PriorityType




// Empty the priority queue
///.Function.clear.param.object.type:Class.PriorityType
///.Function.clear.class:Class.PriorityType
template <typename TValue, typename TLess>
inline void 
clear (PriorityType<TValue,TLess, PriorityHeap> & me)
{
	clear(me.heap); 
}

// true if priority queue is empty 
///.Function.empty.param.object.type:Class.PriorityType
///.Function.empty.class:Class.PriorityType
template <typename TValue, typename TLess>
inline bool 
empty(PriorityType<TValue, TLess, PriorityHeap> const & me) 
{
SEQAN_CHECKPOINT
	return empty(me.heap); 
}

// Number of elements in the priority queue
///.Function.length.param.object.type:Class.PriorityType
///.Function.length.class:Class.PriorityType
template <typename TValue, typename TLess>
inline typename Size<PriorityType<TValue, TLess, PriorityHeap> >::Type
length( PriorityType<TValue, TLess, PriorityHeap> const & me)
{ 
SEQAN_CHECKPOINT
	return length(me.heap);
}






// Return the `best' element
/**
.Function.PriorityType#top:
..summary:Reference to the item with the highest priority.
..cat:Content Manipulation
..signature:top(object)
..class:Class.PriorityType
..param.object:A priority queue.
...type:Class.PriorityType
..remarks:To delete this item and adjust the priority queue use @Function.PriorityType#pop@.
..see:Function.PriorityType#pop
..see:Function.PriorityType#push
..include:seqan/misc.h
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
/**
.Function.adjustTop
..cat:Miscellaneous
..signature:adjustTop(object)
..class:Class.PriorityType
..summary:Adjusts the priority of the first item.
..param.object
...type:Class.PriorityType
..include:seqan/misc.h
*/
template <typename TValue, typename TLess>
inline void 
adjustTop (PriorityType<TValue, TLess, PriorityHeap> & me)	// so könnte man es dann auch nennen
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
/**
.Function.PriorityType#push:
..summary:Inserts a new item and adjusts the priority queue if necessary.
..cat:Content Manipulation
..signature:push(object, element)
..class:Class.PriorityType
..param.object:A priority queue.
...type:Class.PriorityType
..param.element:The item to be inserted in the priority queue.
...metafunction:Metafunction.Value
..remarks:The result of this operation is stored in $object$.
..see:Function.PriorityType#top
..see:Function.PriorityType#pop
..include:seqan/misc.h
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
/**
.Function.PriorityType#pop:
..summary:Deletes item with the highest priority and adjusts the priority queue.
..cat:Content Manipulation
..signature:pop(object)
..class:Class.PriorityType
..param.object:A priority queue.
...type:Class.PriorityType
..remarks:This function only deletes this item, but does not return it. To access the item use @Function.PriorityType#top@.
..see:Function.PriorityType#top
..see:Function.PriorityType#push
..include:seqan/misc.h
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

///.Metafunction.Size.param.T.type:Class.PriorityType
///.Metafunction.Size.class:Class.PriorityType
template < typename TValue, typename TLess>
struct Size<PriorityType<TValue, TLess, PriorityHeap> >
{
	typedef typename Size<typename PriorityType<TValue, TLess, PriorityHeap>::THeap>::Type Type;
};

///.metafunction.value.param.t.type:class.prioritytype
///.metafunction.value.class:class.prioritytype
template < typename TValue, typename TLess>
struct Value<PriorityType<TValue, TLess, PriorityHeap> >
{
	typedef TValue Type;
};




////////////////////////
// debug
//template <typename TValue, typename THeap, typename TLess>
//void check(PriorityType<TValue, THeap, TLess> & me) { // debug
//	typedef typename Size<PriorityType<TValue, THeap, TLess> >::Type TSize;
//	bool okay = true;
//	for ( TSize i = 1; i < length(me.heap)-1; ++i )
//		if ( me.less ( me.heap[(i-1)/2], me.heap[i] ) ) {
//			cout << '\n' << (i-1)/2 << " < " << i << " : "<< (me.heap[(i-1)/2]).value_ << " !< " << (me.heap[i]).value_;
//			okay = false;
//		}
//	if ( okay )
//		cout << " ... seems okay\n";
//	else
//		cout << "\n... there were errors!\n";
//}

	//////////////////////////////////////////////////



//////////////////////////////////////////////////////////////////////////////

} //namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
