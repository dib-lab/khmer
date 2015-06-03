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

// TODO(holtgrew): Remove, nobody implements such simple heaps as trees.

#ifndef SEQAN_HEADER_GRAPH_ALGORITHM_HEAP_TREE_H
#define SEQAN_HEADER_GRAPH_ALGORITHM_HEAP_TREE_H

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////
// HeapTree Specs
//////////////////////////////////////////////////////////////////////////////

template<typename TSpec = Default>
struct KeylessHeap;


template<typename TSpec = Default>
struct KeyedHeap;


//////////////////////////////////////////////////////////////////////////////
// Default HeapTree
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

template<typename TValue = unsigned int, typename TPredicate = std::less<unsigned int>, typename TSpec = KeylessHeap<> >
class HeapTree;

//////////////////////////////////////////////////////////////////////////////

template<typename TValue, typename TPredicate, typename TSpec>
class HeapTree<TValue, TPredicate, KeylessHeap<TSpec> >
{
public:
    typedef typename Size<HeapTree>::Type TSize;
    String<TValue> data_value;
    TSize data_heap_size;
    TPredicate data_predicate;



    //////////////////////////////////////////////////////////////////////////////

    HeapTree() : data_heap_size(0) {
        SEQAN_CHECKPOINT
    }

    //////////////////////////////////////////////////////////////////////////////

    HeapTree(HeapTree const& _other) {
        SEQAN_CHECKPOINT
        data_value = _other.data_value;
        data_heap_size = _other.data_heap_size;
        data_predicate = _other.data_predicate;
    }

    //////////////////////////////////////////////////////////////////////////////

    ~HeapTree() {
        SEQAN_CHECKPOINT
    }


    //////////////////////////////////////////////////////////////////////////////

    HeapTree& operator=(HeapTree const& _other) {
        SEQAN_CHECKPOINT
        if (this == &_other) return *this;
        data_value = _other.data_value;
        data_heap_size = _other.data_heap_size;
        data_predicate = _other.data_predicate;
        return *this;
    }
};


//////////////////////////////////////////////////////////////////////////////

template<typename TValue, typename TPredicate, typename TSpec>
struct Value<HeapTree<TValue, TPredicate, TSpec> >
{
    typedef TValue Type;
};

//////////////////////////////////////////////////////////////////////////////

template<typename TValue, typename TPredicate, typename TSpec>
struct Value<HeapTree<TValue, TPredicate, TSpec> const>
{
    typedef TValue Type;
};


//////////////////////////////////////////////////////////////////////////////

template<typename TSize>
inline TSize _left(TSize i) {
    SEQAN_CHECKPOINT
    return i << 1;
}

//////////////////////////////////////////////////////////////////////////////

template<typename TSize>
inline TSize _right(TSize i) {
    SEQAN_CHECKPOINT
    return ((i << 1) + 1);
}

//////////////////////////////////////////////////////////////////////////////

template<typename TSize>
inline TSize _parent(TSize i) {
    SEQAN_CHECKPOINT
    return (i >> 1);
}

//////////////////////////////////////////////////////////////////////////////

template<typename TValue, typename TPredicate, typename TSpec>
inline bool
empty(HeapTree<TValue, TPredicate, TSpec> const& mHeap)
{
    SEQAN_CHECKPOINT
    return (mHeap.data_heap_size == 0);
}

//////////////////////////////////////////////////////////////////////////////

template<typename TValue, typename TPredicate, typename TSpec>
inline typename Size<HeapTree<TValue, TPredicate, TSpec> >::Type
length(HeapTree<TValue, TPredicate, TSpec> const& mHeap)
{
    SEQAN_CHECKPOINT
    return mHeap.data_heap_size;
}

/////////////////////////////////////////////////////////////////////////////

template<typename TValue, typename TPredicate, typename TSpec>
inline TValue
heapRoot(HeapTree<TValue, TPredicate, TSpec>& mHeap)
{
    return value(mHeap.data_value, 1);
}

/////////////////////////////////////////////////////////////////////////////

template<typename TValue, typename TPredicate, typename TSpec>
inline TValue
heapRoot(HeapTree<TValue, TPredicate, TSpec> const& mHeap)
{
    return value(mHeap.data_value, 1);
}

/////////////////////////////////////////////////////////////////////////////

template<typename TValue, typename TPredicate, typename TSpec, typename TObject>
inline void
heapInsert(HeapTree<TValue, TPredicate, TSpec>& mHeap,
           TObject object)
{
    SEQAN_CHECKPOINT
    typedef HeapTree<TValue, TPredicate, TSpec> THeapTree;
    typedef typename Size<THeapTree>::Type TSize;

    ++mHeap.data_heap_size;
    resize(mHeap.data_value, mHeap.data_heap_size + 1, Generous() );
    TSize i = mHeap.data_heap_size;
    _insertObject(mHeap, i, object);
    while ((i>1) && (!(mHeap.data_predicate(value(mHeap.data_value, _parent(i)), value(mHeap.data_value, i))))) {
        _swapObjects(mHeap, i, _parent(i));
        i = _parent(i);
    }
}

/////////////////////////////////////////////////////////////////////////////

template<typename TValue, typename TPredicate, typename TSpec, typename TSize>
inline void
heapify(HeapTree<TValue, TPredicate, TSpec>& mHeap,
        TSize i)
{
    SEQAN_CHECKPOINT
    TSize l = _left(i);
    TSize r = _right(i);
    TSize largest = i;
    if ((l <= (TSize) mHeap.data_heap_size) && (mHeap.data_predicate(value(mHeap.data_value, l), value(mHeap.data_value, i)))) largest = l;
    if ((r <= (TSize) mHeap.data_heap_size) && (mHeap.data_predicate(value(mHeap.data_value, r), value(mHeap.data_value, largest)))) largest = r;
    if (largest != i) {
        _swapObjects(mHeap, i, largest);
        heapify(mHeap, largest);
    }
}

//////////////////////////////////////////////////////////////////////////////

template<typename TValue, typename TPredicate, typename TSpec, typename TITBegin, typename TITEnd>
inline void
buildHeap(HeapTree<TValue, TPredicate, TSpec>& mHeap,
          TITBegin itBeg,
          TITEnd itEnd)
{
    SEQAN_CHECKPOINT
    typedef HeapTree<TValue, TPredicate, TSpec> THeapTree;
    typedef typename Size<THeapTree>::Type TSize;

    // Leave out the first element
    resize(mHeap.data_value, (itEnd - itBeg) + 1);
    value(mHeap.data_value, 0) = TValue();
    mHeap.data_heap_size = 0;
    for(;itBeg != itEnd; goNext(itBeg), ++mHeap.data_heap_size) {
        _insertObject(mHeap, mHeap.data_heap_size + 1, value(itBeg));
    }
    for(TSize i = (mHeap.data_heap_size / 2); i>0; --i) heapify(mHeap, i);
}

//////////////////////////////////////////////////////////////////////////////

template<typename TValue, typename TPredicate, typename TSpec, typename TIndex1, typename TIndex2>
inline void
_swapObjects(HeapTree<TValue, TPredicate, KeylessHeap<TSpec> >& mHeap,
             TIndex1 i,
             TIndex2 j)
{
    SEQAN_CHECKPOINT
    TValue tmp = value(mHeap.data_value, i);
    value(mHeap.data_value, i) = value(mHeap.data_value, j);
    value(mHeap.data_value, j) = tmp;
}

//////////////////////////////////////////////////////////////////////////////

template<typename TValue, typename TPredicate, typename TSpec, typename TIndex, typename TObject>
inline void
_insertObject(HeapTree<TValue, TPredicate, KeylessHeap<TSpec> >& mHeap,
              TIndex i,
              TObject obj)
{
    SEQAN_CHECKPOINT
    value(mHeap.data_value, i) = obj;
}

/////////////////////////////////////////////////////////////////////////////

template<typename TValue, typename TPredicate, typename TSpec>
inline TValue
heapExtractRoot(HeapTree<TValue, TPredicate, KeylessHeap<TSpec> >& mHeap)
{
    SEQAN_CHECKPOINT
    TValue rootVal = value(mHeap.data_value, 1);
    value(mHeap.data_value, 1) = value(mHeap.data_value, mHeap.data_heap_size);
    --mHeap.data_heap_size;
    heapify(mHeap, 1);
    return rootVal;
}

//////////////////////////////////////////////////////////////////////////////

template<typename TValue, typename TPredicate, typename TSpec>
inline void
clear(HeapTree<TValue, TPredicate, TSpec>& mHeap)
{
    clear(mHeap.data_value);
    mHeap.data_heap_size = 0;
}

//////////////////////////////////////////////////////////////////////////////

template<typename TITBegin, typename TITEnd, typename TPredicate>
inline void
heapSort(TITBegin itBeg,
         TITEnd itEnd,
         TPredicate)
{
    SEQAN_CHECKPOINT
    typedef typename Value<TITBegin>::Type TValue;
    typedef typename Size<TITBegin>::Type TSize;
    HeapTree<TValue, TPredicate> mHeap;
    TITBegin itFill = itBeg;
    buildHeap(mHeap, itBeg, itEnd);
    for(TSize i = mHeap.data_heap_size; i>1; --i) {
        value(itFill) = value(mHeap.data_value, 1);
        goNext(itFill);
        value(mHeap.data_value, 1) = value(mHeap.data_value, i);
        --mHeap.data_heap_size;
        heapify(mHeap, 1);
    }
    value(itFill) = value(mHeap.data_value, 2);
}

//////////////////////////////////////////////////////////////////////////////

template<typename TITBegin, typename TITEnd>
inline void
heapSort(TITBegin itBeg,
         TITEnd itEnd)
{
    SEQAN_CHECKPOINT
    typedef typename Value<TITBegin>::Type TValue;
    heapSort(itBeg, itEnd, std::less<TValue>());
}


//////////////////////////////////////////////////////////////////////////////

template<typename TKey, typename TValue, typename TPredicate>
struct LessPairI2WithFunctor :
    public std::unary_function<Pair<TKey, TValue>, bool>
{
    inline bool
    operator() (Pair<TKey, TValue> const& a1, Pair<TKey, TValue> const& a2) {
        TPredicate private_Predicate;
        return private_Predicate(a1.i2, a2.i2);
    }
};


//////////////////////////////////////////////////////////////////////////////

template<typename TKey, typename TValue, typename TPredicate, typename TSpec>
class HeapTree<Pair<TKey, TValue>, TPredicate, KeyedHeap<TSpec> >
{
public:
    typedef typename Size<HeapTree>::Type TSize;
    String<Pair<TKey, TValue> > data_value;
    std::map<TKey, TSize> data_map;
    TSize data_heap_size;
    LessPairI2WithFunctor<TKey, TValue, TPredicate> data_predicate;


    //////////////////////////////////////////////////////////////////////////////

    HeapTree() : data_heap_size(0) {
        SEQAN_CHECKPOINT
    }

    //////////////////////////////////////////////////////////////////////////////

    HeapTree(HeapTree const& _other) {
        SEQAN_CHECKPOINT
        data_value = _other.data_value;
        data_map = _other.data_map;
        data_heap_size = _other.data_heap_size;
        data_predicate = _other.data_predicate;
    }

    //////////////////////////////////////////////////////////////////////////////

    ~HeapTree() {
        SEQAN_CHECKPOINT
    }


    //////////////////////////////////////////////////////////////////////////////

    HeapTree& operator=(HeapTree const& _other) {
        SEQAN_CHECKPOINT
        data_value = _other.data_value;
        data_map = _other.data_map;
        data_heap_size = _other.data_heap_size;
        data_predicate = _other.data_predicate;
        return *this;
    }
};


//////////////////////////////////////////////////////////////////////////////

template<typename TKey, typename TValue, typename TPredicate, typename TSpec, typename TIndex1, typename TIndex2>
inline void
_swapObjects(HeapTree<Pair<TKey, TValue>, TPredicate, KeyedHeap<TSpec> >& mHeap,
             TIndex1 i,
             TIndex2 j)
{
    SEQAN_CHECKPOINT
    typedef HeapTree<Pair<TKey, TValue>, TPredicate, TSpec> THeapTree;
    typedef typename Size<THeapTree>::Type TSize;
    typedef std::map<TKey, TSize> TMap;
    typedef typename TMap::iterator TMapIter;
    // Swap keys
    TMapIter pos1 = mHeap.data_map.find(value(mHeap.data_value, i).i1);
    TMapIter pos2 = mHeap.data_map.find(value(mHeap.data_value, j).i1);
    TSize tmpVal = pos1->second;
    pos1->second = pos2->second;
    pos2->second = tmpVal;

    // Swap object in heap tree
    Pair<TKey, TValue> tmp = value(mHeap.data_value, i);
    value(mHeap.data_value, i) = value(mHeap.data_value, j);
    value(mHeap.data_value, j) = tmp;
}

//////////////////////////////////////////////////////////////////////////////

template<typename TKey, typename TValue, typename TPredicate, typename TSpec, typename TIndex, typename TObject>
inline void
_insertObject(HeapTree<Pair<TKey, TValue>, TPredicate, KeyedHeap<TSpec> >& mHeap,
              TIndex i,
              TObject obj)
{
    SEQAN_CHECKPOINT
    typedef HeapTree<Pair<TKey, TValue>, TPredicate, TSpec> THeapTree;
    typedef typename Size<THeapTree>::Type TSize;
    value(mHeap.data_value, i) = obj;
    if (!mHeap.data_map.insert(std::pair<TKey, TSize>(obj.i1, i)).second) {
        mHeap.data_map.find(obj.i1)->second = i;
    };
}

/////////////////////////////////////////////////////////////////////////////

template<typename TKey, typename TValue, typename TPredicate, typename TSpec>
inline Pair<TKey, TValue>
heapExtractRoot(HeapTree<Pair<TKey, TValue>, TPredicate, KeyedHeap<TSpec> >& mHeap)
{
    SEQAN_CHECKPOINT
    Pair<TKey, TValue> rootVal = value(mHeap.data_value, 1);
    mHeap.data_map.erase(mHeap.data_map.find(rootVal.i1));
    if (mHeap.data_heap_size != 1) {
        value(mHeap.data_value, 1) = value(mHeap.data_value, mHeap.data_heap_size);
        mHeap.data_map.find(value(mHeap.data_value, 1).i1)->second = 1;
    }
    --mHeap.data_heap_size;
    heapify(mHeap, 1);
    return rootVal;
}

/////////////////////////////////////////////////////////////////////////////

template<typename TKey, typename TValue, typename TPredicate, typename TSpec, typename TKey1, typename TValue2>
inline void
heapChangeValue(HeapTree<Pair<TKey, TValue>, TPredicate, KeyedHeap<TSpec> >& mHeap,
                TKey1 key,
                TValue2 val)
{
    SEQAN_CHECKPOINT
    typedef HeapTree<Pair<TKey, TValue>, TPredicate, TSpec> THeapTree;
    typedef typename Size<THeapTree>::Type TSize;
    TSize i = mHeap.data_map.find(key)->second;
    Pair<TKey, TValue> obj = Pair<TKey, TValue>(key,val);
    if (!(mHeap.data_predicate(obj, value(mHeap.data_value, i)))) {
        _insertObject(mHeap, i, obj);
        heapify(mHeap, i);
    } else {
        _insertObject(mHeap, i, obj);
        while ((i>1) && (!mHeap.data_predicate(value(mHeap.data_value, _parent(i)), value(mHeap.data_value, i)))) {
            _swapObjects(mHeap, i, _parent(i));
            i = _parent(i);
        }
    }
}

/////////////////////////////////////////////////////////////////////////////

template<typename TKey, typename TValue, typename TPredicate, typename TSpec, typename TKey1>
inline TValue
heapGetValue(HeapTree<Pair<TKey, TValue>, TPredicate, KeyedHeap<TSpec> >& mHeap,
             TKey1 key)
{
    SEQAN_CHECKPOINT
    typedef HeapTree<Pair<TKey, TValue>, TPredicate, TSpec> THeapTree;
    typedef typename Size<THeapTree>::Type TSize;
    TSize i = mHeap.data_map.find(key)->second;
    return value(mHeap.data_value, i).i2;
}

//////////////////////////////////////////////////////////////////////////////

template<typename TKey, typename TValue, typename TPredicate, typename TSpec>
inline void
clear(HeapTree<Pair<TKey, TValue>, TPredicate, KeyedHeap<TSpec> >& mHeap)
{
    clear(mHeap.data_value);
    mHeap.data_map.clear();
    mHeap.data_heap_size = 0;
}



}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
