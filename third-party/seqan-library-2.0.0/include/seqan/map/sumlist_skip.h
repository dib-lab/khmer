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

// SEQAN_NO_GENERATED_FORWARDS: no forwards are generated for this file

#ifndef SEQAN_HEADER_SUMLIST_SKIP_H
#define SEQAN_HEADER_SUMLIST_SKIP_H

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////
// Skipsumlist_
//////////////////////////////////////////////////////////////////////////////

//Spec for Map
template <unsigned int DIM, typename TSpec = Default>
struct Skipsumlist_;

//////////////////////////////////////////////////////////////////////////////
// Map<Skipsumlist_>
//////////////////////////////////////////////////////////////////////////////



/*
template <typename TValue, unsigned int DIM, typename TSpec>
struct SkipListBlockSize_<Map<TValue, Skiplist<Skipsumlist_<DIM, TSpec> > > >
{
    typedef SumList<DIM, TValue, MiniSumList< > > TMiniSumList;
    typedef Map<TValue, Skiplist<Skipsumlist_<DIM, TSpec> > > TMap;
    enum
    {
        VALUE = sizeof(TMiniSumList) + TMap::MAX_HEIGHT * sizeof(
    };
};
*/

//////////////////////////////////////////////////////////////////////////////

template <typename TValue, unsigned int DIM, typename TSpec>
class SkiplistNext<TValue, Skipsumlist_<DIM, TSpec> >
{
public:
    typedef Map<TValue, Skiplist<Skipsumlist_<DIM, TSpec> > > TMap;
    typedef SkiplistElement<TValue, Skipsumlist_<DIM, TSpec> > TElement;
    typedef SumListValues<DIM, TValue> TValues;

    TElement * data_element;
    TValues values;

    SkiplistNext()
        : values(MinimalCtor())
    {}
    SkiplistNext(NonMinimalCtor)
        : data_element(0)
    {}
    SkiplistNext(SkiplistNext const & other)
        : data_element(other.data_element)
        , values(other.values)
    {}
    ~SkiplistNext()
    {}
    SkiplistNext const & operator = (SkiplistNext const & other)
    {
        data_element = other.data_element;
        values = other.values;
        return *this;
    }
};

//////////////////////////////////////////////////////////////////////////////


template <typename TValue, unsigned int DIM, typename TSpec>
class SkiplistElement<TValue, Skipsumlist_<DIM, TSpec> >
{
public:
    typedef Map<TValue, Skiplist<Skipsumlist_<DIM, TSpec> > > TSkiplist;
    typedef SkiplistNext<TValue, Skipsumlist_<DIM, TSpec> > TNext;
    typedef SumList<DIM, TValue, MiniSumList< > > TMiniSumList;

    enum
    {
        MAX_HEIGHT = TSkiplist::MAX_HEIGHT
    };

    TMiniSumList minilist;
    TNext data_next[MAX_HEIGHT]; //note: only parts of this array available
};

//////////////////////////////////////////////////////////////////////////////

template <typename TValue, unsigned int DIM, typename TSpec>
class SkiplistPath<TValue, Skipsumlist_<DIM, TSpec> >
{
public:
    typedef Map<TValue, Skiplist<Skipsumlist_<DIM, TSpec> > > TSkiplist;
    typedef SkiplistElement<TValue, Skipsumlist_<DIM, TSpec> > TElement;
    typedef SumListValues<DIM, TValue> TValues;

    enum
    {
        MAX_HEIGHT = TSkiplist::MAX_HEIGHT
    };

    TElement * data_elements[MAX_HEIGHT];
    TValues sums[MAX_HEIGHT];

    SkiplistPath() {}
    SkiplistPath(TSkiplist & sl)
    {
        for (unsigned int i = 0; i <= sl.data_height; ++i)
        {
            data_elements[i] = & (sl.data_border);
        }
    }
};


//____________________________________________________________________________

//make path to a kind of iterator

template <typename TValue, unsigned int DIM, typename TSpec>
bool atEnd(SkiplistPath<TValue, Skipsumlist_<DIM, TSpec> > & path)
{
    return (!path.data_elements[0]);
}

template <typename TValue, unsigned int DIM, typename TSpec, typename TSkiplist>
void goNext(SkiplistPath<TValue, Skipsumlist_<DIM, TSpec> > & path,
            TSkiplist const & skiplist)
{
    typedef SkiplistPath<TValue, Skipsumlist_<DIM, TSpec> > TPath;
    typedef SkiplistElement<TValue, Skipsumlist_<DIM, TSpec> > TElement;

    TElement * next = path.data_elements[0]->data_next[0].data_element;
    for (unsigned int i = 0; i <= skiplist.data_height; ++i)
    {
        if (!path.data_elements[i]) break;
        if (path.data_elements[i]->data_next[i].data_element != next) break;

        path.sums[i] += path.data_elements[i]->data_next[i].values;
        path.data_elements[i] = next;
    }
}

//////////////////////////////////////////////////////////////////////////////

//if true, go right
//if false, go down

//find a value
template <typename TNext, typename THere, typename TValues, typename TValue>
inline bool
_skipsumlistFindGoNext(TNext const & _next,
                       THere const & /*_here*/,
                       unsigned char /*_height*/,
                       TValues const & _sum,
                       unsigned int _dim,
                       TValue const & _value)
{
    return (static_cast<TValue>((_sum[_dim] + _next.values[_dim])) <= _value); // SEARCH SEMANTICS
}
//find an element, given (sum, element_pointer)-pair
template <typename TNext, typename THere, typename TValues, typename TValue, typename TElementPtr>
inline bool
_skipsumlistFindGoNext(TNext const & _next,
                       THere const & _here,
                       unsigned char /*_height*/,
                       TValues const & _sum,
                       unsigned int _dim,
                       Pair<TValue, TElementPtr> const & pair)
{
    return (_sum[_dim] + _next.values[_dim] <= pair.i1) && (pair.i2 != _here);
}
//find the end of the skiplist
template <typename TNext, typename THere, typename TValues>
inline bool
_skipsumlistFindGoNext(TNext const & _next,
                       THere const & /*_here*/,
                       unsigned char /*_height*/,
                       TValues const & /*_sum*/,
                       unsigned int /*_dim*/,
                       GoEnd)
{
    return _next.data_element != 0;
}


template <typename TValue, unsigned int DIM, typename TSpec, typename TFind>
inline void
_skipsumlistFind(Map<TValue, Skiplist<Skipsumlist_<DIM, TSpec> > > & me,
                 TFind const & find,
                 unsigned int dim,
                 /*OUT*/ SkiplistPath<TValue, Skipsumlist_<DIM, TSpec> > & path)
{
    typedef Map<TValue, Skiplist<Skipsumlist_<DIM, TSpec> > > TMap;
    typedef SkiplistElement<TValue, Skipsumlist_<DIM, TSpec> > TElement;
    typedef SkiplistNext<TValue, Skipsumlist_<DIM, TSpec> > TNext;
    typedef SumListValues<DIM, TValue> TValues;

    typedef typename Size<TMap>::Type TSize;

    TElement * here = & me.data_border;

    TValues sum;

    for (int i = me.data_height; i >= 0; --i)
    {
        while (true)
        {
            TNext & next = here->data_next[i];
            if (!next.data_element) break;

            if (!_skipsumlistFindGoNext(next, here, i, sum, dim, find)) break;

            here = next.data_element;
            sum += next.values;
        }
        path.data_elements[i] = here;
        path.sums[i] = sum;
    }
}
template <typename TValue, unsigned int DIM, typename TSpec, typename TFind>
inline void
_skipsumlistFind(Map<TValue, Skiplist<Skipsumlist_<DIM, TSpec> > > const & me,
                 TFind const & find,
                 unsigned int dim,
                 /*OUT*/ SkiplistPath<TValue, Skipsumlist_<DIM, TSpec> > & path)
{//destroy const
    typedef Map<TValue, Skiplist<Skipsumlist_<DIM, TSpec> > > TMap;
    _skipsumlistFind(const_cast<TMap &>(me), find, dim, path);
}


template <typename TValue, unsigned int DIM, typename TSpec>
inline void
assign(Map<TValue, Skiplist< Skipsumlist_<DIM, TSpec> > > & target,
       Map<TValue, Skiplist< Skipsumlist_<DIM, TSpec> > > const & source_)
{
    typedef Map<TValue, Skiplist< Skipsumlist_<DIM, TSpec> > > TSkiplist;
    typedef SkiplistPath<TValue, Skipsumlist_<DIM, TSpec> > TPath;
    typedef SkiplistElement<TValue, Skipsumlist_<DIM, TSpec> > TElement;
    //typedef typename Iterator<TSkiplist>::Type TIterator;
    //typedef typename Value<TSkiplist>::Type TValue2;

    TSkiplist & source = const_cast<TSkiplist &>(source_); //oh, I'm so damned lazy :-P

    //first clear target
    clear(target);

    //init path variable that is used as iterator-like construct to build up target
    //and copy link values outgoing from target border
    TPath target_path(target);
    target.data_height = source.data_height;
    for (unsigned int i = 0; i <= target.data_height; ++i)
    {
        target_path.data_elements[i] = & target.data_border;
        target.data_border.data_next[i].values = source.data_border.data_next[i].values;
    }

    //copy the first minilist
    target.data_border.minilist = source.data_border.minilist;

    TPath source_path(source);
    //copy rest of the values
    while (true)
    {
        goNext(source_path, source);

        if (atEnd(source_path)) break;

        //determine height of current tower
        unsigned char height;
        for (height = 1; height < TSkiplist::MAX_HEIGHT; ++height)
        {
            if (source_path.data_elements[height] != source_path.data_elements[0]) break;
        }
        --height;

        //create new element in target and copy minilist
        TElement & el = _skiplistAllocateElement(target, height);
        el.minilist = source_path.data_elements[0]->minilist;

        //insert element in target
        for (int i = 0; i <= height; ++i)
        {
            el.data_next[i].data_element = 0;
            el.data_next[i].values = source_path.data_elements[i]->data_next[i].values;
            target_path.data_elements[i]->data_next[i].data_element = & el;
            target_path.data_elements[i] = & el;
        }

    }

    //set length
    target.data_length = length(source);
}

//////////////////////////////////////////////////////////////////////////////

template <typename TValue, unsigned int DIM, typename TSpec>
inline void
clear(Map<TValue, Skiplist< Skipsumlist_<DIM, TSpec> > > & me)
{
    typedef Map<TValue, Skiplist< Skipsumlist_<DIM, TSpec> > > TSkiplist;

    me.data_mem_begin = me.data_mem_end = 0;
    me.data_length = 0;
    me.data_height = 0;

    for (unsigned char i = 0; i < TSkiplist::MAX_HEIGHT; ++i)
    {
        me.data_recycle[i] = 0;
        valueConstruct(me.data_border.data_next + i, NonMinimalCtor());
    }
    clear(me.data_border.minilist); //thats new

    clear(value(me.data_allocator));
}

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
// SumList<Skipsumlist_>
//////////////////////////////////////////////////////////////////////////////

template <unsigned int DIM, typename TValue, typename TSpec>
class SumList<DIM, TValue, SkipSumList<TSpec> >
{
public:
    typedef Map<TValue, Skiplist<Skipsumlist_<DIM, TSpec> > > TMap;
    typedef typename Size<SumList>::Type TSize;
    typedef SumListValues<DIM, TValue> TValues;

    TMap map;
    TSize length;
    TValues sum;

    SumList()
        : length(0)
    {}
    SumList(SumList const & other)
        : length(0)
    {
        assign(*this, other);
    }
    ~SumList()
    {}
    SumList const & operator = (SumList const & other)
    {
        assign(*this, other);
        return *this;
    }
};

//////////////////////////////////////////////////////////////////////////////


template <unsigned int DIM, typename TValue, typename TSpec>
inline void
assign(SumList<DIM, TValue, SkipSumList<TSpec> > & target,
       SumList<DIM, TValue, SkipSumList<TSpec> > const & source)
{
    target.map = source.map;
    target.length = source.length;
    target.sum = source.sum;
}


//////////////////////////////////////////////////////////////////////////////

template <unsigned int DIM, typename TValue, typename TSpec>
inline typename Size< SumList<DIM, TValue, SkipSumList<TSpec> > >::Type
length(SumList<DIM, TValue, SkipSumList<TSpec> > const & me)
{
    return me.length;
}

//////////////////////////////////////////////////////////////////////////////

template <unsigned int DIM, typename TValue, typename TSpec>
inline typename Value< SumList<DIM, TValue, SkipSumList<TSpec> > >::Type
getSum(SumList<DIM, TValue, SkipSumList<TSpec> > const & me,
       unsigned int dim)
{
    return me.sum[dim];
}

//////////////////////////////////////////////////////////////////////////////

template <unsigned int DIM, typename TValue, typename TSpec>
inline void
clear(SumList<DIM, TValue, SkipSumList<TSpec> > & me)
{
    clear(me.map);
    me.length = 0;
    clear(me.sum);
}

//////////////////////////////////////////////////////////////////////////////

template <unsigned int DIM, typename TValue, typename TSpec>
inline typename Iterator< SumList<DIM, TValue, SkipSumList<TSpec> > >::Type
begin(SumList<DIM, TValue, SkipSumList<TSpec> > & me)
{
    typedef SumList<DIM, TValue, SkipSumList<TSpec> >  TMe;
    typedef typename Iterator<TMe>::Type TIterator;
    return TIterator(me);
}
template <unsigned int DIM, typename TValue, typename TSpec>
inline typename Iterator< SumList<DIM, TValue, SkipSumList<TSpec> > const>::Type
begin(SumList<DIM, TValue, SkipSumList<TSpec> > const & me)
{
    typedef SumList<DIM, TValue, SkipSumList<TSpec> > const TMe;
    typedef typename Iterator<TMe>::Type TIterator;
    return TIterator(me);
}

//////////////////////////////////////////////////////////////////////////////

template <unsigned int DIM, typename TValue, typename TSpec>
inline typename Iterator< SumList<DIM, TValue, SkipSumList<TSpec> > >::Type
end(SumList<DIM, TValue, SkipSumList<TSpec> > & me)
{
    typedef SumList<DIM, TValue, SkipSumList<TSpec> >  TMe;
    typedef typename Iterator<TMe>::Type TIterator;
    return TIterator(me, GoEnd());
}
template <unsigned int DIM, typename TValue, typename TSpec>
inline typename Iterator< SumList<DIM, TValue, SkipSumList<TSpec> > const>::Type
end(SumList<DIM, TValue, SkipSumList<TSpec> > const & me)
{
    typedef SumList<DIM, TValue, SkipSumList<TSpec> > const TMe;
    typedef typename Iterator<TMe>::Type TIterator;
    return TIterator(me, GoEnd());
}

//////////////////////////////////////////////////////////////////////////////

//NOTE: the first argument is Sumlist, not Map!
template <typename TValue, unsigned int DIM, typename TSpec, typename TSpec2>
inline unsigned char
_skiplistCreateHeight(SumList<DIM, TValue, SkipSumList<TSpec> > & me,
                      SkiplistPath<TValue, TSpec2 > & path) //extend path if height is increased
{
    typedef Map<TValue, Skiplist<Skipsumlist_<DIM, TSpec> > > TSkiplist;

    unsigned char height = geomRand<unsigned char>();
    if (height >= TSkiplist::MAX_HEIGHT) height = TSkiplist::MAX_HEIGHT-1;

    if (height > me.map.data_height)
    {
        for (unsigned char i = me.map.data_height + 1; i <= height; ++i)
        {
            path.data_elements[i] = & me.map.data_border;
            //path.sums[i] = me.sum;
            me.map.data_border.data_next[i].values = me.sum;
        }
        me.map.data_height = height;
    }

    return height;
}

//////////////////////////////////////////////////////////////////////////////

template <unsigned int DIM, typename TValue, typename TSpec, typename TValues>
inline void
appendValues(SumList<DIM, TValue, SkipSumList<TSpec> > & me,
             TValues const & vals)
{
    typedef SkiplistPath<TValue, Skipsumlist_<DIM, TSpec> > TPath;
    typedef SkiplistElement<TValue, Skipsumlist_<DIM, TSpec> > TElement;
    typedef Pair<TValue, TElement *> TPair;
    typedef SumList<DIM, TValue, MiniSumList< > > TMiniSumList;

    //find end of skip list
    TPath path;
    _skipsumlistFind(me.map, GoEnd(), 0, path);

    if (!appendValues(path.data_elements[0]->minilist, vals))
    {//not enough space in the current last element: append new
        //create new element
        unsigned char height = _skiplistCreateHeight(me, path);
        TElement & el = _skiplistAllocateElement(me.map, height); //see map_skiplist.
        valueConstruct(& (el.minilist));

        //append value into minilist of the new element
        appendValues(el.minilist, vals);

        //link new element el into skiplist, adjust link weights
        for (int i = height; i >= 0; --i)
        {
            path.data_elements[i]->data_next[i].data_element = & el;

            el.data_next[i].data_element = 0;
            el.data_next[i].values = vals;
        }
        for (int i = me.map.data_height; i > height ; --i)
        {
            path.data_elements[i]->data_next[i].values += vals;
        }

        //adjust map properties
        ++me.map.data_length;
    }
    else
    {
        //update skiplist path sums
        for (int i = me.map.data_height; i >= 0 ; --i)
        {
            path.data_elements[i]->data_next[i].values += vals;
        }
    }

    //update sumlist variables
    me.sum += vals;
    ++me.length;
}
template <unsigned int DIM, typename TValue, typename TSpec>
inline void
appendValues(SumList<DIM, TValue, SkipSumList<TSpec> > & me,
             TValue * p_vals)
{
    SumListValues<DIM, TValue> vals(p_vals);
    appendValues(me, vals);
}
template <unsigned int DIM, typename TValue, typename TSpec>
inline void
appendValues(SumList<DIM, TValue, SkipSumList<TSpec> > & me,
             TValue const * p_vals)
{
    SumListValues<DIM, TValue> vals(p_vals);
    appendValues(me, vals);
}


//////////////////////////////////////////////////////////////////////////////
// Iterator for SkipSumList
//////////////////////////////////////////////////////////////////////////////

//Helpers for SkipSumListIterator

template <typename TSumList>
struct SumlistSkipListElement_
{//dummy implementation
    typedef SkiplistElement<int, void> Type;
};
template <unsigned int DIM, typename TValue, typename TSpec>
struct SumlistSkipListElement_< SumList<DIM, TValue, SkipSumList<TSpec> > >
{
    typedef SkiplistElement<TValue, Skipsumlist_<DIM, TSpec> > Type;
};
template <unsigned int DIM, typename TValue, typename TSpec>
struct SumlistSkipListElement_< SumList<DIM, TValue, SkipSumList<TSpec> > const >
{
    typedef SkiplistElement<TValue, Skipsumlist_<DIM, TSpec> > const Type;
};

//////////////////////////////////////////////////////////////////////////////

/*
template <typename TSumList>
struct SumListSkiplistMinilist_
{//dummy implementation
    typedef SumList<1, int, MiniSumList< > > Type;
};

template <unsigned int DIM, typename TValue, typename TSpec>
struct SumListSkiplistMinilist_< SumList<DIM, TValue, SkipSumList<TSpec> > >
{
    typedef SumList<DIM, TValue, MiniSumList< > > Type;
};
template <unsigned int DIM, typename TValue, typename TSpec>
struct SumListSkiplistMinilist_< SumList<DIM, TValue, SkipSumList<TSpec> > const >
{
    typedef SumList<DIM, TValue, MiniSumList< > > const Type;
};
*/

//////////////////////////////////////////////////////////////////////////////

struct SkipSumListIterator;

template <typename TSumList>
class Iter<TSumList, SkipSumListIterator>
{
public:
    typedef typename Value<TSumList>::Type TValue;
    typedef SumList<DIMENSION<TSumList>::VALUE, TValue, MiniSumList< > > TMiniSumList_;
    typedef typename CopyConst_<TSumList, TMiniSumList_>::Type TMiniSumList;
//    typedef typename SumListSkiplistMinilist_<TSumList>::Type TMiniSumList;
    typedef typename Iterator<TMiniSumList>::Type TMiniSumListIterator;
    typedef typename SumlistSkipListElement_<TSumList>::Type TElement;

    TSumList * container;
    mutable TMiniSumListIterator iter;
    mutable TElement * element;

    Iter()
    {
    }
    Iter(TSumList & cont)
        : container(& cont)
    {
        goBegin(*this);
    }
    Iter(TSumList & cont, GoEnd)
        : container(& cont)
        , element(0)
    {
    }
    Iter(Iter const & other)
        : container(other.container)
        , iter(other.iter)
        , element(other.element)
    {
    }

    ~Iter()
    {
    }
    inline Iter const &
    operator = (Iter const & other)
    {
        container = other.container;
        iter = other.iter;
        element = other.element;
        return *this;
    }
};


//////////////////////////////////////////////////////////////////////////////

template <unsigned int DIM, typename TValue, typename TSpec, typename TIteratorSpec>
struct Iterator< SumList<DIM, TValue, SkipSumList<TSpec> >, TIteratorSpec>
{
    typedef SumList<DIM, TValue, SkipSumList<TSpec> > TSumList_;
    typedef Iter<TSumList_, SkipSumListIterator> Type;
};
template <unsigned int DIM, typename TValue, typename TSpec, typename TIteratorSpec>
struct Iterator< SumList<DIM, TValue, SkipSumList<TSpec> > const, TIteratorSpec>
{
    typedef SumList<DIM, TValue, SkipSumList<TSpec> > const TSumList_;
    typedef Iter<TSumList_, SkipSumListIterator> Type;
};

//////////////////////////////////////////////////////////////////////////////

template <typename TSumList>
inline void
goNext(Iter<TSumList, SkipSumListIterator> & it)
{
    goNext(it.iter);
    if (atEnd(it.iter))
    {
        it.element = it.element->data_next[0].data_element;
        if (it.element)
        {
            it.iter.container_ = & it.element->minilist;
            it.iter.next_ = it.iter.here_ = it.element->minilist.data_;
            scanValues(it.iter.next_, it.iter.values_);
        }
    }
}

//////////////////////////////////////////////////////////////////////////////

//nicht ganz sauber, aber sollte funzen
template <typename TSumList>
inline void
goPrevious(Iter<TSumList, SkipSumListIterator> & it)
{
    searchSumList(it, getSum(it, 0) - 1, 0);
}

//////////////////////////////////////////////////////////////////////////////

template <typename TSumList>
inline void
goBegin(Iter<TSumList, SkipSumListIterator> & it)
{
    it.element = & (it.container->map.data_border);
    it.iter = begin(it.element->minilist);
}

//////////////////////////////////////////////////////////////////////////////

template <typename TSumList>
inline void
goEnd(Iter<TSumList, SkipSumListIterator> & it)
{
    it.element = 0;
    it.iter.sums_ = it.container->sum;
}

//////////////////////////////////////////////////////////////////////////////

template <unsigned int DIM, typename TValue, typename TSpec>
inline void
goBeforeEnd(Iter<SumList<DIM, TValue, SkipSumList<TSpec> >, SkipSumListIterator > & it)
{
    typedef SkiplistPath<TValue, Skipsumlist_<DIM, TSpec> > TPath;
    typedef SkiplistElement<TValue, Skipsumlist_<DIM, TSpec> > TElement;
    typedef SumList<DIM, TValue, MiniSumList< > > TMiniSumList;

    //find end of skip list
    TPath path;
    _skipsumlistFind(it.container->map, GoEnd(), 0, path);

    it.element = path.data_elements[0];

    //set minilist iterator to end of last minilist
    it.iter.container_ = & it.element->minilist;
    goBeforeEnd(it.iter);
}

//////////////////////////////////////////////////////////////////////////////

template <typename TSumList>
inline bool
atBegin(Iter<TSumList, SkipSumListIterator> & it)
{
    return (it.iter == begin(it.element->minilist));
}

//////////////////////////////////////////////////////////////////////////////

template <typename TSumList>
inline bool
atEnd(Iter<TSumList, SkipSumListIterator> & it)
{
    return (!it.element);
}

//////////////////////////////////////////////////////////////////////////////

template <typename TSumList>
inline typename Value<TSumList>::Type
getValue(Iter<TSumList, SkipSumListIterator > & it,
         int dim)
{
    return it.iter.values_[dim];
}
template <typename TSumList>
inline typename Value<TSumList>::Type
getValue(Iter<TSumList, SkipSumListIterator > const & it,
         int dim)
{
    return it.iter.values_[dim];
}

//////////////////////////////////////////////////////////////////////////////

template <typename TSumList>
inline typename Values<TSumList>::Type
getValues(Iter<TSumList, SkipSumListIterator > & it)
{
    return it.iter.values_;
}
template <typename TSumList>
inline typename Values<TSumList>::Type
getValues(Iter<TSumList, SkipSumListIterator > const & it)
{
    return it.iter.values_;
}

//////////////////////////////////////////////////////////////////////////////

template <typename TSumList>
inline typename Value<TSumList>::Type
getSum(Iter<TSumList, SkipSumListIterator > & it,
       int dim)
{
    return it.iter.sums_[dim];
}
template <typename TSumList>
inline typename Value<TSumList>::Type
getSum(Iter<TSumList, SkipSumListIterator > const & it,
       int dim)
{
    return it.iter.sums_[dim];
}

//////////////////////////////////////////////////////////////////////////////
//split the mini list at it and insert a new mini list in skiplist
//path must be a path to it.element
//path and it are adjusted to the new element if necessary
template <unsigned int DIM, typename TValue, typename TSpec, typename TPath>
inline void
_splitMiniList(Iter<SumList<DIM, TValue, SkipSumList<TSpec> >, SkipSumListIterator > & it,
               TPath & path)
{
    typedef SkiplistElement<TValue, Skipsumlist_<DIM, TSpec> > TElement;
    typedef SumListValues<DIM, TValue> TValues;

    //height for new element, adjust path if necessary (if the height of the complete skiplist increases)
    unsigned char height = _skiplistCreateHeight(* it.container, path);

    //the current element
    TElement & el1 = *(it.element);

    //create new element
    TElement & el2 = _skiplistAllocateElement(it.container->map, height); //see map_skiplist.h
    valueConstruct(& (el2.minilist));

    //move half of the values in el1 to el2
    splitSumList(el1.minilist, el2.minilist);

    if (!length(el2))
    {//el1 was not full enough to create a split: remove el2 from skiplist and exit
     //(this should not happen, since _splitMiniList is only called it el1 is full)
        _skiplistDeallocateElement(it.container->map, el2, height);
        return;
    }

    //must adjust the iterator "it"?
    int it_bytepos = it.iter.here_.data_ptr - el1.minilist.data_;
    bool it_was_moved = (it_bytepos >= el1.minilist.data_size); //true if "it" points to el2

    //insert new element el2 into skiplist
    for (int i = height; i >= 0; --i)
    {
        TElement * predecessor;
        if (path.data_elements[i] == & el1)
        {//el2 is linked behind el1
            predecessor = & el1;

            //adjust values
            el2.data_next[i].values = predecessor->data_next[i].values;
            el2.data_next[i].values -= predecessor->minilist.data_sum;
            predecessor->data_next[i].values = predecessor->minilist.data_sum;

            if (it_was_moved)
            {//adjust path
                path.data_elements[i] = & el2;
                path.sums[i] += predecessor->minilist.data_sum;
            }
        }
        else
        {//el2 is linked behind another element
            predecessor = path.data_elements[i];

            //adjust values
            TValues pred_values = path.sums[0];
            pred_values -= path.sums[i];
            pred_values += el1.minilist.data_sum;

            el2.data_next[i].values = predecessor->data_next[i].values;
            el2.data_next[i].values -= pred_values;
            predecessor->data_next[i].values = pred_values;

            if (it_was_moved)
            {//adjust path
                path.data_elements[i] = & el2;
                path.sums[i] += pred_values;
            }
        }

        //link el2 into skiplist
        el2.data_next[i].data_element = predecessor->data_next[i].data_element;
        predecessor->data_next[i].data_element = & el2;
    }

    //adjust map properties
    ++it.container->map.data_length;

    //adjust it.iter if it was moved
    if (it_was_moved)
    {
        int here_length = it.iter.next_.data_ptr - it.iter.here_.data_ptr; //difference between next_ and here_
        it.iter.here_.data_ptr = el2.minilist.data_ + it_bytepos - el1.minilist.data_size; //redirect here_
        it.iter.next_.data_ptr = it.iter.here_.data_ptr + here_length; //redirect next_
        it.element = & el2;
        it.iter.container_ = & el2.minilist;
    }
}

//////////////////////////////////////////////////////////////////////////////

template <unsigned int DIM, typename TValue, typename TSpec, typename TValue2>
inline void
assignValue(Iter<SumList<DIM, TValue, SkipSumList<TSpec> >, SkipSumListIterator > & it,
            int dim,
            TValue2 val)
{
    typedef SkiplistPath<TValue, Skipsumlist_<DIM, TSpec> > TPath;
    typedef SkiplistElement<TValue, Skipsumlist_<DIM, TSpec> > TElement;
    typedef Pair<TValue, TElement *> TPair;

    TValue val_old = getValue(it, dim);

    TPath path;
    _skipsumlistFind(it.container->map, TPair(getSum(it, 0), it.element), 0, path);

    if (!assignValue(it.iter, dim, val))
    {//split
        _splitMiniList(it, path);
        assignValue(it.iter, dim, val);
    }

    //update skiplist path sums
    for (int i = it.container->map.data_height; i >= 0 ; --i)
    {
        if (path.data_elements[i]->data_next[i].data_element == it.element)
        {// adjust link outgoing from it.element
            it.element->data_next[i].values[dim] += (val - val_old);
        }
        else
        {// adjust link that jumps over it.element
            path.data_elements[i]->data_next[i].values[dim] += (val - val_old);
        }
    }
    it.container->sum[dim] += (val - val_old);
}

//////////////////////////////////////////////////////////////////////////////

template <unsigned int DIM, typename TValue, typename TSpec, typename TValues>
inline void
insertValues(Iter<SumList<DIM, TValue, SkipSumList<TSpec> >, SkipSumListIterator > & it,
             TValues const & vals)
{
    typedef SkiplistPath<TValue, Skipsumlist_<DIM, TSpec> > TPath;
    typedef SkiplistElement<TValue, Skipsumlist_<DIM, TSpec> > TElement;
    typedef Pair<TValue, TElement *> TPair;

    TPath path;
    _skipsumlistFind(it.container->map, TPair(getSum(it, 0), it.element), 0, path);

    if (!insertValues(it.iter, vals))
    {//split
        _splitMiniList(it, path);
        insertValues(it.iter, vals);
    }

    //update skiplist path sums
    for (int i = it.container->map.data_height; i >= 0 ; --i)
    {
        path.data_elements[i]->data_next[i].values += vals;
    }

    //update sumlist variables
    it.container->sum += vals;
    ++it.container->length;
}
template <unsigned int DIM, typename TValue, typename TSpec>
inline void
insertValues(Iter<SumList<DIM, TValue, SkipSumList<TSpec> >, SkipSumListIterator > & it,
             TValue * p_vals)
{
    SumListValues<DIM, TValue > vals(p_vals);
    insertValues(it, vals);
}
template <unsigned int DIM, typename TValue, typename TSpec>
inline void
insertValues(Iter<SumList<DIM, TValue, SkipSumList<TSpec> >, SkipSumListIterator > & it,
             TValue const * p_vals)
{
    SumListValues<DIM, TValue > vals(p_vals);
    insertValues(it, vals);
}

//////////////////////////////////////////////////////////////////////////////

template <unsigned int DIM, typename TValue, typename TSpec>
inline void
removeValues(Iter<SumList<DIM, TValue, SkipSumList<TSpec> >, SkipSumListIterator > & it)
{
    typedef SkiplistPath<TValue, Skipsumlist_<DIM, TSpec> > TPath;
    typedef SkiplistElement<TValue, Skipsumlist_<DIM, TSpec> > TElement;
    typedef Pair<TValue, TElement *> TPair;

    TPath path;
    _skipsumlistFind(it.container->map, TPair(getSum(it, 0), it.element), 0, path);

    //update skiplist path sums
    for (int i = it.container->map.data_height; i >= 0 ; --i)
    {
        path.data_elements[i]->data_next[i].values -= it.iter.values_;
    }

    //update sumlist variables
    it.container->sum -= it.iter.values_;
    --it.container->length;

    //remove values
    removeValues(it.iter);

    if (atEnd(it.iter))
    {//move to next element
        it.element = it.element->data_next[0].data_element;
        if (it.element)
        {
            it.iter.container_ = & it.element->minilist;
            it.iter.here_ = it.iter.next_ = it.iter.container_->data_;
            scanValues(it.iter.next_, it.iter.values_);
        }
    }
}

//////////////////////////////////////////////////////////////////////////////

//template <unsigned int DIM, typename TValue, typename TSpec, typename TValue2>
//inline void
//searchSumList(Iter<SumList<DIM, TValue, SkipSumList<TSpec> >, SkipSumListIterator > & it,
//              TValue2 const & val,
//              int dim)
template <unsigned int DIM, typename TValue, typename TSpec, typename TSumList, typename TValue2>
inline void
_searchSumlistSkip(Iter<TSumList, SkipSumListIterator > & it,
                    TValue2 const & val,
                    int dim)
{
    typedef SkiplistPath<TValue, Skipsumlist_<DIM, TSpec> > TPath;
//    typedef SkiplistElement<TValue, Skipsumlist_<DIM, TSpec> > TElement;
//    typedef Pair<TValue, TElement *> TPair;
//    typedef SumList<DIM, TValue, MiniSumList< > > TMiniSumList;

    TPath path;

  // The assumpiton here is that val is signed and getSum() fits into it.
    if (val >= static_cast<TValue2>(getSum(*it.container, dim)))
    {
        goEnd(it);
        return;
    }

    _skipsumlistFind(it.container->map, val, dim, path);

    it.element = path.data_elements[0];
    it.iter.container_ = & it.element->minilist;
    searchSumList(it.iter, val - path.sums[0][dim], dim);
    it.iter.sums_ += path.sums[0];
}

template <unsigned int DIM, typename TValue, typename TSpec, typename TValue2>
inline void
searchSumList(Iter<SumList<DIM, TValue, SkipSumList<TSpec> >, SkipSumListIterator > & it,
              TValue2 const & val,
              int dim)
{
    _searchSumlistSkip<DIM, TValue, TSpec>(it, val, dim);
}
template <unsigned int DIM, typename TValue, typename TSpec, typename TValue2>
inline void
searchSumList(Iter<SumList<DIM, TValue, SkipSumList<TSpec> > const, SkipSumListIterator > & it,
              TValue2 const & val,
              int dim)
{
    _searchSumlistSkip<DIM, TValue, TSpec>(it, val, dim);
}
//template <unsigned int DIM, typename TValue, typename TSpec, typename TValue2>
//inline void
//searchSumList(Iter<SumList<DIM, TValue, SkipSumList<TSpec> > const, SkipSumListIterator > & it,
//              TValue2 const & val,
//              int dim)
//{
//    _searchSumlistSkip<DIM, TValue, TSpec>(it, val, dim);
//}

//////////////////////////////////////////////////////////////////////////////

template <typename TSumList, typename TSumList2>
inline bool
operator == (Iter<TSumList, SkipSumListIterator> const & left,
             Iter<TSumList2, SkipSumListIterator> const & right)
{
    return ((left.element == 0) && (right.element == 0)) || ((left.element == right.element) && (left.iter == right.iter));
}

template <typename TSumList, typename TSumList2>
inline bool
operator != (Iter<TSumList, SkipSumListIterator> const & left,
             Iter<TSumList2, SkipSumListIterator> const & right)
{
    return !(left == right);
}


//////////////////////////////////////////////////////////////////////////////



//____________________________________________________________________________


//////////////////////////////////////////////////////////////////////////////

}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
