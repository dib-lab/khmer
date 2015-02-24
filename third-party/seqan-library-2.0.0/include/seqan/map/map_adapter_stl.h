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

#ifndef SEQAN_HEADER_MAP_ADAPTER_STL_H
#define SEQAN_HEADER_MAP_ADAPTER_STL_H

//////////////////////////////////////////////////////////////////////////////

namespace SEQAN_NAMESPACE_MAIN
{

template <typename TKey, typename TCargo, typename TCompare, typename TAlloc>
struct Value< std::map<TKey, TCargo, TCompare, TAlloc> >
{
    typedef Pair<TKey,TCargo> Type;
};

template <typename TKey, typename TCompare, typename TAlloc>
struct Value< std::set<TKey, TCompare, TAlloc> >
{
    typedef bool Type;
};

//////////////////////////////////////////////////////////////////////////////

template <typename TKey, typename TCargo, typename TCompare, typename TAlloc>
struct Key< std::map<TKey,TCargo, TCompare, TAlloc>  >
{
    typedef TKey Type;
};

template <typename TKey, typename TCompare, typename TAlloc>
struct Key< std::set<TKey, TCompare, TAlloc>  >
{
    typedef TKey Type;
};

//////////////////////////////////////////////////////////////////////////////

template <typename TKey, typename TCargo, typename TCompare, typename TAlloc>
struct Cargo< std::map<TKey,TCargo, TCompare, TAlloc> >
{
    typedef TCargo Type;
};

//////////////////////////////////////////////////////////////////////////////

template <typename TKey, typename TCompare, typename TAlloc>
struct Size< std::set<TKey, TCompare, TAlloc> >
{
    //typedef std::set<TKey>::size_type Type;
    typedef unsigned Type;
};

template <typename TKey, typename TCargo, typename TCompare, typename TAlloc>
struct Size< std::map<TKey,TCargo, TCompare, TAlloc> >
{
    //typedef std::map<TKey,TCargo>::size_type Type;
    typedef unsigned Type;
};

//////////////////////////////////////////////////////////////////////////////

// traits for stl usage
/*
template <typename TKey, typename TCompare, typename TAlloc>
struct AllocatorType< std::set<TKey, TCompare, TAlloc> >
{
    typedef TAlloc Type;
};

template <typename TKey, typename TCargo, typename TCompare, typename TAlloc>
struct AllocatorType< std::map<TKey, TCargo, TCompare, TAlloc> >
{
    typedef TAlloc Type;
};

template <typename T>
struct StlComparator_;

template <typename TKey, typename TCompare, typename TAlloc>
struct StlComparator_< std::set<TKey, TCompare, TAlloc> >
{
    typedef TCompare Type;
};

template <typename TKey, typename TCargo, typename TCompare, typename TAlloc>
struct StlComparator_< std::map<TKey, TCargo, TCompare, TAlloc> >
{
    typedef TCompare Type;
};

*/

// TODO(rmaerker): Default type should not be int.
// TODO(rmaerker): Metafunction for const iterator is missing.

template <typename T>
struct StlIterator_
{
    typedef int Type;
};

template <typename TKey, typename TCompare, typename TAlloc>
struct StlIterator_< std::set<TKey, TCompare, TAlloc> >
{
    typedef typename std::set<TKey, TCompare, TAlloc>::iterator Type;
};

template <typename TKey, typename TCargo, typename TCompare, typename TAlloc>
struct StlIterator_< std::map<TKey, TCargo, TCompare, TAlloc> >
{
    typedef typename std::map<TKey, TCargo, TCompare, TAlloc>::iterator Type;
};

//////////////////////////////////////////////////////////////////////////////

template <typename TKey, typename TCompare, typename TAlloc>
inline void
assign(std::set<TKey, TCompare, TAlloc> & target,
       std::set<TKey, TCompare, TAlloc> const & source)
{
    target = source;
}

template <typename TKey,typename TCargo, typename TCompare, typename TAlloc>
inline void
assign(std::map<TKey,TCargo, TCompare, TAlloc> & target,
       std::map<TKey,TCargo, TCompare, TAlloc> const & source)
{
    target = source;
}

//////////////////////////////////////////////////////////////////////////////

template <typename TValue, typename TCompare, typename TAlloc, typename TKey2>
inline bool
hasKey(std::set<TValue,TCompare,TAlloc> & me, TKey2 const & _key)
{
SEQAN_CHECKPOINT
    return (me.count(_key) != 0);
}

template <typename TValue, typename TCompare, typename TAlloc, typename TKey2>
inline bool
hasKey(std::set<TValue,TCompare,TAlloc> const & me, TKey2 const & _key)
{
SEQAN_CHECKPOINT
    return (me.count(_key) != 0);
}


template <typename TKey, typename TCargo, typename TCompare, typename TAlloc, typename TKey2>
inline bool
hasKey(std::map<TKey, TCargo, TCompare, TAlloc> & me, TKey2 const & _key)
{
SEQAN_CHECKPOINT
    return (me.count(_key) != 0);
}

template <typename TKey, typename TCargo, typename TCompare, typename TAlloc, typename TKey2>
inline bool
hasKey(std::map<TKey, TCargo, TCompare, TAlloc> const & me, TKey2 const & _key)
{
SEQAN_CHECKPOINT
    return (me.count(_key) != 0);
}

//////////////////////////////////////////////////////////////////////////////

template <typename TValue, typename TCompare, typename TAlloc>
inline typename Size< std::set<TValue, TCompare, TAlloc> >::Type
length(std::set<TValue, TCompare, TAlloc> const & me)
{
SEQAN_CHECKPOINT
    return me.size();
}

template <typename TKey, typename TCargo, typename TCompare, typename TAlloc>
inline typename Size< std::map<TKey,TCargo, TCompare, TAlloc> >::Type
length(std::map<TKey,TCargo, TCompare, TAlloc> const & me)
{
SEQAN_CHECKPOINT
    return me.size();
}

//////////////////////////////////////////////////////////////////////////////

template <typename TValue, typename TCompare, typename TAlloc, typename TValue2>
inline void
insert(std::set<TValue, TCompare, TAlloc> & me,TValue2 const & _value)
{
SEQAN_CHECKPOINT
    me.insert(_value);
}

template <typename TKey, typename TCargo, typename TCompare, typename TAlloc, typename TKey2, typename TCargo2>
inline void
insert(std::map<TKey,TCargo, TCompare, TAlloc> & me,
       TKey2 const & _key,
       TCargo2 const & _cargo)
{
SEQAN_CHECKPOINT

    me[_key] = _cargo;
}
template <typename TKey, typename TCargo, typename TCompare, typename TAlloc, typename TKey2>
inline void
insert(std::map<TKey,TCargo, TCompare, TAlloc> & me,TKey2 const & _key)
{
SEQAN_CHECKPOINT

    insert(me, _key, TCargo());
}
template <typename TKey, typename TCargo, typename TCompare, typename TAlloc, typename TKey2, typename TCargo2, typename TSpec>
inline void
insert(std::map<TKey,TCargo, TCompare, TAlloc> & me, Pair<TKey2,TCargo2,TSpec> const & _value)
{
SEQAN_CHECKPOINT
    insert(me, _value.i1, _value.i2);
}


//////////////////////////////////////////////////////////////////////////////

template <typename TValue, typename TCompare, typename TAlloc>
inline void
clear(std::set<TValue, TCompare, TAlloc> & me)
{
SEQAN_CHECKPOINT
    me.clear();
}

template <typename TKey, typename TCargo, typename TCompare, typename TAlloc>
inline void
clear(std::map<TKey,TCargo, TCompare, TAlloc> & me)
{
SEQAN_CHECKPOINT
    me.clear();
}

//////////////////////////////////////////////////////////////////////////////

//template <typename TKey, typename TCargo, typename TCompare, typename TAlloc, typename TKey2>
//inline typename Value< std::map<TKey,TCargo, TCompare, TAlloc> >::Type &
//value(std::map<TKey,TCargo, TCompare, TAlloc> & me,
//      TKey2 const & _key)
//{
//    return me[_key];
//}

//////////////////////////////////////////////////////////////////////////////

template <typename TKey, typename TCargo, typename TCompare, typename TAlloc, typename TKey2>
inline typename Cargo< std::map<TKey,TCargo, TCompare, TAlloc> >::Type &
cargo(std::map<TKey,TCargo, TCompare, TAlloc> & me,
      TKey2 const & _key)
{
SEQAN_CHECKPOINT
    return me[_key];
}

//////////////////////////////////////////////////////////////////////////////

struct StlSetIterator;
struct StlMapIterator;

//////////////////////////////////////////////////////////////////////////////
//                            MapIterator                                   //
//////////////////////////////////////////////////////////////////////////////

template <typename TKey, typename TCargo, typename TCompare, typename TAlloc, typename TIteratorSpec>
struct Iterator< std::map<TKey,TCargo, TCompare, TAlloc> , TIteratorSpec >
{
    typedef std::map<TKey,TCargo, TCompare, TAlloc> TStlMap_;
    typedef Iter<TStlMap_, StlMapIterator> Type;
};


template <typename TStlMap>
class Iter<TStlMap, StlMapIterator>
{
public:
    //typedef typename std::map<typename Key<TSTLMap>::Type,
    //                            typename Cargo<TSTLMap>::Type,
    //                            typename StlComparator_<TStlMap>::Type,
    //                            typename AllocatorType<TStlMap>::Type >::iterator THostIter;

    typedef typename StlIterator_<TStlMap>::Type THostIter;
    THostIter _iter;
    Holder<TStlMap> host_map_holder;

    Iter()
    {
SEQAN_CHECKPOINT
    }

    Iter(Iter const & other)
        : _iter(other._iter)
    {
SEQAN_CHECKPOINT
        host_map_holder = other.host_map_holder;
    }

    Iter(TStlMap & map)
        : _iter(map.begin())
    {
SEQAN_CHECKPOINT
        setValue(host_map_holder,map);
    }

    ~Iter()
    {
SEQAN_CHECKPOINT
    }

    Iter const &
    operator = (Iter const & other)
    {
SEQAN_CHECKPOINT
        host_map_holder = other.host_map_holder;
        _iter = other._iter;
        return *this;
    }
    operator bool () const
    {
SEQAN_CHECKPOINT
        return (_iter != value(host_map_holder).end());
    }

};

//////////////////////////////////////////////////////////////////////////////

template <typename TStlMap>
inline bool
operator == (Iter<TStlMap, StlMapIterator> const & left,
             Iter<TStlMap, StlMapIterator> const & right)
{
SEQAN_CHECKPOINT
    return left._iter == right._iter;
}

template <typename TStlMap>
inline bool
operator != (Iter<TStlMap, StlMapIterator> const & left,
             Iter<TStlMap, StlMapIterator> const & right)
{
SEQAN_CHECKPOINT
    return left._iter != right._iter;
}

//////////////////////////////////////////////////////////////////////////////

template <typename TKey, typename TCargo, typename TCompare, typename TAlloc, typename TIteratorSpec>
inline typename Iterator< std::map<TKey,TCargo, TCompare, TAlloc>, TIteratorSpec>::Type
begin(std::map<TKey,TCargo, TCompare, TAlloc> & me,
      TIteratorSpec)
{
    typedef std::map<TKey,TCargo, TCompare, TAlloc> TStlMap;
    typedef typename Iterator<TStlMap , TIteratorSpec>::Type TIterator;
    return TIterator(me);
}
template <typename TKey, typename TCargo, typename TCompare, typename TAlloc>
inline typename Iterator< std::map<TKey,TCargo, TCompare, TAlloc> >::Type
begin(std::map<TKey,TCargo, TCompare, TAlloc> & me)
{
    typedef std::map<TKey,TCargo, TCompare, TAlloc> TStlMap;
    typedef typename Iterator<TStlMap>::Type TIterator;
    return TIterator(me);
}

//////////////////////////////////////////////////////////////////////////////

template <typename TKey, typename TCargo, typename TCompare, typename TAlloc, typename TIteratorSpec>
inline typename Iterator< std::map<TKey,TCargo, TCompare, TAlloc> >::Type
end(std::map<TKey,TCargo, TCompare, TAlloc> & me,
      TIteratorSpec)
{
    typedef std::map<TKey,TCargo, TCompare, TAlloc> TStlMap;
    typedef typename Iterator<TStlMap, TIteratorSpec>::Type TIterator;
    TIterator _iter(me);
    _iter._iter = me.end();
    return _iter;
}

template <typename TKey, typename TCargo, typename TCompare, typename TAlloc>
inline typename Iterator< std::map<TKey,TCargo, TCompare, TAlloc> >::Type
end(std::map<TKey,TCargo, TCompare, TAlloc> & me)
{
    typedef std::map<TKey,TCargo, TCompare, TAlloc> TStlMap;
    typedef typename Iterator<TStlMap>::Type TIterator;
    TIterator _iter(me);
    _iter._iter = me.end();
    return _iter;
}

//////////////////////////////////////////////////////////////////////////////

template <typename TStlMap>
inline bool
atEnd(Iter<TStlMap, StlMapIterator> & it)
{
SEQAN_CHECKPOINT
    return (it._iter == value(it.host_map_holder).end());
}

//////////////////////////////////////////////////////////////////////////////

template <typename TStlMap>
inline void
goNext(Iter<TStlMap, StlMapIterator> & it)
{
SEQAN_CHECKPOINT
    it._iter++;
}

//////////////////////////////////////////////////////////////////////////////

template <typename TStlMap>
inline typename Value<TStlMap>::Type &
value(Iter<TStlMap, StlMapIterator> & it)
{
SEQAN_CHECKPOINT
    return it._iter->second;
}
template <typename TStlMap>
inline typename Value<TStlMap>::Type &
value(Iter<TStlMap, StlMapIterator> const & it)
{
SEQAN_CHECKPOINT
    return it._iter->second;
}

//////////////////////////////////////////////////////////////////////////////

template <typename TStlMap>
inline typename Key<TStlMap>::Type const &
key(Iter<TStlMap, StlMapIterator> & it)
{
SEQAN_CHECKPOINT
    return it._iter->first;
}
template <typename TStlMap>
inline typename Key<TStlMap>::Type const &
key(Iter<TStlMap, StlMapIterator> const & it)
{
SEQAN_CHECKPOINT
    return it._iter->first;
}

//////////////////////////////////////////////////////////////////////////////

template <typename TStlMap>
inline typename Cargo<TStlMap>::Type &
cargo(Iter<TStlMap, StlMapIterator> & it)
{
SEQAN_CHECKPOINT
    return it._iter->second;
}
template <typename TStlMap>
inline typename Cargo<TStlMap>::Type &
cargo(Iter<TStlMap, StlMapIterator> const & it)
{
SEQAN_CHECKPOINT
    return it._iter->second;
}

////////////////////////////////////////////////////////////////////////////////
//                             SetIterator                                    //
////////////////////////////////////////////////////////////////////////////////

template <typename TValue, typename TCompare, typename TAlloc, typename TIteratorSpec>
struct Iterator< std::set<TValue, TCompare, TAlloc> , TIteratorSpec >
{
    typedef std::set<TValue, TCompare, TAlloc> TStlSet_;
    typedef Iter<TStlSet_, StlSetIterator> Type;
};

template <typename TStlMap>
class Iter< TStlMap, StlSetIterator>
{
public:
    //typedef typename std::set<typename Key<TSTLMap>::Type,
    //                            typename StlComparator_<TSTLMap>::Type,
    //                            typename AllocatorType<TSTLMap>::Type >::iterator THostIter;
    typedef typename StlIterator_<TStlMap>::Type THostIter;
    THostIter _iter;
    Holder<TStlMap> host_map_holder;

    Iter()
    {
SEQAN_CHECKPOINT
    }

    Iter(Iter const & other)
        : _iter(other._iter)
    {
SEQAN_CHECKPOINT
        host_map_holder = other.host_map_holder;
    }

    Iter(TStlMap & map)
        : _iter(map.begin())
    {
SEQAN_CHECKPOINT
        setValue(host_map_holder,map);
    }

    ~Iter()
    {
SEQAN_CHECKPOINT
    }

    Iter const &
    operator = (Iter const & other)
    {
SEQAN_CHECKPOINT
        host_map_holder = other.host_map_holder;
        _iter = other._iter;
        return *this;
    }
    operator bool () const
    {
SEQAN_CHECKPOINT
        return (_iter != value(host_map_holder).end());
    }

};

//////////////////////////////////////////////////////////////////////////////

template <typename TStlMap>
inline bool
operator == (Iter<TStlMap, StlSetIterator> const & left,
             Iter<TStlMap, StlSetIterator> const & right)
{
SEQAN_CHECKPOINT
    return left._iter == right._iter;
}

template <typename TStlMap>
inline bool
operator != (Iter<TStlMap, StlSetIterator> const & left,
             Iter<TStlMap, StlSetIterator> const & right)
{
SEQAN_CHECKPOINT
    return left._iter != right._iter;
}

//////////////////////////////////////////////////////////////////////////////

template <typename TValue, typename TCompare, typename TAlloc, typename TIteratorSpec>
inline typename Iterator< std::set<TValue,TCompare,TAlloc>, TIteratorSpec>::Type
begin(std::set<TValue,TCompare,TAlloc> & me,
      TIteratorSpec)
{
    typedef std::set<TValue,TCompare,TAlloc> TStlMap;
    typedef typename Iterator<TStlMap , TIteratorSpec>::Type TIterator;
    return TIterator(me);
}
template <typename TValue, typename TCompare, typename TAlloc>
inline typename Iterator< std::set<TValue,TCompare,TAlloc> >::Type
begin(std::set<TValue,TCompare,TAlloc> & me)
{
    typedef std::set<TValue,TCompare,TAlloc> TStlMap;
    typedef typename Iterator<TStlMap>::Type TIterator;
    return TIterator(me);
}

//////////////////////////////////////////////////////////////////////////////

template <typename TValue, typename TCompare, typename TAlloc, typename TIteratorSpec>
inline typename Iterator< std::set<TValue, TCompare, TAlloc>, TIteratorSpec>::Type
end(std::set<TValue, TCompare, TAlloc> & me,
      TIteratorSpec)
{
    typedef std::set<TValue,TCompare,TAlloc> TStlMap;
    typedef typename Iterator<TStlMap, TIteratorSpec>::Type TIterator;
    TIterator _iter(me);
    _iter._iter = me.end();
    return _iter;
}

template <typename TValue, typename TCompare, typename TAlloc>
inline typename Iterator< std::set<TValue,TCompare,TAlloc> >::Type
end(std::set<TValue,TCompare,TAlloc> & me)
{
    typedef std::set<TValue,TCompare,TAlloc> TStlMap;
    typedef typename Iterator<TStlMap>::Type TIterator;
    TIterator _iter(me);
    _iter._iter = me.end();
    return _iter;
}

//////////////////////////////////////////////////////////////////////////////

template <typename TStlMap>
inline bool
atEnd(Iter<TStlMap, StlSetIterator> & it)
{
SEQAN_CHECKPOINT
    return (it._iter == value(it.host_map_holder).end());
}

//////////////////////////////////////////////////////////////////////////////

template <typename TStlMap>
inline void
goNext(Iter<TStlMap, StlSetIterator> & it)
{
SEQAN_CHECKPOINT
    it._iter++;
}

//////////////////////////////////////////////////////////////////////////////
//
//template <typename TStlMap>
//inline typename Value<TStlMap>::Type &
//value(Iter<TStlMap, StlSetIterator> & it)
//{
//    return hasKey(*it);
//}
//template <typename TStlMap>
//inline typename Value<TStlMap>::Type &
//value(Iter<TStlMap, StlSetIterator> const & it)
//{
//    return hasKey(*it);
//}

//////////////////////////////////////////////////////////////////////////////

template <typename TStlMap>
inline typename Key<TStlMap>::Type const &
key(Iter<TStlMap, StlSetIterator> & it)
{
SEQAN_CHECKPOINT
    return (*it._iter);
}
template <typename TStlMap>
inline typename Key<TStlMap>::Type const &
key(Iter<TStlMap, StlSetIterator> const & it)
{
SEQAN_CHECKPOINT
    return (*it._iter);
}

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////

template <typename TValue, typename TCompare, typename TAlloc,typename TFind>
inline typename Iterator< std::set<TValue, TCompare, TAlloc> >::Type
find(std::set<TValue, TCompare, TAlloc> & me,
     TFind const & _find)
{
SEQAN_CHECKPOINT
    typedef std::set<TValue, TCompare, TAlloc> TMap;
    typedef typename Iterator< TMap >::Type TMapIterator;

    TMapIterator _iter(me);
    _iter._iter = me.find(_find);
    if(!_iter){
        _iter._iter = me.upper_bound(_find);
    }
    return _iter;
}

template <typename TKey, typename TCargo, typename TCompare, typename TAlloc, typename TFind>
inline typename Iterator< std::map<TKey,TCargo, TCompare, TAlloc> >::Type
find(std::map<TKey,TCargo, TCompare, TAlloc> & me,
     TFind const & _find)
{
SEQAN_CHECKPOINT
    typedef std::map<TKey,TCargo, TCompare, TAlloc> TMap;
    typedef typename Iterator< TMap >::Type TMapIterator;

    TMapIterator _iter(me);
    _iter._iter = me.find(_find);
    if(!_iter){
        _iter._iter = me.upper_bound(_find);
    }
    return _iter;
}


//////////////////////////////////////////////////////////////////////////////

template <typename TValue, typename TCompare, typename TAlloc, typename TMap2>
inline void
erase(std::set<TValue, TCompare, TAlloc> & me,
      Iter<TMap2, StlSetIterator> const & it)
{
SEQAN_CHECKPOINT
    me.erase(it._iter);
}

template <typename TKey, typename TCargo ,typename TMap2>
inline void
erase(std::map<TKey,TCargo> & me,
      Iter<TMap2, StlMapIterator> const & it)
{
SEQAN_CHECKPOINT
    me.erase(it._iter);
}


//////////////////////////////////////////////////////////////////////////////

template <typename TValue, typename TCompare, typename TAlloc, typename TToRemove>
inline void
erase(std::set<TValue, TCompare, TAlloc> & me,
      TToRemove const & to_remove)
{
SEQAN_CHECKPOINT
    me.erase(to_remove);
}

template <typename TKey, typename TCargo, typename TCompare, typename TAlloc, typename TToRemove>
inline void
erase(std::map<TKey,TCargo, TCompare, TAlloc> & me,
      TToRemove const & to_remove)
{
SEQAN_CHECKPOINT
    me.erase(to_remove);
}

}
#endif // #ifndef SEQAN_HEADER
