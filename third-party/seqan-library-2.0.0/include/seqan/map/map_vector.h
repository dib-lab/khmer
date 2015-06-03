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


#ifndef SEQAN_HEADER_MAP_VECTOR_H
#define SEQAN_HEADER_MAP_VECTOR_H

//////////////////////////////////////////////////////////////////////////////

namespace SEQAN_NAMESPACE_MAIN
{
/*!
 * @class VectorSet
 * @extends Map
 * @headerfile <seqan/map.h>
 * @brief A fast map for small key types.
 *
 * @signature template <typename TValue[, typename TSpec]>
 *            class Map<TValue, VectorSet<TSpec> >;
 *
 * @tparam TValue The type of value stored in the map.
 * @tparam TSpec  The specializing type.  <tt>TSpec</tt> is used as specialization for the String type that is
 *                used to store the values, defaults to <tt>Alloc&lt;&gt;</tt>.
 *
 * The memory needed is linear to the number different values the Key type of <tt>TValue</tt> can
 * get.  So do not use this map for large key types.
 */

//////////////////////////////////////////////////////////////////////////////
//
//    VectorSet
//
//////////////////////////////////////////////////////////////////////////////


template <typename TSpec = Alloc<> >
struct VectorSet;


//////////////////////////////////////////////////////////////////////////////
// internal set meta-functions

template <typename TCargo>
struct VectorSetElement_
{
    bool data_valid;
    TCargo data_cargo;
};
template <>
struct VectorSetElement_<Nothing>
{
    bool data_valid;
};

//////////////////////////////////////////////////////////////////////////////

template <typename T>
struct VectorSetElements_
{
    typedef char * Type; //dummy implementation for VC++
};
template <typename TValue, typename TSpec>
struct VectorSetElements_<Map<TValue, VectorSet<TSpec> > >
{
    typedef Map<TValue, VectorSet<TSpec> > TMap_;
    typedef typename Cargo<TMap_>::Type TCargo_;
    typedef VectorSetElement_<TCargo_> TElement_;
    typedef String<TElement_, TSpec> Type;
};

//////////////////////////////////////////////////////////////////////////////
// VectorSet class

template <typename TValue, typename TSpec>
class Map<TValue, VectorSet<TSpec> >
{
public:
    typedef typename Key<Map>::Type TKey;
    typedef typename Size<Map>::Type TSize;
    typedef typename VectorSetElements_<Map>::Type TElements;

    TElements    data_elements;
    TSize        data_length;

    Map()
    {
        resize(data_elements, (TSize) ValueSize<TKey>::VALUE);
        clear(*this);
    }

    Map(TSize _vectorSize)
    {
        resize(data_elements, (TSize) _vectorSize);
        clear(*this);
    }

    template <typename TKey>
    inline typename MapValue<Map>::Type
    operator [] (TKey const & key)
    {
        return mapValue(*this, key);
    }

};


//////////////////////////////////////////////////////////////////////////////
// set meta-functions

template <typename TValue, typename TSpec>
struct Value< Map<TValue, VectorSet<TSpec> > >
{
    typedef TValue Type;
};
template <typename TValue, typename TSpec>
struct Reference< Map<TValue, VectorSet<TSpec> > >
{
    typedef Map<TValue, VectorSet<TSpec> > TMap_;
    typedef typename Iterator<TMap_>::Type TIterator_;
    typedef Proxy<IteratorProxy<TIterator_> > Type;
};


template <typename TValue, typename TSpec>
struct Size< Map<TValue, VectorSet<TSpec> > >:
    Size<typename VectorSetElements_< Map<TValue, VectorSet<TSpec> > >::Type> {};

template <typename TValue, typename TSpec>
struct Key< Map<TValue, VectorSet<TSpec> > > :
    Key<TValue> {};

template <typename TValue, typename TSpec>
struct Cargo< Map<TValue, VectorSet<TSpec> > > :
    Cargo<TValue> {};


//////////////////////////////////////////////////////////////////////////////

template <typename TValue, typename TSpec>
inline typename Size< Map<TValue, VectorSet<TSpec> > >::Type
length(Map<TValue, VectorSet<TSpec>  > const &set)
{
    return set.data_length;
}

//////////////////////////////////////////////////////////////////////////////
//

struct VectorSetIterator;

template <typename TMap>
class Iter<TMap,  VectorSetIterator>
{
public:
    typedef typename Size<TMap>::Type TSize;
    typedef typename VectorSetElements_<TMap>::Type TElements;
    typedef typename Iterator<TElements, Rooted>::Type TElementsIterator;
    typedef typename Value<TElementsIterator>::Type TValue;

    TElementsIterator data_iterator;

    Iter()
    {    }
    Iter(TElementsIterator it)
        : data_iterator(it)
    {
        if (!atEnd(*this) && !(*data_iterator).data_valid)
        {
            goNext(*this);
        }
    }
    Iter(Iter const & other)
        : data_iterator(other.data_iterator)
    {    }
    ~Iter()
    {    }

    Iter const & operator = (Iter const & other)
    {
        data_iterator = other.data_iterator;
        return *this;
    }

    operator bool()
    {
        return data_iterator != 0;
    }

    TValue *
    operator -> () const
    {
        return & value(data_iterator);
    }

};

//////////////////////////////////////////////////////////////////////////////

template <typename TObject, typename TSpec, typename TIteratorSpec>
struct Iterator< Map<TObject, VectorSet<TSpec> >, TIteratorSpec >
{
    typedef Map<TObject, VectorSet<TSpec> > TMap_;
    typedef Iter<TMap_, VectorSetIterator> Type;
};

//////////////////////////////////////////////////////////////////////////////

template <typename TValue, typename TSpec>
inline void
clear(Map<TValue, VectorSet<TSpec> > & me)
{
    typedef Map<TValue, VectorSet<TSpec> > TMap;
    typedef typename VectorSetElements_<TMap>::Type TElements;
    typedef typename Iterator<TElements>::Type TIterator;
    TIterator it_end = end(me.data_elements);
    for (TIterator it = begin(me.data_elements); it != it_end; ++it)
    {
        (*it).data_valid = false;
    }
    me.data_length = 0;
}

//////////////////////////////////////////////////////////////////////////////

template <typename TCargo>
struct VectorSetInsert_
{
    template <typename TValue, typename TSpec, typename TElement>
    static inline void
    insert_(Map<TValue, VectorSet<TSpec> > & set,
        TElement const &element)
    {
        if (!set.data_elements[ordValue(key(element))].data_valid)
        {
            ++set.data_length;
            set.data_elements[ordValue(key(element))].data_valid = true;
        }
        set.data_elements[ordValue(key(element))].data_cargo = cargo(element);
    }
};
template <>
struct VectorSetInsert_<Nothing>
{
    template <typename TValue, typename TSpec, typename TElement>
    static inline void
    insert_(Map<TValue, VectorSet<TSpec> > & set,
        TElement const &element)
    {
        if (!set.data_elements[ordValue(key(element))].data_valid)
        {
            ++set.data_length;
            set.data_elements[ordValue(key(element))].data_valid = true;
        }
    }
};

template <typename TValue, typename TSpec, typename TElement>
inline void
insert(Map<TValue, VectorSet<TSpec> > & set,
       TElement const &element)
{
    typedef Map<TValue, VectorSet<TSpec> > TMap;
    typedef typename Cargo<TMap>::Type TCargo;
    VectorSetInsert_<TCargo>::insert_(set, element);
}
template <typename TValue, typename TSpec, typename TKey, typename TCargo>
inline void
insert(Map<TValue, VectorSet<TSpec> > &set,
       TKey const & _key,
       TCargo const & _cargo)
{
    typedef Map<TValue, VectorSet<TSpec> > TMap;
    typedef typename Value<TMap>::Type TValue2;
    TValue2 new_val;
    key(new_val) = _key;
    cargo(new_val) = _cargo;
    VectorSetInsert_<TCargo>::insert_(set, new_val);
}

//////////////////////////////////////////////////////////////////////////////

template <typename TValue, typename TSpec, typename T>
inline void
erase(Map<TValue, VectorSet<TSpec> > &set,
      T const & tokill)  //can be an element, a key, or an iterator
{
    if (set.data_elements[ordValue(key(tokill))].data_valid)
    {
        --set.data_length;
        set.data_elements[ordValue(key(tokill))].data_valid = false;
    }
}

//////////////////////////////////////////////////////////////////////////////

template <typename TKey, typename TValue, typename TSpec>
inline bool
hasKey(Map<TValue, VectorSet<TSpec> > const &set, TKey const &key)
{
    return set.data_elements[ordValue(key)].data_valid;
}

//////////////////////////////////////////////////////////////////////////////

template <typename TKey, typename TKey2, typename TSpec>
inline typename Iterator< Map<TKey2, VectorSet<TSpec> > >::Type
find(Map<TKey2, VectorSet<TSpec> > & set,
     TKey const &key)
{
    typedef Map<TKey2, VectorSet<TSpec> > TMap;
    typedef typename Iterator<TMap>::Type TIterator;
    return TIterator(begin(set.data_elements) + ordValue(key));
}

//////////////////////////////////////////////////////////////////////////////

template <typename TKey, typename TSpec, typename TKey2>
inline typename Cargo< Map<TKey, VectorSet<TSpec> > >::Type &
cargo(Map<TKey, VectorSet<TSpec> > & set,
      TKey2 const & _key)
{
    if (!hasKey(set, _key))
    {
        typedef Map<TKey, VectorSet<TSpec> > TMap;
        typedef typename Value<TMap>::Type TValue;
        TValue new_value;
        key(new_value) = _key;
        insert(set, new_value);
    }
    return set.data_elements[ordValue(_key)].data_cargo;
}

//////////////////////////////////////////////////////////////////////////////

template <typename TElement, typename TSpec>
inline typename Iterator< Map<TElement, VectorSet<TSpec> > >::Type
begin(Map<TElement, VectorSet<TSpec> > & set)
{
    typedef Map<TElement, VectorSet<TSpec> > TMap;
    typedef typename Iterator<TMap>::Type TIterator;
    return TIterator(begin(set.data_elements));
}

//////////////////////////////////////////////////////////////////////////////

template <typename TElement, typename TSpec>
inline typename Iterator< Map<TElement, VectorSet<TSpec> > >::Type
end(Map<TElement, VectorSet<TSpec> > &set)
{
    typedef Map<TElement, VectorSet<TSpec> > TMap;
    typedef typename Iterator<TMap>::Type TIterator;
    return TIterator(end(set.data_elements));
}

//////////////////////////////////////////////////////////////////////////////

template <typename TSet>
inline bool
operator==(Iter<TSet, VectorSetIterator> const &a,
           Iter<TSet, VectorSetIterator> const &b)
{
    return a.data_iterator == b.data_iterator;
}
template <typename TSet>
inline bool
operator!=(Iter<TSet, VectorSetIterator> const &a,
           Iter<TSet, VectorSetIterator> const &b)
{
    return a.data_iterator != b.data_iterator;
}

//////////////////////////////////////////////////////////////////////////////

template <typename TSet>
inline bool
atEnd(Iter<TSet, VectorSetIterator> & a)
{
    return atEnd(a.data_iterator);
}
template <typename TSet>
inline bool
atEnd(Iter<TSet, VectorSetIterator> const & a)
{
    return atEnd(a.data_iterator);
}

//////////////////////////////////////////////////////////////////////////////

template <typename TSet>
inline typename Key<TSet>::Type
key(Iter<TSet, VectorSetIterator> &it)
{
    return position(it.data_iterator);
}
template <typename TSet>
inline typename Key<TSet>::Type
key(Iter<TSet, VectorSetIterator> const &it)
{
    return position(it.data_iterator);
}

//////////////////////////////////////////////////////////////////////////////

template <typename TSet>
inline typename Cargo<TSet>::Type &
cargo(Iter<TSet, VectorSetIterator> &it) {
    return (*it.data_iterator).data_cargo;
}
template <typename TSet>
inline typename Cargo<TSet>::Type &
cargo(Iter<TSet, VectorSetIterator> const &it) {
    return (*it.data_iterator).data_cargo;
}

//////////////////////////////////////////////////////////////////////////////

template <typename TSet>
inline typename Reference<TSet>::Type
value(Iter<TSet, VectorSetIterator> &it)
{
    typedef typename Reference<TSet>::Type TReference;
    return TReference(it);
}
template <typename TSet>
inline typename Reference<TSet>::Type
value(Iter<TSet, VectorSetIterator> const &it)
{
    typedef typename Reference<TSet>::Type TReference;
    return TReference(it);
}

//////////////////////////////////////////////////////////////////////////////


template <typename TSet>
inline void
goNext(Iter<TSet, VectorSetIterator> & it)
{
    ++it.data_iterator;
    while (!atEnd(it.data_iterator))
    {
        if ((*it.data_iterator).data_valid)
            break;
        ++it.data_iterator;
    }
}

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////

template <typename TSet>
inline typename Key<TSet>::Type
key(Proxy<IteratorProxy<Iter<TSet, VectorSetIterator> > > & pr)
{
    return key(iter(pr));
}
template <typename TSet>
inline typename Key<TSet>::Type
key(Proxy<IteratorProxy<Iter<TSet, VectorSetIterator> > > const & pr)
{
    return key(iter(pr));
}

//////////////////////////////////////////////////////////////////////////////

template <typename TSet>
inline typename Cargo<TSet>::Type &
cargo(Proxy<IteratorProxy<Iter<TSet, VectorSetIterator> > >& pr)
{
    return cargo(iter(pr));
}
template <typename TSet>
inline typename Cargo<TSet>::Type &
cargo(Proxy<IteratorProxy<Iter<TSet, VectorSetIterator> > > const & pr)
{
    return cargo(iter(pr));
}


//////////////////////////////////////////////////////////////////////////////

}

#endif
