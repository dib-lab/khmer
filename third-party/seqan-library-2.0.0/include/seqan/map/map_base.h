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

#ifndef SEQAN_HEADER_MAP_BASE_H
#define SEQAN_HEADER_MAP_BASE_H


//////////////////////////////////////////////////////////////////////////////

namespace SEQAN_NAMESPACE_MAIN
{

// ==========================================================================
// Forwards
// ==========================================================================

#if !defined(_MSC_VER) || _MSC_VER <= 1600

template <typename TKey, typename TCargo, typename TCompare, typename TAlloc, typename TKey2>
inline typename Cargo< std::map<TKey,TCargo, TCompare, TAlloc> >::Type &
cargo(std::map<TKey,TCargo, TCompare, TAlloc> & me, TKey2 const & _key);

#endif  // #if !defined(_MSC_VER) || _MSC_VER <= 1600

//////////////////////////////////////////////////////////////////////////////
//insertion tags

template <typename TSpec = Default>
struct Skiplist;

/*!
 * @class Map
 * @headerfile <seqan/map.h>
 * @brief Set/dictionary container.
 *
 * @signature template <typename TValue, typename TSpec>
 *            class Map;
 *
 * @tparam TSpec  The specializing type. Default: @link Skiplist @endlink
 * @tparam TValue Type of values that are stored in the map. Use a Pair<Key, Cargo> to implement a dictionary
 *                mapping from <tt>Key</tt> to <tt>Cargo</tt>.
 */

template <typename TElement, typename TSpec = Skiplist<> >
class Map;

/*!
 * @fn Map#length
 * @headerfile <seqan/map.h>
 * @brief Return number of elements in map.
 *
 * @signature TSize length(map);
 *
 * @param[in] map The Map to query for its size.
 * @return TSize The number of elements in the map.
 */

//////////////////////////////////////////////////////////////////////////////
// In SeqAn sets and maps store elements as pairs of (key,cargo)
// the elements of sets without objects are the keys.
//////////////////////////////////////////////////////////////////////////////

/*moved to basic_aggregates.h
template <typename TKey, typename TObject, typename TSpec>
struct Key< Pair<TKey, TObject, TSpec> >
{
    typedef TKey Type;
};

template <typename TKey, typename TCargo, typename TSpec>
struct Cargo< Pair<TKey, TCargo, TSpec> >
{
    typedef TCargo Type;
};
*/

//////////////////////////////////////////////////////////////////////////////
// Type for mapValue function that implements [] for map types

template <typename TMap, typename TCargo>
struct MapValueImpl_
{
    typedef TCargo & Type;
};
template <typename TMap>
struct MapValueImpl_<TMap, Nothing>
{
    typedef bool Type;
};

/*!
 * @mfn Map#MapValue
 * @brief Type of the map value type.
 *
 * @signature MapValue<T>::Type
 * @tparam T A map type. Types: Map
 *
 * @return Type the map value type.
 */

template <typename TMap>
struct MapValue :
    MapValueImpl_< TMap, typename Cargo<TMap>::Type >
{
};



template <typename TCargo>
struct ImplMapValue_
{
    template <typename TMap, typename TKey2>
    static inline TCargo &
    mapValue_(TMap & me,
        TKey2 const & _key)
    {
        return cargo(me, _key);
    }
};

template <>
struct ImplMapValue_<Nothing>
{
    template <typename TMap, typename TKey2>
    static inline bool
    mapValue_(TMap & me,
        TKey2 const & _key)
    {
        return hasKey(me, _key);
    }
};

/*!
 * @fn Map#mapValue
 * @brief Subscript <tt>operator[]</tt> of maps.
 *
 * @signature TMapValue mapValue(map, key);
 *
 * @param[in,out] map A map. Types: Map
 * @param[in]     key A key.
 *
 * @return TMapValue If <tt>map</tt> is a set: The same as Map#hasKey.  If <tt>map</tt> is a dictionary: The same as
 *                   Map#value.
 *
 * @section Remarks
 *
 * Usually, Map#value implements the subscript operator <tt>[ ]</tt>, but for maps, this operator is implemented in
 * <tt>mapValue</tt>. The semantic of this operator depends on the kind of map: If the map has a Cargo.cargo, than
 * <tt>mapValue(map, key)</tt> returns the cargo of the (first) value in the map of the given key. If the map has no
 * Cargo.cargo, than the function returns a <tt>true</tt>, if <tt>key</tt> is in <tt>map</tt>, or <tt>false</tt>
 * otherwise.
 *
 * @section Remarks
 *
 * There is no way to create a set of Pair, since it is always interpreted as a key/value pair.  If you need a key type
 * that holds two members, define your own key type.
 *
 * You may overload Key and Cargo for your own value type in order to define, what part of your value type is used as
 * key and what as cargo.
 */

template <typename TMap, typename TKey>
inline typename MapValue<TMap>::Type
mapValue(TMap & me,
         TKey const & _key)
{
    typedef typename Cargo<TMap>::Type TCargo;
    return ImplMapValue_<TCargo>::mapValue_(me, _key);
}

//////////////////////////////////////////////////////////////////////////////

template <typename TElement>
inline TElement &
key(TElement & element)
{
    return element;
}
template <typename TElement>
inline TElement const &
key(TElement const & element)
{
    return element;
}

template <typename TKey, typename TObject, typename TSpec>
inline TKey &
key(Pair<TKey, TObject, TSpec> & element)
{
    return element.i1;
}
template <typename TKey, typename TObject, typename TSpec>
inline TKey const &
key(Pair<TKey, TObject, TSpec> const & element)
{
    return element.i1;
}

//////////////////////////////////////////////////////////////////////////////

template <typename TElement, typename TSource>
inline void
setKey(TElement & element,
       TSource const & source)
{
    element = source;
}
template <typename TKey, typename TObject, typename TSpec, typename TSource>
inline void
setKey(Pair<TKey, TObject, TSpec> & element,
       TSource const & source)
{
    element.i1 = source;
}

//////////////////////////////////////////////////////////////////////////////
//no default cargo function

template <typename TKey, typename TObject, typename TSpec>
inline TObject &
cargo(Pair<TKey, TObject, TSpec> & element)
{
    return element.i2;
}
template <typename TKey, typename TObject, typename TSpec>
inline TObject const &
cargo(Pair<TKey, TObject, TSpec> const & element)
{
    return element.i2;
}

//////////////////////////////////////////////////////////////////////////////

template <typename TKey, typename TObject, typename TSpec, typename TSource>
inline void
setCargo(Pair<TKey, TObject, TSpec> & element,
       TSource const & source)
{
    element.i2 = source;
}

//////////////////////////////////////////////////////////////////////////////

}

#endif

