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
// Author: Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>
// ==========================================================================
// Implementation of the Union-Find data structure.
// ==========================================================================

#ifndef SEQAN_MISC_UNION_FIND_H_
#define SEQAN_MISC_UNION_FIND_H_

#include <seqan/basic.h>
#include <seqan/sequence.h>
// TODO(holtgrew): Comprehensive tests. Currently, there are *some* tests in the tests for graph_algorithm.

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

template <typename TSpec> class Graph;

// ============================================================================
// Classes, Structs, Enums, Tags
// ============================================================================

/*!
 * @class UnionFind
 * @headerfile <seqan/misc/union_find.h>
 * @brief Union-Find data structure.
 *
 * @signature template <typename T>
 *            class UnionFind;
 *
 * @tparam T The integer type the data structure operates on.
 *
 * The data structure uses union by rank and path compresison to achieve almost linear running time.
 *
 * Note that internally T is used as signed, so not the whole range is available if T is unsigned.
 */

template <typename TValue>
class UnionFind
{
public:
    typedef typename MakeSigned_<TValue>::Type TValue_;
    String<TValue_> _values;

    UnionFind() {}
};

// ============================================================================
// Metafunctions
// ============================================================================

// ----------------------------------------------------------------------------
// Metafunction Value
// ----------------------------------------------------------------------------

/*!
 * @mfn UnionFind#Value
 * @brief Returns the value type for the given UnionFind specialization.
 *
 * @signature Value<TUnionFind>::Type
 *
 * @tparam TUnionFind The UnionFind specialization to query for its value type.
 */

template <typename TValue>
struct Value<UnionFind<TValue> >
{
    typedef typename MakeSigned_<TValue>::Type TValue_;
    typedef String<TValue_> TString_;
    typedef typename Value<TString_>::Type Type;
};

template <typename TValue>
struct Value<UnionFind<TValue> const>
{
    typedef typename MakeSigned_<TValue>::Type TValue_;
    typedef String<TValue_> const TString_;
    typedef typename Value<TString_>::Type Type;
};

// ----------------------------------------------------------------------------
// Metafunction GetValue
// ----------------------------------------------------------------------------

/*!
 * @mfn UnionFind#GetValue
 * @brief Returns the get-value type for the given UnionFind specialization.
 *
 * @signature GetValue<TUnionFind>::Type
 *
 * @tparam TUnionFind The UnionFind specialization to query for its get-value type.
 */

template <typename TValue>
struct GetValue<UnionFind<TValue> >
{
    typedef typename MakeSigned_<TValue>::Type TValue_;
    typedef String<TValue_> TString_;
    typedef typename GetValue<TString_>::Type Type;
};

template <typename TValue>
struct GetValue<UnionFind<TValue> const>
{
    typedef typename MakeSigned_<TValue>::Type TValue_;
    typedef String<TValue_> const TString_;
    typedef typename GetValue<TString_>::Type Type;
};

// ----------------------------------------------------------------------------
// Metafunction Size
// ----------------------------------------------------------------------------

/*!
 * @mfn UnionFind#Size
 * @brief Returns the size type for the given UnionFind specialization.
 *
 * @signature Size<TUnionFind>::Type
 *
 * @tparam TUnionFind The UnionFind specialization to query for its size type.
 */

template <typename TValue>
struct Size<UnionFind<TValue> >
{
    typedef typename MakeSigned_<TValue>::Type TValue_;
    typedef String<TValue_> TString_;
    typedef typename Size<TString_>::Type Type;
};

template <typename TValue>
struct Size<UnionFind<TValue> const>
{
    typedef typename MakeSigned_<TValue>::Type TValue_;
    typedef String<TValue_> const TString_;
    typedef typename Size<TString_>::Type Type;
};

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function clear()
// ----------------------------------------------------------------------------

/*!
 * @fn UnionFind#clear
 * @brief Clear the Union-Find data structure.
 *
 * @signature void clear(uf);
 *
 * @param[in,out] uf The Union-Find object to clear
 */

template <typename TValue>
inline
void
clear(UnionFind<TValue> & unionFind)
{
    clear(unionFind._values);
}

// ----------------------------------------------------------------------------
// Function length()
// ----------------------------------------------------------------------------

/*!
 * @fn UnionFind#length
 * @brief Returns the number of entries in the Union-Find object.
 *
 * @signature TSize length(uf);
 *
 * @param[in] uf The Union-Find object to query.
 *
 * @return TSize The length of the Union-Find object  (Metafunction: @link UnionFind#Size @endlink).
 */

template <typename TValue>
inline
typename Size<UnionFind<TValue> >::Type
length(UnionFind<TValue> const & unionFind)
{
    return length(unionFind._values);
}

// ----------------------------------------------------------------------------
// Function reserve()
// ----------------------------------------------------------------------------

// TODO(holtgrew): Remove resize tag, return value?

/*!
 * @fn UnionFind#reverse
 * @brief Reserve memory for the Union-Find object.
 *
 * @signature void reserve(uf, size, tag);
 *
 * @param[in,out] uf    The Union-Find object to reserve memory in.
 * @param[in]     size  The number of elements to reserve.
 * @param[in]     tag   The tag to use for reserving (defaults to <tt>Generous()</tt>).
 */

template <typename TValue, typename TSize, typename TTag>
inline
typename Size<UnionFind<TValue> >::Type
reserve(UnionFind<TValue> & unionFind,
       TSize const & newSize,
       TTag const & tag)
{
    return reserve(unionFind._values, newSize, tag);
}

// ----------------------------------------------------------------------------
// Function resize()
// ----------------------------------------------------------------------------

/*!
 * @fn UnionFind#resize
 * @brief Allocate number of elements for the Union-Find object.
 *
 * The UF dat structure is resized to the given <tt>size</tt>, the value for each element is set to -1, i.e. they
 * are singletons by default.
 *
 * @signature TSize resize(uf, size, tag);
 *
 * @param[in,out] uf    The Union-Find object to resize.
 * @param[in]     size  The number of elements to reserve.
 * @param[in]     tag   The tag to use for reserving (defaults to <tt>Generous()</tt>).
 *
 * @return TSize The new length of the Union-Find object  (Metafunction: @link UnionFind#Size @endlink).
 */

template <typename TValue, typename TSize, typename TTag>
inline
typename Size<UnionFind<TValue> >::Type
resize(UnionFind<TValue> & unionFind,
       TSize const & newSize,
       TTag const & tag)
{
    return resize(unionFind._values, newSize, -1, tag);
}

// ----------------------------------------------------------------------------
// Function resizeVertexMap()
// ----------------------------------------------------------------------------

// TODO(holtgrew): Change parameter order for resizeVertexMap()!

/*!
 * @fn UnionFind#resizeVertexMap
 * @brief Resize Union-Find data structure to appropriate size for a vertex map.
 *
 * @signature TSize resizeVertexMap(uf, g);
 *
 * @param[in,out] uf The Union-Find object to resize.
 * @param[in]     g  The graph to use for getting the vertex number.
 *
 * @return TSize New size of the vertex map  (Metafunction: @link UnionFind#Size @endlink).
 */

template <typename TSpec, typename TValue>
typename Size<UnionFind<TValue> >::Type
resizeVertexMap(UnionFind<TValue> & unionFind,
                Graph<TSpec> const & g)
{
    clear(unionFind);
    return resize(unionFind, numVertices(g));
}

// ----------------------------------------------------------------------------
// Function findSet()
// ----------------------------------------------------------------------------

/*!
 * @fn UnionFind#findSet
 * @brief Return set identifier, given an element identifier.
 *
 * @signature TValue findSet(uf, q);
 *
 * @param[in] uf  The Union-Find object to query.
 * @param[in] q   The value to get the set identifier for.
 *
 * @return TValue Identifier of the set that <tt>q</tt> is in  (Metafunction: @link UnionFind#Value @endlink).
 *
 * @see UnionFind#joinSets
 */

template <typename TValue, typename TQuery>
inline
TValue
findSet(UnionFind<TValue> & unionFind,
        TQuery const & query)
{
    TValue j = query;
    while (unionFind._values[j] >= 0)
        j = unionFind._values[j];
    SEQAN_ASSERT_LT(unionFind._values[j], static_cast<int>(length(unionFind._values)));
    SEQAN_ASSERT_GEQ(j, static_cast<TValue>(0));
    SEQAN_ASSERT_LT(j, static_cast<TValue>(length(unionFind._values)));

    TValue i = query;
    while (unionFind._values[i] >= 0)
    {
        TValue tmp = i;
        i = unionFind._values[i];
        unionFind._values[tmp] = j;
    }

    return j;
}

// ----------------------------------------------------------------------------
// Function joinSets()
// ----------------------------------------------------------------------------

/*!
 * @fn UnionFind#joinSets
 * @brief UNION() operation for Union-Find data structure.
 *
 * @signature void joinSets(uf, left, right);
 *
 * @param[in,out] uf    The type the data structure operates on.
 * @param[in]     left  Representant of the left set to union.
 * @param[in]     right Representant of the right set to union.
 *
 * This function is called <tt>join</tt> and not <tt>union</tt> since <tt>union</tt> is a reserved keyword in the
 * C and C++ programming languages.
 *
 * Note that you most likely want to put return values of <tt>findSet()</tt> as the values for <tt>left</tt> and
 * <tt>right</tt>.
 *
 * @see UnionFind#findSet
 */

template <typename TValue, typename TLeft, typename TRight>
inline
void
joinSets(UnionFind<TValue> & unionFind,
         TLeft const & left,
         TRight const & right)
{
    if (left == right)
        return;
    TValue sum = unionFind._values[left] + unionFind._values[right];
    if (_abs(unionFind._values[left]) < _abs(unionFind._values[right])) {
        unionFind._values[left] = right;
        unionFind._values[right] = sum;
    } else {
        unionFind._values[right] = left;
        unionFind._values[left] = sum;
    }
}

}  // namespace seqan

#endif // #ifndef SEQAN_MISC_UNION_FIND_H_
