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
// Author: Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>
// ==========================================================================
// Implementation of the Union-Find data structure.
// ==========================================================================

#ifndef SEQAN_MISC_MISC_UNION_FIND_H_
#define SEQAN_MISC_MISC_UNION_FIND_H_

// TODO(holtgrew): Comprehensive tests. Currently, there are *some* tests in the tests for graph_algorithm.

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Classes, Structs, Enums, Tags
// ============================================================================

/**
.Class.UnionFind
..cat:Miscellaneous
..summary:Union-Find data structure.
..signature:UnionFind<T>
..param.T:The type the data structure operates on.
..remarks:The data structure uses union by rank and path compression to achieve almost linear running time.
..remarks:Note that internally, T is used signed so not the whole range might be available.
..include:seqan/misc/misc_union_find.h
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

///.Metafunction.Value.param.T.type:Class.UnionFind
///.Metafunction.Value.class:Class.UnionFind

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

///.Metafunction.GetValue.param.T.type:Class.UnionFind
///.Metafunction.GetValue.class:Class.UnionFind

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

///.Metafunction.Size.param.T.type:Class.UnionFind
///.Metafunction.Size.class:Class.UnionFind

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

///.Function.clear.param.object.type:Class.UnionFind
///.Function.clear.class:Class.UnionFind

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

///.Function.length.param.object.type:Class.UnionFind
///.Function.length.class:Class.UnionFind

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

///.Function.reserve.param.object.type:Class.UnionFind
///.Function.reserve.class:Class.UnionFind

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

///.Function.resize.param.object.type:Class.UnionFind
///.Function.resize.class:Class.UnionFind
///.Function.resize.remarks:If $pm$ is of the @Class.UnionFind@, $value$ will automatically be set to -1.

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

///.Function.resizeVertexMap.param.pm.type:Class.UnionFind
///.Function.resizeVertexMap.class:Class.UnionFind
///.Function.resizeVertexMap.remarks:If $pm$ is of type @Class.UnionFind@ then the $prototype$ parameter is not available.

template <typename TSpec, typename TValue>
inline
typename Size<UnionFind<TValue> >::Type
resizeVertexMap(Graph<TSpec> const & g,
                 UnionFind<TValue> & unionFind)
{
    clear(unionFind);
    return resize(unionFind, numVertices(g));
}

// ----------------------------------------------------------------------------
// Function find()
// ----------------------------------------------------------------------------

/**
.Function.findSet:
..cat:Miscellaneous
..summary:Return set identifier, given an element identifier.
..signature:find(unionFind, query)
..class:Class.UnionFind
..param.unionFind:The Union-Find data structure query.
...type:Class.UnionFind
..param.query:The value to query for.
..returns:The set identifier for the given element id.
..see:Function.joinSets
..include:seqan/misc/misc_union_find.h
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

/**
.Function.joinSets
..cat:Miscellaneous
..summary:UNION() operation for UF data structure.
..signature:joinSets(unionFind, left, right)
..class:Class.UnionFind
..param.unionFind:The type the data structure operates on.
...type:Class.UnionFind
..param.left:Representant of the left set to union.
..param.right:Representant of the right set to union.
..remarks:This function is called $join$ and not $union$ since $union$ is a reserved keyword in the C and C++ programming languages.
..remarks:Note that you most likely want to put return values of $findSet()$ as the values for $left$ and $right$.
..see:Function.findSet
..include:seqan/misc/misc_union_find.h
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

#endif // #ifndef SEQAN_MISC_MISC_UNION_FIND_H_
