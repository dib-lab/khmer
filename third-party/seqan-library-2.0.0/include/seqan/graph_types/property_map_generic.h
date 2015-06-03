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
// Author: Tobias Rausch <rausch@embl.de>
// Author: Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>
// ==========================================================================
// Graph-related property map functions, defaults are used for general types
// implementing the PropertyMapConcept, e.g. containers such as SeqAn
// Strings.
// ==========================================================================

#ifndef INCLUDE_SEQAN_GRAPH_TYPES_PROPERTY_MAP_GENERIC_H_
#define INCLUDE_SEQAN_GRAPH_TYPES_PROPERTY_MAP_GENERIC_H_

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

// --------------------------------------------------------------------------
// Class Graph#resizeVertexMap
// --------------------------------------------------------------------------

/*!
 * @fn Graph#resizeVertexMap
 * @brief Initializes a vertex map.
 * @signature void resizeVertexMap(pm, g[, prototype]);
 *
 * @param[in,out] pm        A @link PropertyMapConcept property map @endlink.
 * @param[in]     g         A Graph.
 * @param[in]     prototype An optional prototype that is used for initializing the property map.
 */

template <typename TPropertyMap, typename TSpec>
SEQAN_FUNC_ENABLE_IF(Is<PropertyMapConcept<TPropertyMap> >, void)
resizeVertexMap(TPropertyMap & pm, Graph<TSpec> const & g)
{
    typedef typename Value<TPropertyMap>::Type TValue;
    resize(pm, getIdUpperBound(_getVertexIdManager(g)), TValue(), Generous());
}

template <typename TPropertyMap, typename TSpec, typename TPrototype>
SEQAN_FUNC_ENABLE_IF(Is<PropertyMapConcept<TPropertyMap> >, void)
resizeVertexMap(TPropertyMap & pm, Graph<TSpec> const & g, TPrototype const & prototype)
{
    resize(pm, getIdUpperBound(_getVertexIdManager(g)), prototype, Generous());
}

// --------------------------------------------------------------------------
// Class Graph#resizeEdgeMap
// --------------------------------------------------------------------------

/*!
 * @fn Graph#resizeEdgeMap
 * @brief Initializes an edge map.
 * @signature void resizeEdgeMap(pm, g[, prototype]);;
 * @param[in,out] pm        A @link PropertyMapConcept property map @endlink.
 * @param[in]     g         A Graph.
 * @param[in]     prototype An optional prototype that is used for initializing the property map.
 */

template <typename TPropertyMap, typename TSpec>
SEQAN_FUNC_ENABLE_IF(Is<PropertyMapConcept<TPropertyMap> >, void)
resizeEdgeMap(TPropertyMap & pm, Graph<TSpec> const & g)
{
    typedef typename Value<TPropertyMap>::Type TValue;
    resize(pm, getIdUpperBound(_getEdgeIdManager(g)), TValue(), Generous());
}

template <typename TPropertyMap, typename TSpec, typename TPrototype>
SEQAN_FUNC_ENABLE_IF(Is<PropertyMapConcept<TPropertyMap> >, void)
resizeEdgeMap(TPropertyMap & pm, Graph<TSpec> const & g, TPrototype const & prototype)
{
    resize(pm, getIdUpperBound(_getEdgeIdManager(g)), prototype, Generous());
}

// --------------------------------------------------------------------------
// Function Graph#assignVertexMap()
// --------------------------------------------------------------------------

/*!
 * @fn Graph#assignVertexMap
 * @brief Initializes a vertex map with values of an array.
 * @signature void assignVertexMap(g, pm, prop);
 *
 * @param[out] pm   A @link PropertyMapConcept property map @endlink.
 * @param[in]  g    A Graph.
 * @param[in]  prop An array with properties that are to be assigned to the items in the property map.
 *
 * For every vertex descriptor there must be an entry in the array.
 */

template <typename TSpec, typename TPropertyMap, typename TProperties>
SEQAN_FUNC_ENABLE_IF(Is<PropertyMapConcept<TPropertyMap> >, void)
assignVertexMap(TPropertyMap & pm,
                Graph<TSpec> const & g,
                TProperties const & prop)
{
    resize(pm, getIdUpperBound(_getVertexIdManager(g)), Generous());
    typedef Graph<TSpec> TGraph;
    typedef typename Iterator<TGraph, VertexIterator>::Type TVertexIterator;
    TVertexIterator it(g);
    for (; !atEnd(it); goNext(it))
        assignProperty(pm, getValue(it), getValue(prop, _getId(value(it))));
}

// --------------------------------------------------------------------------
// Function Graph#assignEdgeMap()
// --------------------------------------------------------------------------

/*!
 * @fn Graph#assignEdgeMap
 * @brief Initializes a vertex map with values of an array.
 * @signature void assignEdgeMap(g, pm, prop);
 *
 * @param[in]  pm   An @link PropertyMapConcept property map @endlink.
 * @param[out] prop An array with properties that are to be assigned to the items in the property map.
 * @param[in]  g    A Graph.
 *
 * For every edge id there must be an entry in the array.
 */

template <typename TSpec, typename TPropertyMap, typename TProperties>
SEQAN_FUNC_ENABLE_IF(Is<PropertyMapConcept<TPropertyMap> >, void)
assignEdgeMap(TPropertyMap & pm,
              Graph<TSpec> const & g,
              TProperties const & prop)
{
    resize(pm, getIdUpperBound(_getEdgeIdManager(g)), Generous());
    typedef Graph<TSpec> TGraph;
    typedef typename Iterator<TGraph, EdgeIterator>::Type TEdgeIterator;
    TEdgeIterator it(g);
    for (; !atEnd(it); goNext(it))
        assignProperty(pm, *it, prop[_getId(*it)]);
}

// --------------------------------------------------------------------------
// Class PropertyMapConcept#assignProperty
// --------------------------------------------------------------------------

template <typename TPropertyMap, typename TDescriptor, typename TValue>
SEQAN_FUNC_ENABLE_IF(Is<PropertyMapConcept<TPropertyMap> >, void)
assignProperty(TPropertyMap & pm, TDescriptor const d, TValue const val)
{
    assignValue(pm, _getId(d), val);
}

// --------------------------------------------------------------------------
// Class PropertyMapConcept#property
// --------------------------------------------------------------------------

template <typename TPropertyMap, typename TDescriptor>
SEQAN_FUNC_ENABLE_IF(Is<PropertyMapConcept<TPropertyMap> >,
                     typename Reference<TPropertyMap>::Type)
property(TPropertyMap & pm, TDescriptor const d)
{
    return value(pm, _getId(d));
}

template <typename TPropertyMap, typename TDescriptor>
SEQAN_FUNC_ENABLE_IF(Is<PropertyMapConcept<TPropertyMap> >,
                     typename Reference<TPropertyMap const>::Type)
property(TPropertyMap const & pm, TDescriptor const d)
{
    return value(pm, _getId(d));
}

// --------------------------------------------------------------------------
// Class PropertyMapConcept#getProperty
// --------------------------------------------------------------------------

template <typename TPropertyMap, typename TDescriptor>
SEQAN_FUNC_ENABLE_IF(Is<PropertyMapConcept<TPropertyMap> >,
                     typename GetValue<TPropertyMap const>::Type)
getProperty(TPropertyMap const & pm, TDescriptor const d)
{
    return getValue(pm, _getId(d));
}

}  // namespace seqan

#endif  // #ifndef INCLUDE_SEQAN_GRAPH_TYPES_PROPERTY_MAP_GENERIC_H_
