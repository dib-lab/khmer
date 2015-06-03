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
// Author: Tobias Rausch <rausch@embl.de>
// ==========================================================================

#ifndef INCLUDE_SEQAN_GRAPH_TYPES_PROPERTY_MAP_INTERNAL_H_
#define INCLUDE_SEQAN_GRAPH_TYPES_PROPERTY_MAP_INTERNAL_H_

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// --------------------------------------------------------------------------
// Class InternalPropertyMap
// --------------------------------------------------------------------------

/*!
 * @class InternalPropertyMap InternalPropertyMap
 * @implements PropertyMapConcept
 * @headerfile <seqan/graph_types.h>
 * @brief An internal property map with direct access to members.
 *
 * @signature template <typename TMember>
 *            class InternalPropertyMap;
 *
 * @tparam TMember The member type.
 *
 * Internal property maps are used to access internal edge cargos.
 */

template <typename TMember>
struct InternalPropertyMap
{};

template <typename TMember>
SEQAN_CONCEPT_IMPL((InternalPropertyMap<TMember>), (PropertyMapConcept));

// ============================================================================
// Metafunctions
// ============================================================================

// --------------------------------------------------------------------------
// Metafunction InternalPropertyMap#Value
// --------------------------------------------------------------------------

template <typename TValue>
struct Value<InternalPropertyMap<TValue> >
{
    typedef TValue Type;
};

template <typename TValue>
struct Value<InternalPropertyMap<TValue> const> : Value<InternalPropertyMap<TValue> >
{};

// --------------------------------------------------------------------------
// Metafunction PropertyMapConcept#GetValue
// --------------------------------------------------------------------------

template <typename TValue>
struct GetValue<InternalPropertyMap<TValue> >
{
    typedef TValue Type;
};

template <typename TValue>
struct GetValue<InternalPropertyMap<TValue> const>
{
    typedef TValue Type;
};

// --------------------------------------------------------------------------
// Metafunction InternalPropertyMap#Reference
// --------------------------------------------------------------------------

template <typename TValue>
struct Reference<InternalPropertyMap<TValue> const>
{
    typedef TValue const & Type;
};

template <typename TValue>
struct Reference<InternalPropertyMap<TValue> >
{
    typedef TValue & Type;
};

// ============================================================================
// Functions
// ============================================================================

// --------------------------------------------------------------------------
// Function Graph#resizeEdgeMap()
// --------------------------------------------------------------------------

// Overload for InternalPropertyMap.

template <typename TPropMap, typename TSpec>
void resizeEdgeMap(InternalPropertyMap<TPropMap> &, Graph<TSpec> const &)
{}

// --------------------------------------------------------------------------
// Function InternalPropertyMap#resize()
// --------------------------------------------------------------------------

template <typename TPropertyMap, typename TSize, typename TPrototype>
void resize(InternalPropertyMap<TPropertyMap> &, TSize, TPrototype)
{}

// --------------------------------------------------------------------------
// Function InternalPropertyMap#assignProperty()
// --------------------------------------------------------------------------

/*!
 * @fn InternalPropertyMap#assignProperty:
 * @brief Assigns a property to an item in the property map.
 * @signature void assignProperty(pm, d, val)
 *
 * @param[in,out] pm  The InternalPropertyMap to assign into.
 * @param[in]     d   A vertex or edge descriptor that identifies the item in the property map.
 *                    Types: @link VertexDescriptor @endlink, @link Graph#EdgeDescriptor @endlink
 * @param[in]     val The new value, where thg type of the new value must match the value type of the property map.
 */

template <typename TValue, typename TEdgeDescriptor>
void assignProperty(InternalPropertyMap<TValue> &,
                    TEdgeDescriptor const e,
                    TValue const val)
{
    cargo(e) = val;
}

// --------------------------------------------------------------------------
// Function InternalPropertyMap#property()
// --------------------------------------------------------------------------

/*!
 * @fn InternalPropertyMap#property:
 * @brief Accesses the property of an item in the property map.
 *
 * @signature TReference property(pm, d)
 *
 * @param[in,out] pm  The property map.
 * @param[in]     d   A vertex or edge descriptor that identifies the item in the property map.
 *                    Types: @link VertexDescriptor @endlink, @link Graph#EdgeDescriptor @endlink
 *
 * @return TReference @link PropertyMapConcept#Reference Reference @endlink to the item in the property map of type
 *                    @link Reference @endlink.
 */

template <typename TValue, typename TEdgeDescriptor>
typename Reference<InternalPropertyMap<TValue> const>::Type
property(InternalPropertyMap<TValue> const &, TEdgeDescriptor const e)
{
    return cargo(e);
}

template <typename TValue, typename TEdgeDescriptor>
typename Reference<InternalPropertyMap<TValue> >::Type
property(InternalPropertyMap<TValue>&, TEdgeDescriptor const e)
{
    return cargo(e);
}

// --------------------------------------------------------------------------
// Function InternalPropertyMap#getProperty()
// --------------------------------------------------------------------------

/*!
 * @fn InternalPropertyMap#getProperty
 * @brief Get method for an item's property.
 *
 * @signature TGetValue getProperty(pm, d)
 *
 * @param[in] pm  The property map.
 * @param[in] d   A vertex or edge descriptor that identifies the item in the property map.
 *               Types: @link VertexDescriptor @endlink, @link Graph#EdgeDescriptor @endlink
 *
 * @return TValue Get-value of the item in the property map of type @link PropertyMapConcept#GetValue GetValue @endlink.
 */

template <typename TValue, typename TEdgeDescriptor>
typename GetValue<InternalPropertyMap<TValue> >::Type
getProperty(InternalPropertyMap<TValue> const &, TEdgeDescriptor const e)
{
    return getCargo(e);
}

}  // namespace seqan

#endif  // #ifndef INCLUDE_SEQAN_GRAPH_TYPES_PROPERTY_MAP_INTERNAL_H_
