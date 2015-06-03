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

#ifndef INCLUDE_SEQAN_GRAPH_TYPES_PROPERTY_MAP_INTERNAL_POINTER_H_
#define INCLUDE_SEQAN_GRAPH_TYPES_PROPERTY_MAP_INTERNAL_POINTER_H_

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// --------------------------------------------------------------------------
// Class InternalPointerPropertyMap
// --------------------------------------------------------------------------

/*!
 * @class InternalPointerPropertyMap InternalPointerPropertyMap
 * @implements PropertyMapConcept
 * @headerfile <seqan/graph_types.h>
 * @brief An internal property map using pointer to members.
 *
 * @signature template <typename TMemberPointer, TMemberPointer const MEMBER_POINTER>
 *            class InternalPointerPropertyMap;
 *
 * @tparam TMemberPointer  A pointer to a member type.
 * @tparam MEMBER_POINTER  A pointer to a type member of type <tt>TMemberPointer</tt>.
 *
 * Internal property maps are used to access internal edge cargos that are structs or classes.
 */

template <typename TPointer, TPointer const MEMBER_POINTER>
struct InternalPointerPropertyMap
{};

template <typename TPointer, TPointer const MEMBER_POINTER>
SEQAN_CONCEPT_IMPL((InternalPointerPropertyMap<TPointer, MEMBER_POINTER>), (PropertyMapConcept));

// ============================================================================
// Metafunctions
// ============================================================================

// --------------------------------------------------------------------------
// Metafunction InternalPointerPropertyMap#Value
// --------------------------------------------------------------------------

template <typename TClass, typename TValue, TValue TClass::* MEMBER_POINTER>
struct Value<InternalPointerPropertyMap<TValue TClass::*, MEMBER_POINTER> >
{
    typedef TValue Type;
};

template <typename TClass, typename TValue, TValue TClass::* MEMBER_POINTER>
struct Value<InternalPointerPropertyMap<TValue TClass::*, MEMBER_POINTER> const> :
    Value<InternalPointerPropertyMap<TValue TClass::*, MEMBER_POINTER> >
{};

// --------------------------------------------------------------------------
// Metafunction InternalPointerPropertyMap#GetValue
// --------------------------------------------------------------------------

template <typename TClass, typename TValue, TValue TClass::* MEMBER_POINTER>
struct GetValue<InternalPointerPropertyMap<TValue TClass::*, MEMBER_POINTER> >
{
    typedef TValue Type;
};

template <typename TClass, typename TValue, TValue TClass::* MEMBER_POINTER>
struct GetValue<InternalPointerPropertyMap<TValue TClass::*, MEMBER_POINTER> const>
{
    typedef TValue Type;
};

// --------------------------------------------------------------------------
// Metafunction InternalPointerPropertyMap#Reference
// --------------------------------------------------------------------------

template <typename TClass, typename TValue, TValue TClass::* MEMBER_POINTER>
struct Reference<InternalPointerPropertyMap<TValue TClass::*, MEMBER_POINTER> const>
{
    typedef TValue const & Type;
};

template <typename TClass, typename TValue, TValue TClass::* MEMBER_POINTER>
struct Reference<InternalPointerPropertyMap<TValue TClass::*, MEMBER_POINTER> >
{
    typedef TValue & Type;
};

// ============================================================================
// Functions
// ============================================================================

// --------------------------------------------------------------------------
// Function Graph#resizeEdgeMap()
// --------------------------------------------------------------------------

// Overload for InternalPointerPropertyMap.

template <typename TPropMap, TPropMap const INSTANCE, typename TSpec>
void resizeEdgeMap(InternalPointerPropertyMap<TPropMap, INSTANCE> &, Graph<TSpec> const &)
{}

// --------------------------------------------------------------------------
// Function InternalPointerPropertyMap#resize()
// --------------------------------------------------------------------------

template <typename TPropertyMap, TPropertyMap INSTANCE, typename TSize, typename TPrototype>
void resize(InternalPointerPropertyMap<TPropertyMap, INSTANCE> &, TSize, TPrototype)
{}

// --------------------------------------------------------------------------
// Function InternalPointerPropertyMap#assignProperty()
// --------------------------------------------------------------------------

/*!
 * @fn InternalPointerPropertyMap#assignProperty:
 * @brief Assigns a property to an item in the property map.
 * @signature void assignProperty(pm, d, val)
 *
 * @param[in,out] pm  The InternalPointerPropertyMap to assign into.
 * @param[in]     d   A vertex or edge descriptor that identifies the item in the property map.
 *                    Types: @link VertexDescriptor @endlink, @link Graph#EdgeDescriptor @endlink
 * @param[in]     val The new value, where thg type of the new value must match the value type of the property map.
 */

template <typename TClass, typename TValue, TValue TClass::* MEMBER_POINTER, typename TEdgeDescriptor>
void assignProperty(InternalPointerPropertyMap<TValue TClass::*, MEMBER_POINTER> &,
                    TEdgeDescriptor const e,
                    TValue const val)
{
    (cargo(e)).*MEMBER_POINTER = val;
}

// --------------------------------------------------------------------------
// Function InternalPointerPropertyMap#property()
// --------------------------------------------------------------------------

/*!
 * @fn InternalPointerPropertyMap#property:
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

template <typename TClass, typename TValue, TValue TClass::* MEMBER_POINTER, typename TEdgeDescriptor>
typename Reference<InternalPointerPropertyMap<TValue TClass::*, MEMBER_POINTER> >::Type
property(InternalPointerPropertyMap<TValue TClass::*, MEMBER_POINTER>&,
         TEdgeDescriptor const e)
{
    return (cargo(e)).*MEMBER_POINTER;
}

template <typename TClass, typename TValue, TValue TClass::* MEMBER_POINTER, typename TEdgeDescriptor>
typename Reference<InternalPointerPropertyMap<TValue TClass::*, MEMBER_POINTER> const>::Type
property(InternalPointerPropertyMap<TValue TClass::*, MEMBER_POINTER> const &,
         TEdgeDescriptor const e)
{
    return (cargo(e)).*MEMBER_POINTER;
}

// --------------------------------------------------------------------------
// Function InternalPointerPropertyMap#getProperty()
// --------------------------------------------------------------------------

/*!
 * @fn InternalPointerPropertyMap#getProperty
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

template <typename TClass, typename TValue, TValue TClass::* MEMBER_POINTER, typename TEdgeDescriptor>
typename GetValue<InternalPointerPropertyMap<TValue TClass::*, MEMBER_POINTER> >::Type
getProperty(InternalPointerPropertyMap<TValue TClass::*, MEMBER_POINTER> const &,
            TEdgeDescriptor const e)
{
    return (getCargo(e)).*MEMBER_POINTER;
}

}  // namespace seqan

#endif  // #ifndef INCLUDE_SEQAN_GRAPH_TYPES_PROPERTY_MAP_INTERNAL_POINTER_H_
