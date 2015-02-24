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

#ifndef INCLUDE_SEQAN_BASIC_PROPERTY_MAP_CONCEPT_H_
#define INCLUDE_SEQAN_BASIC_PROPERTY_MAP_CONCEPT_H_

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// ----------------------------------------------------------------------------
// Concept PropertyMapConcept
// ----------------------------------------------------------------------------

/*!
 * @concept PropertyMapConcept
 * @headerfile <seqan/graph_types.h>
 * @brief Concept for maps from contained elements (such as graph vertices or index nodes) to values.
 *
 * @signature concept PropertyMapConcept;
 */

/*!
 * @mfn PropertyMapConcept#Value
 * @brief Returns the value type of the property map.
 *
 * @signature Value<TPropertyMap>::Type
 *
 * @tparam TPropertyMap The property map to query.
 *
 * @return Type The element type of the container.
 *
 * @see PropertyMapConcept#Value
 */

/*!
 * @mfn PropertyMapConcept#GetValue
 * @brief Returns the get-value type of the property map.
 *
 * @signature GetValue<TPropertyMap>::Type
 *
 * @tparam TPropertyMap The property map to query.
 *
 * @return Type The get-value type of the container.
 *
 * @see PropertyMapConcept#GetValue
 */

/*!
 * @mfn PropertyMapConcept#Reference
 * @brief Returns the reference type of the property map.
 *
 * @signature Reference<TPropertyMap>::Type
 *
 * @tparam TPropertyMap The property map to query.
 *
 * @see PropertyMapConcept#Reference
 */

/*!
 * @fn PropertyMapConcept#assignProperty
 * @brief Assigns a property to an item in the property map.
 *
 * @signature void assignProperty(pm, d, val);
 *
 * @param[in,out] pm  The property map
 * @param[in]     d   A vertex or edge descriptor that identifies the item in the property map.
 * @param[in]     val The new value, where the type of the new value must match the value type of the property map.
 */

/*!
 * @fn PropertyMapConcept#property
 * @brief Accesses the property of an item in the property map.
 *
 * @signature TReference property(pm, d);
 *
 * @param[in,out] pm  The property map.
 * @param[in]     d   A vertex or edge descriptor that identifies the item in the property map.
 *
 * @return TReference Reference to the item in the property map of type @link Reference @endlink.
 */

/*!
 * @fn PropertyMapConcept#getProperty
 * @brief Get method for an item's property.
 *
 * @signature TGetValue getProperty(pm, d);
 *
 * @param[in,out] pm  The property map.
 * @param[in]     d   A vertex or edge descriptor that identifies the item in the property map.
 *
 * @return TGetValue Get-value to the item in the property map of type @link PropertyMapConcept#GetValue
 *                   GetValue @endlink.
 */

/*!
 * @fn PropertyMapConcept#resize
 * @brief Resize a sequence.
 *
 * @signature void resize(pm, len[, val]);
 *
 * @param[in,out] seq Sequence to resize.
 * @param[in]     len Length to resize <tt>seq</tt> to.
 * @param[in]     val When increasing the size, <tt>val</tt> is used to fill new entries.  When omitted,
 *                    <tt>TValue()</tt> is used where <tt>TValue</tt> is the @link ContainerConcept#Value @endlink
 *                    type of the sequence.
 *
 * @see StringConcept#resize
 */

SEQAN_CONCEPT(PropertyMapConcept, (TPropertyMap))
{
public:
    typedef typename Value<TPropertyMap>::Type                TValue;
    typedef typename GetValue<TPropertyMap>::Type             TGetValue;
    typedef typename Reference<TPropertyMap>::Type            TReference;

    SEQAN_CONCEPT_USAGE(PropertyMapConcept)
    {}
};

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

}  // namespace seqan

#endif  // #ifndef INCLUDE_SEQAN_BASIC_PROPERTY_MAP_CONCEPT_H_
