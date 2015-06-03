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
// Author: David Weese <david.weese@fu-berlin.de>
// ==========================================================================

#ifndef SEQAN_HEADER_STORE_READ_H
#define SEQAN_HEADER_STORE_READ_H

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////
// Read Store
//////////////////////////////////////////////////////////////////////////////

/*!
 * @class ReadStoreElement
 * @headerfile <seqan/store.h>
 * @brief Represents a single read (without sequence).
 *
 * @signature template <[typename TSpec]>
 *            struct ReadStoreElement;
 *
 * @tparam TSpec The specializing type.  Default: <tt>void</tt>.
 *
 *
 * @fn ReadStoreElement::ReadStoreElement
 * @signature ReadStoreElement::ReadStoreElement();
 * @brief Constructor.
 *
 * Sets ReadStoreElement::matePairId to ReadStoreElement::INVALID_ID.
 */

/*!
 * @var TId ReadStoreElement::matePairId;
 * @brief Refers to a mate pair in the @link FragmentStore::matePairStore @endlink or is INVALID_ID if the read is
 *        not paired.
 *
 * @var TId ReadStoreElement::INVALID_ID;
 * @brief Constant to represetn an invalid id.
 */

template <typename TSpec = void>
struct ReadStoreElement
{
    typedef typename Id<ReadStoreElement>::Type TId;

    static const TId INVALID_ID;

    TId matePairId;                // refers to the mate-pair, INVALID_ID if not part of a mate-pair

    ReadStoreElement() : matePairId(INVALID_ID) {}

    inline
    bool
    operator==(ReadStoreElement const & other)
    {
        return matePairId == other.matePairId;
    }
};

//////////////////////////////////////////////////////////////////////////////

template <typename TSpec>
const typename Id<ReadStoreElement<TSpec> >::Type
ReadStoreElement<TSpec>::INVALID_ID = MaxValue<typename Id<ReadStoreElement<TSpec> >::Type>::VALUE;

//////////////////////////////////////////////////////////////////////////////

}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
