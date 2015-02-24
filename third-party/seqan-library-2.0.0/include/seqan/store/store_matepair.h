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

#ifndef SEQAN_HEADER_STORE_MATEPAIR_H
#define SEQAN_HEADER_STORE_MATEPAIR_H

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////
// Mate Store
//////////////////////////////////////////////////////////////////////////////

/*!
 * @class MatePairStoreElement
 * @brief Represents a mate pair.
 *
 * @signature template <[typename TSpec]>
 *            struct MatePairStoreElement;
 *
 * @tparam TSpec The specializing type.
 *
 * A mate pair consists of two reads sequences from opposite ends of the same fragment.  The insert size of a mate pair
 * is the inferred size of the fragment.
 *
 *
 * @fn MatePairStoreElement::MatePairStoreElement
 * @brief Constructor.
 *
 * @signature MatePairStoreElement::MatePairStoreElement();
 *
 * Initialize all member values to INVALID_ID.
 */

/*!
 * @var TId MatePairStoreElement::INVALID_ID;
 * @brief Constant to represent an invalid id.
 *
 * @var TId MatePairStoreElement::libId;
 * @brief Refers to a library in the @link FragmentStore::libraryStore @endlink or is INVALID_ID if the mate pair
 *        has no library.
 *
 * @var TId[2] MatePairStoreElement::readId;
 * @brief Refers to two paried reads in @link FragmentStore::readStore @endlink or is INVALID_ID values.
 */

template <typename TSpec = void>
struct MatePairStoreElement
{
    typedef typename Id<MatePairStoreElement>::Type TId;

    static const TId INVALID_ID;

    TId        readId[2];    // refers to the two reads of a mate-pair, INVALID_ID if this is a singleton fragment (e.g. in afg: reads refer to fragments (mate pairs) and these refer to libraries, singletons refer to an empty fragment)
    TId        libId;

    MatePairStoreElement() : libId(INVALID_ID)
    {
        readId[0] = INVALID_ID;
        readId[1] = INVALID_ID;
    }

    inline bool operator==(MatePairStoreElement const & other) const
    {
        return readId[0] == other.readId[0] &&
                readId[1] == other.readId[1] &&
                libId == other.libId;
    }
};

//////////////////////////////////////////////////////////////////////////////

/*!
 * @mfn MatePairStoreElement#Id
 * @headerfile <seqan/store.h>
 * @brief Returns the id type to use for <tt>MatePairStoreElement</tt>.
 *
 * @signature Id<TMatePairStoreElement>::Type;
 *
 * @tparam TMatePairStoreElement The MatePairStoreElement specialization to get the id type for.
 *
 * @return Type The resulting id type.
 */

template <typename TSpec>
const typename Id<MatePairStoreElement<TSpec> >::Type
MatePairStoreElement<TSpec>::INVALID_ID = MaxValue<typename Id<MatePairStoreElement<TSpec> >::Type>::VALUE;

//////////////////////////////////////////////////////////////////////////////

}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
