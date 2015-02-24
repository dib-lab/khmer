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

// TODO(holtgrew): Parse more than just the key/value pair.

#ifndef SEQAN_INCLUDE_SEQAN_VCF_IO_VCF_HEADER_RECORD_H_
#define SEQAN_INCLUDE_SEQAN_VCF_IO_VCF_HEADER_RECORD_H_

namespace seqan {

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// ----------------------------------------------------------------------------
// Class VcfHeaderRecord
// ----------------------------------------------------------------------------

/*!
 * @class VcfHeaderRecord
 * @headerfile <seqan/vcf_io.h>
 * @brief Store key/value pair for VCF header records.
 *
 * @signature class VcfHeaderRecord;
 *
 * @var CharString VcfHeaderRecord::key;
 * @brief Key of the header record.
 *
 * @var CharString VcfHeaderRecord::value;
 * @brief Value of the header record.
 */

/*!
 * @fn VcfHeaderRecord::VcfHeaderRecord
 * @brief Constructor
 *
 * @signature VcfHeaderRecord::VcfHeaderRecord();
 * @signature VcfHeaderRecord::VcfHeaderRecord(key, value);
 *
 * @param[in] key   Key of the header record, @link CharString @endlink.
 * @param[in] value Key of the header record, @link CharString @endlink.
 */

/*!
 * @fn VcfHeaderRecord#clear
 *
 * @brief Clear a VcfHeaderRecord.
 * @signature void clear(record);
 *
 * @param[in,out] record The VcfHeaderRecord to clear.
 */

class VcfHeaderRecord
{
public:
    // Record's key.
    CharString key;
    // Record's value.
    CharString value;

    // Default constructor.
    VcfHeaderRecord()
    {}

    // Construct directly with key/value.
    VcfHeaderRecord(CharString const & key, CharString const & value) :
            key(key), value(value)
    {}
};

// ============================================================================
// Functions
// ============================================================================

inline void clear(VcfHeaderRecord & record)
{
    clear(record.key);
    clear(record.value);
}

}  // namespace seqan

#endif  // #ifndef SEQAN_INCLUDE_SEQAN_VCF_IO_VCF_HEADER_RECORD_H_
