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

// TODO(holtgrew): Get out more than just strings...

#ifndef SEQAN_INCLUDE_SEQAN_VCF_IO_VCF_RECORD_H_
#define SEQAN_INCLUDE_SEQAN_VCF_IO_VCF_RECORD_H_

namespace seqan {

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// ----------------------------------------------------------------------------
// Class VcfRecord
// ----------------------------------------------------------------------------
//
/*!
 * @class VcfRecord
 * @implements FormattedFileRecordConcept
 * @headerfile <seqan/vcf_io.h>
 * @brief Information for one VCF record.
 *
 * @signature class VcfRecord;
 *
 * We store most information as strings and without structure since the VCF format's definition is quite loose.  We plan
 * to provide classes for structured access to these strings later.
 *
 * @section Remarks
 *
 * Although all positions in the VCF text format are 1-based, they are stored 0-based in the @link VcfRecord @endlink.
 * Positions in string members are stored verbatim in the @link VcfRecord @endlink's members, e.g. 1-based.
 *
 * Invalid qualities are stored as a float <tt>NaN</tt> (not a number).  To test a float quality <tt>q</tt> for being
 * <tt>NaN</tt>, test for <tt>q != q</tt>.  Only <tt>NaN</tt> has the property that <tt>NaN != NaN</tt>.
 */

/*!
 * @fn VcfRecord::MISSING_QUAL
 * @brief Return IEEE <tt>NaN</tt> float value.
 *
 * @signature float VcfRecord::MISSING_QUAL();
 *
 * @section Remarks
 *
 * This cannot be implemented portably as a constant in a header-only library.  In C++11, there is a function
 * <tt>std::nanf()</tt> that could be used instead.  For C++98, this is the best we can do.
 */

/*!
 * @fn VcfRecord::VcfRecord
 * @brief Default constructor.
 *
 * @signature VcfRecord::VcfRecord();
 */

/*!
 * @var VariableType VcfRecord::ref
 * @brief Bases in the reference (@link CharString @endlink).
 *
 * @var VariableType VcfRecord::genotypeInfos
 * @brief Genotype information, as in VCF file (@link StringSet @endlink<@link CharString @endlink>).
 *
 * @var VariableType VcfRecord::info
 * @brief Value of the INFO field, empty if "." in VCF file (@link CharString @endlink).
 *
 * @var VariableType VcfRecord::INVALID_REFID
 * @brief Static member as marker for invalid reference (<tt>__int32</tt>)
 *
 * @var VariableType VcfRecord::format
 * @brief Value of the VCF FORMAT field, empty if "." in VCF file (@link CharString @endlink).
 *
 * @var VariableType VcfRecord::alt
 * @brief Alternative bases in the variants, comma-separated if multiple (@link CharString @endlink).
 *
 * @var VariableType VcfRecord::qual
 * @brief Quality, <tt>NaN</tt> if invalid (<tt>float</tt>).
 *
 * @var VariableType VcfRecord::filter
 * @brief Value of FILTER field, empty if "." in VCF file (@link CharString @endlink).
 *
 * @var VariableType VcfRecord::rID
 * @brief Static member as marker for invalid reference (<tt>__int32</tt>)
 *
 * @var VariableType VcfRecord::beginPos
 * @brief Position of the VCF record (<tt>__int32</tt>).
 *
 * @var VariableType VcfRecord::INVALID_POS
 * @brief Static member as marker for invalid position (<tt>__int32</tt>)
 *
 * @var VariableType VcfRecord::id
 * @brief Textual identifier of the variant (@link CharString @endlink).
 */

class VcfRecord
{
public:
    // Constant for invalid reference id.
    static const __int32 INVALID_REFID = -1;
    // Constant for invalid position.
    static const __int32 INVALID_POS = -1;

    // Numeric id of the reference sequence.
    __int32 rID;
    // Position on the reference.
    __int32 beginPos;
    // Textual identifier of the variant.
    CharString id;
    // Bases in the reference.
    CharString ref;
    // Bases in the alternatives, comma-separated.
    CharString alt;
    // Quality
    float qual;
    // Value of FILTER field.
    CharString filter;
    // Value of INFO field.
    CharString info;
    // Value of FORMAT field.
    CharString format;
    // The genotype infos.
    StringSet<CharString> genotypeInfos;

    // Default constructor.
    VcfRecord() : rID(INVALID_REFID), beginPos(INVALID_POS), qual(MISSING_QUAL())
    {}

    // Actually this is IEEE NaN.
    static float MISSING_QUAL()
    {
        union {
            __uint32 u;
            float f;
        } tmp;
        tmp.u = 0x7F800001;
        return tmp.f;
    }
};

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function clear()
// ----------------------------------------------------------------------------

/*!
 * @fn VcfRecord#clear
 * @brief Clear a VcfRecord.
 *
 * @signature void clear(record);
 *
 * @param[in,out] record The VcfRecord to clear.
 */

inline void clear(VcfRecord & record)
{
    clear(record.id);
    clear(record.ref);
    clear(record.alt);
    clear(record.filter);
    clear(record.info);
    clear(record.format);
    clear(record.genotypeInfos);
}

}  // namespace seqan

#endif  // #ifndef SEQAN_INCLUDE_SEQAN_VCF_IO_VCF_RECORD_H_
