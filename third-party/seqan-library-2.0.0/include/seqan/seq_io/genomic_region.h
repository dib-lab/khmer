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
// Author: David Weese <david.weese@fu-berlin.de>
// ==========================================================================
// The GenomicRegion class represents a region on one chromosome/contig in a
// genome, e.g. chr1, chr1:15,000, chr1:100,000-200,000.
// ==========================================================================

#ifndef INCLUDE_SEQAN_SEQ_IO_GENOMIC_REGION_H_
#define INCLUDE_SEQAN_SEQ_IO_GENOMIC_REGION_H_

#include <seqan/sequence.h>
#include <seqan/stream.h>  // for tokenization

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

struct GenomicRegion;
inline void parse(GenomicRegion & region, CharString const & regionString);

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// ---------------------------------------------------------------------------
// Class GenomicRegion
// ---------------------------------------------------------------------------

/*!
 * @class GenomicRegion
 * @headerfile <seqan/seq_io.h>
 * @brief Store information about a genomic region.
 *
 * @signature class GenomicRegion;
 *
 * A genomic region is a range on a chromosome.  The chromosome is identified by its name (as text in @link
 * GenomicRegion::seqName @endlink, optionally also as an <tt>integer</tt> in @link GenomicRegion::rID @endlink).  The
 * range is stored as a half-open interval [@link GenomicRegion::beginPos @endlink, @link GenomicRegion::endPos @endlink).
 * If @link GenomicRegion::beginPos @endlink is set to <tt>-1</tt> then the range spans the whole chromosome.  If @link
 * GenomicRegion::beginPos @endlink is set to a value <tt>&gt;= 0</tt> and @link GenomicRegion::endPos @endlink is set ot
 * <tt>-1</tt>, then the chromosome is selected from @link GenomicRegion::beginPos @endlink to the end.
 *
 * Examples for genomic regions are <tt>chr1</tt>, <tt>chr1:1,000</tt>, <tt>chr1:1,000-2,000</tt>.
 *
 * The textual description of a genomic region has one of the formats <tt>NAME</tt>, <tt>NAME:START</tt>,
 * <tt>NAME:START-END</tt>.  The positions in the textual representation <tt>START</tt> and <tt>END</tt> are one-based.
 * However, the representation in the members of @link GenomicRegion @endlink is zero-based.
 *
 * @section Examples
 *
 * Construct a @link GenomicRegion @endlink object and fill it from different region strings.
 *
 * @code{.cpp}
 * seqan::GenomicRegion genomicRegion;
 *
 * parse(genomicRegion, "chr1");
 * // genomicRegion.seqName == "chr1"
 * // genomicRegion.rID == -1, genomicRegion.beginPos == -1, genomicRegion.beginPos == -1
 *
 * parse(genomicRegion, "chr1:1000");
 * // genomicRegion.seqName == "chr1"
 * // genomicRegion.beginPos == 999
 * // genomicRegion.rID == -1, genomicRegion.beginPos == -1
 *
 * parse(genomicRegion, "chr1:1000-2000");
 * // genomicRegion.seqName == "chr1"
 * // genomicRegion.beginPos == 999
 * // genomicRegion.beginPos == 2000
 * // genomicRegion.rID == -1
 * @endcode
 *
 */

struct GenomicRegion
{
    /*!
     * @var __int32 GenomicRegion::INVALID_POS;
     * @brief Constant for marking a position as invalid (static-const member).
     */
    static __int32 const INVALID_POS = -1;

    /*!
     * @var __int32 GenomicRegion::INVALID_ID;
     * @brief Constant for marking a position as invalid (static-const member).
     */
    static __int32 const INVALID_ID = -1;

    /*!
     * @var CharString GenomicRegion::seqName;
     * @brief Name of the sequence the region lies on, default is the empty string.
     */
    CharString seqName;

    /*!
     * @var __int32 GenomicRegion::rID;
     * @brief An optional field storing an integer.  Default is <tt>-1</tt>.
     */
    __int32 rID;

    /*!
     * @var __int32 GenomicRegion::beginPos;
     * @brief Begin position of the range on the chromosome.  Default is <tt>-1</tt>.
     */
    __int32 beginPos;

    /*!
     * @var __int32 GenomicRegion::endPos;
     * @brief End position of the range on the chromosome.  Default is <tt>-1</tt>.
     */
    __int32 endPos;

    /*!
     * @fn GenomicRegion::GenomicRegion
     * @brief Constructor.
     *
     * @signature GenomicRegion::GenomicRegion();
     * @signature GenomicRegion::GenomicRegion(str);
     *
     * @param[in] str The string to parse region from. Types: CharString
     *
     * The default constructor sets all integer members to <tt>INVALID_POS</tt>, the <tt>seqName</tt> member is left empty.
     */

    GenomicRegion() : rID(INVALID_ID), beginPos(INVALID_POS), endPos(INVALID_POS)
    {}

    GenomicRegion(CharString const & str) : rID(INVALID_ID), beginPos(INVALID_POS), endPos(INVALID_POS)
    {
        parse(*this, str);
    }

    /*!
     * @fn GenomicRegion::toString
     * @brief Write string representation of interval to out.
     *
     * @signature void GenomicRegion::toString(out);
     *
     * @param[in,out] out Target to write textual interval description to.
     */
    template <typename TString>
    void toString(TString & out) const
    {
        clear(out);
        append(out, seqName);
        appendValue(out, ':');
        if (beginPos + 1 == endPos)
        {
            appendNumber(out, beginPos + 1);
        }
        else
        {
            appendNumber(out, beginPos + 1);
            appendValue(out, '-');
            appendNumber(out, endPos);
        }
    }
};

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

// ---------------------------------------------------------------------------
// Function clear()
// ---------------------------------------------------------------------------

/*!
 * @fn GenomicRegion#clear
 * @brief Reset a GenomicRegion object to the same state after default construction.
 *
 * @signature void clear(genomicRegion);
 *
 * @param[in,out] genomicRegion The @link GenomicRegion @endlink object to reset.  Types: GenomicRegion
 */

inline void
clear(GenomicRegion & region)
{
    clear(region.seqName);
    region.rID = region.INVALID_ID;
    region.beginPos = region.INVALID_POS;
    region.endPos = region.INVALID_POS;
}

// ---------------------------------------------------------------------------
// Function parse()
// ---------------------------------------------------------------------------

/*!
 * @fn GenomicRegion#parse
 * @brief Parse genomic region string store results in @link GenomicRegion @endlink.
 *
 * @signature bool parse(genomicRegion, regionString);
 *
 * @param[in]  regionString  The region string to prse.  Types: @link CharString @endlink.
 * @param[out] genomicRegion The @link GenomicRegion @endlink object to write the results to.
 *                           Types: GenomicRegion
 *
 * @return bool true indicates successful parsing.
 */

// Parse regionString and write to region.  region.rID will not be set but
// region.seqName will be.  Return true on success.

inline void
parse(GenomicRegion & region, CharString const & regionString)
{
    DirectionIterator<CharString const, Input>::Type reader = directionIterator(regionString, Input());

    // Parse out sequence name.
    clear(region.seqName);
    readUntil(region.seqName, reader, EqualsChar<':'>());
    if (atEnd(reader))
        return;

    skipOne(reader, EqualsChar<':'>());     // Skip ':'.

    // Parse out begin position.
    CharString buffer;
    readUntil(buffer, reader, EqualsChar<'-'>(), EqualsChar<','>());
    lexicalCastWithException(region.beginPos, buffer);

    if (region.beginPos < 1)
        SEQAN_THROW(ParseError("GenomicRegion: Begin postition less than 1"));

    region.beginPos--;                      // Adjust to 0-based.
    if (atEnd(reader))  // just one position
    {
        region.endPos = region.INVALID_POS;
        return;
    }

    skipOne(reader, EqualsChar<'-'>());     // Skip '-'.

    // Parse out end position.
    clear(buffer);
    readUntil(buffer, reader, False(), EqualsChar<','>());
    lexicalCastWithException(region.endPos, buffer);

    if (region.endPos < 1)
        SEQAN_THROW(ParseError("GenomicRegion: End postition less than 1"));
}

}  // namespace seqan

#endif  // #ifndef INCLUDE_SEQAN_SEQ_IO_GENOMIC_REGION_H_
