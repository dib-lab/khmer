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
// Author: Jochen Singer <jochen.singer@fu-berlin.de>
// Author: Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>
// ==========================================================================

#ifndef SEQAN_INCLUDE_SEQAN_UCSC_RECORD_H_
#define SEQAN_INCLUDE_SEQAN_UCSC_RECORD_H_

namespace seqan {

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// ----------------------------------------------------------------------------
// Class UcscRecord
// ----------------------------------------------------------------------------

/*!
 * @class UcscRecord
 * @implements FormattedFileRecordConcept
 * @headerfile <seqan/ucsc_io.h>
 * @brief Represent the information for one UCSC gene annotation record.
 *
 * @signature class UcscRecord;
 *
 * The UCSC genome browser server allows the download of a set of the known genes and isoforms used in the browser.
 * These files can be downloaded from the UCSC FTP's <tt>goldenPath</tt> folder, e.g. the one for <a
 * href="http://hgdownload.cse.ucsc.edu/goldenpath/hg19/">hg19</a>.
 *
 * To load the annotations, you need to download both the <tt>knownGenes.txt.gz</tt> and the
 * <tt>knownIsoforms.txt.gz</tt> files from the UCSC <tt>goldenPath</tt> database.  These files can then be read record
 * by record into a <tt>UcscRecord</tt>.  This record type is capable of representing records from either file type.
 */

class UcscRecord
{
public:

    /*!
     * @var CharString UcscRecord::transName
     * @brief Name of the transcript.
     */
    CharString transName;

    /*!
     * @var CharString UcscRecord::contigName
     * @brief Name of the contig of the genomic location.
     */
    CharString contigName;

    /*!
     * @var __int32 UcscRecord::cdsBegin
     * @brief Start of the coding region (<tt>0</tt>-based position, defaults to <tt>-1</tt>).
     */
    __int32 cdsBegin;

    /*!
     * @var __int32 UcscRecord::cdsEnd
     * @brief End of the coding region, defaults to <tt>-1</tt>.
     */
    __int32 cdsEnd;

    /*!
     * @var CharString UcscRecord::exonBegin
     * @brief Start of the exon (<tt>0</tt>-based position, defaults to <tt>-1</tt>).
     */
    String<__int32> exonBegin;

    /*!
     * @var CharString UcscRecord::exonEnds
     * @brief End of the exon, defaults to <tt>-1</tt>.
     */
    String<__int32> exonEnds;

    /*!
     * @var CharString UcscRecord::proteinName
     * @brief Name of the coded protein.
     */
    CharString proteinName;

    /*!
     * @var __uint32 UcscRecord::annotationBeginPos
     * @brief Start of the annotation (<tt>0</tt>-based, defaults to <tt>-1</tt>).
     */
    __uint32 annotationBeginPos;

    /*!
     * @var CharString UcscRecord::annotationEndPos
     * @brief End position of the annotation, defaults to <tt>-1</tt>.
     */
    __uint32 annotationEndPos;

    UcscRecord() : cdsBegin(0), cdsEnd(0), annotationBeginPos(0), annotationEndPos(0)
    {}
};

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function clear()
// ----------------------------------------------------------------------------

/*!
 * @fn UcscRecord#clear
 * @brief Clear UcscRecord.
 *
 * @signature void clear(record);
 *
 * @param record The UcscRecord to clear.
 *
 * Clears all strings and resets it to default initialization state.
 */

inline void clear(UcscRecord & record)
{
    clear(record.transName);
    clear(record.contigName);
    record.cdsBegin = -1;
    record.cdsEnd = -1;
    clear(record.exonBegin);
    clear(record.exonEnds);
    clear(record.proteinName);

    record.annotationBeginPos = -1;
    record.annotationEndPos = -1;
}

}  // namespace seqan

#endif  // #ifndef SEQAN_INCLUDE_SEQAN_BAM_IO_BAM_RECORD_H_
