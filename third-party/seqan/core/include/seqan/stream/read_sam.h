// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2013, Knut Reinert, FU Berlin
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
// Record I/O for SAM.  The FragmentStore uses this API for loading whole
// SAM documents.  Only Single-Pass I/O is implemented.
//
// Usage:
//
// // Open file and record reader.
// FILE * f = fopen("example.sam", "rb");
// RecordReader<FILE *, SinglePass<TSpec> > reader(f);
//
// // Read all headers.
// SamHeader header;
// String<SamHeader> headers;
// while (!atEnd(reader) && value(reader) == '@') {
//     readRecord(header, reader, Sam());
//     appendValue(headers, header);
// }
//
// // Read all alignments.
// SamAlignment alignment;
// String<SamAlignment> alignments;
// while (!atEnd(reader)) {
//     readRecord(alignment, reader, Sam());
//     appendValue(alignments, alignment);
// }
// ==========================================================================

// TODO(holtgrew): Delete this, replaced by bam_io module.

#ifndef SEQAN_STREAM_READ_SAM_H_
#define SEQAN_STREAM_READ_SAM_H_

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

struct SamHeader
{
    typedef String<char, Array<2> > TTag;

    TTag recordType;
    String<Pair<TTag, CharString> > dataFields;
};

enum SamFieldType
{
    SAM_CHAR,
    SAM_INT_32,
    SAM_FLOAT,
    SAM_STRING,
    SAM_HEX_STRING
};

struct SamAlignment
{
    typedef String<char, Array<2> > TTag;

    CharString queryName;      //  1 QNAME
    __int32 flag;              //  2 FLAG
    CharString referenceName;  //  3 RNAME
    __int32 pos;               //  4 POS
    __int32 mappingQuality;    //  5 MAPQ
    CharString cigarString;    //  6 CIGAR
    CharString referenceNext;  //  7 RNEXT
    __int32 positionNext;      //  8 PNEXT
    __int32 templateLength;    //  9 TLEN
    Dna5String sequence;       // 10 SEQ
    CharString qualities;      // 11 QUAL
    
    String<Triple<TTag, SamFieldType, CharString> > optionalFields;
};

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function readRecord(SamHeader, );  Single-Pass.
// ----------------------------------------------------------------------------

inline
bool isLineBreak(char c)
{
    return c == '\r' || c == '\n';
}

template <typename TIdString, typename TSeqString, typename TFile, typename TSpec>
int readRecord(SamHeader & header,
               RecordReader<TFile, SinglePass<TSpec> > & reader,
               Sam const & /*tag*/)
{
    // Read record type.
    SEQAN_ASSERT_NOT(atEnd(reader));
    SEQAN_ASSERT_EQ(value(reader), '@');
    goNext(reader);
    SEQAN_ASSERT_NOT(atEnd(reader));
    header.queryName[0] = value(reader);
    goNext(reader);
    SEQAN_ASSERT_NOT(atEnd(reader));
    header.queryName[1] = value(reader);
    goNext(reader);
    SEQAN_ASSERT_NOT(atEnd(reader));

    // Read all fields.
    while (!atEnd(reader) && !isLineBreak(value(reader))) {
        SEQAN_ASSERT_EQ(value(reader), '\t');
        goNext(reader);

        resize(header.dataFields, length(header.dataFields) + 1);
        back(header.dataFields).i1[0] = value(reader);
        goNext(reader);
        back(header.dataFields).i1[1] = value(reader);
        goNext(reader);
        SEQAN_ASSERT_EQ(value(reader), ':');
        goNext(reader);
        for (; value(reader) != '\t' && !isLineBreak(value(reader)); goNext(reader))
            appendValue(back(header.dataFields).i2, value(reader));
    }
    // Skip line break.
    if (!atEnd(reader) && isLineBreak(value(reader))) {
        goNext(reader);
        if (value(reader) == '\n')
            goNext(reader);
    }
    return 0;
}

// ----------------------------------------------------------------------------
// Function readRecord(SamHeader, );  Single-Pass.
// ----------------------------------------------------------------------------

template <typename TIdString, typename TSeqString, typename TFile, typename TSpec>
int readRecord(SamAlignment & record,
               RecordReader<TFile, SinglePass<TSpec> > & reader,
               Sam const & /*tag*/)
{
    SEQAN_ASSERT_NEQ(value(reader), '@');

    // Reading integers.
    CharString buffer;
    tokenizeTo(buffer, reader, '\t');
    record.beginPos = ::seqan::lexical_cast<__int32>(buffer); // We should probably steal/port the boost code, with enforcing english locale.
    goNext(reader);
    
    // ... etc.
}

}  // namespace seqan

#endif  // #ifndef SEQAN_STREAM_READ_SAM_H_
