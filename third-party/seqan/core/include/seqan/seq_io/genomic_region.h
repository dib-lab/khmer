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
// The GenomicRegion class represents a region on one chromosome/contig in a
// genome, e.g. chr1, chr1:15,000, chr1:100,000-200,000.
// ==========================================================================

#ifndef CORE_INCLUDE_SEQAN_SEQ_IO_GENOMIC_REGION_H_
#define CORE_INCLUDE_SEQAN_SEQ_IO_GENOMIC_REGION_H_

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

class GenomicRegion;
inline bool parse(GenomicRegion & region, CharString const & regionString);

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// ---------------------------------------------------------------------------
// Class GenomicRegion
// ---------------------------------------------------------------------------

/**
.Class.GenomicRegion
..cat:Input/Output
..signature:GenomicRegion
..summary:Store information about a genomic region.
..description.text:
A genomic region is a range on a chromosome.
The chromosome is identified by its name (as text in @Memvar.GenomicRegion#seqName@, optionally also as an $integer$ in @Memvar.GenomicRegion#seqId@).
The range is stored as a half-open interval [@Memvar.GenomicRegion#beginPos@, @Memvar.GenomicRegion#endPos@).
If @Memvar.GenomicRegion#beginPos@ is set to $-1$ then the range spans the whole chromosome.
If @Memvar.GenomicRegion#beginPos@ is set to a value $>= 0$ and @Memvar.GenomicRegion#endPos@ is set ot $-1$, then the chromosome is selected from @Memvar.GenomicRegion#beginPos@ to the end.
..description.text:
Examples for genomic regions are $chr1$, $chr1:1,000$, $chr1:1,000-2,000$.
..description.text:
The textual description of a genomic region has one of the formats $NAME$, $NAME:START$, $NAME:START-END$.
The positions in the textual representation $START$ and $END$ are one-based.
However, the representation in the members of @Class.GenomicRegion@ is zero-based.
..example.text:Construct a @Class.GenomicRegion@ object and fill it from different region strings.
..example.code:
seqan::GenomicRegion genomicRegion;

parse(genomicRegion, "chr1");
// genomicRegion.seqName == "chr1"
// genomicRegion.seqId == -1, genomicRegion.beginPos == -1, genomicRegion.beginPos == -1

parse(genomicRegion, "chr1:1000");
// genomicRegion.seqName == "chr1"
// genomicRegion.beginPos == 999
// genomicRegion.seqId == -1, genomicRegion.beginPos == -1

parse(genomicRegion, "chr1:1000-2000");
// genomicRegion.seqName == "chr1"
// genomicRegion.beginPos == 999
// genomicRegion.beginPos == 2000
// genomicRegion.seqId == -1


..include:seqan/seq_io.h

.Memfunc.GenomicRegion#GenomicRegion
..class:Class.GenomicRegion
..summary:Constructor.
..description:
The default constructor sets all integer members to $-1$, the $seqName$ member is left empty.
..signature:GenomicRegion()
..signature:GenomicRegion(str)
..param.str:The string to parse region from.
...type:Shortcut.CharString

.Memvar.GenomicRegion#seqName
..class:Class.GenomicRegion
..summary:Name of the sequence the region lies on, default is the empty string.
..type:Shortcut.CharString

.Memvar.GenomicRegion#seqId
..class:Class.GenomicRegion
..summary:An optional field storing an integer. Default is $-1$.
..type:nolink:$__int32$

.Memvar.GenomicRegion#beginPos
..class:Class.GenomicRegion
..summary:Begin position of the range on the chromosome. Default is $-1$.
..type:nolink:$__int32$

.Memvar.GenomicRegion#endPos
..class:Class.GenomicRegion
..summary:End position of the range on the chromosome. Default is $-1$.
..type:nolink:$__int32$
*/

class GenomicRegion
{
public:
    // Name of sequence.
    CharString seqName;
    // Index of sequence in FASTA file.  -1 if not set.
    __int32 seqId;
    // 0-based begin position.  -1 if not set.
    __int32 beginPos;
    // 0-based, C-style end position.  -1 if not set.
    __int32 endPos;

    GenomicRegion(CharString const & str) :
            seqId(-1), beginPos(-1), endPos(-1)
    {
        parse(*this, str);
    }

    GenomicRegion() :
        seqId(-1), beginPos(-1), endPos(-1)
    {}
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

/**
.Function.GenomicRegion#clear
..cat:Input/Output
..class:Class.GenomicRegion
..summary:Reset a @Class.GenomicRegion@ object to the same state after default construction.
..signature:reset(genomicRegion)
..param.genomicRegion:The @Class.GenomicRegion@ object to reset.
...type:Class.GenomicRegion
..returns:$void$, where $true$ indicates sucess
..include:seqan/seq_io.h
*/

inline void clear(GenomicRegion & region)
{
    clear(region.seqName);
    region.seqId = -1;
    region.beginPos = -1;
    region.endPos = -1;
}

// ---------------------------------------------------------------------------
// Function parse()
// ---------------------------------------------------------------------------

/**
.Function.GenomicRegion#parse
..cat:Input/Output
..class:Class.GenomicRegion
..summary:Parse genomic region string store results in @Class.GenomicRegion@.
..signature:parse(genomicRegion, regionString)
..param.genomicRegion:The @Class.GenomicRegion@ object to write the results to.
...type:Class.GenomicRegion
..param.regionString:The region string to prse.
...type:Shortcut.CharString
..returns:$bool$, where $true$ indicates sucess
..example.text:See the example for parsing in the @Class.GenomicRegion@.
..include:seqan/seq_io.h
*/

// Parse regionString and write to region.  region.seqId will not be set but
// region.seqName will be.  Return true on success.

inline bool parse(GenomicRegion & region, CharString const & regionString)
{
    Stream<CharArray<char const *> > stream(begin(regionString, Standard()),
                                            end(regionString, Standard()));
    RecordReader<Stream<CharArray<char const *> >, SinglePass<> > reader(stream);

    // Parse out sequence name.
    CharString buffer;
    int res = readUntilChar(buffer, reader, ':');
    if (res != 0 && res != EOF_BEFORE_SUCCESS)
        return 1;  // Parse error.

    region.seqName = buffer;
    if (atEnd(reader))
        return true;  // Done after parsing the sequence name.

    goNext(reader);  // Skip ':'.

    // Parse out begin position.
    clear(buffer);
    while (!atEnd(reader) && value(reader) != '-')
    {
        if (!isdigit(value(reader)) && value(reader) != ',')
            return false;  // Error parsing.

        if (isdigit(value(reader)))
            appendValue(buffer, value(reader));
        goNext(reader);
    }
    if (empty(buffer))
        return false;

    if (!lexicalCast2(region.beginPos, buffer))
        return false;

    if (region.beginPos <= 0)
        return false;

    region.beginPos -= 1;  // Adjust to 0-based.
    if (atEnd(reader))
        return true;

    goNext(reader);  // Skip '-'.

    // Parse out end position.
    clear(buffer);
    while (!atEnd(reader))
    {
        if (!isdigit(value(reader)) && value(reader) != ',')
            return false;  // Error parsing.

        if (isdigit(value(reader)))
            appendValue(buffer, value(reader));
        goNext(reader);
    }
    if (empty(buffer))
        return false;

    if (!lexicalCast2(region.endPos, buffer))
        return false;

    if (region.endPos <= 0)
        return false;

    return atEnd(reader);
}

}  // namespace seqan

#endif  // #ifndef CORE_INCLUDE_SEQAN_SEQ_IO_GENOMIC_REGION_H_
