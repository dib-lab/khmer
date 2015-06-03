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
// Code for reading SAM.
// ==========================================================================

#ifndef INCLUDE_SEQAN_BAM_IO_READ_SAM_H_
#define INCLUDE_SEQAN_BAM_IO_READ_SAM_H_

namespace seqan {

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

struct Sam_;
typedef Tag<Sam_> Sam;


template <typename T>
struct FileExtensions<Sam, T>
{
    static char const * VALUE[1];    // default is one extension
};

template <typename T>
char const * FileExtensions<Sam, T>::VALUE[1] =
{
    ".sam"     // default output extension
};


template <typename T>
struct MagicHeader<Sam, T>
{
    static unsigned char const * VALUE;
};

template <typename T>
unsigned char const * MagicHeader<Sam, T>::VALUE = NULL;  // SAM has no magic header


enum SamTokenizeErrors_
{
    SAM_INVALID_RECORD = 2048
};

struct SamHeader_;
typedef Tag<SamHeader_> SamHeader;

struct SamAlignment_;
typedef Tag<SamAlignment_> SamAlignment;

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function nextIs()                                                  SamHeader
// ----------------------------------------------------------------------------

template <typename TForwardIter>
inline bool nextIs(TForwardIter & iter, SamHeader const & /*tag*/)
{
    if (atEnd(iter))
        return false;
    return value(iter) == '@';
}

//// ----------------------------------------------------------------------------
//// Function skipRecord()                                              SamHeader
//// ----------------------------------------------------------------------------
//
//template <typename TForwardIter, typename TPass>
//inline void skipRecord(TForwardIter & iter,
//                      SamHeader const & /*tag*/)
//{
//    skipOne(iter, EqualsChar<'@'>());
//    skipLine(iter);
//    return 0;
//}
//
//// ----------------------------------------------------------------------------
//// Function skipRecord()                                           SamAlignment
//// ----------------------------------------------------------------------------
//
//template <typename TForwardIter, typename TPass>
//inline void skipRecord(TForwardIter & iter,
//                      SamAlignment const & /*tag*/)
//{
//    skipOne(iter, EqualsChar<'@'>());
//    skipLine(iter);
//    return 0;
//}

// ----------------------------------------------------------------------------
// Function readRecord()                                        BamHeaderRecord
// ----------------------------------------------------------------------------

template <typename TForwardIter, typename TNameStore, typename TNameStoreCache, typename TStorageSpec>
inline void
readRecord(BamHeaderRecord & record,
           BamIOContext<TNameStore, TNameStoreCache, TStorageSpec> & context,
           TForwardIter & iter,
           Sam const & /*tag*/)
{
    clear(record);

    // Make sure the first character is '@'.
    skipOne(iter, EqualsChar<'@'>());

    // Read the header tag.
    char c1, c2;
    readOne(c1, iter);
    readOne(c2, iter);

    // Determine header type.
    if (c1 == 'H' && c2 == 'D')
        record.type = BAM_HEADER_FIRST;
    else if (c1 == 'S' && c2 == 'Q')
        record.type = BAM_HEADER_REFERENCE;
    else if (c1 == 'R' && c2 == 'G')
        record.type = BAM_HEADER_READ_GROUP;
    else if (c1 == 'P' && c2 == 'G')
        record.type = BAM_HEADER_PROGRAM;
    else if (c1 == 'C' && c2 == 'O')
        record.type = BAM_HEADER_COMMENT;
    else
        SEQAN_THROW(ParseError("Unknown SAM header type!"));

    CharString &buffer = context.buffer;

    if (record.type == BAM_HEADER_COMMENT)
    {
        skipOne(iter, IsTab());

        appendValue(record.tags, Pair<CharString>());

        clear(buffer);
        readLine(buffer, iter);
        assign(back(record.tags).i2, buffer, Exact());
    }
    else
    {
        // Read the rest of the line into the tag field of record.
        while (!atEnd(iter) && value(iter) == '\t')
        {
            skipOne(iter, IsTab());

            appendValue(record.tags, Pair<CharString>());

            clear(buffer);
            readUntil(buffer, iter, EqualsChar<':'>());
            assign(back(record.tags).i1, buffer, Exact());

            skipOne(iter, EqualsChar<':'>());

            clear(buffer);
            readUntil(buffer, iter, OrFunctor<IsTab, IsNewline>());
            assign(back(record.tags).i2, buffer, Exact());
        }
        // Skip remaining line break
        skipLine(iter);
    }
}

// ----------------------------------------------------------------------------
// Function readRecord()                                              BamHeader
// ----------------------------------------------------------------------------

template <typename TForwardIter, typename TNameStore, typename TNameStoreCache, typename TStorageSpec>
inline void
readHeader(BamHeader & header,
           BamIOContext<TNameStore, TNameStoreCache, TStorageSpec> & context,
           TForwardIter & iter,
           Sam const & tag)
{
    BamHeaderRecord record;
    while (nextIs(iter, SamHeader()))
    {
        clear(record);
        readRecord(record, context, iter, tag);
        appendValue(header, record);

        // Get sequence information from @SQ header.
        if (record.type == BAM_HEADER_REFERENCE)
        {
            CharString name;
            unsigned lRef = 0;
            for (unsigned i = 0; i < length(record.tags); ++i)
            {
                if (record.tags[i].i1 == "SN")
                    name = record.tags[i].i2;
                else if (record.tags[i].i1 == "LN")
                    lexicalCast(lRef, record.tags[i].i2);
            }

            // Add entry to name store and sequenceInfos if necessary.
            size_t globalRefId = nameToId(contigNamesCache(context), name);
            if (length(contigLengths(context)) <= globalRefId)
                resize(contigLengths(context), globalRefId + 1, 0);
            contigLengths(context)[globalRefId] = lRef;
        }
    }
}

// ----------------------------------------------------------------------------
// Function _readBamRecord()
// ----------------------------------------------------------------------------

template <typename TBuffer, typename TForwardIter>
inline void
_readBamRecord(TBuffer & rawRecord, TForwardIter & iter, Sam)
{
    clear(rawRecord);
    readLine(rawRecord, iter);
}

// ----------------------------------------------------------------------------
// Function readRecord()                                     BamAlignmentRecord
// ----------------------------------------------------------------------------

template <typename TForwardIter, typename TNameStore, typename TNameStoreCache, typename TStorageSpec>
inline void
readRecord(BamAlignmentRecord & record,
           BamIOContext<TNameStore, TNameStoreCache, TStorageSpec> & context,
           TForwardIter & iter,
           Sam const & /*tag*/)
{
    // fail, if we read "@" (did you miss to call readRecord(header, bamFile) first?)
    if (nextIs(iter, SamHeader()))
        SEQAN_THROW(ParseError("Unexpected SAM header encountered."));

    OrFunctor<IsTab, AssertFunctor<NotFunctor<IsNewline>, ParseError, Sam> > nextEntry;

    clear(record);
    CharString &buffer = context.buffer;

    // QNAME
    readUntil(record.qName, iter, nextEntry);
    skipOne(iter, IsTab());

    // FLAG
    // TODO(holtgrew): Interpret hex and char as c-samtools -X does?
    clear(buffer);
    readUntil(buffer, iter, nextEntry);
    record.flag = lexicalCast<__uint16>(buffer);
    skipOne(iter, IsTab());

    // RNAME
    clear(buffer);
    readUntil(buffer, iter, nextEntry);
    if (buffer == "*")
        record.rID = BamAlignmentRecord::INVALID_REFID;
    else
        record.rID = nameToId(contigNamesCache(context), buffer);
    skipOne(iter, IsTab());

    // POS
    clear(buffer);
    SEQAN_ASSERT_EQ((__int32)0 - 1, (__int32)BamAlignmentRecord::INVALID_POS);
    readUntil(buffer, iter, nextEntry);
    record.beginPos = (__int32)lexicalCast<__uint32>(buffer) - 1;
    skipOne(iter, IsTab());

    // MAPQ
    clear(buffer);
    if (value(iter) == '*')
    {
        record.mapQ = 255;
        skipOne(iter);
    }
    else
    {
        readUntil(buffer, iter, nextEntry);
        record.mapQ = lexicalCast<__uint16>(buffer);
    }
    skipOne(iter, IsTab());

    // CIGAR
    CigarElement<> element;
    if (value(iter) == '*')
        skipOne(iter);
    else
    {
        do
        {
            clear(buffer);
            readUntil(buffer, iter, OrFunctor<IsAlpha, AssertFunctor<NotFunctor<IsNewline>, ParseError, Sam> >());
            element.count = lexicalCast<__uint32>(buffer);
            element.operation = value(iter);
            skipOne(iter);
            appendValue(record.cigar, element);
        } while (value(iter) != '\t');
    }
    skipOne(iter, IsTab());

    // RNEXT
    clear(buffer);
    readUntil(buffer, iter, nextEntry);
    if (buffer == "*")
        record.rNextId = BamAlignmentRecord::INVALID_REFID;
    else if (buffer == "=")
        record.rNextId = record.rID;
    else
        record.rNextId = nameToId(contigNamesCache(context), buffer);
    skipOne(iter, IsTab());

    // PNEXT
    if (value(iter) == '*')
    {
        record.pNext = BamAlignmentRecord::INVALID_POS;
        skipOne(iter);
    }
    else
    {
        clear(buffer);
        readUntil(buffer, iter, nextEntry);
        record.pNext = (__int32)lexicalCast<__uint32>(buffer) - 1;
    }
    skipOne(iter, IsTab());

    // TLEN
    if (value(iter) == '*')
    {
        record.tLen = MaxValue<__int32>::VALUE;
        skipOne(iter);
    }
    else
    {
        clear(buffer);
        readUntil(buffer, iter, nextEntry);
        record.tLen = lexicalCast<__int32>(buffer);
    }
    skipOne(iter, IsTab());

    // SEQ
    readUntil(record.seq, iter, nextEntry);
    // Handle case of missing sequence:  Clear seq string as documented.
    if (record.seq == "*")
        clear(record.seq);
    skipOne(iter, IsTab());

    // QUAL
    readUntil(record.qual, iter, OrFunctor<IsTab, IsNewline>());

    // Handle case of missing quality:  Clear qual string as documented.
    if (record.qual == "*")
        clear(record.qual);

    // The following list of tags is optional.  A line break or EOF could also follow.
    if (atEnd(iter))
        return;
    if (value(iter) != '\t')
    {
        skipLine(iter);
        return;
    }
    skipOne(iter, IsTab());

    // TAGS
    clear(buffer);
    readLine(buffer, iter);
    appendTagsSamToBam(record.tags, buffer);
}

}  // namespace seqan

#endif  // #ifndef INCLUDE_SEQAN_BAM_IO_READ_SAM_H_

