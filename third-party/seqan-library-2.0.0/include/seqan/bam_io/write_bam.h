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
// Code for writing BAM.
// ==========================================================================

#ifndef INCLUDE_SEQAN_BAM_IO_WRITE_BAM_H_
#define INCLUDE_SEQAN_BAM_IO_WRITE_BAM_H_

namespace seqan {

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function writeRecord()                                             BamHeader
// ----------------------------------------------------------------------------

template <typename TTarget, typename TNameStore, typename TNameStoreCache, typename TStorageSpec>
void write(TTarget & target,
           BamHeader const & header,
           BamIOContext<TNameStore, TNameStoreCache, TStorageSpec> & context,
           Bam const & /*tag*/)
{
    write(target, "BAM\1");
    clear(context.buffer);

    // Create text of header.
    for (unsigned i = 0; i < length(header); ++i)
        write(context.buffer, header[i], context, Sam());

    // Note that we do not write out a null-character to terminate the header.  This would be valid by the SAM standard
    // but the samtools do not expect this and write out the '\0' when converting from BAM to SAM.
    // appendValue(context.buffer, '\0');

    // Write text header.
    appendRawPod(target, (__int32)length(context.buffer));
    write(target, context.buffer);

    // Write references.
    __int32 nRef = _max(length(contigNames(context)), length(contigLengths(context)));
    appendRawPod(target, nRef);

    for (__int32 i = 0; i < nRef; ++i)
    {
        if (i < (__int32)length(contigNames(context)))
        {
            appendRawPod(target, (__int32)(length(contigNames(context)[i]) + 1));
            write(target, contigNames(context)[i]);
        }
        else
        {
            appendRawPod(target, (__int32)1);
        }
        writeValue(target, '\0');
        __int32 lRef = 0;
        if (i < (__int32)length(contigLengths(context)))
            lRef = contigLengths(context)[i];
        appendRawPod(target, lRef);
    }
}

// ----------------------------------------------------------------------------
// Function writeRecord()                                    BamAlignmentRecord
// ----------------------------------------------------------------------------

static inline int _reg2Bin(uint32_t beg, uint32_t end)
{
    --end;
    if (beg >> 14 == end >> 14)
        return 4681 + (beg >> 14);

    if (beg >> 17 == end >> 17)
        return 585 + (beg >> 17);

    if (beg >> 20 == end >> 20)
        return 73 + (beg >> 20);

    if (beg >> 23 == end >> 23)
        return 9 + (beg >> 23);

    if (beg >> 26 == end >> 26)
        return 1 + (beg >> 26);

    return 0;
}

inline __uint32
updateLengths(BamAlignmentRecord const & record)
{
    // update internal lengths.
    record._l_qname = length(record.qName) + 1;
    record._n_cigar = length(record.cigar);
    record._l_qseq = length(record.seq);

    return sizeof(BamAlignmentRecordCore) + record._l_qname +
           record._n_cigar * 4 + (record._l_qseq + 1) / 2 + record._l_qseq +
           length(record.tags);
}


template <typename TTarget>
inline void
_writeBamRecord(TTarget & target,
                BamAlignmentRecord const & record,
                Bam const & /*tag*/)
{
    typedef typename Iterator<String<CigarElement<> > const, Standard>::Type SEQAN_RESTRICT TCigarIter;
    typedef typename Iterator<IupacString const, Standard>::Type SEQAN_RESTRICT             TSeqIter;
    typedef typename Iterator<CharString const, Standard>::Type SEQAN_RESTRICT              TQualIter;

    // bin_mq_nl
    unsigned l = 0;
    _getLengthInRef(l, record.cigar);
    record.bin =_reg2Bin(record.beginPos, record.beginPos + l);

    // Write fixed-size BamAlignmentRecordCore.
    appendRawPod(target, (BamAlignmentRecordCore &)record);

    // read_name
    write(target, record.qName);
    writeValue(target, '\0');

    // cigar
    static unsigned char const MAP[256] =
    {
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 7, 0, 0,
        0, 0, 0, 0, 2, 0, 0, 0, 5, 1, 0, 0, 0, 0, 3, 0,
        6, 0, 0, 4, 0, 0, 0, 0, 8, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
        0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
    };
    TCigarIter citEnd = end(record.cigar, Standard());
    for (TCigarIter cit = begin(record.cigar, Standard()); cit != citEnd; ++cit)
        appendRawPod(target, ((__uint32)cit->count << 4) | MAP[(unsigned char)cit->operation]);

    // seq
    TSeqIter sit = begin(record.seq, Standard());
    TSeqIter sitEnd = sit + (record._l_qseq & ~1);
    while (sit != sitEnd)
    {
        unsigned char x = (ordValue(getValue(sit++)) << 4);
        writeValue(target, x | ordValue(getValue(sit++)));
    }
    if (record._l_qseq & 1)
        writeValue(target, ordValue(getValue(sit++)) << 4);

    // qual
    SEQAN_ASSERT_LEQ(length(record.qual), length(record.seq));
    TQualIter qit = begin(record.qual, Standard());
    TQualIter qitEnd = end(record.qual, Standard());
    TQualIter qitVirtEnd = qit + record._l_qseq;
    while (qit != qitEnd)
        writeValue(target, *qit++ - '!');
    for (; qit != qitVirtEnd; ++qit)
        writeValue(target, '\xff');     // fill with zero qualities

    // tags
    write(target, record.tags);
}

template <typename TTarget>
inline void
_writeBamRecordWrapper(TTarget & target,
                       BamAlignmentRecord const & record,
                       Nothing & /* range */,
                       __uint32 size,
                       Bam const & tag)
{
    appendRawPod(target, size);
    _writeBamRecord(target, record, tag);
}

template <typename TTarget, typename TOValue>
inline void
_writeBamRecordWrapper(TTarget & target,
                       BamAlignmentRecord const & record,
                       Range<TOValue*> & range,
                       __uint32 size,
                       Bam const & tag)
{
    if (SEQAN_LIKELY(size + 4 <= length(range)))
    {
        appendRawPod(range.begin, size);
        _writeBamRecord(range.begin, record, tag);
        advanceChunk(target, size + 4);
    }
    else
    {
        appendRawPod(target, size);
        _writeBamRecord(target, record, tag);
    }
}

template <typename TTarget, typename TNameStore, typename TNameStoreCache, typename TStorageSpec>
void write(TTarget & target,
           BamAlignmentRecord const & record,
           BamIOContext<TNameStore, TNameStoreCache, TStorageSpec> & /* context */,
           Bam const & tag)
{
    // Update internal lengths
    __uint32 size = updateLengths(record);

    // Reserve chunk memory
    reserveChunk(target, 4 + size, Output());

    // Write length and record
    typename Chunk<TTarget>::Type ochunk;
    getChunk(ochunk, target, Output());
    _writeBamRecordWrapper(target, record, ochunk, size, tag);
}

}  // namespace seqan

#endif  // #ifndef INCLUDE_SEQAN_BAM_IO_WRITE_BAM_H_
