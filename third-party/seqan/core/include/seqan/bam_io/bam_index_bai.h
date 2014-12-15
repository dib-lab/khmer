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

// The index building algorithm is based on the samtools implementation which
// is under the MIT License:

/* The MIT License

   Copyright (c) 2008-2009 Genome Research Ltd.

   Permission is hereby granted, free of charge, to any person obtaining a copy
   of this software and associated documentation files (the "Software"), to deal
   in the Software without restriction, including without limitation the rights
   to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
   copies of the Software, and to permit persons to whom the Software is
   furnished to do so, subject to the following conditions:

   The above copyright notice and this permission notice shall be included in
   all copies or substantial portions of the Software.

   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
   IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
   FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
   AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
   LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
   OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
   THE SOFTWARE.
*/

#ifndef CORE_INCLUDE_SEQAN_BAM_IO_BAM_INDEX_BAI_H_
#define CORE_INCLUDE_SEQAN_BAM_IO_BAM_INDEX_BAI_H_

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// ----------------------------------------------------------------------------
// Tag Bai
// ----------------------------------------------------------------------------

struct Bai_;
typedef Tag<Bai_> Bai;

// ----------------------------------------------------------------------------
// Helper Class BaiBamIndexBinData_
// ----------------------------------------------------------------------------

// Store the information of a bin.

struct BaiBamIndexBinData_
{
    String<Pair<__uint64, __uint64> > chunkBegEnds;
};

// ----------------------------------------------------------------------------
// Spec BAI BamIndex
// ----------------------------------------------------------------------------

/**
.Spec.BAI BamIndex
..cat:BAM I/O
..general:Class.BamIndex
..summary:Access to BAI (samtools-style) Indices.
..signature:BamIndex<Bai>
..include:seqan/bam_io.h

.Memfunc.BAI BamIndex#BamIndex
..class:Spec.BAI BamIndex
..signature:BamIndex()
..summary:Constructor.
..remarks:Only the default constructor is provided.
*/

template <>
class BamIndex<Bai>
{
public:
    typedef std::map<__uint32, BaiBamIndexBinData_> TBinIndex_;
    typedef String<__uint64> TLinearIndex_;

    __uint64 _unalignedCount;

    // 1<<14 is the size of the minimum bin.
    static const __int32 BAM_LIDX_SHIFT = 14;
    
    String<TBinIndex_> _binIndices;
    String<TLinearIndex_> _linearIndices;
    
    BamIndex() : _unalignedCount(maxValue<__uint64>())
    {}
};

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function jumpToRegion()
// ----------------------------------------------------------------------------

/**
.Function.BamIndex#jumpToRegion
..class:Class.BamIndex
..cat:BAM I/O
..signature:jumpToRegion(bgzfStream, hasAlignments, bamIOContext, refId, pos, posEnd, bamIndex)
..summary:Seek in BAM BGZF stream using an index.
..remark:Note that because of the structure of BAI indices, you cannot simply jump to a position and you have to jump to region.
..param.bgzfStream:The BGZF Stream to seek in.
...type:Spec.BGZF Stream
..param.refId:Reference ID to seek to.
...type:nolink:$__int32$
..param.hasAlignments:Set to $true$ iff there are alignments at this position.
...type:nolink:$bool$
..param.bamIOContext:Context to use for loading alignments.
...type:Class.BamIOContext
..param.pos:Zero-based begin position in the reference.
...type:nolink:$__int32$
..param.pos:Zero-based (exclusive, C-style) end position in the reference.
...type:nolink:$__int32$
..param.bamIndex:The index to use.
...type:Class.BamIndex
..returns:$bool$ indicating success.
..remarks:This function may fail if the refId/pos is invalid.
..include:seqan/bam_io.h
*/

static inline void
_baiReg2bins(String<__uint16> & list, __uint32 beg, __uint32 end)
{
	unsigned k;
	if (beg >= end) return;
	if (end >= 1u<<29) end = 1u<<29;
	--end;
	appendValue(list, 0);
	for (k =    1 + (beg>>26); k <=    1 + (end>>26); ++k) appendValue(list, k);
	for (k =    9 + (beg>>23); k <=    9 + (end>>23); ++k) appendValue(list, k);
	for (k =   73 + (beg>>20); k <=   73 + (end>>20); ++k) appendValue(list, k);
	for (k =  585 + (beg>>17); k <=  585 + (end>>17); ++k) appendValue(list, k);
    for (k = 4681 + (beg>>14); k <= 4681 + (end>>14); ++k) appendValue(list, k);
}

template <typename TNameStore, typename TNameStoreCache>
inline bool
jumpToRegion(Stream<Bgzf> & stream,
             bool & hasAlignments,
             BamIOContext<TNameStore, TNameStoreCache> /*const*/ & bamIOContext,
             __int32 refId,
             __int32 pos,
             __int32 posEnd,
             BamIndex<Bai> const & index)
{
    hasAlignments = false;
    if (refId < 0)
        return false;  // Cannot seek to invalid reference.
    if (static_cast<unsigned>(refId) >= length(index._binIndices))
        return false;  // Cannot seek to invalid reference.

    // ------------------------------------------------------------------------
    // Compute offset in BGZF file.
    // ------------------------------------------------------------------------
    __uint64 offset = MaxValue<__uint64>::VALUE;

    // Retrieve the candidate bin identifiers for [pos, posEnd).
    String<__uint16> candidateBins;
    _baiReg2bins(candidateBins, pos, posEnd);

    // Retrieve the smallest required offset from the linear index.
    unsigned windowIdx = pos >> 14;  // Linear index consists of 16kb windows.
    __uint64 linearMinOffset = 0;
    if (windowIdx >= length(index._linearIndices[refId]))
    {
        // TODO(holtgrew): Can we simply always take case 1?

        // This is the case were we want to jump in a non-existing window.
        //
        // If there are no linear indices for this reference then we use the linear min offset of the next
        // reference that has an linear index.
        if (empty(index._linearIndices[refId]))
        {
            for (unsigned i = refId; i < length(index._linearIndices); ++i)
            {
                if (!empty(index._linearIndices[i]))
                {
                    linearMinOffset = front(index._linearIndices[i]);
                    if (linearMinOffset != 0u)
                        break;
                    for (unsigned j = 1; j < length(index._linearIndices[i]); ++j)
                    {
                        if (index._linearIndices[i][j] > linearMinOffset)
                        {
                            linearMinOffset = index._linearIndices[i][j];
                            break;
                        }
                    }
                    if (linearMinOffset != 0u)
                        break;
                }
            }
        }
        else
        {
            linearMinOffset = back(index._linearIndices[refId]);
        }
    }
    else
    {
        linearMinOffset = index._linearIndices[refId][windowIdx];
    }

    // Combine candidate bins and smallest required offset from linear index into candidate offset.
    typedef std::set<__uint64> TOffsetCandidates;
    TOffsetCandidates offsetCandidates;
    typedef typename Iterator<String<__uint16>, Rooted>::Type TCandidateIter;
    for (TCandidateIter it = begin(candidateBins, Rooted()); !atEnd(it); goNext(it))
    {
        typedef typename std::map<__uint32, BaiBamIndexBinData_>::const_iterator TMapIter;
        TMapIter mIt = index._binIndices[refId].find(*it);
        if (mIt == index._binIndices[refId].end())
            continue;  // Candidate is not in index!

        typedef typename Iterator<String<Pair<__uint64, __uint64> > const, Rooted>::Type TBegEndIter;
        for (TBegEndIter it2 = begin(mIt->second.chunkBegEnds, Rooted()); !atEnd(it2); goNext(it2))
            if (it2->i2 >= linearMinOffset)
                offsetCandidates.insert(it2->i1);
    }

    // Search through candidate offsets, find smallest with a fitting alignment.
    //
    // Note that it is not necessarily the first.
    //
    // TODO(holtgrew): Can this be optimized similar to how bamtools does it?
    typedef typename TOffsetCandidates::const_iterator TOffsetCandidateIter;
    BamAlignmentRecord record;
    for (TOffsetCandidateIter candIt = offsetCandidates.begin(); candIt != offsetCandidates.end(); ++candIt)
    {
        if (streamSeek(stream, *candIt, SEEK_SET) != 0)
            return false;  // Error while seeking.
        if (readRecord(record, bamIOContext, stream, Bam()) != 0)
            return false;  // Error while reading.

        // std::cerr << "record.beginPos == " << record.beginPos << "\n";
        // __int32 endPos = record.beginPos + getAlignmentLengthInRef(record);
        if (record.rID != refId)
            continue;  // Wrong contig.
        if (record.beginPos >= posEnd)
            continue;  // Cannot overlap with [pos, posEnd).

        // Found an alignment.
        hasAlignments = true;
        offset = *candIt;
        // std::cerr << "offset == " << offset << "\n";
        break;
    }

    if (offset != MaxValue<__uint64>::VALUE)
    {
        if (streamSeek(stream, offset, SEEK_SET) != 0)
            return false;  // Error while seeking.
    }
    // Finding no overlapping alignment is not an error, hasAlignments is false.
    return true;
}

// ----------------------------------------------------------------------------
// Function jumpToOrphans()
// ----------------------------------------------------------------------------

/**
.Function.BamIndex#jumpToOrphans
..class:Class.BamIndex
..cat:BAM I/O
..signature:jumpToOrphans(bgzfStream, bamIOContext, bamIndex)
..summary:Seek to orphans block in BAM BGZF stream using an index.
..param.bgzfStream:The BGZF Stream to seek in.
...type:Spec.BGZF Stream
..param.bamIOContext:Context to use for loading alignments.
...type:Class.BamIOContext
..param.bamIndex:The index to use.
...type:Class.BamIndex
..returns:$bool$ indicating success.
..include:seqan/bam_io.h
*/

template <typename TNameStore, typename TNameStoreCache>
bool jumpToOrphans(Stream<Bgzf> & stream,
                   bool & hasAlignments,
                   BamIOContext<TNameStore, TNameStoreCache> /*const*/ & bamIOContext,
                   BamIndex<Bai> const & index)
{
    hasAlignments = false;

    // Search linear indices for the largest entry of all references.
    __uint64 aliOffset = MaxValue<__uint64>::VALUE;
    for (int i = length(index._linearIndices) - 1; i >= 0; --i)
        if (!empty(index._linearIndices[i]))
        {
            aliOffset = back(index._linearIndices[i]);
            break;
        }
    if (aliOffset == MaxValue<__uint64>::VALUE)
        return false;  // No offset found.

    // Get index of the first orphan alignment by seeking from linear index bucket.
    BamAlignmentRecord record;
    __uint64 offset = MaxValue<__uint64>::VALUE;
    __uint64 result = 0;
    int res = streamSeek(stream, aliOffset, SEEK_SET);
    if (res != 0)
        return false;  // Error while seeking.
    while (!atEnd(stream))
    {
        result = streamTell(stream);
        res = readRecord(record, bamIOContext, stream, Bam());
        if (res != 0)
            return false;  // Error while reading.
        if (record.rID == -1)
        {
            // Found alignment.
            hasAlignments = true;
            offset = result;
            break;
        }
    }

    // Jump back to the first alignment.
    if (offset != MaxValue<__uint64>::VALUE)
    {
        int res = streamSeek(stream, offset, SEEK_SET);
        if (res != 0)
            return false;  // Error while seeking.
    }

    // Finding no orphan alignment is not an error, hasAilgnments is false then.
    return true;
}

// ----------------------------------------------------------------------------
// Function getUnalignedCount()
// ----------------------------------------------------------------------------

/**
.Function.BamIndex#getUnalignedCount
..class:Class.BamIndex
..cat:BAM I/O
..signature:getUnalignedCount(index)
..summary:Query index for number of unaligned reads.
..param.index:Index to query.
...type:Class.BamIndex
..returns:$__uint64$ with number of unaligned reads.
..include:seqan/bam_io.h
*/

inline __uint64
getUnalignedCount(BamIndex<Bai> const & index)
{
    return index._unalignedCount;
}

// ----------------------------------------------------------------------------
// Function read()
// ----------------------------------------------------------------------------

/**
.Function.BamIndex#read
..class:Class.BamIndex
..cat:BAM I/O
..signature:read(index, filename)
..summary:Load a BAM index from a given file name.
..param.index:Target data structure.
...type:Class.BamIndex
..param.filename:Path to file to load.
...type:nolink:$char const *$
..returns:$int$ status code, $0$ indicating success.
..include:seqan/bam_io.h
 */

inline int
read(BamIndex<Bai> & index, char const * filename)
{
    std::fstream fin(filename, std::ios::binary | std::ios::in);
    if (!fin.good())
        return 1;  // Could not open file.

    // Read magic number.
    CharString buffer;
    resize(buffer, 4);
    fin.read(&buffer[0], 4);
    if (!fin.good())
        return 1;
    if (buffer != "BAI\1")
        return 1;  // Magic number is wrong.

    __int32 nRef = 0;
    fin.read(reinterpret_cast<char *>(&nRef), 4);
    if (!fin.good())
        return 1;

    resize(index._linearIndices, nRef);
    resize(index._binIndices, nRef);
    
    for (int i = 0; i < nRef; ++i)  // For each reference.
    {
        // Read bin index.
        __int32 nBin = 0;
        fin.read(reinterpret_cast<char *>(&nBin), 4);
        if (!fin.good())
            return 1;
        index._binIndices[i].clear();
        BaiBamIndexBinData_ data;
        for (int j = 0; j < nBin; ++j)  // For each bin.
        {
            clear(data.chunkBegEnds);

            __uint32 bin = 0;
            fin.read(reinterpret_cast<char *>(&bin), 4);
            if (!fin.good())
                return 1;

            __int32 nChunk = 0;
            fin.read(reinterpret_cast<char *>(&nChunk), 4);
            if (!fin.good())
                return 1;
            reserve(data.chunkBegEnds, nChunk);
            for (int k = 0; k < nChunk; ++k)  // For each chunk;
            {
                __uint64 chunkBeg = 0;
                __uint64 chunkEnd = 0;
                fin.read(reinterpret_cast<char *>(&chunkBeg), 8);
                fin.read(reinterpret_cast<char *>(&chunkEnd), 8);
                if (!fin.good())
                    return 1;
                appendValue(data.chunkBegEnds, Pair<__uint64>(chunkBeg, chunkEnd));
            }

            // Copy bin data into index.
            index._binIndices[i][bin] = data;
        }

        // Read linear index.
        __int32 nIntv = 0;
        fin.read(reinterpret_cast<char *>(&nIntv), 4);
        if (!fin.good())
            return 1;
        clear(index._linearIndices[i]);
        reserve(index._linearIndices[i], nIntv);
        for (int j = 0; j < nIntv; ++j)
        {
            __uint64 ioffset = 0;
            fin.read(reinterpret_cast<char *>(&ioffset), 8);
            if (!fin.good())
                return 1;
            appendValue(index._linearIndices[i], ioffset);
        }
    }

    if (!fin.good())
        return 1;

    // Read (optional) number of alignments without coordinate.
    __uint64 nNoCoord = 0;
    fin.read(reinterpret_cast<char *>(&nNoCoord), 8);
    if (!fin.good())
    {
        fin.clear();
        nNoCoord = 0;
    }
    index._unalignedCount = nNoCoord;

    return 0;
}

// TODO(holtgrew): This is only here because of the read() function with TSequence in old file.h.

inline int
read(BamIndex<Bai> & index, char * filename)
{
    return read(index, static_cast<char const *>(filename));
}

// ----------------------------------------------------------------------------
// Function buildIndex()
// ----------------------------------------------------------------------------

/*DISABLED
.Function.BamIndex#buildIndex
..class:Class.BamIndex
..cat:BAM I/O
..signature:buildIndex(index, filename)
..summary:Build index for BAM file with given filename.
..remarks:This will create an index file named $filename + ".bai"$.
..param.index:Target data structure.
...type:Class.BamIndex
..param.filename:Path to BAM file to load.
...type:nolink:$char const *$
..returns:$bool$ indicating success.
..include:seqan/bam_io.h
 */

inline int _writeIndex(BamIndex<Bai> const & index, char const * filename)
{
    std::cerr << "WRITE INDEX TO " << filename << std::endl;
    // Open output stream.
    std::ofstream out(filename, std::ios::binary | std::ios::out);

    SEQAN_ASSERT_EQ(length(index._binIndices), length(index._linearIndices));
    
    // Write header.
    out.write("BAI\1", 4);
    __int32 numRefSeqs = length(index._binIndices);
    out.write(reinterpret_cast<char *>(&numRefSeqs), 4);

    // Write out indices.
    typedef BamIndex<Bai> const                TBamIndex;
    typedef TBamIndex::TBinIndex_ const        TBinIndex;
    typedef TBinIndex::const_iterator          TBinIndexIter;
    typedef TBamIndex::TLinearIndex_           TLinearIndex;
    for (int i = 0; i < numRefSeqs; ++i)
    {
        TBinIndex const & binIndex = index._binIndices[i];
        TLinearIndex const & linearIndex = index._linearIndices[i];

        // Write out binning index.
        __int32 numBins = binIndex.size();
        out.write(reinterpret_cast<char *>(&numBins), 4);
        for (TBinIndexIter itB = binIndex.begin(), itBEnd = binIndex.end(); itB != itBEnd; ++itB)
        {
            // Write out bin id.
            out.write(reinterpret_cast<char const *>(&itB->first), 4);
            // Write out number of chunks.
            __uint32 numChunks = length(itB->second.chunkBegEnds);
            out.write(reinterpret_cast<char *>(&numChunks), 4);
            // Write out all chunks.
            typedef Iterator<String<Pair<__uint64> > const, Rooted>::Type TChunkIter;
            for (TChunkIter itC = begin(itB->second.chunkBegEnds); !atEnd(itC); goNext(itC))
            {
                out.write(reinterpret_cast<char const *>(&itC->i1), 8);
                out.write(reinterpret_cast<char const *>(&itC->i2), 8);
            }
        }

        // Write out linear index.
        __int32 numIntervals = length(index._linearIndices);
        out.write(reinterpret_cast<char *>(&numIntervals), 4);
        typedef Iterator<String<__uint64> const, Rooted>::Type TLinearIndexIter;
        for (TLinearIndexIter it = begin(linearIndex, Rooted()); !atEnd(it); goNext(it))
            out.write(reinterpret_cast<char const *>(&*it), 8);
    }

    // Write the number of unaligned reads if set.
    std::cerr << "UNALIGNED\t" << index._unalignedCount << std::endl;
    if (index._unalignedCount != maxValue<__uint64>())
        out.write(reinterpret_cast<char const *>(&index._unalignedCount), 8);

    return !out.good();  // 1 on error, 0 on success.
}

inline void _baiAddAlignmentChunkToBin(BamIndex<Bai> & index,
                                       __uint32 currBin,
                                       __uint32 currOffset,
                                       __uint64 prevOffset)
{
    // If this is not the first reference sequence then add previous reference data.
    Pair<__uint64> newChunk(currOffset, prevOffset);

    // If no interest exists yet for this bin, create one and store alignment chunk.
    BamIndex<Bai>::TBinIndex_::iterator binIter = back(index._binIndices).find(currBin);
    if (binIter == back(index._binIndices).end())
    {
        BaiBamIndexBinData_ binData;
        appendValue(binData.chunkBegEnds, newChunk);
        back(index._binIndices).insert(std::make_pair(currBin, binData));
    }
    else
    {
        // Otherwise, just append alignment chunk.
        appendValue(binIter->second.chunkBegEnds, newChunk);
    }
}

inline bool
buildIndex(BamIndex<Bai> & index, char const * filename)
{
    SEQAN_FAIL("This does not work ye!");

    index._unalignedCount = 0;
    clear(index._binIndices);
    clear(index._linearIndices);
    
    // Open BAM file for reading.
    Stream<Bgzf> bamStream;
    if (!open(bamStream, filename, "r"))
        return false;  // Could not open BAM file.

    // Initialize BamIOContext.
    typedef StringSet<CharString>      TNameStore;
    typedef NameStoreCache<TNameStore> TNameStoreCache;
    
    TNameStore refNameStore;
    TNameStoreCache refNameStoreCache(refNameStore);
    BamIOContext<TNameStore> bamIOContext(refNameStore, refNameStoreCache);

    // Read BAM header.
    BamHeader header;
    int res = readRecord(header, bamIOContext, bamStream, Bam());
    if (res != 0)
        return false;  // Could not read BAM header.
    __uint32 numRefSeqs = length(header.sequenceInfos);

    // Scan over BAM file and create index.
    BamAlignmentRecord record;
    __uint32 currBin    = maxValue<__uint32>();
    __uint32 prevBin    = maxValue<__uint32>();
    __int32 currRefId   = BamAlignmentRecord::INVALID_REFID;
    __int32 prevRefId   = BamAlignmentRecord::INVALID_REFID;
    __uint64 currOffset = streamTell(bamStream);
    __uint64 prevOffset = currOffset;
    __int32 prevPos     = minValue<__int32>();

    while (!atEnd(bamStream))
    {
        // Load next record.
        res = readRecord(record, bamIOContext, bamStream, Bam());
        if (res != 0)
            return false;
        
        // Check ordering.
        if (prevRefId == record.rID && prevPos > record.beginPos)
            return false;

        // The reference sequence changed, close bins for previous reference.
        if (prevRefId != record.rID)
        {
            if (prevRefId != BamAlignmentRecord::INVALID_REFID)
            {
                _baiAddAlignmentChunkToBin(index, currBin, currOffset, prevOffset);

                // Add an index for all empty references between prevRefId (excluded) and record.rID (included).
                for (int i = prevRefId + 1; i < record.rID; ++i)
                {
                    BamIndex<Bai>::TBinIndex_ binIndex;
                    appendValue(index._binIndices, binIndex);
                    BamIndex<Bai>::TLinearIndex_ linearIndex;
                    appendValue(index._linearIndices, linearIndex);
                }

                // Update bin book keeping.
                currOffset = prevOffset;
                currBin    = record.bin;
                prevBin    = record.bin;
                currRefId  = record.rID;
            }
            else
            {
                // Otherwise, this is the first pass.  Create an index for all empty references up to and including
                // current refId.
                for (int i = 0; i < record.rID; ++i)
                {
                    BamIndex<Bai>::TBinIndex_ binIndex;
                    appendValue(index._binIndices, binIndex);
                    BamIndex<Bai>::TLinearIndex_ linearIndex;
                    appendValue(index._linearIndices, linearIndex);
                }
            }

            // Update reference book keeping.
            prevRefId = record.rID;
            prevBin = minValue<__int32>();
        }

        // If the alignment's reference id is valid and its bin is not a leaf.
        if (record.rID >= 0 && record.bin < 4681)
        {
            __int32 beginOffset = record.beginPos >> BamIndex<Bai>::BAM_LIDX_SHIFT;
            __int32 endPos      = getAlignmentLengthInRef(record);
            __int32 endOffset   = (endPos - 1) >> BamIndex<Bai>::BAM_LIDX_SHIFT;

            // Resize linear index if necessary.
            unsigned oldSize = length(index._linearIndices);
            unsigned newSize = endOffset + 1;
            if (oldSize < newSize)
                resize(back(index._linearIndices), newSize, 0);

            // Store offset.
            for (int i = beginOffset + 1; i <= endOffset; ++i)
                if (back(index._linearIndices)[i] == 0u)
                    back(index._linearIndices)[i] = prevOffset;
        }

        // Handle the case if we changed to a new BAI bin.
        if (record.bin != prevBin)
        {
            // If not first bin of reference, save previous bin data.
            if (currBin != maxValue<__uint32>())
                _baiAddAlignmentChunkToBin(index, currBin, currOffset, prevOffset);

            // Update markers.
            currOffset = prevOffset;
            currBin    = record.bin;
            prevBin    = record.bin;
            currRefId  = record.rID;

            // If the reference id is invalid then break out.
            if (currRefId < 0)
                break;
        }

        // Make sure that the current file pointer is beyond prevOffset.
        if (streamTell(bamStream) <= static_cast<__int64>(prevOffset))
            return false;  // Calculating offsets failed.

        // Update prevOffset and prevPos.
        prevOffset = streamTell(bamStream);
        prevPos    = record.beginPos;
    }

    // Count remaining unaligned records.
    while (!streamEof(bamStream))
    {
        SEQAN_ASSERT_GT(index._unalignedCount, 0u);

        res = readRecord(record, bamIOContext, bamStream, Bam());
        if (res != 0 || record.rID >= 0)
            return false;  // Could not read record.

        index._unalignedCount += 1;
    }

    // After loading all alignments, if any data was read, perform checks.
    if (currRefId >= 0)
    {
        // Store last alignment chunk to its bin and then write last reference entry with data.
        _baiAddAlignmentChunkToBin(index, currBin, currOffset, prevOffset);

        // Finally, write any empty references remaining at end of file.
        SEQAN_ASSERT_GEQ(numRefSeqs, length(index._binIndices));
        BamIndex<Bai>::TBinIndex_ binIndex;
        resize(index._binIndices, numRefSeqs, binIndex);  // TODO(holtgrew): binIndex is unnecessary if resize used T() as default value as vector.resize() does.
    }

    // Merge small bins if possible.
    SEQAN_FAIL("TODO: Merge bins!");

    // Write out index.
    CharString baiFilename(filename);
    append(baiFilename, ".bai");
    res = _writeIndex(index, toCString(baiFilename));

    return (res == 0);
}

}  // namespace seqan

#endif  // #ifndef CORE_INCLUDE_SEQAN_BAM_IO_BAM_INDEX_BAI_H_
