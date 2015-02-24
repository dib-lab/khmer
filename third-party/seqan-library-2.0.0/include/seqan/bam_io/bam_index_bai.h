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

#ifndef INCLUDE_SEQAN_BAM_IO_BAM_INDEX_BAI_H_
#define INCLUDE_SEQAN_BAM_IO_BAM_INDEX_BAI_H_

namespace seqan {

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// ----------------------------------------------------------------------------
// Class BamIndex
// ----------------------------------------------------------------------------

/*!
 * @class BamIndex
 * @headerfile <seqan/bam_io.h>
 *
 * @brief Access to BAM indices.
 *
 * @signature template <typename TSpec>
 *            class BamIndex;
 *
 * This is an abstract class; don't use it itself but its specializations.
 *
 * @see BamFileIn
 */

template <typename TSpec>
class BamIndex;

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

/*!
 * @class BaiBamIndex
 * @headerfile <seqan/bam_io.h>
 * @extends BamIndex
 * @brief Access to BAI (samtools-style).
 *
 * @signature template <>
 *            class BamIndex<Bai>;
 */

/*!
 * @fn BaiBamIndex::BamIndex
 * @brief Constructor.
 *
 * @signature BamIndex::BamIndex();
 *
 * @section Remarks
 *
 * Only the default constructor is provided.
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
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function jumpToRegion()
// ----------------------------------------------------------------------------

/*!
 * @fn BamFileIn#jumpToRegion
 * @brief Seek in BamFileIn using an index.
 *
 * You provide a region <tt>[pos, posEnd)</tt> on the reference <tt>refID</tt> that you want to jump to and the function
 * jumps to the first alignment in this region, if any.
 *
 * @signature bool jumpToRegion(bamFileIn, hasAlignments, refID, pos, posEnd, index);
 *
 * @param[in,out] bamFileIn     The @link BamFileIn @endlink to jump with.
 * @param[out]    hasAlignments A <tt>bool</tt> that is set true if the region <tt>[pos, posEnd)</tt> has any
 *                              alignments.
 * @param[in]     refID         The reference id to jump to (<tt>__int32</tt>).
 * @param[in]     pos           The begin of the region to jump to (<tt>__int32</tt>).
 * @param[in]     posEnd        The end of the region to jump to (<tt>__int32</tt>).
 * @param[in]     index         The @link BamIndex @endlink to use for the jumping.
 *
 * @return bool true if seeking was successful, false if not.
 *
 * @section Remarks
 *
 * This function fails if <tt>refID</tt>/<tt>pos</tt> are invalid.
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

template <typename TSpec>
inline bool
jumpToRegion(FormattedFile<Bam, Input, TSpec> & bamFile,
             bool & hasAlignments,
             __int32 refId,
             __int32 pos,
             __int32 posEnd,
             BamIndex<Bai> const & index)
{
    if (!isEqual(format(bamFile), Bam()))
        return false;

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

    // Search through candidate offsets, find rightmost possible.
    //
    // Note that it is not necessarily the first.
    //
    // TODO(holtgrew): Can this be optimized similar to how bamtools does it?
    typedef typename TOffsetCandidates::const_iterator TOffsetCandidateIter;
    BamAlignmentRecord record;
    for (TOffsetCandidateIter candIt = offsetCandidates.begin(); candIt != offsetCandidates.end(); ++candIt)
    {
        setPosition(bamFile, *candIt);

        readRecord(record, bamFile);

        // std::cerr << "record.beginPos == " << record.beginPos << "\n";
        // __int32 endPos = record.beginPos + getAlignmentLengthInRef(record);
        if (record.rID != refId)
            continue;  // Wrong contig.
        if (!hasAlignments || record.beginPos <= pos)
        {
            // Found a valid alignment.
            hasAlignments = true;
            offset = *candIt;
        }

        if (record.beginPos >= posEnd)
            break;  // Cannot find overlapping any more.
    }

    if (offset != MaxValue<__uint64>::VALUE)
        setPosition(bamFile, offset);

    // Finding no overlapping alignment is not an error, hasAlignments is false.
    return true;
}

// ----------------------------------------------------------------------------
// Function jumpToOrphans()
// ----------------------------------------------------------------------------

/*!
 * @fn BamFileIn#jumpToOrphans
 * @brief Seek to orphans block in BamFileIn using an index.
 *
 * @signature bool jumpToOrphans(bamFileIn, hasAlignments, index);
 *
 * @param[in,out] bamFileIn      The @link BamFileIn @endlink object to jump with.
 * @param[out]    hasAlignments  A <tt>bool</tt> that is set to true if there are any orphans.
 * @param[in]     index          The @link BamIndex @endlink to use for jumping.
 */

template <typename TSpec, typename TNameStore, typename TNameStoreCache>
bool jumpToOrphans(FormattedFile<Bam, Input, TSpec> & bamFile,
                   bool & hasAlignments,
                   BamIndex<Bai> const & index)
{
    if (!isEqual(format(bamFile), Bam()))
        return false;

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
    if (!setPosition(bamFile, aliOffset))
        return false;  // Error while seeking.
    while (!atEnd(bamFile))
    {
        result = position(bamFile);
        readRecord(record, bamFile);
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
        if (!setPosition(bamFile, offset, SEEK_SET))
            return false;  // Error while seeking.
    }

    // Finding no orphan alignment is not an error, hasAilgnments is false then.
    return true;
}

// ----------------------------------------------------------------------------
// Function getUnalignedCount()
// ----------------------------------------------------------------------------

/*!
 * @fn BamIndex#getUnalignedCount
 * @brief Query index for number of unaligned reads.
 *
 * @signature __uint64 getUnalignedCount(index);
 *
 * @param[in] index     Index to query.
 * @return    __uint64  The number of unaligned reads.
 */

inline __uint64
getUnalignedCount(BamIndex<Bai> const & index)
{
    return index._unalignedCount;
}

// ----------------------------------------------------------------------------
// Function open()
// ----------------------------------------------------------------------------

/*!
 * @fn BamIndex#open
 * @brief Load a BAM index from a given file name.
 * @signature bool open(index, filename);

 * @param[in,out] index    Target data structure.
 * @param[in]     filename Path to file to load. Types: char const *
 *
 * @return        bool     Returns <tt>true</tt> on success, false otherwise.
 */

inline bool
open(BamIndex<Bai> & index, char const * filename)
{
    std::fstream fin(filename, std::ios::binary | std::ios::in);
    if (!fin.good())
        return false;  // Could not open file.

    // Read magic number.
    CharString buffer;
    resize(buffer, 4);
    fin.read(&buffer[0], 4);
    if (!fin.good())
        return false;
    if (buffer != "BAI\1")
        return false;  // Magic number is wrong.

    __int32 nRef = 0;
    fin.read(reinterpret_cast<char *>(&nRef), 4);
    if (!fin.good())
        return false;

    resize(index._linearIndices, nRef);
    resize(index._binIndices, nRef);

    for (int i = 0; i < nRef; ++i)  // For each reference.
    {
        // Read bin index.
        __int32 nBin = 0;
        fin.read(reinterpret_cast<char *>(&nBin), 4);
        if (!fin.good())
            return false;
        index._binIndices[i].clear();
        BaiBamIndexBinData_ data;
        for (int j = 0; j < nBin; ++j)  // For each bin.
        {
            clear(data.chunkBegEnds);

            __uint32 bin = 0;
            fin.read(reinterpret_cast<char *>(&bin), 4);
            if (!fin.good())
                return false;

            __int32 nChunk = 0;
            fin.read(reinterpret_cast<char *>(&nChunk), 4);
            if (!fin.good())
                return false;
            reserve(data.chunkBegEnds, nChunk);
            for (int k = 0; k < nChunk; ++k)  // For each chunk;
            {
                __uint64 chunkBeg = 0;
                __uint64 chunkEnd = 0;
                fin.read(reinterpret_cast<char *>(&chunkBeg), 8);
                fin.read(reinterpret_cast<char *>(&chunkEnd), 8);
                if (!fin.good())
                    return false;
                appendValue(data.chunkBegEnds, Pair<__uint64>(chunkBeg, chunkEnd));
            }

            // Copy bin data into index.
            index._binIndices[i][bin] = data;
        }

        // Read linear index.
        __int32 nIntv = 0;
        fin.read(reinterpret_cast<char *>(&nIntv), 4);
        if (!fin.good())
            return false;
        clear(index._linearIndices[i]);
        reserve(index._linearIndices[i], nIntv);
        for (int j = 0; j < nIntv; ++j)
        {
            __uint64 ioffset = 0;
            fin.read(reinterpret_cast<char *>(&ioffset), 8);
            if (!fin.good())
                return false;
            appendValue(index._linearIndices[i], ioffset);
        }
    }

    if (!fin.good())
        return false;

    // Read (optional) number of alignments without coordinate.
    __uint64 nNoCoord = 0;
    fin.read(reinterpret_cast<char *>(&nNoCoord), 8);
    if (!fin.good())
    {
        fin.clear();
        nNoCoord = 0;
    }
    index._unalignedCount = nNoCoord;

    return true;
}

// TODO(holtgrew): This is only here because of the read() function with TSequence in old file.h.

inline bool
open(BamIndex<Bai> & index, char * filename)
{
    return open(index, static_cast<char const *>(filename));
}

// ----------------------------------------------------------------------------
// Function buildIndex()
// ----------------------------------------------------------------------------

inline bool _saveIndex(BamIndex<Bai> const & index, char const * filename)
{
    std::cerr << "WRITE INDEX TO " << filename << std::endl;
    // Open output file.
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

    return out.good();  // false on error, true on success.
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
    SEQAN_FAIL("This does not work yet!");

    index._unalignedCount = 0;
    clear(index._binIndices);
    clear(index._linearIndices);

    // Open BAM file for reading.
    BamFileIn bamFile;
    if (!open(bamFile, filename))
        return false;  // Could not open BAM file.

    // Read BAM header.
    BamHeader header;
    readHeader(header, bamFile);

    __uint32 numRefSeqs = length(contigNames(context(bamFile)));

    // Scan over BAM file and create index.
    BamAlignmentRecord record;
    __uint32 currBin    = maxValue<__uint32>();
    __uint32 prevBin    = maxValue<__uint32>();
    __int32 currRefId   = BamAlignmentRecord::INVALID_REFID;
    __int32 prevRefId   = BamAlignmentRecord::INVALID_REFID;
    __uint64 currOffset = position(bamFile);
    __uint64 prevOffset = currOffset;
    __int32 prevPos     = minValue<__int32>();

    while (!atEnd(bamFile))
    {
        // Load next record.
        readRecord(record, bamFile);

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
        if (position(bamFile) <= static_cast<__int64>(prevOffset))
            return false;  // Calculating offsets failed.

        // Update prevOffset and prevPos.
        prevOffset = position(bamFile);
        prevPos    = record.beginPos;
    }

    // Count remaining unaligned records.
    while (!atEnd(bamFile))
    {
        SEQAN_ASSERT_GT(index._unalignedCount, 0u);

        readRecord(record, bamFile);
        if (record.rID >= 0)
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
    return _saveIndex(index, toCString(baiFilename));
}

}  // namespace seqan

#endif  // #ifndef INCLUDE_SEQAN_BAM_IO_BAM_INDEX_BAI_H_
