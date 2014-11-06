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
// Code for loading and writing FAI FASTA index files.
// ==========================================================================

#include <iostream>

#include <seqan/sequence.h>
#include <seqan/store.h>
#include <seqan/file.h>
#include <seqan/stream.h>

#ifndef SEQAN_CORE_APPS_RABEMA_FAI_INDEX_H_
#define SEQAN_CORE_APPS_RABEMA_FAI_INDEX_H_

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// ----------------------------------------------------------------------------
// Class FaiIndexEntry_
// ----------------------------------------------------------------------------

// Stores one entry in the FAI file.

class FaiIndexEntry_
{
public:
    // Name of reference sequence.
    CharString name;
    // Number of nucleotides in sequence.
    __uint64 sequenceLength;
    // Offset in the file.
    __uint64 offset;
    // Number of sequence characters per line.
    unsigned lineLength;
    // Number of overall characters per line, including newline character(s).
    unsigned overallLineLength;

    FaiIndexEntry_() :
        sequenceLength(0), offset(0), lineLength(0), overallLineLength(0)
    {}
};

// ----------------------------------------------------------------------------
// Class FaiIndex
// ----------------------------------------------------------------------------

/**
.Class.FaiIndex
..cat:Input/Output
..signature:FaiIndex
..summary:Data type for storing FAI indices.
..wiki:Tutorial/IndexedFastaIO|Tutorial: Indexed FASTA I/O
..example.text:The following example demonstrate the usage of the FAIIndex class.
..example.file:demos/seq_io/fai_index_example.cpp
..include:seqan/seq_io.h

.Memvar.FaiIndex#FaiIndex
..class:Class.FaiIndex
..signature:FaiIndex()
..summary:The @Class.FaiIndex@ class only provides the default constructor.
*/

class FaiIndex
{
public:
    // TODO(holtgrew): The members are not publically documented for now. Keep it this way?

    // The name of the FASTA file.
    CharString fastaFilename;
    // The name of the FAI file.
    CharString faiFilename;

    // The index entries.
    String<FaiIndexEntry_> indexEntryStore;
    // A store for the reference names.
    StringSet<CharString> refNameStore;
    // A cache for fast access to the reference name store.
    NameStoreCache<StringSet<CharString> > refNameStoreCache;

    // We use this memory mapped string (opened read-only) to read from the file.
    String<char, MMap<> > mmapString;
    bool mmapStringOpen;

    FaiIndex() :
        refNameStoreCache(refNameStore), mmapStringOpen(false)
    {}
};

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function clear()
// ----------------------------------------------------------------------------

/**
.Function.FaiIndex#clear
..cat:Input/Output
..class:Class.FaiIndex
..signature:clear(faiIndex)
..param.faiIndex:The @Class.FaiIndex@ to reset.
...type:Class.FaiIndex
..summary:Reset a @Class.FaiIndex@ object to the state after default construction.
..include:seqan/seq_io.h
*/

inline void clear(FaiIndex & index)
{
    clear(index.fastaFilename);
    clear(index.faiFilename);
    clear(index.indexEntryStore);
    clear(index.refNameStore);
    refresh(index.refNameStoreCache);
}

// ----------------------------------------------------------------------------
// Function getIdByName()
// ----------------------------------------------------------------------------

/**
.Function.FaiIndex#getIdByName
..cat:Input/Output
..class:Class.FaiIndex
..signature:getIdByName(faiIndex, name, refId)
..summary:Return id (numeric index in the file) of a sequence in a FAI file.
..param.faiIndex:The @Class.FaiIndex@ to query.
...type:Class.FaiIndex
..param.name:The name of the sequence to get the id for.
...type:Shortcut.CharString
..param.refId:The id of the sequence is written here.
...type:nolink:$unsigned$
..return:$bool$ indicating whether a reference with this name is known in the index.
..include:seqan/seq_io.h
*/

// TODO(holtgrew): Fix parameter order when getIdByName() has good parameter order.

// TODO(holtgrew): Using templates here because of ambiguities hit otherwise.

template <typename TName, typename TId>
inline bool getIdByName(FaiIndex & index, TName const & name, TId & id)
{
    return getIdByName(index.refNameStore, name, id, index.refNameStoreCache);
}

template <typename TName, typename TId>
inline bool getIdByName(FaiIndex const & index, TName const & name, TId & id)
{
    return getIdByName(index.refNameStore, name, id, index.refNameStoreCache);
}

// ----------------------------------------------------------------------------
// Function sequenceLength()
// ----------------------------------------------------------------------------

/**
.Function.FaiIndex#sequenceLength
..cat:Input/Output
..class:Class.FaiIndex
..signature:__uint64 sequenceLength(faiIndex, refId)
..summary:Return length of the sequence with the given id in the @Class.FaiIndex@.
..param.faiIndex:The @Class.FaiIndex@ to query.
...type:Class.FaiIndex
..param.refId:The id of the sequence to get the length of.
...type:nolink:$unsigned$
..return:$__uint64$ with the length of the sequence.
..include:seqan/seq_io.h
*/

template <typename TRefId>
inline __uint64 sequenceLength(FaiIndex const & index, TRefId refId)
{
    return index.indexEntryStore[refId].sequenceLength;
}

// TODO(holtgrew): Wrapper and template only here because sequenceLength in string_set_base.h is weird.

template <typename TRefId>
inline __uint64 sequenceLength(FaiIndex & index, TRefId refId)
{
    return index.indexEntryStore[refId].sequenceLength;
}

// ----------------------------------------------------------------------------
// Function sequenceName()
// ----------------------------------------------------------------------------

/**
.Function.FaiIndex#sequenceName
..cat:Input/Output
..class:Class.FaiIndex
..signature:CharString sequenceName(faiIndex, refId)
..summary:Return the name of the sequence with the given id in the @Class.FaiIndex@.
..param.faiIndex:The @Class.FaiIndex@ to query.
...type:Class.FaiIndex
..param.refId:The id of the sequence to get the name of.
...type:Shortcut.CharString
..return:@Shortcut.CharString@ with the name of the sequence.
..include:seqan/seq_io.h
*/

inline CharString sequenceName(FaiIndex const & index, unsigned refId)
{
    return index.indexEntryStore[refId].name;
}

// ----------------------------------------------------------------------------
// Function numSeqs()
// ----------------------------------------------------------------------------

/**
.Function.FaiIndex#numSeqs
..cat:Input/Output
..class:Class.FaiIndex
..signature:numSeqs(faiIndex)
..summary:Return number of sequences known to an @Class.FaiIndex@.
..param.faiIndex:The @Class.FaiIndex@ to query.
...type:Class.FaiIndex
..return:$__uint64$ with the number of sequences.
..include:seqan/seq_io.h
*/

inline __uint64 numSeqs(FaiIndex const & index)
{
    return length(index.indexEntryStore);
}

// ----------------------------------------------------------------------------
// Function readRegion()
// ----------------------------------------------------------------------------

/**
.Function.FaiIndex#readRegion
..cat:Input/Output
..class:Class.FaiIndex
..signature:readRegion(str, faiIndex, refId, beginPos, endPos)
..signature:readRegion(str, faiIndex, region);
..summary:Load the infix of a sequence from a @Class.FaiIndex@.
..description:The given region is loaded from the indexed FASTA file.
..param.str:The sequence infix is written into this string.
...type:Class.String
..param.faiIndex:The @Class.FaiIndex@ to query.
...type:Class.FaiIndex
..param.refId:The index of the reference in the file.
...type:nolink:$unsigned$
..param.beginPos:The begin position of the infix to write to $str$.
...type:nolink:$unsigned$
..param.endPos:The end position of the infix to write to $str$.
...type:nolink:$unsigned$
..param.region:The @Class.GenomicRegion@ to read.
...type:Class.GenomicRegion
..return:Status code $int$, $0$ indicating success and $1$ an error.
..include:seqan/seq_io.h
*/

template <typename TValue, typename TSpec>
inline int readRegion(String<TValue, TSpec> & str,
                      FaiIndex const & index,
                      unsigned refId,
                      unsigned beginPos,
                      unsigned endPos)
{
    // Limit region to the infix, make sure that beginPos < endPos, compute character to read.
    unsigned seqLen = index.indexEntryStore[refId].sequenceLength;;
    beginPos = std::min(beginPos, seqLen);
    endPos = std::min(std::max(beginPos, endPos), seqLen);
    unsigned toRead = endPos - beginPos;

    typedef typename Iterator<String<char, MMap<> > const, Standard>::Type TSourceIter;
    typedef typename Iterator<String<TValue, TSpec>, Standard>::Type TTargetIter;
    TSourceIter itSource = begin(index.mmapString, Standard());
    __uint64 offset = index.indexEntryStore[refId].offset;
    // First, compute offset of the completely filled lines.
    unsigned numLines = beginPos / index.indexEntryStore[refId].lineLength;
    unsigned numBytes = numLines * index.indexEntryStore[refId].overallLineLength;
    // Then, compute overall offset by adding remaining bytes, too.
    numBytes += beginPos % index.indexEntryStore[refId].lineLength;
    offset += numBytes;
    // Advance iterator in MMap file.
    itSource += offset;

    // Copy out the characters from FASTA file and convert via iterator assignment to target string's type.
    resize(str, toRead, TValue());
    TTargetIter itTarget = begin(str, Standard());
    for (unsigned i = 0; i < toRead; )
    {
        if (isspace(*itSource))
        {
            ++itSource;
            continue;  // Skip spaces.
        }
        *itTarget = *itSource;
        ++itTarget;
        ++itSource;
        ++i;
    }

    return 0;
}

template <typename TValue, typename TSpec>
inline int readRegion(String<TValue, TSpec> & str,
                      FaiIndex const & index,
                      GenomicRegion const & region)
{
    int seqId = region.seqId;
    if (seqId == -1)
    {
        unsigned x = 0;
        if (!getIdByName(index, region.seqName, x))
            return 1;  // Sequence with this name could not be found.
        seqId = x;
    }
    int beginPos = region.beginPos;
    if (beginPos == -1)
        beginPos = 0;
    int endPos = region.endPos;
    if (endPos == -1)
        endPos = sequenceLength(index, seqId);
    return readRegion(str, index, seqId, beginPos, endPos);
}

// ----------------------------------------------------------------------------
// Function readSequence()
// ----------------------------------------------------------------------------

/**
.Function.FaiIndex#readSequence
..cat:Input/Output
..class:Class.FaiIndex
..signature:readSequence(str, faiIndex, refId)
..summary:Load a whole sequence from an @Class.FaiIndex@.
..param.str:The sequence is written into this string.
...type:Class.String
..param.faiIndex:The @Class.FaiIndex@ to use.
...type:Class.FaiIndex
..param.refId:The index of the reference in the file.
...type:nolink:$unsigned$
..return:Status code $int$, $0$ indicating success and $1$ an error.
..include:seqan/seq_io.h
*/

template <typename TValue, typename TSpec>
inline int readSequence(String<TValue, TSpec> & str, FaiIndex const & index, unsigned refId)
{
    if (refId > numSeqs(index))
        return 1;  // Out of bounds.

    return readRegion(str, index, refId, 0, sequenceLength(index, refId));
}

// ----------------------------------------------------------------------------
// Function read()
// ----------------------------------------------------------------------------

/**
.Function.FaiIndex#read
..cat:Input/Output
..class:Class.FaiIndex
..signature:read(faiIndex, fastaFileName[, faiFileName])
..summary:Read a FAI index.
..param.faiIndex:The @Class.FaiIndex@ object to read the file into.
...type:Class.FaiIndex
..param.fastaFileName:The name of the FASTA file to read.
...type:nolink:$char const *$
..param.faiFileName:The name of the FAI file to read.
...default:$fastaFileName + ".fai"$.
...type:nolink:$char const *$
..return:Status code $int$, $0$ indicating success and $1$ an error.
..include:seqan/seq_io.h
*/

inline int read(FaiIndex & index, char const * fastaFilename, char const * faiFilename)
{
    clear(index);  // Also clears filename, thus backup above and restore below.
    index.fastaFilename = fastaFilename;
    index.faiFilename = faiFilename;

    if (index.mmapStringOpen)
        close(index.mmapString);
    if (!open(index.mmapString, toCString(fastaFilename), OPEN_RDONLY))
        return 1;  // Could not open file.
    index.mmapStringOpen = true;

    // Open file.
    std::ifstream faiStream(toCString(index.faiFilename), std::ios::binary | std::ios::in);
    if (!faiStream.good())
        return 1;

    // Read FAI file.
    RecordReader<std::ifstream, SinglePass<> > reader(faiStream);
    CharString buffer;
    while (!atEnd(reader))
    {
        FaiIndexEntry_ entry;

        // Read REF_NAME.
        if (readUntilTabOrLineBreak(entry.name, reader) != 0)
            return 1;

        appendValue(index.refNameStore, entry.name);
        if (atEnd(reader) || value(reader) != '\t')
            return 1;  // Must be on tab.

        skipChar(reader, '\t');  // Must have been on tab, no checking.

        // Read SEQ_LENGTH.
        clear(buffer);
        if (readUntilTabOrLineBreak(buffer, reader) != 0)
            return 1;

        if (!lexicalCast2(entry.sequenceLength, buffer))
            return 1;  // Could not cast to integer.

        if (atEnd(reader) || value(reader) != '\t')
            return 1;  // Must be on tab.

        skipChar(reader, '\t');  // Must have been on tab, no checking.

        // Read OFFSET.
        clear(buffer);
        if (readUntilTabOrLineBreak(buffer, reader) != 0)
            return 1;

        if (!lexicalCast2(entry.offset, buffer))
            return 1;  // Could not cast to integer.

        if (atEnd(reader) || value(reader) != '\t')
            return 1;  // Must be on tab.

        skipChar(reader, '\t');  // Must have been on tab, no checking.

        // Read LINE_LENGTH.
        clear(buffer);
        if (readUntilTabOrLineBreak(buffer, reader) != 0)
            return 1;

        if (!lexicalCast2(entry.lineLength, buffer))
            return 1;  // Could not cast to integer.

        if (atEnd(reader) || value(reader) != '\t')
            return 1;  // Must be on tab.

        skipChar(reader, '\t');  // Must have been on tab, no checking.

        // Read OVERALL_LINE_LENGTH.
        clear(buffer);
        if (readUntilTabOrLineBreak(buffer, reader) != 0)
            return 1;

        if (!lexicalCast2(entry.overallLineLength, buffer))
            return 1;  // Could not cast to integer.

        if (!atEnd(reader) && value(reader) != '\r' && value(reader) != '\n')
            return 1;  // Must be on end of line or file.

        if (!atEnd(reader))
            skipLine(reader);  // Skip over line ending.

        appendValue(index.indexEntryStore, entry);
    }

    // Refresh name store cache.
    refresh(index.refNameStoreCache);

    return 0;
}

// TODO(holtgrew): The wrappers can go when the clash from read() from the file module is gone.

inline int read(FaiIndex & index, char * fastaFilename, char const * faiFilename)
{
    return read(index, static_cast<char const *>(fastaFilename), faiFilename);
}

inline int read(FaiIndex & index, char const * fastaFilename, char * faiFilename)
{
    return read(index, fastaFilename, static_cast<char const *>(faiFilename));
}

inline int read(FaiIndex & index, char * fastaFilename, char * faiFilename)
{
    return read(index, static_cast<char const *>(fastaFilename), static_cast<char const *>(faiFilename));
}

inline int read(FaiIndex & index, char const * fastaFilename)
{
    char buffer[1000];
    snprintf(buffer, 999, "%s.fai", toCString(fastaFilename));
    return read(index, fastaFilename, &buffer[0]);
}

inline int read(FaiIndex & index, char * fastaFilename)
{
    return read(index, static_cast<char const *>(fastaFilename));
}

inline int read(FaiIndex & index)
{
    // Cannot read if FAI filename is empty.
    if (empty(index.faiFilename))
        return 1;

    return read(index, toCString(index.fastaFilename), toCString(index.faiFilename));
}

// ---------------------------------------------------------------------------
// Function write()
// ---------------------------------------------------------------------------

/**
.Function.FaiIndex#write
..cat:Input/Output
..class:Class.FaiIndex
..signature:write(faiIndex[, faiFileName])
..summary:Write out an @Class.FaiIndex@ object.
..param.faiIndex:The @Class.FaiIndex@ object to write out.
...type:Class.FaiIndex
..param.faiFileName:The name of the FAI file to write to.
...remarks:This is optional only if the FAI index knows the FAI file name from a previous @Function.FaiIndex#build@ call.
...default:The FAI file name from the previous call to @Function.FaiIndex#build@, if any.
...type:nolink:$char const *$
..return:Status code $int$, $0$ indicating success and $1$ an error.
..include:seqan/seq_io.h
*/

// TODO(holtgrew): The wrappers can go when the clash from read() from the file module is gone.

inline int write(FaiIndex const & index, char const * faiFilename)
{
    // Open index files.
    std::ofstream indexOut(faiFilename, std::ios::binary | std::ios::out);

    for (unsigned i = 0; i < length(index.indexEntryStore); ++i)
    {
        FaiIndexEntry_ const & entry = index.indexEntryStore[i];
        indexOut << entry.name << '\t' << entry.sequenceLength << '\t' << entry.offset << '\t'
                 << entry.lineLength << '\t' << entry.overallLineLength << '\n';
    }

    return !indexOut.good();  // 1 on errors, 0 on success.
}

inline int write(FaiIndex const & index, char * faiFilename)
{
    return write(index, static_cast<char const *>(faiFilename));
}

inline int write(FaiIndex & index, char const * faiFilename)
{
    return write(static_cast<FaiIndex const &>(index), faiFilename);
}

inline int write(FaiIndex & index, char * faiFilename)
{
    return write(static_cast<FaiIndex const &>(index), static_cast<char const *>(faiFilename));
}

inline int write(FaiIndex const & index)
{
    if (empty(index.faiFilename))
        return 1;  // Cannot write out if faiFilename member is empty.
    return write(index, toCString(index.faiFilename));
}

inline int write(FaiIndex & index)
{
    return write(static_cast<FaiIndex const &>(index));
}

// ---------------------------------------------------------------------------
// Function build()
// ---------------------------------------------------------------------------

/**
.Function.FaiIndex#build
..class:Class.FaiIndex
..summary:Create @Class.FaiIndex@ from FASTA file.
..signature:build(faiIndex, fastaFilename[, faiFilename])
..description:The index can later be written out with @Function.FaiIndex#write@ and be loaded again using @Function.FaiIndex#read@.
..param.faiIndex:@Class.FaiIndex@ to write index to.
...type:Class.FaiIndex
..param.fastaFilename:Name of FASTA file to build an index for.
...type:Shortcut.CharString
..param.faiFilename:Name of FAI index file, stored in $faiIndex$. Optional.
...default:$fastaFilename + ".fai"$
...type:nolink:$char const *$
..returns:$int$, equal to 0 on success, != 0 otherwise.
..include:seqan/stream.h
 */

inline int build(FaiIndex & index, char const * seqFilename, char const * faiFilename)
{
    index.fastaFilename = seqFilename;
    index.faiFilename = faiFilename;
    
    if (index.mmapStringOpen)
        close(index.mmapString);
    if (!open(index.mmapString, toCString(seqFilename), OPEN_RDONLY))
        return 1;  // Could not open file.
    index.mmapStringOpen = true;

    typedef String<char, MMap<> > TMMapString;
    RecordReader<TMMapString, SinglePass<StringReader> > reader(index.mmapString);
    // Get file format, must be FASTA for FAI.
    AutoSeqStreamFormat tagSelector;
    if (!guessStreamFormat(reader, tagSelector))
        return 1;  // Invalid format.

    if (!isEqual(tagSelector, seqan::Fasta()))
        return 1;  // Invalid format, not FASTA.

    // Re-using the FASTA/FASTQ parsing code from read_fasta_fastq is not really feasible here.  We roll our own
    // mini-parser from scratch.
    CharString line;
    CharString seqName;
    __uint32 seqLength = 0;
    __uint64 seqOffset = 0;
    __uint32 lineLength = 0;
    __uint32 lineSize = 0;
    while (!atEnd(reader))
    {
        clear(line);
        clear(seqName);

        if (value(reader) != '>')
            return 1;  // Must be >.

        goNext(reader);

        int res = readUntilWhitespace(seqName, reader);
        if (res != 0)
            return res;  // Error reading.

        res = skipLine(reader);
        if (res != 0)
            return res;  // Error reading.

        seqOffset = reader._current - begin(reader._string, Standard());

        res = readLine(line, reader);
        if (res != 0 && res != EOF_BEFORE_SUCCESS)
            return res;  // Error reading.

        lineSize = reader._current - begin(reader._string, Standard()) - seqOffset;
        lineLength = length(line);
        seqLength = lineLength;

        while (!atEnd(reader))
        {
            char c = value(reader);
            if (c == '>')
                break;
            if (!isspace(c))
                seqLength += 1;
            goNext(reader);
        }

        FaiIndexEntry_ entry;
        entry.name = seqName;
        entry.sequenceLength = seqLength;
        entry.offset = seqOffset;
        entry.lineLength = lineLength;
        entry.overallLineLength = lineSize;
        appendValue(index.indexEntryStore, entry);
        appendValue(index.refNameStore, entry.name);
    }

    // Refresh name store cache.
    refresh(index.refNameStoreCache);

    return 0;
}

inline int build(FaiIndex & index, char const * seqFilename)
{
    CharString faiFilename(seqFilename);
    append(faiFilename, ".fai");
    return build(index, seqFilename, toCString(faiFilename));
}

}  // namespace seqan

#endif  // #ifndef SEQAN_CORE_APPS_RABEMA_FAI_INDEX_H_
