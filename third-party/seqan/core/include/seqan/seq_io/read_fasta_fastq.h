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
// Author: Hannes Hauswedell <hauswedell@mi.fu-berlin.de>
// ==========================================================================
// Record and Document Reading for FASTA and FASTQ files.
// ==========================================================================

//TODO(h4nn3s): double-check if we really want to allow EOF inside meta-line
// and also if a fastq file is legal if a record contains no qualities

#ifndef SEQAN_SEQ_IO_READ_FASTA_FASTQ_H_
#define SEQAN_SEQ_IO_READ_FASTA_FASTQ_H_

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

/**
.Tag.File Format.tag.Fasta:
    FASTA file format for sequences.
..include:seqan/file.h
*/
struct TagFasta_;
typedef Tag<TagFasta_> Fasta;

/**
.Tag.File Format.tag.Fastq:
    FASTQ file format for sequences.
..include:seqan/file.h
*/
struct TagFastq_;
typedef Tag<TagFastq_> Fastq;


// ============================================================================
// Metafunctions
// ============================================================================

// Returns character starting a meta field.

template <typename TTag>
struct MetaFirstChar_;

template <>
struct MetaFirstChar_<Tag<TagFasta_> >
{
    static const char VALUE = '>';
};

template <>
struct MetaFirstChar_<Tag<TagFastq_> >
{
    static const char VALUE = '@';
};

// Returns character starting a field after the sequence field.  In the case of FASTA, this is '>' since a new meta
// field will be created.  In the case of FASTQ, this is '+'.

template <typename TTag>
struct AfterSeqFirstChar_;

template <>
struct AfterSeqFirstChar_<Tag<TagFasta_> >
{
    static const char VALUE = '>';
};

template <>
struct AfterSeqFirstChar_<Tag<TagFastq_> >
{
    static const char VALUE = '+';
};

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function _clearAndReserveMemory() and helpers
// ----------------------------------------------------------------------------

// When not reading char values, we can simply use _countHelper() since the separators '>' and '+' are not part of the
// sequence alphabet.

template <typename TFormatTag, typename TAlph, typename TRecordReader>
inline int
_countSequenceFastAQ(unsigned int & count,
                     TRecordReader & reader,
                     TFormatTag const & /* formatTag */,
                     TAlph const & /* Alphabet type */)
{
    return _countHelper(count, reader, Tag<TAlph>(), Whitespace_(), false);
}

// When reading char values, we have to use a more complicated algorithm:  Read until the separator char '>'/'+' begins
// a line.  Note that when reading text including these characters then there is ambiguity.  When we find such a
// character then we greedily decide that this does not belong to the sequence.

template <typename TFormatTag, typename TRecordReader>
inline int
_countSequenceFastAQ(unsigned int & count,
                     TRecordReader & reader,
                     TFormatTag const & /* formatTag */,
                     char const & /* Alphabet type */)
{
    // The variable afterEol is true if we are after an EOL char and we are not at the beginning of the sequence
    // record.  Empty sequence lines are fine, non-existent ones are not.
    bool afterEol = false;

    while (!atEnd(reader))
    {
        char c = value(reader);
        if (c == '\r' || c == '\n')
        {
            afterEol = true;
            goNext(reader);
            if (resultCode(reader) != 0)
                return resultCode(reader);
            continue;
        }

        if (afterEol && c == AfterSeqFirstChar_<TFormatTag>::VALUE)
        {
            return 0;  // Done, at stop char.
        }

        if (!isspace(c))
            count += 1;
        afterEol = false;
        goNext(reader);
        if (resultCode(reader) != 0)
            return resultCode(reader);
    }

    return EOF_BEFORE_SUCCESS;
}

template <typename TSeqAlph,
          typename TFile,
          typename TSpec,
          typename TTag>
inline int
_countMetaAndSequence(unsigned int & metaLength,
                      unsigned int & seqLength,
                      RecordReader<TFile, DoublePass<TSpec> > & reader,
                      TTag const & formatTag,
                      TSeqAlph const & /* tag*/)
{
    metaLength=0;
    seqLength=0;

    // COUNT META
    if (atEnd(reader) || value(reader) != MetaFirstChar_<TTag>::VALUE)
        return RecordReader<TFile, DoublePass<TSpec> >::INVALID_FORMAT;
    goNext(reader);
    if (resultCode(reader))
        return resultCode(reader);
    if (atEnd(reader)) // empty ID, no sequence, this is legal TODO?
        return 0;

    int res = countLine(metaLength, reader);
    if (res == EOF_BEFORE_SUCCESS)  // EOF half way in ID is legal
        return 0;
    else if (res)
        return res;

    if (atEnd(reader)) // no sequence
        return 0;

    // COUNT SEQUENCE
    res = _countSequenceFastAQ(seqLength, reader, formatTag, TSeqAlph());
    if (res == EOF_BEFORE_SUCCESS)  // EOF half way in sequence is legal
        return 0;
    else if (res)
        return res;

    return 0;
}

// SINGLE-Pass
template <typename TIdString,
          typename TSeqString,
          typename TFile,
          typename TSpec,
          typename TTag>
inline int
_clearAndReserveMemory(TIdString & meta, TSeqString & seq,
                       RecordReader<TFile, SinglePass<TSpec> > & /**/,
                       TTag const & /*tag*/)
{
    clear(meta);
    clear(seq);
    return 0;
}

// DOUBLE-Pass
template <typename TIdString,
          typename TSeqString,
          typename TFile,
          typename TSpec,
          typename TTag>
inline int
_clearAndReserveMemory(TIdString & meta, TSeqString & seq,
                       RecordReader<TFile, DoublePass<TSpec> > & reader,
                       TTag const & /*tag*/)
{
    clear(meta);
    clear(seq);
    startFirstPass(reader);

    unsigned int metaLength=0;
    unsigned int seqLength=0;
    // COUNT
    int res = _countMetaAndSequence(metaLength,
                                    seqLength,
                                    reader,
                                    TTag(),
                                    typename Value<TSeqString>::Type());
    if ((res != 0) && (res != EOF_BEFORE_SUCCESS))
        return res;

    // RESERVE FOR META
    reserve(meta, metaLength, Exact());

    // RESERVE FOR SEQUENCE
    reserve(seq, seqLength, Exact());

    startSecondPass(reader);
    return 0;
}

template <typename TIdString,
          typename TSeqString,
          typename TQualString,
          typename TFile,
          typename TSpec,
          typename TTag>
inline int
_clearAndReserveMemory(TIdString & meta, TSeqString & seq, TQualString & qual,
                       RecordReader<TFile, TSpec> & reader,
                       TTag const & tag)
{
    int res = _clearAndReserveMemory(meta, seq, reader, tag);
    clear(qual);
    reserve(qual, capacity(seq), Exact());
    return res;
}

// ----------------------------------------------------------------------------
// Function _readMetaAndSequence() and helpers
// ----------------------------------------------------------------------------

// This overload of _readSequenceFastAQ() assumes that the separator characters, such as '>'/'+' are not part of
// the alphabet.  This way, we can use _readHelper() and do not have to fall back to the more complicated way
// used in the case of char.

template <typename TAlph, typename TSpec, typename TRecordReader, typename TFormatTag>
inline int
_readSequenceFastAQ(String<TAlph, TSpec> & string,
                    TRecordReader & reader,
                    TFormatTag const & /*formatTag*/)
{
    return _readHelper(string, reader, Tag<TAlph>(), Whitespace_(), false);
}

// If we want to read char values from a FASTA or FASTQ file then we have to fall back to a more complicated
// algorithm, this is similar to the corresponding variant of _countSequenceFastAQ().
//
// TODO(holtgrew): Very similar to the _countSequenceFastAQ() function, can we make improve this using functors?

template <typename TString, typename TRecordReader, typename TFormatTag>
inline int
_readSequenceFastAQCharImpl(TString & string,
                            TRecordReader & reader,
                            TFormatTag const & /*formatTag*/)
{
    // The variable afterEol is true if we are after an EOL char and we are not at the beginning of the sequence
    // record.  Empty sequence lines are fine, non-existent ones are not.
    bool afterEol = false;

    while (!atEnd(reader))
    {
        char c = value(reader);
        if (c == '\r' || c == '\n')
        {
            afterEol = true;
            goNext(reader);
            if (resultCode(reader) != 0)
                return resultCode(reader);
            continue;
        }

        if (afterEol && c == AfterSeqFirstChar_<TFormatTag>::VALUE)
        {
            return 0;  // Done, at stop char.
        }

        if (!isspace(c))
            appendValue(string, c);
        afterEol = false;
        goNext(reader);
        if (resultCode(reader) != 0)
            return resultCode(reader);
    }

    return EOF_BEFORE_SUCCESS;
}

template <typename TSpec, typename TRecordReader, typename TFormatTag>
inline int
_readSequenceFastAQ(String<char, TSpec> & string,
                    TRecordReader & reader,
                    TFormatTag const & formatTag)
{
    return _readSequenceFastAQCharImpl(string, reader, formatTag);
}

template <typename TString, typename TRecordReader, typename TFormatTag>
inline int
_readSequenceFastAQ(TString & string,
                    TRecordReader & reader,
                    TFormatTag const & formatTag)
{
    return _readSequenceFastAQCharImpl(string, reader, formatTag);
}

// This reads Meta and Sequence
template <typename TIdString,
          typename TSeqString,
          typename TFile,
          typename TPass,
          typename TTag>
inline int
_readMetaAndSequence(TIdString & meta, TSeqString & seq,
                     RecordReader<TFile, TPass > & reader,
                     TTag const & formatTag)
{
    // READ META
    if (atEnd(reader) || value(reader) != MetaFirstChar_<TTag>::VALUE)
        return RecordReader<TFile, TPass>::INVALID_FORMAT;
    goNext(reader);
    if (resultCode(reader))
        return resultCode(reader);
    if (atEnd(reader)) // empty ID, no sequence, this is legal
        return 0;

    int res = readLine(meta, reader);
    if (res == EOF_BEFORE_SUCCESS)
        return EOF_BEFORE_SUCCESS;
    else if (res)
        return res;

    if (atEnd(reader)) // no sequence
        return 0;

    // READ SEQUENCE
    res = _readSequenceFastAQ(seq, reader, formatTag);
    if (res && res != EOF_BEFORE_SUCCESS)
        return res;

    return 0;
}

// ----------------------------------------------------------------------------
// Function _skipQualityBlock() and _readQualityBlock()
// ----------------------------------------------------------------------------

template <typename TFile, typename TPass>
inline int
_skipQualityBlock(RecordReader<TFile, TPass > & /**/,
                  unsigned const /**/,
                  Fasta const & /*tag*/)
{
    // NOOP for Fasta
    return 0;
}

template <typename TFile, typename TPass>
inline int
_skipQualityBlock(RecordReader<TFile, TPass > & reader,
                  unsigned const seqLength,
                  Fastq const & /*tag*/)
{
    int res = 0;
    // SKIP QUALITIES' META
    skipLine(reader);
    if (res == EOF_BEFORE_SUCCESS)  // EOF half way in ID is legal
        return 0;
    else if (res)
        return res;

    // SKIP QUALITIES
    res = skipNCharsIgnoringWhitespace(reader, seqLength); // there have to be n qualities
    if (res && res == EOF_BEFORE_SUCCESS)
        return EOF_BEFORE_SUCCESS;
    else if (res)
        return RecordReader<TFile, TPass >::INVALID_FORMAT;
    skipLine(reader); // goto to next line if it exists, result is unimportant

    return 0;
}

// Reading of the quality block.
template <typename TSeqString, typename TFile, typename TPass>
inline int _readQualityBlockHelper(TSeqString & seq,
                                   RecordReader<TFile, TPass > & reader,
                                   unsigned const seqLength,
                                   True const & /*assign qualities to seq*/)
{
    // Copy of readNCharsIgnoringWhitespace but we assign to seq instead of appending
    // to a string.
    typedef typename Iterator<TSeqString, Rooted>::Type TSeqIter;
    TSeqIter it = begin(seq, Rooted());
    
    for (unsigned i = 0; i < seqLength; ++i)
    {
        if (atEnd(reader))
            return EOF_BEFORE_SUCCESS;

        if (_charCompare(value(reader), Whitespace_()))
        {
            --i;
        }
        else
        {
            SEQAN_ASSERT(!atEnd(it));
            assignQualityValue(*it, value(reader));
            goNext(it);
        }

        goNext(reader);
        if (resultCode(reader) != 0)
            return resultCode(reader);
    }
    return 0;
}

template <typename TQualString, typename TFile, typename TPass>
inline int _readQualityBlockHelper(TQualString & qual,
                                   RecordReader<TFile, TPass > & reader,
                                   unsigned const seqLength,
                                   False const & /*qual is character sequence*/)
{
    reserve(qual, seqLength, Exact());
    return readNCharsIgnoringWhitespace(qual, reader, seqLength);
}

template <typename TIdString,
          typename TQualString,
          typename TFile,
          typename TPass>
inline int
_readQualityBlock(TQualString & qual,
                  RecordReader<TFile, TPass > & reader,
                  unsigned const seqLength,
                  TIdString const & meta,
                  Fastq const & /*tag*/)
{
    // READ AND CHECK QUALITIES' META
    if (atEnd(reader))
        return EOF_BEFORE_SUCCESS;
    if (value(reader) != '+')
        return RecordReader<TFile, TPass >::INVALID_FORMAT;
    goNext(reader);
    if (resultCode(reader))
        return resultCode(reader);
    if (atEnd(reader)) // empty ID, no sequence, this is legal? TODO
        return 0;

    CharString qualmeta_buffer;
    int res = readLine(qualmeta_buffer, reader);
    if (res && res == EOF_BEFORE_SUCCESS)
        return EOF_BEFORE_SUCCESS;
    else if (res)
        return RecordReader<TFile, TPass >::INVALID_FORMAT;

    // meta string has to be empty or identical to sequence's meta
    if ((qualmeta_buffer != "") && (qualmeta_buffer != meta))
        return RecordReader<TFile, TPass >::INVALID_FORMAT;

    if (atEnd(reader)) // empty qualities, is this legal?
        return 0;

    // READ QUALITIES
    res = _readQualityBlockHelper(qual, reader, seqLength, typename HasQualities<typename Value<TQualString>::Type>::Type());
    // there have to be n qualities
    if (res && res == EOF_BEFORE_SUCCESS)
        return EOF_BEFORE_SUCCESS;
    else if (res)
        return RecordReader<TFile, TPass >::INVALID_FORMAT;
    skipLine(reader); // goto to next line if it exists, result is unimportant

    return 0;
}

// ----------------------------------------------------------------------------
// Function readRecord()                               [Single and double pass]
// ----------------------------------------------------------------------------
/**
.Function.readRecord
..signature:readRecord(TIdString & meta, TSeqString & seq, TRecordReader & reader, Fasta const &)
..remarks:For FASTA-Files a double-Pass implementation of RecordReader is implemented, which offers better performance. Just pass a Double-Pass reader object (only works with seekable Streams).
*/

/**
.Function.readRecord
..signature:readRecord(TIdString & meta, TSeqString & seq, TRecordReader & reader, Fastq const &)
*/

template <typename TIdString,
          typename TSeqString,
          typename TFile,
          typename TPass>
inline int
readRecord(TIdString & meta,
           TSeqString & seq,
           RecordReader<TFile, TPass > & reader,
           Fasta const & /*tag*/)
{
    int res = _clearAndReserveMemory(meta, seq, reader, Fasta());
    if (res)
        return res;

    res = _readMetaAndSequence(meta, seq, reader, Fasta());
    if (res)
        return res;

    return _skipQualityBlock(reader, length(seq), Fasta());
}

template <typename TIdString,
          typename TSeqString,
          typename TQualString,
          typename TFile,
          typename TPass>
inline int
readRecord(TIdString & meta,
           TSeqString & seq,
           TQualString & quals,
           RecordReader<TFile, TPass > & reader,
           Fasta const & /*tag*/)
{
    int res = _clearAndReserveMemory(meta, seq, reader, Fasta());
    if (res)
        return res;

    clear(quals);

    res = _readMetaAndSequence(meta, seq, reader, Fasta());
    if (res)
        return res;

    return _skipQualityBlock(reader, length(seq), Fasta());
}

template <typename TIdString,
          typename TSeqString,
          typename TFile,
          typename TPass>
inline int
readRecord(TIdString & meta,
           TSeqString & seq,
           RecordReader<TFile, TPass > & reader,
           Fastq const & /*tag*/)
{
    int res = _clearAndReserveMemory(meta, seq, reader, Fastq());
    if (res)
        return res;

    res = _readMetaAndSequence(meta, seq, reader, Fastq());
    if (res)
        return res;

    if (HasQualities<typename Value<TSeqString>::Type>::VALUE)
        return _readQualityBlock(seq, reader, length(seq), meta, Fastq());
    else
        return _skipQualityBlock(reader, length(seq), Fastq());
}


/**
.Function.readRecord
..signature:readRecord(TIdString & meta, TSeqString & seq, TQualString & qual, TRecordReader & reader, Fastq const &)
..remarks:For FASTQ-Files a double-Pass implementation of RecordReader is implemented, which offers better performance. Just pass a Double-Pass reader object (only works with seekable Streams).
*/

// FASTQ and we want explicit qualities
template <typename TIdString,
          typename TSeqString,
          typename TQualString,
          typename TFile,
          typename TPass>
inline int
readRecord(TIdString & meta,
           TSeqString & seq,
           TQualString & qual,
           RecordReader<TFile, TPass > & reader,
           Fastq const & /*tag*/)
{
    int res = _clearAndReserveMemory(meta, seq, qual, reader, Fastq());
    if (res)
        return res;

    res = _readMetaAndSequence(meta, seq, reader, Fastq());
    if (res)
        return res;

    res = _readQualityBlock(qual, reader, length(seq), meta, Fastq());
    if (res == 0 || res == EOF_BEFORE_SUCCESS)
        return 0;
    return res;
}

// ----------------------------------------------------------------------------
// Function read();  Double-Pass.
// ----------------------------------------------------------------------------

// Reads a whole FASTA/FASTQ file into string sets, optimizing memory usage.
// optimized for ConcatDirect StringSets

template <typename TIdString, typename TSeqString, typename TQualString, typename TRecordReader>
inline int
_readFastAQQualityReadDispatcher(StringSet<TIdString, Owner<ConcatDirect<> > > & sequenceIds,
                                 StringSet<TSeqString, Owner<ConcatDirect<> > > & sequences,
                                 StringSet<TQualString, Owner<ConcatDirect<> > > & /*qualities*/,
                                 TRecordReader & reader,
                                 unsigned i,
                                 True const & /*seq has qualities*/)
{
    typedef StringSet<TSeqString, Owner<ConcatDirect<> > > TSequenceSet;
    typedef typename Concatenator<TSequenceSet>::Type TConcat;
    typedef typename Suffix<TConcat>::Type TSuffix;
    TSuffix s(concat(sequences), sequences.limits[i]);
    return _readQualityBlock(s, reader, length(sequences[i]), sequenceIds[i], Fastq());
}

template <typename TIdString, typename TSeqString, typename TQualString, typename TRecordReader>
inline int
_readFastAQQualityReadDispatcher(StringSet<TIdString, Owner<ConcatDirect<> > > & sequenceIds,
                                 StringSet<TSeqString, Owner<ConcatDirect<> > > & sequences,
                                 StringSet<TQualString, Owner<ConcatDirect<> > > & qualities,
                                 TRecordReader & reader,
                                 unsigned i,
                                 False const & /*qualities are explicit*/)
{
    return _readQualityBlock(qualities.concat, reader, length(sequences[i]), sequenceIds[i], Fastq());
}

template <typename TIdString,
          typename TSeqString,
          typename TQualString,
          typename TFile,
          typename TSpec,
          typename TTag>
int _readFastAQ(StringSet<TIdString, Owner<ConcatDirect<> > > & sequenceIds,
                StringSet<TSeqString, Owner<ConcatDirect<> > > & sequences,
                StringSet<TQualString, Owner<ConcatDirect<> > > & qualities,
                RecordReader<TFile, DoublePass<TSpec> > & reader,
                bool const withQual,
                TTag const & /*tag*/)
{
    int res = 0;
    String<unsigned> metaLengths;
    String<unsigned> seqLengths;

    // ------------------------------------------------------------------------
    // First Pass: Compute meta and sequence lengths.
    // ------------------------------------------------------------------------
    startFirstPass(reader);
    size_t sequenceCount = 0;
    while (!atEnd(reader))
    {
        sequenceCount += 1;
        unsigned metaLength = 0;
        unsigned seqLength = 0;
        res = _countMetaAndSequence(metaLength,
                                    seqLength,
                                    reader,
                                    TTag(),
                                    typename Value<TSeqString>::Type());
        if (res)
            return res;

        appendValue(metaLengths, metaLength, Generous());
        appendValue(seqLengths, seqLength, Generous());

        res = _skipQualityBlock(reader, seqLength, TTag());
        if (res)
            return res;
    }

    // ------------------------------------------------------------------------
    // Allocate memory.
    // ------------------------------------------------------------------------
    clear(sequenceIds);
    clear(sequences);
    if (withQual && !HasQualities<typename Value<TSeqString>::Type>::VALUE)
        clear(qualities);

    resize(sequenceIds.limits, sequenceCount + 1, Exact());
    resize(sequences.limits, sequenceCount + 1, Exact());
    sequenceIds.limits[0] = 0;
    sequences.limits[0] = 0;

    unsigned long metaLengthsSum = 0;
    unsigned long seqLengthsSum = 0;

    for (unsigned int i = 0; i < sequenceCount; ++i)
    {
        metaLengthsSum += metaLengths[i];
        sequenceIds.limits[i+1] = metaLengthsSum;
        seqLengthsSum += seqLengths[i];
        sequences.limits[i+1] = seqLengthsSum;
    }

    reserve(sequenceIds.concat, metaLengthsSum + 1, Exact());
    reserve(sequences.concat, seqLengthsSum + 1, Exact());
    if (withQual && !HasQualities<typename Value<TSeqString>::Type>::VALUE)
    {
        assign(qualities.limits, sequences.limits);
        reserve(qualities.concat, seqLengthsSum + 1, Exact());
    }

    // ------------------------------------------------------------------------
    // Second Pass: Actually read data.
    // ------------------------------------------------------------------------
    startSecondPass(reader);
    for (unsigned int i = 0; i < sequenceCount; ++i)
    {
        res = _readMetaAndSequence(sequenceIds.concat,
                                   sequences.concat,
                                   reader,
                                   TTag());
        switch(res)
        {
            case 0:
                break;
            case EOF_BEFORE_SUCCESS:        // file may end without newline
                if (i >= sequenceCount -1)
                    break;
            default:
                return res;
        }
        if (withQual)
        {
            _readFastAQQualityReadDispatcher(sequenceIds, sequences, qualities, reader, i, typename HasQualities<typename Value<TSeqString>::Type>::Type());
        }
        else
        {
            res = _skipQualityBlock(reader, length(sequences[i]), TTag());
        }
        if (res)
            return res;
    }

    return 0;
}

// Reads a whole FASTA/FASTQ file into string sets, optimizing memory usage.
// Generic, non-concat version
template <typename TIdString, typename TIdSpec,
          typename TSeqString, typename TSeqSpec,
          typename TQualString, typename TQualSpec,
          typename TFile,
          typename TSpec,
          typename TTag>
int _readFastAQ(StringSet<TIdString, TIdSpec> & sequenceIds,
                StringSet<TSeqString, TSeqSpec> & sequences,
                StringSet<TQualString, TQualSpec> & qualities,
                RecordReader<TFile, DoublePass<TSpec> > & reader,
                bool const withQual,
                TTag const & /*tag*/)
{
    int res = 0;
    String<unsigned> metaLengths;
    String<unsigned> seqLengths;

    // ------------------------------------------------------------------------
    // First Pass: Compute meta and sequence lengths.
    // ------------------------------------------------------------------------
    startFirstPass(reader);
    size_t sequenceCount = 0;
    while (!atEnd(reader))
    {
        sequenceCount += 1;
        unsigned metaLength = 0;
        unsigned seqLength = 0;
        res = _countMetaAndSequence(metaLength,
                                    seqLength,
                                    reader,
                                    TTag(),
                                    typename Value<TSeqString>::Type());
        if (res)
            return res;

        appendValue(metaLengths, metaLength, Generous());
        appendValue(seqLengths, seqLength, Generous());

        res = _skipQualityBlock(reader, seqLength, TTag());
        if (res)
            return res;
    }

    // ------------------------------------------------------------------------
    // Allocate memory.
    // ------------------------------------------------------------------------
    clear(sequenceIds);
    clear(sequences);
    if (withQual && !HasQualities<typename Value<TSeqString>::Type>::VALUE)
        clear(qualities);

    resize(sequenceIds, sequenceCount, Exact());
    resize(sequences, sequenceCount, Exact());
    if (withQual && !HasQualities<typename Value<TSeqString>::Type>::VALUE)
        resize(qualities, sequenceCount, Exact());

    for (unsigned int i = 0; i < sequenceCount; ++i)
    {
        reserve(sequenceIds[i], metaLengths[i], Exact());
        reserve(sequences[i], seqLengths[i], Exact());
        if (withQual && !HasQualities<typename Value<TSeqString>::Type>::VALUE)
            reserve(qualities[i], seqLengths[i], Exact());
    }

    // ------------------------------------------------------------------------
    // Second Pass: Actually read data.
    // ------------------------------------------------------------------------
    startSecondPass(reader);
    for (unsigned int i = 0; i < sequenceCount; ++i)
    {
        res = _readMetaAndSequence(sequenceIds[i],
                                   sequences[i],
                                   reader,
                                   TTag());
        switch(res)
        {
            case 0:
                break;
            case EOF_BEFORE_SUCCESS:        // file may end without newline
                if (i >= sequenceCount -1)
                    break;
            default:
                return res;
        }
        if (withQual)
        {
            if (HasQualities<typename Value<TSeqString>::Type>::VALUE)
                res = _readQualityBlock(sequences[i], reader, length(sequences[i]), sequenceIds[i], Fastq());
            else
                res = _readQualityBlock(qualities[i], reader, length(sequences[i]), sequenceIds[i], Fastq());
        }
        else
        {
            res = _skipQualityBlock(reader, length(sequences[i]), TTag());
        }
        if (res)
            return res;
    }
    return 0;
}

/**
.Function.read2
..signature:read2(StringSet<TIdString, TIdSpec> & sequenceIds, StringSet<TSeqString, TSeqSpec> & sequences, RecordReader<TFile, DoublePass<TSpec> > & reader, Fasta const &)
*/

// FASTA
template <typename TIdString, typename TIdSpec,
          typename TSeqString, typename TSeqSpec,
          typename TFile,
          typename TSpec>
int read2(StringSet<TIdString, TIdSpec> & sequenceIds,
         StringSet<TSeqString, TSeqSpec> & sequences,
         RecordReader<TFile, DoublePass<TSpec> > & reader,
         Fasta const & /*tag*/)
{
    StringSet<CharString, TSeqSpec> qualities;
    return _readFastAQ(sequenceIds, sequences, qualities, reader, false, Fasta());
}

template <typename TIdString, typename TIdSpec,
          typename TSeqString, typename TSeqSpec,
          typename TQualString, typename TQualSpec,
          typename TFile,
          typename TSpec>
int read2(StringSet<TIdString, TIdSpec> & sequenceIds,
         StringSet<TSeqString, TSeqSpec> & sequences,
         StringSet<TQualString, TQualSpec> & qualities,
         RecordReader<TFile, DoublePass<TSpec> > & reader,
         Fasta const & /*tag*/)
{
    int res = _readFastAQ(sequenceIds, sequences, qualities, reader, false, Fasta());
    clear(qualities);
    resize(qualities, length(sequences));
    return res;
}

/**
.Function.read2
..signature:read2(StringSet<TIdString, TIdSpec> & sequenceIds, StringSet<TSeqString, TSeqSpec> & sequences, RecordReader<TFile, DoublePass<TSpec> > & reader, Fastq const &)
*/

// FASTQ, if we don't have explicit qualities.
template <typename TIdString, typename TIdSpec,
          typename TSeqString, typename TSeqSpec,
          typename TFile,
          typename TSpec>
int read2(StringSet<TIdString, TIdSpec> & sequenceIds,
         StringSet<TSeqString, TSeqSpec> & sequences,
         RecordReader<TFile, DoublePass<TSpec> > & reader,
         Fastq const & /*tag*/)
{
    typedef typename Value<TSeqString>::Type TChar;
    StringSet<CharString, TSeqSpec> qualities;
    return _readFastAQ(sequenceIds, sequences, qualities, reader, HasQualities<TChar>::VALUE, Fastq());
}

/**
.Function.read2
..signature:read2(StringSet<TIdString, TIdSpec> & sequenceIds, StringSet<TSeqString, TSeqSpec> & sequences, StringSet<TQualString, TQualSpec> & qualities, RecordReader<TFile, DoublePass<TSpec> > & reader, Fastq const &)
*/

// FASTQ and we want Qualities
template <typename TIdString, typename TIdSpec,
          typename TSeqString, typename TSeqSpec,
          typename TQualString, typename TQualSpec,
          typename TFile,
          typename TSpec>
int read2(StringSet<TIdString, TIdSpec> & sequenceIds,
         StringSet<TSeqString, TSeqSpec> & sequences,
         StringSet<TQualString, TQualSpec> & qualities,
         RecordReader<TFile, DoublePass<TSpec> > & reader,
         Fastq const & /*tag*/)
{
    return _readFastAQ(sequenceIds, sequences, qualities, reader, true, Fastq());
}

// ----------------------------------------------------------------------------
// Function read();  Single-Pass
// ----------------------------------------------------------------------------

template <typename TIdString, typename TIdSpec,
          typename TSeqString, typename TSeqSpec,
          typename TFile,
          typename TSpec, typename TTag>
typename EnableIf<typename Or<typename IsSameType<typename RemoveConst_<TTag>::Type, typename RemoveConst_<Fasta>::Type >::Type,
                              typename IsSameType<typename RemoveConst_<TTag>::Type, typename RemoveConst_<Fastq>::Type >::Type>::Type,
                   int>::Type
read2(StringSet<TIdString, TIdSpec> & sequenceIds,
      StringSet<TSeqString, TSeqSpec> & sequences,
      RecordReader<TFile, SinglePass<TSpec> > & reader,
      TTag const & tag)
{
    TIdString id;
    TSeqString seq;
    while (!atEnd(reader))
    {
        int res = readRecord(id, seq, reader, tag);
        if (res != 0)
            return res;
        appendValue(sequenceIds, id);
        appendValue(sequences, seq);
    }
    return 0;
}

template <typename TIdString, typename TIdSpec,
          typename TSeqString, typename TSeqSpec,
          typename TQualString, typename TQualSpec,
          typename TFile,
          typename TSpec, typename TTag>
int read2(StringSet<TIdString, TIdSpec> & sequenceIds,
          StringSet<TSeqString, TSeqSpec> & sequences,
          StringSet<TQualString, TQualSpec> & qualities,
          RecordReader<TFile, SinglePass<TSpec> > & reader,
          TTag const & tag)
{
    clear(qualities);
    return read2(sequenceIds, sequences, reader, tag);
}

template <typename TIdString, typename TIdSpec,
          typename TSeqString, typename TSeqSpec,
          typename TQualString, typename TQualSpec,
          typename TFile,
          typename TSpec>
int read2(StringSet<TIdString, TIdSpec> & sequenceIds,
         StringSet<TSeqString, TSeqSpec> & sequences,
         StringSet<TQualString, TQualSpec> & qualities,
         RecordReader<TFile, SinglePass<TSpec> > & reader,
         Fastq const & tag)
{
    TIdString id;
    TSeqString seq;
    TQualString qual;
    while (!atEnd(reader))
    {
        int res = readRecord(id, seq, qual, reader, tag);
        if (res != 0)
            return res;
        appendValue(sequenceIds, id);
        appendValue(sequences, seq);
        appendValue(qualities, qual);
    }
    return 0;
}

// ----------------------------------------------------------------------------
// Function guessStreamFormat()
// ----------------------------------------------------------------------------

/**
.Function.guessStreamFormat
..signature:guessStreamFormat(reader, tag)
..returns:$bool$, indicating whether the file of $reader$ is of the given format.
..param.reader:RecordReader to query.
...type:Class.RecordReader
..param.tag:Format indicator
...type:Tag.File Format
...type:nolink:$Fasta$
...type:nolink:$Fastq$
..include:seqan/stream.h
*/

template <typename TStream, typename TPass>
inline bool
guessStreamFormat(RecordReader<TStream, TPass> & reader, Fasta const & /*tag*/)
{
    LimitRecordReaderInScope<TStream, TPass> limiter(reader);
    while (!atEnd(reader))
    {
        CharString meta;
        CharString seq;
        int r = readRecord(meta, seq, reader, Fasta());
//         std::cout << "Meta: " << toCString(meta) << "\nSeq: " << toCString(seq) << "\n";
        if (r == RecordReader<TStream, TPass>::INVALID_FORMAT)
            return false;
    }
    return true;
}

template <typename TStream, typename TPass>
inline bool
guessStreamFormat(RecordReader<TStream, TPass> & reader, Fastq const & /*tag*/)
{
    LimitRecordReaderInScope<TStream, TPass> limiter(reader);
    while (!atEnd(reader))
    {
        CharString meta;
        CharString seq;
        CharString qual;
        int r = readRecord(meta, seq, qual, reader, Fastq());
        if (r == RecordReader<TStream, TPass>::INVALID_FORMAT)
            return false;
    }
    return true;
}

}  // namespace seqan

#endif  // #ifndef SEQAN_SEQ_IO_READ_FASTA_FASTQ_H_
