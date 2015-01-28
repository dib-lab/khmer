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
// This header contains overloads for read2() and readRecord() that allow
// them to be used for automatic format detection with sequence files.
// ==========================================================================

// (weese:)
// Why did you introduce int i in every function?
// That adds complexity and requires an additional _x() helper function
// for each x() function. Please have a look at assignSeq() in file/file_format_mmap.h

#ifndef SEQAN_STREAM_RECORD_READER_AUTO_FORMAT_H_
#define SEQAN_STREAM_RECORD_READER_AUTO_FORMAT_H_

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

enum ReadAutoFormatErrorCodes_
{
    IOERR_READ_AUTO_FORMAT_UNKNOWN_FORMAT = 2049
};

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function read2
// ----------------------------------------------------------------------------

// For read2(metas, seqs, quals, reader, TTagSelector/AutoSeqStreamFormat).

template <typename TIdString, typename TIdSpec, typename TSeqString, typename TSeqSpec, typename TQualString,
          typename TQualSpec, typename TRecordReader>
int _read2(StringSet<TIdString, TIdSpec> & /*metas*/,
           StringSet<TSeqString, TSeqSpec> & /*seqs*/,
           StringSet<TQualString, TQualSpec> & /*quals*/,
           TRecordReader & /*reader*/,
           TagSelector<void> const & /*tagSelector*/,
           int i)
{
    SEQAN_ASSERT_EQ(i, -1);
    (void)i;  // Only used for assertion.
    return IOERR_READ_AUTO_FORMAT_UNKNOWN_FORMAT;
}

template <typename TIdString, typename TIdSpec, typename TSeqString, typename TSeqSpec, typename TQualString,
          typename TQualSpec, typename TRecordReader, typename TTagList>
int _read2(StringSet<TIdString, TIdSpec> & metas,
           StringSet<TSeqString, TSeqSpec> & seqs,
           StringSet<TQualString, TQualSpec> & quals,
           TRecordReader & reader,
           TagSelector<TTagList> const & tagSelector,
           int i)
{
    if (i == -1)
        return IOERR_READ_AUTO_FORMAT_UNKNOWN_FORMAT;
    if (i != tagSelector.tagId)
        return _read2(metas, seqs, quals, reader, static_cast<typename TagSelector<TTagList>::Base const &>(tagSelector), i - 1);

    typedef typename TTagList::Type TFormat;
    return read2(metas, seqs, quals, reader, TFormat());
}

template <typename TIdString, typename TIdSpec, typename TSeqString, typename TSeqSpec, typename TQualString,
          typename TQualSpec, typename TRecordReader, typename TTagList>
int read2(StringSet<TIdString, TIdSpec> & metas,
          StringSet<TSeqString, TSeqSpec> & seqs,
          StringSet<TQualString, TQualSpec> & quals,
          TRecordReader & reader,
          TagSelector<TTagList> & tagSelector)
{
    if (!guessStreamFormat(reader, tagSelector))
        return IOERR_READ_AUTO_FORMAT_UNKNOWN_FORMAT;

    return _read2(metas, seqs, quals, reader, tagSelector, LENGTH<TTagList>::VALUE - 1);
}

// For read2(metas, seqs, reader, TTagSelector/AutoSeqStreamFormat).

template <typename TIdString, typename TIdSpec, typename TSeqString, typename TSeqSpec,
          typename TRecordReader>
int _read2(StringSet<TIdString, TIdSpec> & /*metas*/,
           StringSet<TSeqString, TSeqSpec> & /*seqs*/,
           TRecordReader & /*reader*/,
           TagSelector<void> const & /*tagSelector*/,
           int i)
{
    SEQAN_ASSERT_EQ(i, -1);
    (void)i;  // Only used for assertion.
    return IOERR_READ_AUTO_FORMAT_UNKNOWN_FORMAT;
}

template <typename TIdString, typename TIdSpec, typename TSeqString, typename TSeqSpec,
          typename TRecordReader, typename TTagList>
int _read2(StringSet<TIdString, TIdSpec> & metas,
           StringSet<TSeqString, TSeqSpec> & seqs,
           TRecordReader & reader,
           TagSelector<TTagList> const & tagSelector,
           int i)
{
    if (i == -1)
        return IOERR_READ_AUTO_FORMAT_UNKNOWN_FORMAT;
    if (i != tagSelector.tagId)
        return _read2(metas, seqs, reader, static_cast<typename TagSelector<TTagList>::Base const &>(tagSelector), i - 1);

    typedef typename TTagList::Type TFormat;
    return read2(metas, seqs, reader, TFormat());
}

template <typename TIdString, typename TIdSpec, typename TSeqString, typename TSeqSpec,
          typename TRecordReader, typename TTagList>
int read2(StringSet<TIdString, TIdSpec> & metas,
          StringSet<TSeqString, TSeqSpec> & seqs,
          TRecordReader & reader,
          TagSelector<TTagList> & tagSelector)
{
    if (!guessStreamFormat(reader, tagSelector))
        return IOERR_READ_AUTO_FORMAT_UNKNOWN_FORMAT;

    return _read2(metas, seqs, reader, tagSelector, LENGTH<TTagList>::VALUE - 1);
}

// ----------------------------------------------------------------------------
// Function readRecord
// ----------------------------------------------------------------------------

// For readRecord(meta, seq, qual, reader, TTagSelector/AutoSeqStreamFormat).

template <typename TIdString, typename TSeqString, typename TQualString, typename TRecordReader>
int _readRecord(TIdString & /*meta*/, TSeqString & /*seq*/, TQualString & /*qual*/, TRecordReader & /*reader*/, TagSelector<void> const & /*tagSelector*/, int i)
{
    SEQAN_ASSERT_EQ(i, -1);
    (void)i;  // only used for assertion.
    return IOERR_READ_AUTO_FORMAT_UNKNOWN_FORMAT;
}

template <typename TIdString, typename TSeqString, typename TQualString, typename TRecordReader, typename TTagList>
int _readRecord(TIdString & meta, TSeqString & seq, TQualString & qual, TRecordReader & reader, TagSelector<TTagList> const & tagSelector, int i)
{
    if (i == -1)
        return IOERR_READ_AUTO_FORMAT_UNKNOWN_FORMAT;
    if (i != tagSelector.tagId)
        return _readRecord(meta, seq, qual, reader, static_cast<typename TagSelector<TTagList>::Base const &>(tagSelector), i - 1);

    typedef typename TTagList::Type TFormat;
    return readRecord(meta, seq, qual, reader, TFormat());
}

template <typename TIdString, typename TSeqString, typename TQualString, typename TRecordReader, typename TTagList>
int readRecord(TIdString & meta, TSeqString & seq, TQualString & qual, TRecordReader & reader, TagSelector<TTagList> & tagSelector)
{
    // TODO(holtgrew): The check should only happen once! Otherwise we waste time :(
    if (!guessStreamFormat(reader, tagSelector))
        return IOERR_READ_AUTO_FORMAT_UNKNOWN_FORMAT;

    return _readRecord(meta, seq, qual, reader, tagSelector, LENGTH<TTagList>::VALUE - 1);
}

// For readRecord(meta, seq, reader, TTagSelector/AutoSeqStreamFormat).

template <typename TIdString, typename TSeqString, typename TRecordReader>
int _readRecord(TIdString & /*meta*/, TSeqString & /*seq*/, TRecordReader & /*reader*/, TagSelector<void> const & /*tagSelector*/, int i)
{
    (void)i;  // Only used in assertion.
    SEQAN_ASSERT_EQ(i, -1);
    return IOERR_READ_AUTO_FORMAT_UNKNOWN_FORMAT;
}

template <typename TIdString, typename TSeqString, typename TRecordReader, typename TTagList>
int _readRecord(TIdString & meta, TSeqString & seq, TRecordReader & reader, TagSelector<TTagList> const & tagSelector, int i)
{
    if (i == -1)
        return IOERR_READ_AUTO_FORMAT_UNKNOWN_FORMAT;
    if (i != tagSelector.tagId)
        return _readRecord(meta, seq, reader, static_cast<typename TagSelector<TTagList>::Base const &>(tagSelector), i - 1);

    typedef typename TTagList::Type TFormat;
    return readRecord(meta, seq, reader, TFormat());
}

// (weese:)
// why is guessStreamFormat called insided readRecord?
// guessStreamFormat should be called outside and only once, directly after opening a stream

template <typename TIdString, typename TSeqString, typename TRecordReader, typename TTagList>
int readRecord(TIdString & meta, TSeqString & seq, TRecordReader & reader, TagSelector<TTagList> & tagSelector)
{
    // TODO(holtgrew): The check should only happen once! Otherwise we waste time :(
    if (!guessStreamFormat(reader, tagSelector))
        return IOERR_READ_AUTO_FORMAT_UNKNOWN_FORMAT;

    return _readRecord(meta, seq, reader, tagSelector, LENGTH<TTagList>::VALUE - 1);
}

}  // namespace seqan

#endif  // #ifndef SEQAN_STREAM_RECORD_READER_AUTO_FORMAT_H_
