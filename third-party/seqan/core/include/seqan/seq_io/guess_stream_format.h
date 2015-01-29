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
// Author: David Weese <david.weese@fu-berlin.de>
// ==========================================================================
// File Format detection with Stream Class.
// ==========================================================================

#ifndef SEQAN_SEQ_IO_GUESSSEQ_IOFORMAT_H_
#define SEQAN_SEQ_IO_GUESSSEQ_IOFORMAT_H_

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// NOTE: Must be kept in sync with Function.getAutoSeqStreamFormatName().

/**
.Shortcut.SeqStreamFormats
..cat:Input/Output
..summary:A Tag list of the currently implemented Sequence-Formats (in RecordReader/Stream-IO)
..signature:SeqStreamFormats
..shortcutfor:Tag.TagList
...signature:TagList<Fastq, TagList<Fasta > >
..include:seqan/stream.h
*/

// TODO(weese): File formats should be independent of the I/O interface?
// Why adding SeqStreamFormats and AutoSeqStreamFormat instead of using
// the existing AutoSeqFormat and SeqFormats?
typedef TagList<Fastq, TagList<Fasta> > SeqStreamFormats;

/**
.Shortcut.AutoSeqStreamFormat
..cat:Input/Output
..summary:A TagSelector for @Shortcut.SeqStreamFormats@, the list of the currently implemented Sequence-Formats (in RecordReader/Stream-IO)
..signature:AutoSeqStreamFormat
..shortcutfor:Class.TagSelector
..shortcutfor:Tag.TagList
..remarks:This shortcut is a typedef of $TagSelector<TagList<Fastq, TagList<Fasta > > >$.
..remarks:Variables of this type an be passed to @Function.guessStreamFormat@ and will offer the index of the detected FileFormat in its member $tagId$.
The values of its member $tagId$ can be as follows:
..remarks.tableheader:value|meaning
..remarks.table:0|Unknown file format.
..remarks.table:1|FASTA
..remarks.table:2|FASTQ
..see:Shortcut.SeqStreamFormats
..see:Function.guessStreamFormat
..include:seqan/stream.h
*/

typedef TagSelector<SeqStreamFormats> AutoSeqStreamFormat;

/**
.Class.LimitRecordReaderInScope
..cat:Input/Output
..summary:manipulates a @Class.RecordReader@ -Object so that it operates only on one buffer
..signature:LimitRecordReaderInScope<TStream, TSpec>
..param.TStream:The @Concept.StreamConcept@ of the @Class.RecordReader@.
...type:Concept.StreamConcept
..param.TSpec:The specialization of the @Class.RecordReader@.
...type:Class.RecordReader
..see:Class.RecordReader
..see:Function.guessStreamFormat
..include:seqan/stream.h
..remarks:This class is intended for situations, where you do not wish the RecordReader to rebuffer and where you wish to return to the original reading position after reading, e.g. when detecting the file format of the stream.
..remarks:It is used by passing the RecordReader-object on construction (this already does the necessary changes in the RecordReader). Upon deconstruction of this object, the RecordReader is reset to its original state, including all iterators.
..remarks:This works on all RecordReader-objects, independent of the underlying stream-object. It also works, if the underlying stream does not support seeking.
..include:seqan/stream.h
*/

// TODO(holtgrew): Rename to LimitRecordReaderRIIA?

template <typename TStream, typename TPass>
class LimitRecordReaderInScope
{
public:
    RecordReader<TStream, TPass> & _recordreader;
    typename RecordReader<TStream, TPass>::TIter _currentBeforeSuspend;
    // the following is only needed for MMapStrings
    typename RecordReader<TStream, TPass>::TIter _endBeforeSuspend;

    LimitRecordReaderInScope(RecordReader<TStream, TPass> & reader)
            : _recordreader(reader),
              _currentBeforeSuspend(reader._current),
              _endBeforeSuspend(reader._end)
    {
        _suspendRefill(TPass());
    }

    ~LimitRecordReaderInScope()
    {
        _resumeRefillAndReset(TPass());
    }
private:
    inline void
    _suspendRefill(SinglePass<void> const & /* tag */)
    {
        _recordreader._stayInOneBuffer = true;
    }

    inline void
    _suspendRefill(SinglePass<StringReader> const & /* tag */)
    {
        if (_recordreader._end - _recordreader._current > BUFSIZ)
            _recordreader._end = _recordreader._current + BUFSIZ;
    }

    template <typename TSpec>
    inline void
    _suspendRefill(DoublePass<TSpec> const & /* tag */)
    {
        _suspendRefill(SinglePass<TSpec>());
    }


    inline void
    _resumeRefillAndReset(SinglePass<void> const & /* tag */)
    {
        _recordreader._stayInOneBuffer = false;
        if (_currentBeforeSuspend == 0) // there had been no Buffer at start
            _recordreader._current = begin(_recordreader._buffer, Standard());
        else
            _recordreader._current = _currentBeforeSuspend;
    }

    inline void
    _resumeRefillAndReset(DoublePass<void> const & /* tag */)
    {
        _recordreader._stayInOneBuffer = false;
        if (_currentBeforeSuspend == 0) // there had been no Buffer at start
            _recordreader._current = begin(*_recordreader._currentBuffer,
                                           Standard());
        else
            _recordreader._current = _currentBeforeSuspend;
    }

    inline void
    _resumeRefillAndReset(SinglePass<StringReader> const & /* tag */)
    {
        _recordreader._end = _endBeforeSuspend;
        _recordreader._current = _currentBeforeSuspend;
    }

    template <typename TSpec>
    inline void
    _resumeRefillAndReset(DoublePass<TSpec> const & /* tag */)
    {
        _resumeRefillAndReset(SinglePass<TSpec>());
    }
};

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function autoSeqStreamFormatName()
// ----------------------------------------------------------------------------

/**
.Function.getAutoSeqStreamFormatName
..summary:Returns
..signature:getAutoSeqStreamFormatName(autoSeqStreamFormat)
..param.autoSeqStreamFormat:The @Shortcut.AutoSeqStreamFormat@ object to get the format name from.
...type:Shortcut.AutoSeqStreamFormat
..returns:$char const *$ with the format name.
..include:seqan/stream.h
..see:Function.guessStreamFormat
..see:Shortcut.AutoSeqStreamFormat
..example.code:
// (Create RecordReader object recordReader with sequence file here).
AutoSeqStreamFormat formatTag;
if (checkFormat(recordReader, formatTag))
    std::cout << "File format is " << getAutoSeqStreamFormatName(formatTag) << '\n';
*/

// NOTE: Must be kept in sync with the typedef SeqStreamFormats.

inline char const *
getAutoSeqStreamFormatName(AutoSeqStreamFormat const & format)
{
    switch (format.tagId)
    {
        case Find<SeqStreamFormats, Fasta>::VALUE: return "FASTA";
        case Find<SeqStreamFormats, Fastq>::VALUE: return "FASTQ";
    }

    return "INVALID FORMAT";
}

// ----------------------------------------------------------------------------
// Function guessStreamFormat()
// ----------------------------------------------------------------------------

/**
.Function.guessStreamFormat
..class:Class.RecordReader
..cat:Input/Output
..summary:check whether the data provided by reader is (one of) the specified format(s).
..signature:guessStreamFormat(TRecordReader & reader, TTag const &)
..param.reader:The @Class.RecordReader@ to read from
..param.TTag:The tag to check against.
..signature:guessStreamFormat(TRecordReader & reader, TagSelector<TTagList> & formats)
..param.formats:A @Class.TagSelector@ object that contains the list of tags to check and provides a tagId member with index of the detected tag.
..returns: $true$ if (one of) the specified Tag(s) tested positive and $False$ otherwise
...type:nolink:$bool$
..remarks:With the help of @Class.LimitRecordReaderInScope@ these functions do not (permanently) alter the position in the stream.
..remarks:The tagId-member of the TagSelector holds the index in inside-to-outside order and begins counting at one. E.g. The Index of FASTQ in TagList<Fastq, TagList<Fasta > > would be 2.
..include:seqan/stream.h
..example.text:
The following example guesses the sequence file format of the already open fstream $in$.
After the call to $guessStreamFormat()$, the $tagSelector.tagId$ contains the 1-based index of the matching tag.
Here, we use the @Shortcut.AutoSeqStreamFormat@ tag selector.
..example.code:
RecordReader<std::fstream, SinglePass<> > reader(in);
AutoSeqStreamFormat tagSelector;
bool b = guessStreamFormat(reader, tagSelector);
// b is true if any format was detected successfully.
if (tagSelector.tagId == 1)
    std::cerr << "Detected FASTA." << std::endl;
else if (tagSelector.tagId == 2)
    std::cerr << "Detected FASTQ." << std::endl;
else
    std::cerr << "Unknown file format!" << std::endl;
..example.text:
Alternatively, we can define your own tag selector.
Note that we reverse the order of FASTA and FASTQ in respect to @Shortcut.AutoSeqStreamFormat@
..example.code:
typedef TagSelector<TagList<Fasta, TagList<Fastq> > > MyTagSelector;

RecordReader<std::fstream, SinglePass<> > reader(in);
AutoSeqStreamFormat tagSelector;
MyTagSelector tagSelector;
bool b = guessStreamFormat(reader, tagSelector);
// b is true if any format was detected successfully.
if (tagSelector.tagId == 1)
    std::cerr << "Detected FASTQ." << std::endl;
else if (tagSelector.tagId == 2)
    std::cerr << "Detected FASTA." << std::endl;
else
    std::cerr << "Unknown file format!" << std::endl;

*/

template < typename TRecordReader >
inline bool
guessStreamFormat(TRecordReader &, TagSelector<> &)
{
    // we get here if the file format could not be determined
    return false;
}

template <typename TRecordReader, typename TTagList >
inline bool
guessStreamFormat(TRecordReader & reader, TagSelector<TTagList> & format)
{
    // (weese:) compare with guessFormat in file_format_mmap.h
    typedef typename TTagList::Type TFormatTag;

    if (value(format) == -1 || value(format) == Find<TTagList, TFormatTag>::VALUE)
    {
        // if tagId is set to -1 (auto-detect) or the current format (TFormatTag) then test for TFormatTag format
        if (guessStreamFormat(reader, TFormatTag()))
        {
            value(format) = Find<TTagList, TFormatTag>::VALUE;
            return true;
        }
    }
    return guessStreamFormat(reader, static_cast<typename TagSelector<TTagList>::Base &>(format));
}

}  // namespace seqan

#endif  // #ifndef SEQAN_SEQ_IO_ADAPT_FSEQ_IO_H_
