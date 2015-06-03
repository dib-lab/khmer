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
// Author: Jochen Singer <jochen.singer@fu-berlin.de>
// ==========================================================================

#ifndef INCLUDE_SEQAN_GFF_IO_GFF_IO_BASE_H_
#define INCLUDE_SEQAN_GFF_IO_GFF_IO_BASE_H_

namespace seqan {

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// ----------------------------------------------------------------------------
// Tag Gff
// ----------------------------------------------------------------------------

/*!
 * @tag FileFormats#Gff
 * @brief Tag for selecting the GFF format.
 *
 * Both the GFF and the GTF file format are represented by @link GffRecord @endlink in SeqAn.
 * Tags and functions in this group can be used for I/O of both formats to and from @link GffRecord @endlink objects.
 *
 * @signature typedef Tag<TagGff_> Gff;
 */
struct TagGff_;
typedef Tag<TagGff_> Gff;

// ----------------------------------------------------------------------------
// Tag Gtf
// ----------------------------------------------------------------------------

/*!
 * @tag FileFormats#Gtf
 * @brief Tag for selecting the GTF format.
 *
 * @signature typedef Tag<TagGtf_> Gtf;
 */
struct TagGtf_;
typedef Tag<TagGtf_> Gtf;

// ----------------------------------------------------------------------------
// Class MagicHeader
// ----------------------------------------------------------------------------

template <typename T>
struct MagicHeader<Gtf, T> :
    public MagicHeader<Nothing, T> {};

template <typename T>
struct MagicHeader<Gff, T> :
    public MagicHeader<Nothing, T> {};

// ----------------------------------------------------------------------------
// Class FileExtensions
// ----------------------------------------------------------------------------

template <typename T>
struct FileExtensions<Gff, T>
{
    static char const * VALUE[2];    // default is one extension
};

template <typename T>
char const * FileExtensions<Gff, T>::VALUE[2] =
{
    ".gff",     // default output extension
    ".gff3"
};

template <typename T>
struct FileExtensions<Gtf, T>
{
    static char const * VALUE[1];    // default is one extension
};

template <typename T>
char const * FileExtensions<Gtf, T>::VALUE[1] =
{
    ".gtf"     // default output extension
};

// ----------------------------------------------------------------------------
// Class GffRecord
// ----------------------------------------------------------------------------

/*!
 * @class GffRecord
 * @implements FormattedFileRecordConcept
 * @headerfile <seqan/gff_io.h>
 * @brief Represent a record from a GFF or GTF file.
 *
 * @signature class GffRecord;
 */
struct GffRecord
{
    /*!
     * @var __int32 GffRecord::INVALID_IDX;
     * @brief Static member with invalid/sentinel rID value.
     */
    static __int32 const INVALID_POS = 2147483647;  // TODO(singer): Should be MaxValue<__int32>::VALUE, but that is not a constant expression :(

    /*!
     * @var CharString GffRecord::ref;
     * @brief The sequence name of the record.
     *
     * The ID of the landmark used to establish the coordinate system for the current feature, most often the
     * contig/chromosome name.
     */
    CharString ref;

    /*!
     * @var CharString GffRecord::source;
     * @brief The source of the record.
     *
     * The source is a free text qualifier intended to describe the algorithm or operating procedure that generated this
     * feature.
     */
    CharString source;

    /*!
     * @var CharString GffRecord::type;
     * @brief The type of the record.
     */
    CharString type;

    /*!
     * @var TCharStringSet GffRecord::tagNames;
     * @brief The names of the attributes of the record, StringSet of CharString.
     *
     * For each value there is a name associated in @link GffRecord::tagNames tagNames @endlink.
     */
    StringSet<CharString> tagNames;

    /*!
     * @var TCharStringSet GffRecord::tagValues;
     * @brief The values of the attributes of the record, StringSet of CharString.
     *
     * @section Remarks
     *
     * For each name there is a value associated in GffRecord::tagValues.
     */
    StringSet<CharString> tagValues;

    /*!
     * @var __int32 GffRecord::beginPos;
     * @brief The begin position of the record.
     */
    __uint32 beginPos;

    /*!
     * @var __int32 GffRecord::endPos;
     * @brief The end position of the record.
     *
     * GFF and GTF use 1-based positions in text, but they are stored as 0-based coordinates.
     */
    __uint32 endPos;

    /*!
     * @var float GffRecord::score;
     * @brief The score of the record.
     */
    float score;

    /*!
     * @var char GffRecord::strand;
     * @brief The strand the record belongs to.
     *
     * The strand of the feature. + for positive strand (relative to the landmark), - for minus strand, and . for
     * features that are not stranded.
     */
    char strand;

    /*!
     * @var char GffRecord::phase;
     * @brief The phase of the record.
     *
     * For features of type "CDS", the phase indicates where the feature begins with reference to the reading frame.
     * The phase is one of the integers 0, 1, or 2, indicating the number of bases that should be removed from the
     * beginning of this feature to reach the first base of the next codon.
     */
    char phase;

    // TODO(holtgrew): C++11 will have a nan() function, use this instead then.
    /*!
     * @fn GffRecord::INVALID_SCORE
     * @signature static float INVALID_SCORE()
     * @brief Returns invalid score (NaN float value).
     *
     * The term <tt>x != x</tt> (for <tt>float x</tt> is only true if <tt>x</tt> is a NaN.
     */
    static float INVALID_SCORE()
    {
        union
        {
            __uint32 u;
            float f;
        } tmp;
        tmp.u = 0x7F800001;
        return tmp.f;
    }

    GffRecord() :
        beginPos(-1), endPos(-1), score(INVALID_SCORE()),
        strand('.'), phase('.')
    {}
};

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function _parseReadGffKeyValue
// ----------------------------------------------------------------------------

template <typename TForwardIter, typename TKeyString, typename TValueString>
inline void
_parseReadGffKeyValue(TValueString & outValue, TKeyString & key, TForwardIter & iter)
{
    //TODO(singer): AssertList functor would be need
    char c = value(iter);
    if (IsWhitespace()(c) || c == '=')
        SEQAN_THROW(ParseError("The key field of an attribute is empty!"));

    for (; !atEnd(iter); goNext(iter))
    {
        c = value(iter);
        if (IsNewline()(c) || c == ' ' || c == '=' || c == ';')
            break;
        appendValue(key, c);
    }
    if (!atEnd(iter) && value(iter) == ';')
    {
        skipOne(iter);
        return;
    }

    if (IsNewline()(value(iter)))
        return;

    skipUntil(iter, NotFunctor<IsWhitespace>());

    if (value(iter) == '=')
    {
        skipOne(iter);
        skipUntil(iter, NotFunctor<IsWhitespace>());
    }

    if (value(iter) == '"')
    {
        // Handle the case of a string literal.
        skipOne(iter);
        skipUntil(iter, NotFunctor<IsWhitespace>());
        readUntil(outValue, iter, OrFunctor<EqualsChar<'"'>, AssertFunctor<NotFunctor<IsNewline>, ParseError, Gff> >());
        skipOne(iter);

        // Go over the trailing semicolon and any trailing space.
        skipUntil(iter, NotFunctor<OrFunctor<EqualsChar<';'>, EqualsChar<' '> > >());
    }
    else
    {
        // Read until the first semicolon, return at whitespace.
        readUntil(outValue, iter, OrFunctor<EqualsChar<';'>, IsNewline>());

        // Skip semicolon and spaces if any.
        skipUntil(iter, NotFunctor<OrFunctor<EqualsChar<';'>, EqualsChar<' '> > >());
    }
    return;
}

// ----------------------------------------------------------------------------
// Function clear
// ----------------------------------------------------------------------------

/*!
 * @fn GffRecord#clear
 * @brief Reset a @link GffRecord @endlink object.
 *
 * @signature void clear(record);
 *
 * @param[in,out] record The GffRecord to reset.
 */
inline void clear(GffRecord & record)
{
    record.beginPos = -1;
    record.endPos = -1;
    record.score = record.INVALID_SCORE();
    record.strand = '.';
    record.phase = '.';

    clear(record.ref);
    clear(record.source);
    clear(record.type);
    clear(record.tagNames);
    clear(record.tagValues);
}

// ----------------------------------------------------------------------------
// Function readRecord
// ----------------------------------------------------------------------------

// NOTE(esiragusa): dox disabled.
/*
 * @fn GffFileIO#readRecord
 * @brief Read one GFF/GTF record from a SinglePassRecordReader.
 *
 * @signature void readRecord(record, context, iter);
 *
 * @param[out]    record  The GffRecord to write the results to.
 * @param[in,out] context A CharString to use for buffers.
 * @param[in,out] iter    A @link ForwardIteratorConcept forward iterator @endlink to use for reading.
 *
 * @throws IOError if something went wrong.
 */
template <typename TFwdIterator>
void readRecord(GffRecord & record, CharString & buffer, TFwdIterator & iter)
{
    IsNewline isNewline;

    clear(record);

    skipUntil(iter, NotFunctor<OrFunctor<EqualsChar<'#'>, IsWhitespace> >());  //skip commments and empty lines

    // read column 1: seqid
    readUntil(record.ref, iter, OrFunctor<IsTab, AssertFunctor<NotFunctor<IsNewline>, ParseError, Gff> >());
    skipOne(iter);

    // read column 2: source
    readUntil(record.source, iter, OrFunctor<IsTab, AssertFunctor<NotFunctor<IsNewline>, ParseError, Gff> >());

    if (record.source == ".")
        clear(record.source);

    skipOne(iter);

    // read column 3: type
    readUntil(record.type, iter, OrFunctor<IsTab, AssertFunctor<NotFunctor<IsNewline>, ParseError, Gff> >());
    skipOne(iter);

    // read column 4: begin position
    clear(buffer);
    readUntil(buffer, iter, OrFunctor<IsTab, AssertFunctor<NotFunctor<IsNewline>, ParseError, Gff> >());
    record.beginPos = lexicalCast<__uint32>(buffer);
    --record.beginPos;  // Translate from 1-based to 0-based.
    skipOne(iter);

    // read column 5: end position
    clear(buffer);
    readUntil(buffer, iter, OrFunctor<IsTab, AssertFunctor<NotFunctor<IsNewline>, ParseError, Gff> >());
    record.endPos = lexicalCast<__uint32>(buffer);
    skipOne(iter);

    //check if end < begin
    if (record.endPos < record.beginPos)
        SEQAN_THROW(ParseError("Begin position of GFF/GTF record is larger than end position!"));

    // read column 6: score
    clear(buffer);
    readUntil(buffer, iter, OrFunctor<IsTab, AssertFunctor<NotFunctor<IsNewline>, ParseError, Gff> >());
    if (buffer != ".")
        record.score = lexicalCast<float>(buffer);
    skipOne(iter, IsTab());

    // read column 7: strand
    readOne(record.strand, iter, OrFunctor<OrFunctor<EqualsChar<'-'>, EqualsChar<'+'> >, EqualsChar<'.'> >());
    skipOne(iter, IsTab());

    // read column 8: phase
    readOne(record.phase, iter, OrFunctor<EqualsChar<'.'>, IsInRange<'0', '2'> >());

    // It's fine if there are no attributes and the line ends here.
    if (atEnd(iter) || isNewline(value(iter)))
    {
        skipLine(iter);
        return;
    }
    skipOne(iter, IsTab());

    // read column 9: attributes
    while (!atEnd(iter))
    {

        CharString _key;
        CharString _value;
        // Read next key/value pair.
        _parseReadGffKeyValue(_value, _key, iter);

        appendValue(record.tagNames, _key);
        appendValue(record.tagValues, _value);

        clear(_key);
        clear(_value);

        // At end of line:  Skip EOL and break.
        if (!atEnd(iter) && isNewline(value(iter)))
        {
            skipOne(iter);
            break;
        }
    }
    return;
}

// ----------------------------------------------------------------------------
// Function _writeSemicolonSensitive()
// ----------------------------------------------------------------------------

// This function checks if the string to be written contains a semicolon. If
// this is the case then quotes are written around the string.
// Returns false on success.

template <typename TTargetStream, typename TString>
inline void
_writeInQuotes(TTargetStream & target, TString & temp)
{
    // TODO(jsinger): What about escaping quote chars '"'?
    writeValue(target, '"');
    write(target, temp);
    writeValue(target, '"');
}

template <typename TTarget, typename TString, typename TMustBeQuotedFunctor>
inline void
_writePossiblyInQuotes(TTarget& target, TString & source, TMustBeQuotedFunctor const &func)
{
    // TODO(jsinger): What about escaping quote chars '"'?
    typedef typename Iterator<TString>::Type TIter;
    TIter itEnd = end(source, Standard());
    for (TIter it = begin(source, Standard()); it != itEnd; ++it)
    {
        // we have a problem if the string contains a '"' or a line break
        if (value(it) =='\n' || value(it) == '"')
            SEQAN_THROW(ParseError("Attribute contains illegal character!"));

        if (func(*it))
        {
            _writeInQuotes(target, source);
            return;
        }
    }
    write(target, source);
}

// ----------------------------------------------------------------------------
// Function writeRecord()
// ----------------------------------------------------------------------------

// NOTE(esiragusa): dox disabled.
/*
 * @fn GffFileIO#writeRecord
 * @brief Writes a @link GffRecord @endlink to a stream as GFF or GTF.
 *
 * @signature void writeRecord(stream, record, tag);
 *
 * @param[in,out] stream  The @link OutputIteratorConcept output iterator @endlink to write to.
 * @param[in]     record  The @link GffRecord @endlink to write out.
 * @param[in]     tag     A tag to select the file format, either @link GffFileIO#Gff @endlink or @link GffFileIO#Gtf
 *                        @endlink.
 *
 * @throws IOError if something went wrong.
 */

template <typename TFormatTag>
struct GffRecordKeyMustBeQuoted_;

template <typename TFormatTag>
struct GffRecordValueMustBeQuoted_;

// GFF quotation rules

template <>
struct GffRecordKeyMustBeQuoted_<Gff>
{
    bool operator() (char c) const
    {
        return c == ';' || c == '=';
    }
};

template <>
struct GffRecordValueMustBeQuoted_<Gff> :
    GffRecordKeyMustBeQuoted_<Gff> {};

// GTF quotation rules

template <>
struct GffRecordKeyMustBeQuoted_<Gtf>
{
    bool operator() (char c) const
    {
        return c == ';' || c == ' ';
    }
};

template <>
struct GffRecordValueMustBeQuoted_<Gtf>
{
    bool operator() (char c) const
    {
//        return c == ';' || c == ' ' || !isdigit(c);
        return !isdigit(c);     // is equivalent to the above, quote everything except integral values
    }
};

template <typename TTarget>
inline void
_writeAdditionalSeperator(TTarget const & /*target*/, Gff)
{
    return;
}

template <typename TTarget>
inline void
_writeAdditionalSeperator(TTarget & target, Gtf)
{
    writeValue(target, ' ');
    return;
}


template <typename TTarget, typename TTag>
inline void
_writeAttributes(TTarget & target, GffRecord const & record, TTag const & tag)
{
    const char separatorBetweenTagAndValue = (IsSameType<TTag, Gff>::VALUE)? '=' : ' ';
    for (unsigned i = 0; i < length(record.tagNames); ++i)
    {
        if (i != 0)
        {
            writeValue(target, ';');

            // In GTF files a space follows the semicolon
            _writeAdditionalSeperator(target, tag);
       }

        _writePossiblyInQuotes(target, record.tagNames[i], GffRecordKeyMustBeQuoted_<TTag>());

        if (!empty(record.tagValues[i]))
        {
            writeValue(target, separatorBetweenTagAndValue);
            _writePossiblyInQuotes(target, record.tagValues[i], GffRecordValueMustBeQuoted_<TTag>());
        }
    }

    // In GTF files each (especially the last) attribute must end with a semi-colon
    if (IsSameType<TTag, Gtf>::VALUE && !empty(record.tagNames))
        writeValue(target, ';');

    return;
}

template <typename TTarget, typename TFormat>
inline void
writeRecord(TTarget & target, GffRecord const & record, Tag<TFormat> const & tag)
{
    // ignore empty annotations, i.e. annotations that are 'guessed' by implicit information from their children (in GFF)
    if (empty(record.ref))
        return;

    // write column 1: seqid
    //typename Iterator<TSeqId const, Rooted>::Type itRef = begin(record.ref);
    write(target, record.ref);
    writeValue(target, '\t');

    // write column 2: source
    if (empty(record.source))
        writeValue(target, '.');
    else
        write(target, record.source);
    writeValue(target, '\t');

    // write column 3: type
    write(target, record.type);
    writeValue(target, '\t');

    // write column 4: begin position
    if (record.beginPos != (unsigned)-1)
        appendNumber(target, record.beginPos + 1);
    else
        SEQAN_THROW(ParseError("No start position!"));
    writeValue(target, '\t');

    // write column 5: end position
    if (record.endPos != (unsigned)-1 && record.beginPos <= record.endPos)
        appendNumber(target, record.endPos);
    else
        SEQAN_THROW(ParseError("No end position!"));
    writeValue(target, '\t');

    // write column 6: score
    if (record.score != record.score)
        writeValue(target, '.');
    else
        appendNumber(target, record.score);
    writeValue(target, '\t');

    // write column 7: strand
    writeValue(target, record.strand);
    writeValue(target, '\t');

    // write column 8: phase
    writeValue(target, record.phase);
    writeValue(target, '\t');

    // write column 9: attributes
    // only until length - 1, because there is no semicolon at the end of the line

    _writeAttributes(target, record, tag);

    writeValue(target, '\n');
    return;
}

}  // namespace seqan

#endif  // INCLUDE_SEQAN_GFF_IO_GFF_IO_BASE_H_

