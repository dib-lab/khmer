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
// Author: Jochen Singer <jochen.singer@fu-berlin.de>
// ==========================================================================

#ifndef CORE_INCLUDE_SEQAN_GFF_IO_GFF_IO_BASE_H_
#define CORE_INCLUDE_SEQAN_GFF_IO_GFF_IO_BASE_H_

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// ----------------------------------------------------------------------------
// Tag Gff
// ----------------------------------------------------------------------------

/**
.Tag.File Format.tag.Gff:
    Gff annotation file.
..include:seqan/gff_io.h
*/

// TODO(singer): const should be non const, but is const elsewhere
struct TagGff_;
typedef Tag<TagGff_> const Gff;

// ----------------------------------------------------------------------------
// Tag Gtf
// ----------------------------------------------------------------------------

/**
.Tag.File Format.tag.Gtf:
    Gtf annotation file.
..include:seqan/gff_io.h
*/

// TODO(singer): const should be non const, but is const elsewhere
struct TagGtf_;
typedef Tag<TagGtf_> const Gtf;

// ----------------------------------------------------------------------------
// Class GffRecord
// ----------------------------------------------------------------------------

/**
.Class.GffRecord
..cat:BAM I/O
..summary:Represent a record from a Gff file.
..include:seqan/gff_io.h

.Memvar.GffRecord#INVALID_POS
..class:Class.GffRecord
..summary:Static member with invalid/sentinel position value.
..type:nolink:$__uint32$

.Memvar.GffRecord#INVALID_SCORE
..class:Class.GffRecord
..summary:Static member with invalid score value.
..type:nolink:$float$

.Memvar.GffRecord#ref
..class:Class.GffRecord
..summary:The sequence id of the record.
..type:Shortcut.CharString

.Memvar.GffRecord#source
..class:Class.GffRecord
..summary:The source of the record.
..type:Shortcut.CharString

.Memvar.GffRecord#type
..class:Class.GffRecord
..summary:The type of the record.
..type:Shortcut.CharString

.Memvar.GffRecord#beginPos
..class:Class.GffRecord
..summary:The begin position of the record.
..type:nolink:$__uint32$

.Memvar.GffRecord#endPos
..class:Class.GffRecord
..summary:The end position of the record.
..type:nolink:$__uint32$

.Memvar.GffRecord#score
..class:Class.GffRecord
..summary:The score of the record.
..type:nolink:$float$

.Memvar.GffRecord#strand
..class:Class.GffRecord
..summary:The strand the record belongs to.
..type:nolink:$char$

.Memvar.GffRecord#phase
..class:Class.GffRecord
..summary:The phase of the record.
..remarks:For features of type "CDS", the phase indicates where the feature begins with reference to the reading frame. The phase is one of the integers 0, 1, or 2, indicating the number of bases that should be removed from the beginning of this feature to reach the first base of the next codon
..type:nolink:$char$

.Memvar.GffRecord#tagName
..class:Class.GffRecord
..summary:The names of the attributes of the record.
..type:Class.StringSet
..remarks:For each name there is a value associated in $Memvar.GffRecord#tagValue$

.Memvar.GffRecord#tagValue
..class:Class.GffRecord
..summary:The values of the attributes of the record.
..type:Class.StringSet
..remarks:For each value there is a name associated in $Memvar.GffRecord#tagName$
*/

struct GffRecord
{
    static __int32 const INVALID_POS = 2147483647;  // TODO(singer): Should be MaxValue<__int32>::VALUE, but that is not a constant expression :(
    static __int32 const INVALID_IDX = -1;

    // The member descriptions are taken from: http://gmod.org/wiki/GFF

    // TODO(singer): Maybe use a I/O context object and store ids as integers
    // The ID of the landmark used to establish the coordinate system for the current feature.
    String<char> ref;
    int rID;

    // The source is a free text qualifier intended to describe the algorithm or operating procedure that generated this feature.
    String<char> source;

    // The type of the feature
    String<char> type;

    // A list of feature attributes in the format tag=value.
    StringSet<String<char> > tagName;
    StringSet<String<char> > tagValue;

    // The start and end of the feature, in 1-based integer coordinates, relative to the landmark given in column 1
    __uint32 beginPos;
    __uint32 endPos;

    // The score of the feature
    float score;

    // The strand of the feature. + for positive strand (relative to the landmark), - for minus strand, and . for features that are not stranded.
    char strand;

    // For features of type "CDS", the phase indicates where the feature begins with reference to the reading frame.
    // The phase is one of the integers 0, 1, or 2, indicating the number of bases that should be removed from the beginning of this feature to reach the first base of the next codon.
    char phase;

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
        rID(INVALID_IDX), beginPos(-1), endPos(-1), score(INVALID_SCORE()),
        strand('.'), phase('.')
    {}
};

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function _parseReadGffKeyValue
// ----------------------------------------------------------------------------

template <typename TReader, typename TKeyString, typename TValueString>
inline int
_parseReadGffKeyValue(TValueString & outValue, TKeyString & key, TReader & reader)
{
    char c = value(reader);
    if (c == ' ' || c == '\t' || c == '\n' || c == '=')
        return 1;  // Key cannot be empty.

    for (; !atEnd(reader); goNext(reader))
    {
        c = value(reader);
        if (c == ' ' || c == '\t' || c == '\n' || c == '=' || c == ';')
            break;
        appendValue(key, c);
    }
    if (!atEnd(reader) && value(reader) == ';')
    {
        goNext(reader);
        return 0;
    }
    if (!atEnd(reader) && (value(reader) == '\r' || value(reader) == '\n'))
        return 0;

    if (skipWhitespaces(reader) != 0)
        return 1;

    if (value(reader) == '=')
    {
        goNext(reader);
        if (skipWhitespaces(reader) != 0)
            return 1;

        if (atEnd(reader))
            return 1;
    }

    if (value(reader) == '"')
    {
        // Handle the case of a string literal.

        goNext(reader);
        // Append all characters in the literal to outValue until the first i
        // line break or the closing '"'.
        for (; !atEnd(reader); goNext(reader))
        {
            if (value(reader) == '\n')
                return 1;

            if (value(reader) == '"')
            {
                goNext(reader);
                break;
            }
            appendValue(outValue, value(reader));
        }
        // Go over the trailing semicolon and any trailing space.
        while (!atEnd(reader) && (value(reader) == ';' || value(reader) == ' '))
            goNext(reader);
    }
    else
    {
        // Handle the non literal case.

        // Read until the first semicolon, return at whitespace.
        for (; !atEnd(reader); goNext(reader))
        {
            if (value(reader) == ';' || value(reader) == '\n' || value(reader) == '\r')
                break;
            appendValue(outValue, value(reader));
        }
        // Skip semicolon and spaces if any.
        while (!atEnd(reader) && (value(reader) == ';' || value(reader) == ' '))
            goNext(reader);
    }
    return 0;
}

// ----------------------------------------------------------------------------
// Function clear
// ----------------------------------------------------------------------------

/**
.Function.GffRecord#clear
..class:Class.GffRecord
..cat:Input/Output
..signature:clear(record)
..param.record:The @Class.GffRecord@ to reset.
...type:Class.GffRecord
..summary:Reset a @Class.GffRecord@ object.
..include:seqan/gff_io.h
*/

inline void clear(GffRecord & record)
{
    clear(record.ref);
    clear(record.source);
    clear(record.type);
    clear(record.tagName);
    clear(record.tagValue);
}

// ----------------------------------------------------------------------------
// Function readRecord
// ----------------------------------------------------------------------------

/**
.Function.GffRecord#readRecord
..class:Class.GffRecord
..cat:Input/Output
..summary:Read one gff record.
..signature:readRecord(record, reader)
..param.record:The gff record.
...type:Class.GffRecord
..param.reader: The record reader object.
...type:Class.RecordReader
..include:seqan/gff_io.h
*/

template <typename TStream, typename TRecordReaderSpec>
inline int
_readGffRecord(GffRecord & record, RecordReader<TStream, TRecordReaderSpec> & reader)
{
    clear(record);

    // read column 1: seqid
    // The letters until the first whitespace will be read.
    // Then, we skip until we hit the first tab character.
    if (readUntilTabOrLineBreak(record.ref, reader))
        return 1;
    record.rID = GffRecord::INVALID_IDX;

    if (!empty(record.ref) && record.ref[0] == '#')
    {
        if (skipLine(reader))
            return 1;

        return 1;
    }
    if (skipWhitespaces(reader))
        return 1;

    // read column 2: source
    if (readUntilTabOrLineBreak(record.source, reader))
        return 1;

    if (record.source == ".")
        clear(record.source);

    if (skipWhitespaces(reader))
        return 1;

    // read column 3: type
    if (readUntilTabOrLineBreak(record.type, reader))
        return 1;

    if (skipWhitespaces(reader))
        return 1;

    // read column 4: begin position
    String<char> temp;
    if (readDigits(temp, reader))
        return 1;

    if (length(temp) > 0u)
    {
        if (!lexicalCast2(record.beginPos, temp))
            return 1;

        --record.beginPos;  // Translate from 1-based to 0-based.
    }
    else
    {
        record.beginPos = GffRecord::INVALID_POS;
        if (skipUntilWhitespace(reader))
            return 1;
    }
    if (skipBlanks(reader))
        return 1;

    // read column 5: end position
    clear(temp);
    if (readDigits(temp, reader))
        return 1;

    if (length(temp) > 0u)
    {
        if (!lexicalCast2(record.endPos, temp))
            return 1;
    }
    else
    {
        record.endPos = GffRecord::INVALID_POS;
        if (skipUntilWhitespace(reader))
            return 1;
    }
    if (skipBlanks(reader))
        return 1;


    // read column 6: score
    clear(temp);
    readFloat(temp, reader);

    if (length(temp) > 0u)
    {
        if (temp != ".")
        {
            if (!lexicalCast2(record.score, temp))
                return 1;
        }
        else
        {
            record.score = GffRecord::INVALID_SCORE();
            if (skipUntilWhitespace(reader))
                return 1;
        }
    }
    else
    {
        return 1;
    }

    if (skipBlanks(reader))
        return 1;

    // read column 7: strand
    clear(temp);
    if (readUntilTabOrLineBreak(temp, reader))
        return 1;

    if (temp[0] != '-' && temp[0] != '+')
    {
        record.strand = '.';
    }
    else
    {
        record.strand = temp[0];
    }

    if (skipBlanks(reader))
        return 1;

    // read column 8: phase
    clear(temp);
    if (readUntilTabOrLineBreak(temp, reader))
        return 1;

    if (temp != "0" && temp != "1" && temp != "2")
    {
        record.phase = '.';
    }
    else
    {
        record.phase = temp[0];
    }

    if (skipBlanks(reader))
        return 1;

    // It's fine if there are no attributes and the line ends here.
    if (atEnd(reader))
        return 0;
    if (value(reader) == '\n' || value(reader) == '\r')
        return skipLine(reader);

    // read column 9: attributes
    while (!atEnd(reader))
    {

        String<char> _key;
        String<char> _value;
        // Read next key/value pair.
        if (_parseReadGffKeyValue(_value, _key, reader) != 0)
            return 1;

        appendValue(record.tagName, _key);
        appendValue(record.tagValue, _value);

        clear(_key);
        clear(_value);

        // At end of line:  Skip EOL and break.
        if (!atEnd(reader) && (value(reader) == '\r' || value(reader) == '\n'))
        {
            if (skipLine(reader) != 0)
                return 1;

            break;
        }
    }
    return 0;
}

template <typename TRecordReader>
inline int
readRecord(GffRecord & record, TRecordReader & reader, Gff /*tag*/)
{
    return _readGffRecord(record, reader);
}

template <typename TRecordReader>
inline int
readRecord(GffRecord & record, TRecordReader & reader, Gtf /*tag*/)
{
    return _readGffRecord(record, reader);
}

// TODO(singer): Needs proper documentation!!! Check the length of the stores!!!
template <typename TRecordReader, typename TContextSpec, typename TContextSpec2>
inline int
_readGffRecord(GffRecord & record, TRecordReader & reader, GffIOContext<TContextSpec, TContextSpec2> & context)
{
    // Read record with string ref from GFF file.
    int res = _readGffRecord(record, reader);
    if (res != 0)
        return res;

    // Translate ref to rID using the context.  If there is no such sequence name in the context yet then we add it.
    unsigned idx = 0;
    if (!getIdByName(nameStore(context), record.ref, idx, nameStoreCache(context)))
    {
        idx = length(nameStore(context));
        appendName(nameStore(context), record.ref, nameStoreCache(context));
    }
    record.rID = idx;

    return 0;
}

template <typename TRecordReader, typename TContextSpec, typename TContextSpec2>
inline int
readRecord(GffRecord & record, TRecordReader & reader, GffIOContext<TContextSpec, TContextSpec2> & context, Gff /*tag*/)
{
    return _readGffRecord(record, reader, context);
}

template <typename TRecordReader, typename TContextSpec, typename TContextSpec2>
inline int
readRecord(GffRecord & record, TRecordReader & reader, GffIOContext<TContextSpec, TContextSpec2> & context, Gtf /*tag*/)
{
    return _readGffRecord(record, reader, context);
}

// ----------------------------------------------------------------------------
// Function _writeSemicolonSensitive()
// ----------------------------------------------------------------------------

// This function checks if the string to be written contains a semicolon. If
// this is the case then quotes are written around the string.
// Returns false on success.

template <typename TTargetStream, typename TString>
inline bool
_writeSemicolon(TTargetStream & target, TString & temp)
{
    // TODO(jsinger): What about escaping quote chars '"'?
    if (streamWriteChar(target, '"') || streamWriteBlock(target, &temp[0], length(temp)) < length(temp) ||
        streamWriteChar(target, '"'))
        return true;

    return false;
}

template <typename TTargetStream, typename TString>
inline bool
_writeSemicolonSensitive(TTargetStream & target, TString & temp)
{
    // TODO(jsinger): What about escaping quote chars '"'?
    if (std::find(begin(temp), end(temp), ';') != end(temp))
    {
        return _writeSemicolon(target, temp);
    }
    else
    {
        if (streamWriteBlock(target, &temp[0], length(temp)) < length(temp))
            return true;
    }

    return false;
}

// ----------------------------------------------------------------------------
// Function writeRecord
// ----------------------------------------------------------------------------

/**
.Function.GffRecord#writeRecord
..class:Class.GffRecord
..cat:Input/Output
..summary:Writes one gff record to a stream.
..signature:writeRecord(TSreamm stream, GffRecord record)
..param.stream:The output stream.
...type:Concept.StreamConcept
..param.record:The gff record.
...type:Class.GffRecord
..include:seqan/gff_io.h
*/

template <typename TStream>
inline int
_writeAttributes(TStream & stream, GffRecord const & record, Gff /*tag*/)
{
    if (empty(record.tagName))
        return 0;

    unsigned i = 0;
    for (; i + 1 < length(record.tagName); ++i)
    {
        if (_writeSemicolonSensitive(stream, record.tagName[i]))
            return 1;

        if (length(record.tagValue[i]) > 0u)
        {
            if (streamWriteChar(stream, '='))
                return 1;

            if (_writeSemicolonSensitive(stream, record.tagValue[i]))
                return 1;
        }
        if (streamWriteChar(stream, ';'))
            return 1;

    }
    if (_writeSemicolonSensitive(stream, record.tagName[i]))
        return 1;

    if (length(record.tagValue[i]) > 0u)
    {
        if (streamWriteChar(stream, '='))
            return 1;

        if (_writeSemicolonSensitive(stream, record.tagValue[i]))
            return 1;
    }
    if (streamWriteChar(stream, '\n'))
        return 1;

    return 0;
}

template <typename TStream>
inline int
_writeAttributes(TStream & stream, GffRecord const & record, Gtf /*tag*/)
{
    unsigned i = 0;
    for (; i < length(record.tagName) - 1; ++i)
    {
        if (_writeSemicolonSensitive(stream, record.tagName[i]))
            return 1;

        if (streamWriteChar(stream, ' '))
            return 1;

        if (length(record.tagValue[i]) > 0u)
        {
            if (_writeSemicolon(stream, record.tagValue[i]))
                return 1;

            if (streamWriteBlock(stream, "; ", 2) != 2u)
                return 1;
        }
    }
    if (_writeSemicolonSensitive(stream, record.tagName[i]))
        return 1;

    if (streamWriteChar(stream, ' '))
        return 1;

    if (length(record.tagValue[i]) > 0u)
    {
        if (_writeSemicolon(stream, record.tagValue[i]))
            return 1;

        if (streamWriteBlock(stream, ";\n", 2) != 2u)
            return 1;
    }
    else
    {
        if (streamWriteChar(stream, '\n'))
            return 1;
    }

    return 0;
}

template <typename TStream, typename TSeqId, typename TTag>
inline int
_writeRecordImpl(TStream & stream, GffRecord const & record, TSeqId const & ref, TTag const tag)
{
    // ignore empty annotations, i.e. annotations that are 'guessed' by implicit information from their children (in GFF)
    if (empty(ref))
        return 0;

    // write column 1: seqid
    if (streamWriteBlock(stream, &ref[0], length(ref)) != length(ref))
        return 1;

    if (streamWriteChar(stream, '\t'))
        return 1;

    // write column 2: source
    if (empty(record.source))
    {
        if (streamWriteChar(stream, '.'))
            return 1;
    }
    else
    {
        if (streamWriteBlock(stream, &record.source[0], length(record.source)) != length(record.source))
            return 1;
    }

    if (streamWriteChar(stream, '\t'))
        return 1;

    // write column 3: type
    if (streamWriteBlock(stream, &record.type[0], length(record.type)) != length(record.type))
        return 1;

    if (streamWriteChar(stream, '\t'))
        return 1;

    // write column 4: begin position
    if (record.beginPos != (unsigned)-1)
    {
        if (streamPut(stream, record.beginPos + 1))
            return 1;
    }
    else
    {
        if (streamWriteChar(stream, '.'))
            return 1;
    }

    if (streamWriteChar(stream, '\t'))
        return 1;


    // write column 5: end position
    if (record.endPos != (unsigned)-1)
    {
        if (streamPut(stream, record.endPos))
            return 1;
    }
    else
    {
        if (streamWriteChar(stream, '.'))
            return 1;
    }

    if (streamWriteChar(stream, '\t'))
        return 1;

    // write column 6: score
    if (record.score != record.score)
    {
        if (streamWriteChar(stream, '.'))
            return 1;
    }
    else
    {
        if (streamPut(stream, record.score))
            return 1;
    }

    if (streamWriteChar(stream, '\t'))
        return 1;

    // write column 7: strand
    if (streamWriteChar(stream, record.strand))
        return 1;

    if (streamWriteChar(stream, '\t'))
        return 1;

    // write column 8: phase
    if (streamWriteChar(stream, record.phase))
        return 1;

    if (streamWriteChar(stream, '\t'))
        return 1;

    // write column 9: attributes
    // only until length - 1, because there is no semicolon at the end of the line

    _writeAttributes(stream, record, tag);

    return 0;
}

template <typename TStream, typename TTag>
inline int
writeRecord(TStream & stream, GffRecord const & record, TTag const tag)
{
    return _writeRecordImpl(stream, record, record.ref, tag);
}

template <typename TStream, typename TContextSpec, typename TContextSpec2, typename TTag>
inline int
writeRecord(TStream & stream, GffRecord const & record, GffIOContext<TContextSpec, TContextSpec2> & context, TTag const tag)
{
    if (record.rID != GffRecord::INVALID_IDX)
    {
        String<char> tempSeqId = nameStore(context)[record.rID];
        return _writeRecordImpl(stream, record, tempSeqId, tag);
    }
    return _writeRecordImpl(stream, record, record.ref, tag);
}

}  // namespace seqan

#endif  // CORE_INCLUDE_SEQAN_GFF_IO_GFF_IO_BASE_H_
