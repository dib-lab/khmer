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
// I/O functionality for score matrices.
// ==========================================================================

// TODO(holtgrew): We need to adapt this to the new I/O scheme.

#ifndef CORE_INCLUDE_SEQAN_SCORE_SCORE_MATRIX_IO_H_
#define CORE_INCLUDE_SEQAN_SCORE_SCORE_MATRIX_IO_H_

#include <seqan/stream.h>

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// ----------------------------------------------------------------------------
// Tag ScoreMatrixFile
// ----------------------------------------------------------------------------

// TODO(holtgrew): What format is this actually?

/**
.Tag.File Format.tag.ScoreMatrixFile:Score matrix file.
..include:seqan/score.h
*/

struct TagScoreMatrixFile_;
typedef Tag<TagScoreMatrixFile_> ScoreMatrixFile;

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Helper Function _sscanfValue()                                 [DEPRECATED!]
// ----------------------------------------------------------------------------

/*
.Function._sscanfValue:
..cat:Input/Output
..summary:Use sscanf to parse a value from a $char *$ buffer.
..signature:_sscanfValue(buffer, value)
..param.buffer:Buffer to parse into.
...type:const char *
..param.value:Variable to parse the value from $buffer$ to.
...type:unsigned int
..include:seqan/score.h
 */

inline void
SEQAN_DEPRECATED_PRE("Use new RecordReader-based parsing functionality instead")
_sscanfValue(const char * buf, unsigned int & val)
SEQAN_DEPRECATED_POST("Use new RecordReader-based parsing functionality instead");

inline void
_sscanfValue(const char * buf, unsigned int & val)
{
    SEQAN_CHECKPOINT;
    std::sscanf(buf, "%u", & val);
}

/*
.Function._sscanfValue.param.value.type:int
..include:seqan/score.h
 */

inline void
SEQAN_DEPRECATED_PRE("Use new RecordReader-based parsing functionality instead")
_sscanfValue(const char * buf, int & val)
SEQAN_DEPRECATED_POST("Use new RecordReader-based parsing functionality instead");

inline void
_sscanfValue(const char * buf, int & val)
{
    SEQAN_CHECKPOINT;
    std::sscanf(buf, "%i", & val);
}

/*
.Function._sscanfValue.param.value.type:float
..include:seqan/score.h
 */

inline void
SEQAN_DEPRECATED_PRE("Use new RecordReader-based parsing functionality instead")
_sscanfValue(const char * buf, float & val)
SEQAN_DEPRECATED_POST("Use new RecordReader-based parsing functionality instead");

inline void
_sscanfValue(const char * buf, float & val)
{
    SEQAN_CHECKPOINT;
    std::sscanf(buf, "%f", & val);
}

/*
.Function._sscanfValue.param.value.type:double
..include:seqan/score.h
 */

inline void
SEQAN_DEPRECATED_PRE("Use new RecordReader-based parsing functionality instead")
_sscanfValue(const char * buf, double & val)
SEQAN_DEPRECATED_POST("Use new RecordReader-based parsing functionality instead");

inline void
_sscanfValue(const char * buf, double & val)
{
    SEQAN_CHECKPOINT;
    std::sscanf(buf, "%lf", & val);
}

// ----------------------------------------------------------------------------
// Function readMeta()                                        [ScoreMatrixFile]
// ----------------------------------------------------------------------------

// TODO(holtgrew): Fix parameter order, adjust to new style I/O system.

template <typename TStream, typename TMeta>
int
readMeta(RecordReader<TStream, SinglePass<> > & reader,
         TMeta & meta,
         ScoreMatrixFile const & /*tag*/)
{
    clear(meta);
    while (!atEnd(reader) && value(reader) == '#')
    {
        goNext(reader);
        int res = readLine(meta, reader);
        if (res != 0 && res != EOF_BEFORE_SUCCESS)
            return 1;
    }

    return 0;
}

// ----------------------------------------------------------------------------
// Function read()                                            [ScoreMatrixFile]
// ----------------------------------------------------------------------------

template <typename TStream, typename TValue, typename TSequenceValue, typename TSpec>
int
read(RecordReader<TStream, SinglePass<> > & reader,
     Score<TValue, ScoreMatrix<TSequenceValue, TSpec> > & sc,
     ScoreMatrixFile const & /*tag*/)
{
    typedef Score<TValue, ScoreMatrix<TSequenceValue, TSpec> > TScore;

    // Clear the matrix.
    std::fill(sc.data_tab, sc.data_tab + TScore::TAB_SIZE, TValue());

    // Skip header.
    while (!atEnd(reader) && value(reader) == '#')
        if (skipLine(reader) != 0)
            return 1;  // Also at EOF!

    // -----------------------------------------------------------------------
    // Read alphabet line.
    // -----------------------------------------------------------------------
    // The first non-header line is the alphabet line.  Read in the characters and build a map from column to character.
    // Such a line could read " A R N D C Q E G H I L K M F P S T W Y V B Z X *".
    String<TSequenceValue> mapping;
    // Skip blanks (TAB and SPACE).
    while (!atEnd(reader) && isblank(value(reader)))
        goNext(reader);
    while (!atEnd(reader) && value(reader) != '\n' && value(reader) != '\r')
    {
        appendValue(mapping, value(reader));  // Conversion char => TSequenceValue
        goNext(reader);
        // Skip blanks (TAB and SPACE).
        while (!atEnd(reader) && isblank(value(reader)))
            goNext(reader);
    }
    if (skipLine(reader) != 0)  // Go to next line.
        return 1;

    // -----------------------------------------------------------------------
    // Read the matrix.
    // -----------------------------------------------------------------------
    // The matrix is read line-wise, we read the score for a (source, dest) replacement.
    for (unsigned row = 0; !atEnd(reader); ++row)
    {
        while (!atEnd(reader) && isblank(value(reader)))  // Skip blanks (TAB and SPACE).
            goNext(reader);
        if (atEnd(reader))
            return 1;  // No label.

        // Read label.
        TSequenceValue source = value(reader);
        unsigned const OFFSET = ordValue(source) * TScore::VALUE_SIZE;
        goNext(reader);
        while (!atEnd(reader) && isblank(value(reader)))  // Skip blanks (TAB and SPACE).
            goNext(reader);
        if (atEnd(reader))
            return 1;  // No entries.

        // Read numbers.
        CharString buffer;
        for (unsigned col = 0; !atEnd(reader) && value(reader) != '\r' && value(reader) != '\n'; ++col)
        {
            if (atEnd(reader))
                return 1;  // No number.

            clear(buffer);
            int res = readUntilWhitespace(buffer, reader);
            if (res != 0 && res != EOF_BEFORE_SUCCESS)
                return 1;  // Error reading.
            TValue val = 0;
            if (!lexicalCast2(val, buffer))
                return 1;  // Invalid number.

            TSequenceValue destChar = mapping[col];
            sc.data_tab[OFFSET + ordValue(destChar)] = val;

            while (!atEnd(reader) && isblank(value(reader)))  // Skip blanks (TAB and SPACE).
                goNext(reader);
        }

        int res = skipLine(reader);  // Go to next line.
        if (res != 0)
            return 1;
    }

    return 0;
}


template <typename TStream, typename TValue, typename TSequenceValue, typename TSpec>
inline void
read(RecordReader<TStream, SinglePass<> > & fl,
     Score<TValue, ScoreMatrix<TSequenceValue, TSpec> > & sc)
{
    read(fl, sc, ScoreMatrixFile());
}

// ----------------------------------------------------------------------------
// Function loadScoreMatrix()
// ----------------------------------------------------------------------------

/**
.Function.loadScoreMatrix
..cat:Input/Output
..summary:Load a score matrix from a file.
..signature:loadScoreMatrix(score, filename)
..param.score:Score matrix object to load into.
...type:Spec.Score Matrix
..param.filename:Path to the file to load.
..include:seqan/score.h
*/
template <typename TValue, typename TSequenceValue, typename TSpec>
inline void
loadScoreMatrix(Score<TValue, ScoreMatrix<TSequenceValue, TSpec> > & sc,
                CharString const & filename)
{
    std::fstream fin(toCString(filename), std::ios::binary | std::ios::in);
    RecordReader<std::fstream, SinglePass<> > reader(fin);
    read(reader, sc);
}

// TODO(holtgrew): Change parameter order?

/**
.Function.loadScoreMatrix
..signature:loadScoreMatrix(score, filename, meta)
..param.meta:Meta information read from the file.
..status:The order of the parameters filename and meta will switch.
..include:seqan/score.h
*/
template <typename TValue, typename TSequenceValue, typename TSpec, typename TMeta>
inline void
loadScoreMatrix(Score<TValue, ScoreMatrix<TSequenceValue, TSpec> > & sc, CharString const & filename, TMeta & meta)
{
    std::fstream fin(toCString(filename), std::ios::binary | std::ios::in);
    RecordReader<std::fstream, SinglePass<> > reader(fin);
    readMeta(reader, meta, ScoreMatrixFile());
    read(reader, sc, ScoreMatrixFile());
}

// ----------------------------------------------------------------------------
// Helper Function _sprintfValue()
// ----------------------------------------------------------------------------

// TODO(holtgrew): Can we solve this more eleganly using std::stringstream's?

/*
.Function._sprintfValue
..cat:Input/Output
..summary:Use sprintf to print a value into a $const char *$ buffer.
..signature:_sprintf(buffer, value)
..param.buffer:Buffer to write to.
...type:char *
...remark:Must be of sufficient size.
..param.value:Variable for which to write the string value to $buffer$.
...type:unsigned int
...type:int
...type:float
...type:double
..includes:seqan/score.h
..include:seqan/score.h
 */

inline void
_sprintfValue(char * buf, unsigned val)
{
    std::sprintf(buf, "%u", val);
}

/*
.Function._sprintfValue.param.value.type:int
..include:seqan/score.h
 */
inline void
_sprintfValue(char * buf, int val)
{
    std::sprintf(buf, "%d", val);
}

/*
.Function._sprintfValue.param.value.type:float
..include:seqan/score.h
 */
inline void
_sprintfValue(char * buf, float val)
{
    double d = val;
    std::sprintf(buf, "%G", d);
}

/*
.Function._sprintfValue.param.value.type:float
..include:seqan/score.h
 */
inline void
_sprintfValue(char * buf, double val)
{
    std::sprintf(buf, "%G", val);
}

// ----------------------------------------------------------------------------
// Function write()                                         [Spec.Score Matrix]
// ----------------------------------------------------------------------------

/**
.Function.Score Matrix#write
..summary:Write a score matrix to a stream.
..cat:Input/Output
..class:Spec.Score Matrix
..signature:write(stream, scoreMatrix[, meta])
..param.stream:A @Concept.StreamConcept@ to write to.
...type:Concept.StreamConcept
..param.scoreMatrix:The matrix to write out.
...type:Spec.Score Matrix
..param.meta:Optional description.
...type:Shortcut.CharString
..include:seqan/score.h
 */

template <typename TStream, typename TValue, typename TSequenceValue, typename TSpec>
inline int
write(TStream & stream,
      Score<TValue, ScoreMatrix<TSequenceValue, TSpec> > const & sc,
      CharString const & meta)
{
    typedef typename ValueSize<TSequenceValue>::Type TSeqValueSize;
    TValue const * tab = sc.data_tab;

    TSeqValueSize const VALUE_SIZE = ValueSize<TSequenceValue>::VALUE;
    __uint32 const TAB_SIZE = VALUE_SIZE * VALUE_SIZE;

    // -----------------------------------------------------------------------
    // Write meta data.
    // -----------------------------------------------------------------------
    if (!empty(meta)) {
        bool lineBegin = true;
        for (unsigned int i = 0; i < length(meta); ++i) {
            if (lineBegin)
            {
                // Escape each line with a starting '#'.
                streamPut(stream, '#');
                lineBegin = false;
            }
            if (meta[i] == '\r')
                continue;
            if (meta[i] == '\n')
                lineBegin = true;
            streamPut(stream, meta[i]);
        }
        if (!lineBegin)
            streamPut(stream, '\n');
    }

    // -----------------------------------------------------------------------
    // Determine maximal column width.
    // -----------------------------------------------------------------------
    unsigned colWidth = 1;
    char buf[100];  // 100 is enough for printing int, unsigned, dobule
    for (unsigned i = 0; i < TAB_SIZE; ++i)
    {
        _sprintfValue(buf, tab[i]);
        unsigned cellWidth = std::strlen(buf);
        if (cellWidth > colWidth)
            colWidth = cellWidth;  // Compute maximum.
    }
    colWidth += 1;  // Increase colWidth for an additional blank.

    // -----------------------------------------------------------------------
    // Write alphabet line.
    // -----------------------------------------------------------------------
    streamPut(stream, ' ');  // A blank for alphabet column.
    for (unsigned j = 0; j < VALUE_SIZE; ++j)
    {
        char val = static_cast<TSequenceValue>(j);  // Conversion integral => TSequenceValue => char.
        // Leading blanks for column j.
        for (unsigned k = 1; k < colWidth; ++k)
            streamPut(stream, ' ');
        streamPut(stream, val);
    }
    streamPut(stream, '\n');

    // -----------------------------------------------------------------------
    // Write rest of matrix.
    // -----------------------------------------------------------------------
    for (unsigned i = 0; i < VALUE_SIZE; ++i)
    {
        // Write alphabet column cell.
        char val = static_cast<TSequenceValue>(i);  // Conversion integral => TSequenceValue => char.
        streamPut(stream, val);

        // Write rest of line i.
        unsigned offset = i * VALUE_SIZE;
        for (unsigned j = 0; j < VALUE_SIZE; ++j)
        {
            _sprintfValue(buf, tab[offset + j]);
            unsigned len = strlen(buf);

            // Leading blanks.
            for (unsigned k = 0; k < colWidth - len; ++k)
                streamPut(stream, ' ');

            // Write cell.
            streamPut(stream, buf);
        }
        streamPut(stream, '\n');
    }

    return 0;
}

template <typename TStream, typename TValue, typename TSequenceValue, typename TSpec>
inline int
write(TStream & stream,
      Score<TValue, ScoreMatrix<TSequenceValue, TSpec> > const & sc)
{
    return write(stream, sc, "");
}

}  // namespace seqan

#endif  // #ifndef CORE_INCLUDE_SEQAN_SCORE_SCORE_MATRIX_IO_H_
