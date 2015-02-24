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
// I/O functionality for score matrices.
// ==========================================================================

// TODO(holtgrew): We need to adapt this to the new I/O scheme.

#ifndef INCLUDE_SEQAN_SCORE_SCORE_MATRIX_IO_H_
#define INCLUDE_SEQAN_SCORE_SCORE_MATRIX_IO_H_

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

struct TagScoreMatrixFile_;
typedef Tag<TagScoreMatrixFile_> ScoreMatrixFile;

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function readMeta()                                        [ScoreMatrixFile]
// ----------------------------------------------------------------------------

template <typename TMeta, typename TSourceIter>
inline void
readMeta(TMeta & meta,
         TSourceIter & reader,
         ScoreMatrixFile const & /*tag*/)
{
    clear(meta);
    while (!atEnd(reader) && value(reader) == '#')
    {
        skipOne(reader);
        readLine(meta, reader);
    }
}

// ----------------------------------------------------------------------------
// Function read()                                            [ScoreMatrixFile]
// ----------------------------------------------------------------------------

template <typename TValue, typename TSequenceValue, typename TSpec, typename TSourceIter>
inline void
read(Score<TValue, ScoreMatrix<TSequenceValue, TSpec> > & sc,
     TSourceIter & reader,
     ScoreMatrixFile const & /*tag*/)
{
    typedef Score<TValue, ScoreMatrix<TSequenceValue, TSpec> > TScore;

    // Clear the matrix.
    arrayFill(sc.data_tab, sc.data_tab + TScore::TAB_SIZE, TValue());

    // Skip header.
    while (!atEnd(reader) && value(reader) == '#')
        skipLine(reader);

    // -----------------------------------------------------------------------
    // Read alphabet line.
    // -----------------------------------------------------------------------
    // The first non-header line is the alphabet line.  Read in the characters and build a map from column to character.
    // Such a line could read " A R N D C Q E G H I L K M F P S T W Y V B Z X *".
    String<TSequenceValue> mapping;

    skipUntil(reader, NotFunctor<IsWhitespace>());
    readUntil(mapping, reader, IsNewline(), IsWhitespace());    // Conversion char => TSequenceValue, Skip blanks
    skipLine(reader);

    // -----------------------------------------------------------------------
    // Read the matrix.
    // -----------------------------------------------------------------------
    // The matrix is read line-wise, we read the score for a (source, dest) replacement.

    CharString buffer;
    for (unsigned row = 0; !atEnd(reader); ++row)
    {
        skipUntil(reader, NotFunctor<IsBlank>());

        // Read label.
        TSequenceValue source;
        readOne(source, reader);
        unsigned const offset = ordValue(source) * TScore::VALUE_SIZE;

        skipUntil(reader, NotFunctor<IsBlank>());

        // Read numbers.
        for (unsigned col = 0; !atEnd(reader) && !IsNewline()(value(reader)); ++col)
        {
            readUntil(buffer, reader, IsWhitespace());
            sc.data_tab[offset + ordValue(mapping[col])] = lexicalCast<TValue>(buffer);
            clear(buffer);
            skipUntil(reader, NotFunctor<IsBlank>());
        }

        skipLine(reader);  // Go to next line.
    }
}


template <typename TValue, typename TSequenceValue, typename TSpec, typename TSourceIter>
inline void
read(Score<TValue, ScoreMatrix<TSequenceValue, TSpec> > & sc,
     TSourceIter & iter)
{
    read(sc, iter, ScoreMatrixFile());
}

// ----------------------------------------------------------------------------
// Function loadScoreMatrix()
// ----------------------------------------------------------------------------

/*!
 * @fn Score#loadScoreMatrix
 * @brief Load a score matrix from a file.
 *
 * @signature void loadScoreMatrix(score, fileName[, meta]);
 *
 * @param[in,out] score    The MatrixScore to load.
 * @param[in]     fileName The path to the file to load, CharString.
 * @param[out]    meta     Meta information is stored here, CharString.
 */

template <typename TValue, typename TSequenceValue, typename TSpec>
inline bool
loadScoreMatrix(Score<TValue, ScoreMatrix<TSequenceValue, TSpec> > & sc,
                const char * filename)
{
    VirtualStream<char, Input> fin;
    if (!open(fin, toCString(filename)))
        return false;

    typename DirectionIterator<VirtualStream<char, Input>, Input>::Type reader = directionIterator(fin, Input());
    read(sc, reader);
    return true;
}

template <typename TMeta, typename TValue, typename TSequenceValue, typename TSpec>
inline bool
loadScoreMatrix(TMeta & meta,
                Score<TValue, ScoreMatrix<TSequenceValue, TSpec> > & sc,
                const char * filename)
{
    VirtualStream<char, Input> fin;
    if (!open(fin, toCString(filename)))
        return false;

    typename DirectionIterator<VirtualStream<char, Input>, Input>::Type reader = directionIterator(fin, Input());
    readMeta(meta, reader, ScoreMatrixFile());
    read(sc, reader, ScoreMatrixFile());
    return true;
}

// ----------------------------------------------------------------------------
// Function write()                                         [Spec.Score Matrix]
// ----------------------------------------------------------------------------

/*!
 * @fn MatrixScore#write
 * @brief Write a score matrix to a stream.
 *
 * @signature int write(stream, scoreMatrix[, meta]);
 *
 * @param[in,out] stream      Stream to write to.
 * @param[out]    scoreMatrix MatrixScore object.
 * @param[out]    meta        Meta information is stored here.
 *
 * @return int Status code, 0 on success, a different value on errors.
 */

template <typename TTarget, typename TValue, typename TSequenceValue, typename TSpec>
inline void
write(TTarget & target,
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
    if (!empty(meta))
    {
        bool lineBegin = true;
        for (unsigned int i = 0; i < length(meta); ++i) {
            if (lineBegin)
            {
                // Escape each line with a starting '#'.
                writeValue(target, '#');
                lineBegin = false;
            }
            if (meta[i] == '\r')
                continue;
            if (meta[i] == '\n')
                lineBegin = true;
            writeValue(target, meta[i]);
        }
        if (!lineBegin)
            writeValue(target, '\n');
    }

    // -----------------------------------------------------------------------
    // Determine maximal column width.
    // -----------------------------------------------------------------------
    size_t colWidth = 1;
    String<char, Array<100> > buf;  // 100 is enough for printing int, unsigned, dobule
    for (unsigned i = 0; i < TAB_SIZE; ++i)
    {
        clear(buf);
        appendNumber(buf, tab[i]);
        colWidth = std::max(colWidth, (size_t)length(buf));
    }
    colWidth++;     // Increase colWidth for an additional blank.

    // -----------------------------------------------------------------------
    // Write alphabet line.
    // -----------------------------------------------------------------------
    writeValue(target, ' ');  // A blank for alphabet column.
    for (unsigned j = 0; j < VALUE_SIZE; ++j)
    {
        char val = static_cast<TSequenceValue>(j);  // Conversion integral => TSequenceValue => char.
        // Leading blanks for column j.
        for (unsigned k = 1; k < colWidth; ++k)
            writeValue(target, ' ');
        writeValue(target, val);
    }
    writeValue(target, '\n');

    // -----------------------------------------------------------------------
    // Write rest of matrix.
    // -----------------------------------------------------------------------
    for (unsigned i = 0; i < VALUE_SIZE; ++i)
    {
        // Write alphabet column cell.
        char val = static_cast<TSequenceValue>(i);  // Conversion integral => TSequenceValue => char.
        writeValue(target, val);

        // Write rest of line i.
        unsigned offset = i * VALUE_SIZE;
        for (unsigned j = 0; j < VALUE_SIZE; ++j)
        {
            clear(buf);
            appendNumber(buf, tab[offset + j]);

            // Leading blanks.
            for (unsigned k = length(buf); k < colWidth; ++k)
                writeValue(target, ' ');

            // Write cell.
            write(target, buf);
        }
        writeValue(target, '\n');
    }
}

template <typename TTarget, typename TValue, typename TSequenceValue, typename TSpec>
inline void
write(TTarget & target,
      Score<TValue, ScoreMatrix<TSequenceValue, TSpec> > const & sc)
{
    write(target, sc, "");
}

template <typename TStream, typename TValue, typename TSequenceValue, typename TSpec>
inline TStream &
operator<<(TStream & target,
           Score<TValue, ScoreMatrix<TSequenceValue, TSpec> > const & sc)
{
    typename DirectionIterator<TStream, Output>::Type it = directionIterator(target, Output());
    write(it, sc);
    return target;
}

}  // namespace seqan

#endif  // #ifndef INCLUDE_SEQAN_SCORE_SCORE_MATRIX_IO_H_
