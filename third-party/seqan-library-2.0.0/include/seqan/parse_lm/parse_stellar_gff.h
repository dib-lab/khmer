// ==========================================================================
//                                  parse_lm
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

// SEQAN_NO_GENERATED_FORWARDS

#ifndef INCLUDE_SEQAN_PARSE_LM_PARSE_STELLAR_GFF_H_
#define INCLUDE_SEQAN_PARSE_LM_PARSE_STELLAR_GFF_H_

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

struct StellarGff_;
typedef Tag<StellarGff_> StellarGff;

// ============================================================================
// Metafunctions
// ============================================================================

// ----------------------------------------------------------------------------
// Metafunction Value
// ----------------------------------------------------------------------------

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function readRecord()
// ----------------------------------------------------------------------------

template <typename TLocalMatchStore, typename TForwardIter>
inline void
readRecord(TLocalMatchStore & store,
           TForwardIter & iter,
           StellarGff const & /*tag*/)
{
    typedef typename TLocalMatchStore::TPosition TPosition;

    // Read line.
    CharString buffer;
    char subjectStrand = 'X';

    CharString subjectName;
    CharString queryName;
    TPosition subjectBeginPos = 0;
    TPosition subjectEndPos = 0;
    TPosition queryBeginPos = 0;
    TPosition queryEndPos = 0;

    // Field: SUBJECT
    readUntil(subjectName, iter, IsTab());

    // Skip TAB.
    skipOne(iter, IsTab());

    // Field: SOURCE
    skipUntil(iter, IsTab());

    // Skip TAB.
    skipOne(iter, IsTab());

    // Field: TYPE
    skipUntil(iter, IsTab());

    // Skip TAB.
    skipOne(iter, IsTab());

    // Field: START
    clear(buffer);
    readUntil(buffer, iter, NotFunctor<IsDigit>());
    subjectBeginPos = lexicalCast<TPosition>(buffer) - 1;

    // Skip TAB.
    skipOne(iter, IsTab());

    // Field: END
    clear(buffer);
    readUntil(buffer, iter, NotFunctor<IsDigit>());
    subjectEndPos = lexicalCast<TPosition>(buffer);

    // Skip TAB.
    skipOne(iter, IsTab());

    // Field: SCORE
    skipUntil(iter, IsTab());

    // Skip TAB.
    skipOne(iter, IsTab());

    // Field: STRAND
    clear(buffer);
    readOne(subjectStrand, iter, OrFunctor<EqualsChar<'+'>, EqualsChar<'-'> >());

    // Skip TAB.
    skipOne(iter, IsTab());

    // Field: FRAME
    skipOne(iter);

    // Skip TAB.
    skipOne(iter, IsTab());

    // The GROUP field contains the information about the query.
    // query sequence
    readUntil(queryName, iter, EqualsChar<';'>());

    // Skip semicolon.
    skipOne(iter, EqualsChar<';'>());

    // "seq2Range="
    clear(buffer);
    readUntil(buffer, iter, EqualsChar<'='>());
    if (buffer != "seq2Range")
    {
        throw ParseError("Format error; Expected 'seq2Range'");
        return;  // FORMAT ERROR, should probably be a constant
    }
    // Skip '='
    skipOne(iter, EqualsChar<'='>());

    // query begin pos
    clear(buffer);
    readUntil(buffer, iter, NotFunctor<IsDigit>());
    queryBeginPos = lexicalCast<TPosition>(buffer) - 1;

    // skip comma
    skipOne(iter, EqualsChar<','>());

    // query end pos
    clear(buffer);
    readUntil(buffer, iter, NotFunctor<IsDigit>());
    queryEndPos = lexicalCast<TPosition>(buffer);

    // Skip semicolon.
    skipOne(iter, EqualsChar<';'>());

    // Skip "eValue=*;"
    clear(buffer);
    readUntil(buffer, iter, EqualsChar<'='>());
    if (buffer != "eValue")
    {
        throw ParseError("Format error; Expected 'eValue'");
        return;  // FORMAT ERROR, should probably be a constant
    }

    // Skip '='
    skipOne(iter, EqualsChar<'='>());

    // Skip until semicolon.
    skipUntil(iter, EqualsChar<';'>());

    // And skip the semicolon.
    skipOne(iter, EqualsChar<';'>());

    // Skip "cigar".
    clear(buffer);
    readUntil(buffer, iter, EqualsChar<'='>());
    if (buffer != "cigar")
    {
        throw ParseError("Format error; Expected 'cigar'");
        return;  // FORMAT ERROR, should probably be a constant
    }

    // Skip '='
    skipOne(iter, EqualsChar<'='>());

    // Read cigar string.
    clear(buffer);
    readUntil(buffer, iter, EqualsChar<';'>());

    // Skip semicolon.
    skipOne(iter, EqualsChar<';'>());

    // ignore rest of the field, skip to next line
    skipLine(iter);

    // Finally, append the local match.
    if (subjectStrand == '-')
        std::swap(subjectBeginPos, subjectEndPos);
    appendLocalMatch(store, subjectName, subjectBeginPos, subjectEndPos, queryName, queryBeginPos, queryEndPos, buffer);
}

}  // namespace seqan

#endif  // INCLUDE_SEQAN_PARSE_LM_PARSE_STELLAR_GFF_H_
