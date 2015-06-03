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
// Parsing for the "general" format of lastz.
// ==========================================================================

// SEQAN_NO_GENERATED_FORWARDS

#ifndef INCLUDE_SEQAN_PARSE_LM_PARSE_LASTZ_GENERAL_H_
#define INCLUDE_SEQAN_PARSE_LM_PARSE_LASTZ_GENERAL_H_

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

struct LastzGeneral_;
typedef Tag<LastzGeneral_> LastzGeneral;

// ============================================================================
// Metafunctions
// ============================================================================

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
           LastzGeneral const & /*tag*/)
{
    typedef typename TLocalMatchStore::TPosition TPosition;

    // Skip any comments.
    while (true)
    {
        if (atEnd(iter))
        {
            throw UnexpectedEnd();
            return;
        }
        if (value(iter) != '#')
            break;
        skipLine(iter);
    }

    SEQAN_ASSERT_NEQ(value(iter), '#');

    // Read line.
    CharString buffer;

    CharString subjectName;
    CharString queryName;
    char subjectStrand = 'X';
    char queryStrand = 'X';
    TPosition subjectBeginPos = 0;
    TPosition subjectEndPos = 0;
    TPosition queryBeginPos = 0;
    TPosition queryEndPos = 0;

    // Field: score
    readUntil(buffer, iter, NotFunctor<IsDigit>());

    // Skip TAB.
    skipOne(iter, IsTab());

    // Field: name1
    readUntil(subjectName, iter, IsTab());

    // Skip TAB.
    skipOne(iter, IsTab());

    // Field: strand1
    readOne(subjectStrand, iter, OrFunctor<EqualsChar<'+'>, EqualsChar<'-'> >());

    // Skip TAB.
    skipOne(iter, IsTab());

    // Field: size1
    clear(buffer);
    readUntil(buffer, iter, NotFunctor<IsDigit>());

    // Skip TAB.
    skipOne(iter, IsTab());

    // Field: zstart1
    clear(buffer);
    readUntil(buffer, iter, NotFunctor<IsDigit>());
    subjectBeginPos = lexicalCast<TPosition>(buffer);

    // Skip TAB.
    skipOne(iter, IsTab());

    // Field: end1
    clear(buffer);
    readUntil(buffer, iter, NotFunctor<IsDigit>());
    subjectEndPos = lexicalCast<TPosition>(buffer);

    // Skip TAB.
    skipOne(iter, IsTab());

    // Field: name2
    readUntil(queryName, iter, IsTab());

    // Skip TAB.
    skipOne(iter, IsTab());

    // Field: strand2
    readOne(queryStrand, iter, OrFunctor<EqualsChar<'+'>, EqualsChar<'-'> >());

    // Skip TAB.
    skipOne(iter, IsTab());

    // Field: size2
    clear(buffer);
    readUntil(buffer, iter, NotFunctor<IsDigit>());

    // Skip TAB.
    skipOne(iter, IsTab());

    // Field: zstart2
    clear(buffer);
    readUntil(buffer, iter, NotFunctor<IsDigit>());
    queryBeginPos = lexicalCast<TPosition>(buffer);

    // Skip TAB.
    skipOne(iter, IsTab());

    // Field: end2
    clear(buffer);
    readUntil(buffer, iter, NotFunctor<IsDigit>());
    queryEndPos = lexicalCast<TPosition>(buffer);

    // Skip TAB.
    skipOne(iter, IsTab());

    // Field: identity
    clear(buffer);
    readUntil(buffer, iter, IsWhitespace());

    // Skip TAB.
    skipOne(iter, IsTab());

    // Field: idPct
    clear(buffer);
    readUntil(buffer, iter, IsWhitespace());

    // Skip TAB.
    skipOne(iter, IsTab());

    // Field: coverage
    clear(buffer);
    readUntil(buffer, iter, IsWhitespace());

    // Skip TAB.
    skipOne(iter, IsTab());

    // Field: covPct
    clear(buffer);
    readUntil(buffer, iter, IsWhitespace());

    // Skip to next line, ignore EOF here.
    skipLine(iter);

    // Finally, append the local match.
    if (subjectStrand == '-')
        std::swap(subjectBeginPos, subjectEndPos);
    if (queryStrand == '-')
        std::swap(queryBeginPos, queryEndPos);
    appendLocalMatch(store, subjectName, subjectBeginPos, subjectEndPos, queryName, queryBeginPos, queryEndPos);
}

}  // namespace seqan

#endif  // INCLUDE_SEQAN_PARSE_LM_PARSE_LASTZ_GENERAL_H_
