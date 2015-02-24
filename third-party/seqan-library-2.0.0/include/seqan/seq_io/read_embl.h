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
// Read support for the EMBL file format.
// ==========================================================================

#ifndef INCLUDE_SEQAN_SEQ_IO_READ_EMBL_H_
#define INCLUDE_SEQAN_SEQ_IO_READ_EMBL_H_

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

struct Embl_;
typedef Tag<Embl_> Embl;

struct EmblHeader_;
typedef Tag<EmblHeader_> EmblHeader;

struct EmblSequence_;
typedef Tag<EmblSequence_> EmblSequence;

// ============================================================================
// Metafunctions
// ============================================================================

// --------------------------------------------------------------------------
// Metafunction MagicHeader
// --------------------------------------------------------------------------

template <typename T>
struct MagicHeader<Embl, T>
{
    static char const VALUE[3];
};
template <typename T>
char const MagicHeader<Embl, T>::VALUE[3] = { 'I','D',' ' };  // typical Embl header

// --------------------------------------------------------------------------
// Metafunction FileExtensions
// --------------------------------------------------------------------------

template <typename T>
struct FileExtensions<Embl, T>
{
    static char const * VALUE[1];
};
template <typename T>
char const * FileExtensions<Embl, T>::VALUE[1] =
{
    ".embl",     // default output extension
};

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function splitEmblHeader()
// ----------------------------------------------------------------------------

/*!
 * @fn readEmblHeader
 * @headerfile <seqan/seq_io.h>
 * @brief Split an EMBL header line.
 *
 * @signature void readEmblHeader(key, value, iter);
 *
 * @param[out] key   A @link ContainerConcept @endlink object to write the key to.
 * @param[out] value A @link ContainerConcept @endlink object to write the value to.
 * @param[in]  iter  Input iterator.
 */

template <typename TKey, typename TValue, typename TFwdIterator>
inline void
readRecord(TKey & key, TValue & val, TFwdIterator & iter, EmblHeader)
{
    clear(key);
    clear(val);

    skipUntil(iter, NotFunctor<IsWhitespace>());
    readUntil(key, iter, IsWhitespace(), AssertFunctor<NotFunctor<CountDownFunctor<True, 2> >, ParseError>());
    if (!atEnd(iter) && IsBlank()(value(iter)))
    {
        skipUntil(iter, NotFunctor<IsBlank>());
        readUntil(val, iter, IsNewline());
    }
    skipLine(iter);
}

// ----------------------------------------------------------------------------
// Function nextIs()
// ----------------------------------------------------------------------------

// nextIs() for EMBL header records.

template <typename TFwdIterator>
inline bool
nextIs(TFwdIterator & iter, EmblHeader)
{
    return !atEnd(iter) && !IsWhitespace()(value(iter));
}

// nextIs() for EMBL sequence records.

template <typename TFwdIterator>
inline bool
nextIs(TFwdIterator & iter, EmblSequence)
{
    return !atEnd(iter) && IsWhitespace()(value(iter));
}


// ----------------------------------------------------------------------------
// Function readRecord()
// ----------------------------------------------------------------------------

// readRecord() for EMBL sequences records.

// Read all sequence, eat/ignore '//' line.

template <typename TSeqString, typename TFwdIterator>
inline void
readRecord(TSeqString & seq, TFwdIterator & iter, EmblSequence)
{
    typedef typename Value<TSeqString>::Type TSeqAlphabet;
    typedef OrFunctor<OrFunctor<IsBlank, IsDigit>,
                      AssertFunctor<IsInAlphabet<TSeqAlphabet>, ParseError, Embl> > TSeqAsserter;

    IsNewline isNewline;
    TSeqAsserter asserter;

    if (!atEnd(iter) && value(iter) == 'I')
        skipLine(iter);

    clear(seq);
    for (; !atEnd(iter); skipLine(iter))
    {
        if (value(iter) == '/')
        {
            skipLine(iter);
            break;
        }

        readUntil(seq, iter, isNewline, asserter);
    }
}

// readRecord() for EMBL id/seq pairs.

template <typename TIdString, typename TSeqString, typename TFwdIterator>
inline void
readRecord(TIdString & meta, TSeqString & seq, TFwdIterator & iter, Embl)
{
    IsBlank isBlank;
    IsWhitespace isWhite;
    String<char, Array<2> > key;

    // extract meta from "ID" line
    clear(meta);
    for (; !atEnd(iter) && !isBlank(value(iter)); skipLine(iter))
    {
        AssertFunctor<NotFunctor<CountDownFunctor<True, 2> >, ParseError, Embl> asserter;

        clear(key);
        readUntil(key, iter, isWhite, asserter);
        if (!atEnd(iter) && isBlank(value(iter)) && key == "ID")
        {
            skipUntil(iter, NotFunctor<IsBlank>());
            readUntil(meta, iter, IsNewline());
        }
    }

    readRecord(seq, iter, EmblSequence());
}

template <typename TIdString, typename TSeqString, typename TQualString, typename TFwdIterator>
inline void
readRecord(TIdString & meta, TSeqString & seq, TQualString & qual, TFwdIterator & iter, Embl)
{
    clear(qual);
    readRecord(meta, seq, iter, Embl());
}

}  // namespace seqan

#endif  // #ifndef INCLUDE_SEQAN_SEQ_IO_READ_EMBL_H_
