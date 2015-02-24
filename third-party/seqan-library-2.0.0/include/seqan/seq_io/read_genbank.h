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
// Read support for the GenBank file format.
//
// Note: We do not want to require a lookahead because then we might need too
// many buffers.  Thus, we heuristically only look at the first character of
// the field name to differentiate between headers and the sequence start
// label "ORIGIN".
//
// The sequence identifier for record-reading is the VERSION field.
// ==========================================================================


#ifndef INCLUDE_SEQAN_SEQ_IO_READ_GENBANK_H_
#define INCLUDE_SEQAN_SEQ_IO_READ_GENBANK_H_

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

struct GenBank_;
typedef Tag<GenBank_> GenBank;

struct GenBankHeader_;
typedef Tag<GenBankHeader_> GenBankHeader;

struct GenBankSequence_;
typedef Tag<GenBankSequence_> GenBankSequence;


// ============================================================================
// Metafunctions
// ============================================================================

// --------------------------------------------------------------------------
// Metafunction MagicHeader
// --------------------------------------------------------------------------

template <typename T>
struct MagicHeader<GenBank, T>
{
    static char const VALUE[6];
};
template <typename T>
char const MagicHeader<GenBank, T>::VALUE[6] = { 'L','O','C','U','S',' ' };  // typical GenBank header

// --------------------------------------------------------------------------
// Metafunction FileExtensions
// --------------------------------------------------------------------------

template <typename T>
struct FileExtensions<GenBank, T>
{
    static char const * VALUE[1];
};
template <typename T>
char const * FileExtensions<GenBank, T>::VALUE[1] =
{
    ".gbk",     // default output extension
};

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function splitGenBankHeader()
// ----------------------------------------------------------------------------

/*!
 * @fn splitGenBankHeader
 * @headerfile <seqan/seq_io.h>
 * @brief Split a GenBank header field/value.
 *
 * @signature void readGenBankHeader(key, value, iter);
 *
 * @param[out] key   A @link ContainerConcept @endlink object to write the key to.
 * @param[out] value A @link ContainerConcept @endlink object to write the value to.
 * @param[in]  iter  Input iterator.
 */

// readGenBankHeader() for GenBank header records.

template <typename TKey, typename TValue, typename TFwdIterator>
inline void
readRecord(TKey & key, TValue & val, TFwdIterator & iter, GenBankHeader)
{
    clear(key);
    clear(val);

    readUntil(key, iter, IsWhitespace());

    if (IsBlank()(value(iter)))
    {
        skipUntil(iter, NotFunctor<IsBlank>());
        readLine(val, iter);
    }

    while (IsBlank()(value(iter)))
    {
        appendValue(val, '\n');
        readLine(val, iter);
    }
}

// ----------------------------------------------------------------------------
// Function nextIs()
// ----------------------------------------------------------------------------

// nextIs() for GenBank header records.

template <typename TFwdIterator>
inline bool
nextIs(TFwdIterator & iter, GenBankHeader)
{
    return !atEnd(iter) && value(iter) != 'O' && !IsWhitespace()(value(iter));
}

// nextIs() for GenBank sequence records.

template <typename TFwdIterator>
inline bool
nextIs(TFwdIterator & iter, GenBankSequence)
{
    return !atEnd(iter) && value(iter) == 'O';
}

// ----------------------------------------------------------------------------
// Function readRecord()
// ----------------------------------------------------------------------------

// readRecord() for GenBank sequences records.

// Read all sequence, eat/ignore '//' line.

template <typename TSeqString, typename TFwdIterator>
inline void
readRecord(TSeqString & seq, TFwdIterator & iter, GenBankSequence)
{
    typedef typename Value<TSeqString>::Type TSeqAlphabet;
    typedef OrFunctor<OrFunctor<IsBlank, IsDigit>,
                      AssertFunctor<IsInAlphabet<TSeqAlphabet>, ParseError, Embl> > TSeqAsserter;

    IsNewline isNewline;
    TSeqAsserter asserter;

    if (!atEnd(iter) && value(iter) == 'O')
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

// readRecord() for GenBank id/seq pairs.

template <typename TIdString, typename TSeqString, typename TFwdIterator>
inline void
readRecord(TIdString & meta, TSeqString & seq, TFwdIterator & iter, GenBank)
{
    IsWhitespace isWhite;
    IsBlank isBlank;
    String<char, Array<20> > key;

    // extract meta from "ID" line
    clear(meta);
    for (; !atEnd(iter); skipLine(iter))
    {
        if (isBlank(value(iter)))
            continue;

        AssertFunctor<NotFunctor<CountDownFunctor<True, 20> >, ParseError, Embl> asserter;

        clear(key);
        readUntil(key, iter, isWhite, asserter);

        if (key == "VERSION")
        {
            skipUntil(iter, NotFunctor<IsBlank>());
            readUntil(meta, iter, IsNewline());
        }
        else if (key == "ORIGIN")
        {
            skipLine(iter);
            break;
        }
    }

    readRecord(seq, iter, GenBankSequence());
}

template <typename TIdString, typename TSeqString, typename TQualString, typename TFwdIterator>
inline void
readRecord(TIdString & meta, TSeqString & seq, TQualString & qual, TFwdIterator & iter, GenBank)
{
    clear(qual);
    readRecord(meta, seq, iter, GenBank());
}

}  // namespace seqan

#endif  // #ifndef INCLUDE_SEQAN_SEQ_IO_READ_GENBANK_H_
