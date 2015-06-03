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
// Author: Stephan Aiche <stephan.aiche@fu-berlin.de>
// ==========================================================================

#ifndef INCLUDE_SEQAN_ARG_PARSE_XML_SUPPORT_H_
#define INCLUDE_SEQAN_ARG_PARSE_XML_SUPPORT_H_

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function _xmlEscape()
// ----------------------------------------------------------------------------

/*
 * make sure that the text we put into the XML does not break the XML
 * candidates are
 *  " -> &quot;
 *  ' -> &apos;
 *  & -> &amp;
 *  < -> &lt;
 *  > -> &gt;
 */
template <typename TSequence>
TSequence xmlEscape(TSequence const & original)
{
    typedef typename Iterator<TSequence const, Standard>::Type TSequenceIterator;
    TSequence escaped;
    TSequenceIterator endIter = end(original, Standard());
    for (TSequenceIterator ch  = begin(original, Standard()); ch != endIter; goNext(ch))
    {
        if (value(ch) == '"')
            append(escaped, "&quot;");
        else if (value(ch) == '\'')
            append(escaped, "&apos;");
        else if (value(ch) == '&')
            append(escaped, "&amp;");
        else if (value(ch) == '<')
            append(escaped, "&lt;");
        else if (value(ch) == '>')
            append(escaped, "&gt;");
        else
            appendValue(escaped, *ch);
    }
    return escaped;
}

}  // namespace seqan

#endif  // #ifndef INCLUDE_SEQAN_ARG_PARSE_XML_SUPPORT_H_
