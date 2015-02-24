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
// Author: David Weese <david.weese@fu-berlin.de>
// ==========================================================================

#ifndef SEQAN_HEADER_INDEX_MANUAL_FORWARDS_H
#define SEQAN_HEADER_INDEX_MANUAL_FORWARDS_H

//SEQAN_NO_GENERATED_FORWARDS: no forwards are generated for this file

//////////////////////////////////////////////////////////////////////////////
// CLASSES
//////////////////////////////////////////////////////////////////////////////

namespace SEQAN_NAMESPACE_MAIN {

    struct FibreText_;        // Original text. Can be a String or a StringSet
    struct FibreRawText_;    // Concatenation of the strings above
    struct FibreSA_;        // suffix array (of raw text with virtual $-delimiters) with Pair entries
    struct FibreRawSA_;    // suffix array with integer entries
    struct FibreSae_;        // suffix array reordered in a b-tree
    struct FibreLcp_;        // lcp table of raw text
    struct FibreLcpe_;        // lcp interval tree
    struct FibreChildtab_;    // childtab (Kurtz et al.) of raw text
    struct FibreBwt_;        // burrows wheeler table of raw text

    typedef Tag<FibreText_> const        FibreText;
    typedef Tag<FibreRawText_> const    FibreRawText;
    typedef Tag<FibreSA_> const        FibreSA;
    typedef Tag<FibreRawSA_> const        FibreRawSA;
    typedef Tag<FibreSae_> const        FibreSae;
    typedef Tag<FibreLcp_> const        FibreLcp;
    typedef Tag<FibreLcpe_> const        FibreLcpe;
    typedef Tag<FibreChildtab_> const    FibreChildtab;
    typedef Tag<FibreBwt_> const        FibreBwt;

}

#endif

