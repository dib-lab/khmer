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

#ifndef SEQAN_HEADER_INDEX_SHAWARMA_H
#define SEQAN_HEADER_INDEX_SHAWARMA_H

// adapt external C libraries
extern "C" {

    void ds_ssort(unsigned char *t, int *sa, int n);
    int init_ds_ssort(int adist, int bs_ratio);

}

namespace SEQAN_NAMESPACE_MAIN {

/*    namespace shawarma {

        extern void ds_ssort(unsigned char *t, int *sa, int n);
        extern int init_ds_ssort(int adist, int bs_ratio);

    }
*/
//////////////////////////////////////////////////////////////////////////////
// SeqAn interface

    template <typename TSpec>
    struct Shawarma {};

    struct MSufSort{};            // MSufSort
    struct DivSufSort{};        // DivSufSort
    struct DeepShallow{};        // Deep-Shallow sort
    struct QSufSort{};            // QSufSort

    // WARNING:
    // 1. text value must have char size
    // 2. SA value must have int size
    // 3. Deep-Shallow sort expects overshoot bytes behind the text
    // 4. SA must be contiguous
    //

    template < typename TSA,
               typename TText >
    void createSuffixArray(
        TSA &SA,
        TText const &s,
        Shawarma<DeepShallow> const)
    {
        typedef typename Value<TText>::Type    TValue   SEQAN_TYPEDEF_FOR_DEBUG;
        typedef typename Value<TSA>::Type    TSAValue SEQAN_TYPEDEF_FOR_DEBUG;

        SEQAN_ASSERT_EQ(sizeof(TValue), sizeof(unsigned char));
        SEQAN_ASSERT_EQ(sizeof(TSAValue), sizeof(int));
//        SEQAN_ASSERT(IsContiguous<TSA>::VALUE);

        int overshoot = init_ds_ssort(500, 2000);

        SEQAN_ASSERT_GT(overshoot, 0);
        reserve(s, length(s) + overshoot);
        ds_ssort(
            (unsigned char*)toCString(s),        // text
            (int*)begin(SA, Standard()),        // SA
            length(s));                            // n
    }


} // namespace SEQAN_NAMESPACE_MAIN

#endif // SEQAN_HEADER_INDEX_PIZZACHILI_H
