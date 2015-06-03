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

#ifndef SEQAN_HEADER_INDEX_SA_MM_H
#define SEQAN_HEADER_INDEX_SA_MM_H

namespace SEQAN_NAMESPACE_MAIN
{

    struct ManberMyers {};


    //////////////////////////////////////////////////////////////////////////////
    // internal Manber Myers algorithm
    //////////////////////////////////////////////////////////////////////////////

    template < typename TSA,
               typename TText >
    void createSuffixArray(
        TSA &SA,
        TText const &s,
        ManberMyers const &,
        unsigned K,
        unsigned maxdepth)
    {
        typedef typename Value<TSA>::Type    TSize;
        typedef typename Value<TText>::Type    TValue;

        typedef String<TSize, Alloc<> >        TString;
        typedef String<bool, Alloc<> >        TBoolString;

        typedef typename Iterator<TString, Standard>::Type    TIter;

        TSize n = length(s);

        TString ISA;        resize(ISA, n, Exact());
        TString count;
        TBoolString Bh;        resize(Bh, n, Exact());
        TBoolString B2h;    resize(B2h, n, Exact());

        // sort the first character (h=1)
        {
            resize(count, _max((TSize)K, n), Exact());
            TIter it = begin(ISA, Standard());
            for(TSize i = 0; i < n; ++i, ++it)
                *it = i;
            SEQAN_PROMARK("Suffix-Array invertiert");

            radixPass(SA, ISA, s, count, length(count));
            arrayFill(begin(Bh, Standard()), end(Bh, Standard()), false);
        }

        // update bucket borders Bh
        {
            TValue c = 0, d;
            for(TSize i = 0; i < n; ++i) {
                if (c != (d = s[SA[i]])) {
                    c = d;
                    Bh[i] = true;
                }
            }
        }

        // extend h-sort to 2h-sort
        {
            resize(count, n, Exact());
            TSize j, k, l, Ti, cd = 0;

            // set necessary number of rounds cd
            for(TSize h = 1; h < n; h <<= 1, ++cd) ;
            if (maxdepth > 0 && maxdepth < cd)
                cd = maxdepth;

            for(TSize h = 1; cd > 0; h <<= 1, --cd) {

                #ifdef SEQAN_DEBUG_INDEX
                    std::cerr << "[" << cd << "] ";
                #endif
                SEQAN_PROADD(SEQAN_PRODEPTH, 1);
                SEQAN_PROMARK("Beginne Durchlauf");
                arrayFill(begin(count, Standard()), end(count, Standard()), 0);
                arrayFill(begin(B2h, Standard()), end(B2h, Standard()), false);

                l = 0;
                for(TSize i = 0; i < n; ++i) {
                    if (Bh[i]) l = i;
                    ISA[SA[i]] = l;
                }

                SEQAN_ASSERT_GEQ(n, h);
                Ti = n - h;
                TIter p = begin(ISA, Standard()) + Ti;

                SEQAN_ASSERT_LT(*p, n);
                j = count[*p]++;
                *p += j;

                SEQAN_ASSERT_LT(*p, n);
                B2h[*p] = true;

                l = 0;
                for(TSize i = 0; i < n; ++i) {

                    if ((Ti = SA[i]) >= h) {
                        Ti -= h;
                        p = begin(ISA) + Ti;

                        SEQAN_ASSERT_LT(*p, n);
                        j = count[*p]++;
                        *p += j;

                        SEQAN_ASSERT_LT(*p, n);
                        B2h[*p] = true;
                    }

                    if (i + 1 == n || Bh[i + 1]) {
                        for(j = l; j <= i; ++j)
                            if ((Ti = SA[j]) >= h) {
                                Ti -= h;
                                if (B2h[k = ISA[Ti]])
                                    while (++k < n && !Bh[k] && B2h[k])
                                        B2h[k] = false;
                            }
                        l = i + 1;
                    }
                }

                for(TSize i = 0; i < n; ++i) {
                    SA[ISA[i]] = i;
                    Bh[i] |= B2h[i];
                }
            }
        }
        SEQAN_PROSET(SEQAN_PRODEPTH, 0);
        #ifdef SEQAN_DEBUG_INDEX
            std::cerr << std::endl;
        #endif
    }

    // creates suffix array sorted by the first maxLCP chars of suffixes
    template < typename TSA,
               typename TText,
               typename TSize >
    inline void createSuffixArrayPart(
        TSA &SA,
        TText &s,
        ManberMyers const &alg,
        TSize maxLCP,
        unsigned K)
    {
        unsigned depth = 0;
        for(TSize i = 1; i < maxLCP; i*=2) ++depth;
        createSuffixArray(SA, s, alg, K, depth);
    }

    template < typename TSA,
               typename TText,
               typename TSize >
    inline void createSuffixArrayPart(
        TSA &SA,
        TText &s,
        ManberMyers const &alg,
        TSize maxLCP)
    {
        SEQAN_CHECKPOINT;
        createSuffixArrayPart(SA, s, alg, maxLCP, ValueSize< typename Value<TText>::Type >::VALUE);
    }
}

#endif //#ifndef SEQAN_HEADER_...
