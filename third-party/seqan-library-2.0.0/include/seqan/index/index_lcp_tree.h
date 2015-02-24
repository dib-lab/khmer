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

#ifndef SEQAN_HEADER_INDEX_LCP_TREE_H
#define SEQAN_HEADER_INDEX_LCP_TREE_H

namespace SEQAN_NAMESPACE_MAIN
{

    template <
        class LCPFwdIt,        // lcp table input iterator
        class FlatOutIt >    // flat tree output iterator
    inline FlatOutIt createLcpBinTree(
        LCPFwdIt First_, LCPFwdIt _Last,
        FlatOutIt Dest_)
    {
        typedef typename Value<LCPFwdIt>::Type  TValue;
        typedef typename Size<LCPFwdIt>::Type   TSize;

        TSize size = difference(First_, _Last);
        if (size <= 1) return Dest_;
        --size;

        // calculate the depth of the lcp tree
        unsigned treeLevels = 1;
        TSize _xSize = 1;
        for(; size > _xSize; _xSize *= 2, ++treeLevels) ;

        // get output iterators for every level in the flat tree
        FlatOutIt *level = new FlatOutIt[treeLevels];
        for(unsigned i = treeLevels - 1; _xSize; --i, _xSize /= 2) {
            level[i] = Dest_;
            goFurther(Dest_, (size + _xSize - 1) / _xSize);
        }

        // fields to keep track of minimum elements and state
        TValue *minVal = new TValue[treeLevels];
        bool *half = new bool[treeLevels];
        for(unsigned i = 0; i < treeLevels; ++i)
            half[i] = false;

        // it works like a binary counter of half[n]...half[1]half[0]
        for(TSize j = 0; j < size; ++j, ++First_) {
            *(level[0]) = minVal[0] = *First_;
            ++(level[0]);
            for(unsigned i = 1; i < treeLevels; ++i) {
                if (half[i]) {
                    if (minVal[i-1] < minVal[i]) minVal[i] = minVal[i-1];
                    *(level[i]) = minVal[i];    // min[i] is the minimum of last 2 values in min[i-1]
                    ++(level[i]);
                    half[i] = false;
                } else {
                    minVal[i] = minVal[i-1];
                    half[i] = true;
                    break;
                }
            }
        }

        // complete half filled nodes
        bool carry = false;
        for(unsigned i = 1; i < treeLevels; ++i)
            if (half[i] || carry) {
                if (half[i]) {
                    if (minVal[i-1] < minVal[i]) minVal[i] = minVal[i-1];
                } else
                    minVal[i] = minVal[i-1];
                *(level[i]) = minVal[i];
                ++(level[i]);
                carry = true;
            }

        // trailing zero
        *Dest_ = 0;
        ++Dest_;

        delete[] half;
        delete[] minVal;
        delete[] level;

        return Dest_;
    }



    template < typename TSize >
    inline TSize sizeofLcpe(TSize n)
    {
        if (n < 2) return n;    // 0 -> 0, 1 -> 1, 2 -> 2, 3 -> 4
        --n;
        TSize size = 2;
        for(TSize _xSize = 1; _xSize < n; _xSize *= 2)
            size += (n + _xSize - 1) / _xSize;
        return size;
    }

    template < typename TSize >
    inline TSize sizeofLcph(TSize n)
    {
        return sizeofLcpe(n);
    }

    template <
        class LCPFwdIt,        // lcp table input iterator
        typename TSize >
    inline void sizeofLcpe(LCPFwdIt First_, LCPFwdIt _Last, TSize &Size_)
    {
        Size_ = sizeofLcpe(difference(First_, _Last));
    }

    template <
        class LCPFwdIt,        // lcp table input iterator
        typename TSize >
    inline void sizeofLcph(LCPFwdIt First_, LCPFwdIt _Last, TSize &Size_)
    {
        sizeofLcpe(First_, _Last, Size_);
        return;
    }


    template < typename TLCPE, typename TLCP >
    inline void createLcpBinTree(TLCPE &lcp_enhanced, TLCP &lcp) {
        createLcpBinTree(begin(lcp, Standard()), end(lcp, Standard()), begin(lcp_enhanced, Standard()));
    }


    template < typename TSize >
    inline unsigned _treeLevels(TSize lcpSize)
    {
        unsigned treeLevels = 1;
        --lcpSize;
        TSize _xSize = 1;
        for (; lcpSize > _xSize; _xSize *= 2, ++treeLevels)
            continue;  // Get smallest 2^k that is <= lcpSize.
        return treeLevels;
    }

    template < typename TValue, typename TConfig, typename TLCP >
    inline void createLcpBinTree(String<TValue, External<TConfig> > &lcp_enhanced, TLCP &lcp) {
        unsigned writeHeads = _treeLevels(length(lcp)) + 1;   // plus 1 write back buffer
        if (lcp_enhanced.cache.size() < writeHeads)
            lcp_enhanced.resizeCache(writeHeads);
        createLcpBinTree(begin(lcp, Standard()), end(lcp, Standard()), begin(lcp_enhanced));
    }

}

#endif
