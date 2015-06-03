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

#ifndef SEQAN_HEADER_INDEX_SA_BTREE_H
#define SEQAN_HEADER_INDEX_SA_BTREE_H

namespace SEQAN_NAMESPACE_MAIN
{

    template <
        class SAFwdIt,        // suffix array input iterator
        class FlatOutIt >    // flat tree output iterator
    inline FlatOutIt createSABTree(
        SAFwdIt First_, SAFwdIt _Last,
        FlatOutIt Dest_, unsigned BlockSize)
    {
        typedef typename Value<SAFwdIt>::Type    TSize;

        TSize size = difference(First_, _Last);
        if (!size) return Dest_;

        // calculate the depth of the sa b-tree
        TSize BlockElements = BlockSize - 1;

        unsigned treeLevels = 1;
        TSize _xSize;
        for(_xSize = 1; _xSize * BlockSize <= size; _xSize *= BlockSize, ++treeLevels) ;

        // get output iterators for every level in the flat tree
        FlatOutIt *level = new FlatOutIt[treeLevels];
        for(int i = treeLevels - 1; _xSize; --i, _xSize /= BlockSize) {
            level[i] = Dest_;
            goFurther(Dest_, ((size / _xSize + BlockSize - 1) / BlockSize) * BlockSize);
        }

        // counter for each b-tree level
        TSize *cnt = new TSize[treeLevels];
        for(unsigned i = 0; i < treeLevels; ++i)
            cnt[i] = 0;

        // distribute to responsible levels
        for(TSize j = 0; j < size; ++j, ++First_)
            for(unsigned i = 0; i < treeLevels; ++i) {
                *(level[i]) = *First_;
                ++(level[i]);
                if (cnt[i] != BlockElements) {
                    ++cnt[i];
                    break;
                } else
                    cnt[i] = 0;
            }

        delete[] cnt;
        delete[] level;

        return Dest_;
    }


    template < typename TSize >
    inline TSize sizeofSABTree(TSize n, unsigned BlockSize)
    {
        TSize size = 0;
        for(TSize _xSize = 1; _xSize <= n; _xSize *= BlockSize)
            size += (n / _xSize + BlockSize - 1) / BlockSize;
        return size * BlockSize;
    }

    template <
        class SAFwdIt,        // suffix array input iterator
        typename TSize >
    inline void sizeofSABTree(SAFwdIt First_, SAFwdIt _Last, TSize &Size_, unsigned BlockSize)
    {
        Size_ = sizeofSABTree(difference(First_, _Last), BlockSize);
    }


    template < typename TSAB, typename TSA >
    inline void createSABTree(TSAB &sa_btree, TSA &sa, unsigned BlockSize) {
        createSABTree(begin(sa), end(sa), begin(sa_btree), BlockSize);
    }


    template < typename TSize >
    inline unsigned treeLevelsSAB(TSize saSize, unsigned BlockSize)
    {
        unsigned treeLevels = 1;
        for (TSize _xSize = 1; _xSize <= saSize; _xSize *= BlockSize, ++treeLevels)
            continue;  // Get smallest 2^k that is > saSize.
        return treeLevels;
    }

    template < typename TValue, typename TConfig, typename TSA >
    inline void createSABTree(String<TValue, External<TConfig> > &sa_btree, TSA &sa, unsigned BlockSize) {
        int writeHeads = treeLevelsSAB(length(sa)) + 1;   // plus 1 write back buffer
        if (sa_btree.cache.size() < writeHeads)
            sa_btree.resizeCache(writeHeads);
        createSABTree(begin(sa), end(sa), begin(sa_btree), BlockSize);
    }

}

#endif
