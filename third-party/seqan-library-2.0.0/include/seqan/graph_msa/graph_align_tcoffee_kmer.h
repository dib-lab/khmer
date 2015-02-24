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

#ifndef SEQAN_HEADER_GRAPH_ALIGN_TCOFFEE_KMER_H
#define SEQAN_HEADER_GRAPH_ALIGN_TCOFFEE_KMER_H

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////
// Simple k-mer counter
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

template<typename TString, typename TTupelString, typename TKTup, typename TAlphabet>
inline void
_getTupelString(TString const& str,
                TTupelString& tupelString,
                TKTup const ktup,
                TAlphabet)
{
    SEQAN_CHECKPOINT
    typedef typename Value<typename Value<TTupelString>::Type>::Type TWord;

    // Alphabet size
    TWord alphabet_size = ValueSize<TAlphabet>::VALUE;

    // Assign a unique number to each k-tupel
    String<TWord> prod;  // Scaling according to position in k-tupel
    resize(prod,ktup);
    for (TWord i=0; i< (TWord) ktup;++i) {
        prod[ktup-i-1] = 1;
        for(TWord j=0;j<i;++j) prod[ktup-i-1] *= alphabet_size;
    }

    TWord len = length(str);
    clear(tupelString);
    if(len < ktup) return;
    resize(tupelString, len-(ktup - 1));
    TWord tupelIndex = 0;
    TWord endTupel = 0;
    tupelString[tupelIndex] = 0;
    for(;endTupel< (TWord) ktup;++endTupel) {
        tupelString[tupelIndex] += (TWord) (ordValue((TAlphabet) str[endTupel])) * prod[endTupel];
    }
    ++tupelIndex;
    for(;endTupel<len;++endTupel) {
        tupelString[tupelIndex] = tupelString[tupelIndex - 1];
        tupelString[tupelIndex] -= (TWord) (ordValue((TAlphabet) str[endTupel - ktup])) * prod[0];
        tupelString[tupelIndex] *= alphabet_size;
        tupelString[tupelIndex] += (TWord) (ordValue((TAlphabet) str[endTupel]));
        ++tupelIndex;
    }
}

//////////////////////////////////////////////////////////////////////////////

template<typename TString, typename TSpec, typename THitMatrix, typename TSize, typename TAlphabet>
inline void
getKmerSimilarityMatrix(StringSet<TString, TSpec> const& strSet,
                        THitMatrix& mat,
                        TSize ktup,
                        TAlphabet)
{
    SEQAN_CHECKPOINT
    typedef TSize TWord;
    typedef String<TWord> TTupelString;
    typedef String<TTupelString> TTupelStringSet;
    typedef typename Value<THitMatrix>::Type TValue;

    // Number of sequences
    TSize nseq = length(strSet);
    TSize alphabet_size = ValueSize<TAlphabet>::VALUE;
    TWord qIndexSize = (TWord) std::pow((double)alphabet_size, (double)ktup);

    // Initialization
    // Matrix for common k-tupels between sequence i and j
    resize(mat, nseq*nseq);

    // Transform the set of strings into a set of strings of k-tupels
    TTupelStringSet tupSet;
    resize(tupSet, length(strSet));
    for(TSize k=0;k<(TSize) length(strSet);++k) _getTupelString(strSet[k], tupSet[k], ktup, TAlphabet());

    // Build for each sequence the q-gram Index and count common hits
    String<TWord> qIndex;
    String<TWord> compareIndex;
    for(TSize k=0;k<nseq;++k) {
        clear(qIndex);
        resize(qIndex, qIndexSize, (TWord) 0, Exact());
        for(TSize i = 0;i < (TSize) length(tupSet[k]);++i) ++qIndex[ tupSet[k][i] ];
        TWord value;
        for (TSize k2=k; k2<nseq; ++k2) {
            clear(compareIndex);
            resize(compareIndex, qIndexSize, (TWord) 0, Exact());
            value = 0;
            for(TSize i = 0;i < (TSize) length(tupSet[k2]);++i) {
                //std::cout << tupSet[k2][i] << "," << compareIndex[ tupSet[k2][i] ] << "," << qIndex[ tupSet[k2][i] ]<< std::endl;
                if (compareIndex[ tupSet[k2][i] ] < qIndex[ tupSet[k2][i] ]) ++value;
                ++compareIndex[ tupSet[k2][i] ];
            }
            mat[k*nseq+k2] = value;
        }
    }

    // Scale counts
    for(TWord row = 0; row < (TWord) nseq; ++row) {
        for(TWord col = row+1; col < (TWord) nseq; ++col) {
            // Fractional common kmer count
            // = Number of common q-grams / Number of possible common q-grams
            TValue minVal = (mat[col*nseq+col] < mat[row*nseq+row]) ? mat[col*nseq+col] : mat[row*nseq+row];
            if (minVal == 0) mat[row*nseq+col] = 0;
            else mat[row*nseq+col] = (mat[row*nseq+col] * SEQAN_DISTANCE_UNITY) / minVal;
            //std::cout << mat[row*nseq+col] << ",";
        }
        //std::cout << std::endl;
    }
}




}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
