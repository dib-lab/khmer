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

#ifndef SEQAN_HEADER_GRAPH_ALIGN_TCOFFEE_REFINEMENT_H
#define SEQAN_HEADER_GRAPH_ALIGN_TCOFFEE_REFINEMENT_H

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////////

template<typename TScore, typename TSc>
inline void
_matchScore(TScore&,
             TSc)
{
    // No operation
}

//////////////////////////////////////////////////////////////////////////////////

template<typename TScore, typename TSc>
inline void
_mismatchScore(TScore&,
                TSc)
{
    // No operation
}

//////////////////////////////////////////////////////////////////////////////////

template<typename TValue, typename TSc>
inline void
_matchScore(Score<TValue, Simple>& sc,
             TSc msc)
{
    sc.data_match = msc;
}

//////////////////////////////////////////////////////////////////////////////////

template<typename TValue, typename TSc>
inline void
_mismatchScore(Score<TValue, Simple>& sc,
                TSc mmsc)
{
    sc.data_mismatch = mmsc;
}


//////////////////////////////////////////////////////////////////////////////


template<typename T>
struct ProfileProfileScore;
//IOREV _notio_

//////////////////////////////////////////////////////////////////////////////

template <typename TValue, typename TScoreMember>
class Score<TValue, ProfileProfileScore<TScoreMember> >
{
public:
    TScoreMember sc;

    Score() {}

    template<typename TSpec>
    Score(Score<TValue, TSpec>& old_sc)
    {
        sc.data_gap_extend = scoreGapExtend(old_sc);
        sc.data_gap_open = scoreGapOpen(old_sc);
    }

    Score(Score<TValue, Simple>& old_sc)
    {
        sc.data_gap_extend = scoreGapExtend(old_sc);
        sc.data_gap_open = scoreGapOpen(old_sc);
        _matchScore(sc, scoreMatch(old_sc));
        _mismatchScore(sc, scoreMismatch(old_sc));
    }

    Score(TValue gap_extend, TValue gap_open)
    {
        sc.data_gap_extend = gap_extend;
        sc.data_gap_open = gap_open;
    }

    Score(TValue match, TValue mismatch, TValue gap_extend, TValue gap_open)
    {
        sc.data_gap_extend = gap_extend;
        sc.data_gap_open = gap_open;
        _matchScore(sc, match);
        _mismatchScore(sc, mismatch);
    }
};


template <typename TValue, typename TScoreMember, typename TPos1, typename TPos2, typename TSeq1, typename TSeq2>
inline TValue
scoreGapExtendHorizontal(
    Score<TValue, ProfileProfileScore<TScoreMember> > const & me,
    TPos1 pos1,
    TPos2,
    TSeq1 const& seq1,
    TSeq2 const&)
{
    typedef typename Size<TSeq1>::Type TSize;
    typedef typename Value<TSeq1>::Type TAlphabet1;
    //typedef typename Value<TSeq2>::Type TAlphabet2;
    typedef typename SourceValue<TAlphabet1>::Type TSourceValue;

    // Last character is a gap
    // std::cout << seq1[pos1] << std::endl;

    TSize alph_size = ValueSize<TSourceValue>::VALUE;
    TValue n1 = 0;
    for(TSize i = 0; i < alph_size; ++i) n1 += seq1[pos1].count[i];
    TValue totalSc = n1 * scoreGapExtend(me.sc);
    n1 += seq1[pos1].count[alph_size];
    //std::cout <<  totalSc / n1 << ',';
    return totalSc / n1;
}

template <typename TValue, typename TScoreMember, typename TPos1, typename TPos2, typename TSeq1, typename TSeq2>
inline TValue
scoreGapOpenHorizontal(
    Score<TValue, ProfileProfileScore<TScoreMember> > const & me,
    TPos1 pos1,
    TPos2,
    TSeq1 const & seq1,
    TSeq2 const &)
{
    typedef typename Size<TSeq1>::Type TSize;
    typedef typename Value<TSeq1>::Type TAlphabet1;
    //typedef typename Value<TSeq2>::Type TAlphabet2;
    typedef typename SourceValue<TAlphabet1>::Type TSourceValue;

    TSize alph_size = ValueSize<TSourceValue>::VALUE;
    TValue n1 = 0;
    for(TSize i = 0; i < alph_size+1; ++i) n1 += seq1[pos1].count[i];
    TSize newGapOpen = n1;
    if ((pos1) && (seq1[pos1 - 1].count[alph_size] < seq1[pos1].count[alph_size])) {
        newGapOpen -= (seq1[pos1].count[alph_size] - seq1[pos1 - 1].count[alph_size]);
    } else if (!pos1) {
        newGapOpen -= seq1[pos1].count[alph_size];
    }
    TValue totalSc = newGapOpen * scoreGapOpen(me.sc);
    return totalSc / n1;
}


template <typename TValue, typename TScoreMember, typename TPos1, typename TPos2, typename TSeq1, typename TSeq2>
inline TValue
scoreGapExtendVertical(
    Score<TValue, ProfileProfileScore<TScoreMember> > const & me,
    TPos1,
    TPos2 pos2,
    TSeq1 const &,
    TSeq2 const & seq2)
{
    typedef typename Size<TSeq1>::Type TSize;
    typedef typename Value<TSeq1>::Type TAlphabet1;
    //typedef typename Value<TSeq2>::Type TAlphabet2;
    typedef typename SourceValue<TAlphabet1>::Type TSourceValue;

    // Last character is a gap
    // std::cout << seq2[pos2] << std::endl;

    TSize alph_size = ValueSize<TSourceValue>::VALUE;
    TValue n2 = 0;
    for(TSize j = 0; j < alph_size; ++j) n2 += seq2[pos2].count[j];
    TValue totalSc = n2 * scoreGapExtend(me.sc);
    n2 += seq2[pos2].count[alph_size];
    //std::cout <<  totalSc / n2 << ',';
    return totalSc / n2;
}

template <typename TValue, typename TScoreMember, typename TPos1, typename TPos2, typename TSeq1, typename TSeq2>
inline TValue
scoreGapOpenVertical(
    Score<TValue, ProfileProfileScore<TScoreMember> > const & me,
    TPos1,
    TPos2 pos2,
    TSeq1 const &,
    TSeq2 const & seq2)
{
    typedef typename Size<TSeq1>::Type TSize;
    typedef typename Value<TSeq1>::Type TAlphabet1;
    //typedef typename Value<TSeq2>::Type TAlphabet2;
    typedef typename SourceValue<TAlphabet1>::Type TSourceValue;

    TSize alph_size = ValueSize<TSourceValue>::VALUE;
    TValue n2 = 0;
    for(TSize j = 0; j < alph_size+1; ++j) n2 += seq2[pos2].count[j];
    TSize newGapOpen = n2;
    if ((pos2) && (seq2[pos2 - 1].count[alph_size] < seq2[pos2].count[alph_size])) {
        newGapOpen -= (seq2[pos2].count[alph_size] - seq2[pos2 - 1].count[alph_size]);
    } else if (!pos2) {
        newGapOpen -= seq2[pos2].count[alph_size];
    }
    TValue totalSc = newGapOpen * scoreGapOpen(me.sc);
    return totalSc / n2;
}


template <typename TValue, typename TScoreMember, typename TPos1, typename TPos2, typename TSeq1, typename TSeq2>
inline TValue
score(Score<TValue, ProfileProfileScore<TScoreMember> > const & me,
      TPos1 pos1,
      TPos2 pos2,
      TSeq1 const &seq1,
      TSeq2 const &seq2)
{
    typedef typename Size<TSeq1>::Type TSize;
    typedef typename Value<TSeq1>::Type TAlphabet1;
    //typedef typename Value<TSeq2>::Type TAlphabet2;
    typedef typename SourceValue<TAlphabet1>::Type TSourceValue;

    // Last character is a gap
    //std::cout << seq1[pos1] << std::endl;
    //std::cout << seq2[pos2] << std::endl;

    TValue totalSc = 0;
    TSize alph_size = ValueSize<TSourceValue>::VALUE;
    TValue n1 = 0;
    TValue n2 = 0;
    for(TSize j = 0; j < alph_size; ++j) n2 += seq2[pos2].count[j];
    for(TSize i = 0; i < alph_size; ++i) {
        if (seq1[pos1].count[i]) {
            n1 += seq1[pos1].count[i];
            for(TSize j = 0; j < alph_size; ++j) {
                if (seq2[pos2].count[j]) {
                    //std::cout << TSourceValue(i) << ',' << TSourceValue(j) << ',' << score(me.sc, TSourceValue(i), TSourceValue(j)) << std::endl;
                    totalSc += seq1[pos1].count[i] * seq2[pos2].count[j] * score(me.sc, TSourceValue(i), TSourceValue(j));
                }
            }
        }
    }
    totalSc += n2 * seq1[pos1].count[alph_size] * scoreGapExtend(me.sc);
    totalSc += n1 * seq2[pos2].count[alph_size] * scoreGapExtend(me.sc);
    n1 += seq1[pos1].count[alph_size];
    n2 += seq2[pos2].count[alph_size];
    //std::cout <<  totalSc / (n1 * n2) << ',';
    return totalSc / (n1 * n2);
}


//////////////////////////////////////////////////////////////////////////////

template<typename TValue, typename TSpec, typename TSize, typename TAlphabet, typename TScore>
inline void
_msaRefinement(String<TValue, TSpec>& mat,
               TSize nseq,
               TSize splitPos,
               TScore& sc,
               TAlphabet)
{
    // Is it a valid split?
    if ((splitPos == 0) || (splitPos > nseq - 1)) return;

    // Initialization
    typedef String<TValue, TSpec> TAlignMat;
    typedef typename Iterator<TAlignMat, Standard>::Type TMatIter;
    TSize alignLen = length(mat) / nseq;

    //// Debug code
    //for(TSize row = 0; row < nseq; ++row) {
    //    for(TSize col = 0; col < alignLen; ++col) {
    //        std::cout << mat[row * alignLen + col];
    //    }
    //    std::cout << std::endl;
    //}

    char gapChar = gapValue<char>();
    TSize nseq1 = splitPos;
    TSize nseq2 = nseq - nseq1;
    TSize gapPos = ValueSize<TAlphabet>::VALUE;
    typedef ProfileChar<TAlphabet> TProfile;
    typedef String<TProfile> TProfileString;
    typedef typename Iterator<TProfileString, Standard>::Type TProfIter;
    TProfileString prof1;
    TProfileString prof2;
    resize(prof1, alignLen, TProfile());
    resize(prof2, alignLen, TProfile());
    TMatIter itMat = begin(mat, Standard());
    TMatIter itMatEnd = begin(mat, Standard());
    itMatEnd += nseq1 * alignLen;
    TProfIter itProf = begin(prof1, Standard());
    for(TSize i = 0;itMat != itMatEnd;++itMat, ++itProf, ++i) {
        if (i % alignLen == 0) itProf = begin(prof1, Standard());
        if (*itMat == gapChar) ++(itProf->count[gapPos]);
        else ++(itProf->count[ordValue((TAlphabet) *itMat)]);
    }
    itProf = begin(prof2, Standard());
    itMatEnd = end(mat, Standard());
    for(TSize i = 0;itMat != itMatEnd;++itMat, ++itProf, ++i) {
        if (i % alignLen == 0) itProf = begin(prof2, Standard());
        if (*itMat == gapChar) ++(itProf->count[gapPos]);
        else ++(itProf->count[ordValue((TAlphabet) *itMat)]);
    }


    TProfileString alignProf1;
    reserve(alignProf1, alignLen);
    itProf = begin(prof1, Standard());
    TProfIter itProfEnd = end(prof1, Standard());
    typedef String<TSize> TDelGap;
    typedef typename Iterator<TDelGap, Standard>::Type TDelGapIter;
    TDelGap delGapCol;
    for(TSize i = 0;itProf != itProfEnd;++itProf, ++i) {
        if (itProf->count[gapPos] != nseq1) appendValue(alignProf1, *itProf);
        else appendValue(delGapCol, i, Generous());
    }
    TDelGapIter itDelGap = begin(delGapCol, Standard());
    TDelGapIter itDelGapEnd = end(delGapCol, Standard());
    TAlignMat mat1;
    TSize alignLen1 = alignLen - length(delGapCol);
    resize(mat1, nseq1 * alignLen1);
    itMat = begin(mat, Standard());
    itMatEnd = begin(mat, Standard());
    itMatEnd += nseq1 * alignLen;
    TMatIter itMat1 = begin(mat1, Standard());
    for(TSize col = 0;itMat != itMatEnd;++itMat, ++col) {
        if (col % alignLen == 0) itDelGap = begin(delGapCol, Standard());
        if ((itDelGap != itDelGapEnd) && (col % alignLen == *itDelGap)) {
            ++itDelGap;
            continue;
        }
        *itMat1 = *itMat;
        ++itMat1;
    }

    //// Debug code
    //for(TSize row = 0; row < nseq1; ++row) {
    //    for(TSize col = 0; col < alignLen1; ++col) {
    //        std::cout << mat1[row * alignLen1 + col];
    //    }
    //    std::cout << std::endl;
    //}


    TProfileString alignProf2;
    reserve(alignProf2, alignLen);
    itProf = begin(prof2, Standard());
    itProfEnd = end(prof2, Standard());
    clear(delGapCol);
    for(TSize i = 0;itProf != itProfEnd;++itProf, ++i) {
        if (itProf->count[gapPos] != nseq2) appendValue(alignProf2, *itProf);
        else appendValue(delGapCol, i, Generous());
    }
    itDelGap = begin(delGapCol, Standard());
    itDelGapEnd = end(delGapCol, Standard());
    TAlignMat mat2;
    TSize alignLen2 = alignLen - length(delGapCol);
    resize(mat2, nseq2 * alignLen2);
    itMatEnd = end(mat, Standard());
    TMatIter itMat2 = begin(mat2, Standard());
    for(TSize col = 0;itMat != itMatEnd;++itMat, ++col) {
        if (col % alignLen == 0) itDelGap = begin(delGapCol, Standard());
        if ((itDelGap != itDelGapEnd) && (col % alignLen == *itDelGap)) {
            ++itDelGap;
            continue;
        }
        *itMat2 = *itMat;
        ++itMat2;
    }

    //// Debug code
    //for(TSize row = 0; row < nseq2; ++row) {
    //    for(TSize col = 0; col < alignLen2; ++col) {
    //        std::cout << mat2[row * alignLen2 + col];
    //    }
    //    std::cout << std::endl;
    //}


    // ReAlign the consensus with the sequence
    typedef StringSet<TProfileString, Dependent<> > TStringSet;
    TStringSet pairSet;
    appendValue(pairSet, alignProf1);
    appendValue(pairSet, alignProf2);

    typedef String<Fragment<> > TFragmentString;
    TFragmentString matches;
    globalAlignment(matches, pairSet, sc, Gotoh());

    typedef typename Iterator<TFragmentString, Standard>::Type TFragIter;
    TFragIter fragIt = end(matches, Standard() );
    TFragIter fragItEnd = begin(matches, Standard() );
    typedef String<TSize> TInsertPattern;
    typedef typename Iterator<TInsertPattern, Standard>::Type TInsertIter;
    TInsertPattern insertPattern1;
    TInsertPattern insertPattern2;
    TSize begin1Pos = 0;
    TSize begin2Pos = 0;
    if (fragIt != fragItEnd) {
        do {
            --fragIt;
            appendValue(insertPattern1, fragIt->begin2 - begin2Pos, Generous());
            appendValue(insertPattern1, fragIt->begin1 - begin1Pos, Generous());
            appendValue(insertPattern1, fragIt->len, Generous());
            appendValue(insertPattern2, fragIt->begin2 - begin2Pos, Generous());
            appendValue(insertPattern2, fragIt->begin1 - begin1Pos, Generous());
            appendValue(insertPattern2, fragIt->len, Generous());
            begin1Pos = fragIt->begin1 + fragIt->len;
            begin2Pos = fragIt->begin2 + fragIt->len;
        } while (fragIt != fragItEnd);
    }
    appendValue(insertPattern1, alignLen2 - begin2Pos, Generous());
    appendValue(insertPattern1, alignLen1 - begin1Pos, Generous());
    appendValue(insertPattern2, alignLen2 - begin2Pos, Generous());
    appendValue(insertPattern2, alignLen1 - begin1Pos, Generous());
    clear(matches);

    TInsertIter itInsert = begin(insertPattern1, Standard());
    TInsertIter itInsertEnd = end(insertPattern1, Standard());
    TSize totalAlignLen = 0;
    for(;itInsert != itInsertEnd; ++itInsert) totalAlignLen += *itInsert;
    itMat1 = begin(mat1, Standard());
    TMatIter itMat1End = end(mat1, Standard());
    clear(mat);
    while(itMat1 != itMat1End) {
        itInsert = begin(insertPattern1, Standard());
        for(TSize i = 0; i < *itInsert; ++i) appendValue(mat, gapChar, Generous());
        ++itInsert;
        for(TSize i = 0; i < *itInsert; ++i, ++itMat1) appendValue(mat, *itMat1, Generous());
        ++itInsert;
        while(itInsert!=itInsertEnd) {
            for(TSize i = 0; i < *itInsert; ++i, ++itMat1) appendValue(mat, *itMat1, Generous());
            ++itInsert;
            for(TSize i = 0; i < *itInsert; ++i) appendValue(mat, gapChar, Generous());
            ++itInsert;
            for(TSize i = 0; i < *itInsert; ++i, ++itMat1) appendValue(mat, *itMat1, Generous());
            ++itInsert;
        }
    }

    itInsertEnd = end(insertPattern2, Standard());
    itMat2 = begin(mat2, Standard());
    TMatIter itMat2End = end(mat2, Standard());
    while(itMat2 != itMat2End) {
        itInsert = begin(insertPattern2, Standard());
        for(TSize i = 0; i < *itInsert; ++i, ++itMat2) appendValue(mat, *itMat2, Generous());
        ++itInsert;
        for(TSize i = 0; i < *itInsert; ++i) appendValue(mat, gapChar, Generous());
        ++itInsert;
        while(itInsert!=itInsertEnd) {
            for(TSize i = 0; i < *itInsert; ++i, ++itMat2) appendValue(mat, *itMat2, Generous());
            ++itInsert;
            for(TSize i = 0; i < *itInsert; ++i, ++itMat2) appendValue(mat, *itMat2, Generous());
            ++itInsert;
            for(TSize i = 0; i < *itInsert; ++i) appendValue(mat, gapChar, Generous());
            ++itInsert;
        }
    }

    //// Debug code
    //for(TSize row = 0; row < nseq; ++row) {
    //    for(TSize col = 0; col < totalAlignLen; ++col) {
    //        std::cout << mat[row * totalAlignLen + col];
    //    }
    //    std::cout << std::endl;
    //}
}



template<typename TStringSet, typename TCargo, typename TSpec, typename TScore>
inline void
msaRefinement(Graph<Alignment<TStringSet, TCargo, TSpec> >& gAlign,
              TScore& sc)
{
    typedef Graph<Alignment<TStringSet, TCargo, TSpec> > TGraph;
    typedef typename Size<TGraph>::Type TSize;
    typedef typename Value<TStringSet>::Type TSequence;
    typedef typename Value<TSequence>::Type TAlphabet;
    TSize nseq = length(stringSet(gAlign));
    String<char> mat;
    convertAlignment(gAlign, mat);
    for(TSize i = 1; i<nseq; ++i) _msaRefinement(mat, nseq, i, sc, TAlphabet());
    convertAlignment(mat, gAlign);
}

}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
