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

#ifndef SEQAN_INCLUDE_SEQAN_ALIGN_GLOBAL_ALIGNMENT_MYERS_HIRSCHBERG_IMPL_H_
#define SEQAN_INCLUDE_SEQAN_ALIGN_GLOBAL_ALIGNMENT_MYERS_HIRSCHBERG_IMPL_H_

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
// Function _writeDebugMatrix()
// ----------------------------------------------------------------------------

#ifdef MYERS_HIRSCHBERG_VERBOSE
    template<typename TSource>
    void _writeDebugMatrix(TSource s1,TSource s2)
    {
//IOREV _notio_ not relevant for iorev
        int l1 = length(s1);
        int l2 = length(s2);

        int i,j,sg,sd;

        String<String<int> > fMatrix,rMatrix,tMatrix;

        resize(fMatrix,l1 + 1);
        resize(rMatrix,l1 + 1);
        resize(tMatrix,l1 + 1);

        for(i = 0;i <= l1;++i)
        {
            resize(fMatrix[i],l2 + 1);
            resize(rMatrix[i],l2 + 1);
            resize(tMatrix[i],l2 + 1);
        }

        for(i = 0;i <= l1;++i)
            fMatrix[i][0] = i * (-1);

        for(i = l1;i >= 0;--i)
            rMatrix[i][l2] = (l1 - i) * (-1);

        // calculate forward matrix
        for(j = 1;j <= l2;++j)
        {
            fMatrix[0][j] = j*(-1);
            for(i = 1;i <= l1;++i)
            {
                sg = -1 + ((fMatrix[i-1][j] > fMatrix[i][j-1]) ? fMatrix[i-1][j] : fMatrix[i][j-1]);
                sd = fMatrix[i-1][j-1] + ((s1[i - 1] == s2[j-1]) ? 0 : -1 );

                fMatrix[i][j] = ((sg > sd) ? sg : sd);
            }
        }

        // calculate reverse matrix
        for(j = l2 - 1;j >= 0;--j)
        {
            rMatrix[l1][j] = (l2 - j)*(-1);
            for(i = l1 - 1;i >= 0;--i)
            {
                sg = -1 + ((rMatrix[i+1][j] > rMatrix[i][j+1]) ? rMatrix[i+1][j] : rMatrix[i][j+1]);
                sd = rMatrix[i+1][j+1] + ((s1[i] == s2[j]) ? 0 : -1 );

                rMatrix[i][j] = ((sg > sd) ? sg : sd);
            }
        }

        // print fMatrix
        std::cout << ";-;";
        for(i = 0;i < l1;++i)
            std::cout << s1[i] << ";";

        std::cout << std::endl << "-;";
        for(j = 0;j <= l2;++j)
        {
            if(j != 0) std::cout << s2[j-1] << ";";
            for(i = 0;i <= l1;++i)
            {
                std::cout << fMatrix[i][j] << ";";
            }
            std::cout << std::endl;
        }
        // print rMatrix
        std::cout << ";";
        for(i = 0;i < l1;++i)
            std::cout << s1[i] << ";";
        std::cout << "-;" << std::endl;

        for(j = 0;j <= l2;++j)
        {
            if(j != l2) std::cout << s2[j] << ";";
            else std::cout << "-;";
            for(i = 0;i <= l1;++i)
            {
                std::cout << rMatrix[i][j] << ";";
            }
            std::cout << std::endl;
        }

        // fill and print target matrix
        std::cout << ";-;";
        for(i = 0;i < l1;++i)
            std::cout << s1[i] << ";";

        std::cout << std::endl << "-;";
        for(j = 0;j <= l2;++j)
        {
            if(j != 0) std::cout << s2[j-1] << ";";
            for(i = 0;i <= l1;++i)
            {
                tMatrix[i][j] = fMatrix[i][j] + rMatrix[i][j];
                std::cout << tMatrix[i][j] << ";";
            }
            std::cout << std::endl;
        }
    }
#endif

// ----------------------------------------------------------------------------
// Function globalAlignment()
// ----------------------------------------------------------------------------

// When using different alphabets, we will internally use the pattern alphabet
// for the comparison.  This means that the text character is converted to the
// pattern alphabet.  Here, the "pattern" is the shorter source sequence.

template <typename TSequenceH, typename TGapsSpecH, typename TSequenceV, typename TGapsSpecV>
int
_globalAlignment(Gaps<TSequenceH, TGapsSpecH> & gapsH,
                 Gaps<TSequenceV, TGapsSpecV> & gapsV,
                 MyersHirschberg const & algorithmTag)
{
    // Switch horizontal and vertical gap roles, gapsV should be the shorter one
    // to fit into less words.
    if (length(source(gapsH)) < length(source(gapsV)))
        return _globalAlignment(gapsV, gapsH, algorithmTag);

    clearGaps(gapsH);
    clearGaps(gapsV);
    clearClipping(gapsH);
    clearClipping(gapsV);

    typedef int TScoreValue;

    // use size of unsigned int as blocksize for bit-vectors
    const unsigned int BLOCK_SIZE = BitsPerValue<unsigned int>::VALUE;

    // saves the score value that will be returned
    TScoreValue score,total_score = 0;

    typedef typename Value<TSequenceV>::Type TPatternAlphabet;
    typedef typename Size<TSequenceH>::Type  TStringSize;

    typedef typename Iterator<TSequenceH const, Rooted>::Type TSequenceHIterator;
    typedef typename Iterator<TSequenceV const, Rooted>::Type TSequenceVIterator;
    typedef Gaps<TSequenceH, TGapsSpecH> TGapsH;
    typedef Gaps<TSequenceV, TGapsSpecV> TGapsV;
    typedef typename Iterator<TGapsH, Rooted>::Type TGapsHIterator;
    typedef typename Iterator<TGapsV, Rooted>::Type TGapsVIterator;

    typedef typename Iterator<Matrix<TScoreValue>, Rooted>::Type TMatrixIterator;

    TGapsHIterator target_0 = begin(gapsH);
    TGapsVIterator target_1 = begin(gapsV);

    TSequenceH const & x = source(gapsH);
    TSequenceV const & y = source(gapsV);

    TStringSize len_x = length(x);
    TStringSize len_y = length(y);

    // string to store the score values for the currently active cell
    String<TScoreValue> c_score;
    resize(c_score, len_x + 1, 0);

    // scoring-scheme specific score values
    TScoreValue score_match = 0;
    TScoreValue score_mismatch = -1;
    TScoreValue score_gap = -1;

    // additional vars
    int i;

    // stack with parts of matrix that have to be processed
    std::stack<HirschbergSet_> to_process;
    HirschbergSet_ target;

    // myers specific vars and preprocessing
    unsigned int patternAlphabetSize = ValueSize<TPatternAlphabet>::VALUE;
    unsigned int blockCount = (len_y + BLOCK_SIZE - 1) / BLOCK_SIZE; // maximal count of blocks

    String<unsigned> VP;
    String<unsigned> VN;
    String<unsigned> forwardBitMask;
    String<unsigned> reverseBitMask;

    resize(VP, blockCount, maxValue<unsigned>());
    resize(VN, blockCount, 0);

    // first bitMask will be constructed from the shorter sequence
    resize(forwardBitMask, patternAlphabetSize * blockCount, 0);
    resize(reverseBitMask, patternAlphabetSize * blockCount, 0);

    // encoding the letters as bit-vectors
    for (unsigned int j = 0; j < len_y; j++){
        forwardBitMask[blockCount * ordValue(getValue(y,j)) + j/BLOCK_SIZE] = forwardBitMask[blockCount * ordValue(getValue(y,j)) + j/BLOCK_SIZE] | 1 << (j%BLOCK_SIZE);
        reverseBitMask[blockCount * ordValue(getValue(y,len_y - j - 1)) + j/BLOCK_SIZE] = reverseBitMask[blockCount * ordValue(getValue(y,len_y - j - 1)) + j/BLOCK_SIZE] | 1 << (j%BLOCK_SIZE);
    }

    HirschbergSet_ hs_complete(0,len_x,0,len_y,1);
    to_process.push(hs_complete);

    while(!to_process.empty())
    {
        target = to_process.top();
        to_process.pop();
        /* if score is zero, the whole part of the sequence can be simply skipped */
        if(_score(target) == 0)
        {
            /* coukd work faster */
            for(i = 0;i < (_end1(target) - _begin1(target));++i)
            {
                ++target_0;
                ++target_1;
            }

#ifdef MYERS_HIRSCHBERG_VERBOSE
            printf("skipped %i to %i in first sequence\n",_begin1(target),_end1(target));
#endif
        }
        else if(_begin1(target) == _end1(target))
        {

#ifdef MYERS_HIRSCHBERG_VERBOSE
            std::cout << "align y " << _begin2(target) << " to " << _end2(target) << std::endl;
            std::cout << "align " << infix(y,_begin2(target),_end2(target)) << std::endl << std::endl;
#endif
            for(i = 0;i < (_end2(target) - _begin2(target));++i)
            {
                insertGap(target_0);
                ++target_0;
                ++target_1;
            }
        }
        else if(_begin2(target) + 1 == _end2(target))
        {
            /* ALIGN */
#ifdef MYERS_HIRSCHBERG_VERBOSE
            std::cout << "align x " << _begin1(target) << " to " << _end1(target) << " and y " << _begin2(target) << " to " << _end2(target) << std::endl;
            std::cout << "align " << infix(x,_begin1(target),_end1(target)) << " and " << infix(y,_begin2(target),_end2(target)) << std::endl << std::endl;
#endif

            TStringSize len_1 = _end1(target) - _begin1(target);
            TStringSize len_2 = _end2(target) - _begin2(target);

            Matrix<TScoreValue> matrix_;

            setDimension(matrix_, 2);
            setLength(matrix_, 0, len_1 + 1);
            setLength(matrix_, 1, len_2 + 1);
            resize(matrix_);

            /* init matrix */
            TSequenceHIterator xs_begin = iter(x,_begin1(target)) - 1;
            TSequenceHIterator xs_end = iter(x,_end1(target)) - 1;
            TSequenceVIterator ys_begin = iter(y,_begin2(target)) - 1;
            TSequenceVIterator ys_end = iter(y,_end2(target)) - 1;

            TSequenceHIterator xs = xs_end;
            TSequenceVIterator ys;

            TMatrixIterator col_ = end(matrix_) - 1;
            TMatrixIterator finger1;
            TMatrixIterator finger2;


            TScoreValue h = 0;
            TScoreValue border_ = score_gap;
            TScoreValue v = border_;


            //-------------------------------------------------------------------------
            // init

            finger1 = col_;
            *finger1 = 0;
            for (xs = xs_end; xs != xs_begin; --xs)
            {
                goPrevious(finger1, 0);
                *finger1 = border_;
                border_ += score_gap;
            }

            //-------------------------------------------------------------------------
            //fill matrix

            border_ = 0;
            for (ys = ys_end; ys != ys_begin; --ys)
            {
                TPatternAlphabet cy = *ys;
                h = border_;
                border_ += score_gap;
                v = border_;

                finger2 = col_;
                goPrevious(col_, 1);
                finger1 = col_;

                *finger1 = v;

                for (xs = xs_end; xs != xs_begin; --xs)
                {
                    goPrevious(finger1, 0);
                    goPrevious(finger2, 0);
                    if (*xs == cy)
                    {
                        v = h + score_match;
                        h = *finger2;
                    }
                    else
                    {
                        TScoreValue s1 = h + score_mismatch;
                        h = *finger2;
                        TScoreValue s2 = score_gap + ((h > v) ? h : v);
                        v = (s1 > s2) ? s1 : s2;
                    }
                    *finger1 = v;
                }
            }

            // if computed the whole matrix last value of v = alignment score
            if(target == hs_complete)   total_score = v;

            /* TRACE BACK */
            finger1 = begin(matrix_);
            xs = iter(x,_begin1(target));
            ys = iter(y,_begin2(target));
            xs_end = iter(x,_end1(target));
            ys_end = iter(y,_end2(target));

            while ((xs != xs_end) && (ys != ys_end))
            {
                bool gv;
                bool gh;

                if (*xs == *ys)
                {
                    gv = gh = true;
                }
                else
                {
                    TMatrixIterator it_ = finger1;

                    goNext(it_, 0);
                    TScoreValue v = *it_;

                    goNext(it_, 1);
                    TScoreValue d = *it_;

                    it_ = finger1;
                    goNext(it_, 1);
                    TScoreValue h = *it_;

                    gv = (v >= h) | (d >= h);
                    gh = (h >= v) | (d >= v);
                }

                if (gv)
                {
                    ++xs;
                    goNext(finger1, 0);
                }
                else
                {
                    insertGap(target_0);
                }

                if (gh)
                {
                    ++ys;
                    goNext(finger1, 1);
                }
                else
                {
                    insertGap(target_1);
                }

                ++target_0;
                ++target_1;
            }

            // if x or y did not reached there end position, fill the rest with gaps
            while(xs != xs_end)
            {
                insertGap(target_1);
                ++target_0;
                ++target_1;
                ++xs;
            }

            while(ys != ys_end)
            {
                insertGap(target_0);
                ++target_0;
                ++target_1;
                ++ys;
            }
            /* END ALIGN */


#ifdef MYERS_HIRSCHBERG_VERBOSE
            std::cout << std::endl << align_ << std::endl << std::endl;
#endif

        }
        else
        {
            /*
                ---------------------------------------------------------------
                Calculate cut position using extended Myers-Bitvector-Algorithm
                ---------------------------------------------------------------
            */

            /* declare variables */
            unsigned int X, D0, HN, HP;

            /* compute cut position */
            int mid = static_cast<int>(floor( static_cast<double>((_begin2(target) + _end2(target))/2) ));

            /* debug infos */
#ifdef MYERS_HIRSCHBERG_VERBOSE
            std::cout << "calculate cut for x " << _begin1(target) << " to " << _end1(target) << " and y " << _begin2(target) << " to " << _end2(target) << std::endl;
            std::cout << "calculate cut for " << infix(x,_begin1(target),_end1(target)) << " and " << infix(y,_begin2(target),_end2(target)) << std::endl;
            std::cout << "cut is in row " << mid << " symbol is " << getValue(x,mid-1) << std::endl << std::endl;

            std::cout << std::endl;
            _writeDebugMatrix(infix(x,_begin1(target),_end1(target)),infix(y,_begin2(target),_end2(target)));
            std::cout << std::endl;
#endif
            /* compute blocks and score masks */
            int fStartBlock = _begin2(target) / BLOCK_SIZE;
            int fEndBlock = (mid - 1) / BLOCK_SIZE;
            int fSpannedBlocks = (fEndBlock - fStartBlock) + 1;

            unsigned int fScoreMask = 1 << ((mid  - 1) % BLOCK_SIZE);

            unsigned int fOffSet = _begin2(target) % BLOCK_SIZE;
            unsigned int fSilencer = ~0;
            fSilencer <<= fOffSet;

            /* reset v-bitvectors */
            std::fill(begin(VP, Standard()) + fStartBlock, end(VP, Standard()) + fEndBlock + 1, maxValue<unsigned>());
            std::fill(begin(VN, Standard()) + fStartBlock, end(VN, Standard()) + fEndBlock + 1, 0);

            /* determine start-position and start-score */
            int pos = _begin1(target);
            score = (mid - _begin2(target)) * score_gap;
            c_score[pos] = score;

            /* compute with myers - forward - begin */
            if(fSpannedBlocks == 1)
            {
                while (pos < _end1(target)) {
                    X = (fSilencer & forwardBitMask[(blockCount * ordValue(static_cast<TPatternAlphabet>(getValue(x,pos)))) + fStartBlock]) | VN[fStartBlock];

                    D0 = ((VP[fStartBlock] + (X & VP[fStartBlock])) ^ VP[fStartBlock]) | X;
                    HN = VP[fStartBlock] & D0;
                    HP = VN[fStartBlock] | ~(VP[fStartBlock] | D0);

                    X = (HP << 1) | (1 << fOffSet);
                    VN[fStartBlock] = X & D0;
                    VP[fStartBlock] = (HN << 1) | ~(X | D0);

                    if (HP & fScoreMask)
                        score--;
                    else if (HN & fScoreMask)
                        score++;

                    c_score[pos + 1] = score;

                    ++pos;
                }
            } /* end - short patten */
            else
            {
                int shift, currentBlock;
                unsigned int temp, carryD0, carryHP, carryHN;

                while (pos < _end1(target))
                {
                    carryD0 = carryHP = carryHN = 0;
                    shift = blockCount * ordValue(static_cast<TPatternAlphabet>(getValue(x,pos)));

                    // computing first the top most block
                    X = (fSilencer & forwardBitMask[shift + fStartBlock]) | VN[fStartBlock];

                    temp = VP[fStartBlock] + (X & VP[fStartBlock]);
                    carryD0 = temp < VP[fStartBlock];

                    D0 = (temp ^ VP[fStartBlock]) | X;
                    HN = VP[fStartBlock] & D0;
                    HP = VN[fStartBlock] | ~(VP[fStartBlock] | D0);

                    X = (HP << 1) | (1 << fOffSet);
                    carryHP = HP >> (BLOCK_SIZE - 1);

                    VN[fStartBlock] = X & D0;

                    temp = (HN << 1);
                    carryHN = HN >> (BLOCK_SIZE - 1);

                     VP[fStartBlock] = temp | ~(X | D0);

                    // compute the remaining blocks
                    for (currentBlock = fStartBlock + 1; currentBlock <= fEndBlock; currentBlock++) {
                        X = forwardBitMask[shift + currentBlock] | VN[currentBlock];

                        temp = VP[currentBlock] + (X & VP[currentBlock]) + carryD0;

                        carryD0 = ((carryD0) ? temp <= VP[currentBlock] : temp < VP[currentBlock]);

                        D0 = (temp ^ VP[currentBlock]) | X;
                        HN = VP[currentBlock] & D0;
                        HP = VN[currentBlock] | ~(VP[currentBlock] | D0);

                        X = (HP << 1) | carryHP;
                        carryHP = HP >> (BLOCK_SIZE-1);

                        VN[currentBlock] = X & D0;

                        temp = (HN << 1) | carryHN;
                        carryHN = HN >> (BLOCK_SIZE - 1);

                         VP[currentBlock] = temp | ~(X | D0);
                    }

                    /* update score */
                    if (HP & fScoreMask)
                        score--;
                    else if (HN & fScoreMask)
                        score++;

                    c_score[pos + 1] = score;

                    ++pos;
                }

            } /* end - long patten */
            /* compute with myers - forward - end */

            /* compute blocks and score masks */
            int rStartBlock = (len_y - _end2(target)) / BLOCK_SIZE;
            int rEndBlock = (len_y - mid - 1) / BLOCK_SIZE;
            int rSpannedBlocks = (rEndBlock - rStartBlock) + 1;

            unsigned int rScoreMask = 1 <<  ((len_y - mid - 1) % BLOCK_SIZE);
            unsigned int rOffSet = (len_y - _end2(target)) % BLOCK_SIZE;
            unsigned int rSilencer = ~0;
            rSilencer <<= rOffSet;

            /* reset v-bitvectors */
            std::fill(begin(VP, Standard()) + rStartBlock, end(VP, Standard()) + rEndBlock + 1, maxValue<unsigned>());
            std::fill(begin(VN, Standard()) + rStartBlock, end(VN, Standard()) + rEndBlock + 1, 0);

            /* determine start-position and start-score */
            pos = _end1(target)-1;
            score = (_end2(target) - mid) * score_gap;

            /* set start score */
            c_score[_end1(target)] += score;

            /* determine optimal cut position -- score extension */
            TScoreValue max = c_score[_end1(target)];
            TScoreValue rmax = score;
            unsigned int pos_max = _end1(target);

            /* compute with myers - reverse - begin */
            if(rSpannedBlocks == 1)
            {
                while (pos >= _begin1(target)) {
                    X = (rSilencer & reverseBitMask[(blockCount * ordValue(static_cast<TPatternAlphabet>(getValue(x,pos)))) + rStartBlock]) | VN[rStartBlock];

                    D0 = ((VP[rStartBlock] + (X & VP[rStartBlock])) ^ VP[rStartBlock]) | X;
                    HN = VP[rStartBlock] & D0;
                    HP = VN[rStartBlock] | ~(VP[rStartBlock] | D0);

                    X = (HP << 1) | (1 << rOffSet);
                    VN[rStartBlock] = X & D0;
                    VP[rStartBlock] = (HN << 1) | ~(X | D0);

                    if (HP & rScoreMask)
                        --score;
                    else if (HN & rScoreMask)
                        ++score;

                    c_score[pos] += score;

                    /* check for optimality -- score extension */
                    if(c_score[pos]> max)
                    {
                        pos_max = pos;
                        max = c_score[pos];
                        rmax =  score;
                    }

                    --pos;
                }
            } /* end - short pattern */
            else
            {
                int shift, currentBlock;
                unsigned int temp, carryD0, carryHP, carryHN;

                while (pos >= _begin1(target))
                {
                    carryD0 = carryHP = carryHN = 0;
                    shift = blockCount * ordValue(static_cast<TPatternAlphabet>(getValue(x,pos)));

                    // compute first the top most block
                    X = (rSilencer & reverseBitMask[shift + rStartBlock]) | VN[rStartBlock];

                    temp = VP[rStartBlock] + (X & VP[rStartBlock]);
                    carryD0 = temp < VP[rStartBlock];

                    D0 = (temp ^ VP[rStartBlock]) | X;
                    HN = VP[rStartBlock] & D0;
                    HP = VN[rStartBlock] | ~(VP[rStartBlock] | D0);

                    X = (HP << 1) | (1 << rOffSet);
                    carryHP = HP >> (BLOCK_SIZE - 1);

                    VN[rStartBlock] = X & D0;

                    temp = (HN << 1);
                    carryHN = HN >> (BLOCK_SIZE - 1);

                     VP[rStartBlock] = temp | ~(X | D0);

                    // compute the remaining blocks
                    for (currentBlock = rStartBlock + 1; currentBlock <= rEndBlock; currentBlock++) {
                        X = reverseBitMask[shift + currentBlock] | VN[currentBlock];

                        temp = VP[currentBlock] + (X & VP[currentBlock]) + carryD0;

                        carryD0 = ((carryD0) ? temp <= VP[currentBlock] : temp < VP[currentBlock]);

                        D0 = (temp ^ VP[currentBlock]) | X;
                        HN = VP[currentBlock] & D0;
                        HP = VN[currentBlock] | ~(VP[currentBlock] | D0);

                        X = (HP << 1) | carryHP;
                        carryHP = HP >> (BLOCK_SIZE-1);

                        VN[currentBlock] = X & D0;

                        temp = (HN << 1) | carryHN;
                        carryHN = HN >> (BLOCK_SIZE - 1);

                         VP[currentBlock] = temp | ~(X | D0);
                    }

                    if (HP & rScoreMask)
                        --score;
                    else if (HN & rScoreMask)
                        ++score;

                    c_score[pos] += score;

                    /* check for optimality -- score extension*/
                    if(c_score[pos] > max)
                    {
                        pos_max = pos;
                        max = c_score[pos];
                        rmax = score;
                    }

                    --pos;
                }

            }  /* end - long pattern */
            /* compute with myers - reverse - end */

            // if computed the whole matrix max = alignment score
            if(target == hs_complete)
                total_score = max;

#ifdef MYERS_HIRSCHBERG_VERBOSE
            printf("Optimal cut is at %i and %i with forward score %i and reverse score %i\n\n",mid,pos_max,(max - rmax),rmax);
#endif
            /* push the two computed parts of the dp-matrix on process stack */
            to_process.push(HirschbergSet_(pos_max,_end1(target),mid,_end2(target),rmax));
            to_process.push(HirschbergSet_(_begin1(target),pos_max,_begin2(target),mid,max - rmax));

        }
        /* END CUT */
    }

    return total_score;
}

}  // namespace seqan

#endif  // #ifndef SEQAN_INCLUDE_SEQAN_ALIGN_GLOBAL_ALIGNMENT_MYERS_HIRSCHBERG_IMPL_H_
