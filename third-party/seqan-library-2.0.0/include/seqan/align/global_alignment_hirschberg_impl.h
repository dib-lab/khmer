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

#include <stack>

// TODO(holtgrew): Get rid of this?
//#define SEQAN_HIRSCHBERG_DEBUG_CUT

#ifndef SEQAN_INCLUDE_SEQAN_ALIGN_GLOBAL_ALIGNMENT_HIRSCHBERG_IMPL_H_
#define SEQAN_INCLUDE_SEQAN_ALIGN_GLOBAL_ALIGNMENT_HIRSCHBERG_IMPL_H_

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

struct Hirschberg_;
typedef Tag<Hirschberg_> Hirschberg;

// ----------------------------------------------------------------------------
// Helper Class HirschbergSet_
// ----------------------------------------------------------------------------

// State for the implementation of Hirschberg's algorithm.

class HirschbergSet_
{
public:
    int x1,x2,y1,y2;
    int score;

    HirschbergSet_()
        : x1(0), x2(0), y1(0), y2(0), score(0)
    {}

    HirschbergSet_(int a1,int a2,int b1,int b2,int sc)
        : x1(a1), x2(a2), y1(b1), y2(b2), score(sc)
    {
        SEQAN_ASSERT_LEQ(a1, a2);
        SEQAN_ASSERT_LEQ(b1, b2);
    }

    HirschbergSet_ &
    operator=(HirschbergSet_ const & other_)
    {
        x1 = other_.x1;
        x2 = other_.x2;
        y1 = other_.y1;
        y2 = other_.y2;
        score = other_.score;
        return *this;
    }

};

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function _begin1()
// ----------------------------------------------------------------------------

inline int&
_begin1(HirschbergSet_ & me) {
    return me.x1;
}


inline int const&
_begin1(HirschbergSet_ const & me) {
    return me.x1;
}

// ----------------------------------------------------------------------------
// Function _setBegin1()
// ----------------------------------------------------------------------------

inline void
_setBegin1(HirschbergSet_ & me, int const & new_begin) {
    me.x1 = new_begin;
}

// ----------------------------------------------------------------------------
// Function _end1()
// ----------------------------------------------------------------------------

inline int&
_end1(HirschbergSet_ & me) {
    return me.x2;
}

inline int const&
_end1(HirschbergSet_ const & me) {
    return me.x2;
}

// ----------------------------------------------------------------------------
// Function _setEnd1()
// ----------------------------------------------------------------------------

inline void
_setEnd1(HirschbergSet_ & me, int const & new_end) {
    me.x2 = new_end;
}

// ----------------------------------------------------------------------------
// Function _begin2()
// ----------------------------------------------------------------------------

inline int&
_begin2(HirschbergSet_ & me) {
    return me.y1;
}

inline int const&
_begin2(HirschbergSet_ const & me) {
    return me.y1;
}

// ----------------------------------------------------------------------------
// Function _setBegin2()
// ----------------------------------------------------------------------------

inline void
_setBegin2(HirschbergSet_ & me, int const & new_begin) {
    me.y1 = new_begin;
}

// ----------------------------------------------------------------------------
// Function _end2()
// ----------------------------------------------------------------------------

inline int&
_end2(HirschbergSet_ & me) {
    return me.y2;
}

inline int const&
_end2(HirschbergSet_ const & me) {
    return me.y2;
}

// ----------------------------------------------------------------------------
// Function _setEnd2()
// ----------------------------------------------------------------------------

inline void
_setEnd2(HirschbergSet_ & me, int const & new_end) {
    me.y2 = new_end;
}

// ----------------------------------------------------------------------------
// Function _score()
// ----------------------------------------------------------------------------

inline int&
_score(HirschbergSet_ & me)    {
    return me.score;
}

// ----------------------------------------------------------------------------
// Function _score()
// ----------------------------------------------------------------------------

inline int const&
_score(HirschbergSet_ const & me)
{
    return me.score;
}

// ----------------------------------------------------------------------------
// Function _setScore()
// ----------------------------------------------------------------------------

inline void
_setScore(HirschbergSet_ & me,int new_score) {
    me.score = new_score;
}

// ----------------------------------------------------------------------------
// Function _setBegin2()
// ----------------------------------------------------------------------------

// //////////////////////////////////////////////////////////////////////////////////////////////
//  Debug Methods
//        functions are only used for debugging or verbose output, therefore they
//      are only active in SEQAN_DEBUG
// //////////////////////////////////////////////////////////////////////////////////////////////
#ifdef SEQAN_DEBUG

inline
void
print(HirschbergSet_ const & me)
{
    std::cout << me.x1 << " " << me.x2 << "\t" << me.y1 << " " << me.y2 << std::endl;
}
#endif

// ----------------------------------------------------------------------------
// Function operator==()
// ----------------------------------------------------------------------------

inline bool
operator==(HirschbergSet_ const & lhs,
           HirschbergSet_ const & rhs)
{
    return ((_begin1(lhs) == _begin1(rhs)) && (_end1(lhs) == _end1(rhs)) &&
            (_begin2(lhs) == _begin2(rhs)) && (_end2(lhs) == _end2(rhs)));
}

// ----------------------------------------------------------------------------
// Function _writeDebugMatrix()
// ----------------------------------------------------------------------------

#ifdef SEQAN_HIRSCHBERG_DEBUG_CUT
    template<typename TSource>
    void _writeDebugMatrix(TSource s1,TSource s2)
    {
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

// debug flag .. define to see where Hirschberg cuts the sequences
//#define SEQAN_HIRSCHBERG_DEBUG_CUT

// ----------------------------------------------------------------------------
// Function globalAlignment()
// ----------------------------------------------------------------------------

template <typename TSequenceH, typename TGapsSpecH, typename TSequenceV, typename TGapsSpecV,
          typename TScoreValue, typename TScoreSpec>
TScoreValue
_globalAlignment(Gaps<TSequenceH, TGapsSpecH> & gapsH,
                 Gaps<TSequenceV, TGapsSpecV> & gapsV,
                 Score<TScoreValue, TScoreSpec> const & score_,
                 Hirschberg const & /*algorithmTag*/)
{
    TSequenceH const & s1 = source(gapsH);
    TSequenceV const & s2 = source(gapsV);

    TScoreValue total_score = 0;

    typedef typename Value<TSequenceV>::Type TValueV;

    typedef typename Size<TSequenceH>::Type TStringSize;

    typedef typename Iterator<TSequenceH const, Standard>::Type TSequenceHIter;
    typedef typename Iterator<TSequenceV const, Standard>::Type TSequenceVIter;

    typedef typename Iterator<Gaps<TSequenceH, TGapsSpecH> >::Type TGapsHIter;
    typedef typename Iterator<Gaps<TSequenceV, TGapsSpecV> >::Type TGapsVIter;

    TGapsHIter target_0 = begin(gapsH);
    TGapsVIter target_1 = begin(gapsV);

    typedef typename Iterator<Matrix<TScoreValue> >::Type TMatrixIterator;

    TValueV v;

    TStringSize len1 = length(s1);
    TStringSize len2 = length(s2);

    // string to store the score values for the currently active cell
    String<TScoreValue> c_score;
    resize(c_score,len2 + 1);
    // string to strore the backpointers
    String<int> pointer;
    resize(pointer,len2 + 1);

    // scoring-scheme specific score values
    TScoreValue score_match = scoreMatch(score_);
    TScoreValue score_mismatch = scoreMismatch(score_);
    TScoreValue score_gap = scoreGapExtend(score_);

    TScoreValue border,s,sg,sd,sg1,sg2;
    int dp;

    std::stack<HirschbergSet_> to_process;
    HirschbergSet_ target;

    int i,j;

    HirschbergSet_ hs_complete(0,len1,0,len2,0);
    to_process.push(hs_complete);

    while(!to_process.empty())
    {
        target = to_process.top();
        to_process.pop();

        if(_begin2(target) == _end2(target))
        {
            for(i = 0;i < (_end1(target) - _begin1(target));++i)
            {
                insertGap(target_1);
                ++target_0;
                ++target_1;
            }
        }
        if(_begin1(target) == _end1(target))
        {
            for(i = 0;i < (_end2(target) - _begin2(target));++i)
            {
                insertGap(target_0);
                ++target_0;
                ++target_1;
            }
        }
        else if(_begin1(target) + 1 == _end1(target) || _begin2(target) + 1 == _end2(target))
        {
            /* ALIGN */
#ifdef SEQAN_HIRSCHBERG_DEBUG_CUT
            std::cout << "align s1 " << _begin1(target) << " to " << _end1(target) << " and s2 " << _begin2(target) << " to " << _end2(target) << std::endl;
            std::cout << "align " << infix(s1,_begin1(target),_end1(target)) << " and " << infix(s2,_begin2(target),_end2(target)) << std::endl << std::endl;
#endif

            TStringSize len_1 = _end1(target) - _begin1(target);
            TStringSize len_2 = _end2(target) - _begin2(target);

            Matrix<TScoreValue> matrix_;

            setDimension(matrix_, 2);
            setLength(matrix_, 0, len_1 + 1);
            setLength(matrix_, 1, len_2 + 1);
            resize(matrix_);

            /* init matrix */
            TSequenceHIter x_begin = iter(s1, _begin1(target), Standard()) - 1;
            TSequenceHIter x_end = iter(s1, _end1(target), Standard()) - 1;
            TSequenceVIter y_begin = iter(s2, _begin2(target), Standard()) - 1;
            TSequenceVIter y_end = iter(s2, _end2(target), Standard()) - 1;

            TSequenceHIter x = x_end;
            TSequenceVIter y;

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
            for (x = x_end; x != x_begin; --x)
            {
                goPrevious(finger1, 0);
                *finger1 = border_;
                border_ += score_gap;
            }

            //-------------------------------------------------------------------------
            //fill matrix
            border_ = 0;
            for (y = y_end; y != y_begin; --y)
            {
                TValueV cy = *y;
                h = border_;
                border_ += score_gap;
                v = border_;

                finger2 = col_;
                goPrevious(col_, 1);
                finger1 = col_;

                *finger1 = v;

                for (x = x_end; x != x_begin; --x)
                {
                    goPrevious(finger1, 0);
                    goPrevious(finger2, 0);
                    if (*x == cy)
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
            total_score += value(matrix_, 0,0);
#ifdef SEQAN_HIRSCHBERG_DEBUG_CUT
            std::cout << "alignment score is " << total_score << std::endl << std::endl;
#endif

            /* TRACE BACK */
            finger1 = begin(matrix_);
            x = iter(s1,_begin1(target));
            y = iter(s2,_begin2(target));
            x_end = iter(s1,_end1(target));
            y_end = iter(s2,_end2(target));

            while ((x != x_end) && (y != y_end))
            {
                bool gv;
                bool gh;

                if (*x == *y)
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
                    ++x;
                    goNext(finger1, 0);
                }
                else
                {
                    insertGap(target_0);
                }

                if (gh)
                {
                    ++y;
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
            while(x != x_end)
            {
                insertGap(target_1);
                ++target_0;
                ++target_1;
                ++x;
            }

            while(y != y_end)
            {
                insertGap(target_0);
                ++target_0;
                ++target_1;
                ++y;
            }
            /* END ALIGN */
        }
        else
        {
            /*
                Calculate cut using the algorithm as proposed in the lecture of Clemens Gröpl
                using a backpointer to remember the position where the optimal alignment passes
                the mid column
            */
            int mid = static_cast<int>(floor( static_cast<double>((_begin1(target) + _end1(target))/2) ));

#ifdef SEQAN_HIRSCHBERG_DEBUG_CUT
            std::cout << "calculate cut for s1 " << _begin1(target) << " to " << _end1(target) << " and s2 " << _begin2(target) << " to " << _end2(target) << std::endl;
            std::cout << "calculate cut for " << infix(s1,_begin1(target),_end1(target)) << " and " << infix(s2,_begin2(target),_end2(target)) << std::endl;
            std::cout << "cut is in row " << mid << " symbol is " << getValue(s1,mid-1) << std::endl << std::endl;


            _writeDebugMatrix(infix(s1,_begin1(target),_end1(target)),infix(s2,_begin2(target),_end2(target)));
#endif

            border = 0;
            for(i = _begin2(target);i <= _end2(target);++i)
            {
                c_score[i] = border;
                border += score_gap;
                pointer[i] = i;
            }

            // iterate over s1 until the mid column is reached
            border = score_gap;
            for(i = _begin1(target) + 1;i <= mid;++i)
            {
                s = c_score[_begin2(target)];
                c_score[_begin2(target)] = border;
                border += score_gap;
                v = getValue(s1,i-1);
                for(j = _begin2(target) + 1;j <= _end2(target);++j)
                {
                    sg = score_gap + ((c_score[j] > c_score[j - 1]) ? c_score[j] : c_score[j - 1]);
                    sd = s + ((v == getValue(s2,j-1)) ? score_match : score_mismatch);

                    s = c_score[j];
                    c_score[j] = (sg > sd) ? sg : sd;
                }
            }

            // from here, rememeber the cell of mid-column, where optimal alignment passed
            for(i = mid + 1;i <= _end1(target);++i)
            {
                s = c_score[_begin2(target)];
                c_score[_begin2(target)] = border;
                border += score_gap;
                v = getValue(s1,i-1);

                dp = _begin2(target);

                for(j = _begin2(target) + 1;j <= _end2(target);++j)
                {
                    sg1 = score_gap + c_score[j];
                    sg2 = score_gap + c_score[j - 1];

                    sd = s + ((v == getValue(s2,j-1)) ? score_match : score_mismatch);

                    s = c_score[j];
                    sg = pointer[j];
                    if(sd >= _max(sg1,sg2))
                    {
                        c_score[j] = sd;
                        pointer[j] = dp;
                    }
                    else
                    {
                        if(sg2 > sg1)
                        {
                            c_score[j] = sg2;
                            pointer[j] = pointer[j-1];
                        }
                        else
                        {
                            // gap introduced from left
                            // no update for the pointer
                            c_score[j] = sg1;
                        }
                    }
                    dp = sg;
                }
            }

#ifdef SEQAN_HIRSCHBERG_DEBUG_CUT
            std::cout << "hirschberg calculates cut in column " << mid << " and row " << pointer[_end2(target)] << std::endl;
            std::cout << "requested position in c_score and pointer is " << _end2(target) << std::endl;
            std::cout << "alignment score is " << c_score[_end2(target)] << std::endl << std::endl;
#endif
            to_process.push(HirschbergSet_(mid,_end1(target),pointer[_end2(target)],_end2(target),0));
            to_process.push(HirschbergSet_(_begin1(target),mid,_begin2(target),pointer[_end2(target)],0));
        }
        /* END CUT */
    }
    return total_score;
}

}  // namespace seqan

#endif  // #ifndef SEQAN_INCLUDE_SEQAN_ALIGN_GLOBAL_ALIGNMENT_HIRSCHBERG_IMPL_H_
