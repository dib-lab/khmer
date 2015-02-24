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
// Author: Andreas Gogol-Doering <andreas.doering@mdc-berlin.de>
// Author: Anne-Katrin Emde <anne-katrin.emde@fu-berlin.de>
// ==========================================================================
// Implementation of the Waterman-Eggert algorithm, sometimes also called
// Smith-Waterman algorithm with declumping.
// ==========================================================================

#ifndef SEQAN_INCLUDE_SEQAN_ALIGN_LOCAL_ALIGNMENT_WATERMAN_EGGERT_IMPL_H_
#define SEQAN_INCLUDE_SEQAN_ALIGN_LOCAL_ALIGNMENT_WATERMAN_EGGERT_IMPL_H_

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// ----------------------------------------------------------------------------
// Helper Class ScoreAndID
// ----------------------------------------------------------------------------

// TODO(holtgrew): Why is Pair<> not enough?
// TODO(holtgrew): Rename to ScoreAndID_ (with trailing underscore)

// Simple class that stores a value with an ID.
template <typename TValue, typename TID>
class ScoreAndID
{
public:
    TValue value_;
    TID id_;

    ScoreAndID() : value_(MinValue<TValue>::VALUE), id_(MaxValue<TValue>::VALUE)
    {}

    ScoreAndID(TValue score, TID id_pos)
    {
        value_ = score;
        id_ = id_pos;
    }

    inline bool operator>(ScoreAndID const & other) const
    {
        return value_ > other.value_;
    }

    inline bool operator<(ScoreAndID const & other) const
    {
        return value_ < other.value_;
    }
};

// ----------------------------------------------------------------------------
// Class LocalAlignmentFinder
// ----------------------------------------------------------------------------

template <typename TScoreValue = int>
class LocalAlignmentFinder
{
public:
    typedef Matrix<TScoreValue> TMatrix;
    typedef typename Position<TMatrix>::Type TMatrixPosition;
    typedef typename Size<TMatrix>::Type TSize;
    typedef ScoreAndID<TScoreValue,TMatrixPosition> TPQEntry;

    typedef Iter<TMatrix,PositionIterator> TMatrixIterator;
    typedef PriorityType<TPQEntry> TPriorityQ;
    typedef String<bool> TBoolMatrix;

    //DP-matrix
    TMatrix matrix;
    //matrix that memorizes the cells from which not to go diagonal
    TBoolMatrix forbidden;
    //priority queue for quickly finding the maximum score in the DP-matrix
    TPriorityQ pQ;
    //position of maximum score (where traceback is started from)
    TMatrixPosition bestEndPos;
    //position where traceback ended and where declumping begins
    TMatrixPosition bestBeginPos;
    //traceback path that is set to forbidden while declumping
    AlignTraceback<TSize> trace;

    bool needReinit; //true: call "smithWaterman", false: call "smithWatermanGetNext"

    LocalAlignmentFinder() : bestEndPos(0), bestBeginPos(0), needReinit(true)
    {}

    // TODO(holtgrew): Remove and replace all occurrences with default constructor.
    template<typename TAlign>
    LocalAlignmentFinder(TAlign const &) : bestEndPos(0), bestBeginPos(0), needReinit(true)
    {}
};

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function _initLocalAlignmentFinder()
// ----------------------------------------------------------------------------

template<typename TSequenceH, typename TSequenceV, typename TScoreValue, typename TTag>
void
_initLocalAlignmentFinder(TSequenceH const & seqH,
                          TSequenceV const & seqV,
                          LocalAlignmentFinder<TScoreValue> & finder,
                          TTag) {
    typedef LocalAlignmentFinder<TScoreValue> TFinder;
    typedef typename TFinder::TMatrix TMatrix;
    typedef typename Size<TMatrix>::Type TSize;

    TSize len0 = length(seqH);
    TSize len1 = length(seqV);

    setDimension(finder.matrix, 2);
    setLength(finder.matrix, 0, len0 + 1);
    setLength(finder.matrix, 1, len1 + 1);
    resize(finder.matrix);

    resize(finder.forbidden, (len0 + 1) * (len1 + 1), false);

    finder.bestEndPos = maxValue<typename TFinder::TMatrixPosition>();
    finder.bestBeginPos = maxValue<typename TFinder::TMatrixPosition>();
}

// ----------------------------------------------------------------------------
// Function clear()
// ----------------------------------------------------------------------------

template <typename TScoreValue>
void clear(LocalAlignmentFinder<TScoreValue> & sw_finder)
{
    sw_finder.needReinit = true;
}

// ----------------------------------------------------------------------------
// Function getScore()
// ----------------------------------------------------------------------------

template <typename TScoreValue>
TScoreValue getScore(LocalAlignmentFinder<TScoreValue> const & sw)
{
    typedef LocalAlignmentFinder<TScoreValue> TFinder;
    if(sw.bestEndPos !=  maxValue<typename TFinder::TMatrixPosition>())
        return getValue(const_cast<typename TFinder::TMatrix &>(sw.matrix), sw.bestEndPos);
    return 0;
}

// ----------------------------------------------------------------------------
// Function _smithWatermanGetMatrix()
// ----------------------------------------------------------------------------

template <typename TScoreValue, typename TStringH, typename TStringV>
TScoreValue
_smithWatermanGetMatrix(LocalAlignmentFinder<TScoreValue> & sw,
                        TStringH const & strH,
                        TStringV const & strV,
                        Score<TScoreValue, Simple> const & score_,
                        TScoreValue cutoff)
{
    // typedefs
    typedef Matrix<TScoreValue> TMatrix;
    typedef typename Position<TMatrix>::Type TMatrixPosition;
    typedef Iter<TMatrix,PositionIterator> TMatrixIterator;

    typedef typename Iterator<TStringH const, Rooted>::Type TStringIteratorH;
    //typedef typename Value<TStringH const>::Type TValueH;
    typedef typename Iterator<TStringV const, Rooted>::Type TStringIteratorV;
    typedef typename Value<TStringV const>::Type TValueV;

    //-------------------------------------------------------------------------
    //define some variables


//    TSize str1_length = length(strH);
//    TSize str2_length = length(strV);
    TStringIteratorH x_begin = begin(strH) - 1;
    TStringIteratorH x_end = end(strH) - 1;
    TStringIteratorV y_begin = begin(strV) - 1;
    TStringIteratorV y_end = end(strV) - 1;

    TStringIteratorH x = x_end;
    TStringIteratorV y;

    TScoreValue score_match = scoreMatch(score_);
    TScoreValue score_mismatch = scoreMismatch(score_);
    TScoreValue score_gap = scoreGapExtend(score_);

    TScoreValue h = 0;
    TScoreValue v = 0;

    TMatrixIterator col_ = end(sw.matrix) - 1;
    TMatrixIterator finger1;
    TMatrixIterator finger2;

    //-------------------------------------------------------------------------
    // init

    finger1 = col_;
    *finger1 = 0;
    //std::cout <<"  ";
    for (x = x_end; x != x_begin; --x)
    {
        goPrevious(finger1, 0);
        *finger1 = 0;
    }

    //-------------------------------------------------------------------------
    //fill matrix
    for (y = y_end; y != y_begin; --y)
    {
        TValueV cy = *y;

        h = 0;
        v = 0;

        finger2 = col_;        //points to last column
        goPrevious(col_, 1);    //points to this column
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
                if (v < 0) v = 0;

            }
            *finger1 = v;
            if (v >= cutoff)
            {
                push(sw.pQ,ScoreAndID<TScoreValue,TMatrixPosition>(v,position(finger1)));
            }
        }
    }

    // check if any scores >= cutoff were found
    if(!empty(sw.pQ))
    {
        ScoreAndID<TScoreValue,TMatrixPosition> best = top(sw.pQ);
        v = getValue(sw.matrix,best.id_);
        sw.bestEndPos = best.id_;
    }
    else
        v=0;

    return v;
}

// ----------------------------------------------------------------------------
// Function _smithWatermanDeclump()
// ----------------------------------------------------------------------------

// declumping
template <typename TScoreValue, typename TSequenceH, typename TGapsSpecH, typename TSequenceV, typename TGapsSpecV>
void
_smithWatermanDeclump(LocalAlignmentFinder<TScoreValue> & sw ,
                      Gaps<TSequenceH, TGapsSpecH> & gapsH,
                      Gaps<TSequenceV, TGapsSpecV> & gapsV,
                      Score<TScoreValue, Simple> const & score_)
{
//-------------------------------------------------------------------------
//typedefs
    //typedef typename LocalAlignmentFinder<TScoreValue>::TMatrixPosition TMatrixPosition;
    typedef typename LocalAlignmentFinder<TScoreValue>::TMatrix         TMatrix;
    typedef Iter<TMatrix, PositionIterator>                             TMatrixIterator;

    typedef Gaps<TSequenceH, TGapsSpecH>                TGapsH;
    typedef typename Iterator<TGapsH, Rooted>::Type     TGapsHIter;
    typedef typename Iterator<TSequenceH, Rooted>::Type TSequenceHIter;
    //typedef typename Value<TSequenceH>::Type            TValueH;

    typedef Gaps<TSequenceV, TGapsSpecV>                TGapsV;
    typedef typename Iterator<TGapsV, Rooted>::Type     TGapsVIter;
    typedef typename Iterator<TSequenceV, Rooted>::Type TSequenceVIter;
    typedef typename Value<TSequenceV>::Type            TValueV;

//-------------------------------------------------------------------------
//variables
    // TRow row0 = row(align_,0);
    // TRow row1 = row(align_,1);

    // beginPosition == # leading gaps
    // endPosition == length of clipped region without trailing gaps
    // clippedEndPosition == source position of clipping end.

    // TAlignIterator ali_it0_stop = iter(row0,beginPosition(row0));
    // TAlignIterator ali_it1_stop = iter(row1,beginPosition(row1));
    TGapsHIter ali_it0_stop = begin(gapsH);
    TGapsVIter ali_it1_stop = begin(gapsV);

    // SEQAN_ASSERT( endPosition(row0)- beginPosition(row0) == endPosition(row1)- beginPosition(row1) );

    // TAlignIterator ali_it0 = iter(row0,endPosition(row0));
    // TAlignIterator ali_it1 = iter(row1,endPosition(row1));
    TGapsHIter ali_it0 = end(gapsH);
    TGapsVIter ali_it1 = end(gapsV);

    // TStringIterator x_begin = begin(source(row0))-1;
    // TStringIterator y_begin = begin(source(row1))-1;
    // TStringIterator x_end = iter(source(row0),clippedEndPosition(row0))-1;
    // TStringIterator y_end = iter(source(row1),clippedEndPosition(row1))-1;
    TSequenceHIter x_begin = begin(source(gapsH))-1;
    TSequenceVIter y_begin = begin(source(gapsV))-1;
    TSequenceHIter x_end = iter(source(gapsH), endPosition(gapsH) - 1);
    TSequenceVIter y_end = iter(source(gapsV), endPosition(gapsV) - 1);

    // TStringIterator x = x_end;
    // TStringIterator y = y_end;
    // TStringIterator x_stop = x_end;
    TSequenceHIter x = x_end;
    TSequenceVIter y = y_end;
    TSequenceHIter x_stop = x_end;


    TScoreValue score_match = scoreMatch(score_);
    TScoreValue score_mismatch = scoreMismatch(score_);
    TScoreValue score_gap = scoreGapExtend(score_);
    TScoreValue h,v;

    TMatrixIterator finger0 = iter(sw.matrix,sw.bestBeginPos);
    TMatrixIterator end_col = finger0;
    TMatrixIterator finger1 = finger0;
    TMatrixIterator forbidden = finger0;

    bool different = true;
    bool forbidden_reached = true;
    bool end_reached = false;
    bool skip_row = false;


/*    int str0_length = length(source(row(align_,0)))+1;
    int str1_length = length(source(row(align_,1)))+1;
    for(int i = 0; i <str1_length; ++i){
         for(int j=0;j<str0_length;++j){
             std::cout << getValue(sw.matrix,(str0_length*i)+j);
             if(sw.forbidden[(str0_length*i)+j]==true)
                 std::cout <<"(1) ";
             else
                 std::cout <<"(0) ";
         }
         std::cout <<"\n";
     }*/

    clearClipping(gapsH);
    clearClipping(gapsV);
    // setClippedBeginPosition(row(align_, 0),0);
    // setClippedBeginPosition(row(align_, 1),0);


    for (y = y_end; (y != y_begin) ; --y)
    {
        different = true;
        //compute next "forbidden" cell (where you are not allowed to go diagonal)
        if(forbidden_reached && !end_reached)
        {

            if(ali_it0==ali_it0_stop && ali_it1==ali_it1_stop)
            {
                end_reached = true;
            }

            if(!end_reached)
            {
                --ali_it0;
                goPrevious(forbidden,1);
                while(isGap(ali_it0)&& ali_it0!=ali_it0_stop)
                {
                    skip_row = true;
                    --ali_it0;
                    --ali_it1;
                    goPrevious(forbidden,1);
                }

                --ali_it1;
                goPrevious(forbidden,0);
                while(isGap(ali_it1)&& ali_it1!=ali_it1_stop)
                {
                    --ali_it1;
                    --ali_it0;
                    goPrevious(forbidden,0);
                }
                // mark the forbidden cell
                sw.forbidden[position(forbidden)]=true;
                forbidden_reached = false;
            }

        }

        TValueV cy = *y;

        h = *end_col;

        finger1 = end_col;        //points to last column
        goPrevious(end_col, 1);    //points to this column
        finger0 = end_col;

        v = *finger0;
        x = x_end;

        // declump the current row until the cells in the declumped and
        // the cells in the original matrix are the same (or the border is reached)
        // indicated by bool different
        while (x != x_begin && different)
        {
            goPrevious(finger0, 0);
            goPrevious(finger1, 0);
            if (*x == cy && !(sw.forbidden[position(finger0)]))
            {
                v = h + score_match;
                h = *finger1;
            }
            else
            {
                TScoreValue s1;

                if(finger0 == forbidden)
                {
                        skip_row = false;
                        forbidden_reached = true;
                        s1 = 0;
                }
                else
                {
                    if(sw.forbidden[position(finger0)]) s1 = 0;
                    else s1 = h + score_mismatch;
                }

                h = *finger1;
                TScoreValue s2 = score_gap + ((h > v) ? h : v);
                v = (s1 > s2) ? s1 : s2;
                if (v < 0) v = 0;

            }

            // value is the same as in the original matrix
            if(*finger0==v)
            {
                //x_stop is as far as we have to go at least
                if(x<x_stop)
                {
                    different=false;
            //        x_stop=x;
                }
                else
                {
                    // x_end indicates where to start in the next row
                    if(x==x_end && ((!forbidden_reached && !skip_row)||end_reached))
                    {
                        --x_end;
                        goPrevious(end_col, 0);
                    }
                }
            }
            if(x<x_stop)// && different)
            {
                x_stop=x;
            }
            *finger0 = v;
            --x;
        }
        if(x_end==x_begin)
            break;
    }

//    cout <<"...declumped.\n";
}

// ----------------------------------------------------------------------------
// Function _smithWatermanTrace()
// ----------------------------------------------------------------------------

// Traceback.
template <typename TSourceH, typename TGapsSpecH, typename TSourceV, typename TGapsSpecV, typename TScoreValue, unsigned DIMENSION>
typename Iterator<Matrix<TScoreValue, DIMENSION>, Standard >::Type
_smithWatermanTrace(Gaps<TSourceH, TGapsSpecH> & gapsH,
                    Gaps<TSourceV, TGapsSpecV> & gapsV,
                    typename LocalAlignmentFinder<TScoreValue>::TBoolMatrix & fb_matrix,
                    Iter< Matrix<TScoreValue, DIMENSION>, PositionIterator > source_,
                    Score<TScoreValue, Simple> const & scoring_) {
    //typedefs
    typedef Iter<Matrix<TScoreValue, DIMENSION>, PositionIterator > TMatrixIterator;
    typedef typename Position<Matrix<TScoreValue, DIMENSION> >::Type TPosition;

//    typedef Segment<TTargetSource, InfixSegment> TTargetSourceSegment;
    typedef typename Iterator<TSourceH, Standard>::Type TSourceIteratorH;
    typedef typename Iterator<TSourceV, Standard>::Type TSourceIteratorV;

    typedef Gaps<TSourceH, TGapsSpecH> TGapsH;
    typedef Gaps<TSourceV, TGapsSpecV> TGapsV;
    typedef typename Iterator<TGapsH, Rooted>::Type TTargetIteratorH;
    typedef typename Iterator<TGapsV, Rooted>::Type TTargetIteratorV;

    //-------------------------------------------------------------------------
    //variables
    TPosition pos_0 = coordinate(source_, 0);
    TPosition pos_1 = coordinate(source_, 1);

    TSourceH strH = source(gapsH);
    TSourceV strV = source(gapsV);

    TTargetIteratorH target_0 = iter(gapsH, pos_0);
    TTargetIteratorV target_1 = iter(gapsV, pos_1);

    TSourceIteratorH it_0 = iter(strH, pos_0, Standard());
    TSourceIteratorH it_0_end = end(strH);

    TSourceIteratorV it_1 = iter(strV, pos_1, Standard());
    TSourceIteratorV it_1_end = end(strV);

    TScoreValue score_mismatch = scoreMismatch(scoring_);
    TScoreValue score_gap = scoreGapExtend(scoring_);

    //-------------------------------------------------------------------------
    //follow the trace until 0 is reached
    while ((*source_!=0) && (it_0 != it_0_end) && (it_1 != it_1_end))
    {
        bool gv;
        bool gh;
        bool forbidden = fb_matrix[position(source_)];

        if (*it_0 == *it_1 && !forbidden)
        {
            gv = gh = true;
        }
        else
        {
            TMatrixIterator it_ = source_;

            goNext(it_, 0);
            TScoreValue v = *it_ + score_gap;

            TScoreValue d;
            if(forbidden)
                d = 0;
            else{
                goNext(it_, 1);
                d = *it_ + score_mismatch;
            }

            it_ = source_;
            goNext(it_, 1);
            TScoreValue h = *it_ + score_gap;

            gv = (v >= h) | (d >= h);
            gh = (h >  v) | (d >= v);
        }

        if (gv)
        {
            ++it_0;
            goNext(source_, 0);
        }
        else
        {
            insertGap(target_0);
        }

        if (gh)
        {
            ++it_1;
            goNext(source_, 1);
        }
        else
        {
            insertGap(target_1);
        }
        ++target_0;
        ++target_1;
    }

    // We have removed all gaps and clippings from gapsH and gapsV in the calling functions, so the following works.
    // Note that we have to set the end position first.
    // TODO(holtgrew): Use setBegin/EndPosition().
    setClippedEndPosition(gapsH, toViewPosition(gapsH, position(it_0, strH)));
    setClippedEndPosition(gapsV, toViewPosition(gapsV, position(it_1, strV)));
    setClippedBeginPosition(gapsH, toViewPosition(gapsH, pos_0));
    setClippedBeginPosition(gapsV, toViewPosition(gapsV, pos_1));

    return source_;
}

// ----------------------------------------------------------------------------
// Function _getNextBestEndPosition()
// ----------------------------------------------------------------------------

// Adjust the priority queue of scores until the true maximum is found.
template <typename TScoreValue>
typename LocalAlignmentFinder<TScoreValue>::TMatrixPosition
_getNextBestEndPosition(LocalAlignmentFinder<TScoreValue> & sw ,
                        TScoreValue cutoff)
{
    // get maximal score from priority queue
    TScoreValue topScore = 0;
    if (!empty(sw.pQ))
        topScore = getValue(sw.matrix, top(sw.pQ).id_);

    // check if matrix entry of topScore did not change while declumping
    if (!empty(sw.pQ)) {
        while (top(sw.pQ).value_ != topScore) {
            if (topScore >= cutoff) {
                ((sw.pQ).heap[0]).value_ = topScore;
                adjustTop(sw.pQ);
            } else {
                pop(sw.pQ);
            }
            if (!empty(sw.pQ)) topScore = getValue(sw.matrix, top(sw.pQ).id_);
            else break;
        }
    }

    // priority queue with top scores is empty
    if(empty(sw.pQ)) {//||top(sw.pQ).value_<cutoff) {
        sw.needReinit = true;
        return 0;
    }

    typename LocalAlignmentFinder<TScoreValue>::TMatrixPosition ret_pos = top(sw.pQ).id_;
    sw.bestEndPos = ret_pos;
    pop(sw.pQ);

    return ret_pos;
}

// ----------------------------------------------------------------------------
// Function _smithWaterman()
// ----------------------------------------------------------------------------

// Wrapper that computes the matrix and does the backtracking for the best alignment
template <typename TSourceH, typename TGapsSpecH, typename TSourceV, typename TGapsSpecV, typename TScoreValue, typename TScoreSpec>
TScoreValue
_smithWaterman(Gaps<TSourceH, TGapsSpecH> & gapsH,
               Gaps<TSourceV, TGapsSpecV> & gapsV,
               LocalAlignmentFinder<TScoreValue> & sw_finder,
               Score<TScoreValue, TScoreSpec> const & score_,
               TScoreValue cutoff)
{
    // TODO(holtgrew): This sourceSegment() stuff is confusing... Do we *really* need this?
    // Clear gaps and clipping from result Gaps structures.
    clearGaps(gapsH);
    clearGaps(gapsV);
    clearClipping(gapsH);
    clearClipping(gapsV);

    _initLocalAlignmentFinder(sourceSegment(gapsH), sourceSegment(gapsV), sw_finder, WatermanEggert());

    TScoreValue ret = _smithWatermanGetMatrix(sw_finder, sourceSegment(gapsH), sourceSegment(gapsV), score_,cutoff);

    if(ret==0)
        return ret;
    sw_finder.needReinit = false;

    typedef Iter<typename LocalAlignmentFinder<TScoreValue>::TMatrix,PositionIterator > TMatrixIterator;
    TMatrixIterator best_begin;

    // TODO(holtgrew): What does the following comment mean?
    // TODO: sw_finder statt kram
    best_begin = _smithWatermanTrace(gapsH, gapsV, sw_finder.forbidden,iter(sw_finder.matrix,(top(sw_finder.pQ)).id_), score_);

    sw_finder.bestBeginPos = position(best_begin);

    pop(sw_finder.pQ);

    return ret;
}

// ----------------------------------------------------------------------------
// Function _smithWatermanGetNext()
// ----------------------------------------------------------------------------

// Wrapper that declumps the matrix and traces back the next best alignment
template <typename TSequenceH, typename TGapsSpecH, typename TSequenceV, typename TGapsSpecV, typename TScoreValue, typename TScoreSpec>
TScoreValue
_smithWatermanGetNext(Gaps<TSequenceH, TGapsSpecH> & gapsH,
                      Gaps<TSequenceV, TGapsSpecV> & gapsV,
                      LocalAlignmentFinder<TScoreValue> & sw_finder ,
                      Score<TScoreValue, TScoreSpec> const & score_,
                      TScoreValue cutoff)
{
    _smithWatermanDeclump(sw_finder, gapsH, gapsV, score_);

    clearGaps(gapsH);
    clearGaps(gapsV);
    clearClipping(gapsH);
    clearClipping(gapsV);

    typename LocalAlignmentFinder<TScoreValue>::TMatrixPosition next_best_end;
    next_best_end = _getNextBestEndPosition(sw_finder,cutoff);
    if(next_best_end==0)
        return 0;
    typename LocalAlignmentFinder<TScoreValue>::TMatrixIterator next_best_begin;
    next_best_begin= _smithWatermanTrace(gapsH, gapsV, sw_finder.forbidden,iter(sw_finder.matrix,next_best_end), score_);
    sw_finder.bestBeginPos = position(next_best_begin);

    return getValue(sw_finder.matrix,next_best_end);
}

}  // namespace seqan

#endif  // #ifndef SEQAN_INCLUDE_SEQAN_ALIGN_LOCAL_ALIGNMENT_WATERMAN_EGGERT_IMPL_H_
