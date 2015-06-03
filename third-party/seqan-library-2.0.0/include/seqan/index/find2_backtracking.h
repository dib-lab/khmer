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
// Author: Enrico Siragusa <enrico.siragusa@fu-berlin.de>
// ==========================================================================
// Approximate string matching via backtracking on two substring indices.
// ==========================================================================

#ifndef SEQAN_FIND_BACKTRACKING_MULTIPLE_H_
#define SEQAN_FIND_BACKTRACKING_MULTIPLE_H_

//#define SEQAN_DEBUG

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

template <typename TDistance, typename TSpec>
struct Backtracking;

struct StageInitial_;
struct StageUpper_;
struct StageDiagonal_;
struct StageLower_;
struct StageFinal_;
struct StageExact_;

// ============================================================================
// Metafunctions
// ============================================================================

// ----------------------------------------------------------------------------
// Metafunction TextIterator_
// ----------------------------------------------------------------------------

template <typename TText, typename TTextIndexSpec, typename TPattern, typename TPatternIndexSpec, typename TDistance, typename TSpec>
struct TextIterator_<Index<TText, TTextIndexSpec>, Index<TPattern, TPatternIndexSpec>, Backtracking<TDistance, TSpec> >
{
    typedef typename Iterator<Index<TText, TTextIndexSpec>, TopDown<> >::Type   Type;
};

// ----------------------------------------------------------------------------
// Metafunction PatternIterator_
// ----------------------------------------------------------------------------

template <typename TText, typename TTextIndexSpec, typename TPattern, typename TPatternIndexSpec, typename TDistance, typename TSpec>
struct PatternIterator_<Index<TText, TTextIndexSpec>, Index<TPattern, TPatternIndexSpec>, Backtracking<TDistance, TSpec> >
{
    typedef typename Iterator<Index<TPattern, TPatternIndexSpec>, TopDown<> >::Type    Type;
};

// ----------------------------------------------------------------------------
// Metafunction VertexScore_
// ----------------------------------------------------------------------------

template <typename TBacktracking>
struct VertexScore_ {};

template <typename TSpec>
struct VertexScore_<Backtracking<HammingDistance, TSpec> >
{
    typedef unsigned char           Type;
};

template <typename TSpec>
struct VertexScore_<Backtracking<EditDistance, TSpec> >
{
    typedef String<unsigned char>                   TString;
    typedef Segment<TString, InfixSegment>          Type;
    typedef Segment<const TString, InfixSegment>    ConstType;
};

// ----------------------------------------------------------------------------
// Metafunction VertexScoreStack_
// ----------------------------------------------------------------------------

template <typename TBacktracking>
struct VertexScoreStack_ {};

template <typename TSpec>
struct VertexScoreStack_<Backtracking<HammingDistance, TSpec> >
{
    typedef String<typename VertexScore_<Backtracking<HammingDistance, TSpec> >::Type> Type;
};

template <typename TSpec>
struct VertexScoreStack_<Backtracking<EditDistance, TSpec> >
{
    typedef StringSet<typename VertexScore_<Backtracking<EditDistance, TSpec> >::TString, Owner<ConcatDirect<> > > Type;
};

// ----------------------------------------------------------------------------
// Metafunction NextStage_
// ----------------------------------------------------------------------------

template <typename TBacktracking, typename TStage>
struct NextStage_
{
    typedef Nothing         Type;
};

// HammingDistance: StageInitial|StageExact -> StageFinal

template <typename TSpec>
struct NextStage_<Backtracking<HammingDistance, TSpec>, StageInitial_>
{
    typedef StageFinal_     Type;
};

template <typename TSpec>
struct NextStage_<Backtracking<HammingDistance, TSpec>, StageExact_>
{
    typedef StageFinal_     Type;
};

// EditDistance: StageInitial -> StageUpper -> StageDiagonal -> StageLower -> StageFinal

template <typename TSpec>
struct NextStage_<Backtracking<EditDistance, TSpec>, StageInitial_>
{
    typedef StageUpper_     Type;
};

template <typename TSpec>
struct NextStage_<Backtracking<EditDistance, TSpec>, StageUpper_>
{
    typedef StageDiagonal_  Type;
};

template <typename TSpec>
struct NextStage_<Backtracking<EditDistance, TSpec>, StageDiagonal_>
{
    typedef StageLower_     Type;
};

template <typename TSpec>
struct NextStage_<Backtracking<EditDistance, TSpec>, StageLower_>
{
    typedef StageFinal_     Type;
};

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// ----------------------------------------------------------------------------
// Tags for backtracking stages
// ----------------------------------------------------------------------------

struct StageInitial_ {};
struct StageUpper_ {};
struct StageDiagonal_ {};
struct StageLower_ {};
struct StageFinal_ {};
struct StageExact_ {};

// ----------------------------------------------------------------------------
// Class Finder_
// ----------------------------------------------------------------------------

template <typename TText, typename TTextIndexSpec, typename TPattern, typename TPatternIndexSpec,
          typename TDistance, typename TSpec>
struct Finder_<Index<TText, TTextIndexSpec>, Index<TPattern, TPatternIndexSpec>, Backtracking<TDistance, TSpec> >
{
    typedef Index<TText, TTextIndexSpec>                                                TTextIndex;
    typedef Index<TPattern, TPatternIndexSpec>                                          TPatternIndex;
    typedef Backtracking<TDistance, TSpec>                                              TBacktracking;
    typedef typename TextIterator_<TTextIndex, TPatternIndex, TBacktracking>::Type      TTextIterator;
    typedef typename PatternIterator_<TTextIndex, TPatternIndex, TBacktracking>::Type   TPatternIterator;
    typedef String<TTextIterator>                                                       TTextStack;
    typedef String<TPatternIterator>                                                    TPatternStack;
    typedef typename VertexScoreStack_<TBacktracking>::Type                             TVertexScoreStack;
    typedef typename Score_<TBacktracking>::Type                                        TScore;

    TTextStack              textStack;
    TPatternStack           patternStack;
    TVertexScoreStack       scoreStack;
    TScore                  maxScore;

    Finder_() :
        textStack(),
        patternStack(),
        scoreStack(),
        maxScore(0)
    {}
};

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function _min3()
// ----------------------------------------------------------------------------

template <typename TValue>
inline TValue
_min3(TValue a, TValue b, TValue c)
{
    TValue m = a;
    if (m > b) m = b;
    if (m > c) m = c;
    return m;
}

// ----------------------------------------------------------------------------
// Function _updateVertexScore()
// ----------------------------------------------------------------------------

template <typename TVertexScore>
inline void
_updateVertexScore(TVertexScore current,
                   TVertexScore const previous,
                   StageInitial_ const & /* tag */)
{
    // Update last cell.
    // C[i,0] = C[i-1,0] + 1 [Upper]
    back(current) = back(previous) + 1;
}

template <typename TVertexScore, typename TTextValue, typename TPatternIterator, typename TStage>
inline void
_updateVertexScore(TVertexScore current,
                   TVertexScore const previous,
                   TTextValue textChar,
                   TPatternIterator patternIt,
                   TStage const & /* tag */)
{
    typedef typename Iterator<TVertexScore, Standard>::Type         TVertexScoreIterator;
    typedef typename Iterator<TVertexScore const, Standard>::Type   TVertexScoreConstIterator;
    typedef typename Value<TVertexScore>::Type                      TScore;

    TVertexScoreIterator currentIt = begin(current, Standard());
    TVertexScoreIterator columnEnd = end(current, Standard());
    TVertexScoreConstIterator previousIt = begin(previous, Standard());

    // Update first cell.
    if (IsSameType<TStage, StageUpper_>::VALUE)
    {
        // C[0,j] = C[0,j-1] + 1 [Left]
        value(currentIt) = value(previousIt) + 1;
    }
    else
    {
        TScore score = ordEqual(textChar, value(patternIt)) ? 0 : 1;

        // C[i,j] = min { C[i-1,j-1] + d(t,p), C[i-1,j] + 1 } [Diagonal, Left]
        value(currentIt) = _min(value(previousIt) + score, value(previousIt + 1) + 1);

        ++previousIt;
        ++patternIt;
    }

    // Update central cells.
    for (++currentIt; currentIt != columnEnd - 1; ++currentIt, ++previousIt, ++patternIt)
    {
        TScore score = ordEqual(textChar, value(patternIt)) ? 0 : 1;

        // C[i,j] = min { C[i-1,j-1] + d(t,p), C[i-1,j] + 1, C[i,j-1] + 1 } [Diagonal, Left, Upper]
        value(currentIt) = _min3(value(previousIt) + score, value(previousIt + 1) + 1, value(currentIt - 1) + 1);
    }

    // Update last cell.
    if (IsSameType<TStage, StageLower_>::VALUE)
    {
        TScore score = ordEqual(textChar, value(patternIt)) ? 0 : 1;

        // C[i,j] = min { C[i-1,j-1] + d(t,p), C[i-1,j] + 1, C[i,j-1] + 1 } [Diagonal, Left, Upper]
        value(currentIt) = _min3(value(previousIt) + score, value(previousIt + 1) + 1, value(currentIt - 1) + 1);

        ++previousIt;
    }
    else
    {
        TScore score = ordEqual(textChar, value(patternIt)) ? 0 : 1;

        // C[i,j] = min { C[i-1,j-1] + d(t,p), C[i,j-1] + 1 } [Diagonal, Upper]
        value(currentIt) = _min(value(previousIt) + score, value(currentIt - 1) + 1);
    }

    // Assert end of columns.
    SEQAN_ASSERT_EQ(currentIt + 1, end(current, Standard()));
    SEQAN_ASSERT_EQ(previousIt + 1, end(previous, Standard()));
}

template <typename TVertexScore, typename TTextValue, typename TPatternValue>
inline void
_updateVertexScore(TVertexScore current,
                   TVertexScore const previous,
                   TTextValue textChar,
                   TPatternValue patternChar,
                   StageFinal_ const & /* tag */)
{
    typedef typename Iterator<TVertexScore, Standard>::Type         TVertexScoreIterator;
    typedef typename Iterator<TVertexScore const, Standard>::Type   TVertexScoreConstIterator;
    typedef typename Value<TVertexScore>::Type                      TScore;

    TVertexScoreIterator currentIt = begin(current, Standard());
    TVertexScoreConstIterator previousIt = begin(previous, Standard());

    // Update last cell.
    TScore score = ordEqual(textChar, patternChar) ? 0 : 1;

    // C[i,j] = min { C[i-1,j-1] + d(t,p), C[i-1,j] + 1 } [Diagonal, Left]
    value(currentIt) = _min(value(previousIt) + score, value(previousIt + 1) + 1);

    // Assert end of columns.
    SEQAN_ASSERT_EQ(currentIt + 1, end(current, Standard()));
    SEQAN_ASSERT_EQ(previousIt + 2, end(previous, Standard()));
}

// ----------------------------------------------------------------------------
// Function copyBackAndResize()                                     [StringSet]
// ----------------------------------------------------------------------------

template <typename TString, typename TSSetSpec, typename TDelta>
inline void
copyBackAndResize(StringSet<TString, TSSetSpec> & stringSet, TDelta delta)
{
    // Copy last element to the back.
    appendValue(stringSet, back(stringSet));

    // Update limits.
    back(stringSet.limits) += delta;

    // Resize concat.
    resize(stringSet.concat, length(stringSet.concat) + delta);
}

// ----------------------------------------------------------------------------
// Function clear()
// ----------------------------------------------------------------------------

template <typename TText, typename TTextIndexSpec, typename TPattern, typename TPatternIndexSpec, typename TDistance, typename TSpec>
inline void
clear(Finder_<Index<TText, TTextIndexSpec>, Index<TPattern, TPatternIndexSpec>, Backtracking<TDistance, TSpec> > & finder)
{
    clear(finder.textStack);
    clear(finder.patternStack);
    clear(finder.scoreStack);
}

// ----------------------------------------------------------------------------
// Function _initState()
// ----------------------------------------------------------------------------

template <typename TText, typename TTextIndexSpec, typename TPattern, typename TPatternIndexSpec, typename TDistance, typename TSpec, typename TTextIterator, typename TPatternIterator>
inline void
_initState(Finder_<Index<TText, TTextIndexSpec>, Index<TPattern, TPatternIndexSpec>, Backtracking<TDistance, TSpec> > & finder,
           TTextIterator const & textIt, TPatternIterator const & patternIt)
{
    // Init iterators.
    appendValue(finder.textStack, textIt);
    appendValue(finder.patternStack, patternIt);

    _initScore(finder);

#ifdef SEQAN_DEBUG
        _printState(finder, StageInitial_());
#endif
}

// ----------------------------------------------------------------------------
// Function _initScore()
// ----------------------------------------------------------------------------

template <typename TText, typename TTextIndexSpec, typename TPattern, typename TPatternIndexSpec, typename TDistance, typename TSpec>
inline void
_initScore(Finder_<Index<TText, TTextIndexSpec>, Index<TPattern, TPatternIndexSpec>, Backtracking<TDistance, TSpec> > & finder)
{
    typedef Backtracking<TDistance, TSpec>              TBacktracking;
    typedef typename Score_<TBacktracking>::Type        TScore;

    // Push zero.
    appendValue(finder.scoreStack, TScore());
}

template <typename TText, typename TTextIndexSpec, typename TPattern, typename TPatternIndexSpec, typename TSpec>
inline void
_initScore(Finder_<Index<TText, TTextIndexSpec>, Index<TPattern, TPatternIndexSpec>, Backtracking<EditDistance, TSpec> > & finder)
{
    typedef Backtracking<EditDistance, TSpec>           TBacktracking;
    typedef typename Score_<TBacktracking>::Type        TScore;

    // Push a column with one zero cell.
    append(finder.scoreStack.limits, 1);
    appendValue(finder.scoreStack.concat, TScore());
}

// ----------------------------------------------------------------------------
// Function _pushState()
// ----------------------------------------------------------------------------

template <typename TText, typename TTextIndexSpec, typename TPattern, typename TPatternIndexSpec, typename TDistance, typename TSpec, typename TStage>
inline bool
_pushState(Finder_<Index<TText, TTextIndexSpec>, Index<TPattern, TPatternIndexSpec>, Backtracking<TDistance, TSpec> > & finder,
           TStage const & /* tag */)
{
    _pushIterators(finder, TStage());
    _pushScore(finder, TStage());

    if (_moveIteratorsDown(finder, TStage()))
    {
        _updateScore(finder, TStage());

#ifdef SEQAN_DEBUG
        _printPush(finder, TStage());
#endif
        return true;
    }

    _popIterators(finder, TStage());
    _popScore(finder, TStage());

    return false;
}

// ----------------------------------------------------------------------------
// Function _pushIterators()
// ----------------------------------------------------------------------------
// TODO(esiragusa): Specialize _pushIterators() for StageInitial_ and StageFinal_ of EditDistance

template <typename TText, typename TTextIndexSpec, typename TPattern, typename TPatternIndexSpec, typename TDistance, typename TSpec, typename TStage>
inline void
_pushIterators(Finder_<Index<TText, TTextIndexSpec>, Index<TPattern, TPatternIndexSpec>, Backtracking<TDistance, TSpec> > & finder,
               TStage const & /* tag */)
{
    appendValue(finder.patternStack, back(finder.patternStack));
    appendValue(finder.textStack, back(finder.textStack));
}

// ----------------------------------------------------------------------------
// Function _pushScore()
// ----------------------------------------------------------------------------

template <typename TText, typename TTextIndexSpec, typename TPattern, typename TPatternIndexSpec, typename TDistance, typename TSpec, typename TStage>
inline void
_pushScore(Finder_<Index<TText, TTextIndexSpec>, Index<TPattern, TPatternIndexSpec>, Backtracking<TDistance, TSpec> > & finder,
           TStage const & /* tag */)
{
    // Copy the last score on top of the stack.
    appendValue(finder.scoreStack, back(finder.scoreStack));
}

template <typename TText, typename TTextIndexSpec, typename TPattern, typename TPatternIndexSpec, typename TSpec>
inline void
_pushScore(Finder_<Index<TText, TTextIndexSpec>, Index<TPattern, TPatternIndexSpec>, Backtracking<HammingDistance, TSpec> > & /* finder */,
           StageExact_ const & /* tag */)
{
    // Do nothing.
}

template <typename TText, typename TTextIndexSpec, typename TPattern, typename TPatternIndexSpec, typename TSpec>
inline void
_pushScore(Finder_<Index<TText, TTextIndexSpec>, Index<TPattern, TPatternIndexSpec>, Backtracking<EditDistance, TSpec> > & finder,
           StageInitial_ const & /* tag */)
{
    // Copy the last column on top of the stack and add one cell.
    copyBackAndResize(finder.scoreStack, 1);
}

template <typename TText, typename TTextIndexSpec, typename TPattern, typename TPatternIndexSpec, typename TSpec>
inline void
_pushScore(Finder_<Index<TText, TTextIndexSpec>, Index<TPattern, TPatternIndexSpec>, Backtracking<EditDistance, TSpec> > & finder,
           StageUpper_ const & /* tag */)
{
    _pushScore(finder, StageInitial_());
}

template <typename TText, typename TTextIndexSpec, typename TPattern, typename TPatternIndexSpec, typename TSpec>
inline void
_pushScore(Finder_<Index<TText, TTextIndexSpec>, Index<TPattern, TPatternIndexSpec>, Backtracking<EditDistance, TSpec> > & finder,
           StageLower_ const & /* tag */)
{
    // Copy the last column on top of the stack and remove one cell.
    copyBackAndResize(finder.scoreStack, -1);
}

template <typename TText, typename TTextIndexSpec, typename TPattern, typename TPatternIndexSpec, typename TSpec>
inline void
_pushScore(Finder_<Index<TText, TTextIndexSpec>, Index<TPattern, TPatternIndexSpec>, Backtracking<EditDistance, TSpec> > & finder,
           StageFinal_ const & /* tag */)
{
    _pushScore(finder, StageLower_());
}

// ----------------------------------------------------------------------------
// Function _moveIteratorsDown()
// ----------------------------------------------------------------------------

template <typename TText, typename TTextIndexSpec, typename TPattern, typename TPatternIndexSpec, typename TDistance, typename TSpec, typename TStage>
inline bool
_moveIteratorsDown(Finder_<Index<TText, TTextIndexSpec>, Index<TPattern, TPatternIndexSpec>, Backtracking<TDistance, TSpec> > & finder,
                   TStage const & /* tag */)
{
    // Go down in text and pattern.
    return goDown(back(finder.textStack)) && goDown(back(finder.patternStack));
}

template <typename TText, typename TTextIndexSpec, typename TPattern, typename TPatternIndexSpec, typename TSpec>
inline bool
_moveIteratorsDown(Finder_<Index<TText, TTextIndexSpec>, Index<TPattern, TPatternIndexSpec>, Backtracking<HammingDistance, TSpec> > & finder,
                   StageExact_ const & /* tag */)
{
    typedef Index<TText, TTextIndexSpec>                                                TTextIndex;
    typedef Index<TPattern, TPatternIndexSpec>                                          TPatternIndex;
    typedef Backtracking<HammingDistance, TSpec>                                        TBacktracking;
    typedef typename TextIterator_<TTextIndex, TPatternIndex, TBacktracking>::Type      TTextIterator;
    typedef typename PatternIterator_<TTextIndex, TPatternIndex, TBacktracking>::Type   TPatternIterator;

    TTextIterator & textIt = back(finder.textStack);
    TPatternIterator & patternIt = back(finder.patternStack);

    // Go down in pattern and search pattern label in text.
    if (goDown(patternIt) && goDown(textIt, parentEdgeLabel(patternIt)))
        return true;

    // Otherwise go right in pattern.
    return _moveIteratorsRight(finder, StageExact_());
}

template <typename TText, typename TTextIndexSpec, typename TPattern, typename TPatternIndexSpec, typename TSpec>
inline bool
_moveIteratorsDown(Finder_<Index<TText, TTextIndexSpec>, Index<TPattern, TPatternIndexSpec>, Backtracking<EditDistance, TSpec> > & finder,
                   StageInitial_ const & /* tag */)
{
    // Go down in pattern.
    return goDown(back(finder.patternStack));
}

template <typename TText, typename TTextIndexSpec, typename TPattern, typename TPatternIndexSpec, typename TSpec>
inline bool
_moveIteratorsDown(Finder_<Index<TText, TTextIndexSpec>, Index<TPattern, TPatternIndexSpec>, Backtracking<EditDistance, TSpec> > & finder,
                   StageLower_ const & /* tag */)
{
    // Go down in text.
    return goDown(back(finder.textStack));
}

template <typename TText, typename TTextIndexSpec, typename TPattern, typename TPatternIndexSpec, typename TSpec>
inline bool
_moveIteratorsDown(Finder_<Index<TText, TTextIndexSpec>, Index<TPattern, TPatternIndexSpec>, Backtracking<EditDistance, TSpec> > & finder,
                   StageFinal_ const & /* tag */)
{
    return _moveIteratorsDown(finder, StageLower_());
}

// ----------------------------------------------------------------------------
// Function _nextState()
// ----------------------------------------------------------------------------

template <typename TText, typename TTextIndexSpec, typename TPattern, typename TPatternIndexSpec, typename TDistance, typename TSpec, typename TStage>
inline bool
_nextState(Finder_<Index<TText, TTextIndexSpec>, Index<TPattern, TPatternIndexSpec>, Backtracking<TDistance, TSpec> > & finder,
           TStage const & /* tag */)
{
    if (_moveIteratorsRight(finder, TStage()))
    {
        _updateScore(finder, TStage());

#ifdef SEQAN_DEBUG
        _printState(finder, TStage());
#endif
        return true;
    }

    return false;
}

// ----------------------------------------------------------------------------
// Function _moveIteratorsRight()
// ----------------------------------------------------------------------------

template <typename TText, typename TTextIndexSpec, typename TPattern, typename TPatternIndexSpec, typename TDistance, typename TSpec, typename TStage>
inline bool
_moveIteratorsRight(Finder_<Index<TText, TTextIndexSpec>, Index<TPattern, TPatternIndexSpec>, Backtracking<TDistance, TSpec> > & finder,
                    TStage const & /* tag */)
{
    // Try to go right in the pattern.
    if (goRight(back(finder.patternStack)))
    {
        return true;
    }
    // Try to go right in the text.
    else if (goRight(back(finder.textStack)))
    {
        // Move to the leftmost pattern.
//        back(finder.patternStack) = finder.patternLeftmost;

        eraseBack(finder.patternStack);
        appendValue(finder.patternStack, back(finder.patternStack));
        goDown(back(finder.patternStack));
        return true;
    }

    return false;
}

template <typename TText, typename TTextIndexSpec, typename TPattern, typename TPatternIndexSpec, typename TSpec>
inline bool
_moveIteratorsRight(Finder_<Index<TText, TTextIndexSpec>, Index<TPattern, TPatternIndexSpec>, Backtracking<HammingDistance, TSpec> > & finder,
                    StageExact_ const & /* tag */)
{
    typedef Index<TText, TTextIndexSpec>                                                TTextIndex;
    typedef Index<TPattern, TPatternIndexSpec>                                          TPatternIndex;
    typedef Backtracking<HammingDistance, TSpec>                                        TBacktracking;
    typedef typename TextIterator_<TTextIndex, TPatternIndex, TBacktracking>::Type      TTextIterator;
    typedef typename PatternIterator_<TTextIndex, TPatternIndex, TBacktracking>::Type   TPatternIterator;

    // TODO(esiragusa): Implement goRight(it, pattern).
//    TTextIterator & textIt = back(finder.textStack);
//    TPatternIterator & patternIt = back(finder.patternStack);
//    return goRight(patternIt) && goRight(textIt, parentEdgeLabel(patternIt));

    TPatternIterator & patternIt = back(finder.patternStack);

    // Try to go right in the pattern.
    while (goRight(patternIt))
    {
        // Move up in the text.
        back(finder.textStack) = finder.textStack[length(finder.textStack) - 2];

        // Search pattern label in text.
        TTextIterator & textIt = back(finder.textStack);
        if (goDown(textIt, parentEdgeLabel(patternIt)))
            return true;
    }

    return false;
}

template <typename TText, typename TTextIndexSpec, typename TPattern, typename TPatternIndexSpec, typename TSpec>
inline bool
_moveIteratorsRight(Finder_<Index<TText, TTextIndexSpec>, Index<TPattern, TPatternIndexSpec>, Backtracking<EditDistance, TSpec> > & finder,
                    StageInitial_ const & /* tag */)
{
    // Try to go right in the pattern.
    return goRight(back(finder.patternStack));
}

template <typename TText, typename TTextIndexSpec, typename TPattern, typename TPatternIndexSpec, typename TSpec>
inline bool
_moveIteratorsRight(Finder_<Index<TText, TTextIndexSpec>, Index<TPattern, TPatternIndexSpec>, Backtracking<EditDistance, TSpec> > & finder,
                    StageLower_ const & /* tag */)
{
    // Try to go right in the text.
    return goRight(back(finder.textStack));
}

// ----------------------------------------------------------------------------
// Function _popState()
// ----------------------------------------------------------------------------

template <typename TText, typename TTextIndexSpec, typename TPattern, typename TPatternIndexSpec, typename TDistance, typename TSpec, typename TStage>
inline void
_popState(Finder_<Index<TText, TTextIndexSpec>, Index<TPattern, TPatternIndexSpec>, Backtracking<TDistance, TSpec> > & finder,
          TStage const & /* tag */)
{
    _popIterators(finder, TStage());
    _popScore(finder, TStage());

#ifdef SEQAN_DEBUG
    _printPop(finder, TStage());
#endif
}

// ----------------------------------------------------------------------------
// Function _popIterators()
// ----------------------------------------------------------------------------
// TODO(esiragusa): Specialize _popIterators() for StageInitial_ and StageFinal_ of EditDistance

template <typename TText, typename TTextIndexSpec, typename TPattern, typename TPatternIndexSpec, typename TDistance, typename TSpec, typename TStage>
inline void
_popIterators(Finder_<Index<TText, TTextIndexSpec>, Index<TPattern, TPatternIndexSpec>, Backtracking<TDistance, TSpec> > & finder,
              TStage const & /* tag */)
{
    eraseBack(finder.textStack);
    eraseBack(finder.patternStack);
}

// ----------------------------------------------------------------------------
// Function _popScore()
// ----------------------------------------------------------------------------

template <typename TText, typename TTextIndexSpec, typename TPattern, typename TPatternIndexSpec, typename TDistance, typename TSpec, typename TStage>
inline void
_popScore(Finder_<Index<TText, TTextIndexSpec>, Index<TPattern, TPatternIndexSpec>, Backtracking<TDistance, TSpec> > & finder,
          TStage const & /* tag */)
{
    eraseBack(finder.scoreStack);
}

template <typename TText, typename TTextIndexSpec, typename TPattern, typename TPatternIndexSpec, typename TSpec>
inline void
_popScore(Finder_<Index<TText, TTextIndexSpec>, Index<TPattern, TPatternIndexSpec>, Backtracking<HammingDistance, TSpec> > & /* finder */,
          StageExact_ const & /* tag */)
{
    // Do nothing.
}

// ----------------------------------------------------------------------------
// Function _updateScore()
// ----------------------------------------------------------------------------

template <typename TText, typename TTextIndexSpec, typename TPattern, typename TPatternIndexSpec, typename TDistance, typename TSpec, typename TStage>
inline void
_updateScore(Finder_<Index<TText, TTextIndexSpec>, Index<TPattern, TPatternIndexSpec>, Backtracking<TDistance, TSpec> > & /* finder */,
             TStage const & /* tag */)
{}

template <typename TText, typename TTextIndexSpec, typename TPattern, typename TPatternIndexSpec, typename TSpec, typename TStage>
inline void
_updateScore(Finder_<Index<TText, TTextIndexSpec>, Index<TPattern, TPatternIndexSpec>, Backtracking<HammingDistance, TSpec> > & finder,
             TStage const & /* tag */)
{
    typedef Backtracking<HammingDistance, TSpec>    TBacktracking;
    typedef typename Score_<TBacktracking>::Type    TScore;

    // Compute score of text and pattern.
    TScore score = ordEqual(parentEdgeLabel(back(finder.textStack)), parentEdgeLabel(back(finder.patternStack))) ? 0 : 1;

    // Add score to previous score.
    back(finder.scoreStack) = value(finder.scoreStack, length(finder.scoreStack) - 2) + score;
}

template <typename TText, typename TTextIndexSpec, typename TPattern, typename TPatternIndexSpec, typename TSpec>
inline void
_updateScore(Finder_<Index<TText, TTextIndexSpec>, Index<TPattern, TPatternIndexSpec>, Backtracking<HammingDistance, TSpec> > & /* finder */,
             StageExact_ const & /* tag */)
{
    // Do nothing.
}

template <typename TText, typename TTextIndexSpec, typename TPattern, typename TPatternIndexSpec, typename TSpec>
inline void
_updateScore(Finder_<Index<TText, TTextIndexSpec>, Index<TPattern, TPatternIndexSpec>, Backtracking<EditDistance, TSpec> > & finder,
             StageInitial_ const & /* tag */)
{
    _updateVertexScore(back(finder.scoreStack), value(finder.scoreStack, length(finder.scoreStack) - 2), StageInitial_());
}

template <typename TText, typename TTextIndexSpec, typename TPattern, typename TPatternIndexSpec, typename TSpec>
inline void
_updateScore(Finder_<Index<TText, TTextIndexSpec>, Index<TPattern, TPatternIndexSpec>, Backtracking<EditDistance, TSpec> > & finder,
             StageUpper_ const & /* tag */)
{
    typedef Index<TPattern, TPatternIndexSpec>                              TPatternIndex;
    typedef typename Fibre<TPatternIndex, FibreText>::Type const            TPatternFibreText;
    typedef typename Infix<TPatternFibreText>::Type                         TPatternRepr;
    typedef typename Iterator<TPatternRepr, Standard>::Type                 TPatternReprIterator;

    // Get pattern read so far.
    TPatternRepr pattern = representative(back(finder.patternStack));
    TPatternReprIterator patternIt = begin(pattern, Standard());

    _updateVertexScore(back(finder.scoreStack),
                       value(finder.scoreStack, length(finder.scoreStack) - 2),
                       parentEdgeLabel(back(finder.textStack)),
                       patternIt,
                       StageUpper_());
}

template <typename TText, typename TTextIndexSpec, typename TPattern, typename TPatternIndexSpec, typename TSpec>
inline void
_updateScore(Finder_<Index<TText, TTextIndexSpec>, Index<TPattern, TPatternIndexSpec>, Backtracking<EditDistance, TSpec> > & finder,
             StageDiagonal_ const & /* tag */)
{
    typedef Index<TPattern, TPatternIndexSpec>                              TPatternIndex;
    typedef typename Fibre<TPatternIndex, FibreText>::Type const            TPatternFibreText;
    typedef typename Infix<TPatternFibreText>::Type                         TPatternRepr;
    typedef typename Iterator<TPatternRepr, Standard>::Type                 TPatternReprIterator;

    // Get last 2k + 1 pattern symbols.
    TPatternRepr pattern = representative(back(finder.patternStack));
    TPatternReprIterator patternIt = end(pattern, Standard()) - (2 * finder.maxScore + 1);

    _updateVertexScore(back(finder.scoreStack),
                       value(finder.scoreStack, length(finder.scoreStack) - 2),
                       parentEdgeLabel(back(finder.textStack)),
                       patternIt,
                       StageDiagonal_());
}

template <typename TText, typename TTextIndexSpec, typename TPattern, typename TPatternIndexSpec, typename TSpec>
inline void
_updateScore(Finder_<Index<TText, TTextIndexSpec>, Index<TPattern, TPatternIndexSpec>, Backtracking<EditDistance, TSpec> > & finder,
             StageLower_ const & /* tag */)
{
    typedef Index<TPattern, TPatternIndexSpec>                              TPatternIndex;
    typedef typename Fibre<TPatternIndex, FibreText>::Type const            TPatternFibreText;
    typedef typename Infix<TPatternFibreText>::Type                         TPatternRepr;
    typedef typename Iterator<TPatternRepr, Standard>::Type                 TPatternReprIterator;

    // Get last pattern symbols.
    TPatternRepr pattern = representative(back(finder.patternStack));
    TPatternReprIterator patternIt = end(pattern, Standard()) - (length(back(finder.scoreStack)));

    _updateVertexScore(back(finder.scoreStack),
                       value(finder.scoreStack, length(finder.scoreStack) - 2),
                       parentEdgeLabel(back(finder.textStack)),
                       patternIt,
                       StageLower_());
}

template <typename TText, typename TTextIndexSpec, typename TPattern, typename TPatternIndexSpec, typename TSpec>
inline void
_updateScore(Finder_<Index<TText, TTextIndexSpec>, Index<TPattern, TPatternIndexSpec>, Backtracking<EditDistance, TSpec> > & finder,
             StageFinal_ const & /* tag */)
{
    _updateVertexScore(back(finder.scoreStack),
                       value(finder.scoreStack, length(finder.scoreStack) - 2),
                       parentEdgeLabel(back(finder.textStack)),
                       parentEdgeLabel(back(finder.patternStack)),
                       StageFinal_());
}

// ----------------------------------------------------------------------------
// Function _setScoreThreshold()
// ----------------------------------------------------------------------------

template <typename TText, typename TTextIndexSpec, typename TPattern, typename TPatternIndexSpec, typename TDistance, typename TSpec, typename TMaxScore>
inline void
_setScoreThreshold(Finder_<Index<TText, TTextIndexSpec>, Index<TPattern, TPatternIndexSpec>, Backtracking<TDistance, TSpec> > & finder,
                  TMaxScore maxScore)
{
    finder.maxScore = maxScore;
}

// ----------------------------------------------------------------------------
// Function _getMinScore()
// ----------------------------------------------------------------------------

template <typename TText, typename TTextIndexSpec, typename TPattern, typename TPatternIndexSpec, typename TSpec>
inline typename Score_<Backtracking<HammingDistance, TSpec> >::Type
_getMinScore(Finder_<Index<TText, TTextIndexSpec>, Index<TPattern, TPatternIndexSpec>, Backtracking<HammingDistance, TSpec> > const & finder)
{
    return _getScore(finder);
}

template <typename TText, typename TTextIndexSpec, typename TPattern, typename TPatternIndexSpec, typename TSpec>
inline typename Score_<Backtracking<EditDistance, TSpec> >::Type
_getMinScore(Finder_<Index<TText, TTextIndexSpec>, Index<TPattern, TPatternIndexSpec>, Backtracking<EditDistance, TSpec> > const & finder)
{
    typedef Backtracking<EditDistance, TSpec>                       TBacktracking;
    typedef typename VertexScore_<TBacktracking>::ConstType         TVertexScore;

    TVertexScore column = back(finder.scoreStack);

    // Return the min value in column.
    return value(std::min_element(begin(column, Standard()), end(column, Standard())));
}

// ----------------------------------------------------------------------------
// Function _getScore()
// ----------------------------------------------------------------------------

template <typename TText, typename TTextIndexSpec, typename TPattern, typename TPatternIndexSpec, typename TSpec>
inline typename Score_<Backtracking<HammingDistance, TSpec> >::Type
_getScore(Finder_<Index<TText, TTextIndexSpec>, Index<TPattern, TPatternIndexSpec>, Backtracking<HammingDistance, TSpec> > const & finder)
{
    // Return the last value.
    return back(finder.scoreStack);
}

template <typename TText, typename TTextIndexSpec, typename TPattern, typename TPatternIndexSpec, typename TSpec>
inline typename Score_<Backtracking<EditDistance, TSpec> >::Type
_getScore(Finder_<Index<TText, TTextIndexSpec>, Index<TPattern, TPatternIndexSpec>, Backtracking<EditDistance, TSpec> > const & finder)
{
    // Return the value of last cell in column.
    return back(back(finder.scoreStack));
}

// ----------------------------------------------------------------------------
// Function _textIterator()
// ----------------------------------------------------------------------------

template <typename TText, typename TTextIndexSpec, typename TPattern, typename TPatternIndexSpec, typename TDistance, typename TSpec>
SEQAN_HOST_DEVICE inline
typename TextIterator_<Index<TText, TTextIndexSpec>, Index<TPattern, TPatternIndexSpec>, Backtracking<TDistance, TSpec> >::Type &
_textIterator(Finder_<Index<TText, TTextIndexSpec>, Index<TPattern, TPatternIndexSpec>, Backtracking<TDistance, TSpec> > & finder)
{
    return back(finder.textStack);
}

template <typename TText, typename TTextIndexSpec, typename TPattern, typename TPatternIndexSpec, typename TDistance, typename TSpec>
SEQAN_HOST_DEVICE inline
typename TextIterator_<Index<TText, TTextIndexSpec>, Index<TPattern, TPatternIndexSpec>, Backtracking<TDistance, TSpec> >::Type const &
_textIterator(Finder_<Index<TText, TTextIndexSpec>, Index<TPattern, TPatternIndexSpec>, Backtracking<TDistance, TSpec> > const & finder)
{
    return back(finder.textStack);
}

// ----------------------------------------------------------------------------
// Function _patternIterator()
// ----------------------------------------------------------------------------

template <typename TText, typename TTextIndexSpec, typename TPattern, typename TPatternIndexSpec, typename TDistance, typename TSpec>
SEQAN_HOST_DEVICE inline
typename PatternIterator_<Index<TText, TTextIndexSpec>, Index<TPattern, TPatternIndexSpec>, Backtracking<TDistance, TSpec> >::Type &
_patternIterator(Finder_<Index<TText, TTextIndexSpec>, Index<TPattern, TPatternIndexSpec>, Backtracking<TDistance, TSpec> > & finder)
{
    return back(finder.patternStack);
}

template <typename TText, typename TTextIndexSpec, typename TPattern, typename TPatternIndexSpec, typename TDistance, typename TSpec>
SEQAN_HOST_DEVICE inline
typename PatternIterator_<Index<TText, TTextIndexSpec>, Index<TPattern, TPatternIndexSpec>, Backtracking<TDistance, TSpec> >::Type const &
_patternIterator(Finder_<Index<TText, TTextIndexSpec>, Index<TPattern, TPatternIndexSpec>, Backtracking<TDistance, TSpec> > const & finder)
{
    return back(finder.patternStack);
}

// ----------------------------------------------------------------------------
// Function _inTerminalState()
// ----------------------------------------------------------------------------

template <typename TText, typename TTextIndexSpec, typename TPattern, typename TPatternIndexSpec, typename TDistance, typename TSpec, typename TStage>
inline bool
_inTerminalState(Finder_<Index<TText, TTextIndexSpec>, Index<TPattern, TPatternIndexSpec>, Backtracking<TDistance, TSpec> > const & /* finder */,
                 TStage const & /* tag */)
{
    // The current state is not terminal by default.
    return false;
}

template <typename TText, typename TTextIndexSpec, typename TPattern, typename TPatternIndexSpec, typename TSpec>
inline bool
_inTerminalState(Finder_<Index<TText, TTextIndexSpec>, Index<TPattern, TPatternIndexSpec>, Backtracking<HammingDistance, TSpec> > const & finder,
                 StageFinal_ const & /* tag */)
{
    // Is the score within the max score?
    return _getScore(finder) <= finder.maxScore;
}

template <typename TText, typename TTextIndexSpec, typename TPattern, typename TPatternIndexSpec, typename TSpec>
inline bool
_inTerminalState(Finder_<Index<TText, TTextIndexSpec>, Index<TPattern, TPatternIndexSpec>, Backtracking<EditDistance, TSpec> > const & finder,
                 StageLower_ const & /* tag */)
{
    // Is the score within the max score?
    return _getScore(finder) <= finder.maxScore;
}

template <typename TText, typename TTextIndexSpec, typename TPattern, typename TPatternIndexSpec, typename TSpec>
inline bool
_inTerminalState(Finder_<Index<TText, TTextIndexSpec>, Index<TPattern, TPatternIndexSpec>, Backtracking<EditDistance, TSpec> > const & finder,
                 StageFinal_ const & /* tag */)
{
    // Is the score within the max score?
    return _getScore(finder) <= finder.maxScore;
}

// ----------------------------------------------------------------------------
// Function _inActiveState()
// ----------------------------------------------------------------------------

template <typename TText, typename TTextIndexSpec, typename TPattern, typename TPatternIndexSpec, typename TDistance, typename TSpec, typename TStage>
inline bool
_inActiveState(Finder_<Index<TText, TTextIndexSpec>, Index<TPattern, TPatternIndexSpec>, Backtracking<TDistance, TSpec> > const & finder,
               TStage const & /* tag */)
{
    // Is the minimum score within the max score?
    return _getMinScore(finder) <= finder.maxScore;
}

template <typename TText, typename TTextIndexSpec, typename TPattern, typename TPatternIndexSpec, typename TSpec>
inline bool
_inActiveState(Finder_<Index<TText, TTextIndexSpec>, Index<TPattern, TPatternIndexSpec>, Backtracking<HammingDistance, TSpec> > const & /* finder */,
               StageExact_ const & /* tag */)
{
    // Exact search only walks through active states.
    return true;
}

// ----------------------------------------------------------------------------
// Function _moveToExactStage()
// ----------------------------------------------------------------------------

template <typename TText, typename TTextIndexSpec, typename TPattern, typename TPatternIndexSpec, typename TDistance, typename TSpec, typename TStage>
inline bool
_moveToStageExact(Finder_<Index<TText, TTextIndexSpec>, Index<TPattern, TPatternIndexSpec>, Backtracking<TDistance, TSpec> > const & /* finder */,
                  TStage const & /* tag */)
{
    // By default there is no exact stage.
    return false;
}

template <typename TText, typename TTextIndexSpec, typename TPattern, typename TPatternIndexSpec, typename TSpec>
inline bool
_moveToStageExact(Finder_<Index<TText, TTextIndexSpec>, Index<TPattern, TPatternIndexSpec>, Backtracking<HammingDistance, TSpec> > const & finder,
                  StageInitial_ const & /* tag */)
{
    // Was the maximum score attained?
    return _getMinScore(finder) == finder.maxScore;
}

//template <typename TText, typename TTextIndexSpec, typename TPattern, typename TPatternIndexSpec, typename TSpec, typename TStage>
//inline bool
//_moveToStageExact(Finder_<Index<TText, TTextIndexSpec>, Index<TPattern, TPatternIndexSpec>, Backtracking<EditDistance, TSpec> > const & finder,
//                  TStage const & /* tag */)
//{
//    // TODO(esiragusa): Implement exact search speedup for EditDistance.
//    return _getMinScore(finder) == finder.maxScore;
//}

// ----------------------------------------------------------------------------
// Function _moveToNextStage()
// ----------------------------------------------------------------------------

template <typename TText, typename TTextIndexSpec, typename TPattern, typename TPatternIndexSpec, typename TDistance, typename TSpec, typename TStage>
inline bool
_moveToNextStage(Finder_<Index<TText, TTextIndexSpec>, Index<TPattern, TPatternIndexSpec>, Backtracking<TDistance, TSpec> > const & /* finder */,
                 TStage const & /* tag */)
{
    // By default there is no next stage.
    return false;
}

template <typename TText, typename TTextIndexSpec, typename TPattern, typename TPatternIndexSpec, typename TSpec>
inline bool
_moveToNextStage(Finder_<Index<TText, TTextIndexSpec>, Index<TPattern, TPatternIndexSpec>, Backtracking<HammingDistance, TSpec> > const & finder,
                 StageInitial_ const & /* tag */)
{
    return isRightTerminal(back(finder.patternStack));
}

template <typename TText, typename TTextIndexSpec, typename TPattern, typename TPatternIndexSpec, typename TSpec>
inline bool
_moveToNextStage(Finder_<Index<TText, TTextIndexSpec>, Index<TPattern, TPatternIndexSpec>, Backtracking<HammingDistance, TSpec> > const & finder,
                 StageExact_ const & /* tag */)
{
    return isRightTerminal(back(finder.patternStack));
}

template <typename TText, typename TTextIndexSpec, typename TPattern, typename TPatternIndexSpec, typename TSpec>
inline bool
_moveToNextStage(Finder_<Index<TText, TTextIndexSpec>, Index<TPattern, TPatternIndexSpec>, Backtracking<EditDistance, TSpec> > const & finder,
                 StageInitial_ const & /* tag */)
{
    // Move to upper stage when k pattern symbols have been consumed.
    return (repLength(back(finder.patternStack)) >= finder.maxScore) || isRightTerminal(back(finder.patternStack));
}

template <typename TText, typename TTextIndexSpec, typename TPattern, typename TPatternIndexSpec, typename TSpec>
inline bool
_moveToNextStage(Finder_<Index<TText, TTextIndexSpec>, Index<TPattern, TPatternIndexSpec>, Backtracking<EditDistance, TSpec> > const & finder,
                 StageUpper_ const & /* tag */)
{
    // Move to diagonal stage when the diagonal has size 2k + 1.
    return (length(back(finder.scoreStack)) >= 2u * finder.maxScore + 1u) || isRightTerminal(back(finder.patternStack));
}

template <typename TText, typename TTextIndexSpec, typename TPattern, typename TPatternIndexSpec, typename TSpec>
inline bool
_moveToNextStage(Finder_<Index<TText, TTextIndexSpec>, Index<TPattern, TPatternIndexSpec>, Backtracking<EditDistance, TSpec> > const & finder,
                 StageDiagonal_ const & /* tag */)
{
    // Move to lower stage when all pattern symbols have been consumed.
    return isRightTerminal(back(finder.patternStack));
}

template <typename TText, typename TTextIndexSpec, typename TPattern, typename TPatternIndexSpec, typename TSpec>
inline bool
_moveToNextStage(Finder_<Index<TText, TTextIndexSpec>, Index<TPattern, TPatternIndexSpec>, Backtracking<EditDistance, TSpec> > const & finder,
                 StageLower_ const & /* tag */)
{
    // Move to final stage when there is only one cell left to compute.
    return length(back(finder.scoreStack)) <= 2;
}

// ----------------------------------------------------------------------------
// Function _stayInCurrentStage()
// ----------------------------------------------------------------------------

template <typename TText, typename TTextIndexSpec, typename TPattern, typename TPatternIndexSpec, typename TDistance, typename TSpec, typename TStage>
inline bool
_stayInCurrentStage(Finder_<Index<TText, TTextIndexSpec>, Index<TPattern, TPatternIndexSpec>, Backtracking<TDistance, TSpec> > const & /* finder */,
                    TStage const & /* tag */)
{
    // Stay in the current stage by default.
    return true;
}

template <typename TText, typename TTextIndexSpec, typename TPattern, typename TPatternIndexSpec, typename TSpec>
inline bool
_stayInCurrentStage(Finder_<Index<TText, TTextIndexSpec>, Index<TPattern, TPatternIndexSpec>, Backtracking<EditDistance, TSpec> > const & finder,
                    StageInitial_ const & /* tag */)
{
    return repLength(back(finder.patternStack)) < finder.maxScore;
}

template <typename TText, typename TTextIndexSpec, typename TPattern, typename TPatternIndexSpec, typename TSpec>
inline bool
_stayInCurrentStage(Finder_<Index<TText, TTextIndexSpec>, Index<TPattern, TPatternIndexSpec>, Backtracking<EditDistance, TSpec> > const & finder,
                    StageUpper_ const & /* tag */)
{
    return length(back(finder.scoreStack)) < 2u * finder.maxScore + 1u;
}

template <typename TText, typename TTextIndexSpec, typename TPattern, typename TPatternIndexSpec, typename TSpec>
inline bool
_stayInCurrentStage(Finder_<Index<TText, TTextIndexSpec>, Index<TPattern, TPatternIndexSpec>, Backtracking<EditDistance, TSpec> > const & finder,
                    StageLower_ const & /* tag */)
{
    return length(back(finder.scoreStack)) > 2;
}

// ----------------------------------------------------------------------------
// Functions _print*()
// ----------------------------------------------------------------------------
// NOTE(esiragusa): Debug functions.

template <typename TText, typename TTextIndexSpec, typename TPattern, typename TPatternIndexSpec, typename TDistance, typename TSpec, typename TStage>
inline void
_printCall(Finder_<Index<TText, TTextIndexSpec>, Index<TPattern, TPatternIndexSpec>, Backtracking<TDistance, TSpec> > const & finder,
           TStage const & /* tag */)
{
    std::cout << "call:           "; _printFindSignature(finder, TStage());
    std::cout << "past text:      " << representative(back(finder.textStack)) << std::endl;
    std::cout << "past pattern:   " << representative(back(finder.patternStack)) << std::endl;
}

template <typename TText, typename TTextIndexSpec, typename TPattern, typename TPatternIndexSpec, typename TDistance, typename TSpec, typename TStage>
inline void
_printReturn(Finder_<Index<TText, TTextIndexSpec>, Index<TPattern, TPatternIndexSpec>, Backtracking<TDistance, TSpec> > const & finder,
             TStage const & /* tag */)
{
    std::cout << "return:         "; _printFindSignature(finder, TStage());
}

template <typename TText, typename TTextIndexSpec, typename TPattern, typename TPatternIndexSpec, typename TDistance, typename TSpec, typename TStage>
inline void
_printPush(Finder_<Index<TText, TTextIndexSpec>, Index<TPattern, TPatternIndexSpec>, Backtracking<TDistance, TSpec> > const & finder,
           TStage const & /* tag */)
{
    std::cout << "push:           "; _printFindSignature(finder, TStage());
    _printState(finder, TStage());
}

template <typename TText, typename TTextIndexSpec, typename TPattern, typename TPatternIndexSpec, typename TDistance, typename TSpec, typename TStage>
inline void
_printPop(Finder_<Index<TText, TTextIndexSpec>, Index<TPattern, TPatternIndexSpec>, Backtracking<TDistance, TSpec> > const & finder,
          TStage const & /* tag */)
{
    std::cout << "pop:            "; _printFindSignature(finder, TStage());
    _printState(finder, TStage());
}

template <typename TText, typename TTextIndexSpec, typename TPattern, typename TPatternIndexSpec, typename TDistance, typename TSpec, typename TStage>
inline void
_printFindSignature(Finder_<Index<TText, TTextIndexSpec>, Index<TPattern, TPatternIndexSpec>, Backtracking<TDistance, TSpec> > const & /* finder */,
                    TStage const & /* tag */)
{
    std::cout << "<TDistance, TStage>" << std::endl;
}

template <typename TText, typename TTextIndexSpec, typename TPattern, typename TPatternIndexSpec, typename TSpec>
inline void
_printFindSignature(Finder_<Index<TText, TTextIndexSpec>, Index<TPattern, TPatternIndexSpec>, Backtracking<HammingDistance, TSpec> > const & /* finder */,
                    StageInitial_ const & /* tag */)
{
    std::cout << "<HammingDistance, StageInitial>" << std::endl;
}

template <typename TText, typename TTextIndexSpec, typename TPattern, typename TPatternIndexSpec, typename TSpec>
inline void
_printFindSignature(Finder_<Index<TText, TTextIndexSpec>, Index<TPattern, TPatternIndexSpec>, Backtracking<HammingDistance, TSpec> > const & /* finder */,
                    StageFinal_ const & /* tag */)
{
    std::cout << "<HammingDistance, StageFinal>" << std::endl;
}

template <typename TText, typename TTextIndexSpec, typename TPattern, typename TPatternIndexSpec, typename TSpec>
inline void
_printFindSignature(Finder_<Index<TText, TTextIndexSpec>, Index<TPattern, TPatternIndexSpec>, Backtracking<HammingDistance, TSpec> > const & /* finder */,
                    StageExact_ const & /* tag */)
{
    std::cout << "<HammingDistance, StageExact>" << std::endl;
}

template <typename TText, typename TTextIndexSpec, typename TPattern, typename TPatternIndexSpec, typename TSpec>
inline void
_printFindSignature(Finder_<Index<TText, TTextIndexSpec>, Index<TPattern, TPatternIndexSpec>, Backtracking<EditDistance, TSpec> > const & /* finder */,
                    StageInitial_ const & /* tag */)
{
    std::cout << "<EditDistance, StageInitial>" << std::endl;
}

template <typename TText, typename TTextIndexSpec, typename TPattern, typename TPatternIndexSpec, typename TSpec>
inline void
_printFindSignature(Finder_<Index<TText, TTextIndexSpec>, Index<TPattern, TPatternIndexSpec>, Backtracking<EditDistance, TSpec> > const & /* finder */,
                    StageUpper_ const & /* tag */)
{
    std::cout << "<EditDistance, StageUpper>" << std::endl;
}

template <typename TText, typename TTextIndexSpec, typename TPattern, typename TPatternIndexSpec, typename TSpec>
inline void
_printFindSignature(Finder_<Index<TText, TTextIndexSpec>, Index<TPattern, TPatternIndexSpec>, Backtracking<EditDistance, TSpec> > const & /* finder */,
                    StageDiagonal_ const & /* tag */)
{
    std::cout << "<EditDistance, StageDiagonal>" << std::endl;
}

template <typename TText, typename TTextIndexSpec, typename TPattern, typename TPatternIndexSpec, typename TSpec>
inline void
_printFindSignature(Finder_<Index<TText, TTextIndexSpec>, Index<TPattern, TPatternIndexSpec>, Backtracking<EditDistance, TSpec> > const & /* finder */,
                    StageLower_ const & /* tag */)
{
    std::cout << "<EditDistance, StageLower>" << std::endl;
}

template <typename TText, typename TTextIndexSpec, typename TPattern, typename TPatternIndexSpec, typename TSpec>
inline void
_printFindSignature(Finder_<Index<TText, TTextIndexSpec>, Index<TPattern, TPatternIndexSpec>, Backtracking<EditDistance, TSpec> > const & /* finder */,
                    StageFinal_ const & /* tag */)
{
    std::cout << "<EditDistance, StageFinal>" << std::endl;
}

template <typename TText, typename TTextIndexSpec, typename TPattern, typename TPatternIndexSpec, typename TSpec, typename TStage>
inline void
_printState(Finder_<Index<TText, TTextIndexSpec>, Index<TPattern, TPatternIndexSpec>, Backtracking<HammingDistance, TSpec> > const & finder,
            TStage const & /* tag */)
{
    std::cout << "text:           " << parentEdgeLabel(back(finder.textStack)) << std::endl;
    std::cout << "pattern:        " << parentEdgeLabel(back(finder.patternStack)) << std::endl;
    std::cout << "errors:         " << static_cast<unsigned>(_getScore(finder)) << std::endl;
}

template <typename TText, typename TTextIndexSpec, typename TPattern, typename TPatternIndexSpec, typename TSpec, typename TStage>
inline void
_printState(Finder_<Index<TText, TTextIndexSpec>, Index<TPattern, TPatternIndexSpec>, Backtracking<EditDistance, TSpec> > const & finder,
            TStage const & /* tag */)
{
    std::cout << "text:           " << parentEdgeLabel(back(finder.textStack)) << std::endl;
    std::cout << "pattern:        " << parentEdgeLabel(back(finder.patternStack)) << std::endl;
    std::cout << "column:         " << "|";
    std::copy(begin(back(finder.scoreStack), Standard()),
              end(back(finder.scoreStack), Standard()),
              std::ostream_iterator<int>(std::cout, "|"));
    std::cout << std::endl;
}

// ----------------------------------------------------------------------------
// Function _find()
// ----------------------------------------------------------------------------

template <typename TText, typename TTextIndexSpec, typename TPattern, typename TPatternIndexSpec, typename TDistance, typename TSpec, typename TStage, typename TDelegate>
inline void
_find(Finder_<Index<TText, TTextIndexSpec>, Index<TPattern, TPatternIndexSpec>, Backtracking<TDistance, TSpec> > & finder,
      TDelegate & delegate,
      TStage const & /* tag */)
{
    typedef Backtracking<TDistance, TSpec>                      TBacktracking;
    typedef typename NextStage_<TBacktracking, TStage>::Type    TNextStage;

#ifdef SEQAN_DEBUG
    _printCall(finder, TStage());
#endif

    if (_moveToStageExact(finder, TStage()))
    {
        _find(finder, delegate, StageExact_());
    }
    else
    {
        if (_moveToNextStage(finder, TStage()))
        {
            _find(finder, delegate, TNextStage());
        }
        if (_stayInCurrentStage(finder, TStage()))
        {
            if (_inTerminalState(finder, TStage()))
            {
                // Inversion of control.
                delegate(finder);
            }
            else if (_inActiveState(finder, TStage()))
            {
                if (_pushState(finder, TStage()))
                {
                    do
                    {
                        _find(finder, delegate, TStage());
                    }
                    while (_nextState(finder, TStage()));

                    _popState(finder, TStage());
                }
            }
        }
    }

#ifdef SEQAN_DEBUG
    _printReturn(finder, TStage());
#endif
}

// ----------------------------------------------------------------------------
// Function _find()
// ----------------------------------------------------------------------------

template <typename TText, typename TTextIndexSpec, typename TPattern, typename TPatternIndexSpec, typename TDistance, typename TSpec, typename TValue, typename TDelegate>
inline void
_find(Finder_<Index<TText, TTextIndexSpec>, Index<TPattern, TPatternIndexSpec>, Backtracking<TDistance, TSpec> > & finder,
      Index<TText, TTextIndexSpec> & text,
      Index<TPattern, TPatternIndexSpec> & pattern,
      TValue maxScore,
      TDelegate & delegate)
{
    typedef Index<TText, TTextIndexSpec>                                                TTextIndex;
    typedef Index<TPattern, TPatternIndexSpec>                                          TPatternIndex;
    typedef Backtracking<TDistance, TSpec>                                              TBacktracking;
    typedef typename TextIterator_<TTextIndex, TPatternIndex, TBacktracking>::Type      TTextIterator;
    typedef typename PatternIterator_<TTextIndex, TPatternIndex, TBacktracking>::Type   TPatternIterator;

    TTextIterator textIt(text);
    TPatternIterator patternIt(pattern);

    _setScoreThreshold(finder, maxScore);
    _initState(finder, textIt, patternIt);
    _find(finder, delegate, StageInitial_());
    _popState(finder, StageInitial_());
}

}

#endif  // #ifndef SEQAN_FIND_BACKTRACKING_MULTIPLE_H_
