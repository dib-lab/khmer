// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2015, Knut Reinert, FU Berlin
// Copyright (c) 2013 NVIDIA Corporation
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
//     * Neither the name of NVIDIA Corporation nor the names of
//       its contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL NVIDIA CORPORATION BE LIABLE
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
// Approximate string matching via backtracking on a substring index.
// ==========================================================================

#ifndef SEQAN_INDEX_FIND_INDEX_H_
#define SEQAN_INDEX_FIND_INDEX_H_

namespace seqan {

// ============================================================================
// Tags
// ============================================================================

// ----------------------------------------------------------------------------
// Tags for backtracking specializations
// ----------------------------------------------------------------------------

struct BacktrackingSemiGlobal_;
typedef Tag<BacktrackingSemiGlobal_>    BacktrackingSemiGlobal;

struct BacktrackingGlobal_;
typedef Tag<BacktrackingGlobal_>        BacktrackingGlobal;

// ----------------------------------------------------------------------------
// Tag Backtracking
// ----------------------------------------------------------------------------

template <typename TDistance = HammingDistance, typename TSpec = BacktrackingSemiGlobal>
struct Backtracking {};

// ============================================================================
// Metafunctions
// ============================================================================

// ----------------------------------------------------------------------------
// Metafunction TextIterator_
// ----------------------------------------------------------------------------

template <typename TText, typename TIndexSpec, typename TPattern, typename TSpec>
struct TextIterator_<Index<TText, TIndexSpec>, TPattern, TSpec>
{
    typedef typename Iterator<Index<TText, TIndexSpec>, TopDown<> >::Type  Type;
};

// ----------------------------------------------------------------------------
// Metafunction TextIterator_
// ----------------------------------------------------------------------------

template <typename TText, typename TIndexSpec, typename TPattern, typename TDistance, typename TSpec>
struct TextIterator_<Index<TText, TIndexSpec>, TPattern, Backtracking<TDistance, TSpec> >
{
    typedef typename Iterator<Index<TText, TIndexSpec>, TopDown<ParentLinks<> > >::Type  Type;
};

// ============================================================================
// Classes
// ============================================================================

// ----------------------------------------------------------------------------
// Class Finder
// ----------------------------------------------------------------------------

template <typename TText, typename TIndexSpec, typename TPattern, typename TSpec>
struct Finder_<Index<TText, TIndexSpec>, TPattern, TSpec>
{
    typedef Index<TText, TIndexSpec>                                    TIndex;
    typedef typename TextIterator_<TIndex, TPattern, TSpec>::Type       TTextIterator;
    typedef typename PatternIterator_<TIndex, TPattern, TSpec>::Type    TPatternIterator;
    typedef typename Score_<TSpec>::Type                                TScore;

    TTextIterator       _textIt;
    TPatternIterator    _patternIt;
    TScore              _scoreThreshold;
    TScore              _score;

    SEQAN_HOST_DEVICE
    Finder_() :
        _scoreThreshold(),
        _score()
    {}

    SEQAN_HOST_DEVICE
    Finder_(TIndex /* const */ & index) :
        _textIt(index),
        _scoreThreshold(),
        _score()
    {}

    SEQAN_HOST_DEVICE
    Finder_(TTextIterator const & textIt) :
        _textIt(textIt),
        _scoreThreshold(),
        _score()
    {}
};

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function clear()
// ----------------------------------------------------------------------------

template <typename TText, typename TIndexSpec, typename TPattern, typename TSpec>
SEQAN_HOST_DEVICE inline void
clear(Finder_<Index<TText, TIndexSpec>, TPattern, TSpec> & finder)
{
    // NOTE(esiragusa): should clear() be called on text/patternIt?
    goRoot(_textIterator(finder));
    // NOTE(esiragusa): if find wasn't called yet, _patternIterator is uninitialized.
//    goBegin(_patternIterator(finder));
    finder._score = typename Score_<TSpec>::Type();
}

// ----------------------------------------------------------------------------
// Function find_()
// ----------------------------------------------------------------------------

template <typename TText, typename TIndexSpec, typename TPattern, typename TDelegate>
SEQAN_HOST_DEVICE inline void
_find(Finder_<Index<TText, TIndexSpec>, TPattern, FinderSTree> & finder,
      TPattern const & pattern,
      TDelegate & delegate)
{
    if (goDown(_textIterator(finder), pattern))
    {
        // TODO(esiragusa): update _patternIterator.
        delegate(finder);
    }
}

// ----------------------------------------------------------------------------
// Function _getVertexScore()
// ----------------------------------------------------------------------------

template <typename TText, typename TPattern, typename TSpec>
SEQAN_HOST_DEVICE inline typename Score_<Backtracking<HammingDistance, TSpec> >::Type
_getVertexScore(Finder_<TText, TPattern, Backtracking<HammingDistance, TSpec> > const & finder)
{
    return !ordEqual(parentEdgeLabel(_textIterator(finder)), value(_patternIterator(finder)));
}

// ----------------------------------------------------------------------------
// Function _printState()
// ----------------------------------------------------------------------------

template <typename TText, typename TIndexSpec, typename TPattern, typename TSpec>
SEQAN_HOST_DEVICE inline void
_printState(Finder_<Index<TText, TIndexSpec>, TPattern, Backtracking<HammingDistance, TSpec> > & finder)
{
    std::cout << "Text:        " << parentEdgeLabel(_textIterator(finder)) << std::endl;
    std::cout << "Pattern:     " << value(_patternIterator(finder)) << std::endl;
    std::cout << "Text Len:    " << repLength(_textIterator(finder)) << std::endl;
    std::cout << "Pattern Len: " << position(_patternIterator(finder)) + 1 << std::endl;
    std::cout << "Errors:      " << static_cast<unsigned>(_getScore(finder)) << std::endl;
    std::cout << "Max errors:  " << static_cast<unsigned>(finder._scoreThreshold) << std::endl;
}

// ----------------------------------------------------------------------------
// Function _find()
// ----------------------------------------------------------------------------

template <typename TText, typename TIndexSpec, typename TPattern, typename TSpec, typename TDelegate>
SEQAN_HOST_DEVICE inline void
_find(Finder_<Index<TText, TIndexSpec>, TPattern, Backtracking<HammingDistance, TSpec> > & finder,
      TPattern const & pattern,
      TDelegate & delegate)
{
    typedef Index<TText, TIndexSpec>                                        TIndex;
    typedef Backtracking<HammingDistance, TSpec>                            TFinderSpec;
    typedef typename TextIterator_<TIndex, TPattern, TFinderSpec>::Type     TTextIterator;
    typedef typename PatternIterator_<TIndex, TPattern, TFinderSpec>::Type  TPatternIterator;

    _patternIterator(finder) = begin(pattern);

    TTextIterator & textIt = _textIterator(finder);
    TPatternIterator & patternIt = _patternIterator(finder);

    do
    {
        // Exact case.
        if (finder._score == finder._scoreThreshold)
        {
            if (goDown(textIt, suffix(pattern, position(patternIt))))
            {
                delegate(finder);
            }

            goUp(textIt);
        }

        // Approximate case.
        else if (finder._score < finder._scoreThreshold)
        {
            // Base case.
            if (atEnd(patternIt))
            {
                delegate(finder);
            }

            // Recursive case.
            else if (goDown(textIt))
            {
                finder._score += _getVertexScore(finder);
                goNext(patternIt);
                continue;
            }
        }

        // Backtrack.
        do
        {
            // Termination.
            if (isRoot(textIt)) break;

            goPrevious(patternIt);
            finder._score -= _getVertexScore(finder);
        }
        while (!goRight(textIt) && goUp(textIt));

        // Termination.
        if (isRoot(textIt)) break;

        finder._score += _getVertexScore(finder);
        goNext(patternIt);
    }
    while (true);
}

}

#endif  // #ifndef SEQAN_INDEX_FIND_INDEX_H_
