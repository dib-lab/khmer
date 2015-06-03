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

#ifndef SEQAN_FIND_FIND_BASE_H_
#define SEQAN_FIND_FIND_BASE_H_

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

template <typename TText, typename TPattern, typename TSpec = void>
struct Finder_;

// ============================================================================
// Metafunctions
// ============================================================================

// ----------------------------------------------------------------------------
// Metafunction View                                                   [Finder]
// ----------------------------------------------------------------------------

template <typename TText, typename TPattern, typename TSpec>
struct View<Finder_<TText, TPattern, TSpec> >
{
    typedef Finder_<typename View<TText>::Type, typename View<TPattern>::Type, TSpec>   Type;
};

// ----------------------------------------------------------------------------
// Metafunction RemoveView                                             [Finder]
// ----------------------------------------------------------------------------

template <typename TText, typename TPattern, typename TSpec>
struct RemoveView<Finder_<TText, TPattern, TSpec> >
{
    typedef Finder_<typename RemoveView<TText>::Type, typename RemoveView<TPattern>::Type, TSpec>   Type;
};

// ----------------------------------------------------------------------------
// Metafunction IsView                                                 [Finder]
// ----------------------------------------------------------------------------

template <typename TText, typename TPattern, typename TSpec>
struct IsView<Finder_<TText, TPattern, TSpec> > : IsView<TText> {};

// ----------------------------------------------------------------------------
// Metafunction IsDevice                                               [Finder]
// ----------------------------------------------------------------------------

template <typename TText, typename TPattern, typename TSpec>
struct IsDevice<Finder_<TText, TPattern, TSpec> > : IsDevice<TText> {};

// ----------------------------------------------------------------------------
// Metafunction View                                                  [Pattern]
// ----------------------------------------------------------------------------

template <typename TNeedle, typename TSpec>
struct View<Pattern<TNeedle, TSpec> >
{
    typedef Pattern<typename View<TNeedle>::Type, TSpec>    Type;
};

// ----------------------------------------------------------------------------
// Metafunction RemoveView                                            [Pattern]
// ----------------------------------------------------------------------------

template <typename TNeedle, typename TSpec>
struct RemoveView<Pattern<TNeedle, TSpec> >
{
    typedef Pattern<typename RemoveView<TNeedle>::Type, TSpec>    Type;
};

// ----------------------------------------------------------------------------
// Metafunction IsView                                                [Pattern]
// ----------------------------------------------------------------------------

template <typename TNeedle, typename TSpec>
struct IsView<Pattern<TNeedle, TSpec> > : IsView<TNeedle> {};

// ----------------------------------------------------------------------------
// Metafunction IsDevice                                              [Pattern]
// ----------------------------------------------------------------------------

template <typename TNeedle, typename TSpec>
struct IsDevice<Pattern<TNeedle, TSpec> > : IsDevice<TNeedle> {};

// ----------------------------------------------------------------------------
// Metafunction TextIterator_
// ----------------------------------------------------------------------------

template <typename TText, typename TPattern, typename TSpec>
struct TextIterator_
{
    typedef typename Iterator<TText>::Type  Type;
};

// ----------------------------------------------------------------------------
// Metafunction PatternIterator_
// ----------------------------------------------------------------------------

template <typename TText, typename TPattern, typename TSpec>
struct PatternIterator_
{
    typedef typename Iterator<TPattern const, Rooted>::Type  Type;
};

// ----------------------------------------------------------------------------
// Metafunction Score_
// ----------------------------------------------------------------------------

template <typename TSpec>
struct Score_
{
    typedef unsigned char   Type;
};

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function _textIterator()
// ----------------------------------------------------------------------------

template <typename TText, typename TPattern, typename TSpec>
SEQAN_HOST_DEVICE inline typename TextIterator_<TText, TPattern, TSpec>::Type &
_textIterator(Finder_<TText, TPattern, TSpec> & finder)
{
    return finder._textIt;
}

template <typename TText, typename TPattern, typename TSpec>
SEQAN_HOST_DEVICE inline typename TextIterator_<TText, TPattern, TSpec>::Type const &
_textIterator(Finder_<TText, TPattern, TSpec> const & finder)
{
    return finder._textIt;
}

// ----------------------------------------------------------------------------
// Function _patternIterator()
// ----------------------------------------------------------------------------

template <typename TText, typename TPattern, typename TSpec>
SEQAN_HOST_DEVICE inline typename PatternIterator_<TText, TPattern, TSpec>::Type &
_patternIterator(Finder_<TText, TPattern, TSpec> & finder)
{
    return finder._patternIt;
}

template <typename TText, typename TPattern, typename TSpec>
SEQAN_HOST_DEVICE inline typename PatternIterator_<TText, TPattern, TSpec>::Type const &
_patternIterator(Finder_<TText, TPattern, TSpec> const & finder)
{
    return finder._patternIt;
}

// ----------------------------------------------------------------------------
// Function _getScore()
// ----------------------------------------------------------------------------

template <typename TText, typename TPattern, typename TSpec>
SEQAN_HOST_DEVICE inline typename Score_<TSpec>::Type
_getScore(Finder_<TText, TPattern, TSpec> const & finder)
{
    return finder._score;
}

// ----------------------------------------------------------------------------
// Function _getScoreThreshold()
// ----------------------------------------------------------------------------

template <typename TText, typename TPattern, typename TSpec>
SEQAN_HOST_DEVICE inline typename Score_<TSpec>::Type
_getScoreThreshold(Finder_<TText, TPattern, TSpec> const & finder)
{
    return finder._scoreThreshold;
}

// ----------------------------------------------------------------------------
// Function _setScoreThreshold()
// ----------------------------------------------------------------------------

template <typename TText, typename TPattern, typename TSpec, typename TScore>
SEQAN_HOST_DEVICE inline void
_setScoreThreshold(Finder_<TText, TPattern, TSpec> & finder, TScore score)
{
    finder._scoreThreshold = score;
}

// ----------------------------------------------------------------------------
// Function clear()
// ----------------------------------------------------------------------------

template <typename TText, typename TIndexSpec, typename TPattern, typename TSpec>
SEQAN_HOST_DEVICE inline void
clear(Finder_<TText, TPattern, TSpec> & finder)
{
    // NOTE(esiragusa): if find wasn't called yet, _textIterator is uninitialized.
//    goBegin(_textIterator(finder));
    // NOTE(esiragusa): if find wasn't called yet, _patternIterator is uninitialized.
//    goBegin(_patternIterator(finder));
    finder._score = typename Score_<TSpec>::Type();
}


}

#endif  // #ifndef SEQAN_FIND_FIND_BASE_H_
