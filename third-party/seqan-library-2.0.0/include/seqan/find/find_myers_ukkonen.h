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

#ifndef SEQAN_HEADER_FIND_MYERS_UKKONEN_H
#define SEQAN_HEADER_FIND_MYERS_UKKONEN_H

#include <seqan/misc/sse2.h>

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////
// MyersUkkonen
//////////////////////////////////////////////////////////////////////////////

/*!
 * @class MyersPattern
 * @extends Pattern
 * @headerfile <seqan/find.h>
 * @brief Provides fast approximate searching of one string in another using Myer's fast bit-parallel algorithm with
 *        application of the Ukkonen- trick.
 *
 * @signature template <typename TNeedle[, typename TSpec[, typename TFindBeginPatternSpec]]>
 *            class Pattern<TNeedle, Myers<TSpec, TFindBeginPatternSpec> >;
 *
 * @tparam TSpec   Specialization tag.  This is @link ApproximateFinderSearchTypeTags#FindInfix @endlink for
 *                 infix search or @link ApproximateFinderSearchTypeTags#FindPrefix @endlink for prefix search.
 *                 Defaults to @linkApproximateFinderSearchTypeTags#FindInfix @endlink.
 * @tparam TFindBeginPatternSpec
 *               Specialization of @link Pattern @endlink used to find the begin of matches.This must be a finder for
 *               prefix search, e.g. @link DPSearchPattern <tt>DPSearch&lt;TScore, FindPrefix&gt;</tt> @endlink or @link
 *               MyersPattern <tt>Myers&lt;FindPrefix&gt;</tt> @endlink. Specify <tt>void</tt> to suppress prefix
 *               searching. Default: @link DefaultFindBeginPatternSpec @endlink
 * @tparam TNeedle The needle type. Types: String
 *
 * The needle-length must be smaller than the highest number that can be stored in an unsigned int.
 */

template <typename TSpec = FindInfix,
          typename THasState = True,
          typename TFindBeginPatternSpec = typename DefaultFindBeginPatternSpec<EditDistanceScore, THasState>::Type>
struct Myers {};

struct NMatchesNone_;
struct NMatchesN_;
struct NMatchesAll_;

//FindInfix and FindPrefix are defined int find_base.h
template <typename TSpec, typename TFinderCharSetPolicy = NMatchesN_, typename TPatternCharSetPolicy = NMatchesN_>
struct AlignTextBanded; // search query in a parallelogram

// TODO(holtgrew): Really deprecated?
//deprecated shortcuts:

/*!
 * @typedef MyersUkkonen
 * @headerfile <seqan/find.h>
 * @brief Semi-global (query-global, text-local) pattern matching without
 *        findBegin() support.
 *
 * @signature typedef Myers<FindInfix, True, void> MyersUkkonen;
 *
 * @deprecated Use <tt>Myers&lt;FindInfix&gt;</tt> instead.
 */

typedef Myers<FindInfix, True, void> MyersUkkonen;

/*!
 * @typedef MyersUkkonenGlobal
 * @headerfile <seqan/find.h>
 * @brief Global (query-global, text-global) pattern matching without findBegin() support.
 *
 * @signature typedef Myers<FindInfix, True, void> MyersUkkonenGlobal;
 */

typedef Myers<FindPrefix, True, void> MyersUkkonenGlobal;

/*!
 * @typedef MyersUkkonenBanded
 * @headerfile <seqan/find.h>
 * @brief Global (query-global, text-local) pattern matching without findBegin() support.
 *
 * @signature Myers<AlignTextBanded<FindInfix, NMatchesN_, NMatchesN_>, True, void> MyersUkkonenBanded;
 */

/*!
 * @typedef MyersUkkonenGlobalBanded
 * @headerfile <seqan/find.h>
 * @brief Global (query-global, text-global) pattern matching without findBegin() support.
 *
 * @signature Myers<AlignTextBanded<FindPrefix, NMatchesN_, NMatchesN_>, True, void> MyersUkkonenGlobalBanded;
 */

typedef Myers<AlignTextBanded<FindInfix, NMatchesN_, NMatchesN_>, True, void> MyersUkkonenBanded;
typedef Myers<AlignTextBanded<FindPrefix, NMatchesN_, NMatchesN_>, True, void> MyersUkkonenGlobalBanded;


//____________________________________________________________________________
// bit 0 of the HP bit-vector
// 0 for begin-gap-free haystack
// 1 for global alignments of haystack

template <typename T>
struct MyersUkkonenHP0_ {
    enum { VALUE = 0 };
};

template <>
struct MyersUkkonenHP0_<FindPrefix> {
    enum { VALUE = 1 };
};

template <typename TValue>
struct MyersSmallAlphabet_:
    public Eval<ValueSize<TValue>::VALUE <= 8> {};


//////////////////////////////////////////////////////////////////////////////
//overwrite FindBegin_ to define host member if find begin is switched on

template <typename TNeedle, typename TSpec, typename THasState, typename TFindBeginPatternSpec, typename TFindBeginPatternSpec2>
struct FindBegin_< Pattern<TNeedle, Myers<TSpec, THasState, TFindBeginPatternSpec> >, TFindBeginPatternSpec2>
{
private:
    typedef Pattern<TNeedle, Myers<TSpec, THasState, TFindBeginPatternSpec> > TPattern;
    typedef typename FindBeginPattern<TPattern>::Type TFindBeginPattern;

public:
    TFindBeginPattern data_findBeginPattern;
//     Holder<TNeedle>    data_host;    //defines the
    typedef False HasHost;
};

template <typename TNeedle, typename TSpec, typename THasState, typename TFindBeginPatternSpec>
struct FindBegin_< Pattern<TNeedle, Myers<TSpec, THasState, TFindBeginPatternSpec> >, void>
{
    typedef False HasHost;
//need no findBegin if FindBeginPatternSpec is void
};


//////////////////////////////////////////////////////////////////////////////
// State Data
//////////////////////////////////////////////////////////////////////////////

// small state
template <typename TNeedle, typename TSpec>
struct MyersSmallState_
{
#ifdef SEQAN_SSE2_INT128
    typedef Sse2Int128 TWord;
#else
    typedef unsigned long TWord;
#endif

    TWord VP0;                    // VP[0] (saves one dereferentiation)
    TWord VN0;                    // VN[0]
    unsigned int errors;        // the current number of errors
    unsigned int maxErrors;        // the maximal number of errors allowed

    MyersSmallState_() : VP0(0), VN0(0), errors(0), maxErrors(0)
    {}
};

template <typename TNeedle, typename TSmallAlphabet>
struct MyersSmallStateBandedShift_ {};
template <typename TNeedle>
struct MyersSmallStateBandedShift_<TNeedle, False> {
    typedef typename Value<TNeedle>::Type TValue;
    unsigned short shift[ValueSize<TValue>::VALUE];
};
template <typename TNeedle, typename TSpec, typename TFinderCSP, typename TPatternCSP>
struct MyersSmallState_<TNeedle, AlignTextBanded<TSpec, TFinderCSP, TPatternCSP> >:
    public MyersSmallStateBandedShift_<TNeedle, typename MyersSmallAlphabet_<typename Value<TNeedle>::Type>::Type>
{
#ifdef SEQAN_SSE2_INT128
    typedef Sse2Int128 TWord;
#else
#ifdef SEQAN_DEBUG_MYERSBITVECTOR
    typedef unsigned char TWord;
#else
    typedef unsigned long TWord;
#endif
#endif
    typedef typename Value<TNeedle>::Type TValue;

    TWord bitMasks[ValueSize<TValue>::VALUE];
    TWord VP0;                    // VP[0] (saves one dereferentiation)
    TWord VN0;                    // VN[0]
    unsigned short errors;      // the current number of errors
    unsigned short maxErrors;   // the maximal number of errors allowed
    unsigned short leftClip;    // clip that many characters from the text begin
//    unsigned short rightClip;   // stop alignment that many characters before the end   <<<< currently unused (autom. determined)

#ifdef SEQAN_DEBUG_MYERSBITVECTOR
    String<int> DPMat;
#endif
    MyersSmallState_() : bitMasks(), VP0(0), VN0(0), errors(0), maxErrors(0), leftClip(0) {}
};

// large state
template <typename TNeedle, typename TSpec>
struct MyersLargeState_
{
#ifdef SEQAN_SSE2_INT128
    typedef Sse2Int128 TWord;
#else
    typedef unsigned long TWord;
#endif
    unsigned lastBlock;            // the block containing the last active cell
    String<TWord> VP;
    String<TWord> VN;
    TWord scoreMask;            // the mask with a bit set at the position of the last active cell

    MyersLargeState_() : lastBlock(0), scoreMask(0) {}
};
template <typename TNeedle, typename TSpec, typename TFinderCSP, typename TPatternCSP>
struct MyersLargeState_<TNeedle, AlignTextBanded<TSpec, TFinderCSP, TPatternCSP> >
{
#ifdef SEQAN_SSE2_INT128
    typedef Sse2Int128 TWord;
#else
#ifdef SEQAN_DEBUG_MYERSBITVECTOR
    typedef unsigned char TWord;
#else
    typedef unsigned long TWord;
#endif
#endif
    unsigned lastBlock;            // the block containing the last active cell
    unsigned blockCount;        // the number of blocks
    String<TWord> VP;
    String<TWord> VN;
};

// TODO: should go elsewhere
// template <typename TNeedle, typename TSpec>
// class Pattern{};
// TODO: should go elsewhere
template <typename TNeedle, typename TSpec>
class PatternState_{};


template <typename TNeedle, typename TSpec, typename TFindBeginPatternSpec>
class PatternState_<TNeedle, Myers<TSpec, False, TFindBeginPatternSpec> > {};

template <typename TNeedle, typename TSpec, typename TFindBeginPatternSpec>
class PatternState_<TNeedle, Myers<TSpec, True, TFindBeginPatternSpec> >:
    public MyersSmallState_<TNeedle, TSpec>
{
public:
    typedef MyersSmallState_<TNeedle, TSpec>    TSmallState;
    typedef MyersLargeState_<TNeedle, TSpec>    TLargeState;
    typedef typename TSmallState::TWord         TWord;

    enum { MACHINE_WORD_SIZE = sizeof(TWord) * 8 };

    TLargeState *largeState;

//____________________________________________________________________________

    PatternState_():
        largeState(NULL) {}

    PatternState_(PatternState_ const & other):
        TSmallState(other),
        largeState(NULL)
    {
        if (other.largeState)
            largeState = new TLargeState(*other.largeState);
    }

    ~PatternState_()
    {
        delete largeState;
    }

    PatternState_ &
    operator = (PatternState_ const & other)
    {
        TSmallState::operator=(other);
        if (other.largeState)
        {
            if (largeState == NULL)
                largeState = new TLargeState;
            (*largeState) = *(other.largeState);
        } else {
            delete largeState;
            largeState = NULL;
        }
        return *this;
    }
};


//////////////////////////////////////////////////////////////////////////////
// Pattern Data
//////////////////////////////////////////////////////////////////////////////

template <typename TNeedle, typename TSpec>
struct MyersSmallPattern_
{
#ifdef SEQAN_SSE2_INT128
    typedef Sse2Int128 TWord;
#else
    typedef unsigned long TWord;
#endif

    String<TWord> bitMasks;        // encode the needle with bitmasks for each alphabet character
    unsigned needleSize;        // needle size

    MyersSmallPattern_():
        needleSize(0) {}
};
template <typename TNeedle, typename TSpec, typename TFinderCSP, typename TPatternCSP>
struct MyersSmallPattern_<TNeedle, AlignTextBanded<TSpec, TFinderCSP, TPatternCSP> >
{
#ifdef SEQAN_SSE2_INT128
    typedef Sse2Int128 TWord;
#else
#ifdef SEQAN_DEBUG_MYERSBITVECTOR
    typedef unsigned char TWord;
#else
    typedef unsigned long TWord;
#endif
#endif

    Holder<TNeedle> data_host;  // needle holder (the banded version needs no preprocessed bitmasks)
};

// large basic pattern
template <typename TNeedle, typename TSpec>
struct MyersLargePattern_
{
#ifdef SEQAN_SSE2_INT128
    typedef Sse2Int128 TWord;
#else
    typedef unsigned long TWord;
#endif

    unsigned blockCount;        // the number of blocks
    TWord finalScoreMask;        // a mask with a bit set on the position of the last row
};
template <typename TNeedle, typename TSpec, typename TFinderCSP, typename TPatternCSP>
struct MyersLargePattern_<TNeedle, AlignTextBanded<TSpec, TFinderCSP, TPatternCSP> > {};


template <typename TNeedle, typename TSpec, typename THasState, typename TFindBeginPatternSpec>
class Pattern<TNeedle, Myers<TSpec, THasState, TFindBeginPatternSpec> >:
    public MyersSmallPattern_<TNeedle, TSpec>,
    public FindBegin_<Pattern<TNeedle, Myers<TSpec, THasState, TFindBeginPatternSpec> > >,
    public PatternState_<TNeedle, Myers<TSpec, THasState, TFindBeginPatternSpec> >
{

public:
    typedef MyersSmallPattern_<TNeedle, TSpec>      TSmallPattern;
    typedef MyersLargePattern_<TNeedle, TSpec>      TLargePattern;
    typedef typename TSmallPattern::TWord           TWord;

    enum { MACHINE_WORD_SIZE = sizeof(TWord) * 8 };

    typedef PatternState_<TNeedle, Myers<TSpec, THasState, TFindBeginPatternSpec> > TPatternState;

    TLargePattern *largePattern;    // extra preprocessing info for large patterns

//____________________________________________________________________________

    Pattern():
        largePattern(NULL) {}

    Pattern(int _limit):
        largePattern(NULL)
    {
        setScoreLimit(*this, _limit);
    }

    Pattern(Pattern const & other) :
        TSmallPattern(other),
        TPatternState(other),
        largePattern(NULL)
    {
        if (other.largePattern)
            largePattern = new TLargePattern(*other.largePattern);
    }

    template <typename TNeedle2>
    Pattern(TNeedle2 const & ndl, int _limit = -1):
        largePattern(NULL)
    {
        setScoreLimit(*this, _limit);
        setHost(*this, ndl);
    }

    ~Pattern()
    {
        delete largePattern;
    }

    Pattern &
    operator = (Pattern const & other)
    {
        TSmallPattern::operator=(other);
        TPatternState::operator=(other);
        if (other.largePattern)
        {
            if (largePattern == NULL)
                largePattern = new TLargePattern;
            (*largePattern) = *(other.largePattern);
        } else {
            delete largePattern;
            largePattern = NULL;
        }
        return *this;
    }
};

//////////////////////////////////////////////////////////////////////////////
// Metafunctions
//////////////////////////////////////////////////////////////////////////////

template <typename TNeedle, typename TSpec, typename THasState, typename TFindBeginPatternSpec>
struct FindBeginPatternSpec <Pattern<TNeedle, Myers<TSpec, THasState, TFindBeginPatternSpec> > >
{
    typedef TFindBeginPatternSpec Type;
};
template <typename TNeedle, typename THasState, typename TFindBeginPatternSpec>
struct FindBeginPatternSpec <Pattern<TNeedle, Myers<FindPrefix, THasState, TFindBeginPatternSpec> > >
{// no find begin for FindPrefix
    typedef void Type;
};


template <typename TPattern>
struct PatternState {};

template <typename TNeedle, typename TSpec, typename THasState, typename TFindBeginPatternSpec>
struct PatternState<Pattern<TNeedle, Myers<TSpec, THasState, TFindBeginPatternSpec> > >
{
    typedef PatternState_<TNeedle, Myers<TSpec, True, TFindBeginPatternSpec> > Type;
};


//////////////////////////////////////////////////////////////////////////////
// Functions
//////////////////////////////////////////////////////////////////////////////


template <typename TNeedle, typename TSpec, typename THasState, typename TFindBeginPatternSpec, typename TNeedle2>
inline void _patternFirstInit(Pattern<TNeedle, Myers<TSpec, THasState, TFindBeginPatternSpec> > & pattern,
                              TNeedle2 & needle)
{
SEQAN_CHECKPOINT

    typedef typename Pattern<TNeedle, Myers<TSpec, TFindBeginPatternSpec> >::TWord TWord;
    typedef typename Value<TNeedle>::Type TValue;

    pattern.needleSize = length(needle);
    unsigned blockCount = (pattern.needleSize + pattern.MACHINE_WORD_SIZE - 1) / pattern.MACHINE_WORD_SIZE;

    if (blockCount > 1)
    {
        if (pattern.largePattern == NULL)
            pattern.largePattern = new MyersLargePattern_<TNeedle, TSpec>();

        pattern.largePattern->blockCount = blockCount;
        pattern.largePattern->finalScoreMask = (TWord)1 << ((pattern.needleSize + pattern.MACHINE_WORD_SIZE - 1) % pattern.MACHINE_WORD_SIZE);
    }
    else
    {
        delete pattern.largePattern;
        pattern.largePattern = NULL;
    }

    clear(pattern.bitMasks);
    resize(pattern.bitMasks, (ValueSize<TValue>::VALUE + 1) * blockCount, 0, Exact());

    // encoding the letters as bit-vectors
    for (unsigned j = 0; j < pattern.needleSize; j++)
        pattern.bitMasks[
            blockCount * ordValue(convert<typename Value<TNeedle>::Type>(getValue(needle, j)))
            + j / pattern.MACHINE_WORD_SIZE
        ] |= (TWord)1 << (j % pattern.MACHINE_WORD_SIZE);
        //pattern.bitMasks[pattern.blockCount * ordValue((CompareType< Value< TNeedle >::Type, Value< Container< THaystack >::Type >::Type >::Type) needle[j]) + j/pattern.MACHINE_WORD_SIZE] = pattern.bitMasks[pattern.blockCount * ordValue((CompareType< Value< TNeedle >::Type, Value< Container< THaystack >::Type >::Type >::Type) needle[j]) + j/MACHINE_WORD_SIZE] | ((TWord)1 << (j%MACHINE_WORD_SIZE));

    _findBeginInit(pattern, needle);
}


template <typename TNeedle, typename TSpec, typename TFinderCSP, typename TPatternCSP, typename THasState, typename TFindBeginPatternSpec, typename TNeedle2>
inline void _patternFirstInit(Pattern<TNeedle, Myers<AlignTextBanded<TSpec, TFinderCSP, TPatternCSP>, THasState, TFindBeginPatternSpec> > & pattern,
                              TNeedle2 & ndl)
{
    _findBeginInit(pattern, ndl);
}


template <typename TNeedle, typename TSpec, typename THasState, typename TFindBeginPatternSpec>
inline void _patternMatchNOfPatternImpl(Pattern<TNeedle, Myers<TSpec, THasState, TFindBeginPatternSpec> > & pattern,
                                        bool match)
{
    typedef typename Pattern<TNeedle, Myers<TSpec, THasState, TFindBeginPatternSpec> >::TWord TWord;
    unsigned blockCount = (pattern.largePattern == NULL)? 1: pattern.largePattern->blockCount;

    // letters are encoded as bit-vectors
    for (unsigned j = 0; j < pattern.needleSize; j++)
    {
        TWord bit = (TWord)1 << (j % pattern.MACHINE_WORD_SIZE);
        bool allNull = true;
        int idx = j / pattern.MACHINE_WORD_SIZE;

        for (int i = 0; i < 4; ++i, idx += blockCount)
            allNull &= (pattern.bitMasks[idx] & bit) == (TWord)0;

        if (allNull)
        {    // all bits are 0 => this letter must be 'N'
            if (match)
            {
                for (; idx >= 0; idx -= blockCount)
                    pattern.bitMasks[idx] |= bit;
            } else
            {
                for (; idx >= 0; idx -= blockCount)
                    pattern.bitMasks[idx] &= ~bit;
            }
        }
    }
}


template <typename TNeedle, typename TSpec, typename THasState, typename TFindBeginPatternSpec>
inline void
_patternMatchNOfPattern(Pattern<TNeedle, Myers<TSpec, THasState, TFindBeginPatternSpec> > & pattern, bool match)
{
    SEQAN_CHECKPOINT;
    _patternMatchNOfPatternImpl(pattern, match);
    _patternMatchNOfPatternImpl(pattern.data_findBeginPattern, match);
}


template <typename TNeedle, typename TSpec, typename THasState>
inline void
_patternMatchNOfPattern(Pattern<TNeedle, Myers<TSpec, THasState, void> > & pattern, bool match)
{
    SEQAN_CHECKPOINT;
    _patternMatchNOfPatternImpl(pattern, match);
}


template <typename TNeedle, typename TSpec, typename THasState, typename TFindBeginPatternSpec>
inline void
_patternMatchNOfFinderImpl(Pattern<TNeedle, Myers<TSpec, THasState, TFindBeginPatternSpec> > & pattern, bool match)
{
    SEQAN_CHECKPOINT;

    typedef typename Pattern<TNeedle, Myers<TSpec, THasState, TFindBeginPatternSpec> >::TWord TWord;
    unsigned blockCount = (pattern.largePattern == NULL)? 1: pattern.largePattern->blockCount;

    // letters are encoded as bit-vectors
    if (match)
    {
        for (unsigned j = 0; j < pattern.needleSize; j++)
            pattern.bitMasks[blockCount * 4 + j / pattern.MACHINE_WORD_SIZE] |= (TWord)1 << (j % pattern.MACHINE_WORD_SIZE);
    } else {
        for (unsigned j = 0; j < pattern.needleSize; j++)
            pattern.bitMasks[blockCount * 4 + j / pattern.MACHINE_WORD_SIZE] &= ~((TWord)1 << (j % pattern.MACHINE_WORD_SIZE));
    }
}


template <typename TNeedle, typename TSpec, typename THasState, typename TFindBeginPatternSpec>
inline void
_patternMatchNOfFinder(Pattern<TNeedle, Myers<TSpec, THasState, TFindBeginPatternSpec> > & pattern, bool match)
{
    SEQAN_CHECKPOINT;
    _patternMatchNOfFinderImpl(pattern, match);
    _patternMatchNOfFinderImpl(pattern.data_findBeginPattern, match);
}


template <typename TNeedle, typename TSpec, typename THasState>
inline void
_patternMatchNOfFinder(Pattern<TNeedle, Myers<TSpec, THasState, void> > & pattern, bool match)
{
    SEQAN_CHECKPOINT;
    _patternMatchNOfFinderImpl(pattern, match);
}


// data_host is not used anymore, the needle can be reconstructed from the bitmasks
template <typename TNeedle, typename TSpec, typename THasState, typename TFindBeginPatternSpec, typename TNeedle2>
inline void
_myersSetHost(Pattern<TNeedle, Myers<TSpec, THasState, TFindBeginPatternSpec> > &, TNeedle2 const &)
{
}

template <typename TNeedle, typename TSpec, typename TFinderCSP, typename TPatternCSP, typename THasState, typename TFindBeginPatternSpec, typename TNeedle2>
inline void
_myersSetHost(Pattern<TNeedle, Myers<AlignTextBanded<TSpec, TFinderCSP, TPatternCSP>, THasState, TFindBeginPatternSpec> > & pattern, TNeedle2 const & ndl)
{
    setValue(pattern.data_host, ndl);
}

template <typename TNeedle, typename TSpec, typename THasState, typename TFindBeginPatternSpec, typename TNeedle2>
inline void
setHost(Pattern<TNeedle, Myers<TSpec, THasState, TFindBeginPatternSpec> > & pattern, TNeedle2 & ndl)
{
    _myersSetHost(pattern, ndl);
    _patternFirstInit(pattern, ndl);
}


template <typename TNeedle, typename TSpec, typename THasState, typename TFindBeginPatternSpec, typename TNeedle2>
inline void
setHost(Pattern<TNeedle, Myers<TSpec, THasState, TFindBeginPatternSpec> > & pattern, TNeedle2 const & ndl)
{
    _myersSetHost(pattern, ndl);
    _patternFirstInit(pattern, ndl);
}


//____________________________________________________________________________


template <typename TNeedle, typename TSpec, typename TFinderCSP, typename TPatternCSP, typename THasState, typename TFindBeginPatternSpec>
inline typename Host<Pattern<TNeedle, Myers<AlignTextBanded<TSpec, TFinderCSP, TPatternCSP>, THasState, TFindBeginPatternSpec> > >::Type &
host(Pattern<TNeedle, Myers<AlignTextBanded<TSpec, TFinderCSP, TPatternCSP>, THasState, TFindBeginPatternSpec> > & pattern)
{
SEQAN_CHECKPOINT
    return value(pattern.data_host);
}


template <typename TNeedle, typename TSpec, typename TFinderCSP, typename TPatternCSP, typename THasState, typename TFindBeginPatternSpec>
inline typename Host<Pattern<TNeedle, Myers<AlignTextBanded<TSpec, TFinderCSP, TPatternCSP>, THasState, TFindBeginPatternSpec> > const>::Type &
host(Pattern<TNeedle, Myers<AlignTextBanded<TSpec, TFinderCSP, TPatternCSP>, THasState, TFindBeginPatternSpec> > const & pattern)
{
SEQAN_CHECKPOINT
    return value(pattern.data_host);
}


template <typename TNeedle, typename TSpec, typename THasState, typename TFindBeginPatternSpec>
inline TNeedle
host(Pattern<TNeedle, Myers<TSpec, THasState, TFindBeginPatternSpec> > const & pattern)
{
SEQAN_CHECKPOINT

    typedef typename Pattern<TNeedle, Myers<TSpec, TFindBeginPatternSpec> >::TWord TWord;
    typedef typename Value<TNeedle>::Type TValue;

    TNeedle temp;
    resize(temp, pattern.needleSize, Exact());

    unsigned blockCount = (pattern.needleSize + pattern.MACHINE_WORD_SIZE - 1) / pattern.MACHINE_WORD_SIZE;
    TValue v = TValue();
    for (unsigned i = 0; i < length(pattern.bitMasks); i += blockCount)
    {
        for (unsigned j = 0; j < pattern.needleSize; j++)
            if ((pattern.bitMasks[i + j / pattern.MACHINE_WORD_SIZE] & (TWord)1 << (j % pattern.MACHINE_WORD_SIZE)) != (TWord)0)
                temp[j] = v;
        ++v;
    }
    return temp;
}

template <typename TNeedle, typename TSpec, typename THasState, typename TFindBeginPatternSpec>
inline TNeedle
host(Pattern<TNeedle, Myers<TSpec, THasState, TFindBeginPatternSpec> > & pattern)
{
SEQAN_CHECKPOINT
    typedef Pattern<TNeedle, Myers<TSpec, THasState, TFindBeginPatternSpec> > TPattern;
    return host(const_cast<TPattern const &>(pattern));
}

//____________________________________________________________________________

template <typename TNeedle, typename TSpec, typename THasState, typename TFindBeginPatternSpec>
inline TNeedle
needle(Pattern<TNeedle, Myers<TSpec, THasState, TFindBeginPatternSpec> > const & pattern)
{
SEQAN_CHECKPOINT
    return host(pattern);
}

template <typename TNeedle, typename TSpec, typename THasState, typename TFindBeginPatternSpec>
inline TNeedle
needle(Pattern<TNeedle, Myers<TSpec, THasState, TFindBeginPatternSpec> > & pattern)
{
SEQAN_CHECKPOINT
    typedef Pattern<TNeedle, Myers<TSpec, THasState, TFindBeginPatternSpec> > TPattern;
    return host(const_cast<TPattern const &>(pattern));
}

//____________________________________________________________________________
/*!
 * @fn MyersPattern#scoreLimit
 * @headerfile <seqan/find.h>
 * @brief The minimal score a match must reach in approximate searching.
 *
 * @signature TScoreValue scoreLimit(pattern);
 *
 * @param[in] pattern The pattern to query.
 *
 * @return TScoreValue The score limit value.
 */

template <typename TNeedle, typename TSpec, typename TFindBeginPatternSpec>
inline int
scoreLimit(PatternState_<TNeedle, Myers<TSpec, True, TFindBeginPatternSpec> > const & state)
{
SEQAN_CHECKPOINT
    return - (int) state.maxErrors;
}

template <typename TNeedle, typename TSpec, typename TFindBeginPatternSpec>
inline int
scoreLimit(Pattern<TNeedle, Myers<TSpec, True, TFindBeginPatternSpec> > const & pattern)
{
SEQAN_CHECKPOINT
    return - (int) pattern.maxErrors;
}


//____________________________________________________________________________
/*!
 * @fn MyersPattern#setSoreLimit
 * @headerfile <seqan/find.h>
 * @brief Set the minimal score a match must reach in approximate serach.
 *
 * @signature void setScoreLimit(pattern, limit);
 *
 * @param[in,out] pattern The pattern to set the limit for.
 * @param[in]     limit   The limit score value to set.
 *
 * @return TScoreValue The score limit value.
 */

template <typename TNeedle, typename TSpec, typename TFindBeginPatternSpec, typename TScoreValue>
inline void
setScoreLimit(PatternState_<TNeedle, Myers<TSpec, True, TFindBeginPatternSpec> > & state,
           TScoreValue minScore)
{
SEQAN_CHECKPOINT
    // we need to convert the minimal score into a maximal penalty
    // that is why minScore is negated
    state.maxErrors = -minScore;
}
template <typename TNeedle, typename TSpec, typename TFindBeginPatternSpec, typename TScoreValue>
inline void
setScoreLimit(Pattern<TNeedle, Myers<TSpec, True, TFindBeginPatternSpec> > & pattern,
           TScoreValue minScore)
{
SEQAN_CHECKPOINT
    // we need to convert the minimal score into a maximal penalty
    // that is why minScore is negated
    pattern.maxErrors = -minScore;
}



//____________________________________________________________________________
/*!
 * @fn MyersPattern#getScore
 * @headerfile <seqan/find.h>
 * @brief Score of the last found match in approximate searching.
 *
 * @signature TScoreValue getScore(pattern);
 *
 * @param[in] pattern A myers pattern that can be used for approximate searching.
 *
 * @return TScoreValue The score of the last match found using <tt>pattern</tt>.  If no match was found then the value
 *                     is undefined.
 */

template <typename TNeedle, typename TSpec, typename TFindBeginPatternSpec>
inline int
getScore(PatternState_<TNeedle, Myers<TSpec, True, TFindBeginPatternSpec> > const & state)
{
    return -(int)state.errors;
}
template<typename TNeedle, typename TSpec, typename TFindBeginPatternSpec>
inline int
getScore(Pattern<TNeedle, Myers<TSpec, True, TFindBeginPatternSpec> > const & state)
{
    return -(int)state.errors;
}


template <typename TNeedle, typename TSpec, typename THasState, typename TFindBeginPatternSpec, typename TFinder>
inline bool
_patternInit(Pattern<TNeedle, Myers<TSpec, THasState, TFindBeginPatternSpec> > const & pattern,
             PatternState_<TNeedle, Myers<TSpec, True, TFindBeginPatternSpec> > & state,
             TFinder &)
{
SEQAN_CHECKPOINT
    typedef MyersLargeState_<TNeedle, TSpec> TLargeState;
    typedef typename TLargeState::TWord TWord;

    if (pattern.largePattern == NULL)
    {
        state.errors = pattern.needleSize;
        state.VP0 = ~(TWord)0;
        state.VN0 = 0;
        delete state.largeState;
        state.largeState = NULL;
    }
    else
    {
        if (state.largeState == NULL)
            state.largeState = new TLargeState;

        TLargeState &largeState = *state.largeState;
        // localMaxErrors either stores the maximal number of errors (me.maxErrors) or the needle size minus one.
        // It is used for the mask computation and setting the initial score (the minus one is there because of the Ukkonen trick).
        int localMaxErrors = _min(state.maxErrors, pattern.needleSize - 1);
        state.errors = localMaxErrors + 1;
        largeState.scoreMask = (TWord)1 << (localMaxErrors % pattern.MACHINE_WORD_SIZE);
        largeState.lastBlock = localMaxErrors / pattern.MACHINE_WORD_SIZE;
        if (largeState.lastBlock >= pattern.largePattern->blockCount)
            largeState.lastBlock = pattern.largePattern->blockCount - 1;

        clear(largeState.VP);
        resize(largeState.VP, pattern.largePattern->blockCount, ~(TWord)0, Exact());

        clear(largeState.VN);
        resize(largeState.VN, pattern.largePattern->blockCount, 0, Exact());
    }
    return true;
}

//____________________________________________________________________________
// bitmask operations - small alphabet

template <typename TNeedle, typename TSpec>
finline void
_myersPreInit(PatternState_<TNeedle, TSpec> &state, True)
{
    typedef typename Value<TNeedle>::Type TValue;
    for (unsigned i = 0; i < ValueSize<TValue>::VALUE; ++i)
        state.bitMasks[i] = 0;
}

template <typename TNeedle, typename TSpec>
finline void
_myersPostInit(PatternState_<TNeedle, TSpec> &state, True)
{
    typedef typename Value<TNeedle>::Type TValue;
    for (unsigned i = 0; i < ValueSize<TValue>::VALUE; ++i)
        state.bitMasks[i] >>= 1;
}

template <typename TNeedle, typename TSpec, typename TFinderCSP, typename TPatternCSP, typename TFindBeginPatternSpec, typename TValue, typename TShift>
finline void
_myersAdjustBitmask(PatternState_<TNeedle, Myers<AlignTextBanded<TSpec, TFinderCSP, TPatternCSP>, True, TFindBeginPatternSpec> > &state, TValue const value, TShift, True)
{
    typedef typename PatternState_<TNeedle, Myers<AlignTextBanded<TSpec, TFinderCSP, TPatternCSP>, True, TFindBeginPatternSpec> >::TWord TWord;

    // compiler will optimize that
    if (IsSameType<TPatternCSP, NMatchesAll_>::VALUE && value == unknownValue<TValue>())
    {
        for (unsigned i = 0; i < ValueSize<TValue>::VALUE; ++i)
            state.bitMasks[i] = (state.bitMasks[i] >> 1) | ((TWord)1 << (BitsPerValue<TWord>::VALUE - 1));
    }
    else
    {
        for (unsigned i = 0; i < ValueSize<TValue>::VALUE; ++i)
            state.bitMasks[i] >>= 1;
        if (!(IsSameType<TPatternCSP, NMatchesNone_>::VALUE && value == unknownValue<TValue>()))
            state.bitMasks[ordValue(value)] |= (TWord)1 << (BitsPerValue<TWord>::VALUE - 1);
    }

    if (IsSameType<TFinderCSP, NMatchesAll_>::VALUE)
        state.bitMasks[ordValue(unknownValue<TValue>())] |= (TWord)1 << (BitsPerValue<TWord>::VALUE - 1);
    if (IsSameType<TFinderCSP, NMatchesNone_>::VALUE)
        state.bitMasks[ordValue(unknownValue<TValue>())] &= ~((TWord)1 << (BitsPerValue<TWord>::VALUE - 1));
}

template <typename TNeedle, typename TSpec, typename TValue, typename TShift>
finline typename PatternState_<TNeedle, TSpec>::TWord
_myersGetBitmask(PatternState_<TNeedle, TSpec> &state, TValue const value, TShift, True)
{
    return state.bitMasks[ordValue(value)];
}



//____________________________________________________________________________
// bitmask operations - large alphabet

template <typename TNeedle, typename TSpec>
finline void
_myersPreInit(PatternState_<TNeedle, TSpec> &state, False)
{
    typedef typename Value<TNeedle>::Type TValue;
    memset(state.bitMasks, 0, (ValueSize<TValue>::VALUE + 1) * sizeof(state.bitMasks[0]));
    memset(state.shift, 0, ValueSize<TValue>::VALUE * sizeof(state.shift[0]));
}

template <typename TNeedle, typename TSpec>
finline void
_myersPostInit(PatternState_<TNeedle, TSpec> &, False)
{
}

template <typename TNeedle, typename TSpec, typename TFinderCSP, typename TPatternCSP, typename TFindBeginPatternSpec, typename TValue, typename TShift>
finline void
_myersAdjustBitmask(PatternState_<TNeedle, Myers<AlignTextBanded<TSpec, TFinderCSP, TPatternCSP>, True, TFindBeginPatternSpec> > &state, TValue const value, TShift const shift, False)
{
    typedef typename PatternState_<TNeedle, Myers<AlignTextBanded<TSpec, TFinderCSP, TPatternCSP>, True, TFindBeginPatternSpec> >::TWord TWord;

    if (IsSameType<TPatternCSP, NMatchesNone_>::VALUE && value == unknownValue<TValue>())
        return;

    unsigned ord = ordValue(value);
    unsigned short x = shift - state.shift[ord];
    if (x < BitsPerValue<TWord>::VALUE)
        state.bitMasks[ord] = (state.bitMasks[ord] >> x) | ((TWord)1 << (BitsPerValue<TWord>::VALUE - 1));
    else
        state.bitMasks[ord] = (TWord)1 << (BitsPerValue<TWord>::VALUE - 1);
    state.shift[ord] = shift;
}

template <typename TNeedle, typename TSpec, typename TFinderCSP, typename TPatternCSP, typename TFindBeginPatternSpec, typename TValue, typename TShift>
finline typename PatternState_<TNeedle, Myers<AlignTextBanded<TSpec, TFinderCSP, TPatternCSP>, True, TFindBeginPatternSpec> >::TWord
_myersGetBitmask(PatternState_<TNeedle, Myers<AlignTextBanded<TSpec, TFinderCSP, TPatternCSP>, True, TFindBeginPatternSpec> > &state, TValue const value, TShift const shift, False)
{
    typedef typename PatternState_<TNeedle, Myers<AlignTextBanded<TSpec, TFinderCSP, TPatternCSP>, True, TFindBeginPatternSpec> >::TWord TWord;

    if (IsSameType<TFinderCSP, NMatchesNone_>::VALUE && value == unknownValue<TValue>())
        return 0;

    if (IsSameType<TFinderCSP, NMatchesAll_>::VALUE && value == unknownValue<TValue>())
        return (shift < BitsPerValue<TWord>::VALUE)? -1 << shift: -1;

    unsigned ord = ordValue(value);
    TWord res;
    TShift x = shift - state.shift[ord];
    if (x < BitsPerValue<TWord>::VALUE)
        res = state.bitMasks[ord] >> x;
    else
        res = 0;

    if (IsSameType<TPatternCSP, NMatchesAll_>::VALUE)
    {
        ord = ordValue(unknownValue<TValue>());
        x = shift - state.shift[ord];
        if (x < BitsPerValue<TWord>::VALUE)
            res |= state.bitMasks[ord] >> x;
    }
    return res;
}


template <typename TFinder, typename TNeedle, typename TNeedle2, typename TSpec, typename TFinderCSP, typename TPatternCSP, typename TFindBeginPatternSpec>
inline bool
_patternInitSmallStateBanded(
    TFinder &finder,
    TNeedle2 const & needle,
    PatternState_<TNeedle, Myers<AlignTextBanded<TSpec, TFinderCSP, TPatternCSP>, True, TFindBeginPatternSpec> > & state)
{
    typedef Pattern<TNeedle, Myers<AlignTextBanded<TSpec, TFinderCSP, TPatternCSP>, True, TFindBeginPatternSpec> > TPattern;
    typedef typename TPattern::TWord TWord;
    typedef typename Iterator<TNeedle2 const, Standard>::Type TIter;
    typedef typename Value<TNeedle>::Type TValue;

#ifdef SEQAN_DEBUG_MYERSBITVECTOR
    int col = state.leftClip + 1;
#ifdef SEQAN_DEBUG_MYERSBITVECTOR_DUMP
    std::cerr << "     ";
    for (int i = length(needle); i != 0; --i)
        std::cerr << std::setw(5) << needle[i - 1];
    std::cerr << std::endl;
    std::cerr << "     ";
    for (int i = length(needle); i >= 0; --i)
        std::cerr << std::setw(5) << state.DPMat[i];
    std::cerr << std::endl;
#endif
#endif

    _myersPreInit(state, typename MyersSmallAlphabet_<TValue>::Type());

    typename Size<TNeedle>::Type const ndlLength = length(needle);

    // Initialize row 0 either with zeros or increasing numbers
    // This can be realized using the following DP pattern and
    // assuming character mismatches at rows -1, -2,...
    // Thus we initialize the bitmasks and VN with 0.
    // VP depends on global/local alignment
    //
    //  0  1  2  3  4   -2 -2 -2 -2 -2   (row -2)
    //  0  1  2  3  4   -1 -1 -1 -1 -1   (row -1)
    //  0  1  2  3  4    0  0  0  0  0   (row  0)
    //  1                1
    //        global           local
    //
    //  VP = 100...      VP = 111...
    //

    TWord VP = (MyersUkkonenHP0_<TSpec>::VALUE == 1)? (TWord)1 << ((int)BitsPerValue<TWord>::VALUE-1): maxValue<TWord>(); // HP[0]==1 <-> global, HP[0]==0 <-> local
    TWord VN = 0;

    // Errors are counted along the lowest diagonal and the
    // lowest row of the band.
    // The errors in the top-left corner are 0.
    //
    // 0 * * * *
    //   x * * * *
    //     x * * * *
    //       x x x x x
    //
    //       |-------|
    //     diagWidth + 1 = 5
    //
    // diagWidth = length(container(finder)) + state.leftClip + state.rightClip - length(needle)

    unsigned errors = 0;
    TIter ndlIter = begin(needle, Standard());
    TIter ndlEnd;

    // The errors along the diagonal can only increase or stay the same.
    // There is only the last row of length diagWidth where errors can decrease.
    // If errors exceeds cutOff it cannot reach maxErrors again.


    typename Size<TFinder>::Type const columns = length(container(finder)) + state.leftClip;
    unsigned cutOff = state.maxErrors;
    if (columns > ndlLength)
    {
        cutOff += columns - ndlLength;        // clipping case *0
        ndlEnd = end(needle, Standard());
    } else {
        errors += ndlLength - columns;
        ndlEnd = ndlIter + columns;            // clipping case *1
    }

//    std::cerr<<std::hex<<"\t  "<<std::setw(17)<<' '<<"\tVN"<<std::setw(17)<<VN<<"\tVP"<<std::setw(17)<<VP<<std::dec<<std::endl;

    unsigned short shift = 0;

    if (state.leftClip != 0)
    {
        //////////////////////////////////////////////////////////////////
        // PART 0: go down the parallelogram in a empty (clipped) area
        //////////////////////////////////////////////////////////////////

        errors += state.leftClip;
        if (errors > ndlLength) errors = ndlLength;
        if (errors > cutOff) return false;

    // leftClip = 2
    //   |-|
    //
    //   . . * * *
    //     . * * * *
    //       * * * * *
    //         * * * * .
    //           * * * . .
    //             * * . . .
    //
    //                 |---|
    //               rightClip = 3
    //
    // We divide the parallelogam into 3 sections:
    //
    //   A A A A
    //     A A A B
    //       A A B B
    //         A B B C
    //           B B C C
    //             B C C C
    //               C C C C
    //
    // Depending on where the clipping ends we identify 4 different clipping cases:
    //
    //     case 00            case 10            case 01            case 11
    //   . . * *            . . . .            . . * *            . . . .
    //     . * * *            . . . *            . * * *            . . . *
    //       * * * *            . . * *            * * * .            . . * .
    //         * * * *            . * * *            * * . .            . * . .
    //           * * * .            * * * .            * . . .            * . . .
    //             * * . .            * * . .            . . . .            . . . .
    //

        // adjust bitmasks (errors = number of needle chars to preprocess)
        for (; shift < errors; ++ndlIter, ++shift)
            _myersAdjustBitmask(state, getValue(ndlIter), shift, typename MyersSmallAlphabet_<TValue>::Type());

        // initialise left column with
        //
        //  0  1  2  3  4   -2 -2 -2 -2 -2
        //  0  1  2  3  4   -1 -1 -1 -1 -1
        //  0  1  2  3  4    0  0  0  0  0
        //  1                1
        //  2   global       2   local
        //  3                3
        //  4                4
        //
        //  VP = 111100...   VP = 111111...
        if (errors < (unsigned)BitsPerValue<TWord>::VALUE-1)
            VP |= ((TWord) -1) << ((unsigned)BitsPerValue<TWord>::VALUE-1 - errors);
        else
            VP = (TWord)-1;
    }

    for (; ndlIter != ndlEnd; ++ndlIter, goNext(finder), ++shift)
    {
        //////////////////////////////////////////////////////////////////
        // PART 1: go down the parallelogram
        //////////////////////////////////////////////////////////////////

        // adjust bitmask
        _myersAdjustBitmask(state, getValue(ndlIter), shift, typename MyersSmallAlphabet_<TValue>::Type());

        /////////////////////////
        // DIAGONAL MYERS CORE

        // VP/VN --> D0  (original Myers)
        TWord X = _myersGetBitmask(state, ordValue(*finder), shift, typename MyersSmallAlphabet_<TValue>::Type()) | VN;
        TWord D0 = ((VP + (X & VP)) ^ VP) | X;

        // adjust errors corresponding to rightmost bit of D0
        errors += (~D0 >> (BitsPerValue<TWord>::VALUE - 1)) & 1;
        if (errors > cutOff) return false;

        // D0 --> HP/HN  (original Myers)
        TWord HN = VP & D0;
        TWord HP = VN | ~(VP | D0);
    //    const int PADDING = sizeof(TWord)*2 + 1;
    //    std::cerr << std::hex;
    //    std::cerr << "\tD0"<<std::setw(PADDING)<<(__uint64)D0<<"\tHN"<<std::setw(PADDING)<<(__uint64)HN<<"\tHP"<<std::setw(PADDING)<<(__uint64)HP << std::endl;

        // moving register down corresponds to shifting HP/HN up (right shift)
        // HP/HN --> shift --> VP/VN (modified Myers)
        X = D0 >> 1;
        VN = X & HP;
        VP = HN | ~(X | HP);
    //    std::cerr << "\t  "<<std::setw(PADDING)<<' '<<"\tVN"<<std::setw(PADDING)<<(__uint64)VN<<"\tVP"<<std::setw(PADDING)<<(__uint64)VP << std::endl;
    //    std::cerr << std::dec;

#ifdef SEQAN_DEBUG_MYERSBITVECTOR
#ifdef SEQAN_DEBUG_MYERSBITVECTOR_DUMP
        std::cerr << "diag ";
#endif
        int val = errors;
        state.DPMat[(col-state.leftClip)*(length(needle)+1)+col] = val;
        for (int i = length(needle); i >=0; --i)
        {
            if (i > col)
            {
#ifdef SEQAN_DEBUG_MYERSBITVECTOR_DUMP
                std::cerr << "     ";
#endif
            } else
            {
                int shft = (int)BitsPerValue<TWord>::VALUE-1 - (col-i);
                if (shft >= 0)
                {
                    if (i < col)
                    {
                        TWord mask = (TWord)1 << (shft);
                        val -= ((VP & mask) != (TWord)0)? 1:0;
                        val += ((VN & mask) != (TWord)0)? 1:0;
                    }
                    state.DPMat[(col-state.leftClip)*(length(needle)+1)+i] = val;
#ifdef SEQAN_DEBUG_MYERSBITVECTOR_DUMP
                    std::cerr << std::setw(5) << val;
                } else
                {
                    std::cerr << "     ";
#endif
                }
            }
        }
#ifdef SEQAN_DEBUG_MYERSBITVECTOR_DUMP
        std::cerr << std::setw(5) << *finder;
        std::cerr << std::setw(5) << errors << std::endl;
#endif
        ++col;
#endif
    }
    state.VP0 = VP;
    state.VN0 = VN;
    state.errors = errors;
    _myersPostInit(state, typename MyersSmallAlphabet_<TValue>::Type());
    return true;
}


template <typename TFinder, typename TNeedle, typename TNeedle2, typename TSpec, typename TFinderCSP, typename TPatternCSP, typename TFindBeginPatternSpec>
bool _stateInit(
    TFinder &finder,
    TNeedle const & needle,
    PatternState_<TNeedle2, Myers<AlignTextBanded<TSpec, TFinderCSP, TPatternCSP>, True, TFindBeginPatternSpec> > & state)
{
    typedef PatternState_<TNeedle2, Myers<AlignTextBanded<TSpec, TFinderCSP, TPatternCSP>, True, TFindBeginPatternSpec> > TState;
    typedef typename TState::TLargeState TLargeState;

//    unsigned diagWidth = length(container(finder)) + state.leftClip - length(needle);
//    unsigned blockCount = diagWidth / state.MACHINE_WORD_SIZE + 1;
    unsigned blockCount = 1;

//    SEQAN_ASSERT_GEQ(length(container(finder)), length(needle));

#ifdef SEQAN_DEBUG_MYERSBITVECTOR
    clear(state.DPMat);
    resize(state.DPMat, (length(container(finder)) + 1) * (length(needle) + 1), -9);
    for (unsigned i = 0; i <= length(needle); ++i)
        state.DPMat[i] = i;
    for (unsigned i = 0; i <= length(container(finder)); ++i)
        state.DPMat[i * (length(needle) + 1)] = 0;
#endif

    if (blockCount <= 1)
    {
        delete state.largeState;
        state.largeState = NULL;
        return _patternInitSmallStateBanded(finder, needle, state);
    }
    else
    {
        // TODO: is that good here?
        if (state.largeState == NULL)
            state.largeState = new TLargeState;

        TLargeState &largeState = *state.largeState;
        largeState.blockCount = blockCount;

        clear(largeState.VP);
        resize(largeState.VP, blockCount, ~0, Exact());

        clear(largeState.VN);
        resize(largeState.VN, blockCount, 0, Exact());

        state.errors = 0;

        return true;
    }
}

template <typename TNeedle, typename TSpec, typename TFindBeginPatternSpec, typename TFinder>
inline bool
_patternInit(Pattern<TNeedle, Myers<TSpec, True, TFindBeginPatternSpec> > & pattern, TFinder & finder)
{
    SEQAN_CHECKPOINT;
    return _patternInit(pattern, pattern, finder);
}



//////////////////////////////////////////////////////////////////////////////
// Myers-Ukkonen for semi-global edit-distance-alignments
// (version for needles longer than one machineword)
//////////////////////////////////////////////////////////////////////////////


template <typename TFinder, typename TNeedle, typename TSpec, typename THasState, typename THasState2, typename TFindBeginPatternSpec, typename TSize>
inline bool _findMyersLargePatterns (TFinder & finder,
                                     Pattern<TNeedle, Myers<TSpec, THasState, TFindBeginPatternSpec> > const & pattern,
                                     PatternState_<TNeedle, Myers<TSpec, THasState2, TFindBeginPatternSpec> > & state,
                                     TSize haystack_length)
{
SEQAN_CHECKPOINT
    typedef MyersLargePattern_<TNeedle, TSpec> TLargePattern;
    typedef MyersLargeState_<TNeedle, TSpec> TLargeState;
    typedef typename TLargeState::TWord TWord;

    TWord X, D0, HN, HP, temp;
    TWord carryD0, carryHP, carryHN;
    unsigned shift, limit, currentBlock;

    TLargePattern &largePattern = *pattern.largePattern;
    TLargeState &largeState = *state.largeState;

    while (position(finder) < haystack_length)
    {
        carryD0 = carryHN = 0;
        carryHP = (int)MyersUkkonenHP0_<TSpec>::VALUE; // FIXME: replace Noting with TSpec

        // if the active cell is the last of it's block, one additional block has to be calculated
        limit = largeState.lastBlock + (unsigned)(largeState.scoreMask >> (pattern.MACHINE_WORD_SIZE - 1));

        if (limit == largePattern.blockCount)
            limit--;

        shift = largePattern.blockCount * ordValue((typename Value< TNeedle >::Type) *finder);

        // computing the necessary blocks, carries between blocks following one another are stored
        for (currentBlock = 0; currentBlock <= limit; currentBlock++)
        {
            X = pattern.bitMasks[shift + currentBlock] | largeState.VN[currentBlock];

            temp = largeState.VP[currentBlock] + (X & largeState.VP[currentBlock]) + carryD0;
            if (carryD0 != (TWord)0)
                carryD0 = temp <= largeState.VP[currentBlock];
            else
                carryD0 = temp < largeState.VP[currentBlock];

            D0 = (temp ^ largeState.VP[currentBlock]) | X;
            HN = largeState.VP[currentBlock] & D0;
            HP = largeState.VN[currentBlock] | ~(largeState.VP[currentBlock] | D0);

            X = (HP << 1) | carryHP;
            carryHP = HP >> (pattern.MACHINE_WORD_SIZE - 1);

            largeState.VN[currentBlock] = X & D0;

            temp = (HN << 1) | carryHN;
            carryHN = HN >> (pattern.MACHINE_WORD_SIZE - 1);

             largeState.VP[currentBlock] = temp | ~(X | D0);

            // if the current block is the one containing the last active cell
            // the new number of errors is computed
            if (currentBlock == largeState.lastBlock) {
                if ((HP & largeState.scoreMask) != (TWord)0)
                    state.errors++;
                else if ((HN & largeState.scoreMask) != (TWord)0)
                    state.errors--;
            }
        }

        // updating the last active cell
        while (state.errors > state.maxErrors) {
            if ((largeState.VP[largeState.lastBlock] & largeState.scoreMask) != (TWord)0)
                state.errors--;
            else if ((largeState.VN[largeState.lastBlock] & largeState.scoreMask) != (TWord)0)
                state.errors++;

            largeState.scoreMask >>= 1;
            if (largeState.scoreMask == (TWord)0)
            {
                largeState.lastBlock--;
                if (IsSameType<TSpec, FindPrefix>::VALUE && largeState.lastBlock == (unsigned)-1)
                    break;
                largeState.scoreMask = (TWord)1 << (pattern.MACHINE_WORD_SIZE - 1);
            }
        }

        if ((largeState.scoreMask == largePattern.finalScoreMask) && (largeState.lastBlock == largePattern.blockCount - 1))
        {
            _setFinderEnd(finder);
            if (IsSameType<TSpec, FindPrefix>::VALUE)
            {
                _setFinderLength(finder, endPosition(finder));
            }
            return true;
        }
        else {
            largeState.scoreMask <<= 1;
            if (!largeState.scoreMask) {
                largeState.scoreMask = 1;
                largeState.lastBlock++;
            }

            if ((largeState.VP[largeState.lastBlock] & largeState.scoreMask) != (TWord)0)
                state.errors++;
            else if ((largeState.VN[largeState.lastBlock] & largeState.scoreMask) != (TWord)0)
                state.errors--;
        }

//        SEQAN_ASSERT (state.errors >= 0);

        goNext(finder);
    }

    return false;
}


template <typename TFinder, typename TNeedle, typename TSpec, typename THasState, typename THasState2, typename TFindBeginPatternSpec, typename TSize>
inline bool
_findMyersSmallPatterns(
    TFinder & finder,
    Pattern<TNeedle, Myers<TSpec, THasState, TFindBeginPatternSpec> > const & pattern,
    PatternState_<TNeedle, Myers<TSpec, THasState2, TFindBeginPatternSpec> > & state,
    TSize haystack_length)
{
SEQAN_CHECKPOINT

    typedef typename Pattern<TNeedle, Myers<TSpec, THasState, TFindBeginPatternSpec> >::TWord TWord;

    TWord X, D0, HN, HP;
    TWord lastBit = (TWord)1 << (pattern.needleSize - 1);

    // computing the blocks
    while (position(finder) < haystack_length)
    {
        X = pattern.bitMasks[ordValue((typename Value<TNeedle>::Type) *finder)] | state.VN0;

        D0 = ((state.VP0 + (X & state.VP0)) ^ state.VP0) | X;
        HN = state.VP0 & D0;
        HP = state.VN0 | ~(state.VP0 | D0);
        X = (HP << 1) | (TWord)(int)MyersUkkonenHP0_<TSpec>::VALUE; // FIXME: replace Nothing by TSpec
        state.VN0 = X & D0;
        state.VP0 = (HN << 1) | ~(X | D0);

        if ((HP & lastBit) != (TWord)0)
            state.errors++;
        else if ((HN & lastBit) != (TWord)0)
            state.errors--;

        if (state.errors <= state.maxErrors)
        {
            _setFinderEnd(finder);
            if (IsSameType<TSpec, FindPrefix>::VALUE)
            {
                _setFinderLength(finder, endPosition(finder));
            }
            return true;
        }
        //
        // if (IsSameType<TSpec, FindPrefix>::VALUE)
        // {//limit haystack length during prefix search
        //
        // }

        goNext(finder);
    }

    return false;
}


//////////////////////////////////////////////////////////////////////////////
// Myers-Ukkonen as a banded alignment
// the band includes the main diagonal and the diagonals above
// the band width is (blockCount * MACHINE_WORD_SIZE)
//////////////////////////////////////////////////////////////////////////////

template <
    typename TFinder,
    typename TNeedle,
    typename TNeedle2,
    typename TSpec,
    typename TFinderCSP,
    typename TPatternCSP,
    typename TFindBeginPatternSpec,
    typename TDoPatternSearch
>
inline bool
_findMyersSmallPatternsBanded(
    TFinder & finder,
    TNeedle const & needle,
    PatternState_<TNeedle2, Myers<AlignTextBanded<TSpec, TFinderCSP, TPatternCSP>, True, TFindBeginPatternSpec> > & state,
    TDoPatternSearch const)
{
    typedef PatternState_<TNeedle, Myers<AlignTextBanded<TSpec, TFinderCSP, TPatternCSP>, True, TFindBeginPatternSpec> > TState;
    typedef typename TState::TWord TWord;
    typedef typename Value<TNeedle>::Type TValue;

    TWord VP = state.VP0;
    TWord VN = state.VN0;
    TWord errors = state.errors;
    TWord const maxErrors = state.maxErrors;
    unsigned short const shift = length(needle);

#ifdef SEQAN_DEBUG_MYERSBITVECTOR
    unsigned col = position(finder) + 1;
#endif

    for (; !atEnd(finder); goNext(finder))
    {
        // PART 2: go right

        // normal Myers
        TWord X = _myersGetBitmask(state, ordValue(*finder), shift, typename MyersSmallAlphabet_<TValue>::Type()) | VN;
        TWord D0 = ((VP + (X & VP)) ^ VP) | X;
        TWord HN = VP & D0;
        TWord HP = VN | ~(VP | D0);
    //    const int PADDING = sizeof(TWord)*2 + 1;
    //    std::cerr << std::hex;
    //    std::cerr << "\tD0"<<std::setw(PADDING)<<(__uint64)D0<<"\tHN"<<std::setw(PADDING)<<(__uint64)HN<<"\tHP"<<std::setw(PADDING)<<(__uint64)HP<<std::endl;
        X = (HP << 1) | 1;
        VN = X & D0;
        VP = (HN << 1) | ~(X | D0);
    //    std::cerr << "\t  "<<std::setw(PADDING)<<' '<<"\tVN"<<std::setw(PADDING)<<(__uint64)VN<<"\tVP"<<std::setw(PADDING)<<(__uint64)VP<<std::endl;
    //    std::cerr << std::dec;
        errors += (HP >> (BitsPerValue<TWord>::VALUE - 2)) & 1;
        errors -= (HN >> (BitsPerValue<TWord>::VALUE - 2)) & 1;

        // shift bitmasks and states
#ifdef SEQAN_DEBUG_MYERSBITVECTOR
#ifdef SEQAN_DEBUG_MYERSBITVECTOR_DUMP
        std::cerr << "horiz";
#endif
        int val = errors;
        state.DPMat[col*(length(needle)+1)+length(needle)] = val;
        for (int i = length(needle); i >= 0; --i)
        {
            int shft = (int)BitsPerValue<TWord>::VALUE-1 - (length(needle)-i);
            if (shft >= 0)
            {
                if (i < (int)length(needle))
                {
                    TWord mask = (TWord)1 << (shft);
                    val -= ((VP & mask) != (TWord)0)? 1:0;
                    val += ((VN & mask) != (TWord)0)? 1:0;
                }
                state.DPMat[col*(length(needle)+1)+i] = val;
#ifdef SEQAN_DEBUG_MYERSBITVECTOR_DUMP
                std::cerr << std::setw(5) << val;
            } else {
                std::cerr << "     ";
#endif
            }
        }
#ifdef SEQAN_DEBUG_MYERSBITVECTOR_DUMP
        std::cerr << std::setw(5) << *finder;
        std::cerr << std::setw(5) << errors << std::endl;
#endif
        ++col;
#endif

        if (TDoPatternSearch::VALUE)
        {
            // pattern search
            if (errors <= maxErrors)
            {
                state.VP0 = VP;
                state.VN0 = VN;
                state.errors = errors;
                _setFinderEnd(finder);
                if (IsSameType<TSpec, FindPrefix>::VALUE)
                {
                    _setFinderLength(finder, endPosition(finder));
                }
                return true;
            }
        }
        else
        {
            // edit distance
        }
    }

    if (!TDoPatternSearch::VALUE)
    {
        // edit distance
        state.errors = errors;
    }
    return false;
}

template <typename TSeq1, typename TSeq2>
inline unsigned
_computeEditDistanceBanded(
    TSeq1 const &seq1,
    TSeq2 const &seq2,
    unsigned maxErrors)
{
    PatternState_<TSeq2, MyersUkkonenGlobalBanded> state;
    state.maxErrors = maxErrors;
    state.leftClip = (length(seq2) - length(seq1) + maxErrors) / 2;
    typename Iterator<TSeq1 const, Rooted>::Type seq1Iter = begin(seq1, Rooted());

    if (!_patternInitSmallStateBanded(seq1Iter, seq2, state))
        return maxErrors + 1;
    _findMyersSmallPatternsBanded(seq1Iter, seq2, state, False());
    return state.errors;
}


//////////////////////////////////////////////////////////////////////////////
// find

template <typename TFinder, typename TNeedle, typename TNeedle2, typename TSpec, typename TFinderCSP, typename TPatternCSP, typename TFindBeginPatternSpec>
inline bool find (TFinder & finder,
                  TNeedle const & needle,
                  PatternState_<TNeedle2, Myers<AlignTextBanded<TSpec, TFinderCSP, TPatternCSP>, True, TFindBeginPatternSpec> > & state)
{
    if (empty(finder))
    {
        _finderSetNonEmpty(finder);
        if (!_stateInit(finder, needle, state))
        {
            goEnd(finder);
            return false;
        }

        if (state.errors <= state.maxErrors)
        {
            goPrevious(finder);
            _setFinderEnd(finder);
            if (IsSameType<TSpec, FindPrefix>::VALUE)
            {
                _setFinderLength(finder, endPosition(finder));
            }
            return true;
        }
        //TODO: adapt myers-ukkonnen to dynamically change maxErrors
    }
    else
    {
        if (atEnd(finder)) return false;
        goNext(finder);
    }

    // distinguish between the version for needles not longer than one machineword and the version for longer needles
    if (state.largeState == NULL)
        return _findMyersSmallPatternsBanded(finder, needle, state, True());
//    else
//        return _findMyersLargePatterns(finder, needle, state);
    return false;
}

// First two for AlignTextBanded
template <typename TFinder, typename TNeedle, typename TSpec, typename TFinderCSP, typename TPatternCSP, typename THasState, typename TFindBeginPatternSpec>
inline bool find (TFinder & finder,
                  Pattern<TNeedle, Myers<AlignTextBanded<TSpec, TFinderCSP, TPatternCSP>, THasState, TFindBeginPatternSpec> > const & pattern,
                  PatternState_<TNeedle, Myers<AlignTextBanded<TSpec, TFinderCSP, TPatternCSP>, True, TFindBeginPatternSpec> > & state)
{
SEQAN_CHECKPOINT
    return find(finder, host(pattern), state);
}

template <typename TFinder, typename TNeedle, typename TSpec, typename THasState, typename TFindBeginPatternSpec>
inline bool find (TFinder & finder,
                  Pattern<TNeedle, Myers<TSpec, THasState, TFindBeginPatternSpec> > const & pattern,
                  PatternState_<TNeedle, Myers<TSpec, True, TFindBeginPatternSpec> > & state)
{
SEQAN_CHECKPOINT
    typedef typename Haystack<TFinder>::Type THaystack;
    typedef typename Size<THaystack>::Type TSize;

    TSize prefix_begin_position; //for prefix search: the position where the prefix begins

    if (empty(finder))
    {
        _patternInit(pattern, state, finder);
        _finderSetNonEmpty(finder);

        prefix_begin_position = position(finder);

        //TODO: adapt myers-ukkonnen to dynamically change maxErrors
    }
    else
    {
        if (atEnd(finder)) return false;
        goNext(finder);

        prefix_begin_position = beginPosition(finder);
    }

    TSize haystack_length = length(container(finder));
    // limit search width for prefix search
    if (IsSameType<TSpec, FindPrefix>::VALUE)
    {
        TSize maxlen = prefix_begin_position + pattern.needleSize - scoreLimit(state) + 1;
        if (haystack_length > maxlen)
            haystack_length = maxlen;
    }

    // distinguish between the version for needles not longer than one machineword and the version for longer needles
    if (pattern.largePattern == NULL)
        return _findMyersSmallPatterns(finder, pattern, state, haystack_length);
    else
        return _findMyersLargePatterns(finder, pattern, state, haystack_length);
}


template <typename TFinder, typename TNeedle, typename TSpec, typename TFindBeginPatternSpec>
inline bool find (TFinder & finder,
                  Pattern<TNeedle, Myers<TSpec, True, TFindBeginPatternSpec> > & pattern)
{
    typedef typename Pattern<TNeedle, Myers<TSpec, True, TFindBeginPatternSpec> >::TPatternState TPatternState;
    return find(finder, pattern, static_cast<TPatternState&>(pattern));
}


template <typename TFinder, typename TNeedle, typename TSpec, typename THasState, typename TFindBeginPatternSpec>
inline bool find (TFinder & finder,
                  Pattern<TNeedle, Myers<TSpec, THasState, TFindBeginPatternSpec> > const & pattern,
                  PatternState_<TNeedle, Myers<TSpec, True, TFindBeginPatternSpec> > & state,
                  int const minScore)
{
    setScoreLimit(state, minScore);
    return find(finder, pattern, state);
}

template <typename TFinder, typename TNeedle, typename TNeedle2, typename TSpec, typename TFindBeginPatternSpec>
inline bool find (TFinder & finder,
                  TNeedle const & needle,
                  PatternState_<TNeedle2, Myers<TSpec, True, TFindBeginPatternSpec> > & state,
                  int const minScore)
{
    setScoreLimit(state, minScore);
    return find(finder, needle, state);
}

template <typename TFinder, typename TNeedle, typename TSpec, typename TFindBeginPatternSpec>
inline bool find (TFinder & finder,
                  Pattern<TNeedle, Myers<TSpec, True, TFindBeginPatternSpec> > & pattern,
                  int const minScore)
{
    return find(finder, pattern, pattern, minScore); //static cast
}

}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...

