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
// Status: Testing.  Should work, some performance improvements possible,
// see TODOs.
// ==========================================================================
// Author: Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>
// ==========================================================================

#ifndef SEQAN_FIND_FIND_SIMPLE_H_
#define SEQAN_FIND_FIND_SIMPLE_H_

#include <algorithm>

namespace SEQAN_NAMESPACE_MAIN {

/*!
 * @class HammingSimplePattern
 * @extends Pattern
 * @headerfile <seqan/find.h>
 * @brief A brute force online searching algorithm for approximate string matching with hamming distance.
 *
 * @signature template <typename TNeedle>
 *            class Pattern<TNeedle, HammingSimple>;
 *
 * @tparam TNeedle The needle type. Types: String
 *
 * This specialization should only be used if no other is applicable or for verification purposes.
 */

struct HammingSimple_;
typedef Tag<HammingSimple_> HammingSimple;


template <typename TNeedle>
class Pattern<TNeedle, HammingSimple> {
public:
    // The holder for the needle.
    Holder<TNeedle> data_host;

    // The maximal distance.  Must be >= 0, i.e. -score.
    int maxDistance;

    // The current distance, >= 0, i.e. -current score.
    int distance;

    // TODO(holtgrew): This is extremely hacky, works only with Dna5.
    // Flags -- ..1xx activate, ..01 match pattern, ..10 match finder.
    unsigned matchNFlags;

    Pattern() : maxDistance(-1), distance(0), matchNFlags(0) {}

    template <typename TNeedle2>
    Pattern(const TNeedle2 &ndl, int k = -1) : maxDistance(-1), distance(0), matchNFlags(0) {
        SEQAN_CHECKPOINT;
        setHost(*this, ndl, k);
    }
};


template <typename TNeedle>
inline void _patternMatchNOfPattern(Pattern<TNeedle, HammingSimple> & pattern, bool matchN)
{
    if (matchN)
        pattern.matchNFlags |= 1;  // |= 01b
    else
        pattern.matchNFlags &= 2;  // &= 10b
    pattern.matchNFlags |= 4;
}


template <typename TNeedle>
inline void _patternMatchNOfFinder(Pattern<TNeedle, HammingSimple> & pattern, bool matchN)
{
    if (matchN)
        pattern.matchNFlags |= 2;  // |= 10b
    else
        pattern.matchNFlags &= 1;  // &= 01b
    pattern.matchNFlags |= 4;
}


template <typename TNeedle, typename TNeedle2>
void setHost (Pattern<TNeedle, HammingSimple> & me,
              const TNeedle2 & needle, int k) {
    SEQAN_CHECKPOINT;

    SEQAN_ASSERT_NOT(empty(needle));
    SEQAN_ASSERT_LEQ_MSG(k, 0, "Are you confusing distances and scores?");

    setValue(me.data_host, needle);
    me.maxDistance = -k;
}


template <typename TNeedle, typename TNeedle2>
void
setHost(Pattern<TNeedle, HammingSimple> &horsp, TNeedle2 &ndl, int k) {
    SEQAN_CHECKPOINT;
    setHost(horsp, reinterpret_cast<const TNeedle2&>(ndl), k);
}


template <typename TNeedle>
inline void _finderInit(Pattern<TNeedle, HammingSimple> & me) {
    SEQAN_CHECKPOINT;
    (void) me;  // Suppress unused variable warning.
}


template <typename TNeedle>
inline int score(const Pattern<TNeedle, HammingSimple> &me) {
    SEQAN_CHECKPOINT;
    return -me.distance;
}


template <typename TNeedle>
inline int getScore(const Pattern<TNeedle, HammingSimple> &me) {
    SEQAN_CHECKPOINT;
    return -me.distance;
}


template <typename TNeedle>
inline void setScoreLimit(Pattern<TNeedle, HammingSimple> & me, int _limit) {
    SEQAN_CHECKPOINT;
    SEQAN_ASSERT_LEQ(_limit, 0);
    me.maxDistance = -_limit;
}


template <typename TAlphabet, typename TNeedle>
inline bool _findHammingSimpleCharsEqual(TAlphabet const & a, TAlphabet const & b, Pattern<TNeedle, HammingSimple> const &) {
    return a == b;
}


template <typename TNeedle>
inline bool _findHammingSimpleCharsEqual(Dna5 const & ndlChar, Dna5 const & hstckChar, Pattern<TNeedle, HammingSimple> const & pattern) {
    if (ndlChar == Dna5('N') && pattern.matchNFlags & 1) {
        return true;
    } else if (hstckChar == Dna5('N') && pattern.matchNFlags & 2) {
        return true;
    } else {
        return ndlChar == hstckChar;
    }
}


template <typename TNeedle>
inline bool _findHammingSimpleCharsEqual(Dna5Q const & a, Dna5 const & b, Pattern<TNeedle, HammingSimple> const & pattern) {
    return _findHammingSimpleCharsEqual(convert<Dna5>(a), b, pattern);
}


template <typename TNeedle>
inline bool _findHammingSimpleCharsEqual(Dna5 const & a, Dna5Q const & b, Pattern<TNeedle, HammingSimple> const & pattern) {
    return _findHammingSimpleCharsEqual(a, convert<Dna5>(b), pattern);
}


template <typename TNeedle>
inline bool _findHammingSimpleCharsEqual(Dna5Q const & a, Dna5Q const & b, Pattern<TNeedle, HammingSimple> const & pattern) {
    return _findHammingSimpleCharsEqual(convert<Dna5>(a), convert<Dna5>(b), pattern);
}


template <typename TFinder, typename TNeedle>
inline bool find(TFinder &finder,
                 Pattern<TNeedle, HammingSimple> &me,
                 int minScore) {
    SEQAN_CHECKPOINT;

    typedef typename Haystack<TFinder>::Type THaystack;
    typedef typename Size<THaystack>::Type TSize;

    // Shortcuts to haystack and needle.
    const THaystack &hstk = haystack(finder);
    const TNeedle &ndl = needle(me);

    // If the needle is longer than the haystack then we cannot find anything.
    if (length(hstk) < length(ndl))
        return false;

    // Initialize or advance finder, depending whether it has been
    // initialized before.
    if (empty(finder)) {
        _finderInit(me);
        _setFinderLength(finder, length(needle(me)));
        _finderSetNonEmpty(finder);
    } else {
        finder += 1;
    }

    // Check whether we are beyond the last possible match position.
    if (position(finder) > length(hstk) - length(ndl))
        return false;

    // TODO(holtgrew): Switch from indices to iterators to improve performance.

    // Perform a naive search for the needle in the haystack such that
    // the difference is <= -minScore.
    //
    // If a special behaviour is enabled for N then we use a different case.
    TSize i;
    if (!(me.matchNFlags & 6u)) {
        // No special behaviour for N.
        for (i = position(finder); i <= length(hstk) - length(ndl); ++i) {
            me.distance = 0;  // Reset mismatch count.
            for (TSize j = 0; j < length(ndl); ++j) {
                me.distance += (ndl[j] != hstk[i + j]);
                if (me.distance > -minScore)
                    break;
            }
            if (me.distance <= -minScore)
                break;
        }
    } else {
        // Special behaviour for N enabled.
        for (i = position(finder); i <= length(hstk) - length(ndl); ++i) {
            me.distance = 0;  // Reset mismatch count.
            for (TSize j = 0; j < length(ndl); ++j) {
                me.distance += !(_findHammingSimpleCharsEqual(ndl[j], hstk[i + j], me));
                if (me.distance > -minScore)
                    break;
            }
            if (me.distance <= -minScore)
                break;
        }
    }

    // Return false if we did not break out of the for-loop but it
    // stopped normally.
    if (i > length(hstk) - length(ndl))
        return false;

    _setFinderEnd(finder, i + length(ndl));
    setPosition(finder, beginPosition(finder));
    return true;
}

template <typename TFinder, typename TNeedle>
inline bool find(TFinder &finder,
                 Pattern<TNeedle, HammingSimple> &me)
{
    return find(finder, me, -me.maxDistance);
}
}  // namespace SEQAN_NAMESPACE_MAIN

#endif  // SEQAN_FIND_FIND_SIMPLE_H_
