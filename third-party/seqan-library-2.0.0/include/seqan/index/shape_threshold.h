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

#ifndef SEQAN_HEADER_SHAPE_THRESHOLD_H
#define SEQAN_HEADER_SHAPE_THRESHOLD_H

namespace SEQAN_NAMESPACE_MAIN
{

struct ThreshQGramLemma_;
struct ThreshExact_;
struct ThreshHeuristic_;

typedef Tag<ThreshQGramLemma_> const    ThreshQGramLemma;
typedef Tag<ThreshHeuristic_> const        ThreshHeuristic;
typedef Tag<ThreshExact_> const            ThreshExact;


//////////////////////////////////////////////////////////////////////////////
// q-gram lemma
//
// - exact for ungapped shapes or errors <= 1
// - lower bound gapped shapes
//////////////////////////////////////////////////////////////////////////////

template <typename TShape, typename TPatternLength, typename TErrors, typename TDistance>
inline int qgramThreshold(TShape const & shape, TPatternLength patternLength, TErrors errors, TDistance const, ThreshQGramLemma const)
{
    int t = (int)patternLength - (int)length(shape) + 1 - (int)errors * (int)weight(shape);
    return (t > 0)? t: 0;
}


//////////////////////////////////////////////////////////////////////////////
// q-gram heuristic
//
// - exact for errors <= 1
// - upper bound
//////////////////////////////////////////////////////////////////////////////

template <typename TShape, typename TPatternSize, typename TErrors, typename TDistance>
int qgramThreshold(TShape const & shape, TPatternSize patternLength, TErrors errors, TDistance const, ThreshHeuristic const)
{
    String<unsigned char> coverage;
    String<bool> preserved;
    String<unsigned> ones;
    CharString bitString;

    // initialize coverage map and bitmap of preserved q-grams
    resize(preserved, patternLength - length(shape) + 1, true);
    resize(coverage, patternLength, 0);

    shapeToString(bitString, shape);
    for (unsigned i = 0; i < length(bitString); ++i)
        if (bitString[i] == '1')
        {
            appendValue(ones, i);
            for (unsigned j = 0; j < length(preserved); ++j)
                ++coverage[i + j];
        }

    // greedily destroy a maximum number of q-grams
    for (; errors > 0; --errors)
    {
        // find position that destroys a maximum number of q-grams
        unsigned maxCoverage = 0;
        unsigned maxCoveragePos = 0;
        for (unsigned i = 0; i < length(coverage); ++i)
            if (maxCoverage < coverage[i])
            {
                maxCoverage = coverage[i];
                maxCoveragePos = i;
            }

        // destroy q-grams
        for (unsigned k = 0; k < length(ones); ++k)
            if (ones[k] <= maxCoveragePos)
            {
                unsigned startPos = maxCoveragePos - ones[k];
                if (startPos < length(preserved) && preserved[startPos])
                {
                    preserved[startPos] = false;
                    for (unsigned l = 0; l < length(ones); ++l)
                        --coverage[startPos + ones[l]];
                }
            }
    }

    unsigned thresh = 0;
    for (unsigned i = 0; i < length(preserved); ++i)
        if (preserved[i])
            ++thresh;

    return thresh;
}





//____________________________________________________________________________
// Extensions to SeqAn

    struct ErrorAlphabet_ {};
    typedef SimpleType<unsigned char, ErrorAlphabet_> ErrorAlphabet;

    template <> struct ValueSize< ErrorAlphabet >    { enum { VALUE = 4 }; };
    template <> struct BitsPerValue< ErrorAlphabet > { enum { VALUE = 2 }; };

    template <typename T = void>
    struct TranslateTableErrorToChar_
    {
        static char const VALUE[4];
    };
    template <typename T>
    char const TranslateTableErrorToChar_<T>::VALUE[4] = {'.', 'M', 'I', 'D'};

    inline void assign(char & c_target,
                       ErrorAlphabet const & source)
    {
    SEQAN_CHECKPOINT
        c_target = TranslateTableErrorToChar_<>::VALUE[source.value];
    }


    struct ErrorPackedString;

    template <typename TValue>
    struct Host<String<TValue, Packed<ErrorPackedString> > >
    {
        typedef String<__int64, Array<1> > Type;
    };

    template <typename TValue>
    struct Host<String<TValue, Packed<ErrorPackedString> > const >
    {
        typedef String<__int64, Array<1> > const Type;
    };


//____________________________________________________________________________

    enum ErrorType {
        SEQAN_MATCH    = 0,
        SEQAN_MISMATCH = 1,
        SEQAN_INSERT   = 2,
        SEQAN_DELETE   = 3
    };

    template <typename TDistance>
    struct ErrorTypes {
        enum { VALUE = 4 };
    };

    template <>
    struct ErrorTypes<HammingDistance> {
        enum { VALUE = 2 };
    };

    // descriptor of the modification pattern
    // in the recursion it modifies the last q-gram of a read sequence
    template <typename TDistance, typename TFloat>
    struct SensitivityDPState_
    {
        enum { TRANSITIONS = ErrorTypes<TDistance>::VALUE };
        TFloat prob;                // probability of this state
        int transition[ErrorTypes<TDistance>::VALUE];    // returns previous state
        unsigned char len;            // length of this pattern (shapeSpan-errors <= this value <= shapeSpan+errors)
        unsigned char errors:4;        // errors in this state
        bool skipFirst:1;            // skip this pattern if it is the first
        bool skipLast:1;            // skip this pattern if it is the last
        bool intermediate:1;        // this is an intermediate result (beginning with INSERT)
        bool qgramHit:1;            // is this a q-gram hit? (result of the former delta function)
    }
#ifndef PLATFORM_WINDOWS
    __attribute__((packed))
#endif
    ;

    // descriptor of the modification pattern
    // in the recursion it modifies the last q-gram of a read sequence
    template <typename TDistance>
    struct ThreshDPState_
    {
        enum { TRANSITIONS = ErrorTypes<TDistance>::VALUE };
        int transition[ErrorTypes<TDistance>::VALUE];    // returns previous state
        unsigned char len;            // length of this pattern (shapeSpan-errors <= this value <= shapeSpan+errors)
        unsigned char errors:4;        // errors in this state
        bool skipFirst:1;            // skip this pattern if it is the first
        bool skipLast:1;            // skip this pattern if it is the last
        bool intermediate:1;        // this is an intermediate result (beginning with INSERT)
        bool qgramHit:1;            // is this a q-gram hit? (result of the former delta function)
    }
#ifndef PLATFORM_WINDOWS
    __attribute__((packed))
#endif
    ;

#if defined(PLATFORM_WINDOWS) | defined(PLATFORM_CUDA)

    template<typename TValue>
    inline bool isNan(TValue value)
    {
        return value != value;
    }

    template<typename TValue>
    inline bool isInf(TValue value)
    {
        return value == log(0.0);
    }

#else

    template<typename TValue>
    inline bool isNan(TValue value)
    {
        return std::isnan(value);
    }

    template<typename TValue>
    inline bool isInf(TValue value)
    {
        return std::isinf(value);
    }

#endif


    template <typename TValue>
    inline long double
    _transform(TValue a)
    {
#ifdef USE_LOGVALUES
        return log(a);
#else
        return a;
#endif
    }

    template <typename TValue>
    inline long double
    _transformBack(TValue a)
    {
#ifdef USE_LOGVALUES
        return exp(a);
#else
        return a;
#endif
    }

    //////////////////////////////////////////////////////////////////////////////
    // Returns the sum of two probability values in log space
    template <typename TValue>
    inline void
    _probAdd(TValue &a, TValue b)
    {
#ifdef USE_LOGVALUES
        if (isInf(a)) {
            a = b;
            return;
        }
        if (isInf(b)) return;
        if (isNan(a + log(1 + exp(b - a)))) return;
        a += log(1 + exp(b - a));
#else
        a += b;
#endif
    }

    template <typename TValue>
    inline TValue
    _probMul(TValue a, TValue b)
    {
#ifdef USE_LOGVALUES
        return a + b;
#else
        return a * b;
#endif
    }

    template <typename TValue>
    inline TValue
    _probDiv(TValue a, TValue b)
    {
#ifdef USE_LOGVALUES
        return a - b;
#else
        return a / b;
#endif
    }


struct ErrorPatternLess
{
    template <typename TPattern>
    bool operator() (TPattern const &a, TPattern const &b) const
    {
        typedef typename Iterator<TPattern const>::Type TIter;
        TIter itA = end(a, Standard());
        TIter itB = end(b, Standard());
        TIter itEnd;
        if (length(a) <= length(b))
        {
            itEnd = begin(a, Standard());
            for (; itA != itEnd;)
            {
                --itA;
                --itB;
                if (*itA < *itB) return true;
                if (*itA > *itB) return false;
            }
            return false;
        } else
        {
            itEnd = begin(b, Standard());
            for (; itB != itEnd;)
            {
                --itA;
                --itB;
                if (*itA < *itB) return true;
                if (*itA > *itB) return false;
            }
            return true;
        }
    }
};

template <typename TPatternStore, typename TPattern>
inline int
_getErrorPatternIndex(TPatternStore const &patternStore, TPattern const &pattern)
{
    typedef typename Iterator<TPatternStore const>::Type TIter;
    TIter lb = std::lower_bound(begin(patternStore, Standard()), end(patternStore, Standard()), pattern, ErrorPatternLess());
    TIter invalid = end(patternStore, Standard());
    if (lb != invalid && *lb == pattern) {
//        std::cout << pattern;
        return lb - begin(patternStore, Standard());
    } else {
/*        std::cerr << "  !Pattern Not Found! " << pattern;
        if (lb != invalid) std::cerr << "\tnext is " << *lb;
        std::cerr << std::endl;
*/        return -1;
    }
}

// Cut 1 read character and trailing INSERTs of the pattern
template <typename TPattern>
inline int
_cutErrorPattern(TPattern &_pattern)
{
    typedef typename Iterator<TPattern const, Standard>::Type TIter;
    TPattern const & pattern = const_cast<TPattern const&>(_pattern);
    TIter it = end(pattern, Standard());
    int cuttedErrors = -2;

    // cut trailing INSERTs
    do {
        --it;
        ++cuttedErrors;
    } while ((int)getValue(it) == SEQAN_INSERT);

    // cut non INSERT
    if ((int)getValue(it) != SEQAN_MATCH)
        ++cuttedErrors;

    //  and all adjacent INSERTs
    do {
        --it;
        ++cuttedErrors;
    } while ((int)getValue(it) == SEQAN_INSERT);

    resize(_pattern, 1 + (it - begin(pattern, Standard())));
    return cuttedErrors;
}

template < typename TLogErrorDistr >
typename Value<TLogErrorDistr>::Type
_getProb(TLogErrorDistr const &logError, int errorType, int readPos)
{
    int maxN = length(logError) / 4;
    SEQAN_ASSERT(readPos >= 0 && readPos < maxN);
    return logError[maxN * (int)errorType + readPos];
}

//////////////////////////////////////////////////////////////////////////////
// Returns log probability of q-gram-configuration q ending at position pos in sequence
template < typename TState, typename TLogErrorDistr, typename TPattern >
inline void
_getLastPatternProb(TState &state, TLogErrorDistr const &logError, TPattern const &pattern, int span)
{
    int maxN = length(logError) / 4;
    typename Value<TLogErrorDistr>::Type prob = _transform(1.0);
    for (int i = 0, j = 0; j < (int)length(pattern); ++j)
    {
        prob = _probMul(prob, _getProb(logError, getValue(pattern, j), maxN - span + i));
        if ((int)getValue(pattern, j) != SEQAN_INSERT)
            ++i;
    }
    state.prob = prob;
}

template < typename TState, typename TPattern >
inline void
_getLastPatternProb(TState &, Nothing const &, TPattern const &, int)
{
}


//////////////////////////////////////////////////////////////////////////////
// Initialize states-string for edit/hamming-distance filters
template <
    typename TStateString,
    typename TShape,
    typename TLogErrorDistr,
    typename TDistance >
void initPatterns(
    TStateString &states,                // resulting states-string
    TShape const &bitShape,                // bit-string of the shape
    int maxErrors,                        // allowed errors per pattern
    TLogErrorDistr const &logError,        // error distribution (Nothing or string of 4*patternLen floats)
    TDistance,                            // enumerate hamming or edit distance patterns
    bool optionMinOutput)                // omit output
{
#ifndef DEBUG_RECOG_DP
//    typedef String<ErrorAlphabet, Packed<ErrorPackedString> >    TPattern;
    typedef String<ErrorAlphabet>                                TPattern;
#endif

    typedef typename Iterator<TPattern, Standard>::Type            TIter;
    typedef typename Value<TStateString>::Type                    TState;

    ErrorType lastErrorType = (IsSameType<TDistance, HammingDistance>::VALUE)? SEQAN_MISMATCH: SEQAN_DELETE;

    SEQAN_ASSERT_EQ(SEQAN_MATCH, 0);
    SEQAN_ASSERT((IsSameType<TLogErrorDistr,Nothing>::VALUE || (length(logError) % 4) == 0u));

#ifndef DEBUG_RECOG_DP
    String<TPattern> patternStore;
#endif

    // a modifier is a pair of position and error type
    String<Pair<int, ErrorType> > mods;
    resize(mods, maxErrors, Pair<int, ErrorType> (0, SEQAN_MATCH));

    TPattern pattern;
    int span = length(bitShape);

    //////////////////////////////////////////////////////////////////////////////
    // Enumerate all edit-modification patterns with up to k errors
    if (maxErrors == 0)
    {
        resize(pattern, span, (ErrorAlphabet)SEQAN_MATCH);
        appendValue(patternStore, pattern, Generous());
    }
    else
    do
    {
        clear(pattern);
        resize(pattern, span, (ErrorAlphabet)SEQAN_MATCH);

        // place errors in the pattern
        bool skip = false;
        for (int i = 0; (i < maxErrors) && !skip; ++i)
        {
//            std::cout << mods[i].i1 << " " << (ErrorAlphabet)mods[i].i2 << "\t";
            switch (mods[i].i2)
            {
            case SEQAN_MISMATCH:
            case SEQAN_DELETE:
                if (pattern[mods[i].i1] != (ErrorAlphabet)SEQAN_MATCH)
                {
                    skip = true;
                    break;
                }
                pattern[mods[i].i1] = (ErrorAlphabet)mods[i].i2;
                break;

            case SEQAN_INSERT:
                insertValue(pattern, mods[i].i1, (ErrorAlphabet)SEQAN_INSERT);
                break;

            case SEQAN_MATCH:
                break;
            }
        }

        // remove redundant patterns
        if (!skip)
        {
            TIter it = begin(pattern, Standard());
            TIter itEnd = end(pattern, Standard());
            int left = getValue(it);
            int right;
            for (++it; (it != itEnd) && !skip; ++it, left = right)
            {
                right = getValue(it);

#ifdef NON_REDUNDANT
                if (left == SEQAN_MISMATCH && right == SEQAN_DELETE)
                    skip = true;    // MISMATCH before DELETE is DELETE before MISMATCH (already enumerated)

                if (left == SEQAN_MISMATCH && right == SEQAN_INSERT)
                    skip = true;    // MISMATCH before INSERT is INSERT before MISMATCH (already enumerated)

                if (left == SEQAN_INSERT && right == SEQAN_DELETE)
                    skip = true;    // INSERT before DELETE is one MISMATCH (already enumerated)

                if (left == SEQAN_DELETE && right == SEQAN_INSERT)
                    skip = true;    // DELETE before INSERT is one MISMATCH (already enumerated)
#endif
            }
            if (left == SEQAN_INSERT)
                skip = true;        // no trailing INSERT allowed
        }

        if (!skip)
        {
            appendValue(patternStore, pattern, Generous());
//            std::cout << pattern << std::endl;
        }

        // reposition modifiers
        int i = 0;
        for (; i < maxErrors; ++i)
        {
            if (mods[i].i2 == SEQAN_MATCH) continue;
            int endPos = (mods[i].i2 == SEQAN_INSERT)? span + 1: span;
            if (++mods[i].i1 < endPos)
            {
                for(--i; i >= 0; --i)
                    mods[i].i1 = mods[i + 1].i1;
                break;
            }
        }

        if (i < maxErrors) continue;

        for (i = 0; i < maxErrors; ++i)
            mods[i].i1 = 0;

        // next state combination
        for (i = 0; i < maxErrors; ++i)
        {
            if (mods[i].i2 == lastErrorType) continue;
            mods[i].i2 = (ErrorType)(mods[i].i2 + 1);
            for(--i; i >= 0; --i)
                mods[i].i2 = SEQAN_MISMATCH;
            break;
        }

        if (i == maxErrors) break;

    } while (true);

    if (!optionMinOutput)
        std::cout << "Stored " << length(patternStore) << " modification patterns" << std::flush;

    reserve(patternStore, length(patternStore), Exact());
    std::sort(begin(patternStore, Standard()), end(patternStore, Standard()), ErrorPatternLess());
    for (int p = 1; p < (int)length(patternStore); ++p)
    {
        if (patternStore[p-1] == patternStore[p])
            std::cerr << "  !Found duplicate! " << patternStore[p] << std::endl;
    }

    if (!optionMinOutput)
        std::cout << " and sorted them." << std::endl;

    //////////////////////////////////////////////////////////////////////////////
    // Calculate transitions
    resize(states, length(patternStore));
    for (int p = 0; p < (int)length(patternStore); ++p)
    {
        pattern = patternStore[p];
        TState &state = states[p];

//        std::cout << pattern << "\t";

        // count errors of current pattern
        int errors = 0;
        for (int i = 0; i < (int)length(pattern); ++i)
            if ((int)getValue(pattern, i) != SEQAN_MATCH)
                ++errors;

        state.len = length(pattern);
        state.errors = errors;
        state.intermediate = (int)getValue(pattern, 0) == SEQAN_INSERT;
        _getLastPatternProb(state, logError, pattern, span);
//        std::cout << pattern << "\t";

        state.skipFirst = false;
        state.skipLast = false;

#ifdef NON_REDUNDANT
        int err = 0, del = 0;
        for (int j = 0; j < (int)length(pattern); ++j)
        {
            switch ((int)getValue(pattern, j)) {
                case SEQAN_MATCH:
                    ++del;
                    break;

                case SEQAN_DELETE:
                    ++del;

                case SEQAN_INSERT:
                    ++err;
                    break;

                default:;
            }
            if (del > 0 && del <= err)
                state.skipFirst = true;
        }
        err = del = 0;
        for (int j = (int)length(pattern) - 1; j >= 0; --j)
        {
            switch ((int)getValue(pattern, j)) {
                case SEQAN_MATCH:
                    ++del;
                    break;

                case SEQAN_DELETE:
                    ++del;

                case SEQAN_INSERT:
                    ++err;
                    break;

                default:;
            }
            if (del > 0 && del <= err)
                state.skipLast = true;
        }
#else
        state.skipFirst = (int)getValue(pattern, 0) == SEQAN_INSERT;
#endif
        // apply pattern to read q-gram
        // and check if shape is recognized in the genome
        state.qgramHit = false;
        int delta = 0;
        for (int j = 0, readPos = 0, genomePos = 0; j < (int)length(pattern); ++j)
        {
            switch ((int)getValue(pattern, j))
            {
                case SEQAN_MATCH:
                    if (readPos == 0) {
                        // SEQAN_ASSERT_EQ(bitShape[0], '1')
                        delta = genomePos;
                        state.qgramHit = true;
                    } else
                        if (bitShape[readPos] == '1')
                            state.qgramHit &= (readPos + delta == genomePos);
//                    std::cout << readPos;
                    ++readPos; ++genomePos;
                    break;
                case SEQAN_MISMATCH:
                    // was it a relevant read position?
                    if (bitShape[readPos] == '1')
                        state.qgramHit = false;
//                    std::cout << 'x';
                    ++readPos; ++genomePos;
                    break;
                case SEQAN_DELETE:
                    // was it a relevant read position?
                    if (bitShape[readPos] == '1')
                        state.qgramHit = false;
                    ++readPos;
                    break;
                case SEQAN_INSERT:
                    ++genomePos;
//                    std::cout << 'x';
            }
        }
//        std::cout << std::endl;

        // prepend INSERT
        ++errors;
        insertValue(pattern, 0, SEQAN_INSERT);
        if ((int)SEQAN_INSERT < (int)state.TRANSITIONS)
        {
            if (errors <= maxErrors)
                state.transition[SEQAN_INSERT] = _getErrorPatternIndex(patternStore, pattern);
            else
                state.transition[SEQAN_INSERT] = -1;
        }

        // prepend MISMATCH and cut INSERTS
        errors -= _cutErrorPattern(pattern);
        if ((int)SEQAN_MISMATCH < (int)state.TRANSITIONS)
        {
            pattern[0] = SEQAN_MISMATCH;
            if (errors <= maxErrors)
                state.transition[SEQAN_MISMATCH] = _getErrorPatternIndex(patternStore, pattern);
            else
                state.transition[SEQAN_MISMATCH] = -1;
        }

        // prepend DELETE
        if ((int)SEQAN_DELETE < (int)state.TRANSITIONS)
        {
            pattern[0] = SEQAN_DELETE;
            if (errors <= maxErrors)
                state.transition[SEQAN_DELETE] = _getErrorPatternIndex(patternStore, pattern);
            else
                state.transition[SEQAN_DELETE] = -1;
        }

        // prepend MATCH
        if ((int)SEQAN_MATCH < (int)state.TRANSITIONS)
        {
            --errors;
            pattern[0] = SEQAN_MATCH;
            if (errors <= maxErrors)
                state.transition[SEQAN_MATCH] = _getErrorPatternIndex(patternStore, pattern);
            else
                state.transition[SEQAN_MATCH] = -1;
        }
/*
        std::cout << "\t" << state.errors;
        std::cout << "\t" << state.qgramHit;
        std::cout << "\t" << state.leftError;
        std::cout << "\t" << state.rightError;
        std::cout << "\t" << state.transition[0];
        std::cout << "\t" << state.transition[1];
        std::cout << "\t" << state.transition[2];
        std::cout << "\t" << state.transition[3];
        std::cout << std::endl;
*/    }
    if (!optionMinOutput)
        std::cout << "Preprocessing finished." << std::endl;
}

//////////////////////////////////////////////////////////////////////////////
// Compute filtering loss of any q-gram filter (given a states-string)
template <
    typename TThreshString,
    typename TStateString >
void computeExactQGramThreshold(
    TThreshString &treshPerError,
    TStateString const &states,
    int span,
    int maxErrors,
    int maxN,
    bool optionMinOutput)
{
    typedef typename Value<TStateString>::Type        TState;
    typedef unsigned                                TThresh;
    typedef String<TThresh>                            TMatrixCol;

    int statesCount = length(states);
//    int span = length(bitShape);

    // columns n-1 and n for recursion
    TMatrixCol col0;    // addressing is colx[errors * statesCount + state]
    TMatrixCol col1;
    const TThresh infty = MaxValue<TThresh>::VALUE >> 1;

    resize(col0, maxErrors * statesCount, infty);
    resize(col1, maxErrors * statesCount);

    // RECURSION BEGIN
    for (int s = 0; s < statesCount; ++s)
    {
        TState const &state = states[s];
        if (state.skipFirst) continue;

        // threshold is 1 iff we have a q-gram hit at the end
        col0[s] = (state.qgramHit)? 1: 0;
    }

    // iterate over sequence length n
    TMatrixCol *col = &col1;
    TMatrixCol *colPrev = &col0;

#ifdef DEBUG_RECOG_DP
    std::cout << span << ":0";
    dump(col0, 0,statesCount);
    std::cout << " :1";
    dump(col0, 1,statesCount);
#endif


    // RECURSION
    //
    // thresh(n,q,e) = min(thresh(n-1,0|(q>>1),e),              delta=1/0 <-> q hat 0/>0 error
    //                     thresh(n-1,1|(q>>1),e-1)) + delta

    for (int n = span; n < maxN; ++n)
    {
        for (int e = 0; e < maxErrors * statesCount; e += statesCount)
        {
            for (int s = 0; s < statesCount; ++s)
            {
                TState const &state = states[s];

                // MATCH
                TThresh t = (*colPrev)[e + state.transition[SEQAN_MATCH]];

                // MISMATCH, INSERT, DELETE
                if (e > 0)
                    for (int m = SEQAN_MISMATCH; m < TState::TRANSITIONS; ++m)
                    {
                        int prevState = state.transition[m];
                        if (prevState >= 0)
                        {
                            if (m == SEQAN_INSERT)
                                t = _min(t, (*col)[(e - statesCount) + prevState]);
                            else
                                t = _min(t, (*colPrev)[(e - statesCount) + prevState]);
                        }
                    }

                (*col)[e + s] = t + state.qgramHit;
            }
            if (!optionMinOutput)
                std::cout << '.' << std::flush;
        }

        TMatrixCol *tmp = col;
        col = colPrev;
        colPrev = tmp;

#ifdef DEBUG_RECOG_DP
        std::cout << n+1 << ":0";
        dump(*colPrev, 0,statesCount);
        std::cout << " :1";
        dump(*colPrev, 1,statesCount);
        std::cout << " :2";
        dump(*colPrev, 2,statesCount);
#endif
    }

    if (!optionMinOutput)
        std::cout << std::endl;

    resize(treshPerError, maxErrors);

    // RECURSION END
    for (int eSum = 0; eSum < maxErrors; ++eSum)
    {
        TThresh t = infty;
        for (int s = 0; s < statesCount; ++s)
        {
            TState const &state = states[s];

            // skip intermediate results
            if (state.intermediate || state.skipLast) continue;
            if (state.errors <= eSum)
            {
                int e = eSum - state.errors;
                // multiply probability for the trailing pattern
                t = _min(t, (*colPrev)[e * statesCount + s]);
            }
        }

        if (t >= infty) t = 0;
        treshPerError[eSum] = t;
    }
}


//////////////////////////////////////////////////////////////////////////////
// Compute filtering loss of any q-gram filter (given a states-string)
template <
    typename TLossString,
    typename TLogErrorDistr,
    typename TStateString >
void computeQGramFilteringSensitivity(
    TLossString &found,
    TStateString const &states,
    int span,
    int maxT,
    int maxErrors,
    TLogErrorDistr const &logError,
//    bool optionAbsolute = false,
    bool optionMinOutput)
{
    typedef typename Value<TLossString>::Type        TFloat;
    //typedef typename Value<TLogErrorDistr>::Type    TProbValue;
    typedef typename Value<TStateString>::Type        TState;

    typedef String<TFloat>                            TMatrixCol;
    //typedef String<int>                                TIntCol;

    SEQAN_ASSERT_EQ((length(logError) % 4), 0u);

    int maxN = length(logError) / 4;
    int statesCount = length(states);
    const bool optionAbsolute = false;
//    int span = length(bitShape);

    // columns n-1 and n for recursion
    TMatrixCol col0;
    TMatrixCol col1;
    resize(col0, maxErrors * statesCount * maxT, (TFloat)_transform(0.0));
    resize(col1, maxErrors * statesCount * maxT);

#ifdef COUNT_LOSSES
    TFloat positive = _transform(0.0);
    TFloat negative = _transform(1.0);
#else
    TFloat positive = _transform(1.0);
    TFloat negative = _transform(0.0);
#endif

    // RECURSION BEGIN
    for (int s = 0; s < statesCount; ++s)
    {
        TState const &state = states[s];

        if (state.skipFirst) continue;

        // we miss no match if threshold t is 0
        col0[s*maxT] = positive;

        // for n==0
        if (state.qgramHit)
        {
            // we miss no match if read q-gram is recognized
            // --> probability of finding this MMP is 1, if t=1
            col0[s*maxT+1] = positive;
            // --> probability of finding this MMP is 0, if t>1
            for (int t = 2; t < maxT; ++t)
                col0[s*maxT+t] = negative;
        } else
        {
            // we miss 1 match if t>0 and read q-gram is not recognized
            // --> probability of finding this MMP is 0, if t>=1
            for (int t = 1; t < maxT; ++t)
                col0[s*maxT+t] = negative;
        }
    }

    // iterate over sequence length n
    TMatrixCol *col = &col1;
    TMatrixCol *colPrev = &col0;

#ifdef DEBUG_RECOG_DP
    std::cout << span << ":0";
    dump(col0, 0,statesCount);
    std::cout << " :1";
    dump(col0, 1,statesCount);
#endif


    // RECURSION
    //
    // found(n,q,t,e) = (1-errorProb[n-span]) * found(n-1,0|(q>>1),t-delta,e) delta=1/0 <-> q hat 0/>0 fehler
    //               + errorProb[n-span] * found(n-1,1|(q>>1),t-delta,e-1)

    // rekursion (fuer q-gram matches <=1 fehler)
    // found(n,q,t,e) = (1-errorProb[n-span]) * found(n-1,0|(q>>1),t-delta,e) delta=1/0 <-> q hat <=1/>1 fehler
    //               + errorProb[n-span] * found(n-1,1|(q>>1),t-delta,e-1)

    for (int n = span; n < maxN; ++n)
    {
        for (int e = 0; e < maxErrors * statesCount; e += statesCount)
        {
            for (int s = 0; s < statesCount; ++s)
            {
                TState const &state = states[s];
                for (int t = 0; t < maxT; ++t)
                {
                    int _t = t;
                    if (_t > 0 && state.qgramHit) --_t;

                    // MATCH
                    TFloat recovered = _probMul(
                        _getProb(logError, SEQAN_MATCH, n-span),
                        (*colPrev)[(e+state.transition[SEQAN_MATCH])*maxT+_t]);

                    // MISMATCH, INSERT, DELETE
                    for (int m = SEQAN_MISMATCH; m < (int)state.TRANSITIONS; ++m)
                        if (e > 0)
                        {
                            int prevState = state.transition[m];
                            if (prevState >= 0)
                            {
                                if (m == SEQAN_INSERT)
                                    _probAdd(recovered, _probMul(_getProb(logError,m,n-span), (*col)[((e-statesCount)+prevState)*maxT+t]));
                                else
                                    _probAdd(recovered, _probMul(_getProb(logError,m,n-span), (*colPrev)[((e-statesCount)+prevState)*maxT+_t]));
                            }
                        }
                    (*col)[(e+s)*maxT+t] = recovered;
                }
            }
            if (!optionMinOutput)
                std::cout << '.' << std::flush;
        }

        TMatrixCol *tmp = col;
        col = colPrev;
        colPrev = tmp;

#ifdef DEBUG_RECOG_DP
        std::cout << n+1 << ":0";
        dump(*colPrev, 0,statesCount);
        std::cout << " :1";
        dump(*colPrev, 1,statesCount);
        std::cout << " :2";
        dump(*colPrev, 2,statesCount);
#endif
    }

    if (!optionMinOutput)
        std::cout << std::endl;

    // RECURSION END
    for (int eSum = 0; eSum < maxErrors; ++eSum)
        for (int t = 0; t < maxT; ++t)
        {
            TFloat recovered = _transform(0.0);
            for (int s = 0; s < statesCount; ++s)
            {
                TState const &state = states[s];

                // skip intermediate results
                if (state.intermediate || state.skipLast) continue;
                if (state.errors <= eSum)
                {
                    int e = eSum - state.errors;
                    // multiply probability for the trailing pattern
                    _probAdd(recovered, _probMul(state.prob, (*colPrev)[(e*statesCount+s)*maxT+t]));
                }
            }

#ifndef COUNT_LOSSES
            // we can only normalize probs if t==0 contains all k-pattern probs
            if (t > 0 && !optionAbsolute)
                recovered = _probDiv(recovered, found[eSum*maxT]);
#endif

            found[eSum*maxT+t] = recovered;
        }
}


//////////////////////////////////////////////////////////////////////////////
// q-gram threshold DP algorithm
//
// - exact threshold
//////////////////////////////////////////////////////////////////////////////

template <typename TShape, typename TPatternSize, typename TErrors, typename TDistance>
int qgramThreshold(TShape const & shape, TPatternSize patternLength, TErrors errors, TDistance const dist, ThreshExact const)
{
    String<ThreshDPState_<TDistance> > states;
    String<unsigned> thresh;
    String<char> bitString;

    shapeToString(bitString, shape);
    initPatterns(states, bitString, errors, Nothing(), dist, true);
    computeExactQGramThreshold(thresh, states, length(bitString), errors + 1, patternLength, true);

    return thresh[errors];
}

//////////////////////////////////////////////////////////////////////////////
// q-gram filter sensitivity DP algorithm
//
// - exact threshold
//////////////////////////////////////////////////////////////////////////////

template <typename TSensitivityMatrix, typename TShape, typename TPatternSize, typename TErrors, typename TThresh, typename TDistance, typename TErrorDist>
void qgramFilteringSensitivity(
    TSensitivityMatrix & sensMat,
    TShape const & shape,
    TPatternSize patternLength,
    TErrors errors,
    TThresh maxThresh,
    TDistance const dist,
    ThreshExact const,
    TErrorDist const & logErrorDistribution)
{
    typedef typename Value<TSensitivityMatrix>::Type TFloat;
    String<SensitivityDPState_<TDistance, TFloat> > states;
    String<unsigned> thresh;
    String<char> bitString;

    maxThresh = _min(maxThresh, patternLength - length(shape) + 1);
    resize(sensMat, (maxThresh + 1) * (errors + 1));
    shapeToString(bitString, shape);

    initPatterns(states, bitString, errors, logErrorDistribution, dist, true);
    computeQGramFilteringSensitivity(sensMat, states, length(bitString), maxThresh + 1, errors + 1, logErrorDistribution, true);
}

}    // namespace seqan

#endif
