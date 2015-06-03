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

#ifndef SEQAN_HEADER_FIND_SWIFT_H
#define SEQAN_HEADER_FIND_SWIFT_H

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////
// SWIFT to search a text for
// - semi-global alignments of one/multiple short sequences
// - local epsilon matches of one/multiple short sequences
//////////////////////////////////////////////////////////////////////////////

// TODO(bkehr): Is this documentatin right? Should the Specializations be Tags?
// weese: of minimal length?
/*!
 * @class SwiftPattern
 * @extends Pattern
 * @headerfile <seqan/index.h>
 * @brief Pattern for SWIFT search, must be used together with @link SwiftFinder @endlink.
 *
 * @signature template <typename TIndex, typename TSpec>
 *            class Pattern<TIndex, Swift<TSpec> >;
 *
 * @tparam TSpec  Specifies the type of Swift filter.
 * @tparam TIndex A q-gram index of needle(s).  Types: @link IndexQGram @endlink
 *
 * The @link Pattern @endlink must be a @link IndexQGram @endlink over multiple patterns.  The allowed error rate must
 * be given when @link Finder#find @endlink or @link SwiftFinder#windowFindBegin @endlink is called.
 *
 * Provides a fast filter alogrithm that guarantees to find all regions overlapping with potential &epsilon;-matches.
 * An &epsilon;-match is a matching region of minimal length and an error rate of at most &epsilon;.
 *
 * @see SwiftFinder
 *
 * @var SwiftParameters SwiftPattern::parameters;
 * @brief The SWIFT parameters to use for configuration.
 */

/*!
 * @class SwiftFinder
 * @extends Finder
 * @headerfile <seqan/index.h>
 * @brief Finder for SWIFT search, must be used together with @link SwiftPattern @endlink.
 *
 * @signature template <typename THaystack, typename TSpec>
 *            class Finder<THaystack, Swift<TSpec> >;
 *
 * @tparam TSpec     Specifies the type of Swift filter.
 * @tparam THaystack A haystack type.  Types: @link Index @endlink, @link String @endlink, @link StringSet @endlink.
 *
 * See @link SwiftPattern @endlink for an overview of the SWIFT algorithm.
 *
 * @see SwiftPattern
 */

/*!
 * @class SwiftLocalPattern
 * @extends SwiftPattern
 * @headerfile <seqan/index.h>
 *
 * @brief The specialization for the general SWIFT filter that finds &epsilon;-matches between haystack and needle.
 *
 * @signature template <typename THaystack>
 *            class Finder<THaystack, Swift<SwiftLocal> >;
 *
 * @tparam THaystack A haystack type. Types: @link Index @endlink, @link String @endlink, @link StringSet @endlink.
 */

/*!
 * @class SwiftLocalFinder
 * @extends SwiftFinder
 * @headerfile <seqan/index.h>
 * @brief The specialization for the general SWIFT filter that finds &epsilon;-matches between haystack and needle.
 *
 * @signature template <typename THaystack>
 *            class Finder<THaystack, Swift<SwiftLocal> >;
 *
 * @tparam THaystack A haystack type. Types: @link Index @endlink, @link String @endlink, @link StringSet @endlink.
 */

/*!
 * @class SwiftSemiGlobalPattern
 * @extends SwiftPattern
 * @headerfile <seqan/index.h>
 * @brief The specialization for the semi-global swift filter that finds regions of the haystack where a needle matches
 *        with an error rate less than \epsilon.
 *
 * @signature template <typename TIndex>
 *            class Pattern<TIndex, Swift<SwiftSemiGlobal> >;
 *
 * @tparam TIndex A @link IndexQGram q-gram index @endlink of needle(s).
 */

/*!
 * @class SwiftSemiGlobalFinder
 * @extends SwiftFinder
 * @headerfile <seqan/index.h>
 * @brief The specialization for the semi-global swift filter that finds regions of the haystack where a needle matches
 *        with an error rate less than \epsilon.
 *
 * @signature template <typename THaystack>
 *            class Finder<THaystack, Swift<SwiftSemiGlobal> >;
 *
 * @tparam THaystack A haystack type. Types: @link Index @endlink, @link String @endlink, @link StringSet @endlink
 */

template < typename TObject, typename TSpec > class Index;
template < typename TObject > struct SAValue;

struct SwiftLocal_;
typedef Tag<SwiftLocal_> SwiftLocal;

template <typename TSpec = void>
struct SwiftSemiGlobal_;
typedef Tag<SwiftSemiGlobal_<void> > SwiftSemiGlobal;

struct Hamming_;
typedef Tag<SwiftSemiGlobal_<Hamming_> > SwiftSemiGlobalHamming;


template <typename TSpec = SwiftSemiGlobal>
struct Swift;

template <>
struct Swift<SwiftSemiGlobal> {
    enum { SEMIGLOBAL = 1 };        // 0..match eps-match of min.length n0; 1..match the whole read
    enum { DIAGONAL = 1 };          // 0..use rectangular buckets (QUASAR); 1..use diagonal buckets (SWIFT)
    enum { QGRAM_ERRORS = 0 };      // q-gram must match exactly
    enum { HAMMING_ONLY = 0 };      // 0..no indels; 1..allow indels
    enum { PARAMS_BY_LENGTH = 1 };  // params are determined only by seq.length
};

template <>
struct Swift<SwiftSemiGlobalHamming> {
    enum { SEMIGLOBAL = 1 };        // 0..match eps-match of min.length n0; 1..match the whole read
    enum { DIAGONAL = 1 };          // 0..use rectangular buckets (QUASAR); 1..use diagonal buckets (SWIFT)
    enum { QGRAM_ERRORS = 0 };      // q-gram must match exactly
    enum { HAMMING_ONLY = 1 };      // 0..no indels; 1..allow indels
    enum { PARAMS_BY_LENGTH = 1 };  // params are determined only by seq.length
};

template <>
struct Swift<SwiftLocal> {
    enum { SEMIGLOBAL = 0 };        // 0..match eps-match of min.length n0; 1..match the whole read
    enum { DIAGONAL = 1 };          // 0..use rectangular buckets (QUASAR); 1..use diagonal buckets (SWIFT)
    enum { QGRAM_ERRORS = 0 };      // allow 0 errors per q-gram
    enum { HAMMING_ONLY = 0 };      // 0..no indels; 1..allow indels
    enum { PARAMS_BY_LENGTH = 0 };  // params are determined by seq.no.
};

/*!
 * @class SwiftParameters
 * @headerfile <seqan/index.h>
 * @brief Parameters for the SWIFT algorithm.
 *
 * @signature struct SwiftParameters;
 *
 * @see SwiftPattern
 *
 *
 * @fn SwiftParameters::SwiftParameters
 * @brief Default constructor.
 *
 * @signature SwiftParameters::SwiftParameters();
 *
 * Initializes all member variables, default values are given in the variable documentation.
 *
 *
 * @var int SwiftParameters::minThreshold;
 * @brief Minimal threshold to use initially, defaults to <tt>1</tt>.
 *
 * @var int SwiftParameters::minLog2Delta;
 * @brief Minimal delta to use, defaults to <tt>16</tt>.
 *
 * @var int SwiftParameters::tabooLength;
 * @brief Taboo length to use, defaults to <tt>1</tt>.
 *
 * @var bool SwiftParameters::printDots;
 * @brief Whether to print dots for the scanning progress, defaults to <tt>false</tt>.
 *
 * @var bool SwiftParameters::debug;
 * @brief Whether to print debug information, defaults to <tt>false</tt>.
 */

struct SwiftParameters
{
    int minThreshold;
    int minLog2Delta;
    int tabooLength;
    bool printDots;
    bool debug;

    SwiftParameters() :
        minThreshold(1),        // set minimal threshold to 1
        minLog2Delta(4),        // set minimal delta to 16
        tabooLength(1),         // minimal genomic distance between q-gram hits
        printDots(false),       // print a . for every 100kbps mapped genome
        debug(false)
    {}
};

//////////////////////////////////////////////////////////////////////////////

    template <typename TSpec, typename TSize_, typename TShortSize_ = unsigned short>
    struct SwiftBucket_
    {
        typedef TSize_          TSize;
        typedef TShortSize_ TShortSize;

        TSize                   firstIncrement;
        TSize                   lastIncrement;
        TShortSize              counter;        // q-gram hits
        TShortSize              threshold;      // at least threshold q-gram hits induce an approx match
        bool                    notListed;         // true if bucket is not listed in the patterns' verify list
#ifdef SEQAN_DEBUG_SWIFT
        TSize                   _lastIncDiag;
#endif
    };

    template <typename TSpec_, typename TSize_, typename TShortSize_>
    struct SwiftBucket_<SwiftSemiGlobal_<TSpec_>, TSize_, TShortSize_>
    {
        typedef TSize_          TSize;
        typedef TShortSize_     TShortSize;

        TSize                   lastIncrement;
        TShortSize              counter;        // q-gram hits
        TShortSize              threshold;      // at least threshold q-gram hits induce an approx match
#ifdef SEQAN_DEBUG_SWIFT
        int                     _lastIncDiag;
#endif
    };

    template <typename TSpec, typename TSize_, typename TShortSize_ = unsigned short>
    struct SwiftBucketParams_
    {
        typedef TSize_          TSize;
        typedef TShortSize_     TShortSize;

        TSize           firstBucket;    // first SwiftBucket_ entry in pattern.buckets
        TSize           reuseMask;      // 2^ceil(log2(x)) reuse every x-th bucket)
        TShortSize      threshold;      // at least threshold q-gram hits induce an approx match
        TShortSize      distanceCut;    // if lastIncrement is this far or farer away, threshold can't be reached
        TShortSize      delta;          // buckets begin at multiples of delta
        TShortSize      overlap;        // number of diagonals/columns a bucket shares with its neighbor
        TShortSize      tabooLength;    // minimal genomic distance between q-gram hits
        unsigned char   logDelta;       // log2(delta)
    };

    template <typename TSpec_, typename TSize_, typename TShortSize_>
    struct SwiftBucketParams_< Swift<Tag<SwiftSemiGlobal_<TSpec_> > >, TSize_, TShortSize_>
    {
        typedef TSize_          TSize;
        typedef TShortSize_     TShortSize;

        TSize           firstBucket;    // first SwiftBucket_ entry in pattern.buckets
        TSize           reuseMask;      // 2^ceil(log2(x)) reuse every x-th bucket)
        TShortSize      threshold;      // at least threshold q-gram hits induce an approx match
        TShortSize      delta;          // buckets begin at multiples of delta
        TShortSize      overlap;        // number of diagonals/columns a bucket shares with its neighbor
        TShortSize      tabooLength;    // minimal genomic distance between q-gram hits
        unsigned char   logDelta;       // log2(delta)
    };

//____________________________________________________________________________


    template <typename THstkPos>
    struct SwiftHit_
    {
        THstkPos    hstkPos;            // parallelogram begin in haystack
        unsigned    ndlSeqNo;           // needle sequence number
        THstkPos    ndlPos;             // begin position of hit in needle
        unsigned    bucketWidth;        // (non-diagonal) bucket width (hitLengthNeedle + delta + overlap (for diagonals))
        unsigned    hitLengthNeedle;    // length of the hit in needle
    };

    template <typename THstkPos>
    struct SwiftHitSemiGlobal_
    {
        THstkPos    hstkPos;            // parallelogram begin in haystack
        unsigned    ndlSeqNo;           // needle sequence number
        unsigned    bucketWidth;        // (non-diagonal) bucket width (bktHeight + delta + overlap (for diagonals))
    };


    template <typename TFinder, typename TPattern = void>
    struct FindResult
    {
    };
    template <typename TFinder, typename TPattern = void>
    struct WindowFindResult
    {
        typedef String<typename FindResult<TFinder, TPattern>::Type> Type;
    };



    template <typename THaystack, typename TPattern, typename TSwiftSpec>
    struct FindResult<Finder<THaystack, Swift<TSwiftSpec> >, TPattern>
    {
        typedef SwiftHit_<__int64> Type;
    };

    template <typename THaystack, typename TPattern, typename TSpec>
    struct FindResult<Finder<THaystack, Swift<Tag<SwiftSemiGlobal_<TSpec> > > >, TPattern>
    {
        typedef SwiftHitSemiGlobal_<__int64> Type;
    };


//____________________________________________________________________________


/*!
 * @typedef SwiftFinder::THitString
 * @brief The type to use for the hit string returned in @link SwiftFinder#getWindowFindHits @endlink.
 *
 * @signature typedef (...) SwiftFinder::THitString;
 */

    template <typename THaystack, typename TSpec>
    class Finder<THaystack, Swift<TSpec> >
    {
    public:
        typedef typename Iterator<THaystack, Rooted>::Type          TIterator;
        typedef typename Position<THaystack>::Type                  THstkPos;
        typedef typename FindResult<Finder>::Type                   TSwiftHit;
        typedef typename WindowFindResult<Finder>::Type             THitString;
        typedef typename Iterator<THitString, Standard>::Type       THitIterator;
        typedef typename SAValue<THaystack>::Type                   TSAValue;
        typedef Repeat<TSAValue, unsigned>                          TRepeat;
        // TODO(holtgrew): Make this a holder so we can share the actual string when copying.
        typedef String<TRepeat>                                     TRepeatString;
        typedef typename Iterator<TRepeatString, Standard>::Type    TRepeatIterator;

        TIterator       data_iterator;
        TIterator       haystackEnd;
        bool            _needReinit;    // if true, the Pattern needs to be reinitialized
        THitString      hits;
        THitIterator    curHit, endHit;
        THstkPos        startPos, curPos, endPos;
        THstkPos        windowStart;
        THstkPos        dotPos, dotPos2;
        TRepeatString   data_repeats;
        TRepeatIterator curRepeat, endRepeat;

        Finder() :
            data_iterator(), haystackEnd(), _needReinit(true), curHit(), endHit(),
            startPos(), curPos(), endPos(), windowStart(), dotPos(), dotPos2(),
            curRepeat(), endRepeat()
        {}

        Finder(THaystack & haystack) :
            data_iterator(begin(haystack, Rooted())), haystackEnd(), _needReinit(true), curHit(), endHit(),
            startPos(), curPos(), endPos(), windowStart(), dotPos(), dotPos2(), curRepeat(), endRepeat()
        {}

        template <typename TRepeatSize, typename TPeriodSize>
        Finder(THaystack &haystack, TRepeatSize minRepeatLen, TPeriodSize maxPeriod) :
            data_iterator(begin(haystack, Rooted())), haystackEnd(), _needReinit(true), curHit(), endHit(),
            startPos(), curPos(), endPos(), windowStart(), dotPos(), dotPos2(), curRepeat(), endRepeat()
        {
            findRepeats(data_repeats, haystack, minRepeatLen, maxPeriod);
        }

        Finder(TIterator &iter) :
            data_iterator(iter), haystackEnd(), _needReinit(true), curHit(), endHit(),
            startPos(), curPos(), endPos(), windowStart(), dotPos(), dotPos2(), curRepeat(), endRepeat()
        {}

        Finder(TIterator const & iter):
            data_iterator(iter), haystackEnd(), _needReinit(true), curHit(), endHit(),
            startPos(), curPos(), endPos(), windowStart(), dotPos(), dotPos2(), curRepeat(), endRepeat()
        {}

        Finder(Finder const & orig):
            data_iterator(orig.data_iterator),
            haystackEnd(orig.haystackEnd),
            _needReinit(orig._needReinit),
            hits(orig.hits),
            startPos(orig.startPos),
            curPos(orig.curPos),
            endPos(orig.endPos),
            windowStart(orig.windowStart),
            dotPos(orig.dotPos),
            dotPos2(orig.dotPos2),
            data_repeats(orig.data_repeats)
        {
            curHit = begin(hits, Standard()) + (orig.curHit - begin(orig.hits, Standard()));
            endHit = end(hits, Standard());
            curRepeat = begin(data_repeats, Standard()) + (orig.curRepeat - begin(orig.data_repeats, Standard()));
            endRepeat = end(data_repeats, Standard());
        }

        inline typename Reference<TIterator>::Type
        operator* () { return value(hostIterator(*this)); }

        inline typename Reference<TIterator const>::Type
        operator* () const { return value(hostIterator(*this)); }

        operator TIterator () const { return data_iterator; }

        Finder & operator = (Finder const &orig)
        {
            data_iterator = orig.data_iterator;
            haystackEnd = orig.haystackEnd;
            _needReinit = orig._needReinit;
            hits = orig.hits;
            startPos = orig.startPos;
            windowStart = orig.windowStart;
            curPos = orig.curPos;
            endPos = orig.endPos;
            dotPos = orig.dotPos;
            dotPos2 = orig.dotPos2;
            data_repeats = orig.data_repeats;
            curHit = begin(hits, Standard()) + (orig.curHit - begin(orig.hits, Standard()));
            endHit = end(hits, Standard());
            curRepeat = begin(data_repeats, Standard()) + (orig.curRepeat - begin(orig.data_repeats, Standard()));
            endRepeat = end(data_repeats, Standard());
            return *this;
        }
    };


//____________________________________________________________________________

    // forward
    template < typename TInput, typename TSpec >
    struct Pipe;

    template < typename TTuples, typename TPipeSpec, typename TSpec >
    class Finder< Pipe<TTuples, TPipeSpec>, Swift<TSpec> >
    {
    public:
        typedef Pipe<TTuples, TPipeSpec>                        TInput;
        typedef typename Size<TInput>::Type                     THstkPos;
        typedef typename FindResult<Finder>::Type               TSwiftHit;
        typedef typename WindowFindResult<Finder>::Type         THitString;
        typedef typename Iterator<THitString, Standard>::Type   THitIterator;

        TInput          &in;
        bool            _needReinit;    // if true, the Pattern needs to be reinitialized
        THitString      hits;
        THitIterator    curHit, endHit;
        THstkPos        curPos, dotPos, dotPos2;

        Finder(TInput &_in):
            in(_in),
            _needReinit(true) {}

        Finder(Finder const &orig):
            in(orig.in),
            hits(orig.hits),
            _needReinit(orig._needReinit)
        {
            curHit = begin(hits, Standard()) + (orig.curHit - begin(orig.hits, Standard()));
            endHit = end(hits, Standard());
        }
    };


//____________________________________________________________________________


    template <typename THaystack, typename TSpec>
    inline bool
    atEnd(Finder<THaystack, Swift<TSpec> > & me)
    {
        return hostIterator(hostIterator(me)) == hostIterator(me.haystackEnd);
    }

    template <typename THaystack, typename TSpec>
    inline void
    goEnd(Finder<THaystack, Swift<TSpec> > & me)
    {
        hostIterator(me) = me.haystackEnd;
    }


//____________________________________________________________________________


    template <typename TIndex, typename TSpec>
    class Pattern<TIndex, Swift<TSpec> >
    {
    public:
        typedef typename Size<TIndex>::Type                             TSize;
        typedef unsigned                                                TShortSize;
        typedef typename Fibre<TIndex, Tag<FibreSA_> const >::Type      TSA;
        typedef typename Fibre<TIndex, Tag<Fibre_Shape_> const >::Type  TShape;
        typedef typename Iterator<TSA const, Standard>::Type            TIterator;

        typedef SwiftBucket_<TSpec, TSize, TShortSize>                  TBucket;
        typedef String<TBucket>                                         TBucketString;
        typedef SwiftBucketParams_<TSpec, TSize, TShortSize>            TBucketParams;
        typedef String<TBucketParams>                                   TBucketParamsString;
        typedef String<Pair<unsigned> >                                 TVerifyList;

        TShape                  shape;
        TBucketString           buckets;
        TBucketParamsString     bucketParams;
        TVerifyList             verifyList;                             // numbers of buckets that need to be verified
        SwiftParameters         params;
        unsigned                curSeqNo;
        __int64                 curBeginPos, curEndPos;
        TSize                   finderPosOffset;                        // these must be of
        TSize                   finderPosNextOffset;                    // type TSize of TBucket
        __int64                 finderLength;
        __int64                 maxPatternLength;

        double                  _currentErrorRate;
        int                     _currentMinLengthForAll;

        Holder<TIndex>  data_host;

        Pattern()
        {
            clear(*this);
        }

        Pattern(TIndex &_index): data_host(_index)
        {
            clear(*this);
        }

        Pattern(TIndex const &_index): data_host(_index)
        {
            clear(*this);
        }
    };

//____________________________________________________________________________


template <typename TSpec_, typename TSize, typename TShortSize>
inline void _printSwiftParams(SwiftBucketParams_<TSpec_, TSize, TShortSize > &bucketParams)
{
    std::cout << "  firstBucket: " << bucketParams.firstBucket << std::endl;
    std::cout << "  reuseMask:   " << bucketParams.reuseMask << std::endl;
    std::cout << "  distanceCut: " << bucketParams.distanceCut << std::endl;
    std::cout << "  delta:       " << bucketParams.delta << std::endl;
    std::cout << "  threshold:   " << bucketParams.threshold << std::endl;
    std::cout << "  overlap:     " << bucketParams.overlap << std::endl;
    std::cout << "  logDelta:    " << (int)bucketParams.logDelta << std::endl << std::endl;
}

template <typename TSpec_, typename TSize, typename TShortSize>
inline void _printSwiftParams(SwiftBucketParams_<Tag<SwiftSemiGlobal_<TSpec_> >, TSize, TShortSize > &bucketParams)
{
    std::cout << "  firstBucket: " << bucketParams.firstBucket << std::endl;
    std::cout << "  reuseMask:   " << bucketParams.reuseMask << std::endl;
    std::cout << "  delta:       " << bucketParams.delta << std::endl;
    std::cout << "  threshold:   " << bucketParams.threshold << std::endl;
    std::cout << "  overlap:     " << bucketParams.overlap << std::endl;
    std::cout << "  logDelta:    " << (int)bucketParams.logDelta << std::endl << std::endl;
}

template <typename TIndex, typename TSpec>
inline void _printSwiftBuckets(Pattern< TIndex, Swift<TSpec> > &p)
{
    typedef typename Pattern<TIndex, Swift<TSpec> >::TBucketParams TParams;

    unsigned j = 0;
    TParams *bucketParams = &_swiftBucketParams(p, 0);

    for(unsigned i=0; i<length(p.buckets) && i<10; ++i)
    {
        if ((i & bucketParams->reuseMask) == 0)
        {
            std::cout << std::endl << "ReadBucket #" << j << "    " << '"';
            std::cout << indexText(host(p))[j] << '"' << std::endl;
            std::cout << "  length:      " << sequenceLength(j, host(p)) << std::endl;
            bucketParams = &_swiftBucketParams(p, j++);
            _printSwiftParams(*bucketParams);
        }

        std::cout << "    lastInc: " << (int)p.buckets[i].lastIncrement;
        std::cout << "  \tCounter: " << p.buckets[i].counter << std::endl;
    }
}

template <typename TIndex, typename TSpec, typename TSize>
inline typename Pattern<TIndex, Swift<TSpec> >::TBucketParams &
_swiftBucketParams(Pattern<TIndex, Swift<TSpec> > & pattern, TSize seqNo)
{
    if (Swift<TSpec>::PARAMS_BY_LENGTH)
        return pattern.bucketParams[sequenceLength(seqNo, host(pattern))];
    else
        return pattern.bucketParams[seqNo];
}

template <typename TIndex, typename TSpec, typename TParams, typename TSize>
inline unsigned
_swiftBucketNo(Pattern<TIndex, Swift<TSpec> > const &, TParams &bucketParams, TSize seqNo)
{
    if (Swift<TSpec>::PARAMS_BY_LENGTH)
        return (bucketParams.reuseMask + 1) * seqNo;    // assumes the same reuseMask for all reads
    else
        return bucketParams.firstBucket;
}

template <typename TIndex, typename TSpec, typename TSeqNo>
inline int
_qgramLemma(Pattern<TIndex, Swift<TSpec> > const & pattern, TSeqNo seqNo, int errors)
{
    // q-gram lemma: How many conserved q-grams we see at least?
    // each error destroys at most <weight> many (gapped) q-grams
    return qgramThreshold(indexShape(host(pattern)), sequenceLength(seqNo, host(pattern)), errors, EditDistance(), ThreshQGramLemma());
}

/*!
 * @fn SwiftPattern#setMinThreshold
 * @headerfile <seqan/index.h>
 * @brief Set minimal threshold (<i>t</i>) for a given needle.
 *
 * This function can be used to make the filter allow less errors for a given needle.
 *
 * @signature void setMinThreshold(pattern, seqNo, thresh);
 *
 * @param[in,out] pattern The SwiftPattern to modify.
 * @param[in]     seqNo   The number to set the threshold for.
 * @param[in]     tresh   The threshold to use.
 */


template <typename TIndex, typename TSpec, typename TSeqNo, typename TThreshold>
inline void
setMinThreshold(Pattern<TIndex, Swift<TSpec> > & pattern, TSeqNo seqNo, TThreshold thresh)
{
    typedef Pattern<TIndex, Swift<TSpec> >                      TPattern;
    typedef typename TPattern::TBucketParams                    TBucketParams;
    typedef typename TPattern::TBucketString                    TBucketString;
    typedef typename Iterator<TBucketString, Standard>::Type    TBucketIterator;

    TBucketParams &bucketParams = _swiftBucketParams(pattern, seqNo);
    TBucketIterator it = begin(pattern.buckets, Standard()) + _swiftBucketNo(pattern, bucketParams, seqNo);
    TBucketIterator itEnd = it + (bucketParams.reuseMask + 1);

    for (; it != itEnd; ++it)
    {

        // increase the threshold if it is less than the minimal threshold
        if ((*it).threshold < thresh)
        {
            // increase the counter once it has reached the threshold
            // otherwise we could output the same hit multiple times
            if ((*it).counter >= (*it).threshold)
                (*it).counter = thresh;

            (*it).threshold = thresh;
        }
    }
}

template <typename TSpec, typename TSize, typename TShortSize, typename TPos>
inline void
_resetBucket(SwiftBucket_<TSpec, TSize, TShortSize> & bkt, TPos lastIncrement)
{
    bkt.lastIncrement = lastIncrement;
    bkt.counter = 0;
    bkt.notListed = true;
}

template <typename TSpec, typename TSize, typename TShortSize, typename TPos, typename TThresh>
inline void
_resetBucket(SwiftBucket_<TSpec, TSize, TShortSize> & bkt, TPos lastIncrement, TThresh threshold)
{
    bkt.lastIncrement = lastIncrement;
    bkt.counter = 0;
    bkt.threshold = threshold;
    bkt.notListed = true;
}

template <typename TSpec_, typename TSize, typename TShortSize, typename TPos>
inline void
_resetBucket(SwiftBucket_<SwiftSemiGlobal_<TSpec_>, TSize, TShortSize> & bkt, TPos lastIncrement)
{
    bkt.lastIncrement = lastIncrement;
    bkt.counter = 0;
}

template <typename TSpec_, typename TSize, typename TShortSize, typename TPos, typename TThresh>
inline void
_resetBucket(SwiftBucket_<SwiftSemiGlobal_<TSpec_>, TSize, TShortSize> & bkt, TPos lastIncrement, TThresh threshold)
{
    bkt.lastIncrement = lastIncrement;
    bkt.counter = 0;
    bkt.threshold = threshold;
}

template <typename TIndex, typename TFloat, typename TSize_, typename TSpec>
inline void _patternInit(Pattern<TIndex, Swift<TSpec> > &pattern, TFloat errorRate, TSize_ minLengthForAll)
{
    typedef Pattern<TIndex, Swift<TSpec> >                      TPattern;
    typedef typename Size<TIndex>::Type                         TSize;
    //typedef typename Fibre<TIndex, QGramSA>::Type               TSA;
    //typedef typename Iterator<TSA, Standard>::Type              TSAIter;
    typedef typename TPattern::TBucket                          TBucket;
    typedef typename TBucket::TSize                             TBucketSize;
    typedef typename TPattern::TBucketParams                    TBucketParams;
    typedef typename TPattern::TBucketString                    TBucketString;
    typedef typename Iterator<TBucketString, Standard>::Type    TBucketIterator;

    double _newErrorRate = errorRate;
    TSize seqCount = countSequences(host(pattern));

    clear(pattern.verifyList);

    if (pattern._currentErrorRate != _newErrorRate || pattern._currentMinLengthForAll != minLengthForAll)
    {
        // settings have been changed -> initialize bucket parameters

        pattern._currentErrorRate = _newErrorRate;
        pattern._currentMinLengthForAll = minLengthForAll;

        indexRequire(host(pattern), QGramSADir());
        pattern.shape = indexShape(host(pattern));

        TSize span = length(pattern.shape);
        TSize count = 0;
        TSize bucketsPerCol2Max = 0;
        TSize maxLength = 0;

        if (Swift<TSpec>::PARAMS_BY_LENGTH) {
            for(unsigned seqNo = 0; seqNo < seqCount; ++seqNo) {
                TSize length = sequenceLength(seqNo, host(pattern));
                if (maxLength < length)
                    maxLength = length;
            }
            resize(pattern.bucketParams, maxLength + 1);
        } else
            resize(pattern.bucketParams, seqCount);

        pattern.maxPatternLength = maxLength;
        pattern.finderPosOffset = 0;
        pattern.finderPosNextOffset = pattern.finderLength + pattern.maxPatternLength;

        if (Swift<TSpec>::SEMIGLOBAL == 0)
        {
            // global matches
            TSize minLength = minLengthForAll;
            for(unsigned seqNo = 0; seqNo < seqCount; ++seqNo)
            {
                // swift q-gram lemma
                TBucketParams &bucketParams = _swiftBucketParams(pattern, seqNo);
                // n..next length that could decrease threshold
                TSize n = (TSize) ceil((floor(errorRate * minLength) + 1) / errorRate);
                // minimal threshold is minimum errors of minLength and n
                int threshold = (TSize) _min(
                    (n + 1) - span * (floor(errorRate * n) + 1),
                    (minLength + 1) - span * (floor(errorRate * minLength) + 1));

                if (threshold > pattern.params.minThreshold)
                    bucketParams.threshold = threshold;
                else
                    bucketParams.threshold = pattern.params.minThreshold;

                SEQAN_ASSERT_GT_MSG((1 / errorRate), span, "SWIFT only works if span < 1 / error rate!");
                TSize errors = (TSize) floor((2 * bucketParams.threshold + span - 3) / (1 / errorRate - span));


                // a bucket has distanceCut different positions of q-grams
                // if a q-gram is this far or further away it can't belong to the
                // same bucket
                bucketParams.distanceCut = (bucketParams.threshold - 1) + span * errors;

                // from now, errors is the maximal number of indels
                if (Swift<TSpec>::HAMMING_ONLY != 0)
                    errors = 0;

                TSize bucketsPerCol2;
                if(Swift<TSpec>::DIAGONAL == 1)
                {
                    // Use overlapping parallelograms
                    bucketParams.overlap = errors;

                    // delta must be a power of 2 and greater than errors
                    bucketParams.logDelta = (TSize) ceil(log((double)errors + 1) / log(2.0));
                    if (bucketParams.logDelta < pattern.params.minLog2Delta)
                        bucketParams.logDelta = pattern.params.minLog2Delta;
                    bucketParams.delta = 1 << bucketParams.logDelta;
                    bucketParams.tabooLength = pattern.params.tabooLength;

                    // maximal number of buckets in one column
                    TSize bucketsPerCol = (sequenceLength(seqNo, host(pattern)) - span + 2 * bucketParams.delta + errors - 1) / bucketParams.delta;
                    bucketsPerCol2 = 1 << (TSize) ceil(log((double)bucketsPerCol) / log(2.0)); // next greater or equal power of 2
                }
                else
                {
                    // TODO: classical swift for rectangular buckets
                    // Use overlapping rectangles
                    //bucketParams.overlap = ;

                    // delta must be a power of 2 greater than seq.length + errors (define a minimal delta of 32)
                    //bucketParams.logDelta = ;
                    //if (bucketParams.logDelta < pattern.params.minLog2Delta)
                    //  bucketParams.logDelta = pattern.params.minLog2Delta;
                    //bucketParams.delta = 1 << bucketParams.logDelta;
                    bucketsPerCol2 = 1;
                }

                bucketParams.firstBucket = count; // firstBucket is only used if Swift<TSpec>::PARAMS_BY_LENGTH == 0
                bucketParams.reuseMask = bucketsPerCol2 - 1;
                bucketParams.tabooLength = pattern.params.tabooLength;

                if (Swift<TSpec>::PARAMS_BY_LENGTH) {
                    ++count;
                    if (bucketsPerCol2Max < bucketsPerCol2)
                        bucketsPerCol2Max = bucketsPerCol2;
                } else
                    count += bucketsPerCol2;

                /*if (seqNo<3)
                    _printSwiftParams(bucketParams);*/
            }
        } else
            for(unsigned seqNo = 0; seqNo < seqCount; ++seqNo)
            {
                // get pattern length and max. allowed errors
                TBucketParams &bucketParams = _swiftBucketParams(pattern, seqNo);
                TSize length = 0;
                if (minLengthForAll != static_cast<TSize_>(0))
                    length = minLengthForAll;
                else
                    length = sequenceLength(seqNo, host(pattern));
                TSize errors = (TSize) floor(errorRate * length);
                TSize errorsWC = errors / (1 + Swift<TSpec>::QGRAM_ERRORS);

                // q-gram lemma: How many conserved q-grams we see at least?
                // (define a minimal threshold of 1)
                int threshold = length - span + 1 - errorsWC * weight(pattern.shape);
                if (threshold > pattern.params.minThreshold)
                    bucketParams.threshold = threshold;
                else
                    bucketParams.threshold = pattern.params.minThreshold;

                // from now, errors is the maximal number of indels
                if (Swift<TSpec>::HAMMING_ONLY != 0)
                    errors = 0;

                // a bucket has distanceCut different positions of q-grams
                // if a q-gram is this far or farer away it can't belong to the
                // same bucket
    //          bucketParams.distanceCut = length - (span - 1) + errors;

                TSize bucketsPerCol2;
                if (Swift<TSpec>::DIAGONAL == 1)
                {
                    // Use overlapping parallelograms
                    bucketParams.overlap = errors;

                    // delta must be a power of 2 greater then errors (define a minimal delta of 8)
                    bucketParams.logDelta = (TSize) ceil(log((double)(errors + 1)) / log(2.0));
                    if (bucketParams.logDelta < pattern.params.minLog2Delta)
                        bucketParams.logDelta = pattern.params.minLog2Delta;
                    bucketParams.delta = 1 << bucketParams.logDelta;

                    // the formula for bucketsPerCol is (worst-case):
                    // (height-(q-1) - 1 - (delta+1-e))/delta + 3
                    //    ^-- full paral. in the middle --^     ^-- 2 at the bottom, 1 at the top
                    TSize bucketsPerCol = (sequenceLength(seqNo, host(pattern)) - span + 2 * bucketParams.delta + errors - 1) / bucketParams.delta;
                    bucketsPerCol2 = 1 << (TSize) ceil(log((double)bucketsPerCol) / log(2.0));
                }
                else
                {
                    // Use overlapping rectangles
                    bucketParams.overlap = length - span + errors;

                    // delta must be a power of 2 greater then seq.length + errors (define a minimal delta of 32)
                    bucketParams.logDelta = (TSize) ceil(log((double)(length - span + 1 + errors)) / log(2.0));
                    if (bucketParams.logDelta < pattern.params.minLog2Delta)
                        bucketParams.logDelta = pattern.params.minLog2Delta;
                    bucketParams.delta = 1 << bucketParams.logDelta;

                    bucketsPerCol2 = 2;
                }

    //          SEQAN_ASSERT_LEQ(distanceCut, bucketsPerCol * (TSize) delta);

                bucketParams.firstBucket = count;
                bucketParams.reuseMask = bucketsPerCol2 - 1;
                bucketParams.tabooLength = pattern.params.tabooLength;

                if (Swift<TSpec>::PARAMS_BY_LENGTH) {
                    ++count;
                    if (bucketsPerCol2Max < bucketsPerCol2)
                        bucketsPerCol2Max = bucketsPerCol2;
                } else
                    count += bucketsPerCol2;

/*              if (seqNo<3)
                    _printSwiftParams(bucketParams);
*/          }

        if (Swift<TSpec>::PARAMS_BY_LENGTH) {
            count *= bucketsPerCol2Max;
            for(unsigned i = 0; i < length(pattern.bucketParams); ++i)
                pattern.bucketParams[i].reuseMask = bucketsPerCol2Max - 1;
        }
        resize(pattern.buckets, count);

        TBucketIterator bktEnd, bkt = begin(pattern.buckets, Standard());
        for(unsigned seqNo = 0; seqNo < seqCount; ++seqNo)
        {
            TBucketParams &bucketParams = _swiftBucketParams(pattern, seqNo);
            TBucketSize lastIncrement = (TBucketSize)0 - (TBucketSize)bucketParams.tabooLength;
            for(bktEnd = bkt + bucketParams.reuseMask + 1; bkt != bktEnd; ++bkt)
            {
                _resetBucket(*bkt, lastIncrement, bucketParams.threshold);
            }
        }
    }
    else
    {
        // settings are unchanged -> reset buckets

        // finderPosOffset is used to circumvent expensive resetting of all buckets
        __int64 clearance = pattern.finderLength + pattern.maxPatternLength;
        pattern.finderPosOffset = pattern.finderPosNextOffset;
        pattern.finderPosNextOffset += clearance;

        // reset buckets only if an overflow of the finder position would occur
        // we would not detect an overflow if clearance is larger than the largest TBucketSize value
        if (pattern.finderPosNextOffset <= pattern.finderPosOffset || (__int64)(TBucketSize)clearance < (__int64)clearance)
        {
            pattern.finderPosOffset = 0;
            pattern.finderPosNextOffset = pattern.finderLength + pattern.maxPatternLength;

            TBucketIterator bktEnd, bkt = begin(pattern.buckets, Standard());
            for (TSize ndlSeqNo = 0; ndlSeqNo < seqCount; ++ndlSeqNo)
            {
                TBucketParams &bucketParams = _swiftBucketParams(pattern, ndlSeqNo);
                TBucketSize lastIncrement = (TBucketSize)0 - (TBucketSize)bucketParams.tabooLength;
                for (bktEnd = bkt + (bucketParams.reuseMask + 1); bkt != bktEnd; ++bkt)
                {
                    _resetBucket(*bkt, lastIncrement);
                }
            }
        }
    }

/*
    std::cerr << "Swift bucket params: " << length(pattern.bucketParams) << std::endl;
    std::cerr << "Swift buckets:       " << length(pattern.buckets) << std::endl;
    std::cerr << "Buckets per read:    " << bucketsPerCol2Max << std::endl;
*/
}


/////////////////////////////////////////////////////////////
// Creates a new hit and appends it to the finders hit list
template <
    typename THaystack,
    typename TIndex,
    typename TSpec,
    typename TBucket,
    typename TBucketParams,
    typename TSize
>
inline void _createHit(
    Finder<THaystack, Swift<TSpec> > & finder,
    Pattern<TIndex, Swift<TSpec> > & pattern,
    TBucket & bkt,
    TBucketParams & bucketParams,
    __int64 diag,
    TSize ndlSeqNo)
{
    typedef typename FindResult<Finder<THaystack, Swift<TSpec> >, Pattern<TIndex, Swift<TSpec> > >::Type THit;
    __int64 lastInc = (__int64)bkt.lastIncrement - pattern.finderPosOffset;
    __int64 firstInc = (__int64)bkt.firstIncrement - pattern.finderPosOffset;

    if(diag > lastInc)
    {
        // bucket is reused since last increment
        TSize reusePos = (bucketParams.reuseMask + 1) << bucketParams.logDelta;
        diag -= (__int64)ceil((diag-lastInc)/(double)reusePos) * reusePos;
    }

    // determine width, height, and begin position in needle
    TSize width = lastInc - firstInc + length(pattern.shape);
    TSize height = width + bucketParams.delta + bucketParams.overlap;
    __int64 ndlBegin = lastInc + length(pattern.shape) - diag - height;

    // create the hit
    THit hit = {                //                              *
        firstInc,               // bucket begin in haystack     * *
        ndlSeqNo,               // needle seq. number           *   *
        ndlBegin,               // bucket begin in needle       *     *
        width,                  // bucket width (non-diagonal)    *   *
        height                  // bucket height                    * *
    };                          //                                    *

    // append the hit to the finders hit list
    appendValue(finder.hits, hit);
}

//////////////////////////////////////////////////////////////////////
// Updates the counters of the buckets in which the q-gram with hash value hash occurs.
// Assures that those updated bucket counters are set to one
//    - that exceeded the reuse mask since last increment
//    - for which the last increment lies more than distanceCut away.
// If a bucket counter reaches threshold a hit is appended to the hit-list of the finder.
// Returns true if the hit-list of the finder is not empty after this call.
template <
    typename TFinder,
    typename TIndex,
    typename TSpec,
    typename THashValue
>
inline bool _swiftMultiProcessQGram(
    TFinder & finder,
    Pattern<TIndex, Swift<TSpec> > & pattern,
    THashValue hash)
{
    typedef Pattern<TIndex, Swift<TSpec> >                      TPattern;

    //typedef typename Size<TIndex>::Type                         TSize;
    typedef typename Fibre<TIndex const, QGramSA>::Type         TSA;
    typedef typename Iterator<TSA, Standard>::Type              TSAIter;
    typedef typename TPattern::TBucketString                    TBucketString;
    typedef typename Iterator<TBucketString, Standard>::Type    TBucketIter;
    typedef typename Value<TBucketString>::Type                 TBucket;
    typedef typename TBucket::TShortSize                        TShortSize;
    typedef typename TPattern::TBucketParams                    TBucketParams;
    //typedef typename FindResult<TFinder, TPattern>::Type        THit;

    TIndex const &index = host(pattern);

    // create an iterator over the positions of the q-gram occurrences in pattern
    TSAIter saBegin = begin(indexSA(index), Standard());
    TSAIter occ = saBegin + indexDir(index)[getBucket(index.bucketMap, hash)];
    TSAIter occEnd = saBegin + indexDir(index)[getBucket(index.bucketMap, hash) + 1];
    TBucketIter bktBegin = begin(pattern.buckets, Standard());
    Pair<unsigned> ndlPos;

/*  std::cerr<<"\t["<<(occEnd-occ)<<"]"<< std::flush;

    if ((occEnd-occ)>100)
    {
        std::cerr<<" ";
        for(int i=0;i<length(indexShape(host(pattern)));++i)
            std::cerr<<*(hostIterator(hostIterator(finder))+i);
    }
*/

    // iterate over all q-gram occurrences and do the processing
    __int64 curPos = finder.curPos + pattern.finderPosOffset;
    for(; occ != occEnd; ++occ)
    {
        posLocalize(ndlPos, *occ, stringSetLimits(index)); // get pair of SeqNo and Pos in needle
        TBucketParams &bucketParams = _swiftBucketParams(pattern, getSeqNo(ndlPos));

        // begin position of the diagonal of q-gram occurrence in haystack (possibly negative)
        __int64 diag = finder.curPos;
        if (Swift<TSpec>::DIAGONAL == 1) diag -= getSeqOffset(ndlPos);

        unsigned bktNo = (diag >> bucketParams.logDelta) & bucketParams.reuseMask; // bucket no of diagonal
        unsigned bktOfs = diag & (bucketParams.delta - 1); // offset of diagonal to bucket begin
        __int64  bktBeginHstk = diag & ~(__int64)(bucketParams.delta - 1); // haystack position of bucket begin diagonal

        // global (over all pattern sequences) number of current bucket
        TBucketIter bkt = bktBegin + (_swiftBucketNo(pattern, bucketParams, getSeqNo(ndlPos)) + bktNo);

        TShortSize hitCount;

        do {
            if ((__int64)(*bkt).lastIncrement < bktBeginHstk + (__int64)pattern.finderPosOffset
                || (__int64)((*bkt).lastIncrement + bucketParams.distanceCut) < curPos)
            {
                // last increment was before the beginning of the current bucket => bucket is reused
                // Or last increment was in the same bucket but lies more than distanceCut away

                if ((*bkt).counter >= (*bkt).threshold)
                {
                    // create a new hit and append it to the finders hit list
                    _createHit(finder, pattern, *bkt, bucketParams, bktBeginHstk, getSeqNo(ndlPos));
                }

                // reuse bucket
                hitCount = 1;
                (*bkt).firstIncrement = curPos;
            }
            else if((__int64)((*bkt).lastIncrement + bucketParams.tabooLength) > curPos)
            {
                // bkt counter was already incremented for another q-gram at
                //   a haystack position that is closer than tabooLength
                // we jump directly to
                //   where we check whether the q-gram falls into another overlapping bucket or not
                goto checkOverlap;
            }
            else
            {
                if((*bkt).counter == 0) (*bkt).firstIncrement = curPos;
                hitCount = (*bkt).counter + 1;
            }

            (*bkt).lastIncrement = curPos;
            (*bkt).counter = hitCount;
#ifdef SEQAN_DEBUG_SWIFT
            (*bkt)._lastIncDiag = diag;
#endif

            if (hitCount == (*bkt).threshold && (*bkt).notListed) {
                // append bkt no to patterns verify list
                appendValue(pattern.verifyList, Pair<unsigned>(getSeqNo(ndlPos), bktNo));
                (*bkt).notListed = false;
            }

checkOverlap:
            // check if q-gram falls into another overlapping bucket
            if(bktOfs >= bucketParams.overlap) break;

            // set to previous overlapping bucket for next iteration
            bktBeginHstk -= bucketParams.delta;
            bktOfs += bucketParams.delta;
            if(bktNo) {
                --bktNo;
                --bkt;
            } else {
                bktNo = bucketParams.reuseMask;
                bkt += bktNo;
            }
        }
        while(true);
    }

    finder.curHit = begin(finder.hits, Standard());
    finder.endHit = end(finder.hits, Standard());

    return !empty(finder.hits);
}

///////////////////////////////////////////////////////////////////
// Updates the counters of the buckets in which the q-gram with hash value hash occurs.
// Assures that those updated bucket counters that exceeded the reuse mask since last increment are set to one.
// If a bucket counter reaches threshold a hit is appended to the hit-list of the finder.
// Returns true if the hit-list of the finder is not empty after this call.
template <
    typename TFinder,
    typename TIndex,
    typename TSpec_,
    typename THValue
>
inline bool _swiftMultiProcessQGram(
    TFinder &finder,
    Pattern<TIndex, Swift<Tag<SwiftSemiGlobal_<TSpec_> > > > &pattern,
    THValue hash)
{
    typedef Pattern<TIndex, Swift<Tag<SwiftSemiGlobal_<TSpec_> > > >    TPattern;

    typedef typename Size<TIndex>::Type                         TSize;
    typedef typename Fibre<TIndex, QGramSA>::Type               TSA;
    typedef typename Iterator<TSA const, Standard>::Type        TSAIter;
    typedef typename TPattern::TBucketString                    TBucketString;
    typedef typename Iterator<TBucketString, Standard>::Type    TBucketIter;
    typedef typename Value<TBucketString>::Type                 TBucket;
    typedef typename TBucket::TShortSize                        TShortSize;
    typedef typename TPattern::TBucketParams                    TBucketParams;
    typedef typename FindResult<TFinder, TPattern>::Type        THit;

    TIndex const & index = host(pattern);

    // create an iterator over the positions of the q-gram occurrences in pattern
    TSAIter saBegin = begin(indexSA(index), Standard());
    TSAIter occ = saBegin + indexDir(index)[getBucket(index.bucketMap, hash)];
    TSAIter occEnd = saBegin + indexDir(index)[getBucket(index.bucketMap, hash) + 1];
    TBucketIter bktBegin = begin(pattern.buckets, Standard());
    Pair<unsigned> ndlPos;

/*  std::cerr<<"\t["<<(occEnd-occ)<<"]"<< std::flush;

    if ((occEnd-occ)>100)
    {
        std::cerr<<" ";
        for(int i=0;i<length(indexShape(host(pattern)));++i)
            std::cerr<<*(hostIterator(hostIterator(finder))+i);
    }
*/
    // iterate over all q-gram occurrences and do the processing
    __int64 curPos = finder.curPos + pattern.finderPosOffset;
//    bool dbg=pattern.params.debug && (finder.curPos > 496 && finder.curPos < 591);
    for(; occ != occEnd; ++occ)
    {
        posLocalize(ndlPos, *occ, stringSetLimits(index));
        TBucketParams &bucketParams = _swiftBucketParams(pattern, getSeqNo(ndlPos));

        __int64 diag = finder.curPos;
        if (Swift<Tag<SwiftSemiGlobal_<TSpec_> > >::DIAGONAL == 1) diag -= getSeqOffset(ndlPos);

        unsigned bktNo = (diag >> bucketParams.logDelta) & bucketParams.reuseMask;
        unsigned bktOfs = diag & (bucketParams.delta - 1);
        __int64  bktBeginHstk = diag & ~(__int64)(bucketParams.delta - 1);

        TBucketIter bkt = bktBegin + (_swiftBucketNo(pattern, bucketParams, getSeqNo(ndlPos)) + bktNo);
        TShortSize hitCount;

//        if (dbg && (*occ).i1==15)
//            std::cout << finder.curPos << '\t' << *occ << '\t' << (*bkt).counter<<std::endl;
//
        do
        {
            if ((__int64)(*bkt).lastIncrement < bktBeginHstk + (__int64)pattern.finderPosOffset)
            {
                // last increment was before the beginning of the current bucket
                // (we must ensure that bucketIdx doesn't collide)
                hitCount = 1;
            }
            else
            {
                if ((__int64)((*bkt).lastIncrement + bucketParams.tabooLength) > curPos)
                    goto checkOverlap;  // increment only once per sequence
                hitCount = (*bkt).counter + 1;
            }

            (*bkt).lastIncrement = curPos;
            (*bkt).counter = hitCount;
#ifdef SEQAN_DEBUG_SWIFT
            (*bkt)._lastIncDiag = diag;
#endif

            if (hitCount == (*bkt).threshold)
            {

                TSize height = 0;
                if (Swift<Tag<SwiftSemiGlobal_<TSpec_> > >::DIAGONAL == 1)
                    height = sequenceLength(getSeqNo(ndlPos), host(pattern)) - 1;

#ifdef SEQAN_DEBUG_SWIFT
                // upper bucket no. of lastIncr. q-gram
                __int64 upperBktNo = ((*bkt).lastIncrement - pattern.finderPosOffset) >> bucketParams.logDelta;

                // we must decrement bucket no. until (no. mod reuse == bktNo)
                __int64 _bktBeginHstk =
                     (upperBktNo - ((upperBktNo - bktNo) & bucketParams.reuseMask)) << bucketParams.logDelta;

                if ((*bkt)._lastIncDiag - _bktBeginHstk >= bucketParams.delta + bucketParams.overlap || (*bkt)._lastIncDiag < _bktBeginHstk) {
                    std::cerr << "qgram stored in wrong bucket (diag:" << (*bkt)._lastIncDiag << ", begin:" << _bktBeginHstk;
                    std::cerr << ", delta:" << bucketParams.delta << ", overlap:" << bucketParams.overlap << ")" << std::endl;
                }
#endif
//              if (bktBeginHstk >= 0)
//              {
                    THit hit = {
                        bktBeginHstk,                                       // bucket begin in haystack
                        getSeqNo(ndlPos),                                   // needle seq. number
                        static_cast<unsigned>(height + bucketParams.delta +
                                              bucketParams.overlap)         // bucket width (non-diagonal)
                    };
                    appendValue(finder.hits, hit);
//              } else {
//                  // match begins left of haystack begin
//                  THit hit = {
//                      0,                                                  // bucket begin in haystack
//                      getSeqNo(ndlPos),                                   // needle seq. number
//                      height + bucketParams.delta + bucketParams.overlap  // bucket width (non-diagonal)
//                      + (diag & ~(__int64)(bucketParams.delta - 1))
//                  };
//                  appendValue(finder.hits, hit);
//              }
            }

        checkOverlap:
            if (bktOfs >= bucketParams.overlap) break;

            // repeat with the previous overlapping bucket
            bktBeginHstk -= bucketParams.delta;
            bktOfs += bucketParams.delta;
            if (bktNo) {
                --bktNo;
                --bkt;
            } else {
                bktNo = bucketParams.reuseMask;
                bkt += bktNo;
            }
        } while (true);
    }

    finder.curHit = begin(finder.hits, Standard());
    finder.endHit = end(finder.hits, Standard());

    return !empty(finder.hits);
}

////////////////////////////////////////////////////////////////////////////////////
// resets counter and lastIncrement of all buckets listed in patterns verify list
template <
    typename TFinder,
    typename TIndex,
    typename TSpec
>
inline bool _swiftMultiFlushBuckets(
    TFinder & finder,
    Pattern<TIndex, Swift<TSpec> > & pattern
    )
{
    typedef Pattern<TIndex, Swift<TSpec> >                      TPattern;

    typedef typename TPattern::TBucket                          TBucket;
    typedef typename TBucket::TSize                             TBucketSize;
    typedef typename TPattern::TBucketString                    TBucketString;
    typedef typename Iterator<TBucketString, Standard>::Type    TBucketIterator;
    typedef typename TPattern::TBucketParams                    TBucketParams;

    typedef typename TPattern::TVerifyList                      TVerifyList;
    typedef typename Iterator<TVerifyList, Standard>::Type      TListIterator;

    typedef typename Size<TIndex>::Type                         TSize;

    __int64 hstkLength = length(haystack(finder));

    TListIterator verifyBkt = begin(pattern.verifyList, Standard());
    TListIterator verifyListEnd = end(pattern.verifyList, Standard());
    for (; verifyBkt < verifyListEnd; ++verifyBkt)
    {
        unsigned bktNo = (*verifyBkt).i2;
        unsigned ndlSeqNo = (*verifyBkt).i1;
        TBucketParams &bucketParams = _swiftBucketParams(pattern, ndlSeqNo);

        TBucketIterator bkt = begin(pattern.buckets, Standard()) + _swiftBucketNo(pattern, bucketParams, ndlSeqNo) + bktNo;
        if ((*bkt).counter >= (*bkt).threshold)
        {
            // hstkPos / delta: gives the number of the bucket that is at the top of this column (modulo reuseMask missing)
            TSize topBucket = (TSize)(hstkLength >> bucketParams.logDelta);
            // number of buckets in last column above the bucket with the number bktNo
            TSize bucketNoInCol = (topBucket + bucketParams.reuseMask + 1 - bktNo) & bucketParams.reuseMask;
            // begin position of lower diagonal of this bucket in haystack (possibly negative)
            __int64 diag = (hstkLength & ~(__int64)(bucketParams.delta - 1)) - (bucketNoInCol << bucketParams.logDelta);

            // create a new hit and append it to the finders hit list
            _createHit(finder, pattern, *bkt, bucketParams, diag, ndlSeqNo);
        }
        _resetBucket(*bkt, (TBucketSize)0 - (TBucketSize)bucketParams.tabooLength);
    }

    clear(pattern.verifyList);

    finder.curHit = begin(finder.hits, Standard());
    finder.endHit = end(finder.hits, Standard());

    return !empty(finder.hits);
}

//////////////////////////////////////////////////////
// no resetting is needed for the semiglobal version
template <
    typename TFinder,
    typename TIndex,
    typename TSpec_
>
inline bool _swiftMultiFlushBuckets(
    TFinder &,
    Pattern<TIndex, Swift<Tag<SwiftSemiGlobal_<TSpec_> > > > &)
{
    // there is nothing to be done here as we dump matches immediately after reaching the threshold
    return false;
}

template <typename TIndex, typename TSpec>
inline bool
empty(Pattern<TIndex, Swift<TSpec> > & me)
{
    return empty(me.bucketParams);
}

template <typename TIndex, typename TSpec>
inline void
clear(Pattern<TIndex, Swift<TSpec> > & me)
{
    me.finderPosOffset = 0;
    me.finderPosNextOffset = 0;
    me.finderLength = 0;
    me.maxPatternLength = 0;
    me._currentErrorRate = -1;
    me._currentMinLengthForAll = -1;
    me.curSeqNo = 0;
    me.curBeginPos = 0;
    me.curEndPos = 0;
    clear(me.bucketParams);
    clear(me.buckets);
}

//____________________________________________________________________________

template <typename THaystack, typename TSpec>
inline typename Position<Finder<THaystack, Swift<TSpec> > >::Type
position(Finder<THaystack, Swift<TSpec> > const & finder)
{
    typename Finder<THaystack, Swift<TSpec> >::TSwiftHit &hit = *finder.curHit;
    return hit.hstkPos + hit.bucketWidth;
}

template <typename THaystack, typename TSpec>
inline typename Position<Finder<THaystack, Swift<TSpec> > >::Type
position(Finder<THaystack, Swift<TSpec> > & finder)
{
    return position(const_cast<Finder<THaystack, Swift<TSpec> > const &>(finder));
}

//____________________________________________________________________________

template <typename TIndex, typename TSpec>
inline typename SAValue<TIndex>::Type
position(Pattern<TIndex, Swift<TSpec> > const & pattern)
{
    __int64 hitEnd = pattern.curEndPos;
    __int64 textLength = sequenceLength(pattern.curSeqNo, needle(pattern));
    if(hitEnd > textLength) hitEnd = textLength;

    typename SAValue<TIndex >::Type pos;
    posLocalToX(pos, Pair<unsigned, __int64>(pattern.curSeqNo, hitEnd), stringSetLimits(host(pattern)));
    return pos;
}

template <typename TIndex, typename TSpec>
inline typename SAValue<TIndex>::Type
position(Pattern<TIndex, Swift<Tag<SwiftSemiGlobal_<TSpec> > > > const & pattern)
{
    typedef typename Size<TIndex>::Type TSize;
    typename SAValue<TIndex >::Type pos;
    posLocalToX(pos, Pair<unsigned, TSize>(pattern.curSeqNo, length(needle(pattern))), stringSetLimits(host(pattern)));
    return pos;
}

template <typename TIndex, typename TSpec>
inline typename SAValue<TIndex>::Type
position(Pattern<TIndex, Swift<TSpec> > & pattern)
{
    return position(const_cast<Pattern<TIndex, Swift<TSpec> > const &>(pattern));
}

//____________________________________________________________________________

template <typename THaystack, typename TSpec>
inline __int64
beginPosition(Finder<THaystack, Swift<TSpec> > const & finder)
{
    return (*finder.curHit).hstkPos;
}

template <typename THaystack, typename TSpec>
inline __int64
beginPosition(Finder<THaystack, Swift<TSpec> > & finder)
{
    return beginPosition(const_cast<Finder<THaystack, Swift<TSpec> > const &>(finder));
}

//____________________________________________________________________________

/*!
 * @fn SwiftPattern#beginPosition
 * @headerfile <seqan/index.h>
 * @brief Returns the begin position of the local match in the pattern.
 *
 * @signature TPosition beginPosition(pattern);
 *
 * @param[in] pattern   The SwiftPattern to query.
 * @return    TPosition The @link TextConcept#SAValue SA value @endlink of the pattern text.
 */

template <typename TIndex, typename TSpec>
inline typename SAValue<TIndex>::Type
beginPosition(Pattern<TIndex, Swift<TSpec> > const & pattern)
{
    __int64 hitBegin = pattern.curBeginPos;
    if (hitBegin < 0) hitBegin = 0;

    typename SAValue<TIndex >::Type pos;
    posLocalToX(pos, Pair<unsigned, __int64>(pattern.curSeqNo, hitBegin), stringSetLimits(host(pattern)));
    return pos;
}

template <typename TIndex, typename TSpec>
inline typename SAValue<TIndex >::Type
beginPosition(Pattern<TIndex, Swift<Tag<SwiftSemiGlobal_<TSpec> > > > const & pattern)
{
    typename SAValue<TIndex >::Type pos;
    posLocalToX(pos, Pair<unsigned>(pattern.curSeqNo, 0), stringSetLimits(host(pattern)));
    return pos;
}

template <typename TIndex, typename TSpec>
inline typename SAValue<TIndex>::Type
beginPosition(Pattern<TIndex, Swift<TSpec> > & pattern)
{
    return beginPosition(const_cast<Pattern<TIndex, Swift<TSpec> > const &>(pattern));
}

//____________________________________________________________________________

/*!
 * @fn SwiftPattern#endPosition
 * @headerfile <seqan/index.h>
 * @brief Returns the end position of the local match in the pattern.
 *
 * @signature TPosition endPosition(pattern);
 *
 * @param[in] pattern   The SwiftPattern to query.
 * @return    TPosition The @link TextConcept#SAValue SA value @endlink of the pattern text.
 */

template <typename THaystack, typename TSpec>
inline typename Position<Finder<THaystack, Swift<TSpec> > >::Type
endPosition(Finder<THaystack, Swift<TSpec> > const & finder)
{
    typename Finder<THaystack, Swift<TSpec> >::TSwiftHit &hit = *finder.curHit;
    return hit.hstkPos + hit.bucketWidth;
}

template <typename THaystack, typename TSpec>
inline typename Position<Finder<THaystack, Swift<TSpec> > >::Type
endPosition(Finder<THaystack, Swift<TSpec> > & finder)
{
    return endPosition(const_cast<Finder<THaystack, Swift<TSpec> > const &>(finder));
}

//____________________________________________________________________________

template <typename TIndex, typename TSpec>
inline typename SAValue<TIndex>::Type
endPosition(Pattern<TIndex, Swift<TSpec> > const & pattern)
{
    __int64 hitEnd = pattern.curEndPos;
    __int64 textLength = sequenceLength(pattern.curSeqNo, needle(pattern));
    if(hitEnd > textLength) hitEnd = textLength;

    typename SAValue<TIndex >::Type pos;
    posLocalToX(pos, Pair<unsigned, __int64>(pattern.curSeqNo, hitEnd), stringSetLimits(host(pattern)));
    return pos;
}

template <typename TIndex, typename TSpec>
inline typename SAValue<TIndex >::Type
endPosition(Pattern<TIndex, Swift<Tag<SwiftSemiGlobal_<TSpec> > > > const & pattern)
{
    typedef typename Size<TIndex>::Type TSize;
    typename SAValue<TIndex >::Type pos;
    posLocalToX(pos, Pair<unsigned, TSize>(pattern.curSeqNo, length(needle(pattern))), stringSetLimits(host(pattern)));
    return pos;
}

template <typename TIndex, typename TSpec>
inline typename SAValue<TIndex>::Type
endPosition(Pattern<TIndex, Swift<TSpec> > & pattern)
{
    return endPosition(const_cast<Pattern<TIndex, Swift<TSpec> > const &>(pattern));
}

//____________________________________________________________________________
/*!
 * @fn SwiftFinder#positionRangeNoClip
 * @headerfile <seqan/index.h>
 * @brief Returns a pair of the begin and end position in or beyond the haystack or needle for the last hit found.
 *
 * @signature TPair positionRangeNoClip(finder);
 *
 * @param[in] finder A SwiftFinder object to query.
 *
 * @return TPair A pair of the begin and end position in the haystack or needle for the last hit found.  These positions
 *               could be negative or beyond the end of <tt>finder</tt> or <tt>pattern</tt> when using filter
 *               algorithms.  The return type is <tt>Pair&lt;typename SAValue&lt;THost&gt;::Type&gt;</tt> if
 *               <tt>THost</tt> is the type of haystack or needle.
 *
 * @see SwiftFinder#positionRange
 * @see SwiftPattern#positionRangeNoClip
 */

/*!
 * @fn SwiftPattern#positionRangeNoClip
 * @headerfile <seqan/index.h>
 * @brief Returns a pair of the begin and end position in or beyond the haystack or needle for the last hit found.
 *
 * @signature TPair positionRangeNoClip(pattern)
 *
 * @param[in] pattern SwiftPattern object to query.
 *
 * @return TPair A pair of the begin and end position in the haystack or needle for the last hit found.  These positions
 *               could be negative or beyond the end of <tt>finder</tt> or <tt>pattern</tt> when using filter
 *               algorithms.  The return type is <tt>Pair&lt;typename SAValue&lt;THost&gt;::Type&gt;</tt> if
 *               <tt>THost</tt> is the type of haystack or needle.
 *
 * @see SwiftPattern#positionRange
 * @see SwiftFinder#positionRangeNoClip
 */

template <typename THaystack, typename TSpec>
inline Pair<typename Position<Finder<THaystack, Swift<TSpec> > >::Type>
positionRangeNoClip(Finder<THaystack, Swift<TSpec> > const & finder)
{
    typedef typename Position<Finder<THaystack, Swift<TSpec> > >::Type TPosition;
    typedef Pair<TPosition> TPair;
    typename Finder<THaystack, Swift<TSpec> >::TSwiftHit &hit = *finder.curHit;
    return TPair((TPosition)hit.hstkPos, (TPosition)(hit.hstkPos + hit.bucketWidth));
}

template <typename THaystack, typename TSpec>
inline Pair<typename Position<Finder<THaystack, Swift<TSpec> > >::Type>
positionRangeNoClip(Finder<THaystack, Swift<TSpec> > & finder)
{
    return positionRangeNoClip(const_cast<Finder<THaystack, Swift<TSpec> > const &>(finder));
}

//____________________________________________________________________________
/*!
 * @fn SwiftFinder#positionRange
 * @headerfile <seqan/index.h>
 * @brief Returns a pair of the begin and end position in the haystack or needle for the last hit found.
 *
 * @signature TPair positionRange(finder);
 *
 * @param[in] finder A @link Finder @endlink object.
 *
 * @return TPair A @link Pair @endlink of the begin and end position in the haystack or needle for the last
 *               hit found.  The return type is <tt>Pair&lt;typename SAValue&ltTHost&g;::Type&gt;</tt> if
 *               <tt>THost</tt> is the type of haystack or needle.
 *
 * @see Finder#beginPosition
 * @see Finder#endPosition
 * @see SwiftFinder#positionRangeNoClip
 */

/*!
 * @fn SwiftPattern#positionRange
 * @headerfile <seqan/index.h>
 * @brief Returns a pair of the begin and end position in the haystack or needle for the last hit found.
 *
 * @signature TPair positionRange(pattern);
 *
 * @param[in] pattern A @link Pattern @endlink object.
 *
 * @return TPair A @link Pair @endlink of the begin and end position in the haystack or needle for the last
 *               hit found.  The return type is <tt>Pair&lt;typename SAValue&ltTHost&g;::Type&gt;</tt> if
 *               <tt>THost</tt> is the type of haystack or needle.
 *
 * @see SwiftPattern#positionRangeNoClip
 */

template <typename THaystack, typename TSpec>
inline Pair<typename Position<Finder<THaystack, Swift<TSpec> > >::Type>
positionRange(Finder<THaystack, Swift<TSpec> > const & finder)
{
    typedef typename Position<Finder<THaystack, Swift<TSpec> > >::Type TPosition;
    typedef Pair<TPosition> TPair;
    typename Finder<THaystack, Swift<TSpec> >::TSwiftHit &hit = *finder.curHit;

    __int64 hitBegin = hit.hstkPos;
    __int64 hitEnd = hit.hstkPos + hit.bucketWidth;
    __int64 textEnd = length(haystack(finder));

    if (hitBegin < 0) hitBegin = 0;
    if (hitEnd > textEnd) hitEnd = textEnd;
    return TPair((TPosition)hitBegin, (TPosition)hitEnd);
}

template <typename THaystack, typename TSpec>
inline Pair<typename Position<Finder<THaystack, Swift<TSpec> > >::Type>
positionRange(Finder<THaystack, Swift<TSpec> > & finder)
{
    return positionRange(const_cast<Finder<THaystack, Swift<TSpec> > const &>(finder));
}

//____________________________________________________________________________

template <typename TIndex, typename TSpec>
inline Pair<typename SAValue<TIndex>::Type>
positionRange(Pattern<TIndex, Swift<TSpec> > & pattern)
{
    return Pair<typename SAValue<TIndex>::Type> (beginPosition(pattern), endPosition(pattern));
}

//____________________________________________________________________________

// TODO(holtgrew): Document this properly.

template <typename TSwiftHit, typename TText>
inline typename Infix<TText>::Type
swiftInfixNoClip(TSwiftHit const &hit, TText &text)
{
    return infix(text, hit.hstkPos, hit.hstkPos + hit.bucketWidth);
}

template <typename TSwiftHit, typename TText>
inline typename Infix<TText>::Type
swiftInfix(TSwiftHit const & hit, TText & text)
{
    __int64 hitBegin = hit.hstkPos;
    __int64 hitEnd = hit.hstkPos + hit.bucketWidth;
    __int64 textEnd = length(text);

    if (hitBegin < 0) hitBegin = 0;
    if (hitEnd > textEnd) hitEnd = textEnd;
    SEQAN_ASSERT_LEQ(hitBegin, hitEnd);
    return infix(text, hitBegin, hitEnd);
}

//____________________________________________________________________________

// now in find_base.h
template <typename THaystack, typename TSpec>
inline typename Infix<THaystack>::Type
infix(Finder<THaystack, Swift<TSpec> > &finder)
{
    typename Parameter_<THaystack>::Type tmpHaystack = haystack(finder);
    return swiftInfix(*finder.curHit, tmpHaystack);
}

template <typename THaystack, typename TSpec, typename TText>
inline typename Infix<TText>::Type
infix(Finder<THaystack, Swift<TSpec> > &finder, TText &text)
{
    return swiftInfix(*finder.curHit, text);
}

//____________________________________________________________________________

template <typename THaystack, typename TSpec>
inline typename Infix<THaystack>::Type
infixNoClip(Finder<THaystack, Swift<TSpec> > &finder)
{
    return swiftInfixNoClip(*finder.curHit, haystack(finder));
}

template <typename THaystack, typename TSpec, typename TText>
inline typename Infix<TText>::Type
infixNoClip(Finder<THaystack, Swift<TSpec> > &finder, TText &text)
{
    return swiftInfixNoClip(*finder.curHit, text);
}

//____________________________________________________________________________

template <typename TIndex, typename TSpec, typename TText>
inline typename Infix<TText>::Type
infix(Pattern<TIndex, Swift<TSpec> > const & pattern, TText &text)
{
    __int64 hitBegin = pattern.curBeginPos;
    __int64 hitEnd = pattern.curEndPos;
    __int64 textLength = sequenceLength(pattern.curSeqNo, needle(pattern));

    if (hitEnd > textLength) hitEnd = textLength;
    if (hitBegin < 0) hitBegin = 0;

    return infix(text, hitBegin, hitEnd);
}

template <typename TIndex, typename TSpec>
inline typename Infix< typename GetSequenceByNo< TIndex const >::Type >::Type
infix(Pattern<TIndex, Swift<TSpec> > const & pattern)
{
    return infix(pattern, getSequenceByNo(pattern.curSeqNo, needle(pattern)));
}

template <typename TIndex, typename TSpec>
inline typename Infix< typename GetSequenceByNo< TIndex const >::Type >::Type
infix(Pattern<TIndex, Swift<Tag<SwiftSemiGlobal_<TSpec> > > > const & pattern)
{
    return infix(getSequenceByNo(pattern.curSeqNo, needle(pattern)), 0, sequenceLength(pattern.curSeqNo, needle(pattern)));
}

template <typename TIndex, typename TSpec>
inline typename Infix< typename GetSequenceByNo< TIndex const >::Type >::Type
infix(Pattern<TIndex, Swift<TSpec> > & pattern)
{
    return infix(const_cast<Pattern<TIndex, Swift<TSpec> > const &>(pattern));
}

//____________________________________________________________________________

template <typename THaystack, typename TSpec>
inline void
_printDots(Finder<THaystack, Swift<TSpec> > &finder)
{
    while (finder.curPos >= finder.dotPos)
    {
        finder.dotPos += 100000;
        if (finder.dotPos >= finder.dotPos2)
        {
            std::cerr << (finder.dotPos2 / 1000000) << "M" << std::flush;
            finder.dotPos2 += 1000000;
        } else
            std::cerr << "." << std::flush;
    }
}

template <typename TFinder, typename TIndex, typename TSpec>
inline bool
_nextNonRepeatRange(
    TFinder &finder,
    Pattern<TIndex, Swift<TSpec> > &pattern)
{
    //typedef typename TFinder::TRepeat       TRepeat;

    if (finder.curRepeat == finder.endRepeat) return false;

    do
    {
        finder.startPos = (*finder.curRepeat).endPosition;
        if (++finder.curRepeat == finder.endRepeat)
        {
            finder.endPos = length(host(finder));
            if (finder.startPos + length(pattern.shape) > finder.endPos)
                return false;
            else
                break;
        } else
            finder.endPos = (*finder.curRepeat).beginPosition;
        // repeat until the shape fits in non-repeat range
    } while (finder.startPos + length(pattern.shape) > finder.endPos);

    finder.curPos = finder.startPos;
    hostIterator(finder) = begin(host(finder)) + finder.startPos;
    finder.haystackEnd = begin(host(finder)) + (finder.endPos - length(pattern.shape) + 1);

//  if (pattern.params.printDots)
//      std::cerr << std::endl << "  scan range (" << finder.startPos << ", " << finder.endPos << ") " << std::flush;

    return true;
}

template <typename TFinder, typename TIndex, typename TSpec>
inline bool
_firstNonRepeatRange(
    TFinder &finder,
    Pattern<TIndex, Swift<TSpec> > &pattern)
{
    //typedef typename TFinder::TRepeat       TRepeat;

    finder.curRepeat = begin(finder.data_repeats, Standard());
    finder.endRepeat = end(finder.data_repeats, Standard());

    if (finder.curRepeat == finder.endRepeat)
        finder.endPos = length(host(finder));
    else
        finder.endPos = (*finder.curRepeat).beginPosition;

    if (length(pattern.shape) > finder.endPos)
        return _nextNonRepeatRange(finder, pattern);

    finder.curPos = finder.startPos = 0;
    typename Parameter_<typename Host<TFinder>::Type>::Type tmpHost = host(finder);
    hostIterator(finder) = begin(tmpHost);
    finder.haystackEnd = begin(tmpHost) + (finder.endPos - length(pattern.shape) + 1);

//  if (pattern.params.printDots)
//      std::cerr << std::endl << "  scan range (" << finder.startPos << ", " << finder.endPos << ") " << std::flush;

    return true;
}

template <typename TFinder, typename TIndex, typename TSpec>
inline void
_copySwiftHit(
    TFinder &finder,
    Pattern<TIndex, Swift<TSpec> > &pattern)
{
    pattern.curSeqNo = (*finder.curHit).ndlSeqNo;
    pattern.curBeginPos = (*finder.curHit).ndlPos;
    pattern.curEndPos = (*finder.curHit).ndlPos + (*finder.curHit).hitLengthNeedle;
}

template <typename TFinder, typename TIndex, typename TSpec>
inline void
_copySwiftHit(
    TFinder &finder,
    Pattern<TIndex, Swift<Tag<SwiftSemiGlobal_<TSpec> > > > &pattern)
{
    pattern.curSeqNo = (*finder.curHit).ndlSeqNo;
    pattern.curBeginPos = 0;
    pattern.curEndPos = length(indexText(needle(pattern))[pattern.curSeqNo]);
}

template <typename TFinder, typename TIndex, typename TSpec>
inline bool
find(
    TFinder &finder,
    Pattern<TIndex, Swift<Tag<SwiftSemiGlobal_<TSpec> > > > &pattern,
    double errorRate)
{
    return find(finder, pattern, errorRate, 0);
}

template <typename THaystack, typename TIndex, typename TSpec, typename TSize>
inline bool
find(
    Finder<THaystack, Swift<TSpec> > &finder,
    Pattern<TIndex, Swift<TSpec> > &pattern,
    double errorRate,
    TSize minLength)
{
    //typedef typename Fibre<TIndex, QGramShape>::Type    TShape;

    if (empty(finder))
    {
        pattern.finderLength = pattern.params.tabooLength + length(container(finder));
        _patternInit(pattern, errorRate, minLength);
        _finderSetNonEmpty(finder);
        finder.dotPos = 100000;
        finder.dotPos2 = 10 * finder.dotPos;

        if (!_firstNonRepeatRange(finder, pattern)) return false;
        if (_swiftMultiProcessQGram(finder, pattern, hash(pattern.shape, hostIterator(hostIterator(finder)))))
        {
            _copySwiftHit(finder, pattern);
            return true;
        }
    }
    else
    {
        if (++finder.curHit < finder.endHit)
        {
            _copySwiftHit(finder, pattern);
            return true;
        }
    }

    // all previous matches reported -> search new ones
    clear(finder.hits);

    // are we at the end of the text?
    if (atEnd(finder) && finder.curRepeat == finder.endRepeat)
    {
        finder.curHit = finder.endHit;
        return false;
    }

    do
    {
        if (pattern.params.printDots) _printDots(finder);
        if (atEnd(++finder))
        {
            if (!_nextNonRepeatRange(finder, pattern))
            {
                if(_swiftMultiFlushBuckets(finder, pattern))
                {
                    _copySwiftHit(finder, pattern);
                    return true;
                }
                else
                    return false;
            }
            hash(pattern.shape, hostIterator(hostIterator(finder)));
        }
        else
        {
            ++finder.curPos;
            hashNext(pattern.shape, hostIterator(hostIterator(finder)));
        }

        if (_swiftMultiProcessQGram(finder, pattern, value(pattern.shape)))
        {
            _copySwiftHit(finder, pattern);
            return true;
        }

    } while (true);
}

template <typename THashes, typename TPipeSpec, typename TIndex, typename TSpec>
inline bool
find(
    Finder<Pipe<THashes, TPipeSpec>, Swift<TSpec> > &finder,
    Pattern<TIndex, Swift<TSpec> > &pattern,
    double errorRate)
{
    if (empty(finder))
    {
        pattern.finderLength = 0;
        _patternInit(pattern, errorRate, 0);
        _finderSetNonEmpty(finder);
        finder.dotPos = 100000;
        finder.dotPos2 = 10 * finder.dotPos;

        beginRead(finder.in);
        if (eof(finder.in))
        {
            endRead(finder.in);
            return false;
        }
        finder.curPos = (*finder.in).i1;
        if (_swiftMultiProcessQGram(finder, pattern, hash(pattern.shape, (*finder.in).i2)))
        {
            _copySwiftHit(finder, pattern);
            return true;
        }
    } else
        if (++finder.curHit != finder.endHit)
        {
            _copySwiftHit(finder, pattern);
            return true;
        }

    clear(finder.hits);
    if (eof(finder.in)) return false;

    do
    {
        ++finder.in;
        if (eof(finder.in))
        {
            endRead(finder.in);
#ifdef SEQAN_DEBUG_SWIFT
            _printSwiftBuckets(pattern);
#endif
            if(_swiftMultiFlushBuckets(finder, pattern))
            {
                _copySwiftHit(finder, pattern);
                return true;
            }
            else
                return false;
        }
        finder.curPos = (*finder.in).i1;
        if (pattern.params.printDots) _printDots(finder);

    } while (!_swiftMultiProcessQGram(finder, pattern, hash(pattern.shape, (*finder.in).i2)));

    _copySwiftHit(finder, pattern);
    return true;
}

/*!
 * @fn SwiftFinder#windowFindBegin
 * @headerfile <seqan/index.h>
 * @brief Initializes the pattern. Sets the finder on the begin position.  Gets the first non-repeat range and sets it
 *        in the finder.  Used together with @link SwiftFinder#windowFindEnd @endlink.
 *
 * @signature bool windowFindBegin(finder, pattern, errorRate);
 *
 * @param[in,out] finder    A finder with window interface.
 * @param[in,out] pattern   A pattern with window interface. Types: @link SwiftPattern @endlink.
 * @param[in]     errorRate Error rate that is allowed between reads and reference.  Should be the same in as in
 *                          @link SwiftFinder#windowFindNext @endlink.  Types: <tt>double</tt>
 *
 * @return bool <tt>true</tt> if there were bases to be scanned and <tt>false</tt> otherwise.
 *
 * @see SwiftFinder#windowFindNext
 * @see SwiftFinder#windowFindEnd
 */
template <typename THaystack, typename TIndex, typename TSpec>
inline bool
windowFindBegin(
    Finder<THaystack, Swift<TSpec> > &finder,
    Pattern<TIndex, Swift<TSpec> > &pattern,
    double errorRate)
{
    SEQAN_CHECKPOINT

    pattern.finderLength = pattern.params.tabooLength + length(container(finder));
    _patternInit(pattern, errorRate, 0);
    _finderSetNonEmpty(finder);
    finder.dotPos = 100000;
    finder.dotPos2 = 10 * finder.dotPos;
    finder.windowStart = 0;

    if (!_firstNonRepeatRange(finder, pattern)) return false;

    return true;
}


/*!
 * @fn SwiftFinder#windowFindNext
 * @headerfile <seqan/index.h>
 * @brief Searches over the next window with the finder. The found hits can be retrieved with @link
 *        SwiftFinder#getWindowFindHits @endlink.  Used together with @link SwiftFinder#windowFindBegin @endlink and
 *        @link SwiftFinder#windowFindEnd @endlink.
 *
 * @signature bool windowFindNext(finder, pattern, finderWindowLength)
 *
 * @param[in,out] finder             A SwiftFinder.
 * @param[in,out] pattern            A SwiftPattern.
 * @param[in] finderWindowLength     Number of bases that are scanned beginning from the position the finder is at.
 *                                   Including bases that are marked as repeats and that are skipped. Types:
 *                                   nolink:unsigned int
 *
 * @return bool <tt>true</tt>, if there are bases that can be scanned, <tt>false</tt> otherwise.
 *
 * @see SwiftFinder#windowFindBegin
 * @see SwiftFinder#windowFindEnd
 * @see SwiftFinder#getWindowFindHits
 */

template <typename THaystack, typename TIndex, typename TSpec, typename TSize>
inline bool
windowFindNext(
    Finder<THaystack, Swift<TSpec> > &finder,
    Pattern<TIndex, Swift<TSpec> > &pattern,
    TSize finderWindowLength
    )
{
    SEQAN_CHECKPOINT

    typedef typename Fibre<TIndex, QGramShape>::Type    TShape;

    typedef Finder<THaystack, Swift<TSpec> >            TFinder;
    typedef typename TFinder::THstkPos                  THstkPos;

    // all previous matches reported -> search new ones
    clear(finder.hits);

    THstkPos windowEnd = finder.windowStart + finderWindowLength;

    // iterate over all non-repeat regions within the window
    for (; finder.curPos < windowEnd; )
    {
        THstkPos nonRepeatEnd = finder.endPos - length(pattern.shape) + 1;
        THstkPos localEnd = _min(windowEnd, nonRepeatEnd);

        // filter a non-repeat region within the window
        if (finder.curPos < localEnd)
        {
            TShape &shape = pattern.shape;
            _swiftMultiProcessQGram(finder, pattern, hash(shape, hostIterator(hostIterator(finder))));

            for (++finder.curPos, ++finder; finder.curPos < localEnd; ++finder.curPos, ++finder){
                _swiftMultiProcessQGram(finder, pattern, hashNext(shape, hostIterator(hostIterator(finder))));
            }
        }

        if (pattern.params.printDots) _printDots(finder);

        if (finder.curPos >= nonRepeatEnd)
            if (!_nextNonRepeatRange(finder, pattern))
            {
                finder.windowStart = windowEnd;
                return false;
            }
    }
    finder.windowStart = windowEnd;
    return true;
}

/*!
 * @fn SwiftFinder#windowFindEnd
 * @headerfile <seqan/index.h>
 * @brief Flushes the finder.  Used together with @link SwiftFinder#windowFindBegin @endlink and @link
 *        SwiftFinder#windowFindNext @endlink.
 *
 * @signature void windowFindEnd(finder, pattern);
 *
 * @param[in,out] finder  A finder with window interface.
 * @param[in,out] pattern A pattern with window interface. Types: @link SwiftPattern @endlink
 *
 * @see SwiftFinder#windowFindBegin
 * @see SwiftFinder#windowFindNext
 */

template <typename THaystack, typename TIndex, typename TSpec>
inline void
windowFindEnd(
    Finder<THaystack, Swift<TSpec> > & finder,
    Pattern<TIndex, Swift<TSpec> > &pattern)
{
    SEQAN_CHECKPOINT
    _swiftMultiFlushBuckets(finder, pattern);
}

/*!
 * @fn SwiftFinder#getWindowFindHits
 * @headerfile <seqan/index.h>
 * @brief Returns the string of hits from the finder.
 *
 * @signature THitString getWindowFindHits(finder);

 * @param[in] finder     A finder with window interface.
 * @return    THitString @link String @endlink of hits, see @link SwiftFinder::THitString @endlink.
 *
 * @see SwiftFinder#windowFindNext
 */

template <typename THaystack, typename TSpec>
inline typename WindowFindResult<Finder<THaystack, Swift<TSpec> >, void>::Type &
getWindowFindHits(Finder<THaystack, Swift<TSpec> > &finder)
{
    SEQAN_CHECKPOINT

    return finder.hits;
}

/*!
 * @fn SwiftPattern#getMaxDeviationOfOrder
 * @headerfile <seqan/index.h>
 * @brief Returns the maximal out-of-order distance of adjacent hits.
 *
 * @signature TSize getMaxDeviationOfOrder(pattern);
 *
 * @param[in] pattern A SwiftPattern.
 * @return    TSize   Returns the maximal distance two adjacent hits can have which are not in increasing order.
 *                    <tt>TSize</tt> is the size type of the underlying index.
 */

template <typename TIndex, typename TSpec>
inline typename Size<TIndex>::Type
getMaxDeviationOfOrder(Pattern<TIndex, Swift<TSpec> > &pattern)
{
    // TODO(holtgrew): Can pattern be made const?
    return back(pattern.bucketParams).delta + back(pattern.bucketParams).overlap + length(pattern.bucketParams) - 2;
}


}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_FIND_SHIFTAND_H

