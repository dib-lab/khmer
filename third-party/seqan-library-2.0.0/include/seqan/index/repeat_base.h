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

#ifndef SEQAN_HEADER_REPEAT_BASE_H
#define SEQAN_HEADER_REPEAT_BASE_H

#if SEQAN_ENABLE_PARALLELISM
#include <seqan/parallel.h>
#endif  // #if SEQAN_ENABLE_PARALLELISM

namespace seqan {

/*!
 * @class Repeat
 * @headerfile <seqan/index.h>
 * @brief Store information about a repeat.
 *
 * @signature template <typename TPos, typename TPeriod>
 *            struct Repeat;
 *
 * @tparam TPeriod Type to use for storing the repeat period. Default: 1
 * @tparam TPos Type to use for storing positions.
 *
 * @see findRepeats
 *
 * @var TPos Repeat::endPosition;
 * @brief The end position of the repeat of type <tt>TPos</tt>.
 *
 * @var TPos Repeat::beginPosition;
 * @brief The begin position of the repeat of type <tt>TPos</tt>.
 *
 * @var TPeriod Repeat::period;
 * @brief The period of the repeat of type <tt>TPeriod</tt>.
 */

    template <typename TPos, typename TPeriod>
    struct Repeat {
        TPos        beginPosition;
        TPos        endPosition;
        TPeriod        period;
    };

    template <typename TPos, typename TPeriod>
    struct Value< Repeat<TPos, TPeriod> > {
        typedef TPos Type;
    };

    template <typename TPos, typename TPeriod>
    struct Size< Repeat<TPos, TPeriod> > {
        typedef TPeriod Type;
    };


    template <typename TSize>
    struct RepeatFinderParams {
        TSize minRepeatLen;
        TSize maxPeriod;
    };

    // custom TSpec for our customized wotd-Index
    struct TRepeatFinder;

    template <typename TText>
    struct Cargo<Index<TText, IndexWotd<TRepeatFinder> > >
    {
        typedef Index<TText, IndexWotd<TRepeatFinder> >    TIndex;
        typedef typename Size<TIndex>::Type                    TSize;
        typedef RepeatFinderParams<TSize>                    Type;
    };


    // node predicate
    template <typename TText, typename TSpec>
    bool nodePredicate(Iter<Index<TText, IndexWotd<TRepeatFinder> >, TSpec> &it)
    {
//        return countOccurrences(it) * nodeDepth(it) >= cargo(container(it)).minRepeatLen;
        return countOccurrences(it) * repLength(it) >= cargo(container(it)).minRepeatLen;
    }

    // monotonic hull
    template <typename TText, typename TSpec>
    bool nodeHullPredicate(Iter<Index<TText, IndexWotd<TRepeatFinder> >, TSpec> &it)
    {
//        return nodeDepth(it) <= cargo(container(it)).maxPeriod;
        return repLength(it) <= cargo(container(it)).maxPeriod;
    }

    template <typename TPos>
    struct RepeatLess_ : public std::binary_function<TPos, TPos, bool>
    {
        // key less
        inline bool operator() (TPos const &a, TPos const &b) const {
            return posLess(a, b);
        }
    };

    template <typename TValue>
    inline bool _repeatMaskValue(TValue const &)
    {
        // TODO(holtgrew): Maybe use unknownValue<TValue>() instead of specializing for all alphabets, especially since we have Rna5 now and might want Rna5Q later.
        return false;
    }

    template <>
    inline bool _repeatMaskValue(Dna5 const &val)
    {
        return val == unknownValue<Dna5>(); // 'N'
    }

    template <>
    inline bool _repeatMaskValue(Dna5Q const &val)
    {
        return val == unknownValue<Dna5Q>(); // 'N'
    }

    template <>
    inline bool _repeatMaskValue(Iupac const &val)
    {
        return val == unknownValue<Iupac>(); // 'N'
    }
/*
    template <>
    inline bool _repeatMaskValue(AminoAcid val)
    {
        return val == 'X';
    }
*/
/*!
 * @fn findRepeats
 * @headerfile <seqan/index.h>
 * @brief Search for repeats in a text.
 *
 * @signature void findRepeats(repeatString, text, minRepeatLength[, maxPeriod]);
 *
 * @param[out] repeatString    A @link String @endlink of @link Repeat @endlink objects.
 * @param[in]  text            The text to search repeats in.  Types: @link ContainerConcept @endlink
 * @param[in]  minRepeatLength The minimum length each reported repeat must have.
 * @param[in]  maxPeriod       Optionally, the maximal period that reported repeats can have. Default: 1
 *
 * Subsequences of undefined values/<tt>N</tt>s will always be reported.
 *
 * @section Examples
 *
 * The following demonstrates finding repeats of period 3.
 *
 * @include demos/index/find_repeats.cpp
 *
 * @code{.console}
 * # of repeats: 15
 * i == 0, beginPosition = 3, endPosition = 7, period = 1
 * i == 1, beginPosition = 46, endPosition = 53, period = 1
 * i == 2, beginPosition = 101, endPosition = 105, period = 1
 * i == 3, beginPosition = 105, endPosition = 109, period = 1
 * i == 4, beginPosition = 164, endPosition = 169, period = 1
 * i == 5, beginPosition = 291, endPosition = 297, period = 1
 * i == 6, beginPosition = 319, endPosition = 327, period = 1
 * i == 7, beginPosition = 400, endPosition = 404, period = 1
 * i == 8, beginPosition = 442, endPosition = 446, period = 1
 * i == 9, beginPosition = 468, endPosition = 473, period = 1
 * i == 10, beginPosition = 476, endPosition = 480, period = 1
 * i == 11, beginPosition = 507, endPosition = 513, period = 1
 * i == 12, beginPosition = 561, endPosition = 566, period = 1
 * i == 13, beginPosition = 623, endPosition = 627, period = 1
 * i == 14, beginPosition = 655, endPosition = 659, period = 1
 * @endcode
 *
 * @see AlphabetWithUnknownValueConcept#unknownValue
 * @see Repeat
 */
// TODO(holtgrew): minRepeatLength is 1-off.

    // period-1 optimization
    template <typename TRepeatStore, typename TString, typename TRepeatSize>
    inline void findRepeats(TRepeatStore &repString, TString const &text, TRepeatSize minRepeatLen)
    {
        typedef typename Value<TRepeatStore>::Type    TRepeat;
        typedef typename Iterator<TString const>::Type    TIterator;
        typedef typename Size<TString>::Type        TSize;

#if SEQAN_ENABLE_PARALLELISM
        typedef typename Value<TString>::Type        TValue;

        if (length(text) > (TSize)(omp_get_max_threads() * 2 * minRepeatLen)) {
            // std::cerr << ">>> PARALLEL WABOOGIE!" << std::endl;
            // std::cerr << "omp_get_max_threads() == " << omp_get_max_threads() << std::endl;
            // Parallel case.

            // NOTE(holtgrew): The minimum text length check above makes it impossible that more than two chunks are
            // required to form an otherwise too short repeat.

            // TODO(holtgrew): Load balancing? Probably not worth it.
            String<TSize> splitters;
            String<TRepeatStore> threadLocalStores;

            // Each threads finds repeats on its chunk in parallel.
            #pragma omp parallel
            {
                // We have to determine the number of available threads at this point.  We will use the number of thread
                // local stores to determin the number of available threads later on.
                #pragma omp master
                {
                    // std::cerr << "omp_get_num_threads() == " << omp_get_num_threads() << std::endl;
                    computeSplitters(splitters, length(text), omp_get_num_threads());
                    resize(threadLocalStores, omp_get_num_threads());
                }  // end of #pragma omp master
                #pragma omp barrier

                int const t = omp_get_thread_num();
                TRepeatStore & store = threadLocalStores[t];

                TRepeat rep;
                rep.beginPosition = 0;
                rep.endPosition = 0;
                rep.period = 1;

                // Flags used for force-adding repeats for the chunks that have a left/right neighbour.
                bool forceFirst = t > 0;
                bool forceLast = (t + 1) < omp_get_num_threads();
                // #pragma omp critical
                // std::cerr << "omp_get_num_threads() == " << omp_get_num_threads() << std::endl;

                TIterator it = iter(text, splitters[t], Standard());
                TIterator itEnd = iter(text, splitters[t + 1], Standard());
                if (it != itEnd)
                {
                    TValue last = *it;
                    TSize repLeft = 0;
                    TSize repRight = 1;

                    for (++it; it != itEnd; ++it, ++repRight)
                    {
                        if (*it != last)
                        {
                            // #pragma omp critical
                            // std::cerr << "t == " << t << ", last == " << last << ", repRight = " << repRight << ", repLeft == " << repLeft << ", minRepeatLen = " << minRepeatLen << ", forceFirst = " << forceFirst << std::endl;
                            if (_repeatMaskValue(last) || (TRepeatSize)(repRight - repLeft) > minRepeatLen || forceFirst)
                            {
                                forceFirst = false;
                                // insert repeat
                                rep.beginPosition = splitters[t] + repLeft;
                                rep.endPosition = splitters[t] + repRight;
                                // #pragma omp critical
                                // std::cerr << "  t == " << t << ", append" << std::endl;
                                appendValue(store, rep);
                            }
                            repLeft = repRight;
                            last = *it;
                        }
                    }
                    // #pragma omp critical
                    // std::cerr << "t == " << t << ", last == " << last << ", repRight = " << repRight << ", repLeft == " << repLeft << ", minRepeatLen = " << minRepeatLen << ", forceLast = " << forceLast << std::endl;
                    if (_repeatMaskValue(last) || (TRepeatSize)(repRight - repLeft) > minRepeatLen || forceLast)
                    {
                        // Insert repeat but only if it is not already in there.
                        if (empty(store) || (back(store).beginPosition != repLeft && back(store).endPosition != repRight))
                        {
                            rep.beginPosition = splitters[t] + repLeft;
                            rep.endPosition = splitters[t] + repRight;
                            // #pragma omp critical
                            // std::cerr << "  t == " << t << ", append" << std::endl;
                            appendValue(store, rep);
                        }
                    }
                }
            }  // end of #pragma omp parallel

            // std::cerr << ",-- REPEATS BEFORE MENDING\n";
            // for (unsigned i = 0; i < length(threadLocalStores); ++i)
            // {
            //     std::cerr << "| i = " << i << std::endl;
            //     for (unsigned j = 0; j < length(threadLocalStores[i]); ++j)
            //         std::cerr << "| threadLocalStores[" << i << "][" << j << "] == {" << threadLocalStores[i][j].beginPosition << ", " << threadLocalStores[i][j].endPosition << "}" << std::endl;
            // }
            // std::cerr << "`--" << std::endl;

            // Mend the splice points.
            //
            // We will copy out infixes described by fromPositions.
            String<Pair<TSize> > fromPositions;
            resize(fromPositions, length(threadLocalStores));
            for (unsigned i = 0; i < length(fromPositions); ++i)
            {
                fromPositions[i].i1 = 0;
                fromPositions[i].i2 = length(threadLocalStores[i]);
            }
            // First, merge repeats spanning blocks.  Do this iteratively until all has been merged.
            bool anyChange;
            do
            {
                anyChange = false;
                int lastNonEmpty = -1;
                for (unsigned i = 0; i < length(threadLocalStores); ++i)
                {
                    if (fromPositions[i].i1 == fromPositions[i].i2)
                        continue;  // Skip empty buckets.

                    if (lastNonEmpty != -1)
                    {
                        bool const adjacent = back(threadLocalStores[lastNonEmpty]).endPosition == front(threadLocalStores[i]).beginPosition;
                        bool const charsEqual = text[back(threadLocalStores[lastNonEmpty]).beginPosition] == text[front(threadLocalStores[i]).beginPosition];
                        if (adjacent && charsEqual)
                        {
                            anyChange = true;
                            back(threadLocalStores[lastNonEmpty]).endPosition = front(threadLocalStores[i]).endPosition;
                            fromPositions[i].i1 += 1;
                        }
                    }

                    if (fromPositions[i].i1 != fromPositions[i].i2)
                        lastNonEmpty = i;
                }
            }
            while (anyChange);
            // Then, remove any repeats in the beginning and end of blocks that are too short.
            for (unsigned i = 0; i < length(threadLocalStores); ++i)
            {
                if (fromPositions[i].i1 == fromPositions[i].i2)
                    continue;
                unsigned j = fromPositions[i].i1;
                TRepeatSize len = threadLocalStores[i][j].endPosition - threadLocalStores[i][j].beginPosition;
                if (!_repeatMaskValue(text[threadLocalStores[i][j].beginPosition]) &&  // Never remove mask value.
                    len <= minRepeatLen)
                    fromPositions[i].i1 += 1;
                if (fromPositions[i].i1 == fromPositions[i].i2)
                    continue;
                j = fromPositions[i].i2 - 1;
                len = threadLocalStores[i][j].endPosition - threadLocalStores[i][j].beginPosition;
                if (!_repeatMaskValue(text[threadLocalStores[i][j].beginPosition]) &&  // Never remove mask value.
                    len <= minRepeatLen)
                    fromPositions[i].i2 -= 1;
            }
            // Last, build splitters for output in parallel.
            String<unsigned> outSplitters;
            appendValue(outSplitters, 0);
            for (unsigned i = 0; i < length(threadLocalStores); ++i)
                appendValue(outSplitters, back(outSplitters) + fromPositions[i].i2 - fromPositions[i].i1);

            // std::cerr << ",-- REPEATS AFTER MENDING\n";
            // for (unsigned i = 0; i < length(threadLocalStores); ++i)
            // {
            //     std::cerr << "| i = " << i << std::endl;
            //     std::cerr << "`--, fromPositions[" << i << "] = (" << fromPositions[i].i1 << ", " << fromPositions[i].i2 << std::endl;
            //     for (unsigned j = 0; j < length(threadLocalStores[i]); ++j)
            //         std::cerr << "   | threadLocalStores[" << i << "][" << j << "] == {" << threadLocalStores[i][j].beginPosition << ", " << threadLocalStores[i][j].endPosition << "}" << std::endl;
            // }
            // std::cerr << "    `--" << std::endl;

            // Allocate memory.
            clear(repString);
            resize(repString, back(outSplitters));

            // Copy back the repeats in parallel.
            unsigned nt = length(threadLocalStores);
            (void) nt;  // Otherwise, GCC 4.6 warns, does not see it used in pragma clause below.
            #pragma omp parallel num_threads(nt)
            {
                int const t = omp_get_thread_num();
                arrayCopy(iter(threadLocalStores[t], fromPositions[t].i1, Standard()),
                          iter(threadLocalStores[t], fromPositions[t].i2, Standard()),
                          iter(repString, outSplitters[t], Standard()));
            }  // end of #pragma omp parallel
        } else {
#endif  // #if SEQAN_ENABLE_PARALLELISM
            // Sequential case.
            TRepeat rep;
            rep.period = 1;
            clear(repString);

            TIterator it = begin(text, Standard());
            TIterator itEnd = end(text, Standard());
            if (it == itEnd) return;

            TSize repLen = 1;
            for (++it; it != itEnd; ++it)
            {
                if (*it != *(it-1))
                {
                    if (_repeatMaskValue(*(it-1)) || repLen > (TSize)minRepeatLen)
                    {
                        // insert repeat
                        rep.endPosition = it - begin(text, Standard());
                        rep.beginPosition = rep.endPosition - repLen;
                        //                    std::cerr<<"left:"<<rep.beginPosition<<"  right:"<<rep.endPosition<<"  length:"<<posSub(rep.endPosition,rep.beginPosition)<<"  period:"<<rep.period<<std::endl;
                        appendValue(repString, rep);
                    }
                    repLen = 1;
                } else
                    ++repLen;
            }
            if (_repeatMaskValue(*(it-1)) || repLen > (TSize)minRepeatLen)
            {
                // insert repeat
                rep.endPosition = length(text);
                rep.beginPosition = rep.endPosition - repLen;
                //            std::cerr<<"left:"<<rep.beginPosition<<"  right:"<<rep.endPosition<<"  length:"<<posSub(rep.endPosition,rep.beginPosition)<<"  period:"<<rep.period<<std::endl;
                appendValue(repString, rep);
            }
#if SEQAN_ENABLE_PARALLELISM
        }
#endif  // #if SEQAN_ENABLE_PARALLELISM
        // #pragma omp critical
        // {
        //     std::cerr << "thread #" << omp_get_thread_num() << " REPEATS:";
        //     for (unsigned i = 0; i < length(repString); ++i) {
        //         std::cerr << " (" << repString[i].beginPosition << ", " << repString[i].endPosition << ", " << repString[i].period << ")";
        //     }
        //     std::cerr << std::endl;
        // }
    }

    // TODO(holtgrew): Why for TString const and StringSet<> const?
    template <typename TRepeatStore, typename TString, typename TSpec, typename TRepeatSize>
    inline void findRepeats(TRepeatStore &repString, StringSet<TString, TSpec> const &text, TRepeatSize minRepeatLen)
    {
        typedef typename Value<TRepeatStore>::Type    TRepeat;
        typedef typename Iterator<TString>::Type    TIterator;
        typedef typename Value<TString>::Type        TValue;
        typedef typename Size<TString>::Type        TSize;

        TRepeat rep;
        rep.period = 1;
        clear(repString);

        for (unsigned i = 0; i < length(text); ++i)
        {
            TIterator it = begin(text[i], Standard());
            TIterator itEnd = end(text[i], Standard());
            if (it == itEnd) continue;

            TValue last = *it;
            TSize repLeft = 0;
            TSize repRight = 1;
            rep.beginPosition.i1 = i;
            rep.endPosition.i1 = i;

            for (++it; it != itEnd; ++it, ++repRight)
            {
                if (last != *it)
                {
                    if (_repeatMaskValue(last) || (TRepeatSize)(repRight - repLeft) > minRepeatLen)
                    {
                        // insert repeat
                        rep.beginPosition.i2 = repLeft;
                        rep.endPosition.i2 = repRight;
//                        std::cerr<<"left:"<<rep.beginPosition<<"  right:"<<rep.endPosition<<"  length:"<<posSub(rep.endPosition,rep.beginPosition)<<"  period:"<<rep.period<<std::endl;
                        appendValue(repString, rep);
                    }
                    repLeft = repRight;
                    last = *it;
                }
            }
            if (_repeatMaskValue(last) || (TRepeatSize)(repRight - repLeft) > minRepeatLen)
            {
                // insert repeat
                rep.beginPosition.i2 = repLeft;
                rep.endPosition.i2 = repRight;
//                std::cerr<<"left:"<<rep.beginPosition<<"  right:"<<rep.endPosition<<"  length:"<<posSub(rep.endPosition,rep.beginPosition)<<"  period:"<<rep.period<<std::endl;
                appendValue(repString, rep);
            }
        }
    }

    // main function
    template <typename TRepeatStore, typename TText, typename TRepeatSize, typename TPeriodSize>
    void findRepeats(TRepeatStore &repString, TText const &text, TRepeatSize minRepeatLen, TPeriodSize maxPeriod)
    {
        typedef Index<TText, IndexWotd<TRepeatFinder> >                    TIndex;
        typedef typename Size<TIndex>::Type                                    TSize;
        typedef typename Iterator<TIndex, TopDown<ParentLinks<> > >::Type    TNodeIterator;
        typedef typename Fibre<TIndex, FibreSA>::Type const                TSA;
        typedef typename Infix<TSA>::Type                                    TOccString;
        typedef typename Iterator<TOccString>::Type                            TOccIterator;

        typedef typename Value<TRepeatStore>::Type                            TRepeat;
        typedef typename Value<TOccString>::Type                            TOcc;

        typedef std::map<TOcc,TRepeat,RepeatLess_<TOcc> >                    TRepeatList;

        if (maxPeriod < 1) return;
        if (maxPeriod == 1)
        {
            findRepeats(repString, text, minRepeatLen);
            return;
        }

        TIndex        index(text);
        TRepeatList list;

        // set repeat finder parameters
        cargo(index).minRepeatLen = minRepeatLen;
        cargo(index).maxPeriod = maxPeriod;

        TNodeIterator nodeIt(index);
        TOccIterator itA, itB, itRepBegin, itEnd;
        TRepeat rep;
        for (; !atEnd(nodeIt); goNext(nodeIt))
        {
            if (isRoot(nodeIt)) continue;

            // get occurrences
            TOccString occ = getOccurrences(nodeIt);
            itA = begin(occ, Standard());
            itEnd = end(occ, Standard());
            itRepBegin = itB = itA;

            TSize repLen = repLength(nodeIt);        // representative length
            if ((TSize)minRepeatLen <= repLen) continue;

            TSize diff, period = 0;                    // period of current repeat
            TSize repeatLen = 0;                    // overall length of current repeat
            TSize minLen = minRepeatLen - repLen;    // minimum repeat length minus length of representative

            for (++itB; itB != itEnd; ++itB)
            {
                diff = posSub(*itB, *itA);
                if (diff != period || getSeqNo(*itA) != getSeqNo(*itB))
                {
                    // is the repeat long enough?
                    if (repeatLen >= minLen)
                        // is the repeat self overlapping or connected?
                        if (parentRepLength(nodeIt) < period && period <= repLen)
                        {
                            // insert repeat
                            rep.beginPosition = *itRepBegin;
                            rep.endPosition = posAdd(*itA, period);
                            rep.period = period;
//                            std::cerr<<"left:"<<rep.beginPosition<<"  right:"<<rep.endPosition<<"  length:"<<posSub(rep.endPosition,rep.beginPosition)<<"  period:"<<rep.period<<std::endl;
                            list.insert(std::pair<TOcc,TRepeat>(rep.beginPosition, rep));
                        }
                    itRepBegin = itA;
                    period = diff;
                    repeatLen = 0;
                }
                repeatLen += period;
                itA = itB;
            }

            // is the last repeat long enough?
            if (repeatLen >= minLen)
                // is the repeat self overlapping or connected?
                if (parentRepLength(nodeIt) < period && period <= repLen)
                {
                    // insert repeat
                    rep.beginPosition = *itRepBegin;
                    rep.endPosition = posAdd(*itA, period);
                    rep.period = period;
//                    std::cerr<<"left:"<<rep.beginPosition<<"  right:"<<rep.endPosition<<"  length:"<<posSub(rep.endPosition,rep.beginPosition)<<"  period:"<<rep.period<<std::endl;
                    list.insert(std::pair<TOcc,TRepeat>(rep.beginPosition, rep));
                }
        }

        // copy low-complex regions to result string
        clear(repString);
        reserve(repString, list.size(), Exact());
        typename TRepeatList::const_iterator lit = list.begin();
        typename TRepeatList::const_iterator litEnd = list.end();
        for (TSize i = 0; lit != litEnd; ++lit, ++i)
            appendValue(repString, (*lit).second);
    }

}    // namespace seqan

#endif
