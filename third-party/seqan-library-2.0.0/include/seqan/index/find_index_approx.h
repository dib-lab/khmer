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

#ifndef SEQAN_HEADER_FIND_INDEX_APPROX_H
#define SEQAN_HEADER_FIND_INDEX_APPROX_H

namespace SEQAN_NAMESPACE_MAIN
{

template <typename TSAValue1, typename TSAValue2>
struct AddResultsFunctor
{
    String<Pair<TSAValue1, String<TSAValue2> > > results;

    template <typename TIter1, typename TIter2>
    void operator() (TIter1 &iter1, TIter2 &iter2)
    {
        unsigned ofs = length(results);
        resize(results, ofs + countOccurrences(iter1));
        for (unsigned i = 0; i < countOccurrences(iter1); ++i)
        {
            results[ofs + i].i1 = getOccurrences(iter1)[i];
            results[ofs + i].i2 = getOccurrences(iter2);
        }
    }
};

template <
    bool enumerateA,
    bool enumerateB,
    typename TOnFoundFunctor,
    typename TTreeIteratorA,
    typename TIterPosA,
    typename TTreeIteratorB,
    typename TIterPosB,
    typename TErrors >
inline void
_exactTreeSearch(
    TOnFoundFunctor    &onFoundFunctor,    // functor to store matches
    TTreeIteratorA    iterA,              // pattern tree iterator
    TIterPosA        iterPosA,           // position in representative
    TTreeIteratorB    iterB,              // text suffix tree iterator
    TIterPosB        iterPosB)           // position in representative
{
    while (true)
    {
        TIterPosA repLenA = repLength(iterA);
        TIterPosB repLenB = repLength(iterB);
        if (iterPosA < repLenA)
        {
            // search patterns in the suffix tree
            if (!goDown(iterB, suffix(representative(iterA), iterPosA)))
                return;
            // end of pattern found?
            if (isLeaf(iterA))
            {
                onFoundFunctor(iterA, iterB);
                return;
            }
            iterPosB += repLenA - iterPosA;
            iterPosA = repLenA;
        }
        else if (iterPosB < repLenB)
        {
            // search text in the pattern tree
            TIterPosA lcp = 0;
            if (!_goDownString(iterA, suffix(representative(iterB), iterPosB), lcp))
            {
                // couldn't we proceed to go down because of a leaf in the pattern tree?
                if (iterPosA + lcp == repLength(iterA) && isLeaf(iterA))
                    onFoundFunctor(iterA, iterB);
                return;
            }
            // end of pattern found?
            if (isLeaf(iterA))
            {
                onFoundFunctor(iterA, iterB);
                return;
            }
            iterPosA += repLenB - iterPosB;
            iterPosB = repLenB;
        }
        else
        {
            if (!goDown(iterA))
            {
                onFoundFunctor(iterA, iterB);
                return;
            }
            while (emptyParentEdge(iterA))
            {
                onFoundFunctor(iterA, iterB);
                if (!goRight(iterA))
                    return;
            }

            if (!goDown(iterB))
                return;
            while (emptyParentEdge(iterB))
            {
                if (!goRight(iterA))
                    return;
            }

            // search pairs of edges with the same character
            while (true)
            {
                if (parentEdgeFirstChar(iterA) < parentEdgeFirstChar(iterB))
                {
                    if (!goRight(iterA))
                        return;
                }
                else if (parentEdgeFirstChar(iterB) < parentEdgeFirstChar(iterA))
                {
                    if (!goRight(iterB))
                        return;
                }
                else
                {
                    _exactTreeSearch(onFoundFunctor, iterA, iterPosA, iterB, iterPosB);
                    if (!goRight(iterA) || !goRight(iterB))
                        return;
                }
            }
        }
    }
}

template <
    bool enumerateA,
    bool enumerateB,
    typename TOnFoundFunctor,
    typename TTreeIteratorA,
    typename TIterPosA,
    typename TTreeIteratorB,
    typename TIterPosB,
    typename TErrors >
inline void
_approximateTreeSearch(
    TOnFoundFunctor    &onFoundFunctor,    // functor to store matches
    TTreeIteratorA    iterA,              // pattern tree iterator
    TIterPosA        iterPosA,           // position in representative
    TTreeIteratorB    iterB_,             // text suffix tree iterator
    TIterPosB        iterPosB,           // position in representative
    TErrors            errorsLeft)         // tolerated errors left
{
    if (enumerateA && !goDown(iterA))
    {
        onFoundFunctor(iterA, iterB_);
        return;
    }
    if (enumerateB && !goDown(iterB_)) return;

    do
    {
        TTreeIteratorB iterB = iterB_;
        do
        {
            TErrors e = errorsLeft;
            TIterPosA ipA = iterPosA;
            TIterPosB ipB = iterPosB;

            if (ipB == repLength(iterB)) continue;

            while (true)
            {
                if (ipA == repLength(iterA))
                {
                    if (ipB == repLength(iterB))
                        _approximateTreeSearch<true,true>(onFoundFunctor, iterA, ipA, iterB, ipB, e);
                    else
                        _approximateTreeSearch<true,false>(onFoundFunctor, iterA, ipA, iterB, ipB, e);
                    break;
                }
                else if (ipB == repLength(iterB))
                {
                    _approximateTreeSearch<false,true>(onFoundFunctor, iterA, ipA, iterB, ipB, e);
                    break;
                }

                if (representative(iterA)[ipA] != representative(iterB)[ipB])
                    if (e-- == 0) break;

                ++ipA;
                ++ipB;
            }
        } while (enumerateB && goRight(iterB));
    } while (enumerateA && goRight(iterA));
}


template <typename TSAValue>
struct AddSingleResultsFunctor
{
    String<TSAValue> results;

    template <typename TPattern, typename TIter>
    void operator() (TPattern & /*pattern*/, TIter &iter)
    {
        append(results, getOccurrences(iter));
    }
};


template <
    typename TOnFoundFunctor,
    typename TString,
    typename TStringPos,
    typename TTreeIterator,
    typename TIterPos,
    typename TErrors >
inline void
_approximateStringSearch(
    TOnFoundFunctor    &onFoundFunctor,    // functor to store matches
    TString const &string,              // search pattern
    TStringPos stringPos,               // position in pattern
    TTreeIterator iter,                 // text suffix tree iterator
    TIterPos iterPos,                   // position in representative
    TErrors errorsLeft)                 // tolerated errors left
{
    if (errorsLeft == 0)
    {
        if (goDown(iter, suffix(string, stringPos)))
            onFoundFunctor(string, iter);
        return;
    }

    if (!goDown(iter)) return;
    do
    {
        TErrors e = errorsLeft;
        TStringPos sp = stringPos;
        TIterPos ip = iterPos;

        if (ip == repLength(iter)) continue;

        while (true)
        {
            if (representative(iter)[ip] != string[sp])
                if (e-- == 0) break;

            if (++sp == length(string))
            {
                onFoundFunctor(string, iter);
                break;
            }

            if (++ip == repLength(iter))
            {
                _approximateStringSearch(onFoundFunctor, string, sp, iter, ip, e);
                break;
            }
        }
    } while (goRight(iter));
}

template <
    typename TOnFoundFunctor,
    typename TString,
    typename TTreeIterator,
    typename TErrors>
inline void
approximateStringSearch(
    TOnFoundFunctor &onFoundFunctor,
    TString const &string,
    TTreeIterator &iter,
    TErrors errorsLeft)
{
    if (length(string) <= errorsLeft)
    {
        onFoundFunctor(string, iter);
        return;
    }
    _approximateStringSearch(onFoundFunctor, string, 0u, iter, repLength(iter), errorsLeft);
}

template <
    typename TOnFoundFunctor,
    typename TTreeIteratorA,
    typename TTreeIteratorB,
    typename TErrors>
inline void
approximateTreeSearch(
    TOnFoundFunctor &onFoundFunctor,
    TTreeIteratorA const &iterA,
    TTreeIteratorB const &iterB,
    TErrors errorsLeft)
{
    _approximateTreeSearch<true,true>(
        onFoundFunctor,
        iterA,
        repLength(iterA),
        iterB,
        repLength(iterB),
        errorsLeft);
}

}
#endif

