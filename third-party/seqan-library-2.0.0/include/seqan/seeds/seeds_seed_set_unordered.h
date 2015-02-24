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
// Author: Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>
// ==========================================================================
// The Unordered specialization of the class SeedSet.  Seeds are stored
// in a string and not kept in a particular order.
// ==========================================================================

#ifndef SEQAN_SEEDS_SEEDS_SEED_SET_UNORDERED_H_
#define SEQAN_SEEDS_SEEDS_SEED_SET_UNORDERED_H_

#include <cmath>

namespace seqan {

// ===========================================================================
// Enums, Tags, Classes, Specializations
// ===========================================================================

// ---------------------------------------------------------------------------
// Class Unordered SeedSet
// ---------------------------------------------------------------------------

struct Unordered_;
typedef Tag<Unordered_> Unordered;

// TODO(holtgrew): Maybe allow iterating over seeds that have reached a certain quality (length/score).

template <typename TSeed>
class SeedSet<TSeed, Unordered>
{
public:
    typedef typename SeedScore<TSeed>::Type TScoreValue_;
    typedef typename Size<TSeed>::Type TSize_;
    typedef LessBeginDiagonal<TSeed> TSeedCmp_;
    typedef std::multiset<TSeed, TSeedCmp_> TSet_;

    TSet_ _seeds;

    TScoreValue_ _minScore;
    TSize_ _minSeedSize;

    SeedSet() : _minScore(0), _minSeedSize(0)
    {}
};

// ===========================================================================
// Metafunctions
// ===========================================================================

// ---------------------------------------------------------------------------
// Metafunction Position
// ---------------------------------------------------------------------------

template <typename TSeed>
struct Position<SeedSet<TSeed, Unordered> >
{
    typedef SeedSet<TSeed, Unordered> TSeedSet_;
    typedef String<typename TSeedSet_::TSeed> TSeedString_;
    typedef typename Position<TSeedString_>::Type Type;
};

template <typename TSeed>
struct Position<SeedSet<TSeed, Unordered> const> : Position<SeedSet<TSeed, Unordered> >
{};

// ---------------------------------------------------------------------------
// Metafunction Size
// ---------------------------------------------------------------------------

template <typename TSeed>
struct Size<SeedSet<TSeed, Unordered> >
{
    typedef SeedSet<TSeed, Unordered> TSeedSet_;
    typedef String<TSeed> TSeedString_;
    typedef typename Size<TSeedString_>::Type Type;
};

template <typename TSeed>
struct Size<SeedSet<TSeed, Unordered> const> : Size<SeedSet<TSeed, Unordered> >
{};

// ---------------------------------------------------------------------------
// Metafunction Value
// ---------------------------------------------------------------------------

template <typename TSeed>
struct Value<SeedSet<TSeed, Unordered> >
{
    typedef TSeed Type;
};

template <typename TSeed>
struct Value<SeedSet<TSeed, Unordered> const>
{
    typedef typename Value<SeedSet<TSeed, Unordered> >::Type const Type;
};

// ---------------------------------------------------------------------------
// Metafunction Reference
// ---------------------------------------------------------------------------

template <typename TSeed>
struct Reference<SeedSet<TSeed, Unordered> >
{
    typedef typename Value<SeedSet<TSeed, Unordered> >::Type & Type;
};

// ---------------------------------------------------------------------------
// Metafunction Iterator
// ---------------------------------------------------------------------------

template <typename TSeed>
struct Iterator<SeedSet<TSeed, Unordered>, Standard>
{
    typedef SeedSet<TSeed, Unordered> TSeedSet_;
    typedef typename TSeedSet_::TSet_ TMultiSet_;
    typedef Iter<TMultiSet_, StdIteratorAdaptor> Type;
};

template <typename TSeed>
struct Iterator<SeedSet<TSeed, Unordered> const, Standard>
{
    typedef SeedSet<TSeed, Unordered> const TSeedSet_;
    typedef typename TSeedSet_::TSet_ const TMultiSet_;
    typedef Iter<TMultiSet_, StdIteratorAdaptor> Type;
};

// ===========================================================================
// Functions
// ===========================================================================

// ---------------------------------------------------------------------------
// Function length()
// ---------------------------------------------------------------------------

// Standard Container Functions

template <typename TSeed>
inline typename Size<SeedSet<TSeed, Unordered> >::Type
length(SeedSet<TSeed, Unordered> & seedSet)
{
    return seedSet._seeds.size();
}

template <typename TSeed>
inline typename Size<SeedSet<TSeed, Unordered> const>::Type
length(SeedSet<TSeed, Unordered> const & seedSet)
{
    return seedSet._seeds.size();
}

// ---------------------------------------------------------------------------
// Function begin()
// ---------------------------------------------------------------------------

template <typename TSeed>
inline typename Iterator<SeedSet<TSeed, Unordered> >::Type
begin(SeedSet<TSeed, Unordered> & seedSet, Standard const &)
{
    return seedSet._seeds.begin();
}

template <typename TSeed>
inline typename Iterator<SeedSet<TSeed, Unordered> const>::Type
begin(SeedSet<TSeed, Unordered> const & seedSet, Standard const &)
{
    return seedSet._seeds.begin();
}

// ---------------------------------------------------------------------------
// Function end()
// ---------------------------------------------------------------------------

template <typename TSeed>
inline typename Iterator<SeedSet<TSeed, Unordered> >::Type
end(SeedSet<TSeed, Unordered> & seedSet, Standard const &)
{
    return seedSet._seeds.end();
}

template <typename TSeed>
inline typename Iterator<SeedSet<TSeed, Unordered> const>::Type
end(SeedSet<TSeed, Unordered> const & seedSet, Standard const &)
{
    return seedSet._seeds.end();
}

// ---------------------------------------------------------------------------
// Function front()
// ---------------------------------------------------------------------------

template <typename TSeed>
inline typename Reference<typename Iterator<SeedSet<TSeed, Unordered>, Standard>::Type >::Type
front(SeedSet<TSeed, Unordered> & seedSet)
{
    return *seedSet._seeds.begin();
}

template <typename TSeed>
inline typename Reference<typename Iterator<SeedSet<TSeed, Unordered> const, Standard>::Type >::Type
front(SeedSet<TSeed, Unordered> const & seedSet)
{
    return *seedSet._seeds.begin();
}

// ---------------------------------------------------------------------------
// Function back()
// ---------------------------------------------------------------------------

template <typename TSeed>
inline typename Reference<typename Iterator<SeedSet<TSeed, Unordered> , Standard>::Type >::Type
back(SeedSet<TSeed, Unordered> & seedSet)
{
    return *seedSet._seeds.rbegin();
}

template <typename TSeed>
inline typename Reference<typename Iterator<SeedSet<TSeed, Unordered> const, Standard>::Type >::Type
back(SeedSet<TSeed, Unordered> const & seedSet)
{
    return *seedSet._seeds.rbegin();
}

// SeedSet Functions

// ---------------------------------------------------------------------------
// Function _findSeedForCombination()
// ---------------------------------------------------------------------------

// TODO(holtgrew): Add bulk-addSeeds functions.

template <typename TSeedIter, typename TSeed, typename TDistanceThreshold, typename TBandwidth, typename TCombination>
bool
_findSeedForCombination(
        TSeedIter & mergePartner,
        bool & seedIsOnTheLeft,
        SeedSet<TSeed, Unordered> & seedSet,
        typename Value<SeedSet<TSeed, Unordered> >::Type const & seed,
        TDistanceThreshold const & maxDistance,
        TBandwidth const & bandwidth,
        TCombination const & tag)
{
    // Iterate over all seeds and search for the first one in this
    // arbitrary order that is combineable with parameter seed within
    // a maximal diagonal distance maxDistance.  We allow either seed
    // to be the left one.
    //
    // TODO(holtgrew): Search for *closest* overlapping one instead!
    for (TSeedIter it = seedSet._seeds.begin(); it != seedSet._seeds.end(); ++it)
    {
        if (_seedsCombineable(*it, seed, maxDistance, bandwidth, tag))
        {
//            std::cout << "Combineable: " << (*value(it)) << " and " << seed << std::endl;
            // seed is to be merged into *it.
            mergePartner = it;
            seedIsOnTheLeft = false;
            return true;
        }
        else if (_seedsCombineable(seed, *it, maxDistance, bandwidth, tag))
        {
//            std::cout << "Combineable: " << seed << " and " << (*value(it)) << std::endl;
            // *it is to be merged into seed.
            mergePartner = it;
            seedIsOnTheLeft = true;
            return true;
        }
    }

    // Found no seed to combine with.
    return false;
}

// ---------------------------------------------------------------------------
// Function addSeed()
// ---------------------------------------------------------------------------

template <typename TSeed, typename TDistanceThreshold, typename TBandwidth, typename TScoreValue, typename TSequence0, typename TSequence1, typename TCombination>
inline bool
addSeed(SeedSet<TSeed, Unordered> & seedSet,
        TSeed const & seed,
        TDistanceThreshold const & maxDiagDist,
        TBandwidth const & bandwidth,
        Score<TScoreValue, Simple> const & scoringScheme,
        TSequence0 const & sequence0,
        TSequence1 const & sequence1,
        TCombination const & tag)
{
    SEQAN_CHECKPOINT;

    typedef SeedSet<TSeed, Unordered> TSeedSet;
    typedef typename TSeedSet::TSet_ TSet;
    typedef typename TSet::iterator TSeedIterator;

    // Try to find a seed for recombination.
    TSeedIterator it;
    bool seedIsOnTheLeft = false;
    bool foundSeed = _findSeedForCombination(it, seedIsOnTheLeft, seedSet, seed, maxDiagDist, bandwidth, tag);

    // If we could find a seed: Combine them.
    if (foundSeed)
    {
        TSeed left;
        if (!seedIsOnTheLeft)
        {
            left = *it;
            _combineSeeds(left, seed, scoringScheme, sequence0, sequence1, tag);
        }
        else
        {
            left = seed;
            _combineSeeds(left, *it, scoringScheme, sequence0, sequence1, tag);
        }

        seedSet._seeds.erase(it);
        seedSet._seeds.insert(left);
        return true;
    }
    return false;
}

// TODO(holtgrew): Score not needed for Merge!
template <typename TSeed, typename TDistanceThreshold>
inline bool
addSeed(SeedSet<TSeed, Unordered> & seedSet,
        TSeed const & seed,
        TDistanceThreshold const & maxDiagDist,
        Merge const &)
{
    return addSeed(seedSet, seed, maxDiagDist, 0, Score<int, Simple>(), Nothing(), Nothing(), Merge());
}

template <typename TSeed, typename TDistanceThreshold, typename TScoreValue>
inline bool
addSeed(SeedSet<TSeed, Unordered> & seedSet,
        TSeed const & seed,
        TDistanceThreshold const & maxDiagDist,
        Score<TScoreValue, Simple> const & scoringScheme,
        SimpleChain const &)
{
    return addSeed(seedSet, seed, maxDiagDist, 0, scoringScheme, Nothing(), Nothing(), SimpleChain());
}

template <typename TSeed>
inline bool
addSeed(SeedSet<TSeed, Unordered> & seedSet,
        TSeed const & seed,
        Single const &)
{
    seedSet._seeds.insert(seed);
    return true;    // Always returns true.
}

}  // namespace seqan

#endif  // SEQAN_SEEDS_SEEDS_SEED_SET_UNORDERED_H_

