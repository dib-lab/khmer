// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
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

#ifndef INDEX_FM_RANK_DICTIONARY_NAIVE_H_
#define INDEX_FM_RANK_DICTIONARY_NAIVE_H_

namespace seqan {

// ============================================================================
// Tags
// ============================================================================

// ----------------------------------------------------------------------------
// Tag Naive
// ----------------------------------------------------------------------------

template <typename TSpec = void, typename TConfig = RDConfig<> >
struct Naive {};

// ============================================================================
// Metafunctions
// ============================================================================

// ----------------------------------------------------------------------------
// Metafunction Fibre
// ----------------------------------------------------------------------------

template <typename TValue, typename TSpec, typename TConfig>
struct Fibre<RankDictionary<TValue, Naive<TSpec, TConfig> >, FibreRanks>
{
    typedef RankDictionary<TValue, Naive<TSpec, TConfig> >          TRankDictionary_;
    typedef typename Size<TRankDictionary_>::Type                   TSize_;
    typedef typename DefaultIndexStringSpec<TRankDictionary_>::Type TFibreSpec_;

    typedef String<TSize_, TFibreSpec_>                             Type;
};

// ============================================================================
// Classes
// ============================================================================

// ----------------------------------------------------------------------------
// Class Naive RankDictionary
// ----------------------------------------------------------------------------

template <typename TValue, typename TSpec, typename TConfig>
struct RankDictionary<TValue, Naive<TSpec, TConfig> >
{
    // ------------------------------------------------------------------------
    // Fibres
    // ------------------------------------------------------------------------

    typename Fibre<RankDictionary, FibreRanks>::Type    ranks;

    // ------------------------------------------------------------------------
    // Constructors
    // ------------------------------------------------------------------------

    RankDictionary() {}

    template <typename TText>
    RankDictionary(TText const & text)
    {
        createRankDictionary(*this, text);
    }
};

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function getRank()
// ----------------------------------------------------------------------------

template <typename TSpec, typename TConfig, typename TPos>
SEQAN_HOST_DEVICE inline typename Size<RankDictionary<bool, Naive<TSpec, TConfig> > const>::Type
getRank(RankDictionary<bool, Naive<TSpec, TConfig> > const & dict, TPos pos, bool c = true)
{
    typedef RankDictionary<bool, Naive<TSpec, TConfig> > const                       TRankDictionary;
    typedef typename Fibre<TRankDictionary, FibreRanks>::Type               TFibreRanks;
    typedef typename Iterator<TFibreRanks, Standard>::Type                  TFibreRanksIterator;

    TFibreRanksIterator ranksBegin = begin(dict.ranks, Standard());
    TFibreRanksIterator ranksEnd   = end(dict.ranks, Standard());
    TFibreRanksIterator ranksIt    = ranksBegin;

    for (; ranksIt != ranksEnd && value(ranksIt) <= pos; ++ranksIt) ;

    return c ? ranksIt - ranksBegin : pos + 1 - (ranksIt - ranksBegin);
}

// ----------------------------------------------------------------------------
// Function getValue()
// ----------------------------------------------------------------------------

template <typename TSpec, typename TConfig, typename TPos>
SEQAN_HOST_DEVICE inline typename Value<RankDictionary<bool, Naive<TSpec, TConfig> > >::Type
getValue(RankDictionary<bool, Naive<TSpec, TConfig> > & dict, TPos pos)
{
    typedef RankDictionary<bool, Naive<TSpec, TConfig> > const                       TRankDictionary;
    typedef typename Fibre<TRankDictionary, FibreRanks>::Type               TFibreRanks;
    typedef typename Iterator<TFibreRanks, Standard>::Type                  TFibreRanksIterator;

    TFibreRanksIterator ranksBegin = begin(dict.ranks, Standard());
    TFibreRanksIterator ranksEnd   = end(dict.ranks, Standard());
    TFibreRanksIterator ranksIt    = ranksBegin;

    for (; ranksIt != ranksEnd && value(ranksIt) < pos; ++ranksIt) ;

    return ranksIt != ranksEnd && value(ranksIt) == pos;
}

template <typename TSpec, typename TConfig, typename TPos>
SEQAN_HOST_DEVICE inline typename Value<RankDictionary<bool, Naive<TSpec, TConfig> > const>::Type
getValue(RankDictionary<bool, Naive<TSpec, TConfig> > const & dict, TPos pos)
{
    typedef RankDictionary<bool, Naive<TSpec, TConfig> > const                       TRankDictionary;
    typedef typename Fibre<TRankDictionary, FibreRanks>::Type               TFibreRanks;
    typedef typename Iterator<TFibreRanks, Standard>::Type                  TFibreRanksIterator;

    TFibreRanksIterator ranksBegin = begin(dict.ranks, Standard());
    TFibreRanksIterator ranksEnd   = end(dict.ranks, Standard());
    TFibreRanksIterator ranksIt    = ranksBegin;

    for (; ranksIt != ranksEnd && value(ranksIt) < pos; ++ranksIt) ;

    return ranksIt != ranksEnd && value(ranksIt) == pos;
}

// ----------------------------------------------------------------------------
// Function setValue()
// ----------------------------------------------------------------------------

template <typename TValue, typename TSpec, typename TConfig, typename TPos, typename TChar>
inline void setValue(RankDictionary<TValue, Naive<TSpec, TConfig> > & dict, TPos pos, TChar c)
{
//    SEQAN_ASSERT_GT(pos, (TPos)back(dict.ranks));

    if (c == false) return;

    appendValue(dict.ranks, pos);
}

// ----------------------------------------------------------------------------
// Function appendValue()
// ----------------------------------------------------------------------------
// NOTE(esiragusa): Better not to have appendValue() - it is not efficient - and thus neither length().

template <typename TSpec, typename TConfig, typename TChar, typename TExpand>
inline void appendValue(RankDictionary<bool, Naive<TSpec, TConfig> > & dict, TChar c, Tag<TExpand> const tag)
{
    if (c == false) return;

// NOTE(esiragusa): RankDictionary's resize() is desabled.
//    resize(dict, length(dict) + 1, tag);
    resize(dict.ranks, length(dict) + 1, tag);
    setValue(dict, length(dict) - 1, c);
}

// ----------------------------------------------------------------------------
// Function updateRanks()
// ----------------------------------------------------------------------------

template <typename TValue, typename TSpec, typename TConfig>
inline void updateRanks(RankDictionary<TValue, Naive<TSpec, TConfig> > & /* dict */) {}

// ----------------------------------------------------------------------------
// Function length()
// ----------------------------------------------------------------------------

template <typename TValue, typename TSpec, typename TConfig>
inline typename Size<RankDictionary<TValue, Naive<TSpec, TConfig> > >::Type
length(RankDictionary<TValue, Naive<TSpec, TConfig> > const & dict)
{
    return length(dict.ranks);
}

// ----------------------------------------------------------------------------
// Function reserve()
// ----------------------------------------------------------------------------

template <typename TValue, typename TSpec, typename TConfig, typename TSize, typename TExpand>
inline typename Size<RankDictionary<TValue, Naive<TSpec, TConfig> > >::Type
reserve(RankDictionary<TValue, Naive<TSpec, TConfig> > & dict, TSize newCapacity, Tag<TExpand> const tag)
{
   return reserve(dict.ranks, newCapacity, tag);
}

// ----------------------------------------------------------------------------
// Function resize()
// ----------------------------------------------------------------------------
// NOTE(esiragusa): disabled because LF::_createBwt() resizes the rank dict to the bwt length.

template <typename TValue, typename TSpec, typename TConfig, typename TSize, typename TExpand>
inline typename Size<RankDictionary<TValue, Naive<TSpec, TConfig> > >::Type
resize(RankDictionary<TValue, Naive<TSpec, TConfig> > & dict, TSize /* newLength */, Tag<TExpand> const /* tag */)
{
    return length(dict);
//    return resize(dict.ranks, newLength, tag);
}

}

#endif  // INDEX_FM_RANK_DICTIONARY_NAIVE_H_
