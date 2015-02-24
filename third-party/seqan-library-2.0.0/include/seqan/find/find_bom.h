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
// Author: Tobias Rausch <rausch@embl.de>
// ==========================================================================

#ifndef SEQAN_HEADER_FIND_BOM_H
#define SEQAN_HEADER_FIND_BOM_H

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////
// BomAlgo
//////////////////////////////////////////////////////////////////////////////

/*!
 * @class BfamPattern
 * @extends Pattern
 * @headerfile <seqan/find.h>
 * @brief Backward Factor Automaton Matching algorithm.
 *
 * @signature template <typename TNeedle, typename TAutomaton>
 *            class Pattern<TNeedle, Bfam<TAutomaton> >;
 *
 * @tparam TAutomaton A tag that specifies the used automaton. Default: <tt>Bfam&lt;Oracle&gt;.
 * @tparam TNeedle    The needle type. Types: String
 *
 * To be used in combination with the default specialization of @link Finder @endlink.
 *
 * @see MultiBfamPattern
 */

/*!
 * @class OracleBfamPattern
 * @extends BfamPattern
 * @headerfile <seqan/find.h>
 * @brief Backward Oracle Matching algorithm.
 *
 * @signature template <typename TNeedle>
 *            class Pattern<TNeedle, Bfam<Oracle> >;
 *
 * @tparam TNeedle The needle type. Types: String
 *
 * To be used in combination with the default specialization of @link Finder @endlink.
 *
 * @see TrieBfamPattern
 * @see OracleMultiBfamPattern
 */

/*!
 * @class TrieBfamPattern
 * @extends BfamPattern
 * @headerfile <seqan/find.h>
 * @brief Backward Suffix Trie Matching algorithm.
 *
 * @signature template <typename TNeedle>
 *            class Pattern<TNeedle, Bfam<Trie> >;
 *
 * @tparam TNeedle The needle type. Types: String
 *
 * To be used in combination with the default specialization of @link Finder @endlink.
 *
 * @see OracleBfamPattern
 */

struct Oracle; //Oracle Tag => "BOM"
struct Trie {}; //Trie Tag => "BTM"

template <typename TSpec = Oracle>
struct Bfam; //backward factor automaton searching

typedef Bfam<Oracle> BomAlgo; //deprecated, still there for compatibility reasons

//////////////////////////////////////////////////////////////////////////////

template <typename TNeedle, typename TSpec>
class Pattern<TNeedle, Bfam<TSpec> > {
//____________________________________________________________________________
public:
    typedef typename Value<TNeedle>::Type TAlphabet;
    typedef typename Size<TNeedle>::Type TSize;
    Holder<TNeedle> data_host;
    TSize needleLength;
    TSize haystackLength;
    TSize step;
    Graph<Automaton<TAlphabet, void, WithoutEdgeId> > automaton;

//____________________________________________________________________________

    Pattern() {
    }

    template <typename TNeedle2>
    Pattern(TNeedle2 const & ndl)
    {
SEQAN_CHECKPOINT
        setHost(*this, ndl);
    }

    ~Pattern() {
        SEQAN_CHECKPOINT
    }
//____________________________________________________________________________
};

//////////////////////////////////////////////////////////////////////////////
// Host Metafunctions
//////////////////////////////////////////////////////////////////////////////

template <typename TNeedle>
struct Host< Pattern<TNeedle, BomAlgo> >
{
    typedef TNeedle Type;
};

template <typename TNeedle>
struct Host< Pattern<TNeedle, BomAlgo> const>
{
    typedef TNeedle const Type;
};


//////////////////////////////////////////////////////////////////////////////
// Functions
//////////////////////////////////////////////////////////////////////////////

//Bfam<Oracle>: BOM Algorithm
template <typename TNeedle, typename TNeedle2>
inline void
setHost (Pattern<TNeedle, Bfam<Oracle> > & me, TNeedle2 const& needle)
{
    SEQAN_CHECKPOINT
    me.needleLength = length(needle);
    clear(me.automaton);
    createOracleOnReverse(me.automaton,needle);
    setValue(me.data_host, needle);
}

//Bfam<Trie>: BTM Algorithm (the same as BOM, but with an trie)
template <typename TNeedle, typename TNeedle2>
inline void
setHost (Pattern<TNeedle, Bfam<Trie> > & me, TNeedle2 const& needle)
{
    SEQAN_CHECKPOINT;
    typedef typename Position<TNeedle>::Type TPosition;
    me.needleLength = length(needle);
    clear(me.automaton);

    String<String<TPosition> > terminal_state_map; //dummy
    typedef typename Value<TNeedle2 const>::Type TValue;
    String<TValue> reverse_string = needle;
    reverse(reverse_string);

    createSuffixTrie(me.automaton, terminal_state_map, reverse_string);

    setValue(me.data_host, needle);
}

template <typename TNeedle, typename TNeedle2, typename TSpec>
inline void
setHost (Pattern<TNeedle, Bfam<TSpec> > & me, TNeedle2 & needle)
{
    setHost(me, reinterpret_cast<TNeedle2 const &>(needle));
}

//____________________________________________________________________________


template <typename TNeedle, typename TSpec>
inline void _patternInit (Pattern<TNeedle, Bfam<TSpec> > & me)
{
SEQAN_CHECKPOINT
    me.step = 0;
}


//____________________________________________________________________________


template <typename TFinder, typename TNeedle, typename TSpec>
inline bool
find(TFinder & finder, Pattern<TNeedle, Bfam<TSpec> > & me)
{
    if (empty(finder)) {
        _patternInit(me);
        _setFinderLength(finder, length(needle(me)));
        _finderSetNonEmpty(finder);
        me.haystackLength = length(container(finder));
    } else
        finder+=me.step;

    if (me.haystackLength < me.needleLength) return false;
    typedef typename Value<TNeedle>::Type TAlphabet;
    typedef Graph<Automaton<TAlphabet> > TOracle;
    typedef typename Size<TNeedle>::Type TSize;
    typedef typename VertexDescriptor<TOracle>::Type TVertexDescriptor;
    TVertexDescriptor nilVal = getNil<TVertexDescriptor>();
    while (position(finder) <= me.haystackLength - me.needleLength) {
        TVertexDescriptor current = getRoot(me.automaton);
        TSize j = me.needleLength;
        while ((j>0) &&    (current != nilVal))
        {
            TAlphabet c = *(finder+(j-1));
            current = targetVertex(me.automaton, findEdge(me.automaton, current, c));
            --j;
        }
        if (current != nilVal) {
            me.step = j + 1;
            _setFinderEnd(finder, position(finder) + me.needleLength);
            return true;
        }
        finder += j + 1;
    }
    return false;
}

}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_FIND_SHIFTAND_H
