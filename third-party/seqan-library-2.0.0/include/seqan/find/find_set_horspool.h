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

#ifndef SEQAN_HEADER_FIND_SETHORSPOOL_H
#define SEQAN_HEADER_FIND_SETHORSPOOL_H

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////
// Set Horspool Algorithm
//////////////////////////////////////////////////////////////////////////////

/*!
 * @class SetHorspoolPattern
 * @extends Pattern
 * @headerfile <seqan/find.h>
 *
 * @brief Multiple exact string matching using set horspool algorithm.
 *
 * @signature template <typename TNeedle>
 *            class Pattern<TNeedle, SetHorspool>;
 *
 * @tparam TNeedle The needle type, a string of keywords.  Types: @link ContainerConcept @endlink.
 *
 * The types of all keywords in the needle and the haystack have to match.
 */

struct SetHorspool_;
typedef Tag<SetHorspool_> SetHorspool;

//////////////////////////////////////////////////////////////////////////////

template <typename TNeedle>
class Pattern<TNeedle, SetHorspool> {
//____________________________________________________________________________
private:
    Pattern(Pattern const& other);
    Pattern const& operator=(Pattern const & other);

//____________________________________________________________________________
public:
    typedef typename Size<TNeedle>::Type TSize;
    typedef typename Value<TNeedle>::Type TValue;
    typedef typename Value<TValue>::Type TAlphabet;
    typedef Graph<Automaton<TAlphabet> > TGraph;
    typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;

    Holder<TNeedle> data_host;
    Graph<Automaton<TAlphabet> > data_reverseTrie;  // Search trie
    String<String<TSize> > data_terminalStateMap;
    String<TSize> data_dMap;    // Jump table
    TSize data_lmin;
    String<TSize> data_endPositions;    // All remaining keyword indices
    TSize data_keywordIndex;            // Current keyword that produced a hit
    TSize data_needleLength;            // Last length of needle to reposition finder
    TVertexDescriptor data_lastState;   // Last state in the trie

//____________________________________________________________________________

    Pattern() {
    }

    template <typename TNeedle2>
    Pattern(TNeedle2 const & ndl)
    {
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
struct Host< Pattern<TNeedle, SetHorspool> >
{
    typedef TNeedle Type;
};

template <typename TNeedle>
struct Host< Pattern<TNeedle, SetHorspool> const>
{
    typedef TNeedle const Type;
};

//////////////////////////////////////////////////////////////////////////////
// Functions
//////////////////////////////////////////////////////////////////////////////

template <typename TNeedle, typename TNeedle2>
void setHost (Pattern<TNeedle, SetHorspool> & me, TNeedle2 const & needle) {
    SEQAN_CHECKPOINT
    typedef typename Value<TNeedle>::Type TKeyword;
    typedef typename Size<TKeyword>::Type TSize;
    typedef typename Value<TKeyword>::Type TAlphabet;

    // clean-up
    clear(me.data_reverseTrie);
    clear(me.data_terminalStateMap);
    clear(me.data_endPositions);
    clear(me.data_dMap);
    me.data_lmin=0;

    // Create Trie
    createTrieOnReverse(me.data_reverseTrie,me.data_terminalStateMap,needle);
    assignRoot(me.data_reverseTrie,0);
    setValue(me.data_host, needle);

    // Create jump map
    TSize alphabet_size = ValueSize<TAlphabet>::VALUE;
    resize(me.data_dMap, alphabet_size);
    me.data_lmin = _getInfinity<TSize>();
    typename Iterator<TNeedle2 const, Rooted>::Type it = begin(needle);
    for(;!atEnd(it);goNext(it)) {
        TSize tmp = length(*it);
        if (tmp<me.data_lmin) me.data_lmin = tmp;
    }
    for(TSize i=0;i<alphabet_size;++i) {
        me.data_dMap[i]=me.data_lmin;
    }
    goBegin(it);
    for(;!atEnd(it);goNext(it)) {
        for(TSize pos = 0;pos < length(*it) - 1; ++pos) {
            TSize ind = ordValue((TAlphabet)(*it)[pos]);
            if ((length(*it)- 1 - pos) < me.data_dMap[ind]) {
                me.data_dMap[ind] = (length(*it) - 1 - pos);
            }
        }
    }

    /*
    fstream strm;
    strm.open(TEST_PATH "my_trie.dot", ios_base::out | ios_base::trunc);
    String<String<char> > nodeMap;
    _createTrieNodeNames(me.data_reverseTrie, me.data_terminalStateMap, nodeMap);
    String<String<char> > edgeMap;
    _createEdgeNames(me.data_reverseTrie,edgeMap);
    writeRecords(strm,me.data_reverseTrie,nodeMap,edgeMap,DotDrawing());
    strm.close();
    */
}

template <typename TNeedle, typename TNeedle2>
void setHost (Pattern<TNeedle, SetHorspool> & me, TNeedle2 & needle)
{
    setHost(me, reinterpret_cast<TNeedle2 const &>(needle));
}

//____________________________________________________________________________


template <typename TNeedle>
inline void _patternInit (Pattern<TNeedle, SetHorspool> & me)
{
SEQAN_CHECKPOINT
    clear(me.data_endPositions);
    me.data_keywordIndex = 0;
    me.data_lastState = getRoot(me.data_reverseTrie);
}


//____________________________________________________________________________

template <typename TNeedle>
inline typename Host<Pattern<TNeedle, SetHorspool>const>::Type &
host(Pattern<TNeedle, SetHorspool> & me)
{
SEQAN_CHECKPOINT
    return value(me.data_host);
}

template <typename TNeedle>
inline typename Host<Pattern<TNeedle, SetHorspool>const>::Type &
host(Pattern<TNeedle, SetHorspool> const & me)
{
SEQAN_CHECKPOINT
    return value(me.data_host);
}

//____________________________________________________________________________


template <typename TNeedle>
inline typename Size<TNeedle>::Type
position(Pattern<TNeedle, SetHorspool> & me)
{
    return me.data_keywordIndex;
}


template <typename TFinder, typename TNeedle>
inline bool find(TFinder & finder, Pattern<TNeedle, SetHorspool> & me) {
    SEQAN_CHECKPOINT
    typedef typename Value<TNeedle>::Type TKeyword;
    typedef typename Size<TKeyword>::Type TSize;
    typedef typename Value<TKeyword>::Type TAlphabet;
    typedef Graph<Automaton<TAlphabet> > TGraph;
    typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;

    TVertexDescriptor current = getRoot(me.data_reverseTrie);

    // Process left-over hits
    if ((!empty(finder)) &&
        (!empty(me.data_endPositions))) {
        finder += me.data_needleLength;
        current = me.data_lastState;
        me.data_keywordIndex = me.data_endPositions[length(me.data_endPositions)-1];
        me.data_needleLength = length(getValue(host(me), me.data_keywordIndex))-1;
        if (length(me.data_endPositions) > 1) resize(me.data_endPositions, (length(me.data_endPositions)-1));
        else clear(me.data_endPositions);
        me.data_lastState = current;
        finder -= me.data_needleLength;
        return true;
    }

    TVertexDescriptor nilVal = getNil<TVertexDescriptor>();
    TSize j = 0;
    if (empty(finder)) {
        _patternInit(me);
        _finderSetNonEmpty(finder);
        finder += me.data_lmin - 1;
    } else {
        finder += me.data_needleLength;
        j = me.data_needleLength + 1;
        current = me.data_lastState;
    }

    TSize haystackLength = length(container(finder));
    bool oldMatch = true;
    // Do not change to !atEnd(finder) because of jump map!
    while(position(finder) < haystackLength) {
        while ((position(finder)>=j) &&
                (getSuccessor(me.data_reverseTrie, current, *(finder-j))!= nilVal))
        {
            me.data_endPositions = getProperty(me.data_terminalStateMap,current);
            if ((!oldMatch) && (!empty(me.data_endPositions))) break;
            current = getSuccessor(me.data_reverseTrie, current, *(finder-j));
            if (current == nilVal) break;
            ++j;
            oldMatch = false;
        }
        me.data_endPositions = getProperty(me.data_terminalStateMap,current);
        if ((!oldMatch) &&
            (!empty(me.data_endPositions)))
        {
            me.data_keywordIndex = me.data_endPositions[length(me.data_endPositions)-1];
            me.data_needleLength = length(getValue(host(me), me.data_keywordIndex))-1;
            if (length(me.data_endPositions) > 1) resize(me.data_endPositions, length(me.data_endPositions)-1);
            else clear(me.data_endPositions);
            me.data_lastState = current;
            finder -= me.data_needleLength;
            return true;
        }
        oldMatch = false;
        TSize ind = ordValue(*finder);
        setPosition(finder, position(finder) + getValue(me.data_dMap, ind));
        j = 0;
        current = getRoot(me.data_reverseTrie);
    }
    return false;
}

}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_FIND_SETHORSPOOL_H
