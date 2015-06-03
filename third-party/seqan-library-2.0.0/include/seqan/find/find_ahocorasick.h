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

#ifndef SEQAN_HEADER_FIND_AHOCORASICK_H
#define SEQAN_HEADER_FIND_AHOCORASICK_H

// TODO(holtgrew): Needles should be a StringSet<CharString>!

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////
// AhoCorasick
//////////////////////////////////////////////////////////////////////////////

/*!
 * @class AhoCorasickPattern
 * @extends Pattern
 * @headerfile <seqan/find.h>
 * @brief Multiple exact string matching using Aho-Corasick.
 *
 * @signature template <typename TNeedle>
 *            class Pattern<TNeedle, AhoCorasick>;
 *
 * @tparam TNeedle The needle type, a string of keywords.
 *
 * The types of the keywords in the needle container and the haystack have to match.
 *
 * Matching positions do not come in order because we report beginning positions of matches.
 *
 * Likewise, if multiple keywords match at a given position no pre-specified order is guaranteed.
 *
 * @section Examples
 *
 * The following example program searches for three needles (<tt>queries</tt>) in two haystack sequences (<tt>db</tt>)
 * using the Aho-Corasick algorithm.
 *
 * @include demos/find/finder_aho_corasick.cpp
 *
 * When executed, this program will create the following output.
 *
 * @code{.console}
 * DB      POS     ENDPOS  TEXT
 * 0       0       4       MARD
 * 0       3       7       DPLY
 * 1       1       6       VGGGG
 * 1       6       9       AAA
 * @endcode
 */

struct AhoCorasick_;
typedef Tag<AhoCorasick_> AhoCorasick;

//////////////////////////////////////////////////////////////////////////////

template <typename TNeedle>
class Pattern<TNeedle, AhoCorasick> {
//____________________________________________________________________________
private:
    Pattern(Pattern const& other);
    Pattern const& operator=(Pattern const & other);

//____________________________________________________________________________
public:
    typedef typename Size<TNeedle>::Type TSize;
    typedef typename Value<TNeedle>::Type TKeyword;
    typedef typename Value<TKeyword>::Type TAlphabet;
    typedef Graph<Automaton<TAlphabet> > TGraph;
    typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;

    Holder<TNeedle> data_host;
    String<TVertexDescriptor> data_supplyMap;
    String<String<TSize> > data_terminalStateMap;
    TGraph data_graph;

    // To restore the automaton after a hit
    String<TSize> data_endPositions;    // All remaining keyword indices
    TSize data_keywordIndex;            // Current keyword that produced a hit
    TSize data_needleLength;            // Last length of needle to reposition finder
    TVertexDescriptor data_lastState;   // Last state in the trie

//____________________________________________________________________________

    Pattern() : data_keywordIndex(0), data_needleLength(0)
    {}

    template <typename TNeedle2>
    Pattern(TNeedle2 const & ndl) : data_keywordIndex(0), data_needleLength(0)
    {
        setHost(*this, ndl);
    }

//____________________________________________________________________________
};


//////////////////////////////////////////////////////////////////////////////
// Host Metafunctions
//////////////////////////////////////////////////////////////////////////////

template <typename TNeedle>
struct Host< Pattern<TNeedle, AhoCorasick> >
{
    typedef TNeedle Type;
};

template <typename TNeedle>
struct Host< Pattern<TNeedle, AhoCorasick> const>
{
    typedef TNeedle const Type;
};


//////////////////////////////////////////////////////////////////////////////
// Functions
//////////////////////////////////////////////////////////////////////////////

template <typename TNeedle>
inline void
_createAcTrie(Pattern<TNeedle, AhoCorasick> & me)
{
    typedef typename Position<TNeedle>::Type TPosition;
    typedef typename Value<TNeedle>::Type TKeyword;
    typedef typename Value<TKeyword>::Type TAlphabet;
    typedef Graph<Automaton<TAlphabet> > TGraph;
    typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
    TVertexDescriptor nilVal = getNil<TVertexDescriptor>();

    // Create regular trie
    createTrie(me.data_graph,me.data_terminalStateMap, host(me));

    // Create parent map
    String<TVertexDescriptor> parentMap;
    String<TAlphabet> parentCharMap;
    resizeVertexMap(parentMap, me.data_graph);
    resizeVertexMap(parentCharMap, me.data_graph);
    for(TPosition i = 0;i<length(parentMap);++i) {
        assignProperty(parentMap, i, nilVal);
    }
    typedef typename Iterator<TGraph, EdgeIterator>::Type TEdgeIterator;
    TEdgeIterator itEd(me.data_graph);
    for(;!atEnd(itEd);goNext(itEd)) {
        assignProperty(parentMap, targetVertex(itEd), sourceVertex(itEd));
        assignProperty(parentCharMap, targetVertex(itEd), label(itEd));
    }

    // Build AC
    TVertexDescriptor root = getRoot(me.data_graph);
    resizeVertexMap(me.data_supplyMap, me.data_graph);
    assignProperty(me.data_supplyMap, root, nilVal);

    // Bfs Traversal
    typedef typename Iterator<TGraph, BfsIterator>::Type TBfsIterator;
    TBfsIterator it(me.data_graph,root);
    for(;!atEnd(it);goNext(it)) {
        if (atBegin(it)) continue;
        TVertexDescriptor parent = getProperty(parentMap, *it);
        TAlphabet sigma = getProperty(parentCharMap, *it);
        TVertexDescriptor down = getProperty(me.data_supplyMap, parent);
        while ((down != nilVal) &&
            (getSuccessor(me.data_graph, down, sigma) == nilVal))
        {
            down = getProperty(me.data_supplyMap, down);
        }
        if (down != nilVal) {
            assignProperty(me.data_supplyMap, *it, getSuccessor(me.data_graph, down, sigma));
            String<TPosition> endPositions = getProperty(me.data_terminalStateMap, getProperty(me.data_supplyMap, *it));
            if (!empty(endPositions)) {
                String<TPosition> endPositionsCurrent = getProperty(me.data_terminalStateMap, *it);
                typedef typename Iterator<String<TPosition>, Rooted >::Type TStringIterator;
                TStringIterator sit = begin(endPositions);
                for(;!atEnd(sit);goNext(sit)) {
                    appendValue(endPositionsCurrent, *sit);
                }
                assignProperty(me.data_terminalStateMap, *it, endPositionsCurrent);
            }
        } else {
            assignProperty(me.data_supplyMap, *it, root);
        }

    }
}


template <typename TNeedle, typename TNeedle2>
void setHost (Pattern<TNeedle, AhoCorasick> & me, TNeedle2 const & needle) {
    SEQAN_CHECKPOINT;
    SEQAN_ASSERT_NOT(empty(needle));
    setValue(me.data_host, needle);
    clear(me.data_graph);
    clear(me.data_supplyMap);
    clear(me.data_endPositions);
    clear(me.data_terminalStateMap);
    _createAcTrie(me);
    me.data_needleLength = 0;


    //fstream strm;
    //strm.open(TEST_PATH "my_trie.dot", ios_base::out | ios_base::trunc);
    //String<String<char> > nodeMap;
    //_createTrieNodeAttributes(me.data_graph, me.data_terminalStateMap, nodeMap);
    //String<String<char> > edgeMap;
    //_createEdgeAttributes(me.data_graph,edgeMap);
    //writeRecords(strm,me.data_graph,nodeMap,edgeMap,DotDrawing());
    //strm.close();
    // Supply links
    //for(unsigned int i=0;i<length(me.data_supplyMap);++i) {
    //    std::cout << i << "->" << getProperty(me.data_supplyMap,i) << std::endl;
    //}
}

template <typename TNeedle, typename TNeedle2>
inline void
setHost (Pattern<TNeedle, AhoCorasick> & me, TNeedle2 & needle)
{
    setHost(me, reinterpret_cast<TNeedle2 const &>(needle));
}

//____________________________________________________________________________


template <typename TNeedle>
inline void _patternInit (Pattern<TNeedle, AhoCorasick> & me)
{
SEQAN_CHECKPOINT
    clear(me.data_endPositions);
    me.data_keywordIndex = 0;
    me.data_lastState = getRoot(me.data_graph);
}


//____________________________________________________________________________


template <typename TNeedle>
inline typename Size<TNeedle>::Type
position(Pattern<TNeedle, AhoCorasick> & me)
{
    return me.data_keywordIndex;
}


template <typename TFinder, typename TNeedle>
inline bool find(TFinder & finder, Pattern<TNeedle, AhoCorasick> & me) {
    typedef typename Value<TNeedle>::Type TKeyword;
    typedef typename Value<TKeyword>::Type TAlphabet;
    typedef Graph<Automaton<TAlphabet> > TGraph;
    typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;

    if (empty(finder)) {
        _patternInit(me);
        _finderSetNonEmpty(finder);
    } else {
        finder += me.data_needleLength;
        ++finder; // Set forward the finder
    }

    // Process left-over hits
    if (!empty(me.data_endPositions)) {
        --finder; // Set back the finder
        me.data_keywordIndex = me.data_endPositions[length(me.data_endPositions)-1];
        me.data_needleLength = length(value(host(me), me.data_keywordIndex))-1;
        if (length(me.data_endPositions) > 1) resize(me.data_endPositions, (length(me.data_endPositions)-1));
        else clear(me.data_endPositions);
        finder -= me.data_needleLength;
        _setFinderLength(finder, me.data_needleLength+1);
        _setFinderEnd(finder, position(finder)+length(finder));
        return true;
    }

    TVertexDescriptor current = me.data_lastState;
    TVertexDescriptor nilVal = getNil<TVertexDescriptor>();
    while (!atEnd(finder)) {
        while ((getSuccessor(me.data_graph, current, *finder) == nilVal) &&
            (getProperty(me.data_supplyMap, current) != nilVal))
        {
            current = getProperty(me.data_supplyMap,current);
        }
        if (getSuccessor(me.data_graph, current, *finder) != nilVal) {
            current = getSuccessor(me.data_graph, current, *finder);
        }
        else {
            current = getRoot(me.data_graph);
        }
        me.data_endPositions = getProperty(me.data_terminalStateMap,current);
        if (!empty(me.data_endPositions)) {
            me.data_keywordIndex = me.data_endPositions[length(me.data_endPositions)-1];
            me.data_needleLength = length(value(host(me), me.data_keywordIndex))-1;
            if (length(me.data_endPositions) > 1) resize(me.data_endPositions, length(me.data_endPositions)-1);
            else clear(me.data_endPositions);
            me.data_lastState = current;
            finder -= me.data_needleLength;
            _setFinderLength(finder, me.data_needleLength+1);
            _setFinderEnd(finder, position(finder)+length(finder));
            return true;
        }
        ++finder;
    }
    return false;
}

}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_FIND_AHOCORASICK_H
