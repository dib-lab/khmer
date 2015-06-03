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

// TODO(holtgrew): We should probably specialize Automaton for this.

#ifndef SEQAN_HEADER_GRAPH_IMPL_TRIE_H
#define SEQAN_HEADER_GRAPH_IMPL_TRIE_H

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////
// Graph - Trie
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////
// INTERNAL FUNCTIONS
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

template<typename TAlphabet, typename TCargo, typename TSpec, typename TTerminalStateMap, typename TKeyword, typename TPos>
inline void
_addStringToTrie(Graph<Automaton<TAlphabet, TCargo, TSpec> >& g,
                 TTerminalStateMap& terminalStateMap,
                 TKeyword const& str,
                 TPos const& keywordIndex)
{
    typedef Graph<Automaton<TAlphabet, TCargo, TSpec> > TGraph;
    typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;

    TVertexDescriptor current = getRoot(g);
    TVertexDescriptor nilVal = getNil<TVertexDescriptor>();
    typename Iterator<TKeyword const, Rooted>::Type sIt = begin(str);
    for(;!atEnd(sIt);goNext(sIt)) {
        if (getSuccessor(g, current, *sIt) == nilVal) break;
        current = getSuccessor(g, current, *sIt);
    }
    for(;!atEnd(sIt);goNext(sIt)) {
        TVertexDescriptor newState = addVertex(g);
        resize(terminalStateMap, numVertices(g), Generous());
        assignProperty(terminalStateMap,newState,String<TPos>());
        addEdge(g,current,newState,*sIt);
        current = newState;
    }
    String<TPos> tmp = getProperty(terminalStateMap,current);
    appendValue(tmp, keywordIndex);
    assignProperty(terminalStateMap,current,tmp);
}

//////////////////////////////////////////////////////////////////////////////
// FUNCTIONS
//////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////

/*!
 * @fn Automaton#createTrie
 * @brief Creates a trie in an Automaton.
 *
 * @signature createTrie(g, terminalStateMap, keywords);
 *
 * @param[in,out] g        Automaton object to create the trie in.
 * @param[in]     terminalStateMap
 *                         An external property map.  The type must be <tt>String&lt;String&lt;unsigned&gt;&gt;</tt>
 *                         since a number of keywords can end in each vertex of the trie.  This is the case in the
 *                         Aho-Corasick algorithm if one pattern is a suffix of another pattern.  Hence, we must
 *                         associate with every vertex a set of indices that correspond to keywords.
 * @param[in]     keywords A set of strings.
 *
 * @see Automaton#createTrieOnReverse
 */

template<typename TAlphabet, typename TCargo, typename TSpec, typename TTerminalStateMap, typename TKeywords>
inline void
createTrie(Graph<Automaton<TAlphabet, TCargo, TSpec> >& g,
           TTerminalStateMap& terminalStateMap,
           TKeywords const& keywords)
{
    typedef Graph<Automaton<TAlphabet, TCargo, TSpec> > TGraph;
    typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
    typedef typename Position<TKeywords>::Type TPos;
    TVertexDescriptor root = addVertex(g);
    assignRoot(g,root);
    resize(terminalStateMap, numVertices(g), Generous());
    assignProperty(terminalStateMap,root,String<TPos>());
    typename Iterator<TKeywords const, Rooted>::Type it = begin(keywords);
    for(;!atEnd(it);goNext(it)) _addStringToTrie(g,terminalStateMap,*it,position(it));
}

//////////////////////////////////////////////////////////////////////////////

/*!
 * @fn Automaton#createTrieOnReverse
 * @brief Creates a trie for all reversed keywords.
 *
 * @signature createTrie(g, terminalStateMap, keywords);
 *
 * @param[in,out] g        Automaton object to create the trie in.
 * @param[in]     terminalStateMap
 *                         An external property map.  The type must be <tt>String&lt;String&lt;unsigned&gt;&gt;</tt>
 *                         since a number of keywords can end in each vertex of the trie.  This is the case in the
 *                         Aho-Corasick algorithm if one pattern is a suffix of another pattern.  Hence, we must
 *                         associate with every vertex a set of indices that correspond to keywords.
 * @param[in]     keywords A set of strings.
 *
 * @see Automaton#createTrie
 */

template<typename TAlphabet, typename TCargo, typename TSpec, typename TTerminalStateMap, typename TKeywords>
inline void
createTrieOnReverse(Graph<Automaton<TAlphabet, TCargo, TSpec> >& g,
                    TTerminalStateMap& terminalStateMap,
                    TKeywords const& keywords)
{
    typedef Graph<Automaton<TAlphabet, TCargo, TSpec> > TGraph;
    typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
    typedef typename Position<TKeywords>::Type TPos;
    TVertexDescriptor root = addVertex(g);
    assignRoot(g,root);
    resize(terminalStateMap, numVertices(g), Generous());
    assignProperty(terminalStateMap,root,String<TPos>());
    typename Iterator<TKeywords const, Rooted>::Type it = begin(keywords);
    for(;!atEnd(it);goNext(it)) {
        typedef typename Value<TKeywords>::Type TKeyword;
        TKeyword tmp;
        typename Iterator<TKeyword const, Rooted>::Type sIt = end(*it);
        while(!atBegin(sIt)) {
            goPrevious(sIt);
            appendValue(tmp,getValue(sIt));
        }
        _addStringToTrie(g,terminalStateMap,tmp,position(it));
    }
}

/*!
 * @fn Automaton#createSuffixTrie
 * @brief Creates a trie of all suffixes of a text.
 *
 * @signature void createSuffixTrie(g, terminalStateMap, text);
 *
 * @param[out] g                The Automaton to create the suffix trie in.
 * @param[out] terminalStateMap An external property map; of type <tt>String&lt;String&lt;unsigned&gt; &gt;</tt>.
 * @param[in]  text             A @link TextConcept @endlink object.
 */

template <typename TAlphabet, typename TCargo, typename TSpec, typename TTerminalStateMap, typename TText>
inline void
createSuffixTrie(Graph<Automaton<TAlphabet, TCargo, TSpec> >& g,
                 TTerminalStateMap& terminalStateMap,
                 TText const& text)
{
    SEQAN_CHECKPOINT
    typedef Graph<Automaton<TAlphabet, TCargo, TSpec> > TGraph;
    typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
    typedef typename Position<TText const>::Type TPosition;
    TVertexDescriptor root = addVertex(g);
    assignRoot(g,root);
    resize(terminalStateMap, numVertices(g), Generous());
    assignProperty(terminalStateMap,root,String<TPosition>());

    for (TPosition i = 0; i < length(text); ++i)
    {
        _addStringToTrie(g,terminalStateMap,suffix(text, i),i);
    }
}

//////////////////////////////////////////////////////////////////////////////


template <typename TAlphabet, typename TCargo, typename TSpec, typename TTerminalStateMap, typename TTexts>
inline void
createSetSuffixTrie(Graph<Automaton<TAlphabet, TCargo, TSpec> >& g,
                    TTerminalStateMap& terminalStateMap,
                    TTexts const& texts)
{
    SEQAN_CHECKPOINT
    typedef Graph<Automaton<TAlphabet, TCargo, TSpec> > TGraph;
    typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
    typedef typename Position<TTexts const>::Type TTextsPosition;
    typedef typename Value<TTexts const>::Type TText;
    typedef typename Position<TText const>::Type TPosition;

    TVertexDescriptor root = addVertex(g);
    assignRoot(g,root);
    resize(terminalStateMap, numVertices(g), Generous());
    assignProperty(terminalStateMap,root,String<TPosition>());

    for (TTextsPosition j = 0; j < length(texts); ++j)
    {
        TText const & text = texts[j];
        for (TPosition i = 0; i < length(text); ++i)
        {
            _addStringToTrie(g,terminalStateMap,suffix(text, i),j);
        }
    }
}

}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
