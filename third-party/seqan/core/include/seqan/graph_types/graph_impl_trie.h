// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2013, Knut Reinert, FU Berlin
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

#ifndef SEQAN_HEADER_GRAPH_IMPL_TRIE_H
#define SEQAN_HEADER_GRAPH_IMPL_TRIE_H

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////
// Graph - Trie
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

/**
.Spec.Trie:
..cat:Graph
..general:Class.Graph
..summary:A keyword trie.
..description:
...image:trieGraph|A trie for the words announce, annual, and annually.
..remarks:A keyword trie is a special automaton and thus, it is not implemented in its own class.
It solely provides create functions where based upon a set of strings a keyword trie is created.
..signature:Graph<Automaton<TAlphabet, TCargo, TSpec> > 
..param.TAlphabet:The alphabet type that is used for the transition labels.
...metafunction:Metafunction.Alphabet
...remarks:Use @Metafunction.Alphabet@ to get the type of the labels in an automaton.
...default:$char$
..param.TCargo:The cargo type that can be attached to the edges.
...metafunction:Metafunction.Cargo
...remarks:Use @Metafunction.Cargo@ to get the cargo type of an undirected graph.
...default:$void$
..param.TSpec:The specializing type for the graph.
...metafunction:Metafunction.Spec
...remarks:Use WithoutEdgeId here to omit edge ids.
Note: If edges do not store ids external property maps do not work.
...default:$Default$, see @Tag.Default@.
..include:seqan/graph_types.h
*/


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

/**
.Function.createTrie
..class:Spec.Trie
..cat:Graph
..summary:Creates a trie.
..signature:createTrie(g, terminalStateMap, keywords)
..param.g:Out-parameter: An automaton.
...type:Spec.Trie
..param.terminalStateMap:Out-parameter: An external property map.
...type:Class.External Property Map
...remarks:The external property map must be a String<String<unsigned int> > because
in every vertex of the trie a number of keywords can end. This is the case in the Aho-Corasick
algorithm if one pattern is a suffix of another pattern! Hence, we must associate with every vertex a set of indices that correspond to keywords.
..param.keywords:In-parameter: A set of strings.
...type:Class.String
..returns:void
..see:Function.createTrieOnReverse
..include:seqan/graph_types.h
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

/**
.Function.createTrieOnReverse
..class:Spec.Trie
..cat:Graph
..summary:Creates a trie for all reversed keywords.
..signature:createTrieOnReverse(g, terminalStateMap, keywords)
..returns.param.g:Out-parameter: An automaton.
...type:Spec.Trie
..returns.param.terminalStateMap:Out-parameter: An external property map.
...type:Class.External Property Map
...remarks:The external property map must be a String<String<unsigned int> > because
in every vertex of the trie a number of keywords can end. This is the case in the Aho-Corasick
algorithm if one pattern is a suffix of another pattern! Hence, we must associate with every vertex a set of indices that correspond to keywords.
..param.keywords:In-parameter: A set of strings.
...type:Class.String
..returns:void
..see:Function.createTrie
..include:seqan/graph_types.h
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



/**
.Function.createSuffixTrie
..class:Spec.Trie
..cat:Graph
..summary:Creates a trie of all suffixes of a text.
..signature:createSuffixTrie(g, terminalStateMap, text)
..param.g:Out-parameter: An automaton.
...type:Spec.Trie
..param.terminalStateMap:Out-parameter: An external property map.
...type:Class.External Property Map
...remarks:The external property map must be a String<String<unsigned int> >.
..param.text:In-parameter: A text.
...type:Class.String
..returns:void
..see:Function.createTrie
..see:Function.createTrieOnReverse
..include:seqan/graph_types.h
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
