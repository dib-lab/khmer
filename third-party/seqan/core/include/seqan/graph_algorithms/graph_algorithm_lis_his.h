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

#ifndef SEQAN_HEADER_GRAPH_ALGORITHM_LIS_HIS_H
#define SEQAN_HEADER_GRAPH_ALGORITHM_LIS_HIS_H

namespace SEQAN_NAMESPACE_MAIN
{

struct Lcs_;
typedef Tag<Lcs_> Lcs;

//////////////////////////////////////////////////////////////////////////////
// LIS: Longest Increasing Subsequence
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

template<typename TSortedSequence, typename TKey>
inline typename TSortedSequence::const_iterator
_previousInSortedSequence(TSortedSequence const& list, TKey const key) {
	SEQAN_CHECKPOINT
	typedef typename TSortedSequence::const_iterator TSortedSequenceIter;

	TSortedSequenceIter a_k_it = list.lower_bound(key);
	// Now we need to move one to the front

	if (a_k_it != list.end()) {
		// If we are at the beginning, no predecessor
		if (a_k_it == list.begin()) a_k_it = list.end();
		else --a_k_it;
	} else {
		// If we are at the end, the predecessor is the last element of the list
		TSortedSequenceIter tmp = list.begin();
		if (tmp != list.end()) {
			do {
				a_k_it = tmp;
			} while(++tmp != list.end());
		}
	}
	return a_k_it;
}


//////////////////////////////////////////////////////////////////////////////


template<typename TSortedSequence, typename TIterator>
inline typename TSortedSequence::const_iterator
_nextInSortedSequence(TSortedSequence const& list, TIterator const& prev) {
	SEQAN_CHECKPOINT
	typedef typename TSortedSequence::const_iterator TSortedSequenceIter;
		
	TSortedSequenceIter b_l_it;
	if (prev == list.end()) b_l_it = list.begin();
	else b_l_it = list.upper_bound(*prev);

	return b_l_it;
}


//////////////////////////////////////////////////////////////////////////////

/**
.Function.longestIncreasingSubsequence:
..summary:Computes the longest increasing subsequence.
..cat:Alignments
..signature:longestIncreasingSubsequence(str, pos)
..param.str:In-parameter: An arbitrary string.
...type:Class.String
..param.pos:Out-parameter: A String with the positions that belong to the longest increasing subsequence.
...remarks:
The last position in pos indicates the first element in the longest increasing subsequence.
That's why pos should be a Block-String (Stack).
..include:seqan/graph_algorithms.h
*/
template<typename TString, typename TPositions>
inline void
longestIncreasingSubsequence(TString const& str, TPositions& pos) {
	SEQAN_CHECKPOINT

	// The list of decreasing covers, only the smallest number must be remembered
	// See Gusfield
	typedef std::pair<typename Value<TString>::Type, typename Position<TPositions>::Type> TKey;
	typedef std::set<TKey, std::less<TKey> > TSortedSequence;
	typedef typename TSortedSequence::const_iterator TSortedSequenceIter;
	TSortedSequence list;

	// The trace-back graph
	typedef Graph<Directed<void, WithoutEdgeId> > TGraph;
	typedef VertexDescriptor<TGraph>::Type TVertexDescriptor;
	TGraph g;

	// Walk through the sequence and build the decreasing covers
	typedef typename Iterator<TString const, Rooted>::Type TStringIter;
	TStringIter endIt = end(str);
	for(TStringIter it = begin(str); it != endIt; ++it) {
		// Get previous element
		TSortedSequenceIter a_k_it = _previousInSortedSequence(list, std::make_pair(*it, 0)); 
		
		// Get next element
		TSortedSequenceIter b_l_it = _nextInSortedSequence(list, a_k_it);
		
		// Delete from list
		if (b_l_it != list.end()) list.erase(*b_l_it);

		// Insert new list element
		list.insert(std::make_pair(*it, position(it)));

		// Create the corresponding node
		// Note: The VertexDescriptor == position(it)
		addVertex(g);


		// Connect to predecessor
		if (a_k_it != list.end()) addEdge(g, (TVertexDescriptor) position(it), (TVertexDescriptor) a_k_it->second);
	}

	// Trace-back
	if (list.rbegin() == list.rend()) return;
	else {
		// Start with the maximal position in the list == Vertex Descriptor
		TVertexDescriptor v = list.rbegin()->second;
		while (true) {
			appendValue(pos, v, Generous());
			if (g.data_vertex[v]) v = (*g.data_vertex[v]).data_target;
			else break;
		}
	}
}


//////////////////////////////////////////////////////////////////////////////
// LCS: Longest Common Subsequence
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

/**
.Function.longestCommonSubsequence:
..summary:Computes the longest common subsequence.
..cat:Alignments
..signature:longestCommonSubsequence(str1, str2, nSize, pos)
..param.str1:In-parameter: An arbitrary string.
...type:Class.String
..param.str2:In-parameter: An arbitrary string.
...type:Class.String
..param.pos:Out-parameter: A String with pairs of positions that indicate the longest common subsequence.
...remarks:
..include:seqan/graph_algorithms.h
*/
template<typename TString1, typename TString2, typename TNeighborhoodSize, typename TFinalPos>
inline void
longestCommonSubsequence(TString1 const& str1,
						 TString2 const& str2,
						 TNeighborhoodSize nSize,
						 TFinalPos& pos) 
{
	SEQAN_CHECKPOINT
	typedef typename Value<TString1>::Type TValue;
	typedef typename Size<TString1>::Type TSize;
	typedef typename Position<TString1>::Type TPos;
	TSize alphabet_size = ValueSize<TValue>::VALUE;

	// The occurrences of each letter in the second string
	typedef String<TPos> TPositions;
	String<TPositions> occ;
	resize(occ, alphabet_size, TPositions());
	typedef typename Iterator<TString2 const, Standard>::Type TStringIter;
	TStringIter itStr2 = begin(str2, Standard());
	TStringIter endItStr2 = end(str2, Standard());
	TPos current_pos = 0;
	for(; itStr2 != endItStr2; ++itStr2, ++current_pos) appendValue(occ[ordValue(*itStr2)], current_pos, Generous());

	// Build the combined string
	String<TPos> finalSeq;
	String<TPos> mapping;
	TStringIter itStr1 = begin(str1, Standard());
	TStringIter endItStr1 = end(str1, Standard());
	current_pos = 0;
	TPos diff = 0;
	for(; itStr1 != endItStr1; ++itStr1, ++current_pos) {
		TPositions& current_occ = occ[ordValue(*itStr1)];
		for(int i = length(current_occ)-1; i>=0; --i) {
			// Do we have a neighborhood
			diff = (current_pos < current_occ[i]) ? current_occ[i] - current_pos : current_pos - current_occ[i];
			if (diff > (TPos) nSize) continue;
			appendValue(finalSeq, current_occ[i], Generous());
			appendValue(mapping, current_pos, Generous());
		}
	}

	// Call longest increasing subsequence
	typedef String<TSize> TResult;
	TResult result;
	longestIncreasingSubsequence(finalSeq, result);
	
	// Insert the common pairs
	typedef typename Iterator<TResult, Standard>::Type TResultIter;
	TResultIter itResult = begin(result, Standard());
	TResultIter endResult = end(result, Standard());
	for(; itResult != endResult; ++itResult) 
		appendValue(pos, std::make_pair(mapping[*itResult], finalSeq[*itResult]), Generous());
}

//////////////////////////////////////////////////////////////////////////////

template<typename TAlign, typename TStringSet>
inline int
globalAlignment(TAlign& align,
				TStringSet const& str,
				Lcs)
{
	SEQAN_CHECKPOINT
	typedef typename Id<TStringSet>::Type TId;
	typedef typename Size<TStringSet>::Type TSize;
	TId id1 = positionToId(str, 0);
	TId id2 = positionToId(str, 1);
		
	// Lcs between first and second string
	String<std::pair<TSize, TSize> > pos1;
	longestCommonSubsequence(str[0], str[1], 100, pos1);
	//longestCommonSubsequence(str[0], str[1], 800, pos1);

	// Extend the matches as long as possible
	TSize oldI = 0;
	TSize oldJ = 0;
	TSize totalLen = 0;
	if (length(pos1)) {
		TSize lenMatch = 1;				
		int last = length(pos1)-1;		
		TSize iBegin = pos1[last].first;
		TSize jBegin = pos1[last].second;
		for(int z = last - 1; z>=0; --z) {
			if ((pos1[z].first == pos1[z+1].first + 1) &&
				(pos1[z].second == pos1[z+1].second + 1)) 
			{
				++lenMatch;
			} else {
				if (oldI < iBegin) _alignTracePrint(align, str[0], str[1], id1, oldI, id2, (TSize) 0, (TSize) iBegin - oldI, 1);
				if (oldJ < jBegin) _alignTracePrint(align, str[0], str[1], id1, (TSize) 0, id2, oldJ, (TSize) jBegin - oldJ, 2);
				oldI = iBegin + lenMatch;
				oldJ = jBegin + lenMatch;
			
				_alignTracePrint(align, str[0], str[1], id1, iBegin, id2, jBegin, lenMatch, 0);
				totalLen += lenMatch;
				lenMatch = 1;
				iBegin = pos1[z].first;
				jBegin = pos1[z].second;
			}
		}
		// Process last match
		if (oldI < iBegin) _alignTracePrint(align, str[0], str[1], id1, oldI, id2, (TSize) 0, (TSize) iBegin - oldI, 1);
		if (oldJ < jBegin) _alignTracePrint(align, str[0], str[1], id1, (TSize) 0, id2, oldJ, (TSize) jBegin - oldJ, 2);
		oldI = iBegin + lenMatch;
		oldJ = jBegin + lenMatch;
		_alignTracePrint(align, str[0], str[1], id1, iBegin, id2, jBegin, lenMatch, 0);
		totalLen += lenMatch;
	}
	// Process left overs
	if (oldI < length(str[0])) _alignTracePrint(align, str[0], str[1], id1, oldI, id2, (TSize) 0, (TSize) length(str[0]) - oldI,  1);
	if (oldJ < length(str[1])) _alignTracePrint(align, str[0], str[1], id1, (TSize) 0, id2, oldJ, (TSize) length(str[1]) - oldJ, 2);
	
	return (int) totalLen;
}


//////////////////////////////////////////////////////////////////////////////
// HIS: Heaviest Increasing Subsequence
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

/**
.Function.heaviestIncreasingSubsequence:
..summary:Computes the heaviest increasing subsequence.
..cat:Alignments
..signature:heaviestIncreasingSubsequence(str, weights, pos)
..param.str:In-parameter: An arbitrary string.
...type:Class.String
..param.weights:In-parameter: A weight for each position in the string.
..param.pos:Out-parameter: A String of positions that indicate the members of the heaviest increasing subsequence.
...remarks:
The last position in pos indicates the first member of the heaviest increasing subsequence.
That's why pos should be a Block-String (Stack).
Note that only members that contribute a weight are selected, that is, positions with associated weight=0 are ignored.
..include:seqan/graph_algorithms.h
*/
template<typename TString, typename TWeightMap, typename TPositions>
inline typename Value<TWeightMap>::Type
heaviestIncreasingSubsequence(TString const& str, 
							  TWeightMap const& weights, 
							  TPositions& pos) 
{
	SEQAN_CHECKPOINT
	typedef typename Size<TString>::Type TSize;
	typedef typename Value<TString>::Type TValue;
	typedef typename Value<TPositions>::Type TPos;
	typedef typename Value<TWeightMap>::Type TWeight;

	// The list of decreasing covers, only the smallest element of each member must be remembered
	typedef std::pair<TValue, std::pair<TWeight, TPos> > TKey;
	typedef std::set<TKey, std::less<TKey> > TSortedSequence;
	typedef typename TSortedSequence::const_iterator TSortedSequenceIter;
	TSortedSequence list;
	
	// The trace-back graph
	typedef Graph<Directed<void, WithoutEdgeId> > TGraph;
	typedef VertexDescriptor<TGraph>::Type TVertexDescriptor;
	TGraph g;

	// Walk through the sequence and build the decreasing covers
	typedef typename Iterator<TString const, Standard>::Type TStringIter;
	TStringIter it = begin(str, Standard());
	TStringIter endIt = end(str, Standard());
	TSize pos_of_iterator = 0;
	TWeight w = 0;
	for(; it != endIt; ++it, ++pos_of_iterator) {
		w = weights[pos_of_iterator];
		// Letters that do not contribute a weight (e.g., w = 0) are excluded!
		// Weights must increase!
		if (w == 0) {
			addVertex(g);  // Note: The vertex id corresponds to the position
			continue;
		}


		// Get previous element
		TSortedSequenceIter a_k_it = _previousInSortedSequence(list, std::make_pair(*it, std::make_pair(0, 0))); 
		
		// Get next element
		TSortedSequenceIter b_l_it = _nextInSortedSequence(list, a_k_it);

		// Determine new weight
		if (a_k_it != list.end()) w += a_k_it->second.first;

		// Delete from list
		while ((b_l_it != list.end()) && 
				(w >= b_l_it->second.first)) {
					TSortedSequenceIter tmp = b_l_it;
					b_l_it = _nextInSortedSequence(list, b_l_it);
					list.erase(*tmp);
		}

		// Insert new list element
		if ((b_l_it == list.end()) ||
			(*it < b_l_it->first)) {
				list.insert(std::make_pair(*it, std::make_pair(w, pos_of_iterator)));
		}

		// Create the corresponding node, pos_of_iterator == Vertex Descriptor
		addVertex(g);

		// Connect to predecessor
		if (a_k_it != list.end()) addEdge(g, (TVertexDescriptor) pos_of_iterator, (TVertexDescriptor) a_k_it->second.second);
	}

	// Trace-back
	w = 0;
	if (list.rbegin() == list.rend()) return 0;
	else {
		// Last vertex is end of heaviest increasing subsequence
		TVertexDescriptor v = list.rbegin()->second.second;
		while (true) {
			appendValue(pos, v, Generous());
			w+=weights[v];
			if (g.data_vertex[v]) v = (*g.data_vertex[v]).data_target;
			else break;
		}
	}
	return w;
}



}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
