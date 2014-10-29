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
// Author: Tobias Rausch <rausch@embl.de>
// Author: Anne-Katrin Emde <anne-katrin.emde@fu-berlin.de>
// ==========================================================================
// Iterative instead of recursive implementation of segment match refinement.
// This is required for some large inputs to circumvent stack overflows.
// ==========================================================================

#ifndef SEQAN_CORE_INCLUDE_SEQAN_GRAPH_ALGORITHM_REFINE_EXACT_ITERATIVE_H_
#define SEQAN_CORE_INCLUDE_SEQAN_GRAPH_ALGORITHM_REFINE_EXACT_ITERATIVE_H_

namespace seqan {
	
struct TagExactRefinement_;
typedef Tag<TagExactRefinement_> const ExactRefinement;

// exact method, every cut is made (unless it already exists)
template<typename TValue, typename TValue2, typename TSize>
inline bool
_cutIsValid(String<std::set<TValue> > & all_nodes,
		TValue2 seq_i_pos,
		TSize,
		typename std::set<TValue>::iterator iter,
		TSize,
		Tag<TagExactRefinement_> const)
{
SEQAN_CHECKPOINT
	//cut already exists
	if(iter != all_nodes[seq_i_pos].end())
		return false;
	return true;
}


// necessary for reversed fragments: projected position pos_j is shifted one to the left --> ++pos_j if fragement reversed
template<typename TSize, typename TSpec,typename TPos>
inline void
_updateCutPosition(Fragment<TSize, ExactReversableFragment<TSpec> > & f, TPos & pos_j)
{
	if(f.reversed)
		++pos_j;
}
//template<typename TSize, typename TSpec,typename TPos>
//inline void
//_updateCutPosition(Fragment<TSize, ExactReversableFragment<TSpec> > const& f, TPos & pos_j)
//{
//	if(f.reversed)
//		++pos_j;
//}

// for all other fragment types --> no shifting necessary 
template<typename TFrag,typename TPos>
inline void
_updateCutPosition(TFrag &, TPos &)
{
	return;
}

///////////////////////////////////////////////////////////////////////////////////////////////////////	
//Recursive Refinement
//refine position node_i on sequence seq_i
template<typename TValue, typename TAlignmentString, typename TStringSet,typename TGraph, typename TPropertyMap,typename TSeqMap, typename TTagSpec>
inline void
_refine(TValue node_i, 
	 TValue seq_i_id, 
	 TStringSet & seqs,
	 TSeqMap & seq_map,
	 TAlignmentString & alis, 
	 String<TGraph> & gs, 
	 String<TPropertyMap> & pms, 
     String<std::set<TValue> > & all_nodes, 
	 TValue min_len,
	 Tag<TTagSpec> tag)
{
SEQAN_CHECKPOINT
	typedef typename Cargo<typename Value<TPropertyMap>::Type>::Type TAlignmentPointer;
	typedef typename Iterator<String<TAlignmentPointer>, Rooted>::Type TSegmentIterator;
	//find all segment matches that contain the current position (node_i)
	String<TAlignmentPointer> relevant_segments;
	TValue seq_i_pos = idToPosition(seqs,seq_i_id);
	findIntervalsExcludeTouching(gs[seq_i_pos],pms[seq_i_pos],node_i,relevant_segments);

	
	TSegmentIterator segment_it = begin(relevant_segments);
	TSegmentIterator segment_end = end(relevant_segments);
	//foreach of those segments
	while(segment_it != segment_end)
	{
		TValue match_id = (*segment_it).i1; // segment match
		TValue seg_num = (*segment_it).i2; // first or second segment in segment match?

		//get the sequence that node_i needs to be projected onto (seq_j)
		//and get the projected position (pos_j)
		TValue seq_j_id, node_j;
		_getOtherSequenceAndProject(alis[match_id],seg_num,seq_map,seq_i_id,node_i,seq_j_id,node_j);
		TValue seq_j_pos = idToPosition(seqs,seq_j_id);
		_updateCutPosition(alis[match_id],node_j);

		typename std::set<TValue>::iterator iter;
		iter = all_nodes[seq_j_pos].find(node_j);
		
		//if node does not exist yet ---> insert and continue cutting
		if(_cutIsValid(all_nodes,seq_j_pos,node_j,iter,min_len,tag))
		{
			all_nodes[seq_j_pos].insert(node_j);
			_refine(node_j,seq_j_id,seqs,seq_map,alis,gs,pms,all_nodes,min_len,tag);
		}
		// TODO: else //verschmelzen, abschneiden und ergehen, erst sp�ter... 	
		// do nothing or resolve problems  
	
		++segment_it;
	}
}

//template<typename TFragSize, typename TFragSpec>
//void
//printMatch(Fragment<TFragSize,TFragSpec> & f)
//{
//	::std::cout << "FRAGMENT:" << " f.len = "<< f.len <<std::endl;
//	::std::cout << "f.seqId1 = "<< f.seqId1 << " f.begin1 = " << f.begin1 << std::endl;
//	::std::cout << "f.seqId2 = "<< f.seqId2 << " f.begin2 = " << f.begin2 << std::endl;
//}

//template<typename TAlign>
//void
//printMatch(TAlign & f)
//{
//	::std::cout << f;
//}

///////////////////////////////////////////////////////////////////////////////////////////////////////	
// Construct interval trees 
////////////////////////////////////////////////////////////////////////////////////////////////////


//construct intervals from allignments for each sequence (other Alignment types)
template<typename TInterval, typename TStringSet, typename TAlignmentString, typename TSeqMap>
void
_buildIntervalsForAllSequences(TAlignmentString & alis, 
							   String<String<TInterval> > & intervals, 
	   						   TStringSet & seqs,
							   TSeqMap & seq_map)
{
SEQAN_CHECKPOINT
	
	typedef typename Value<TInterval>::Type TValue;
	typedef typename Cargo<TInterval>::Type TCargo;
	typedef typename Iterator<TAlignmentString,Standard>::Type TAliIterator;
	TAliIterator ali_it = begin(alis,Standard());
	TAliIterator ali_end = end(alis,Standard());
	TValue ali_counter = 0;
	//foreach alignment
	while(ali_it != ali_end)
	{
		TValue seq_i_id,begin_,end_;
		//printMatch(*ali_it);
		//get the first sequence (and its begin and end) that takes part in the alignment (seq_i)
		_getSeqBeginAndEnd(*ali_it,seq_map,seq_i_id,begin_,end_,0);
		TValue seq_i_pos = idToPosition(seqs, seq_i_id);
		//and append the interval (ali_begin, ali_end) with cargo ali* to the list of intervals of seq_i
		appendValue(intervals[seq_i_pos],IntervalAndCargo<TValue,TCargo>(begin_,end_,TCargo(ali_counter,0))); 
	
		//get the second sequence (and its begin and end) that takes part in the alignment (seq_i)
		_getSeqBeginAndEnd(*ali_it,seq_map,seq_i_id,begin_,end_,1);
		seq_i_pos = idToPosition(seqs, seq_i_id);
		//and again append the interval (ali_begin, ali_end) with cargo ali* to the list of intervals of seq_i
		appendValue(intervals[seq_i_pos],IntervalAndCargo<TValue,TCargo>(begin_,end_,TCargo(ali_counter,1))); 
	
		++ali_counter;
		++ali_it;
	}

}


//get all intervals from the alignments and construct an interval tree for each sequence
template<typename TGraph, typename TPropertyMap, typename TAlignmentString, typename TSequence, typename TSetSpec, typename TValue, typename TSeqMap>
void
_createTreesForAllSequences(String<TGraph> & gs, 
						   String<TPropertyMap> & pms, 
						   TAlignmentString & alis, 
						   StringSet<TSequence,TSetSpec> & seqs,
                           TSeqMap & seq_map,
						   TValue numSequences)
{
	//typedef typename Value<TAlignmentString>::Type TAlignment;
//	typedef TAlignment* TCargo;
	typedef Pair<unsigned,unsigned,BitPacked<31,1> > TCargo;
	typedef IntervalAndCargo<int,TCargo> TInterval;
	//typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;

	//std::cout <<"create interval trees...";
	// clock_t start, finish1;
	// double duration;
	// start = clock();
	//one tree for each sequence
	resize(gs,numSequences);
	resize(pms,numSequences);
	
	// and one string of intervals for each sequence
	String<String<TInterval> > intervals;
	resize(intervals,numSequences);
	// fill intervals
	_buildIntervalsForAllSequences(alis,intervals,seqs,seq_map);
	
	TValue i = 0;
	
	while(i < numSequences)
	{
		//std::cout << (numSequences-i) <<" more ("<<length(intervals[i])<<" intervals)... "<<std::flush;
		TValue center = length(seqs[i])/2; // center raus, hat hier nix zu suchen
		//create interval tree!
		createIntervalTree(gs[i], pms[i], intervals[i], center);
		
		//intervals for sequence i are not needed anymore
		clear(intervals[i]);
		++i;
	}
	// finish1 = clock();
	// duration = (double)(finish1 - start) / CLOCKS_PER_SEC;
	// std::cout << "\ntook " << duration << " seconds.\n";
}


///////////////////////////////////////////////////////////////////////////////////////////////////////	
// Construct refined alignment graph 
////////////////////////////////////////////////////////////////////////////////////////////////////

// step 1 of constructing the refined alignment graph: create all nodes
template<typename TStringSet,typename TValue,typename TAliGraph>
void
_makeRefinedGraphNodes(String<std::set<TValue> > & all_nodes,
					  TStringSet & seqs,
					  TAliGraph & ali_g)
{
SEQAN_CHECKPOINT
	typedef typename std::set<TValue>::iterator TSetIterator;
	//for each sequence look at all cut positions and create nodes between them
	for(unsigned int seq_i_pos = 0; seq_i_pos < length(seqs); ++seq_i_pos)
	{
		TValue seq_i_id = positionToId(stringSet(ali_g), seq_i_pos);
		TSetIterator it = all_nodes[seq_i_pos].begin();
		TSetIterator end_it = all_nodes[seq_i_pos].end();
		TSetIterator next_it = it;
		if(next_it != end_it)
			++next_it;
		else
			addVertex(ali_g, seq_i_id, 0, length(seqs[seq_i_pos]));
		
		//first unaligned node
		if(it != end_it && *it != 0)
			addVertex(ali_g, seq_i_id, 0, *it);
		//a new node for each interval
		while(next_it != end_it)
		{
			TValue pos_i = *it;
			addVertex(ali_g, seq_i_id, pos_i, *next_it - pos_i); 
			++it;
			++next_it;
		}
		//last unaligned node
		if(it !=end_it && *it<length(seqs[seq_i_pos]))
			addVertex(ali_g, seq_i_id, *it, (length(seqs[seq_i_pos])) - *it);
		all_nodes[seq_i_pos].clear();
	}
}


// step 2 of constructing the refined alignment graph: add all edges    
// version for exact refinement
template<typename TAlignmentString,typename TStringSet,typename TSeqMap, typename TPropertyMap,typename TScore,typename TAliGraph > 
void
_makeRefinedGraphEdges(TAlignmentString & alis,
					   TPropertyMap & pm,
					  TStringSet & seqs,
				      TSeqMap & seq_map,
				      TScore & score_type,
					  TAliGraph & ali_g,
					  Tag<TagExactRefinement_> const)
{
SEQAN_CHECKPOINT
	typedef typename Value<TAlignmentString>::Type TAlign;
	typedef typename Size<TAlign>::Type TValue;
	typedef typename Iterator<TAlignmentString, Rooted>::Type TAliIterator;
	typedef typename VertexDescriptor<TAliGraph>::Type TVertexDescriptor;
	typedef typename EdgeDescriptor<TAliGraph>::Type TEdgeDescriptor;
	typedef typename Cargo<TAliGraph>::Type TCargo;
	//make edges
	TAliIterator ali_it = begin(alis);
	TAliIterator ali_end = end(alis);
	//for each segment/fragment/alignment
	while(ali_it != ali_end)
	{
		//get sequence, begin position and end position
		TValue seq_id,begin_pos,end_pos;
		_getSeqBeginAndEnd(*ali_it,seq_map,seq_id,begin_pos,end_pos,(TValue)0);
        SEQAN_ASSERT_LEQ(end_pos, length(seqs[idToPosition(seqs, seq_id)]));
        SEQAN_ASSERT(ali_it.data_container == ali_end.data_container);
        SEQAN_ASSERT(ali_it.data_iterator != ali_end.data_iterator);
		
		//get the node represents the current interval (begin_pos until next_cut_pos or end_pos)
		TVertexDescriptor act_knot = findVertex(ali_g,seq_id,begin_pos);
		TValue act_pos = begin_pos;
		TValue seq_j_id_temp,pos_j_begin;
		_getOtherSequenceAndProject(*ali_it,(TValue)0,seq_map,seq_id,act_pos,seq_j_id_temp,pos_j_begin);
	
		//for each interval that lies within the current segment/fragement/alignment
		while(act_pos < end_pos)
		{
			//get other sequence and projected position
			TValue seq_j_id,pos_j;
			_getOtherSequenceAndProject(*ali_it,(TValue)0,seq_map,seq_id,act_pos,seq_j_id,pos_j);
			SEQAN_ASSERT_NEQ(pos_j, static_cast<TValue>(-1));
			//find node that contains the projected position (pos_j)
			TVertexDescriptor vd = findVertex(ali_g, seq_j_id, pos_j);
			bool doAddEdge = true;

//			if(doAddEdge && fragmentBegin(ali_g,vd)!=pos_j) // check if edge makes sense
			if (vd == getNil<TVertexDescriptor>())
				doAddEdge = false;
			else
			{ 
				TValue temp_seq_i_id,temp_act_pos;
				_getOtherSequenceAndProject(*ali_it,(TValue)1,seq_map,seq_j_id,static_cast<TValue>(fragmentBegin(ali_g,vd)),temp_seq_i_id,temp_act_pos);
				if(temp_act_pos == static_cast<TValue>(-1))
					doAddEdge = false;
				else 
				{
					TVertexDescriptor temp_act_knot = findVertex(ali_g, temp_seq_i_id, temp_act_pos);
					if(act_knot!=temp_act_knot)
						doAddEdge = false;
				}
			}
			if(doAddEdge)
			{
				typename Value<TScore>::Type score = _getRefinedMatchScore(score_type,seqs,*ali_it,act_pos,pos_j,fragmentLength(ali_g,act_knot),fragmentLength(ali_g,vd));//,fragmentLength(ali_g,vd));
		//		typename Value<TScore>::Type score = fragmentLength(ali_g,vd);
				score *= _getRefinedAnnoScore(ali_g,pm,vd,act_knot,score_type);
			//this needs to be generalized (makes sense for positive scores only)
				if(score <= 0) score = 1;
				if(score > 0)
				{
					if (findEdge(ali_g, act_knot, vd) == 0) {
						//if(abs((double)fragmentLength(ali_g, act_knot) - (double)fragmentLength(ali_g, vd)) > 20) {
						//	std::cerr << "added edge: " << fragmentLength(ali_g, act_knot) << "  " <<  fragmentLength(ali_g, vd) << std::endl;
						//	std::cerr << *ali_it;
						//	std::cerr << "act_pos=" << act_pos-begin_pos << " pos_j=" << pos_j-pos_j_begin << std::endl; 
						//} else {
							addEdge(ali_g,act_knot,vd,(TCargo)score);
						//}
					}
					else {
						TEdgeDescriptor ed = findEdge(ali_g, act_knot, vd);
						//if((TCargo)score > getCargo(ed))
							//assignCargo(ed, score);
						assignCargo(ed, getCargo(ed)+score);
					}
				}
			}
			//prepare for next interval
			act_pos += fragmentLength(ali_g,act_knot);
			act_knot = findVertex(ali_g,seq_id,act_pos);
		}
		++ali_it;
	}
}





//build refined alignment graph, nodes are numbered ascendingly:
//seq1   0  1  2  3  4 
//seq2   5  6  7  8  9 10
//seq3  11 12 13 14 15 
template<typename TValue,typename TAlignmentString,typename TScore,typename TSequence, typename TSetSpec,typename TAliGraph,typename TSeqMap,typename TTagSpec>
void
_makeAlignmentGraphFromRefinedSegments(String<std::set<TValue> > & all_nodes,
				   TAlignmentString & alis,
				   TScore & score_type,
				   StringSet<TSequence, TSetSpec> & seqs,
				   TSeqMap & seq_map,
				   TAliGraph & ali_g,
			   	   Tag<TTagSpec> const tag, 
				   bool)
{
SEQAN_CHECKPOINT
	//std::cout << "making refined alignment graph...";
	//clock_t start, finish1;
	//double duration;
	//start = clock();
	
	//make nodes (same function for inexact and exact refinement)
	_makeRefinedGraphNodes(all_nodes,seqs,ali_g);

	bool pm = false;
	//add edges (different functions depending on exact/inexact refinement)
	_makeRefinedGraphEdges(alis,pm,seqs,seq_map,score_type,ali_g,tag);
	
	//std::cout << "check\n";
	//finish1 = clock();
	//duration = (double)(finish1 - start) / CLOCKS_PER_SEC;
	//std::cout << "\ntook " << duration << " seconds.\n";
}


//build refined alignment graph as above, but with additional annotation information
template<typename TValue,typename TAlignmentString,typename TScore,typename TSequence, typename TSetSpec,typename TAliGraph,typename TSeqMap,typename TAnnoString,typename TTagSpec>
void
_makeAlignmentGraphFromRefinedSegments(String<std::set<TValue> > & all_nodes,
				   TAlignmentString & alis,
				   TScore & score_type,
				   StringSet<TSequence, TSetSpec> & seqs,
				   TSeqMap & seq_map,
				   TAliGraph & ali_g,
			   	   Tag<TTagSpec> const tag,
				   TAnnoString & annotation)
{
SEQAN_CHECKPOINT
	//std::cout << "making refined alignment graph...";
	//clock_t start, finish1;
	//double duration;
	//start = clock();
	
	//make nodes (same function for inexact and exact refinement)
	_makeRefinedGraphNodes(all_nodes,seqs,ali_g);

	//add annotation to nodes
	//typedef typename Value<TAnnoString>::Type TAnnotation;
	//typedef typename Value<TAnnotation>::Type TLabel;
	typedef char TLabel;
	String<String<TLabel> > pm;
	_addNodeAnnotation(seqs,seq_map,annotation,pm,ali_g,tag);

	//add edges (different functions depending on exact/inexact refinement)
	_makeRefinedGraphEdges(alis,pm,seqs,seq_map,score_type,ali_g,tag);
	
	//std::cout << "check\n";
	//finish1 = clock();
	//duration = (double)(finish1 - start) / CLOCKS_PER_SEC;
	//std::cout << "\ntook " << duration << " seconds.\n";
}





////////////////////////////////////////////////////////////////////////////////////////
//The big matchRefinement function that does everything: build interval trees, do the 
//refinement and construct a refined alignment graph
////////////////////////////////////////////////////////////////////////////////////////
template<typename TAlignmentString, typename TAnnotation, typename TOutGraph, typename TSequence, typename TSetSpec, typename TScore,typename TTagSpec>
void
matchRefinement(TAlignmentString & alis,
				StringSet<TSequence, TSetSpec> & seq, 
				TScore & score_type,
				TOutGraph & ali_graph,
				typename Size<typename Value<TAlignmentString>::Type>::Type min_fragment_len,
				TAnnotation & annotation,
				Tag<TTagSpec> const tag)
{
SEQAN_CHECKPOINT
	////////////////////////////////////////////////////////////////
	//typedefs
	typedef typename Value<TAlignmentString>::Type TAlign;
	typedef typename Iterator<TAlignmentString, Rooted>::Type TAliIterator;
	typedef typename Size<TAlign>::Type TValue;
//	typedef TValue TCargo;
//	typedef Pair<unsigned,unsigned,BitPacked<31,1> > TCargo;
	typedef Pair<unsigned,unsigned,BitPacked<31,1> > TCargo;
	typedef IntervalAndCargo<int,TCargo> TInterval;
	typedef Graph<Directed<void,WithoutEdgeId> > TGraph;
	typedef IntervalTreeNode<TInterval> TNode;
	typedef String<TNode> TPropertyMap;
	typedef typename std::set<TValue>::iterator TSetIterator;
	typedef typename Cargo<typename Value<TPropertyMap>::Type>::Type TAlignmentPointer;
	typedef typename Iterator<String<TAlignmentPointer>, Rooted>::Type TSegmentIterator;


	////////////////////////////////////////////////////////////////
	TValue numSequences = length(seq);
	//weird ID --> good ID map
	std::map<const void * ,int> seq_map;
	for(int i = 0; i < (int) numSequences; ++i)
		seq_map[getObjectId(seq[i])] = i;
	////////////////////////////////////////////////////////////////
	//build interval trees
	String<TGraph> gs;
	String<TPropertyMap> pms;
	_createTreesForAllSequences(gs, pms, alis, seq, seq_map, numSequences);

	////////////////////////////////////////////////////////////////
	//do refinement
	//std::cout <<"refining..."<<std::flush;
	// clock_t start, finish1;
	// double duration;
	// start = clock();
	
	//all_nodes = set of all cut positions
	String<std::set<TValue> > all_nodes;
	resize(all_nodes,numSequences);

	//all_nodes that need to be processed set of all cut positions
	String<std::set<TValue> > all_node_queues;
	resize(all_node_queues,numSequences);

	//call function _refine for each startknoten
	TAliIterator ali_it = begin(alis);
	TAliIterator ali_end = end(alis);
	//for each segment/fragement/alignment
	while(ali_it != ali_end)
	{
		//for each of the two sequences
		for(TValue i = 0; i < 2; ++i)
		{
			TValue seq_i_id,begin_i,end_i;
			_getSeqBeginAndEnd(*ali_it,seq_map,seq_i_id,begin_i,end_i,i);
			TValue seq_i_pos = idToPosition(seq,seq_i_id);
			
			all_node_queues[seq_i_pos].insert(begin_i);
			all_node_queues[seq_i_pos].insert(end_i);
		}	
		++ali_it;
	}


	TSetIterator queueIt;
	bool done = false;
	while(!done)
	{
		for(unsigned seq_i_pos = 0; seq_i_pos < numSequences; ++seq_i_pos)
		{
			queueIt = all_node_queues[seq_i_pos].begin();
			while (queueIt != all_node_queues[seq_i_pos].end())
			{
				TValue node_i = *queueIt;
				TSetIterator iter = all_nodes[seq_i_pos].find(node_i);		
		//		TSetIterator qiter = all_node_queues[seq_i_pos].find(node_i);		
                if(_cutIsValid(all_nodes,seq_i_pos,node_i,iter,min_fragment_len,tag))
                   //&& _cutIsValid(all_node_queues,seq_i_pos,node_i,qiter,min_fragment_len,tag))
//				if(iter == all_nodes[seq_i_pos].end())
				{
					TValue seq_i_id = positionToId(seq, seq_i_pos);
					all_nodes[seq_i_pos].insert(node_i);
					String<TAlignmentPointer> relevant_segments;
					findIntervalsExcludeTouching(gs[seq_i_pos],pms[seq_i_pos],node_i,relevant_segments);
					
					TSegmentIterator segment_it = begin(relevant_segments);
					TSegmentIterator segment_end = end(relevant_segments);
					//foreach of those segments
					while(segment_it != segment_end)
					{
						TValue match_id = (*segment_it).i1;
						TValue seg_num = (*segment_it).i2;						//get the sequence that node_i needs to be projected onto (seq_j)
						//and get the projected position (pos_j)
						TValue seq_j_id, node_j;
						_getOtherSequenceAndProject(alis[match_id],seg_num,seq_map,seq_i_id,node_i,seq_j_id,node_j);
						TValue seq_j_pos = idToPosition(seq,seq_j_id);
						_updateCutPosition(alis[match_id],node_j);

						typename std::set<TValue>::iterator iter_j, qiter_j;
						iter_j = all_nodes[seq_j_pos].find(node_j);
						qiter_j = all_node_queues[seq_j_pos].find(node_j);
						
						//if node does not exist yet ---> insert and continue cutting
                        if(_cutIsValid(all_nodes,seq_j_pos,node_j,iter_j,min_fragment_len,tag)
                            && _cutIsValid(all_node_queues,seq_j_pos,node_j,qiter_j,min_fragment_len,tag))
                        //if(iter_j == all_nodes[seq_j_pos].end())   
						{
							all_node_queues[seq_j_pos].insert(node_j);
						}
					
						++segment_it;
					}
						
				}
				++queueIt;
			}
			all_node_queues[seq_i_pos].clear();
		}
		unsigned i;
		for(i = 0; i < numSequences; ++i)
		{
			queueIt = all_node_queues[i].begin();
			if (queueIt != all_node_queues[i].end())
				break;
		}
		if(i==numSequences)
			done=true;
	}
	_addAnnotationCuts(all_nodes,alis,gs,pms,seq,seq_map,annotation,min_fragment_len,tag);

	// finish1 = clock();
	// duration = (double)(finish1 - start) / CLOCKS_PER_SEC;
	//std::cout << "\ntook " << duration << " seconds.\n";
	//for(int seq_i = 0; seq_i < length(seq); ++seq_i)
	//{
	//	typename std::set<TValue>::iterator it = all_nodes[seq_i].begin();
	//	typename std::set<TValue>::iterator end_it = all_nodes[seq_i].end();
	//
	//	while(it != end_it)
	//	{
	//		std::cout << *it << ",";
	//		++it;
	//	}
	//	std::cout << "\n";
	//}
	//std::cout <<"building tree..."<<std::flush;
	
	////////////////////////////////////////////////////////////////
	//build refined alignment graph
	_makeAlignmentGraphFromRefinedSegments(all_nodes,alis,score_type,seq,seq_map,ali_graph,tag,annotation);
}


///////WRAPPERS

/**
.Function.matchRefinement:
..signature:matchRefinement(matches,stringSet,scoringScheme,refinedGraph)
..param.matches:The set of matches.
...type:Class.Fragment
...type:Class.Align
...type:Spec.Alignment Graph
..param.scoringScheme:The scoring scheme used to score the refined matches (scores are attached to 
edges in the refined Alignment Graph).
...remarks:If no scoring scheme is given, all edges get weight 1.
...type:Class.Score
..include:seqan/graph_align.h
*/
//exact refinement, score type given
template<typename TAlignmentString, typename TScoreValue,typename TScoreSpec,typename TOutGraph, typename TSequence, typename TSetSpec>
void
matchRefinement(TAlignmentString & alis,
				StringSet<TSequence, TSetSpec> & seq, 
				Score<TScoreValue,TScoreSpec> & score_type,
				TOutGraph & ali_graph)
{
SEQAN_CHECKPOINT
	//min_fragment_len = 1   ==> Exact cutting
	bool anno = false;
	matchRefinement(alis,seq,score_type,ali_graph,1,anno,ExactRefinement());
}



/**
.Function.matchRefinement:
..cat:Alignments
..summary:Refines (i.e. cuts into smaller parts) a set of pairwise segment 
matches in such a way that none of the segments partly overlap. They are either 
identical (fully overlapping) or non-overlapping.
..signature:matchRefinement(matches,stringSet,refinedGraph)
..param.stringSet:The StringSet containing the sequences which the matches lie on.
...type:Class.StringSet
..param.refinedGraph:The resulting refined set of matches stored in a graph.
...type:Spec.Alignment Graph
..include:seqan/graph_align.h
*/
//exact refinement, score type not given
template<typename TFragmentString, typename TOutGraph, typename TSequence, typename TSetSpec>
void
matchRefinement(TFragmentString & matches,
				StringSet<TSequence, TSetSpec> & strSet, 
				TOutGraph & ali_graph)
{
	SEQAN_CHECKPOINT
	typename Cargo<TOutGraph>::Type fake_score = 1;
	bool anno = false;
	matchRefinement(matches,strSet,fake_score,ali_graph,1,anno,ExactRefinement());
}

}  // namespace seqan

#endif  // #ifndef SEQAN_CORE_INCLUDE_SEQAN_GRAPH_ALGORITHM_REFINE_EXACT_ITERATIVE_H_
