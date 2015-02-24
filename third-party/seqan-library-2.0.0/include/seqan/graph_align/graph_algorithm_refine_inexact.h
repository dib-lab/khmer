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
// Author: Anne-Katrin Emde <anne-katrin.emde@fu-berlin.de>
// ==========================================================================

#ifndef SEQAN_INCLUDE_SEQAN_GRAPH_ALGORITHM_REFINE_INEXACT_H_
#define SEQAN_INCLUDE_SEQAN_GRAPH_ALGORITHM_REFINE_INEXACT_H_

//SEQAN_NO_DDDOC: do not generate documentation for this file

namespace seqan {


struct TagInexactRefinement_;
typedef Tag<TagInexactRefinement_> const InexactRefinement;


///inexact refinement (cuts that would produce segments shorter than min_len are not made)
template<typename TValue, typename TValue2, typename TSize>
inline bool
_cutIsValid(String<std::set<TValue> > & all_nodes,
        TValue2 seq_i_pos,
        TSize pos_i,
        typename std::set<TValue>::iterator iter,
        TSize min_len,
        Tag<TagInexactRefinement_> const)
{
SEQAN_CHECKPOINT

    //cut already exists
    if(iter != all_nodes[seq_i_pos].end())
        return false;
    typename std::set<TValue>::iterator tmp_iter = all_nodes[seq_i_pos].upper_bound(pos_i);
    if(tmp_iter != all_nodes[seq_i_pos].end())
        if((*tmp_iter - pos_i) < min_len)
            return false;
    if(tmp_iter != all_nodes[seq_i_pos].begin())
    {
        --tmp_iter;
        if((pos_i - *tmp_iter) < min_len)
            return false;
    }
    return true;
}


// returns node begin or end position closest to pos
// i.e. closest refined position
template<typename TAliGraph, typename TVertexDescriptor, typename TId, typename TPosition>
TPosition
_getClosestRefinedNeighbor(TAliGraph & ali_g,
                           TVertexDescriptor & vd,
                           TId /*seq*/,
                           TPosition pos)
{
SEQAN_CHECKPOINT
    if(pos-fragmentBegin(ali_g,vd) < fragmentBegin(ali_g,vd)+fragmentLength(ali_g,vd)-pos)
        return fragmentBegin(ali_g,vd);
    else
        return fragmentBegin(ali_g,vd) + fragmentLength(ali_g,vd);
}


// get closest refined position to end position of fragment (cut_end_pos)
// and corresponding node (end_knot)
template<typename TAliGraph, typename TId, typename TPosition>
void
_getCutEndPos(TAliGraph & ali_g,
              typename VertexDescriptor<TAliGraph>::Type & end_knot,
              TId seq,
              TPosition act_begin_pos,
              TPosition end_pos,
              TPosition & cut_end_pos)
{
SEQAN_CHECKPOINT
    end_knot = findVertex(ali_g,seq,end_pos-1);//end_pos1 is the first position of the next node
    if(end_pos == fragmentBegin(ali_g,end_knot) + fragmentBegin(ali_g,end_knot))
        cut_end_pos = end_pos;
    else
    {
        cut_end_pos = _getClosestRefinedNeighbor(ali_g,end_knot,seq,end_pos);
        if(cut_end_pos <= act_begin_pos) end_knot = getNil<typename VertexDescriptor<TAliGraph>::Type>();
        else
        {
            end_knot =  findVertex(ali_g,seq,cut_end_pos-1);
            SEQAN_ASSERT(cut_end_pos == fragmentBegin(ali_g,end_knot)+fragmentLength(ali_g,end_knot));
        }
    }
}


// get closest refined position to begin position of fragment (cut_act_pos)
// and corresponding node (act_knot)
template<typename TAliGraph, typename TId, typename TPosition>
void
_getCutBeginPos(TAliGraph & ali_g,
              typename VertexDescriptor<TAliGraph>::Type & act_knot,
              TId seq,
              TPosition act_end_pos,
              TPosition act_pos,
              TPosition & cut_act_pos)
{
SEQAN_CHECKPOINT

    act_knot = findVertex(ali_g,seq,act_pos);
    //if completely refined
    if(act_pos == fragmentBegin(ali_g,act_knot))
        cut_act_pos = act_pos;
    else //if incompletely refined
    {
        cut_act_pos = _getClosestRefinedNeighbor(ali_g,act_knot,seq,act_pos);
        if(cut_act_pos > act_end_pos) act_knot = getNil<typename VertexDescriptor<TAliGraph>::Type>();
        else
        {
            act_knot =  findVertex(ali_g,seq,cut_act_pos); // have to watch out with cut_act_pos==seqLength!
            SEQAN_ASSERT(act_knot == getNil<typename VertexDescriptor<TAliGraph>::Type>() ||
                         cut_act_pos == fragmentBegin(ali_g, act_knot));
        }
    }
}



//step 2 of constructing the refined alignment graph: add all edges
//version for inexact refinement
template<typename TAlignmentString,typename TPropertyMap,typename TStringSet,typename TSeqMap, typename TScore,typename TAliGraph>
void
_makeRefinedGraphEdges(TAlignmentString & alis,
                       TPropertyMap & , //pm,
                      TStringSet & seqs,
                      TSeqMap & seq_map,
                      TScore & score_type,
                      TAliGraph & ali_g,
                      Tag<TagInexactRefinement_> const)
{
SEQAN_CHECKPOINT
    typedef typename Value<TAlignmentString>::Type TAlign;
    typedef typename Position<TAlign>::Type TPosition;
    typedef typename Id<TAlign>::Type TId;
    typedef typename Iterator<TAlignmentString, Rooted>::Type TAliIterator;
    typedef typename VertexDescriptor<TAliGraph>::Type TVertexDescriptor;
    typedef typename EdgeDescriptor<TAliGraph>::Type TEdgeDescriptor;
    //typedef typename Cargo<TAliGraph>::Type TCargo;

        TVertexDescriptor nilVertex = getNil<TVertexDescriptor>();
    //make edges
    TAliIterator ali_it = begin(alis);
    TAliIterator ali_end = end(alis);
    //for each segment/fragment/alignment
    while(ali_it != ali_end)
    {
        //get first sequence that takes part in the alignment + boundaries of the ali
        TId seq1;
        TPosition begin_pos1,end_pos1;
        _getSeqBeginAndEnd(*ali_it,seq_map,seq1,begin_pos1,end_pos1,(TId)0);

        //get second sequence that takes part in the alignment + boundaries of the ali
        TId seq2;
        TPosition begin_pos2,end_pos2;
        _getSeqBeginAndEnd(*ali_it,seq_map,seq2,begin_pos2,end_pos2,(TId)1);

        //get the last node that is within the current ali
        TVertexDescriptor end_knot1;
        TPosition cut_end_pos1;
        _getCutEndPos(ali_g,end_knot1,seq1,begin_pos1,end_pos1,cut_end_pos1);
        if(end_knot1 == nilVertex) // there is no node --> fragment disappeared in min_frag_len heuristic
            continue;

        //get the node that represents the current interval (begin_pos until next_cut_pos or end_pos)
        TVertexDescriptor act_knot1;
        TPosition cut_act_pos1,act_pos1;
        act_pos1 = begin_pos1;
        _getCutBeginPos(ali_g,act_knot1,seq1,end_pos1,act_pos1,cut_act_pos1);
        if(act_knot1 == nilVertex) // there is no node, can this happen here?
            continue;
        TPosition act_end_pos1 = cut_act_pos1 + fragmentLength(ali_g,act_knot1);
        //walk through cuts on the first sequence
//        while (act_end_pos1 <= cut_end_pos1)
        while (true)
        {
            //get other sequence and projected position
            //TId seq2;
            TPosition act_pos2;
            _getOtherSequenceAndProject(*ali_it,0,seq_map,seq1,act_pos1,seq2,act_pos2);

            //get node that corresponds to that position
            TVertexDescriptor act_knot2;
            TPosition cut_act_pos2;
            _getCutBeginPos(ali_g,act_knot2,seq2,end_pos2,act_pos2,cut_act_pos2);
            if(act_knot2 == nilVertex) // there is no corresponding node in second sequence
                break;
            //corresponding end on seq2 (there might be more than one node on seq2 that corresponds
            //to the same interval (=node) on seq1)
            TPosition act_end_pos2;
             _getOtherSequenceAndProject(*ali_it,0,seq_map,seq1,_min(act_end_pos1,end_pos1)-1,seq2,act_end_pos2);
            ++act_end_pos2;
            TVertexDescriptor act_end_knot2;
            TPosition cut_act_end_pos2;
            _getCutEndPos(ali_g,act_end_knot2,seq2,begin_pos2,act_end_pos2,cut_act_end_pos2);
            if(act_end_knot2 == nilVertex) // there is no node at all in second sequence
                break;

            if(cut_act_pos2 == cut_act_end_pos2)
                break;
            while(true)
            {
                //should at the moment return score for:
                //
                //seq1 = ....cr...rc....
                //            ||||||
                //seq2 = ...c.r...rc....
                //bzw
                //seq1 = ..cr.....x....   man will aber nur    ..cr......x....
                //          |||||||-                             ---||||||
                //seq2 = ...r.c...rc...                        ...r.c...rc....
                typename Value<TScore>::Type score = 0;
                score = _getRefinedMatchScore(score_type,seqs,*ali_it,act_pos1,act_pos2,act_end_pos1-act_pos1,cut_act_end_pos2);
                //score *= _getRefinedAnnoScore(ali_g,pm,act_knot1,act_knot2,score_type);
                //add score for
                //
                //seq1 = ...-cr....x....
                //          ||
                //seq2 = ...c.r...rc....
//                    score += getLeftRestScore(score_type,seqs,seq1,seq2,act_pos1,cut_act_pos1,act_pos2,cut_act_pos2);
                if(score > 0)
                {    if(findEdge(ali_g,act_knot1,act_knot2)==0)
                        addEdge(ali_g,act_knot1,act_knot2,score);
                    else
                    {
                        TEdgeDescriptor ed = findEdge(ali_g, act_knot1, act_knot2);
                        //if((TCargo)score > getCargo(ed))
                            //assignCargo(ed, score);
                        assignCargo(ed, getCargo(ed)+score);
                    }
                }
                if(act_knot2==act_end_knot2)
                    break;
                act_pos2 = cut_act_pos2 + fragmentLength(ali_g,act_knot2);
                _getCutBeginPos(ali_g,act_knot2,seq2,end_pos2,act_pos2,cut_act_pos2);
            }
            if(act_knot1 == end_knot1)
                break;
            act_pos1 = act_end_pos1;
            act_knot1 = findVertex(ali_g,seq1,act_pos1);
            cut_act_pos1 = act_pos1;
            act_end_pos1 = cut_act_pos1 + fragmentLength(ali_g,act_knot1);
        }
        ++ali_it;
    }
}

// TODO(holtgrew): Documentation is incomplete.

/*!
 * @fn matchRefinement
 * @headerfile <seqan/graph_align.h>
 * @brief Refines (i.e. cuts into smaller parts) a set of pairwise segment matches in such a way that none of the
 *        segments partly overlap. They are either identical (fully overlapping) or non-overlapping.
 *
 * @signature void matchRefinement(matches, stringSet[, scoringScheme], refinedGraph);
 *
 * @param[out] matches       The set of matches. Types: Fragment, Align, Alignment Graph
 * @param[out] refinedGraph  The resulting refined set of matches stored in a graph.  Types: Alignment Graph
 * @param[out] stringSet     The StringSet containing the sequences which the matches lie on. Types: StringSet
 * @param[in]  scoringScheme The scoring scheme used to score the refined matches (scores are attached to edges
 *                           in the refined AlignmentGraph).  If no scoring scheme is given, all edges get weight 1.
 *                           Types: Score
 */

//score type given, min fragment length given, if > 1 ==> inexact refinement
template<typename TAlignmentString, typename TScoreValue,typename TScoreSpec, typename TOutGraph, typename TSequence, typename TSetSpec>
void
matchRefinement(TAlignmentString & alis,
                StringSet<TSequence, TSetSpec> & seq,
                Score<TScoreValue,TScoreSpec> & score_type,
                TOutGraph & ali_graph,
                unsigned int min_frag_len)
{
SEQAN_CHECKPOINT
    bool anno = false;
    if(min_frag_len > 1)
        matchRefinement(alis,seq,score_type,ali_graph,min_frag_len,anno,InexactRefinement());
    else
        matchRefinement(alis,seq,score_type,ali_graph,min_frag_len,anno,ExactRefinement());
}


//score type not given, min fragment length given, if > 1 ==> inexact refinement
template<typename TAlignmentString, typename TOutGraph, typename TSequence, typename TSetSpec>
void
matchRefinement(TAlignmentString & alis,
                StringSet<TSequence, TSetSpec> & seq,
                TOutGraph & ali_graph,
                unsigned int min_frag_len)
{
SEQAN_CHECKPOINT
//    Score<int,FakeScore > fake_score;
    typename Cargo<TOutGraph>::Type fake_score = 1;
    bool anno = false;
    if(min_frag_len > 1)
        matchRefinement(alis,seq,fake_score,ali_graph,min_frag_len,anno,InexactRefinement());
    else
        matchRefinement(alis,seq,fake_score,ali_graph,min_frag_len,anno,ExactRefinement());
}

}  // namespace seqan

#endif  // #ifndef SEQAN_INCLUDE_SEQAN_GRAPH_ALGORITHM_REFINE_INEXACT_H_
