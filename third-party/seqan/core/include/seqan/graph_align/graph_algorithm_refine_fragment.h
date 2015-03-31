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

#ifndef SEQAN_CORE_INCLUDE_SEQAN_GRAPH_ALGORITHM_REFINE_FRAGMENT_H_
#define SEQAN_CORE_INCLUDE_SEQAN_GRAPH_ALGORITHM_REFINE_FRAGMENT_H_

namespace seqan {

	
///////////////////////////////////////////////////////////////////////////////////////////////////////	
//Functions for Fragments
//project onto other sequence for Graph<Alignment>
template<typename TFragSize, typename TFragSpec,typename TValue,typename TId1, typename TPos1, typename TId2, typename TPos2, typename TMap>
void
_getOtherSequenceAndProject(Fragment<TFragSize,TFragSpec> & segment,
						TValue seg_num,
						   TMap &,
						   TId1 seq_i_id,
						   TPos1 pos_i,
						   TId2 & seq_j_id,
						   TPos2 & pos_j)
{
SEQAN_CHECKPOINT
	getProjectedPosition(segment,seg_num, seq_i_id, pos_i,seq_j_id,pos_j);
	
	//if(seq_i_id == sequenceId(segment,0))
	//	seq_j_id = sequenceId(segment,1);
	//else
	//	seq_j_id = sequenceId(segment,0);
}


//given seq and segment, get the sequenceId (seq_i) and its begin and end
//if seq = 0 get first sequence (that takes part in the segment match)
//if seq = 1 get second sequence
template<typename TFragSize, typename TFragSpec, typename TId, typename TPosition, typename TId2>
void
_getSeqBeginAndEnd(Fragment<TFragSize,TFragSpec> & segment,
				  std::map<const void * ,int> &, 
				  TId & seq_i_id, 
				  TPosition & begin_i, 
				  TPosition & end_i,
				  TId2 seq)
{
SEQAN_CHECKPOINT
	seq_i_id = sequenceId(segment,seq);
	if(seq==0)
		begin_i = segment.begin1; // fragmentBegin(segment,seq_i_id);
	else
		begin_i = segment.begin2; // fragmentBegin(segment,seq_i_id);
	end_i = begin_i + segment.len; //fragmentLength(segment,seq_i_id);
}


////////////////////////////////////////////////////////////////////////////////////////
// 50000 _getRefinedMatchScore Functions
////////////////////////////////////////////////////////////////////////////////////////
//get score for alignment starting at pos_i on one sequence (first sequence if i_am_first==true)
//and pos_j on other sequence (second sequence if i_am_first==true), if len1!=len2 then the refinement
//process was stopped (the cut is not exact)
//template<typename TScore,typename TStringSet,typename TFragId,typename TFragPos,typename TFragSize, typename TFragSpec,typename TValue>
//typename Value<TScore>::Type
//_getRefinedMatchScore(TScore & score_type, 
//		 TStringSet & seqs,
//		 Fragment<TFragId,TFragPos,TFragSize,TFragSpec> & segment, 
//		 TValue pos_i, 
//		 TValue pos_j,
//		 TValue len1, 
//		 TValue len2)
//{
//SEQAN_CHECKPOINT
//	typename Infix<typename Value<TStringSet>::Type>::Type label0 = label(segment,seqs,sequenceId(segment,0));
//	typename Infix<typename Value<TStringSet>::Type>::Type label1 = label(segment,seqs,sequenceId(segment,1));
//	int i = 0;
//	typename Value<TScore>::Type ret_score = 0;
//	TValue len = (len1 < len2) ? len1 : len2;
//	while(i < len)
//	{
//		ret_score += score(score_type,label0[i],label1[i]);
//		++i;
//	}
//	len = (len1 > len2) ? len1 : len2;
//	ret_score += (len - i) * scoreGapExtend(score_type);
//	return ret_score;
//}				
//
//
///////////////////////////////////////////////////////////////////////////////////////////////
////die n�chsten beiden funktionen: f�r Fragmente und Score vom Typ Simple
////f�r den fall dass es keine mismatches innerhalb der segmente gibt und Score vom typ Simple ist
////TODO: m�sste f�r einen bestimmten TFragSpec sein (Exact oder noMismatches)
////get score for alignment starting at pos_i on one sequence (first sequence if i_am_first==true)
////and pos_j on other sequence (second sequence if i_am_first==true), if len1!=len2 then the refinement
////process was stopped (the cut is not exact)
//template<typename TScoreValue,typename TStringSet,typename TFragId,typename TFragPos,typename TFragSize, typename TFragSpec>
//TScoreValue
//_getRefinedMatchScore(Score<TScoreValue, Simple> & score_type,
//		 TStringSet & seqs, 
//		 Fragment<TFragId,TFragPos,TFragSize,TFragSpec> & segment, 
//		 TFragPos pos_i, 
//		 TFragPos pos_j, 
//		 TFragSize len1, 
//		 TFragSize len2)
//{
//SEQAN_CHECKPOINT
//	typename Infix<typename Value<TStringSet>::Type>::Type label0 = label(segment,seqs,sequenceId(segment,0));
//	typename Infix<typename Value<TStringSet>::Type>::Type label1 = label(segment,seqs,sequenceId(segment,1));
//	TScoreValue ret_score = 0;
//	TFragSize len;
//	if (len1 < len2) len = len1;
//	else len = len2;
//	if(len1 <= len2)
//	{
//		ret_score += len1 * scoreMatch(score_type);
//		ret_score += (len2 - len1) * scoreGapExtend(score_type);
//	}
//	else{
//		ret_score += len2 * scoreMatch(score_type);
//		ret_score += (len1 - len2) * scoreGapExtend(score_type);
//	}
//	return ret_score;
//}				


//get score for alignment of length len starting at pos_i on one sequence (first sequence if i_am_first==true)
//and pos_j on other sequence (second sequence if i_am_first==true)
template<typename TScoreValue,typename TScoreSpec,typename TStringSet,typename TFragment,typename TFragPos,typename TFragSize>
TScoreValue
_getRefinedMatchScore(Score<TScoreValue,TScoreSpec> & score_type,
		 TStringSet & seqs,
		 TFragment& segment,
		 TFragPos pos_i,
		 TFragPos pos_j,
		 TFragSize len,
		 TFragSize)
{
SEQAN_CHECKPOINT
	typedef typename Infix<typename Value<TStringSet>::Type>::Type TSegmentLabel;
	TSegmentLabel label0 = label(segment,seqs, sequenceId(segment, 0));
	TSegmentLabel label1 = label(segment,seqs, sequenceId(segment, 1));
	typename Iterator<TSegmentLabel, Rooted>::Type label_it0 = begin(label0) + (pos_i - fragmentBegin(segment,sequenceId(segment,0)));
	typename Iterator<TSegmentLabel, Rooted>::Type label_it1 = begin(label1) + (pos_j - fragmentBegin(segment,sequenceId(segment,1)));
	int i = 0;
	TScoreValue ret_score = 0;
	while(i < (int) len)
	{
		ret_score += score(score_type,*label_it0,*label_it1);
		++label_it0;
		++label_it1;
		++i;
	}
	return ret_score;
}				


//get score for alignment of length len starting at pos_i on one sequence (first sequence if i_am_first==true)
//and pos_j on other sequence (second sequence if i_am_first==true)
template<typename TScoreValue,typename TStringSet,typename TFragPos,typename TFragSize, typename TSpec>
TScoreValue
_getRefinedMatchScore(Score<TScoreValue, Simple> & score_type,
		 TStringSet &,
		 Fragment<TFragSize,ExactFragment<TSpec> > &,
		 TFragPos,
		 TFragPos,
		 TFragSize len,
		 TFragSize)
{
SEQAN_CHECKPOINT
	return len*scoreMatch(score_type);
}				


}  // namespace

#endif  // #ifndef SEQAN_CORE_INCLUDE_SEQAN_GRAPH_ALGORITHM_REFINE_FRAGMENT_H_
