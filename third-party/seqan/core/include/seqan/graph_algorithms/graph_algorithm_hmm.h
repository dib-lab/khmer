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

#ifndef SEQAN_HEADER_GRAPH_ALGORITHM_HMM_H
#define SEQAN_HEADER_GRAPH_ALGORITHM_HMM_H

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////
// Basic HMM algorithms
//////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////

/**
.Function.viterbiAlgorithm:
..cat:Graph
..summary:Implements the Viterbi algorithm
..description:
The Viterbi algorithm computes the most likely sequence of hidden states of the Hidden Markov Model $hmm$ given the sequence $seq$ using dynamic programming.
The result is the most likely sequence of hidden states and returned in $path$.
..signature:TProbability viterbiAlgorithm(hmm, seq, path);
..param.hmm:In-parameter:Input HMM.
...type:Spec.Hmm
..param.seq:In-parameter:Input sequence.
..param.path:Out-parameter:State path.
..returns:Probability of the path, the type parameter $TCargo$ from type of $hmm$.
..see:Function.forwardAlgorithm
..see:Function.backwardAlgorithm
..remarks:
See the @http://en.wikipedia.org/wiki/Viterbi_algorithm|Wikipedia article on the Viterbi algorithm@ for an introduction to the algorithm itself.
..include:seqan/graph_algorithms.h
*/
template<typename TAlphabet, typename TProbability, typename TSpec, typename TSequence, typename TPath>
inline TProbability
viterbiAlgorithm(Graph<Hmm<TAlphabet, TProbability, TSpec> > const& hmm,
				 TSequence const& seq,
				 TPath& path)
{
	SEQAN_CHECKPOINT
	typedef Graph<Hmm<TAlphabet, TProbability, TSpec> > TGraph;
	typedef typename Size<TGraph>::Type TSize;
	typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
	typedef typename Iterator<TGraph, VertexIterator>::Type TVertexIterator;

	// Initialization
	String<TProbability> vMat;
	String<TSize> traceback;
	TSize numCols = length(seq) + 2;
	TSize numRows = getIdUpperBound(_getVertexIdManager(hmm));
	resize(vMat, numCols * numRows, 0.0);
	resize(traceback, numCols * numRows);
	value(vMat, getBeginState(hmm)) = 1.0;
	TVertexDescriptor bState = getBeginState(hmm);
	TVertexDescriptor eState = getEndState(hmm);
	TVertexDescriptor nilVertex = getNil<TVertexDescriptor>();

	// Distinct between silent states and real states
	typedef String<TVertexDescriptor> TStateSet;
	typedef typename Iterator<TStateSet>::Type TStateIter;
	TStateSet silentStates;
	TStateSet realStates;
	TVertexIterator itVertex(hmm);
	for(;!atEnd(itVertex);goNext(itVertex)) {
		if (isSilent(hmm, value(itVertex))) { 
			appendValue(silentStates, value(itVertex));
		} else {
			appendValue(realStates, value(itVertex));
		}
	}

	// Initialization for silent states connected to the begin state
	TStateIter itSilentStateEnd = end(silentStates);
	for(TStateIter itSilentState = begin(silentStates); itSilentState != itSilentStateEnd; goNext(itSilentState)) {
		if ((value(itSilentState) == bState) || (value(itSilentState) == eState)) continue;
		TProbability maxValue = 0.0;
		TVertexDescriptor maxVertex = nilVertex;

		for(TStateIter itBelow = begin(silentStates); itBelow != itSilentState; goNext(itBelow)) {
			TProbability local = value(vMat, value(itBelow)) * getTransitionProbability(hmm, value(itBelow), value(itSilentState));
			if (local > maxValue) {
				maxValue = local;
				maxVertex = value(itBelow);
			}
		}

		// Set traceback vertex
		if (maxVertex != nilVertex) {
			value(vMat, value(itSilentState)) = maxValue;		
			value(traceback, value(itSilentState)) = maxVertex;
		}
	}

	// Recurrence
	TSize len = length(seq);
	for(TSize i=1; i<=len; ++i) {
		// Iterate over real states
		TStateIter itRealStateEnd = end(realStates);
		for(TStateIter itRealState = begin(realStates); itRealState != itRealStateEnd; goNext(itRealState)) {
			// Find maximum
			TProbability maxValue = 0.0;
			TVertexDescriptor maxVertex = nilVertex;
			TVertexIterator itMax(hmm);
			for(;!atEnd(itMax);++itMax) {
				TProbability local = value(vMat, (i-1) * numRows + value(itMax)) * getTransitionProbability(hmm, value(itMax), value(itRealState));
				if (local > maxValue) {
					maxValue = local;
					maxVertex = *itMax;
				}
			}
			// Set traceback vertex
			if (maxVertex != nilVertex) {
				value(vMat, i * numRows + value(itRealState)) = maxValue * getEmissionProbability(hmm, value(itRealState), value(seq, i-1));
				value(traceback, i * numRows + value(itRealState)) = maxVertex;
			}
		}

		// Iterate over silent states
		for(TStateIter itSilentState = begin(silentStates); itSilentState != itSilentStateEnd; goNext(itSilentState)) {
			if ((value(itSilentState) == bState) || (value(itSilentState) == eState)) continue;
			// Find maximum
			TProbability maxValue = 0.0;
			TVertexDescriptor maxVertex = nilVertex;
	
			// Iterate over real states
			for(TStateIter itRealState = begin(realStates); itRealState != itRealStateEnd; goNext(itRealState)) {
				TProbability local = value(vMat, i * numRows + value(itRealState)) * getTransitionProbability(hmm, value(itRealState), value(itSilentState));
				if (local > maxValue) {
					maxValue = local;
					maxVertex = value(itRealState);
				}
			}
			// Iterate over silent states in increasing order
			for(TStateIter itBelow = begin(silentStates); itBelow != itSilentState; goNext(itBelow)) {
				TProbability local = value(vMat, i * numRows + value(itBelow)) * getTransitionProbability(hmm, value(itBelow), value(itSilentState));
				if (local > maxValue) {
					maxValue = local;
					maxVertex = value(itBelow);
				}
			}
			// Set traceback vertex
			if (maxVertex != nilVertex) {
				value(traceback, i * numRows + value(itSilentState)) = maxVertex;
				value(vMat, i * numRows + value(itSilentState)) = maxValue;		
			}
		}
	}

	// Termination
	TProbability maxValue = 0.0;
	TVertexDescriptor maxVertex = 0;
	TVertexIterator itMax(hmm);
	for(;!atEnd(itMax);++itMax) {
		TProbability local = value(vMat, len * numRows + *itMax) * getTransitionProbability(hmm, value(itMax), eState);
		if (local > maxValue) {
			maxValue = local;
			maxVertex = value(itMax);
		}
	}
	value(traceback, (len + 1) * numRows + eState) = maxVertex;
	if (maxVertex != nilVertex) value(vMat, (len+1) * numRows + eState) = maxValue;

	// Traceback
	if (maxValue > 0.0) {
		clear(path);
		TVertexDescriptor oldState = eState;
		appendValue(path, oldState);
		for(TSize i = len + 1; i>=1; --i) {
			do {
				if ((!isSilent(hmm, oldState)) || (oldState == eState)) oldState = value(traceback, i * numRows + oldState);
				else oldState = value(traceback, (i - 1) * numRows + oldState);
				appendValue(path, oldState);
			} while ((isSilent(hmm, oldState)) && (oldState != bState));
		}
		std::reverse(begin(path), end(path));
	}
	
	//// Debug code
	//for(TSize i = 0; i<numRows; ++i) {
	//	for(TSize j=0; j<numCols; ++j) {
	//		std::cout << value(vMat, j*numRows + i) << ',';
	//		//std::cout << value(traceback, j*numRows + i) << ',';
	//	}
	//	std::cout << std::endl;
	//}
	//for(TSize i = 0; i<length(path); ++i) {
	//	std::cout << path[i] << ',';
	//}
	//std::cout << std::endl;

	return value(vMat, (len+1) * numRows + eState);
}


//////////////////////////////////////////////////////////////////////////////

template<typename TAlphabet, typename TProbability, typename TSpec, typename TSequence, typename TForwardMatrix>
inline TProbability
_forwardAlgorithm(Graph<Hmm<TAlphabet, TProbability, TSpec> > const& hmm,
				   TSequence const& seq,
				   TForwardMatrix& fMat)
{
	SEQAN_CHECKPOINT
	typedef Graph<Hmm<TAlphabet, TProbability, TSpec> > TGraph;
	typedef typename Size<TGraph>::Type TSize;
	typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
	typedef typename Iterator<TGraph, VertexIterator>::Type TVertexIterator;

	// Initialization
	TSize numCols = length(seq) + 2;
	TSize numRows = getIdUpperBound(_getVertexIdManager(hmm));
	resize(fMat, numCols * numRows, 0.0);
	value(fMat, getBeginState(hmm)) = 1.0;
	TVertexDescriptor bState = getBeginState(hmm);
	TVertexDescriptor eState = getEndState(hmm);

	// Distinct between silent states and real states
	typedef String<TVertexDescriptor> TStateSet;
	typedef typename Iterator<TStateSet>::Type TStateIter;
	TStateSet silentStates;
	TStateSet realStates;
	TVertexIterator itVertex(hmm);
	for(;!atEnd(itVertex);goNext(itVertex)) {
		if (isSilent(hmm, value(itVertex))) { 
			appendValue(silentStates, value(itVertex));
		} else {
			appendValue(realStates, value(itVertex));
		}
	}

	// Initialization for silent states connected to the begin state
	TStateIter itSilentStateEnd = end(silentStates);
	for(TStateIter itSilentState = begin(silentStates); itSilentState != itSilentStateEnd; goNext(itSilentState)) {
		if ((value(itSilentState) == bState) || (value(itSilentState) == eState)) continue;
		TProbability sumValue = 0.0;
		for(TStateIter itBelow = begin(silentStates); itBelow != itSilentState; goNext(itBelow)) {
			sumValue += value(fMat, value(itBelow)) * getTransitionProbability(hmm, value(itBelow), value(itSilentState));
		}
		value(fMat, value(itSilentState)) = sumValue;		
	}

	// Recurrence
	TSize len = length(seq);
	for(TSize i=1; i<=len; ++i) {
		// Iterate over real states
		TStateIter itRealStateEnd = end(realStates);
		for(TStateIter itRealState = begin(realStates); itRealState != itRealStateEnd; goNext(itRealState)) {
			TProbability sum = 0.0;
			TVertexIterator itAll(hmm);
			for(;!atEnd(itAll);++itAll) sum += value(fMat, (i-1) * numRows + value(itAll)) * getTransitionProbability(hmm, value(itAll), value(itRealState));
			value(fMat, i * numRows + value(itRealState)) = getEmissionProbability(hmm, value(itRealState), seq[i-1]) * sum;
		}

		// Iterate over silent states
		for(TStateIter itSilentState = begin(silentStates); itSilentState != itSilentStateEnd; goNext(itSilentState)) {
			if ((value(itSilentState) == bState) || (value(itSilentState) == eState)) continue;
			TProbability sumValue = 0.0;
			// Iterate over real states
			for(TStateIter itRealState = begin(realStates); itRealState != itRealStateEnd; goNext(itRealState)) {
				sumValue += value(fMat, i * numRows + value(itRealState)) * getTransitionProbability(hmm, value(itRealState), value(itSilentState));
			}
			// Iterate over silent states in increasing order
			for(TStateIter itBelow = begin(silentStates); itBelow != itSilentState; goNext(itBelow)) {
				sumValue += value(fMat, i * numRows + value(itBelow)) * getTransitionProbability(hmm, value(itBelow), value(itSilentState));
			}
			value(fMat, i * numRows + value(itSilentState)) = sumValue;
		}
	}

	// Termination
	TProbability sum = 0.0;
	TVertexIterator itAll(hmm);
	for(;!atEnd(itAll);++itAll) {
		sum += value(fMat, len * numRows + value(itAll)) * getTransitionProbability(hmm, value(itAll), eState);
	}
	value(fMat, (len+1) * numRows + eState) = sum;

	//// Debug code
	//for(TSize i = 0; i<numRows; ++i) {
	//	for(TSize j=0; j<numCols; ++j) {
	//		std::cout << value(fMat, j*numRows + i) << ',';
	//	}
	//	std::cout << std::endl;
	//}

	return value(fMat, (len+1) * numRows + eState);
}

//////////////////////////////////////////////////////////////////////////////

/**
.Function.forwardAlgorithm
..cat:Graph
..summary:Implements the forward algorithm.
..description:Given a Hidden Markov Model $hmm$, the forward algorithm computes the probability of the sequence $seq$.
..signature:TProbability forwardAlgorithm(hmm, seq);
..param.hmm:In-parameter:Input HMM.
...type:Spec.Hmm
..param.seq:In-parameter:Input sequence.
..returns:Probability of the sequence $seq$, $TProbability$ is the type parameter $TCargo$ of the type of $hmm$.
..see:Function.viterbiAlgorithm
..see:Function.backwardAlgorithm
..remarks:See the @http://en.wikipedia.org/wiki/Forward_algorithm|Wikipedia article on the Foward algorithm@ for an introduction to the algorithm itself.
..include:seqan/graph_algorithms.h
*/
template<typename TAlphabet, typename TProbability, typename TSpec, typename TSequence>
inline TProbability
forwardAlgorithm(Graph<Hmm<TAlphabet, TProbability, TSpec> > const& hmm,
				 TSequence const& seq)
{
	SEQAN_CHECKPOINT
	String<TProbability> fMat;
	return _forwardAlgorithm(hmm, seq, fMat);
}

//////////////////////////////////////////////////////////////////////////////

template<typename TAlphabet, typename TProbability, typename TSpec, typename TSequence, typename TBackwardMatrix>
inline TProbability
_backwardAlgorithm(Graph<Hmm<TAlphabet, TProbability, TSpec> > const& hmm,
					TSequence const& seq,
					TBackwardMatrix& bMat)
{
	SEQAN_CHECKPOINT
	typedef Graph<Hmm<TAlphabet, TProbability, TSpec> > TGraph;
	typedef typename Size<TGraph>::Type TSize;
	typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
	typedef typename Iterator<TGraph, VertexIterator>::Type TVertexIterator;

	// Initialization
	TSize numCols = length(seq) + 2;
	TSize numRows = getIdUpperBound(_getVertexIdManager(hmm));
	resize(bMat, numCols * numRows, 0.0);
	TVertexDescriptor bState = getBeginState(hmm);
	TVertexDescriptor eState = getEndState(hmm);
	TSize len = length(seq);
	value(bMat, (len + 1) * numRows + eState) = 1.0;
	
	// Distinct between silent states and real states
	typedef String<TVertexDescriptor> TStateSet;
	typedef typename Iterator<TStateSet>::Type TStateIter;
	TStateSet silentStates;
	TStateSet realStates;
	TVertexIterator itVertex(hmm);
	for(;!atEnd(itVertex);goNext(itVertex)) {
		if (isSilent(hmm, value(itVertex))) { 
			appendValue(silentStates, value(itVertex));
		} else {
			appendValue(realStates, value(itVertex));
		}
	}
	// Reverse the silent states order
	std::reverse(begin(silentStates), end(silentStates));

	// Initialization for silent states connected to the end state
	TStateIter itSilentStateEnd = end(silentStates);
	for(TStateIter itSilentState = begin(silentStates); itSilentState != itSilentStateEnd; goNext(itSilentState)) {
		if (value(itSilentState) == eState) continue;
		TProbability sumValue = getTransitionProbability(hmm, value(itSilentState), eState);
		for(TStateIter itAbove = begin(silentStates); itAbove != itSilentState; goNext(itAbove)) {
			sumValue += value(bMat, len * numRows + value(itAbove)) * getTransitionProbability(hmm, value(itSilentState), value(itAbove));
		}
		value(bMat, len * numRows + value(itSilentState)) = sumValue;
	}

	// Initialization for real states
	TStateIter itRealStateEnd = end(realStates);
	for(TStateIter itRealState = begin(realStates); itRealState != itRealStateEnd; goNext(itRealState)) {
		TProbability sumValue = getTransitionProbability(hmm, value(itRealState), eState);
		for(TStateIter itSilentState = begin(silentStates); itSilentState != itSilentStateEnd; goNext(itSilentState)) {
			sumValue += value(bMat, len * numRows + value(itSilentState)) * getTransitionProbability(hmm, value(itRealState), value(itSilentState));
		}
		value(bMat, len * numRows + value(itRealState)) = sumValue;
	}

	
	// Recurrence
	if (len > 0) {
		for(TSize i=len - 1; i>0; --i) {
			// Iterate over silent states
			for(TStateIter itSilentState = begin(silentStates); itSilentState != itSilentStateEnd; goNext(itSilentState)) {
				if ((value(itSilentState) == bState) || (value(itSilentState) == eState)) continue;
				TProbability sumValue = 0.0;
				// Iterate over real states
				for(TStateIter itRealState = begin(realStates); itRealState != itRealStateEnd; goNext(itRealState)) {
					sumValue += value(bMat, (i+1) * numRows + value(itRealState)) * getTransitionProbability(hmm, value(itSilentState), value(itRealState))* getEmissionProbability(hmm, value(itRealState), value(seq, i));
				}
				// Iterate over silent states in decreasing order
				for(TStateIter itAbove = begin(silentStates); itAbove != itSilentState; goNext(itAbove)) {
					if ((value(itAbove) == bState) || (value(itAbove) == eState)) continue;
					sumValue += value(bMat, i * numRows + value(itAbove)) * getTransitionProbability(hmm, value(itSilentState), value(itAbove));
				}
				value(bMat, i * numRows + value(itSilentState)) = sumValue;
			}

			// Iteration over real states
			for(TStateIter itRealState = begin(realStates); itRealState != itRealStateEnd; goNext(itRealState)) {
				TProbability sumValue = 0.0;
				// Iterate over silent states
				for(TStateIter itSilentState = begin(silentStates); itSilentState != itSilentStateEnd; goNext(itSilentState)) {
					if ((value(itSilentState) == bState) || (value(itSilentState) == eState)) continue;
					sumValue += value(bMat, i * numRows + value(itSilentState)) * getTransitionProbability(hmm, value(itRealState), value(itSilentState));
				}
				// Iterate over real states
				for(TStateIter itR = begin(realStates); itR != itRealStateEnd; goNext(itR)) {
					sumValue += value(bMat, (i+1) * numRows + value(itR)) * getTransitionProbability(hmm, value(itRealState), value(itR)) * getEmissionProbability(hmm, value(itR), value(seq, i));
				}
				value(bMat, i * numRows + value(itRealState)) =  sumValue;
			}
		}
	
		// Termination
		// Iterate over silent states
		for(TStateIter itSilentState = begin(silentStates); itSilentState != itSilentStateEnd; goNext(itSilentState)) {
			if ((value(itSilentState) == bState) || (value(itSilentState) == eState)) continue;
			TProbability sumValue = 0.0;
			// Iterate over real states
			for(TStateIter itRealState = begin(realStates); itRealState != itRealStateEnd; goNext(itRealState)) {
				sumValue += value(bMat, 1 * numRows + value(itRealState)) * getTransitionProbability(hmm, value(itSilentState), value(itRealState))* getEmissionProbability(hmm, value(itRealState), value(seq, 0));
			}
			// Iterate over silent states in decreasing order
			for(TStateIter itAbove = begin(silentStates); itAbove != itSilentState; goNext(itAbove)) {
				if ((value(itAbove) == bState) || (value(itAbove) == eState)) continue;
				sumValue += value(bMat, value(itAbove)) * getTransitionProbability(hmm, value(itSilentState), value(itAbove));
			}
			value(bMat, value(itSilentState)) = sumValue;
		}
		// Sum up all values
		TProbability sumValue = 0.0;
		for(TStateIter itSilentState = begin(silentStates); itSilentState != itSilentStateEnd; goNext(itSilentState)) {
			if ((value(itSilentState) == bState) || (value(itSilentState) == eState)) continue;
			sumValue += value(bMat, value(itSilentState)) * getTransitionProbability(hmm, bState, value(itSilentState));
		}
		for(TStateIter itRealState = begin(realStates); itRealState != itRealStateEnd; goNext(itRealState)) {
			sumValue += value(bMat, 1 * numRows + value(itRealState)) * getTransitionProbability(hmm, bState, value(itRealState)) * getEmissionProbability(hmm, value(itRealState), value(seq, 0));
		}
		value(bMat, bState) = sumValue;
	}

	//// Debug code
	//for(TSize i = 0; i<numRows; ++i) {
	//	for(TSize j=0; j<numCols; ++j) {
	//		std::cout << value(bMat, j*numRows + i) << ',';
	//	}
	//	std::cout << std::endl;
	//}

	return value(bMat, bState);
}
	


//////////////////////////////////////////////////////////////////////////////

/**
.Function.backwardAlgorithm
..cat:Graph
..summary:Implements the backward algorithm.
..description:Execute the backward algorithm on the Hidden Markov Model $hmm$ to get the probability of $seq$.
..signature:TProbability backwardAlgorithm(hmm, seq);
..param.hmm:In-parameter:Input HMM.
...type:Spec.Hmm
..param.seq:In-parameter:Input sequence.
..returns:Probability of the sequence $seq$.  $TProbability$ is the type parameter $TCargo$ of the type of $hmm$.
..see:Function.viterbiAlgorithm
..see:Function.forwardAlgorithm
..remarks:See the @http://en.wikipedia.org/wiki/Forward-backward_algorithm|Wikipedia article on the Forward-backward algorithm@ for an introduction to the algorithm.
..include:seqan/graph_algorithms.h
*/
template<typename TAlphabet, typename TProbability, typename TSpec, typename TSequence>
inline TProbability
backwardAlgorithm(Graph<Hmm<TAlphabet, TProbability, TSpec> > const& hmm,
				  TSequence const& seq)
{
	SEQAN_CHECKPOINT
	String<TProbability> bMat;
	return _backwardAlgorithm(hmm, seq, bMat);
}

///////////////////////////////////////////////////////////////////////////////////////

/**
.Function.generateSequence:
..cat:Graph
..summary:Generates random state and alphabet sequences of a given HMM.
..signature:generateSequence(hmm, sequences, states, numSeq, maxLength)
..param.hmm:In-parameter:Input HMM.
...type:Spec.Hmm
..param.sequences:The StringSet of alphabet sequences.
...type:Class.StringSet
..param.sequences:The StringSet of state sequences.
...type:Class.StringSet
..param.numSeq:The number of sequences to generate.
..param.maxLength:The maximum length of the sequences.
...remarks:Sequences might be shorter if the end state is reached prior to maxLength. 
..remarks: Because of silent states, generated alphabet and state sequences might have different length.
..returns:void
..include:seqan/graph_algorithms.h
*/
template<typename TAlphabet, typename TProbability, typename TSpec,typename TSequenceSet, typename TStateSeqSet, typename TSize>
inline void
generateSequence(Graph<Hmm<TAlphabet, TProbability, TSpec> > const& hmm,
				 TSequenceSet& sequences,
				 TStateSeqSet& states,
				 TSize numSeq,
				 TSize maxLength) 
{
	SEQAN_CHECKPOINT
	typedef Graph<Hmm<TAlphabet, TProbability, TSpec> > TGraph;
	typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
	typedef typename Iterator<TGraph, OutEdgeIterator>::Type TOutEdgeIterator;
	typedef typename Value<TSequenceSet>::Type TSequence;
	typedef typename Value<TStateSeqSet>::Type TStateSeq;

    // TODO(holtgrew): It should be possible to use a per-call RNG.
	typedef typename GetDefaultRng<TGraph>::Type TRng;
	typedef Pdf<Uniform<double> > TPdf;
	TPdf pdf(0.0, 1.0);
	TRng & rng = defaultRng(TGraph());
	
	// Initialization
	clear(sequences);
	clear(states);
	TSize alphSize = ValueSize<TAlphabet>::VALUE;
	
	// Simulate sequences
	TVertexDescriptor currentState;
	TVertexDescriptor endState = getEndState(hmm);
	for(TSize i=0;i<numSeq;++i){
		currentState = getBeginState(hmm);
		TSequence seq;
		TStateSeq stat;
		appendValue(stat, getBeginState(hmm));
		bool stop = false;
		TSize pos = 0; 
		while (pos < maxLength) {
			TProbability prob = pickRandomNumber(rng, pdf);
			TProbability compareProb = 0.0;
			TOutEdgeIterator itState(hmm, currentState);
			// Determine the next state
			for(;!atEnd(itState);++itState) {
				// Probability of the next transition
				compareProb += getTransitionProbability(hmm, value(itState));
				// Compare with random probability
				if (prob <= compareProb){
					TVertexDescriptor nextState = targetVertex(hmm, value(itState));
					if (nextState == endState) {
						stop = true;
						break;
					}
					appendValue(stat, nextState);
					if (!isSilent(hmm, nextState)) {
						compareProb =0.0;
						prob = pickRandomNumber(rng, pdf);
						for (TSize c=0;c<alphSize;++c){
							compareProb += getEmissionProbability(hmm,targetVertex(hmm, value(itState)), TAlphabet(c));
							if (prob <= compareProb) {
								appendValue(seq, TAlphabet(c));
								++pos;
								break;
							}
						}
					}
					currentState = nextState;
					break;
				}
			}
			if (stop==true) break;
		}
		appendValue(stat, getEndState(hmm));
		appendValue(sequences, seq);
		appendValue(states, stat);
	}
}

///////////////////////////////////////////////////////////////////////////////////////

template<typename TAlphabet, typename TProbability, typename TSpec,typename TSequenceSet, typename TSize>
inline void
generateSequence(Graph<Hmm<TAlphabet, TProbability, TSpec> > const& hmm,
				 TSequenceSet& sequences,
				 TSize numSeq,
				 TSize maxLength) 
{
	SEQAN_CHECKPOINT
	typedef Graph<Hmm<TAlphabet, TProbability, TSpec> > TGraph;
	typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
	StringSet<String<TVertexDescriptor> > states;
	generateSequence(hmm, sequences, states, numSeq, maxLength);
}



//////////////////////////////////////////////////////////////////////////////
// Training algorithms
//////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////

template<typename TAlphabet, typename TProbability, typename TSpec,typename TEmissionCounter, typename TTransitionCounter>
inline void
_parameterEstimator(Graph<Hmm<TAlphabet, TProbability, TSpec> >& hmm,
					 TEmissionCounter const& emission,
					 TTransitionCounter const& transition)
{
	typedef Graph<Hmm<TAlphabet, TProbability, TSpec> > TGraph;
	typedef typename Size<TGraph>::Type TSize;
	typedef typename Iterator<TGraph, VertexIterator>::Type TVertexIterator;
	typedef typename Iterator<TGraph, OutEdgeIterator>::Type TOutEdgeIterator;
	typedef typename Value<TTransitionCounter>::Type TCounterValue;

	// Initialization
	TSize alphSize = ValueSize<TAlphabet>::VALUE;
	TCounterValue pseudoCount = (TCounterValue) 1.0;

	// Estimate the parameters from the counter values
	TVertexIterator itAll(hmm);
	for(;!atEnd(itAll);++itAll){
		if (!isSilent(hmm, value(itAll))) {
			TCounterValue summedCount = (TCounterValue) (0.0);
			for(TSize i=0; i<alphSize;++i) summedCount += (pseudoCount + value(value(emission, value(itAll)),i));
			for(TSize i=0; i<alphSize;++i) emissionProbability(hmm,value(itAll),TAlphabet(i)) = (pseudoCount + value(value(emission, value(itAll)),i)) / summedCount;
		}
		TCounterValue summedCount = (TCounterValue) (0.0);
		TOutEdgeIterator itOutSum(hmm,value(itAll));
		for(;!atEnd(itOutSum);++itOutSum) summedCount += (pseudoCount + getProperty(transition, value(itOutSum)));
		TOutEdgeIterator itOut(hmm, value(itAll));
		for(;!atEnd(itOut);++itOut) transitionProbability(hmm, value(itOut)) = (pseudoCount + getProperty(transition, value(itOut))) / summedCount;
	}
}

///////////////////////////////////////////////////////////////////////////////////////

template<typename TAlphabet, typename TProbability, typename TSpec, typename TSequenceSet, typename TStateSeqSet>
inline void 
estimationWithStates(Graph<Hmm<TAlphabet, TProbability, TSpec> >& hmm,
					 TSequenceSet& sequences,
					 TStateSeqSet& states)
{
	typedef Graph<Hmm<TAlphabet, TProbability, TSpec> > TGraph;
	typedef typename Size<TGraph>::Type TSize;
	typedef typename Value<TStateSeqSet>::Type TStateSeq;
	typedef typename Value<TStateSeq>::Type TState;
	typedef typename EdgeDescriptor<TGraph>::Type TEdgeDescriptor;

	// Initialization
	TSize alphSize = ValueSize<TAlphabet>::VALUE;
	typedef String<TSize> TCountString;
	TCountString transitionCounter;
	resize(transitionCounter, getIdUpperBound(_getEdgeIdManager(hmm)), 0);
	StringSet<TCountString> emissionCounter;
	TSize numRows = getIdUpperBound(_getVertexIdManager(hmm));
	for(TSize i =0; i<numRows;++i) {
		TCountString emisCount;
		resize(emisCount,alphSize,0);
		appendValue(emissionCounter,emisCount);
	}
	
	// Iterate over all sequences
	for (TSize j=0;j<length(states);++j) {
		TSize posSeq = 0;
		for (TSize i = 0;i<length(value(states,j));++i) {
			TState s = value(value(states,j),i);
			if (!isSilent(hmm, s)) {
				TAlphabet c = value(value(sequences,j),posSeq);
				value(value(emissionCounter, s), ordValue(c)) += 1;
				++posSeq;
			}
			if (i<(length(value(states,j))-1)) {
				TEdgeDescriptor currEdge = findEdge(hmm, s, value(value(states,j), (i+1)));
				property(transitionCounter, currEdge) += 1;
			}
		}
	}

	//// Debug Code
	//typedef typename Iterator<TGraph, EdgeIterator>::Type TEdgeIterator;
	//TEdgeIterator itE(hmm);
	//for(;!atEnd(itE);goNext(itE)) std::cout << sourceVertex(itE) << ',' << targetVertex(itE) << ':' << property(transitionCounter, value(itE)) << std::endl;
	//for(TSize j=0; j<length(emissionCounter); ++j) {
	//	for(TSize i=0;i<length(value(emissionCounter, j));++i) {
	//		std::cout << j << ":" << TAlphabet(i) << '=' << value(value(emissionCounter, j), i) << std::endl;
	//	}
	//}

	// Estimate Parameters
	_parameterEstimator(hmm,emissionCounter, transitionCounter);
}

///////////////////////////////////////////////////////////////////////////////////////

template<typename TAlphabet, typename TProbability, typename TSpec>
inline void
_fillHmmUniform(Graph<Hmm<TAlphabet, TProbability, TSpec> >& hmm)
{
	SEQAN_CHECKPOINT
	typedef Graph<Hmm<TAlphabet, TProbability, TSpec> > TGraph;
	typedef typename Size<TGraph>::Type TSize;
	typedef typename Iterator<TGraph, VertexIterator>::Type TVertexIterator;
	typedef typename Iterator<TGraph, OutEdgeIterator>::Type TOutEdgeIterator;
	
	// Initialization
	TSize alphSize = ValueSize<TAlphabet>::VALUE;
	
	// Iterate over all states
	TVertexIterator itState(hmm);
	for(;!atEnd(itState);goNext(itState)) {		//pass through the states of the hmm	
		TSize oD = outDegree(hmm, value(itState));
		TOutEdgeIterator itOut(hmm, value(itState));
		for (;!atEnd(itOut);goNext(itOut)) transitionProbability(hmm, value(itOut)) = (TProbability) (1.0 / (double) (oD));
		if (!isSilent(hmm, value(itState))) {
			for(TSize i=0;i<alphSize;++i){
				emissionProbability(hmm,value(itState),TAlphabet(i)) = (TProbability) (1.0 / (double) (alphSize));
			}
		}
	}
}

///////////////////////////////////////////////////////////////////////////////////////

template<typename TAlphabet, typename TProbability, typename TSpec, typename TRNG>
inline void
_fillHmmRandom(Graph<Hmm<TAlphabet, TProbability, TSpec> >& hmm,
                TRNG & rng)
{
	SEQAN_CHECKPOINT
	typedef Graph<Hmm<TAlphabet, TProbability, TSpec> > TGraph;
	typedef typename Size<TGraph>::Type TSize;
	typedef typename Iterator<TGraph, VertexIterator>::Type TVertexIterator;
	typedef typename Iterator<TGraph, OutEdgeIterator>::Type TOutEdgeIterator;
	
	// Initialization
    Pdf<Uniform<TSize> > pdf(20, 99);
	TSize alphSize = ValueSize<TAlphabet>::VALUE;
	
	// Iterate over all states
	TVertexIterator itState(hmm);
	for(;!atEnd(itState);goNext(itState)) {		//pass through the states of the hmm	
		TSize oD = outDegree(hmm, value(itState));
		if (oD > 0) {
			String<TSize> counts;
			TSize sum = 0;
			for(TSize i = 0;i<oD;++i) {
				TSize rd = pickRandomNumber(rng, pdf);
				sum += rd;
				appendValue(counts, rd);
			}
			TSize pos = 0;
			TOutEdgeIterator itOut(hmm, value(itState));
			for (;!atEnd(itOut);goNext(itOut), ++pos) transitionProbability(hmm, value(itOut)) = (TProbability) ((double) (value(counts, pos)) / (double) (sum));
		}	
		if (!isSilent(hmm, value(itState))) {
			String<TSize> counts;
			TSize sum = 0;
			for(TSize i = 0;i<alphSize;++i) {
				TSize rd = pickRandomNumber(rng, pdf);
				sum += rd;
				appendValue(counts, rd);
			}
			for(TSize i=0;i<alphSize;++i) emissionProbability(hmm,value(itState),TAlphabet(i)) = (TProbability) ((double) (value(counts, i)) / (double) (sum));
		}
	}
}

///////////////////////////////////////////////////////////////////////////////////////

template<typename TAlphabet, typename TProbability, typename TSpec, typename TRNG>
inline void
randomizeHmm(Graph<Hmm<TAlphabet, TProbability, TSpec> >& hmm,
             TRNG & /*rng*/)
{
	// TODO(holtgrew): Why uniform and not random?
	//_fillHmmRandom(hmm, rng);
	_fillHmmUniform(hmm);
}

///////////////////////////////////////////////////////////////////////////////////////

template <typename TAlphabet, typename TProbability, typename TSpec, typename TSequence, typename TSize>
inline TProbability
_baumWelchAlgorithm(Graph<Hmm<TAlphabet, TProbability, TSpec > >& hmm,
					 StringSet<TSequence> const& seqSet,
					 TSize maxIter,
					 TProbability epsilon)
{
	typedef Graph<Hmm<TAlphabet, TProbability, TSpec> > TGraph;
	typedef typename Iterator<TGraph, VertexIterator>::Type TVertexIterator;
	typedef typename Iterator<TGraph, OutEdgeIterator>::Type TOutEdgeIterator;
	typedef String<TProbability> TCountString;

	// Initialization
	TSize alphSize = ValueSize<TAlphabet>::VALUE;
	TProbability lastTotalModelProb = 0.0;

	// Randomize current HMM so only topology is preserved
	randomizeHmm(hmm);

	// Iterative Optimization
	for(TSize iter=0; iter<maxIter; ++iter){
		std::cout << "Iteration: "<< iter << std::endl;
		std::cout << hmm << std::endl;
		TCountString transitionCounter;
		resize(transitionCounter, getIdUpperBound(_getEdgeIdManager(hmm)), 0.0);
		StringSet<TCountString> emissionCounter;
		TSize numRows = getIdUpperBound(_getVertexIdManager(hmm));
		for(TSize i =0; i<numRows;++i) {
			TCountString emisCount;
			resize(emisCount,alphSize, 0.0);
			appendValue(emissionCounter,emisCount);
		}
						
		// MaximationStep
		// Determine new contributions
		TProbability totalModelLogProb = 0.0;
		for (TSize i=0; i <length(seqSet);++i){		//sequences
			TSize len = length(value(seqSet,i));

			// Forward algorithm
			String<TProbability> fMat;
			TProbability modelLogProb = _forwardAlgorithm(hmm, value(seqSet,i), fMat);
			totalModelLogProb += modelLogProb;
			
			// Backward algorithm
			String<TProbability> bMat;
			_backwardAlgorithm(hmm, value(seqSet,i), bMat);
			
			// Use the posterior probabilities to estimate counter values
			for (TSize j=0;j<len;++j){  
				TAlphabet c = value(value(seqSet,i),j);
				TAlphabet nextC = c;
				if (j < (len - 1))  nextC = value(value(seqSet,i),(j+1));
				
				// Iterate over all states
				TVertexIterator itAll(hmm);		
				for(;!atEnd(itAll);++itAll) {
					// Handle begin state
					if (value(itAll) == beginState(hmm)) {
						if (j == 0) {
							TOutEdgeIterator itOut(hmm, value(itAll));
							for(;!atEnd(itOut); goNext(itOut)) {
								if (!isSilent(hmm, targetVertex(itOut))) {
									property(transitionCounter, value(itOut)) += (value(fMat, beginState(hmm)) * getTransitionProbability(hmm, value(itOut)) * getEmissionProbability(hmm, targetVertex(itOut), c) * value(bMat, numRows + targetVertex(itOut)) / modelLogProb);
								} else {
									property(transitionCounter, value(itOut)) += (value(fMat, beginState(hmm)) * getTransitionProbability(hmm, value(itOut)) * value(bMat, targetVertex(itOut)) / modelLogProb);
								}
							}
						
						}
						continue;
					}
					// Ignore the end state
					if (value(itAll) == endState(hmm)) continue;

					// Determine emission expectation values
					if (!isSilent(hmm, value(itAll))) value(value(emissionCounter, value(itAll)),ordValue(c)) += ((value(fMat, (j+1) * numRows + value(itAll)) * value(bMat, (j+1) * numRows + value(itAll))) / modelLogProb);
					
					// Determine transition expectation values
					TOutEdgeIterator itOut(hmm, value(itAll));
					for(;!atEnd(itOut); goNext(itOut)) {
						// Handle the end state
						if ((j == (len - 1)) && (targetVertex(itOut)==endState(hmm))) {
							property(transitionCounter, value(itOut)) += (value(fMat, (j+1) * numRows + value(itAll)) * getTransitionProbability(hmm, value(itAll),endState(hmm)) / modelLogProb);
						} else {
							if (!isSilent(hmm, targetVertex(itOut))) {
								if (j < (len - 1)) property(transitionCounter, value(itOut)) += (value(fMat, (j+1) * numRows + value(itAll)) * getTransitionProbability(hmm, value(itOut)) * getEmissionProbability(hmm, targetVertex(itOut), nextC) * value(bMat, (j+2) * numRows + targetVertex(itOut)) / modelLogProb);
							} else {
								property(transitionCounter, value(itOut)) += (value(fMat, (j+1) * numRows + value(itAll)) * getTransitionProbability(hmm, value(itOut)) * value(bMat, (j+1) * numRows + targetVertex(itOut)) / modelLogProb);
							}
						}
					}
				}
			}
		}
		// Expectation step
		_parameterEstimator(hmm,emissionCounter, transitionCounter);
		

		// Termination?
		if ((iter > 5) && ((totalModelLogProb - lastTotalModelProb) < epsilon)) {
			lastTotalModelProb = totalModelLogProb;
			break;
		} else {
			lastTotalModelProb = totalModelLogProb;
		}
	}
	return lastTotalModelProb;
}

///////////////////////////////////////////////////////////////////////////////////////

template <typename TAlphabet, typename TProbability, typename TSpec, typename TSequence>
inline TProbability 
baumWelchAlgorithm(Graph<Hmm<TAlphabet, TProbability, TSpec > >& hmm,
				   StringSet<TSequence> const& seqSet)
{
	SEQAN_CHECKPOINT
	typedef Graph<Hmm<TAlphabet, TProbability, TSpec> > TGraph;
	typedef typename Size<TGraph>::Type TSize;
	TSize maxIter = 100;
	TProbability epsilon = 0.00001;
	return _baumWelchAlgorithm(hmm, seqSet, maxIter, epsilon);
}

/*
//////////////////////////////////////////////////////////////////////////////
// Profile HMMs
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////


template < typename TAlphabet, typename TProbability, typename TSpec, typename TMat, typename TConsensus, typename TSize>
inline void
_profileHmmCounter(Graph<Hmm<TAlphabet, TProbability, TSpec> >& pHmm,
					TMat const& matr,
					TConsensus const& consensus,
					TSize const& numRows,
					TSize const& numCols)
{
	SEQAN_CHECKPOINT
	typedef Graph<Hmm<TAlphabet, TProbability, TSpec> > TGraph;
	typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
	typedef typename EdgeDescriptor<TGraph>::Type TEdgeDescriptor;
	typedef typename Iterator<TGraph, VertexIterator>::Type TVertexIterator;
	typedef typename Iterator<TGraph, OutEdgeIterator>::Type TOutEdgeIterator;
	typedef typename Value<TConsensus>::Type TValue;

	// Initialization
	TSize alphSize = ValueSize<TAlphabet>::VALUE;
	TValue gapChar = gapValue<TValue>();
	TVertexDescriptor begState = beginState(pHmm);
	TVertexDescriptor eState = endState(pHmm);
	TVertexDescriptor matchState = 1;
	TVertexDescriptor insertState = 2;
	TVertexDescriptor deleteState = 3;

	typedef String<TProbability> TCountString;
	TCountString transitionCounter;
	resize(transitionCounter, getIdUpperBound(_getEdgeIdManager(pHmm)), 0.0);
	StringSet<TCountString> emissionCounter;
	TSize nR = getIdUpperBound(_getVertexIdManager(pHmm));
	for(TSize i =0; i<nR;++i) {
		TCountString emisCount;
		resize(emisCount,alphSize, 0.0);
		appendValue(emissionCounter,emisCount);
	}


	String<TVertexDescriptor> oldCol;
	String<TVertexDescriptor> currCol;
	resize(oldCol, numRows, begState);
	resize(currCol, numRows, 0);
	TEdgeDescriptor currEdge;
	
	for(TSize i = 0; i<length(consensus); ++i) {
		//transitionvalues
		
		//being in insertState
		if (value(consensus, i)==gapChar){
			for(TSize j = 0; j<numRows; ++j){
				if ((value(matr, j * numCols + i)!=gapChar)){
					value(currCol,j) = insertState; 
					currEdge = findEdge(pHmm, value(oldCol,j), value(currCol,j));
					property(transitionCounter, currEdge) += 1;
					value(oldCol,j) = value(currCol,j);
				}
				else if ((value(matr, j * numCols + i)==gapChar) ){
					value(currCol,j) = value(oldCol,j);
				}
			}
		}
		//being in normal State
		else{
			for(TSize j = 0; j<numRows; ++j){
				if(value(matr, j * numCols + i)!=gapChar) value(currCol,j) = matchState; 
				else value(currCol,j) = deleteState;
				currEdge = findEdge(pHmm, value(oldCol,j), value(currCol,j));
				property(transitionCounter, currEdge) += 1;
				value(oldCol,j)=value(currCol,j);
			}
		}
		//emissionvalues
		if (value(consensus, i)==gapChar) {
			for(TSize j = 0; j<numRows; ++j) 
				if(value(matr, j * numCols + i)!=gapChar) 
					value(value(emissionCounter, insertState), ordValue( (TAlphabet)  value(matr, j * numCols + i) ) ) += 1;
			
			continue;
		}
		else 
			for(TSize j = 0; j<numRows; ++j) 
				if(value(matr, j * numCols + i)!=gapChar) 
					value(value(emissionCounter, matchState), ordValue( (TAlphabet)  value(matr, j * numCols + i) ) ) += 1;
		
		matchState+=3;
		if ((insertState+3)==eState) insertState+=2;
		else insertState+=3;
		deleteState+=3;
	}
	
	//transition in endState
	for(TSize j = 0; j<numRows; ++j){
		currEdge = findEdge(pHmm,value(currCol,j),eState);
		property(transitionCounter, currEdge) += 1;
	}
	
	_parameterEstimator(pHmm, emissionCounter, transitionCounter);
}


//////////////////////////////////////////////////////////////////////////////

template<typename TAlphabet, typename TCargo, typename TSpec, typename TConsensus>
inline void
_createProfileHmm(Graph<Hmm<TAlphabet, TCargo, TSpec> >& pHmm,
				   TConsensus const& consensus)
{
	SEQAN_CHECKPOINT
	typedef Graph<Hmm<TAlphabet, TCargo, TSpec> > TGraph;
	typedef typename Size<TGraph>::Type TSize;
	typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
	typedef typename EdgeDescriptor<TGraph>::Type TEdgeDescriptor;
	typedef typename Value<TConsensus>::Type TValue;
	
	// Initialization
	TSize alphSize = ValueSize<TAlphabet>::VALUE;
	TValue gapChar = gapValue<TValue>();
	clear(pHmm);
	
	// Add begin state
	TVertexDescriptor begState = addVertex(pHmm);
	assignBeginState(pHmm, begState);
	
	// Add for each consensus letter 3 states
	for (TSize i=0;i<length(consensus);++i){
		if (value(consensus,i) == gapChar) continue;
		addVertex(pHmm); addVertex(pHmm); addVertex(pHmm, true);
	}
	
	// Add last insertion state
	TVertexDescriptor lastIState = addVertex(pHmm);
	
	// Add end state
	TVertexDescriptor endState = addVertex(pHmm);
	assignEndState(pHmm, endState);
	
	// Is there no consensus letter?
	if (lastIState == 1) {
		clear(pHmm);
		return;
	}

	// Remember the kind of state
	TVertexDescriptor mState = 1;
	TVertexDescriptor iState = 2;
	TVertexDescriptor dState = 3;

	// Add tranistions from begin state
	addEdge(pHmm, begState, mState);
	addEdge(pHmm, begState, iState);
	addEdge(pHmm, begState, dState);

	// Add all remaining transitions
	for (TSize i=0;i<length(consensus);++i){	
		if (value(consensus,i) == gapChar) continue;
		else if ((mState + 3) == lastIState) {
			addEdge(pHmm, iState, mState);	
			addEdge(pHmm, iState, iState);
			addEdge(pHmm, iState, dState);
			break;
		}
		else{
			addEdge(pHmm, mState, (mState+3));
			addEdge(pHmm, iState, mState);
			addEdge(pHmm, dState, (mState+3));
			addEdge(pHmm, mState, (iState+3));
			addEdge(pHmm, iState, iState);
			addEdge(pHmm, dState, (iState+3));
			addEdge(pHmm, mState, (dState+3));
			addEdge(pHmm, iState, dState);
			addEdge(pHmm, dState, (dState+3));
			mState+=3;
			iState+=3;
			dState+=3;
		}
	}

	// Transitions to the endState and the last I-state
	addEdge(pHmm, mState, endState);
	addEdge(pHmm, mState, lastIState);
	addEdge(pHmm, lastIState, endState);
	addEdge(pHmm, lastIState, lastIState);
	addEdge(pHmm, dState, endState);
	addEdge(pHmm, dState, lastIState);
}

//////////////////////////////////////////////////////////////////////////////

//e.g.
//String<char> matr = "-AT---GAG-G-AG-CT-C--A--GT-G-CT---G";
//msaToProfileHmm(matr, hmm, 5);

template<typename TAlignmentChar, typename TAlphabet, typename TProbability, typename TSpec, typename TSize>
inline void 
msaToProfileHmm(String<TAlignmentChar> const& matr,
				Graph<Hmm<TAlphabet, TProbability, TSpec> >& pHmm,
				TSize nSeq)
{
	SEQAN_CHECKPOINT
	typedef Graph<Hmm<TAlphabet, TProbability, TSpec> > THmm;
	
	// Consensus 
	String<unsigned int> coverage;
	String<char> gappedConsensus;
	String<Dna> consensusSequence;
	consensusCalling(matr, consensusSequence, gappedConsensus, coverage, nSeq, MajorityVote() );

	// Build the HMM topology
	_createProfileHmm(pHmm,gappedConsensus);

	// Parameterize the pHmm
	TSize numCols = length(matr) / nSeq;
	_profileHmmCounter(pHmm, matr, gappedConsensus, nSeq, numCols);
}

*/


}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
