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
// Author: David Weese <david.weese@fu-berlin.de>
// ==========================================================================

#ifndef SEQAN_HEADER_STORE_ALIGN_INTERVALS_H
#define SEQAN_HEADER_STORE_ALIGN_INTERVALS_H

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////
// MatchedIntervals Store
//////////////////////////////////////////////////////////////////////////////

template<typename TValue = int >
struct Interval
{
	TValue i1;
	TValue i2;
	
	Interval() : i1(0), i2(0) {}
};

template<typename TValue>
inline bool
operator == (Interval<TValue> const & interval1, Interval<TValue> const & interval2)
{
	return (interval1.i1 == interval2.i1 && interval1.i2 == interval2.i2);
}

template<typename TValue>
inline bool
operator < (Interval<TValue> const & interval1, Interval<TValue> const & interval2)
{
    if (interval1.i1 == interval2.i1)
        return (interval1.i2 < interval2.i2);
    else
        return (interval1.i1 < interval2.i1);
}

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////
template <typename TInterval = Interval<unsigned>, typename TSpec = void>         
struct AlignIntervalsStoreElement
{
	typedef  String<TInterval> 	TIntervals;
	
	TIntervals intervals; 
};


///////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////
template<typename TAlignIntervalsStore, typename TSpec, typename TConfig >
inline void
buildAlignIntervalsStore(TAlignIntervalsStore & alignIntervalsStore,
                         FragmentStore<TSpec, TConfig> & me, const unsigned & thresholdGaps)
{
	typedef typename FragmentStore<TSpec, TConfig>::TAlignedReadStore 	TAlignedReadStore;

	typedef typename Iterator<TAlignedReadStore, Standard>::Type 		TAlignIter;
	typedef typename Iterator<TAlignIntervalsStore>::Type 			TAlignIntervalsStoreIter;

	if (!empty(me.alignedReadStore))
	{
		resize(alignIntervalsStore, length(me.alignedReadStore));

		TAlignIntervalsStoreIter it = begin(alignIntervalsStore);

		TAlignIter itAlign = begin(me.alignedReadStore);
		TAlignIter itAlignEnd = end(me.alignedReadStore);

		// get matched intervals for each aligned read
		for ( ; itAlign != itAlignEnd; goNext(itAlign), goNext(it))
		{
			extractAlignIntervals(value(it).intervals, getValue(itAlign), me, thresholdGaps);
		}
	}
}


////////////////////////////////////////////////////////////////////////////////////////////////////////
/////// extract contig intervals of aligned read
////////////////////////////////////////////////////////////////////////////////////////////////////////

template<typename TIntervals, typename TAlignedReadStoreElement, typename TSpec, typename TConfig>
inline void 
extractAlignIntervals(TIntervals & contigIntervals, TAlignedReadStoreElement & align, FragmentStore<TSpec, TConfig> & me, const unsigned & thresholdGaps)
{

	////////////////////////////////////////////////////////////////////////////////////////////////////////
	///////// get read intervals in gapped contig sequence:		
	
	typedef typename FragmentStore<TSpec, TConfig>::TContigPos 		TContigPos;
	typedef typename FragmentStore<TSpec, TConfig>::TReadGapAnchor 		TReadGapAnchor;
	typedef 	 String<TReadGapAnchor> 				TReadGaps;
	typedef typename Position<TReadGaps >::Type 				TGapPos;
	typedef typename Value<TIntervals>::Type 				TInterval;
	typedef typename Iterator<TReadGaps>::Type 				TReadGapsIter;
	
	clear(contigIntervals);
	
	TIntervals readIntervals;
	
	TContigPos beginPos = align.beginPos;
	TContigPos endPos = align.endPos;
	if (beginPos > endPos)
	{
		TContigPos help = beginPos;
		beginPos = endPos;
		endPos = help;
	}
	TInterval interval;
	if (length(align.gaps) == 0)
	{
		resize(readIntervals, 1, Generous());
		interval.i1 = beginPos;
		interval.i2 = endPos -1;
		assignValue(readIntervals, 0, interval);
	}
	else 
	{
		TReadGapsIter itGap = begin(align.gaps);
		TReadGapsIter itGapEnd = end(align.gaps);
	
		TContigPos firstSeqPos = getValue(itGap).seqPos + beginPos;   // seqPos from Read projected to pos in ungapped contig sequence
		TContigPos firstGapPos = getValue(itGap).gapPos + beginPos;   // gapPos from Read projected to pos in ungapped contig sequence
		TGapPos i = 0;
		// without gaps at the beginnig
		if (beginPos < firstGapPos) 		// end != 0 ??
		{
			resize(readIntervals, length(align.gaps), Generous());
			interval.i1 = beginPos; 
			interval.i2 = firstSeqPos - 1;
			assignValue(readIntervals, i, interval);
			++i;
		}
		// with gaps at the beginning
		else resize(readIntervals, length(align.gaps) -1, Generous());
		
		TContigPos seqPos1 = 0;
		TContigPos gapPos1 = 0;
		TContigPos seqPos2 = firstSeqPos;
		TContigPos gapPos2 = firstGapPos;
		goNext(itGap);
		// calculate interval between 2 gaps:
		for (; itGap != itGapEnd; goNext(itGap), ++i)
		{
			seqPos1 = seqPos2;
			gapPos1 = gapPos2;
			seqPos2 = getValue(itGap).seqPos + beginPos;
			gapPos2 = getValue(itGap).gapPos + beginPos;
			interval.i1 = gapPos1;					// ==  gapPos of Gap-Anchor1  + beginPos 
			interval.i2 = seqPos2 + gapPos1 - seqPos1 - 1;   	// ==  gapPos2 - ((gapPos2-seqPos2) - (gapPos1-seqPos1)) - 1
			assignValue(readIntervals, i, interval);
		}
		
		// without gaps at the end
		if (endPos != gapPos2)
		{
			interval.i1 = gapPos2;
			interval.i2 = endPos - 1;
			appendValue(readIntervals, interval, Generous());
		}
		
		// with gaps at the end
		else if (endPos == gapPos2)
		{
			interval.i1 = gapPos1;
			interval.i2 = endPos - 1;
			assignValue(readIntervals, i - 1, interval);
		}
	}
	
	////////////////////////////////////////////////////////////////////////////////////////////////////////
	///////// get contig intervals in ungapped sequence:
	
	typedef typename TAlignedReadStoreElement::TId 				TId;
	typedef typename FragmentStore<TSpec, TConfig>::TContigGapAnchor 	TContigGapAnchor;
	typedef		 String<TContigGapAnchor> 				TContigGaps;
	typedef typename Iterator<TContigGaps const>::Type 			TContigGapsIter;
	typedef typename Iterator<TIntervals>::Type 				TReadIntervalIter;
	
	
	TId contigId = align.contigId;
	
	if (empty(getValue(me.contigStore, contigId).gaps))
	{
		contigIntervals = readIntervals;
	}
	else
	{
		if (front(getValue(me.contigStore, contigId).gaps).seqPos != 0) 
		{
			TContigGapAnchor pseudoGapAnchor;
			pseudoGapAnchor.seqPos = 0;
			pseudoGapAnchor.gapPos = 0;
			insertValue(value(me.contigStore, contigId).gaps, 0, pseudoGapAnchor, Generous());
		}
		
		TReadIntervalIter itI = begin(readIntervals);
		TReadIntervalIter itIEnd = end(readIntervals);
		
		TContigGapsIter it1Gap = begin(getValue(me.contigStore, contigId).gaps);
		TContigGapsIter it2Gap = begin(getValue(me.contigStore, contigId).gaps);
		TContigGapsIter it2GapEnd = end(getValue(me.contigStore, contigId).gaps);
		goNext(it2Gap);
		TInterval contigInterval;
		TInterval readInterval;
		TContigPos lastGapPosBeforeGap1;
		TContigPos lastGapPosBeforeGap2;
		for ( ; itI != itIEnd; goNext(itI))
		{
			readInterval = getValue(itI);

			if (static_cast<TContigPos>(readInterval.i2) < getValue(it1Gap).gapPos) continue;    		// interval-positions are smaller than first-position of contig: read-interval mappt in contig gaps
		
			// interval i1:
			while ( (it2Gap != it2GapEnd) && (static_cast<TContigPos>(readInterval.i1) >= getValue(it2Gap).gapPos) )   // iterate over 2 contig gap-anchors, until gapPos2 is behind start of read-interval
			{
				goNext(it1Gap);
				goNext(it2Gap);
			}
		
			if (it2Gap != it2GapEnd)
			{
				lastGapPosBeforeGap1 = getValue(it2Gap).seqPos + getValue(it1Gap).gapPos - getValue(it1Gap).seqPos -1;  // last position in gapped sequence before gap 
			
				if (static_cast<TContigPos>(readInterval.i1) < getValue(it1Gap).gapPos)  				// occurs, if interval.i1 position is smaller than start-postion of contig
				{
					contigInterval.i1 = getValue(it1Gap).gapPos;
				}
				else if (static_cast<TContigPos>(readInterval.i1) <= lastGapPosBeforeGap1)				// read-interval starts in contig-interval between the 2 gaps		
				{
					contigInterval.i1 = readInterval.i1 - (getValue(it1Gap).gapPos - getValue(it1Gap).seqPos);  // project position onto ungapped sequence
				}
				else if (static_cast<TContigPos>(readInterval.i2) < getValue(it2Gap).gapPos)				// whole read-interval lies in gaps of contig
				{
					continue;
				}
				else 										// read-interval starts in contig gap -> alignment-interval starts in seqPos2
				{
					contigInterval.i1 = getValue(it2Gap).seqPos;
				}
			}
			else		// read-interval lies behind last gap-anchor
			{
				goPrevious(it1Gap);
				goPrevious(it2Gap);
				contigInterval.i1 = readInterval.i1 - (getValue(it2Gap).gapPos - getValue(it2Gap).seqPos); // not goPrev -> it1Gap
			}
		
			// interval i2:	
			while ( (it2Gap != it2GapEnd) && (static_cast<TContigPos>(readInterval.i2) >= getValue(it2Gap).gapPos) )		// iterate over 2 contig gap-anchors, until gapPos2 is behind end of read-interval
			{
				goNext(it1Gap);
				goNext(it2Gap);
			}
	
			if (it2Gap != it2GapEnd)
			{
				lastGapPosBeforeGap2 = getValue(it2Gap).seqPos + getValue(it1Gap).gapPos - getValue(it1Gap).seqPos -1;  // last position in gapped sequence before gap 
				
				if (static_cast<TContigPos>(readInterval.i2) <= lastGapPosBeforeGap2)					// read-interval ends in contig-interval between the 2 gaps	
				{
					contigInterval.i2 = readInterval.i2 - (getValue(it1Gap).gapPos - getValue(it1Gap).seqPos);    // project position onto ungapped sequence
				}
				else 										// read-interval ends in gap -> end of alignment-interval is last position before gap
				{
					contigInterval.i2 = lastGapPosBeforeGap2 - (getValue(it1Gap).gapPos - getValue(it1Gap).seqPos);
				}
			}
			else if (it2Gap == it2GapEnd)				// read-intervall ends behind the last gap-anchor
			{
				goPrevious(it1Gap);
				goPrevious(it2Gap);
				contigInterval.i2 = readInterval.i2 - (getValue(it2Gap).gapPos - getValue(it2Gap).seqPos);
			}	
		
			appendValue(contigIntervals, contigInterval, Generous());
		}
	}
	mergeIntervals(contigIntervals, thresholdGaps);
	
}


////////////////////////////////////////////////////////////////////////////////////////////////////////
//// merge intervals (for sortet intervals; i1 <= i2 <= (i+1)1 <= (i+1)2 )
////////////////////////////////////////////////////////////////////////////////////////////////////////
template<typename TIntervals>
inline void
mergeIntervals(TIntervals & intervals, const unsigned & thresholdGaps)
{
	typedef typename Position<TIntervals>::Type 	TPos;
	typedef typename Value<TIntervals>::Type 	TInterval;
	
	TPos j;
	TInterval newInterval;
	for (TPos i = 0; i < length(intervals) - 1; ++i)
	{
		j = i;
		while ( (j < length(intervals) - 1) && (getValue(intervals, j).i2 + thresholdGaps >= getValue(intervals, j + 1).i1 - 1) )		// merges intervals, if the no. of gaps inbetween is smaller than threshold 
		{
			++j;
		}
		
		if ( j != i)
		{
			newInterval.i1 = getValue(intervals, i).i1;
			newInterval.i2 = getValue(intervals, j).i2;
			replace(intervals, i, j + 1, newInterval);       
		}
		
	}
}


////////////////////////////////////////////////////////////////////////////////////////////////////////

}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
