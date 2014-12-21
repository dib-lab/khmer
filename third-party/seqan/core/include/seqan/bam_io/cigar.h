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

#ifndef CORE_INCLUDE_SEQAN_BAM_IO_CIGAR_H_
#define CORE_INCLUDE_SEQAN_BAM_IO_CIGAR_H_

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// ----------------------------------------------------------------------------
// struct CigarElement
// ----------------------------------------------------------------------------

/**
.Class.CigarElement
..cat:Fragment Store
..summary:One entry of a CIGAR string.
..signature:CigarElement<TOperation, TCount>
..param.TOperation:Type to use for storing operations.
...default:nolink:$char$
..param.TCount:Type to use for storing counts.
...default:nolink:$unsigned$
..include:seqan/store.h

.Memfunc.CigarElement#CigarElement
..class:Class.CigarElement
..summary:Constructor
..signature:CigarElement()
..signature:CigarElement(operation, count)
..param.operation:The operation to use.
...type:nolink:$TOperation$, typically $char$.
..param.count:The operation count.
...type:nolink:$Count$, typically $unsigned$.
..remarks:The default constructor initialized both @Memvar.CigarElement#operation@ and @Memvar.CigarElement#count@ with $0$.

.Memvar.CigarElement#operation
..class:Class.CigarElement
..summary:The described operation.
..type:nolink:$TOperation$

.Memvar.CigarElement#count
..class:Class.CigarElement
..summary:The number of operations.
..type:nolink:$TCount$
*/

template <typename TOperation_ = char, typename TCount_ = unsigned>
struct CigarElement
{
    typedef TOperation_ TOperation;
    typedef TCount_     TCount;

    TOperation          operation;
    TCount              count;

    CigarElement() : operation(0), count(0) {}
    
    CigarElement(TOperation o, TCount c):
        operation(o),
        count(c) {}
};

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

template <typename TOperation, typename TCount>
inline bool operator>(CigarElement<TOperation, TCount> const & lhs,
                      CigarElement<TOperation, TCount> const & rhs)
{
    return lhs.operation > rhs.operation || (lhs.operation == rhs.operation && lhs.count > rhs.count);
}

template <typename TOperation, typename TCount>
inline bool operator<(CigarElement<TOperation, TCount> const & lhs,
                      CigarElement<TOperation, TCount> const & rhs)
{
    return lhs.operation < rhs.operation || (lhs.operation == rhs.operation && lhs.count < rhs.count);
}
    
template <typename TOperation, typename TCount>
inline bool operator==(CigarElement<TOperation, TCount> const & lhs,
                       CigarElement<TOperation, TCount> const & rhs)
{
    return lhs.operation == rhs.operation && lhs.count == rhs.count;
}

// ----------------------------------------------------------------------------
// toBamCigarElement()
// ----------------------------------------------------------------------------

template <typename TOperation, typename TCount>
__uint32 toBamCigarElement(CigarElement<TOperation, TCount> const & cigarElement)
{
    char operation = 0;
    switch (cigarElement.operation) {
        case 'X': operation += 1;
        case '=': operation += 1;
        case 'P': operation += 1;
        case 'H': operation += 1;
        case 'S': operation += 1;
        case 'N': operation += 1;
        case 'D': operation += 1;
        case 'I': operation += 1;
        case 'M': break;
    }
    return (cigarElement.count << 4) | operation;
}

// ----------------------------------------------------------------------------
// getMDString()
// ----------------------------------------------------------------------------

template <
    typename TMDString,
    typename TGaps1,
    typename TGaps2>
inline void
getMDString(
    TMDString &md,
    TGaps1 &gaps1,
    TGaps2 &gaps2)
{
    typedef typename Value<TMDString>::Type TMDChar;
	typename Iterator<TGaps1>::Type it1 = begin(gaps1);
	typename Iterator<TGaps2>::Type it2 = begin(gaps2);
	char op, lastOp = ' ';
	unsigned numOps = 0;

    clear(md);
    for (; !atEnd(it1) && !atEnd(it2); goNext(it1), goNext(it2))
    {
        if (isGap(it1)) continue;
        if (isGap(it2))
        {
            op = 'D';
        } 
        else
            op = (*it1 == *it2)? 'M': 'R';
        
        // append match run
        if (lastOp != op)
        {
            if (lastOp == 'M')
            {
                std::stringstream num;
                num << numOps;
                append(md, num.str());
            }
            numOps = 0;
            lastOp = op;
        }

        // append deleted/replaced reference character
        if (op != 'M')
        {
            // add ^ from non-deletion to deletion
            if (op == 'D' && lastOp != 'D')
                appendValue(md, '^');
            // add 0 from deletion to replacement
            if (op == 'R' && lastOp == 'D')
                appendValue(md, '0');
            appendValue(md, convert<TMDChar>(*it1));
        }

        ++numOps;
    }
    SEQAN_ASSERT_EQ(atEnd(it1), atEnd(it2));
    if (lastOp == 'M')
    {
        std::stringstream num;
        num << numOps;
        append(md, num.str());
    }
}

// ----------------------------------------------------------------------------
// getCigarString()
// ----------------------------------------------------------------------------

template <
    typename TCigar,
    typename TGaps1,
    typename TGaps2,
    typename TThresh>
inline void
getCigarString(
    TCigar &cigar,
    TGaps1 &gaps1,
    TGaps2 &gaps2,
    TThresh splicedGapThresh)
{
	typename Iterator<TGaps1>::Type it1 = begin(gaps1);
	typename Iterator<TGaps2>::Type it2 = begin(gaps2);
	clear(cigar);
	char op, lastOp = ' ';
	unsigned numOps = 0;

	// std::cout << "gaps1\t" << gaps1 << std::endl;
	// std::cout << "gaps2\t" << gaps2 << "\t" << clippedBeginPosition(gaps2) << std::endl;
	for (; !atEnd(it1) && !atEnd(it2); goNext(it1), goNext(it2))
	{
		if (isGap(it1))
		{
			if (isGap(it2))
				op = 'P';
			else if (isClipped(it2))
				op = '?';
			else
				op = 'I';
		} 
		else if (isClipped(it1))
		{
			op = '?';
		}
		else 
		{
			if (isGap(it2))
				op = 'D';
			else if (isClipped(it2))
				op = 'S';
			else
				op = 'M';
		}
        
        // append CIGAR operation
        if (lastOp != op)
        {
            if (lastOp == 'D' && numOps >= (unsigned)splicedGapThresh)
                lastOp = 'N';
            if (numOps > 0)
            {
                std::stringstream num;
                num << numOps;
                append(cigar, num.str());
                appendValue(cigar, lastOp);
            }
            numOps = 0;
            lastOp = op;
        }
        ++numOps;
    }
//  if (atEnd(it1) != atEnd(it2))
//        std::cerr << "Invalid pairwise alignment:" << std::endl << gaps1 << std::endl << gaps2 << std::endl;
	SEQAN_CHECK(atEnd(it1) == atEnd(it2), "Cannot get CIGAR from invalid pairwise alignment!");
	if (lastOp == 'D' && numOps >= (unsigned)splicedGapThresh)
		lastOp = 'N';
	if (numOps > 0)
	{
		std::stringstream num;
		num << numOps;
		append(cigar, num.str());
		appendValue(cigar, lastOp);
	}
}

template <
    typename TCigar,
    typename TGaps1,
    typename TGaps2>
inline void
getCigarString(
    TCigar &cigar,
    TGaps1 &gaps1,
    TGaps2 &gaps2)
{
    return getCigarString(cigar, gaps1, gaps2, 20);
}

template <
    typename TOperation,
    typename TCount,
    typename TSpec,
    typename TGaps1,
    typename TGaps2,
    typename TThresh>
inline void
getCigarString(
        String<CigarElement<TOperation, TCount>, TSpec> &cigar,
        TGaps1 &gaps1,
        TGaps2 &gaps2,
        TThresh splicedGapThresh)
{
	typename Iterator<TGaps1>::Type it1 = begin(gaps1);
	typename Iterator<TGaps2>::Type it2 = begin(gaps2);
	clear(cigar);
	char op = '?', lastOp = ' ';
	unsigned numOps = 0;

//  std::cout << gaps1 << std::endl;
//  std::cout << gaps2 << std::endl;
    for (; !atEnd(it1) && !atEnd(it2); goNext(it1), goNext(it2))
    {
        if (isGap(it1))
        {
            if (isGap(it2))
                op = 'P';
            else if (isClipped(it2))
                op = '?';
            else
                op = 'I';
        } 
        else if (isClipped(it1))
        {
            op = '?';
        }
        else 
        {
            if (isGap(it2))
                op = 'D';
            else if (isClipped(it2))
                op = 'S';
            else
                op = 'M';
        }
        if (lastOp != op)
        {
            if (lastOp == 'D' && numOps >= (unsigned)splicedGapThresh)
                lastOp = 'N';
            if (numOps > 0)
                appendValue(cigar, CigarElement<>(lastOp, numOps));
            numOps = 0;
            lastOp = op;
        }
        ++numOps;
    }
    SEQAN_ASSERT_EQ(atEnd(it1), atEnd(it2));
    if (lastOp == 'D' && numOps >= (unsigned)splicedGapThresh)
        lastOp = 'N';
    if (numOps > 0)
        appendValue(cigar, CigarElement<>(op, numOps));
}

// ----------------------------------------------------------------------------
// alignAndGetCigarString()
// ----------------------------------------------------------------------------

template <
    typename TCigar, typename TMDString, typename TContig, typename TReadSeq,
    typename TAlignedRead, typename TErrors, typename TAlignFunctor>
inline void
alignAndGetCigarString(
    TCigar &cigar, TMDString &md, TContig &contig, TReadSeq &readSeq,
    TAlignedRead &alignedRead, TErrors &errors, TAlignFunctor const & functor)
{
    typedef Align<TReadSeq, ArrayGaps> TAlign;

    TAlign align;
    resize(rows(align), 2);

    if (alignedRead.beginPos <= alignedRead.endPos)
        assignSource(row(align, 0), infix(contig.seq, alignedRead.beginPos, alignedRead.endPos));
    else
        assignSource(row(align, 0), infix(contig.seq, alignedRead.endPos, alignedRead.beginPos));

    assignSource(row(align, 1), readSeq);
    
    if (!(errors == 0 || (errors == 1 && length(readSeq) == length(source(row(align, 0))))))
        errors = functor.align(align);

    getCigarString(cigar, row(align, 0), row(align, 1));
    getMDString(md, row(align, 0), row(align, 1));
}

template <typename TCigar, typename TMDString, typename TContig, typename TReadSeq, typename TErrors, typename TAlignedRead>
inline void
alignAndGetCigarString(TCigar &cigar, TMDString &md, TContig &contig, TReadSeq &readSeq, TAlignedRead &alignedRead, TErrors &, Nothing const &)
{
    typedef typename TContig::TContigSeq                                    TContigSeq;
    typedef Gaps<TContigSeq, AnchorGaps<typename TContig::TGapAnchors> >    TContigGaps;
    typedef Gaps<TReadSeq, AnchorGaps<typename TAlignedRead::TGapAnchors> > TReadGaps;

    TContigGaps contigGaps(contig.seq, contig.gaps);
    
    if (alignedRead.beginPos <= alignedRead.endPos) 
    {
        setClippedBeginPosition(contigGaps, alignedRead.beginPos);
        setClippedEndPosition(contigGaps, alignedRead.endPos);
    } else
    {
        setClippedBeginPosition(contigGaps, alignedRead.endPos);
        setClippedEndPosition(contigGaps, alignedRead.beginPos);
    }

    TReadGaps readGaps(readSeq, alignedRead.gaps);
    // TContigGaps  contigGaps2(contig.seq, contig.gaps);
    // if (i == 4)
    //     printf("It's it!\n");
    // std::cerr << "read gaps:  " << readGaps << std::endl;
    // std::cerr << "contig gaps:" << contigGaps << std::endl;
    
    getCigarString(cigar, contigGaps, readGaps);
    getMDString(md, contigGaps, readGaps);
}

// ----------------------------------------------------------------------------
// _getClippedLength()
// ----------------------------------------------------------------------------

template <typename TCigarString, typename TNum>
inline void _getClippedLength(TCigarString const & cigar, TNum & sum)
{
    typedef typename Iterator<TCigarString, Standard>::Type TCigarIter;
    
    TCigarIter it = begin(cigar, Standard());
    TCigarIter itEnd = end(cigar, Standard());
    
    sum = 0;
    for (; it != itEnd; ++it)
        if (getValue(it).operation != 'S' && getValue(it).operation != 'H')
            sum += getValue(it).count;
}

// ----------------------------------------------------------------------------
// _getLengthInRef()
// ----------------------------------------------------------------------------

template <typename TCigarString, typename TNum>
inline void _getLengthInRef(TCigarString const & cigar, TNum & sum)
{
    typedef typename Iterator<TCigarString const, Standard>::Type TCigarIter;
    
    TCigarIter it = begin(cigar, Standard());
    TCigarIter itEnd = end(cigar, Standard());
    
    sum = 0;
    for (; it != itEnd; ++it)
        if (getValue(it).operation != 'S' && getValue(it).operation != 'H' && getValue(it).operation != 'I')
            sum += getValue(it).count;
}

// ----------------------------------------------------------------------------
// cigarToGapAnchorRead()
// ----------------------------------------------------------------------------

template<typename TCigarString, typename TGaps>
inline unsigned
cigarToGapAnchorRead(TCigarString const & cigar, TGaps & gaps)
{
    typename Iterator<TGaps>::Type it = begin(gaps);
    bool atBegin = true;
    unsigned beginGaps = 0;
    for (unsigned i = 0; i < length(cigar); ++i)
    {
        switch (cigar[i].operation)
        {
            case 'D':
            case 'N':
            case 'P':
                if (atBegin)
                    beginGaps += cigar[i].count;
                insertGaps(it, cigar[i].count);
            case 'I':
            case 'M':
            case 'S':
                it += cigar[i].count;
                atBegin = false;
        }
    }
    return beginGaps;
}

// ----------------------------------------------------------------------------
// cigarToGapAnchorContig()
// ----------------------------------------------------------------------------

template<typename TCigarString, typename TGaps>
inline unsigned
cigarToGapAnchorContig(TCigarString const & cigar, TGaps & gaps)
{
    typename Iterator<TGaps>::Type it = begin(gaps);
    bool atBegin = true;
    unsigned beginGaps = 0;
    for (unsigned i = 0; i < length(cigar); ++i)
    {
        switch (cigar[i].operation)
        {
            case 'I':
            case 'P':
                if (atBegin)
                    beginGaps += cigar[i].count;
                insertGaps(it, cigar[i].count);
            case 'D':
            case 'M':
            case 'N':
            case 'S':
                it += cigar[i].count;
                atBegin = false;
        }
    }
    return beginGaps;
}

}  // namespace seqan

#endif  // #ifndef CORE_INCLUDE_SEQAN_BAM_IO_CIGAR_H_
