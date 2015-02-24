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
// Author: David Weese <david.weese@fu-berlin.de>
// ==========================================================================

#ifndef INCLUDE_SEQAN_BAM_IO_CIGAR_H_
#define INCLUDE_SEQAN_BAM_IO_CIGAR_H_

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

/*!
 * @class CigarElement
 * @headerfile <seqan/bam_io.h>
 * @brief One entry of a CIGAR string.
 *
 * @signature template <[typename TOperation[, typename TCount]]>
 *            class CigarElement;
 *
 * @tparam TOperation Type to use for storing operations, defaults to <tt>char</tt>.
 * @tparam TCount Type to use for storing counts, defaults to <tt>unsigned</tt>.
 */

/*!
 * @fn CigarElement::CigarElement
 * @brief Constructor
 *
 * @signature CigarElement::CigarElement();
 * @signature CigarElement::CigarElement(operation, count);
 *
 * @param[in] operation The operation to use, of type <tt>TOperation</tt>.
 * @param[in] count     The operation count, of type <tt>TCount</tt>.
 *
 * @section Remarks
 *
 * The default constructor initialized both @link CigarElement::operation @endlink and @link CigarElement::count
 * @endlink with <tt>0</tt>.
 */

/*!
 * @var TCount CigarElement::count;
 *
 * @brief The number of operations.
 */

/*!
 * @var TOperation CigarElement::operation;
 *
 * @brief The described operation.
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

template <typename TOperation, typename TCount>
struct Size<CigarElement<TOperation, TCount> >
{
    typedef TCount Type;
};

// ============================================================================
// Functions
// ============================================================================

template <typename TOperation, typename TCount>
inline bool operator>(CigarElement<TOperation, TCount> const & lhs,
                      CigarElement<TOperation, TCount> const & rhs)
{
    return lhs.operation > rhs.operation || (lhs.operation == rhs.operation && (lhs.count) > (rhs.count));
}

template <typename TOperation, typename TCount>
inline bool operator<(CigarElement<TOperation, TCount> const & lhs,
                      CigarElement<TOperation, TCount> const & rhs)
{
    return lhs.operation < rhs.operation || (lhs.operation == rhs.operation && (lhs.count) < (rhs.count));
}

template <typename TOperation, typename TCount>
inline bool operator==(CigarElement<TOperation, TCount> const & lhs,
                       CigarElement<TOperation, TCount> const & rhs)
{
    return lhs.operation == rhs.operation && (lhs.count) == (rhs.count);
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
inline unsigned
getMDString(
    TMDString &md,
    TGaps1 &gaps1,  // typically reference
    TGaps2 &gaps2)  // typically read
{
    typedef typename Value<TMDString>::Type TMDChar;
    typedef typename Value<typename Host<TGaps1>::Type>::Type TVal1;
    typedef typename Value<typename Host<TGaps2>::Type>::Type TVal2;

    typename Iterator<TGaps1>::Type it1 = begin(gaps1);
    typename Iterator<TGaps2>::Type it2 = begin(gaps2);
    char op, lastOp = ' ';
    unsigned numOps = 0;
    unsigned errors = 0;

    clear(md);
    for (; !atEnd(it1) && !atEnd(it2); goNext(it1), goNext(it2))
    {
        if (isGap(it1))
        {
            if (!isGap(it2))
                ++errors;
            continue;       // insertion to the reference (gaps1)
//            op = 'I';     // ignore insertions completely
        }
        if (isGap(it2))
        {
            ++errors;
//            if (op == 'I')  // ignore paddings
//                continue;
            op = 'D';       // deletion from the reference (gaps1)
        }
        else
        {
            if ((TVal1)*it1 == (TVal2)*it2)
            {
                op = 'M';
            }
            else
            {
                op = 'R';
                ++errors;
            }
        }

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
        }

        // append deleted/replaced reference character
        if (op != 'M')
        {
            // add ^ for deleted reference bases (from non-deletion to deletion)
            if (op == 'D' && lastOp != 'D')
                appendValue(md, '^');
            // add 0 for each replaced base that doesn't follow a match (for samtools/BWA compatibility)
            else if (op == 'R' && lastOp != 'M')
                appendValue(md, '0');
            appendValue(md, convert<TMDChar>(*it1));
        }

        lastOp = op;
        ++numOps;
    }
    SEQAN_ASSERT_EQ(atEnd(it1), atEnd(it2));
    if (lastOp == 'M')
    {
        std::stringstream num;
        num << numOps;
        append(md, num.str());
    }
    return errors;
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
    TGaps1 &gaps1,  // typically reference
    TGaps2 &gaps2,  // typically read
    TThresh splicedGapThresh)
{
    typename Iterator<TGaps1>::Type it1 = begin(gaps1);
    typename Iterator<TGaps2>::Type it2 = begin(gaps2);
//    typedef typename Value<typename Host<TGaps1>::Type>::Type TVal1;
//    typedef typename Value<typename Host<TGaps2>::Type>::Type TVal2;

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
//                op = ((TVal1)*it1 == (TVal2)*it2)? '=': 'X';
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
    TGaps1 &gaps1,  // typically reference
    TGaps2 &gaps2)  // typically read
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
//    typedef typename Value<typename Host<TGaps1>::Type>::Type TVal1;
//    typedef typename Value<typename Host<TGaps2>::Type>::Type TVal2;

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
//                op = ((TVal1)*it1 == (TVal2)*it2)? '=': 'X';
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
    typename TAlignedRead, typename TErrors >
inline void
alignAndGetCigarString(
    TCigar &cigar, TMDString &md, TContig &, TReadSeq &,
    TAlignedRead &, TErrors &, Nothing const &)
{
    cigar = "*";
    clear(md);
}

struct BamAlignFunctorEditDistance
{
    typedef String<GapAnchor<int> > TGapAnchors;

    TGapAnchors contigAnchors, readAnchors;

    template <typename TGaps1, typename TGaps2, typename TErrors>
    inline int
    align(TGaps1 &gaps1, TGaps2 &gaps2, TErrors maxErrors)
    {
        return -globalAlignment(
            gaps1, gaps2,
            Score<short, EditDistance>(),
            -(int)maxErrors, (int)maxErrors
        );
    }
};

struct BamAlignFunctorSemiGlobalGotoh
{
    typedef String<GapAnchor<int> > TGapAnchors;

    Score<int> score;
    TGapAnchors contigAnchors, readAnchors;

    BamAlignFunctorSemiGlobalGotoh(Score<int> score_) :
        score(score_)
    {}

    template <typename TGaps1, typename TGaps2, typename TErrors>
    inline int
    align(TGaps1 &gaps1, TGaps2 &gaps2, TErrors maxErrors)
    {
        return globalAlignment(
            gaps1, gaps2, score,
            AlignConfig<true, false, false, true>(),
            -(int)maxErrors, (int)maxErrors,
            Gotoh()
        ) / scoreMismatch(score);
    }
};

struct BamAlignFunctorDefault
{
};

template <
    typename TCigar, typename TMDString, typename TContigInfix, typename TReadSeq,
    typename TAlignedRead, typename TErrors, typename TAlignFunctor>
inline void
_alignAndGetCigarString(
    TCigar &cigar, TMDString &md, TContigInfix const &contigInfix, TReadSeq const &fwdReadSeq,
    TAlignedRead &, TErrors &errors, TAlignFunctor &functor)
{
    typedef Gaps<TContigInfix, AnchorGaps<typename TAlignFunctor::TGapAnchors> >    TContigGaps;
    typedef Gaps<TReadSeq, AnchorGaps<typename TAlignFunctor::TGapAnchors> >        TReadGaps;

    clear(functor.contigAnchors);
    clear(functor.readAnchors);

    TContigGaps contigGaps(contigInfix, functor.contigAnchors);
    TReadGaps readGaps(fwdReadSeq, functor.readAnchors);

    // if there is already an alignment between contigInfix and fwdReadSeq with 0 or 1 error then
    // we don't to realign as it contains no gaps
    if (!(errors == 0 || (errors == 1 && length(contigInfix) == length(fwdReadSeq))))
        errors = functor.align(contigGaps, readGaps, errors);

    getCigarString(cigar, contigGaps, readGaps);
    TErrors mdErrors = getMDString(md, contigGaps, readGaps);

    ignoreUnusedVariableWarning(mdErrors);
    SEQAN_ASSERT_EQ(errors, mdErrors);
}

template <
    typename TCigar, typename TMDString, typename TContig, typename TReadSeq,
    typename TAlignedRead, typename TErrors, typename TAlignFunctor>
inline void
alignAndGetCigarString(
    TCigar &cigar, TMDString &md, TContig const &contig, TReadSeq const &readSeq,
    TAlignedRead &alignedRead, TErrors &errors, TAlignFunctor &functor)
{
    typedef typename TContig::TContigSeq            TContigSeq;
    typedef typename Infix<TContigSeq const>::Type  TContigInfix;

    TContigInfix contigInfix;

    if (alignedRead.beginPos <= alignedRead.endPos)
    {
        contigInfix = infix(contig.seq, alignedRead.beginPos, alignedRead.endPos);
        _alignAndGetCigarString(cigar, md, contigInfix, readSeq, alignedRead, errors, functor);
    }
    else
    {
        contigInfix = infix(contig.seq, alignedRead.endPos, alignedRead.beginPos);
        _alignAndGetCigarString(cigar, md, contigInfix, reverseComplementString(readSeq), alignedRead, errors, functor);
    }
}

template <
    typename TCigar, typename TMDString, typename TContig, typename TReadSeq,
    typename TAlignedRead, typename TErrors>
inline void
alignAndGetCigarString(
    TCigar &cigar, TMDString &md, TContig const &contig, TReadSeq const &readSeq,
    TAlignedRead &alignedRead, TErrors &errors, BamAlignFunctorDefault &)
{
    typedef typename TContig::TContigSeq                                            TContigSeq;
    typedef Gaps<TContigSeq, AnchorGaps<typename TContig::TGapAnchors> >            TContigGaps;
    typedef typename ReverseComplementString<TReadSeq const>::Type                  TRefCompReadSeq;
    typedef Gaps<TReadSeq const, AnchorGaps<typename TAlignedRead::TGapAnchors> >   TReadGaps;
    typedef Gaps<TRefCompReadSeq, AnchorGaps<typename TAlignedRead::TGapAnchors> >  TRCReadGaps;

    TContigGaps contigGaps(contig.seq, contig.gaps);

    if (alignedRead.beginPos <= alignedRead.endPos)
    {
        setClippedBeginPosition(contigGaps, alignedRead.beginPos);
        setClippedEndPosition(contigGaps, alignedRead.endPos);

        TReadGaps readGaps(readSeq, alignedRead.gaps);

        getCigarString(cigar, contigGaps, readGaps);
        errors = getMDString(md, contigGaps, readGaps);
    }
    else
    {
        setClippedBeginPosition(contigGaps, alignedRead.endPos);
        setClippedEndPosition(contigGaps, alignedRead.beginPos);

        TRCReadGaps readGaps(reverseComplementString(readSeq), alignedRead.gaps);

        getCigarString(cigar, contigGaps, readGaps);
        errors = getMDString(md, contigGaps, readGaps);
    }
}

// ----------------------------------------------------------------------------
// _getClippedLength()
// ----------------------------------------------------------------------------

template <typename TCigarString, typename TNum>
inline void _getClippedLength(TNum & sum, TCigarString const & cigar)
{
    typedef typename Iterator<TCigarString const, Standard>::Type TCigarIter;

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
inline void _getLengthInRef(TNum & sum, TCigarString const & cigar)
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
// _getQueryLength()
// ----------------------------------------------------------------------------

template <typename TCigarString>
inline typename Size<typename Value<TCigarString>::Type>::Type
_getQueryLength(TCigarString const & cigar)
{
    typedef typename Iterator<TCigarString const, Standard>::Type TCigarIter;
    typedef typename Size<typename Value<TCigarString>::Type>::Type TSize;
    TCigarIter it = begin(cigar, Standard());
    TCigarIter itEnd = end(cigar, Standard());

    TSize len = 0;
    for (; it != itEnd; ++it)
        if (getValue(it).operation != 'D' && getValue(it).operation != 'H' && getValue(it).operation != 'N' && getValue(it).operation != 'P')
            len += getValue(it).count;
    return len;
}

// ----------------------------------------------------------------------------
// cigarToGapAnchorRead()
// ----------------------------------------------------------------------------

template <typename TGaps, typename TCigarString>
unsigned cigarToGapAnchorRead(TGaps & gaps, TCigarString const & cigar)
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
unsigned cigarToGapAnchorContig(TGaps & gaps, TCigarString const & cigar)
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

#endif  // #ifndef INCLUDE_SEQAN_BAM_IO_CIGAR_H_
