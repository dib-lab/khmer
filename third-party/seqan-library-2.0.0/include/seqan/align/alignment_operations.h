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
// Author: Birte Kehr <birte.kehr@fu-berlin.de>
// ==========================================================================
// Operations on alignments such as integration
// ==========================================================================

#ifndef SEQAN_INCLUDE_SEQAN_ALIGN_ALIGNMENT_OPERATIONS_H_
#define SEQAN_INCLUDE_SEQAN_ALIGN_ALIGNMENT_OPERATIONS_H_

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function integrateAlign()
// ----------------------------------------------------------------------------

/*!
 * @fn integrateAlign
 * @headerfile <seqan/align.h>
 * @brief Integrates an alignment into another by copying the gaps.
 *
 * @signature void integrateAlign(align1, align2[, positions]);
 *
 * @param[in,out] align1    Target Alignment object into which align2 is to be integrated.
 * @param[in]     align2    Alignment object that is to be integrated into align1.
 * @param[in]     positions The integration positions in align1 for all rows (view positions), String of positions.
 *
 * @section Examples
 *
 * @include demos/align/integrate_align.cpp
 *
 * The output is as follows:
 *
 * @include demos/align/integrate_align.cpp.stdout
 */

template <typename TSource1, typename TSpec1, typename TSource2, typename TSpec2, typename TPos>
void integrateAlign(Align<TSource1, TSpec1> & align,
                    Align<TSource2, TSpec2> const & infixAlign,
                    String<TPos> const & viewPos)
{
    typedef Align<TSource1, TSpec1> TAlign;
    typedef Align<TSource2, TSpec2> TInfixAlign;
    typedef typename Size<TAlign>::Type TSize;

    typedef typename Row<TAlign>::Type TRow;
    typedef typename Row<TInfixAlign const>::Type TInfixRow;

    // Iterators on align and infixAlign.
    typename Iterator<TRow>::Type it;
    typedef typename Iterator<TInfixRow, Standard>::Type TInfixRowIt;


    for (TSize i = 0; i < length(rows(align)); ++i)
    {
        TInfixRow const & infixRow = row(infixAlign, i);
        // This assertion ensures that the number of sequence characters after viewPos[i] is greater than or equal to
        // the number of source characters in the clipped infix row.
        SEQAN_ASSERT_GEQ(endPosition(row(align, i)) - toSourcePosition(row(align, i), viewPos[i]),
                endPosition(infixRow) - beginPosition(infixRow));

        // init iterators
        it = iter(row(align, i), value(viewPos, i));

        // walk through Gaps containers and copy gaps
        for (TInfixRowIt infixIt = begin(infixRow, Standard()); !atEnd(infixIt);)
        {
            TSize gapSize = countGaps(infixIt);
            insertGaps(it, gapSize);
            goFurther(it, gapSize+1);
            goFurther(infixIt, gapSize+1);
        }
    }
}

template <typename TSource, typename TSpec1, typename TSpec2>
void integrateAlign(Align<TSource, TSpec1> & align,
                    Align<typename Infix<TSource>::Type, TSpec2> const & infixAlign)
{
    typedef typename Size<TSource>::Type TSize;
    typedef typename Position<typename Row<Align<TSource, TSpec1> >::Type>::Type TPos;

    String<TPos> viewPos;
    TPos pos;
    for (TSize i = 0; i < length(rows(infixAlign)); ++i)
    {
        pos = beginPosition(source(row(infixAlign, i)))
            - beginPosition(source(row(align, i)))
            + beginPosition(row(infixAlign, i));
        appendValue(viewPos, toViewPosition(row(align, i), pos));
    }

    integrateAlign(align, infixAlign, viewPos);
}

}  // namespace seqan

#endif  // #ifndef SEQAN_INCLUDE_SEQAN_ALIGN_ALIGNMENT_OPERATIONS_H_
