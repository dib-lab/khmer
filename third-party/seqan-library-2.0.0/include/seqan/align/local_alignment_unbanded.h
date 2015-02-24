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
// Author: Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>
// ==========================================================================
// Interface functions for unbanded local alignment.
// ==========================================================================

#ifndef SEQAN_INCLUDE_SEQAN_ALIGN_LOCAL_ALIGNMENT_UNBANDED_H_
#define SEQAN_INCLUDE_SEQAN_ALIGN_LOCAL_ALIGNMENT_UNBANDED_H_

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
// Function localAlignment()
// ----------------------------------------------------------------------------

/*!
 * @fn localAlignment
 * @headerfile <seqan/align.h>
 * @brief Computes the best pairwise local alignment using the Smith-Waterman algorithm.
 *
 * @signature TScoreVal localAlignment(align,          scoringScheme, [lowerDiag, upperDiag]);
 * @signature TScoreVal localAlignment(gapsH, gapsV,   scoringScheme, [lowerDiag, upperDiag]);
 * @signature TScoreVal localAlignment(fragmentString, scoringScheme, [lowerDiag, upperDiag]);
 *
 * @param[in,out] gapsH Horizontal gapped sequence in alignment matrix. Types: @link Gaps @endlink
 * @param[in,out] gapsV Vertical gapped sequence in alignment matrix. Types: @link Gaps @endlink
 * @param[in,out] align An @link Align @endlink object that stores the alignment. The
 *                      number of rows must be 2 and the sequences must have already
 *                      been set. <tt>align[0]</tt> is the horizontal one in the
 *                      alignment matrix alignment, <tt>align[1]</tt> is the vertical
 *                      one.
 * @param[in,out] fragmentString
 *                      String of @link Fragment @endlink objects. The sequence
 *                      with id <tt>0</tt> is the horizontal one, the sequence
 *                      with id <tt>1</tt> is the vertical one.
 * @param[in] scoringScheme
 *                      The @link Score scoring scheme @endlink to use for the alignment.
 * @param[in] lowerDiag Optional lower diagonal (<tt>int</tt>).
 * @param[in] upperDiag Optional upper diagonal (<tt>int</tt>).
 *
 * @return TScoreVal Score value of the resulting alignment  (Metafunction @link Score#Value @endlink of the type of
 *                   <tt>scoringScheme</tt>).
 *
 * The Waterman-Eggert algorithm (local alignment with declumping) is available through the @link
 * LocalAlignmentEnumerator @endlink class.
 *
 * When using @link Gaps @endlink and @link Align @endlink objects, only parts (i.e. one infix) of each sequence will be
 * aligned.  This will be presented to the user by setting the clipping begin and end position of the gaps (the rows in
 * the case of @link Align @endlink objects).  When using @link Fragment @endlink strings, these parts of the sequences
 * will not appear in any fragment.
 *
 * There exist multiple overloads for this function with two configuration dimensions.
 *
 * First, you can select the type of the target storing the alignment. This can be either an @link Align @endlink
 * object, two @link Gaps @endlink objects, or a string of @link Fragment @endlink objects. @link Align @endlink objects
 * provide an interface to tabular alignments with the restriction of all rows having the same type. Using two @link
 * Gaps @endlink objects has the advantage that you an align sequences with different types, for example @link DnaString
 * @endlink and @link Dna5String @endlink. Using @link Fragment @endlink strings is useful for collecting many pairwise
 * alignments, for example in the construction of @link AlignmentGraph Alignment Graphs @endlink for multiple- sequence
 * alignments (MSA).
 *
 * Second, you can optionally give a band for the alignment using <tt>lowerDiag</tt> and <tt>upperDiag</tt>. The center
 * diagonal has index <tt>0</tt>, the <tt>i</tt>th diagonal below has index <tt>-i</tt>, the <tt>i</tt>th above has
 * index <tt>i</tt>.
 *
 * The examples below show some common use cases.
 *
 * @section Examples
 *
 * Local alignment of two sequences using an @link Align @endlink object.
 *
 * @code{.cpp}
 * Dna5String seqH = "CGATT";
 * Dna5String seqV = "CGAAATT";
 *
 * Align<Dna5String> align;
 * resize(rows(align), 2);
 * assignSource(row(align, 0), seqH);
 * assignSource(row(align, 0), seqV);
 * Score<int, Simple> scoringScheme(2, -1, -2);
 *
 * int result = localAlignment(align, scoringScheme);
 * @endcode
 *
 * Local banded alignment of two sequences using two @link Gaps @endlink objects.
 *
 * @code{.cpp}
 * Dna5String seqH = "CGATT";
 * Gaps<Dna5String, ArrayGaps> gapsH(seqH);
 * DnaString seqV = "CGAAATT";
 * Gaps<Dna5String, AnchorGaps<> > gapsV(seqV);
 *
 * Score<int, Simple> scoringScheme(5, -3, -1, -5);
 *
 * int result = localAlignment(gapsH, gapsV, scoringScheme, -2, 2);
 * @endcode
 *
 * http://seqan.readthedocs.org/en/develop/Tutorial/PairwiseSequenceAlignment.html
 *
 * @section References
 *
 * <ul>
 *   <li>Smith TF, Waterman, MS: Identification of Common Molecular Subsequences. J Mol Biol 1981, 147(1):195-7.</li>
 * </ul>
 *
 * @see globalAlignment
 * @see LocalAlignmentEnumerator
 * @see PairwiseLocalAlignmentAlgorithms
 */

// ----------------------------------------------------------------------------
// Function localAlignment()                                  [unbanded, Align]
// ----------------------------------------------------------------------------

template <typename TSequence, typename TAlignSpec, typename TScoreValue, typename TScoreSpec, typename TTag>
TScoreValue localAlignment(Align<TSequence, TAlignSpec> & align,
                           Score<TScoreValue, TScoreSpec> const & scoringScheme,
                           TTag const & tag)
{
    SEQAN_ASSERT_EQ(length(rows(align)), 2u);
    typedef Align<TSequence, TAlignSpec> TAlign;
    typedef typename Size<TAlign>::Type TSize;
    typedef typename Position<TAlign>::Type TPosition;
    typedef TraceSegment_<TPosition, TSize> TTraceSegment;
    typedef AlignConfig2<DPLocal, DPBandConfig<BandOff>, FreeEndGaps_<> > TAlignConfig2;

    String<TTraceSegment> trace;
    DPScoutState_<Default> dpScoutState;
    TScoreValue res = _setUpAndRunAlignment(trace, dpScoutState, source(row(align, 0)), source(row(align, 1)),
                                            scoringScheme, TAlignConfig2(), tag);

    _adaptTraceSegmentsTo(row(align, 0), row(align, 1), trace);
    return res;
}

 template <typename TSequence, typename TAlignSpec,
           typename TScoreValue, typename TScoreSpec>
 TScoreValue localAlignment(Align<TSequence, TAlignSpec> & align,
                            Score<TScoreValue, TScoreSpec> const & scoringScheme)
 {
     SEQAN_ASSERT(length(rows(align)) == 2u);
     if (_usesAffineGaps(scoringScheme, source(row(align, 0)), source(row(align, 1))))
         return localAlignment(align, scoringScheme, AffineGaps());
     else
         return localAlignment(align, scoringScheme, LinearGaps());
 }

// ----------------------------------------------------------------------------
// Function localAlignment()                                   [unbanded, Gaps]
// ----------------------------------------------------------------------------

 template <typename TSequenceH, typename TGapsSpecH, typename TSequenceV, typename TGapsSpecV, typename TScoreValue,
           typename TScoreSpec, typename TTag>
 TScoreValue localAlignment(Gaps<TSequenceH, TGapsSpecH> & gapsH,
                            Gaps<TSequenceV, TGapsSpecV> & gapsV,
                            Score<TScoreValue, TScoreSpec> const & scoringScheme,
                            TTag const & tag)
 {
     typedef typename Size<TSequenceH>::Type TSize;
     typedef typename Position<TSequenceH>::Type TPosition;
     typedef TraceSegment_<TPosition, TSize> TTraceSegment;
     typedef AlignConfig2<DPLocal, DPBandConfig<BandOff>, FreeEndGaps_<> > TAlignConfig2;

     String<TTraceSegment> trace;
     DPScoutState_<Default> dpScoutState;
     TScoreValue res = _setUpAndRunAlignment(trace, dpScoutState, source(gapsH), source(gapsV), scoringScheme,
                                             TAlignConfig2(), tag);
     _adaptTraceSegmentsTo(gapsH, gapsV, trace);
     return res;
 }

 template <typename TSequenceH, typename TGapsSpecH,
          typename TSequenceV, typename TGapsSpecV,
          typename TScoreValue, typename TScoreSpec>
TScoreValue localAlignment(Gaps<TSequenceH, TGapsSpecH> & gapsH,
                           Gaps<TSequenceV, TGapsSpecV> & gapsV,
                           Score<TScoreValue, TScoreSpec> const & scoringScheme)
{
     if (_usesAffineGaps(scoringScheme, source(gapsH), source(gapsV)))
         return localAlignment(gapsH, gapsV, scoringScheme, AffineGaps());
     else
         return localAlignment(gapsH, gapsV, scoringScheme, LinearGaps());
}

// ----------------------------------------------------------------------------
// Function localAlignment()                     [unbanded, Graph<Alignment<>>]
// ----------------------------------------------------------------------------

// Full interface.

template <typename TStringSet, typename TCargo, typename TGraphSpec,
          typename TScoreValue, typename TScoreSpec, typename TTag>
TScoreValue localAlignment(Graph<Alignment<TStringSet, TCargo, TGraphSpec> > & alignmentGraph,
                           Score<TScoreValue, TScoreSpec> const & scoringScheme,
                           TTag const & tag)
{
    typedef Graph<Alignment<TStringSet, TCargo, TGraphSpec> > TGraph;
    typedef typename Size<TGraph>::Type TSize;
    typedef typename Position<TGraph>::Type TPosition;
    typedef TraceSegment_<TPosition, TSize> TTraceSegment;
    typedef AlignConfig2<DPLocal, DPBandConfig<BandOff>, FreeEndGaps_<> > TAlignConfig2;

    String<TTraceSegment> trace;
    DPScoutState_<Default> dpScoutState;
    TScoreValue res = _setUpAndRunAlignment(trace, dpScoutState, value(stringSet(alignmentGraph), 0),
                                            value(stringSet(alignmentGraph), 1), scoringScheme, TAlignConfig2(), tag);

    _adaptTraceSegmentsTo(alignmentGraph, positionToId(stringSet(alignmentGraph), 0),
                          positionToId(stringSet(alignmentGraph), 1), trace);
    return res;
}

template <typename TStringSet, typename TCargo, typename TGraphSpec,
          typename TScoreValue, typename TScoreSpec>
TScoreValue localAlignment(Graph<Alignment<TStringSet, TCargo, TGraphSpec> > & alignmentGraph,
                           Score<TScoreValue, TScoreSpec> const & scoringScheme)
{
    SEQAN_ASSERT(length(stringSet(alignmentGraph)) == 2u);

    if (_usesAffineGaps(scoringScheme, stringSet(alignmentGraph)[0], stringSet(alignmentGraph)[1]))
        return localAlignment(alignmentGraph, scoringScheme, AffineGaps());
    else
        return localAlignment(alignmentGraph, scoringScheme, LinearGaps());
}

// ----------------------------------------------------------------------------
// Function localAlignment()                    [unbanded, String<Fragment<> >]
// ----------------------------------------------------------------------------

// Full interface.

template <typename TSize, typename TFragmentSpec, typename TStringSpec,
          typename TSequence, typename TStringSetSpec,
          typename TScoreValue, typename TScoreSpec, typename TTag>
TScoreValue localAlignment(String<Fragment<TSize, TFragmentSpec>, TStringSpec> & fragmentString,
                           StringSet<TSequence, TStringSetSpec> const & strings,
                           Score<TScoreValue, TScoreSpec> const & scoringScheme,
                           TTag const & tag)
{
    typedef String<Fragment<TSize, TFragmentSpec>, TStringSpec> TFragments;
    typedef typename Position<TFragments>::Type TPosition;
    typedef TraceSegment_<TPosition, TSize> TTraceSegment;
    typedef AlignConfig2<DPLocal, DPBandConfig<BandOff>, FreeEndGaps_<> > TAlignConfig2;

    String<TTraceSegment> trace;
    DPScoutState_<Default> dpScoutState;
    TScoreValue res = _setUpAndRunAlignment(trace, dpScoutState, value(strings, 0), value(strings, 1), scoringScheme,
                                            TAlignConfig2(), tag);

    _adaptTraceSegmentsTo(fragmentString, positionToId(strings, 0), positionToId(strings, 1), trace);
    return res;
}

template <typename TSize, typename TFragmentSpec, typename TStringSpec,
          typename TSequence, typename TStringSetSpec,
          typename TScoreValue, typename TScoreSpec>
TScoreValue localAlignment(String<Fragment<TSize, TFragmentSpec>, TStringSpec> & fragmentString,
                           StringSet<TSequence, TStringSetSpec> const & strings,
                           Score<TScoreValue, TScoreSpec> const & scoringScheme)
{
    SEQAN_ASSERT(length(strings) == 2u);

    if (_usesAffineGaps(scoringScheme, strings[0], strings[1]))
        return localAlignment(fragmentString, strings, scoringScheme, AffineGaps());
    else
        return localAlignment(fragmentString, strings, scoringScheme, LinearGaps());
}

}  // namespace seqan

#endif  // #ifndef SEQAN_INCLUDE_SEQAN_ALIGN_LOCAL_ALIGNMENT_UNBANDED_H_
