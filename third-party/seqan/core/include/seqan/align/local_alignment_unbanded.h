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
// Author: Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>
// ==========================================================================
// Interface functions for unbanded local alignment.
// ==========================================================================

#ifndef SEQAN_CORE_INCLUDE_SEQAN_ALIGN_LOCAL_ALIGNMENT_UNBANDED_H_
#define SEQAN_CORE_INCLUDE_SEQAN_ALIGN_LOCAL_ALIGNMENT_UNBANDED_H_

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
 * @param lowerDiag Optional lower diagonal (<tt>int</tt>).
 * @param lowerDiag Optional upper diagonal (<tt>int</tt>).
 *
 * @param gapsH Horizontal gapped sequence in alignment matrix. Types: Gaps
 * @param align An @link Align @endlink object that stores the alignment. The
 *              number of rows must be 2 and the sequences must have already
 *              been set. <tt>align[0]</tt> is the horizontal one in the
 *              alignment matrix alignment, <tt>align[1]</tt> is the vertical
 *              one. Types: Align
 * @param fragmentString String of @link Fragment @endlink objects. The sequence
 *                       with id <tt>0</tt> is the horizontal one, the sequence
 *                       with id <tt>1</tt> is the vertical one.
 * @param gapsV Vertical gapped sequence in alignment matrix. Types: Gaps
 * @param scoringScheme The scoring scheme to use for the alignment. Note that
 *                      the user is responsible for ensuring that the scoring
 *                      scheme is compatible with <tt>algorithmTag</tt>. Types:
 *                      Score
 * 
 * @return TScoreVal The score value of the alignmetn.
 * 
 * @section Remarks
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
 * alignments, for example in the construction of @link Alignment Graph Alignment Graphs @endlink for multiple- sequence
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
 * http://trac.seqan.de/wiki/Tutorial/PairwiseSequenceAlignment
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

/**
.Function.localAlignment
..summary:Computes the best pairwise local alignment using the Smith-Waterman algorithm.
..cat:Alignments
..signature:localAlignment(align,          scoringScheme, [lowerDiag, upperDiag])
..signature:localAlignment(gapsH, gapsV,   scoringScheme, [lowerDiag, upperDiag])
..signature:localAlignment(fragmentString, scoringScheme, [lowerDiag, upperDiag])
..param.align:
An @Class.Align@ object that stores the alignment.
The number of rows must be 2 and the sequences must have already been set.
$align[0]$ is the horizontal one in the alignment matrix alignment, $align[1]$ is the vertical one.
...type:Class.Align
..param.gapsH:Horizontal gapped sequence in alignment matrix.
...type:Class.Gaps
..param.gapsV:Vertical gapped sequence in alignment matrix.
...type:Class.Gaps
..param.fragmentString:
String of @Class.Fragment@ objects.
The sequence with id $0$ is the horizontal one, the sequence with id $1$ is the vertical one.
..param.scoringScheme:
The scoring scheme to use for the alignment.
Note that the user is responsible for ensuring that the scoring scheme is compatible with $algorithmTag$.
...type:Class.Score
..param.lowerDiag:Optional lower diagonal.
...type:nolink:$int$
..param.upperDiag:Optional upper diagonal.
...type:nolink:$int$
..returns:An integer with the alignment score, as given by the @Metafunction.Value@ metafunction of the @Class.Score@ type.
..remarks:The Waterman-Eggert algorithm (local alignment with declumping) is available through the @Class.LocalAlignmentEnumerator@ class.
..remarks:
When using @Class.Gaps@ and @Class.Align@ objects, only parts (i.e. one infix) of each sequence will be aligned.
This will be presented to the user by setting the clipping begin and end position of the gaps (the rows in the case of @Class.Align@ objects).
When using @Class.Fragment@ strings, these parts of the sequences will not appear in any fragment.
..remarks:
There exist multiple overloads for this function with two configuration dimensions.
..remarks:
First, you can select the type of the target storing the alignment.
This can be either an @Class.Align@ object, two @Class.Gaps@ objects, or a string of @Class.Fragment@ objects.
@Class.Align@ objects provide an interface to tabular alignments with the restriction of all rows having the same type.
Using two @Class.Gaps@ objects has the advantage that you an align sequences with different types, for example @Shortcut.DnaString@ and @Shortcut.Dna5String@.
Using @Class.Fragment@ strings is useful for collecting many pairwise alignments, for example in the construction of @Spec.Alignment Graph|Alignment Graphs@ for multiple-sequence alignments (MSA).
..remarks:
Second, you can optionally give a band for the alignment using $lowerDiag$ and $upperDiag$.
The center diagonal has index $0$, the $i$th diagonal below has index $-i$, the $i$th above has index $i$.
..remarks:
The examples below show some common use cases.
..example.text:Local alignment of two sequences using an @Class.Align@ object.
..example.code:
Dna5String seqH = "CGATT";
Dna5String seqV = "CGAAATT";

Align<Dna5String> align;
resize(rows(align), 2);
assignSource(row(align, 0), seqH);
assignSource(row(align, 0), seqV);
Score<int, Simple> scoringScheme(2, -1, -2);

int result = localAlignment(align, scoringScheme);
..example.text:Local banded alignment of two sequences using two @Class.Gaps@ objects.
..example.code:
Dna5String seqH = "CGATT";
Gaps<Dna5String, ArrayGaps> gapsH(seqH);
DnaString seqV = "CGAAATT";
Gaps<Dna5String, AnchorGaps<> > gapsV(seqV);

Score<int, Simple> scoringScheme(5, -3, -1, -5);

int result = localAlignment(gapsH, gapsV, scoringScheme, -2, 2);
..see:Function.globalAlignment
..see:Class.LocalAlignmentEnumerator
..include:seqan/align.h
..wiki:Tutorial/PairwiseSequenceAlignment
..cite:Smith TF, Waterman, MS: Identification of Common Molecular Subsequences. J Mol Biol 1981, 147(1):195-7.
.
*/

// ----------------------------------------------------------------------------
// Function localAlignment()                                  [unbanded, Align]
// ----------------------------------------------------------------------------

template <typename TSequence, typename TAlignSpec,
          typename TScoreValue, typename TScoreSpec>
TScoreValue localAlignment(Align<TSequence, TAlignSpec> & align,
                           Score<TScoreValue, TScoreSpec> const & scoringScheme)
{
    SEQAN_ASSERT_EQ(length(rows(align)), 2u);
    typedef Align<TSequence, TAlignSpec> TAlign;
    typedef typename Size<TAlign>::Type TSize;
    typedef typename Position<TAlign>::Type TPosition;
    typedef TraceSegment_<TPosition, TSize> TTraceSegment;

    String<TTraceSegment> traceSegments;
    TScoreValue score = _setUpAndRunAlignment(traceSegments, source(row(align, 0)), source(row(align, 1)),
                                              scoringScheme, SmithWaterman());
    _adaptTraceSegmentsTo(row(align, 0), row(align, 1), traceSegments);
    return score;
}

// ----------------------------------------------------------------------------
// Function localAlignment()                                   [unbanded, Gaps]
// ----------------------------------------------------------------------------

template <typename TSequenceH, typename TGapsSpecH,
          typename TSequenceV, typename TGapsSpecV,
          typename TScoreValue, typename TScoreSpec>
TScoreValue localAlignment(Gaps<TSequenceH, TGapsSpecH> & gapsH,
                           Gaps<TSequenceV, TGapsSpecV> & gapsV,
                           Score<TScoreValue, TScoreSpec> const & scoringScheme)
{
    typedef typename Size<TSequenceH>::Type TSize;
    typedef typename Position<TSequenceH>::Type TPosition;
    typedef TraceSegment_<TPosition, TSize> TTraceSegment;

    String<TTraceSegment> traceSegments;
    TScoreValue score = _setUpAndRunAlignment(traceSegments, source(gapsH), source(gapsV), scoringScheme,
                                              SmithWaterman());
    _adaptTraceSegmentsTo(gapsH, gapsV, traceSegments);
    return score;
}

// ----------------------------------------------------------------------------
// Function localAlignment()                     [unbanded, Graph<Alignment<>>]
// ----------------------------------------------------------------------------

// Full interface.

template <typename TStringSet, typename TCargo, typename TGraphSpec,
          typename TScoreValue, typename TScoreSpec>
TScoreValue localAlignment(Graph<Alignment<TStringSet, TCargo, TGraphSpec> > & alignmentGraph,
                           Score<TScoreValue, TScoreSpec> const & scoringScheme)
{
    typedef Graph<Alignment<TStringSet, TCargo, TGraphSpec> > TGraph;
    typedef typename Size<TGraph>::Type TSize;
    typedef typename Position<TGraph>::Type TPosition;
    typedef TraceSegment_<TPosition, TSize> TTraceSegment;

    String<TTraceSegment> traceSegments;
    TScoreValue score = _setUpAndRunAlignment(traceSegments, value(stringSet(alignmentGraph), 0),
                                              value(stringSet(alignmentGraph), 1), scoringScheme, SmithWaterman());
    _adaptTraceSegmentsTo(alignmentGraph, positionToId(stringSet(alignmentGraph), 0),
                          positionToId(stringSet(alignmentGraph), 1), traceSegments);
    return score;
}

// ----------------------------------------------------------------------------
// Function localAlignment()                    [unbanded, String<Fragment<> >]
// ----------------------------------------------------------------------------

// Full interface.

template <typename TSize, typename TFragmentSpec, typename TStringSpec,
          typename TSequence, typename TStringSetSpec,
          typename TScoreValue, typename TScoreSpec>
TScoreValue localAlignment(String<Fragment<TSize, TFragmentSpec>, TStringSpec> & fragmentString,
                           StringSet<TSequence, TStringSetSpec> const & strings,
                           Score<TScoreValue, TScoreSpec> const & scoringScheme)
{
    typedef String<Fragment<TSize, TFragmentSpec>, TStringSpec> TFragments;
    typedef typename Position<TFragments>::Type TPosition;
    typedef TraceSegment_<TPosition, TSize> TTraceSegment;

    String<TTraceSegment> traceSegments;
    TScoreValue score = _setUpAndRunAlignment(traceSegments, value(strings, 0), value(strings, 1), scoringScheme,
                                              SmithWaterman());
    _adaptTraceSegmentsTo(fragmentString, positionToId(strings, 0), positionToId(strings, 1), traceSegments);
    return score;
}

}  // namespace seqan

#endif  // #ifndef SEQAN_CORE_INCLUDE_SEQAN_ALIGN_LOCAL_ALIGNMENT_UNBANDED_H_
