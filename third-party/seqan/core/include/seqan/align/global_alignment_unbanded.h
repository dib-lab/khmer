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
// Global alignment interface for the unbanded Needleman-Wunsch and Gotoh
// algorithms.
//
// We define the interface functions pretty explicitely (versus just TAlign,
// TFragments etc.) so the candidates the compiler gives when resolution to
// the globalFunction() fails is actually meaningful.
// ==========================================================================

#ifndef SEQAN_CORE_INCLUDE_SEQAN_ALIGN_GLOBAL_ALIGNMENT_UNBANDED_H_
#define SEQAN_CORE_INCLUDE_SEQAN_ALIGN_GLOBAL_ALIGNMENT_UNBANDED_H_

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

template <typename TScoreValue, typename TSpec>
class Score;
template <typename TSpec>
class Graph;
template <typename TStringSet, typename TCargo, typename TGraphSpec>
struct Alignment;
template <typename TSize, typename TFragmentSpec>
class Fragment;

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
// Function globalAlignment()
// ----------------------------------------------------------------------------

/*!
 * @fn globalAlignment
 * 
 * @headerfile seqan/align.h
 * 
 * @brief Computes the best global pairwise alignment.
 * 
 * @signature TScoreVal globalAlignment(align,          scoringScheme, [alignConfig,] [lowerDiag, upperDiag,] [algorithmTag]);
 * @signature TScoreVal globalAlignment(gapsH, gapsV,   scoringScheme, [alignConfig,] [lowerDiag, upperDiag,] [algorithmTag]);
 * @signature TScoreVal globalAlignment(frags, strings, scoringScheme, [alignConfig,] [lowerDiag, upperDiag,] [algorithmTag]);
 * @signature TScoreVal globalAlignment(alignGraph,     scoringScheme, [alignConfig,] [lowerDiag, upperDiag,] [algorithmTag]);
 * 
 * @param align        The @link Align @endlink object to use for storing the pairwise alignment.
 * @param gapsH        The @link Gaps @endlink object for the first row (horizontal in the DP matrix).
 * @param gapsV        The @link Gaps @endlink object for the second row (vertical in the DP matrix).
 * @param frags        String of @link Fragment @endlink objects to store alignment in.
 * @param strings      StringSet of length two with the strings to align.
 * @param alignGraph   Alignment Graph for the resulting alignment.  Must be initialized with two strings.
 * @param scoringScheme The @link Score scoring scheme @endlink to use for the alignment.  Note that
 *                      the user is responsible for ensuring that the scoring scheme is compatible with <tt>algorithmTag</tt>.
 * @param alignConfig  @link AlignConfig @endlink instance to use for the alignment configuration.
 * @param lowerDiag    Optional lower diagonal (<tt>int</tt>).
 * @param upperDiag    Optional upper diagonal (<tt>int</tt>).
 * @param algorithmTag Tag to select the alignment algorithm (see @link AlignmentAlgorithmTags @endlink).
 *
 * @return TScoreVal Score value of the resulting alignment.  Of type <tt>Value<TScore>::Type</tt> where
 *                   <tt>TScore</tt> is the type of <tt>scoringScheme</tt>.
 * 
 * There exist multiple overloads for this function with four configuration dimensions.
 * 
 * First, you can select whether begin and end gaps are free in either sequence using <tt>alignConfig</tt>.
 * 
 * Second, you can select the type of the target storing the alignment. This can be either an @link Align @endlink
 * object, two @link Gaps @endlink objects, a @link Alignment Graph @endlink, or a string of @link Fragment @endlink
 * objects. @link Align @endlink objects provide an interface to tabular alignments with the restriction of all rows
 * having the same type. Using two @link Gaps @endlink objects has the advantage that you an align sequences with
 * different types, for example @link DnaString @endlink and @link Dna5String @endlink. @link Alignment Graph Alignment
 * Graphs @endlink provide a graph-based representation of segment-based colinear alignments. Using @link Fragment
 * @endlink strings is useful for collecting many pairwise alignments, for example in the construction of @link
 * Alignment Graph Alignment Graphs @endlink for multiple-sequence alignments (MSA).
 * 
 * Third, you can optionally give a band for the alignment using <tt>lowerDiag</tt> and <tt>upperDiag</tt>. The center
 * diagonal has index <tt>0</tt>, the <tt>i</tt>th diagonal below has index <tt>-i</tt>, the <tt>i</tt>th above has
 * index <tt>i</tt>.
 * 
 * Fourth, you can select the algorithm to use with <tt>algorithmTag</tt>. This can be one of @link
 * AlignmentAlgorithmTags @endlink and @link Pairwise Global Alignment Algorithms.value.Gotoh @endlink.  The
 * Needleman-Wunsch algorithm supports scoring schemes with linear gap costs only while Gotoh's algorithm also allows
 * affine gap costs.
 * 
 * The available alignment algorithms all have some restrictions.  Gotoh's algorithm can handle arbitrary substitution
 * and affine gap scores.  Needleman-Wunsch is limited to linear gap scores.  The implementation of Hirschberg's
 * algorithm is further limited that it does not support <tt>alignConfig</tt> objects or banding.  The implementation of
 * the Myers-Hirschberg algorithm further limits this to only support edit distance (as scores, matches are scored with
 * 0, mismatches are scored with -1).
 * 
 * The examples below show some common use cases.
 * 
 * @section Examples
 * 
 * Global alignment of two sequences using an @link Align @endlink object and
 * the Needleman-Wunsch algorithm.
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
 * AlignConfig<> alignConfig;
 *  
 * int result = globalAlignment(align, scoringScheme, alignConfig,
 *                              NeedlemanWunsch());
 * @endcode
 *
 * Global banded alignment of two sequences using two @link Gaps @endlink objects and the Gotoh algorithm.
 * 
 * @code{.cpp}
 * Dna5String seqH = "CGATT";
 * Gaps<Dna5String, ArrayGaps> gapsH(seqH);
 * DnaString seqV = "CGAAATT";
 * Gaps<Dna5String, AnchorGaps<> > gapsV(seqV);
 *  
 * Score<int, Simple> scoringScheme(5, -3, -1, -5);
 * AlignConfig<> alignConfig;
 *  
 * int result = globalAlignment(gapsH, gapsV, scoringScheme, alignConfig, -2, 2);
 * @endcode
 * 
 * http://trac.seqan.de/wiki/Tutorial/PairwiseSequenceAlignment
 * 
 * @section References
 *
 * <ul>
 *   <li>Needleman SB, Wunsch CD: A general method applicable to the search for similarities in the amino acid sequence
 *       of two proteins. J Mol Biol 1970, 48(3): 443-53.</li>
 *   <li>Gotoh O: An improved algorithm for matching biological sequences. J Mol Biol 1982, 162(3):705-8</li>
 * </ul>
 * 
 * @see localAlignment
 * @see globalAlignmentScore
 * @see AlignmentAlgorithmTags
 */

/**
.Function.globalAlignment
..summary:Computes the best global pairwise alignment.
..cat:Alignments
..signature:globalAlignment(align,          scoringScheme, [alignConfig,] [lowerDiag, upperDiag,] [algorithmTag])
..signature:globalAlignment(gapsH, gapsV,   scoringScheme, [alignConfig,] [lowerDiag, upperDiag,] [algorithmTag])
..signature:globalAlignment(frags, strings, scoringScheme, [alignConfig,] [lowerDiag, upperDiag,] [algorithmTag])
..signature:globalAlignment(alignmentGraph, scoringScheme, [alignConfig,] [lowerDiag, upperDiag,] [algorithmTag])
..param.align:
An @Class.Align@ object that stores the alignment.
The number of rows must be 2 and the sequences must have already been set.
$row(align, 0)$ is the horizontal one in the alignment matrix alignment, $row(align, 1)$ is the vertical one.
...type:Class.Align
..param.gapsH:Horizontal gapped sequence in alignment matrix.
...type:Class.Gaps
..param.gapsV:Vertical gapped sequence in alignment matrix.
...type:Class.Gaps
..param.frags:
String of @Class.Fragment@ objects.
The sequence with id $0$ is the horizontal one, the sequence with id $1$ is the vertical one.
..param.alignmentGraph:
@Spec.Alignment Graph@ object to store the alignment in.
...type:Spec.Alignment Graph
...remarks:The underlying @Class.StringSet@ must be an @Spec.Owner|Owner StringSet@.
..param.strings:A @Class.StringSet@ containing two sequences.
...type:Class.StringSet
..param.scoringScheme:
The scoring scheme to use for the alignment.
Note that the user is responsible for ensuring that the scoring scheme is compatible with $algorithmTag$.
...type:Class.Score
..param.alignConfig:The @Class.AlignConfig@ to use for the alignment.
...type:Class.AlignConfig
..param.lowerDiag:Optional lower diagonal.
...type:nolink:$int$
..param.upperDiag:Optional upper diagonal.
...type:nolink:$int$
..param.algorithmTag:The Tag for picking the alignment algorithm.
...type:Tag.Pairwise Global Alignment Algorithms.tag.Gotoh
...type:Tag.Pairwise Global Alignment Algorithms.tag.NeedlemanWunsch
...type:Tag.Pairwise Global Alignment Algorithms.tag.Hirschberg
...type:Tag.Pairwise Global Alignment Algorithms.tag.MyersHirschberg
..returns:An integer with the alignment score, as given by the @Metafunction.Value@ metafunction of the @Class.Score@ type.
..remarks:
There exist multiple overloads for this function with four configuration dimensions.
..remarks:
First, you can select whether begin and end gaps are free in either sequence using $alignConfig$.
..remarks:
Second, you can select the type of the target storing the alignment.
This can be either an @Class.Align@ object, two @Class.Gaps@ objects, a @Spec.Alignment Graph@, or a string of @Class.Fragment@ objects.
@Class.Align@ objects provide an interface to tabular alignments with the restriction of all rows having the same type.
Using two @Class.Gaps@ objects has the advantage that you an align sequences with different types, for example @Shortcut.DnaString@ and @Shortcut.Dna5String@.
@Spec.Alignment Graph|Alignment Graphs@ provide a graph-based representation of segment-based colinear alignments.
Using @Class.Fragment@ strings is useful for collecting many pairwise alignments, for example in the construction of @Spec.Alignment Graph|Alignment Graphs@ for multiple-sequence alignments (MSA).
..remarks:
Third, you can optionally give a band for the alignment using $lowerDiag$ and $upperDiag$.
The center diagonal has index $0$, the $i$th diagonal below has index $-i$, the $i$th above has index $i$.
..remarks:
Fourth, you can select the algorithm to use with $algorithmTag$.
This can be one of @Tag.Pairwise Global Alignment Algorithms.value.NeedlemanWunsch@ and @Tag.Pairwise Global Alignment Algorithms.value.Gotoh@.
The Needleman-Wunsch algorithm supports scoring schemes with linear gap costs only while Gotoh's algorithm also allows affine gap costs.
..remarks:
The available alignment algorithms all have some restrictions.
Gotoh's algorithm can handle arbitrary substitution and affine gap scores.
Needleman-Wunsch is limited to linear gap scores.
The implementation of Hirschberg's algorithm is further limited that it does not support $alignConfig$ objects or banding.
The implementation of the Myers-Hirschberg algorithm further limits this to only support edit distance (as scores, matches are scored with 0, mismatches are scored with -1).
..remarks:
The examples below show some common use cases.
..example.text:Global alignment of two sequences using an @Class.Align@ object and the Needleman-Wunsch algorithm. The Needleman-Wunsch algorithm is automatically selected since the scoring scheme uses linear gap costs.
..example.file:demos/align/global_alignment_unbanded.cpp
..example.text:Global banded alignment of two sequences using two @Class.Gaps@ objects and the Gotoh algorithm. The Gotoh algorithm is automatically selected since the scoring scheme uses affine gap costs.
..example.file:demos/align/global_alignment_banded.cpp
..see:Function.localAlignment
..see:Function.globalAlignmentScore
..include:seqan/align.h
..wiki:Tutorial/PairwiseSequenceAlignment
..cite:Needleman SB, Wunsch CD: A general method applicable to the search for similarities in the amino acid sequence of two proteins. J Mol Biol 1970, 48(3): 443-53.
..cite:Gotoh O: An improved algorithm for matching biological sequences. J Mol Biol 1982, 162(3):705-8
.
*/

// ----------------------------------------------------------------------------
// Function globalAlignment()                                 [unbanded, Align]
// ----------------------------------------------------------------------------

template <typename TSequence, typename TAlignSpec,
          typename TScoreValue, typename TScoreSpec,
          bool TOP, bool LEFT, bool RIGHT, bool BOTTOM, typename TACSpec,
          typename TAlgoTag>
TScoreValue globalAlignment(Align<TSequence, TAlignSpec> & align,
                            Score<TScoreValue, TScoreSpec> const & scoringScheme,
                            AlignConfig<TOP, LEFT, RIGHT, BOTTOM, TACSpec> const & alignConfig,
                            TAlgoTag const & algoTag)
{
    typedef Align<TSequence, TAlignSpec> TAlign;
    typedef typename Size<TAlign>::Type TSize;
    typedef typename Position<TAlign>::Type TPosition;
    typedef TraceSegment_<TPosition, TSize> TTraceSegment;

    String<TTraceSegment> trace;

    // We do not need string ids for this variant and set them to 0u.  They are
    // only required for the Fragment String and the Alignment Graph variant.
    TScoreValue res = _setUpAndRunAlignment(trace, source(row(align, 0)), source(row(align, 1)), scoringScheme,
                                            alignConfig, algoTag);
    _adaptTraceSegmentsTo(row(align, 0), row(align, 1), trace);
    return res;
}

// Interface without AlignConfig<>.

template <typename TSequence, typename TAlignSpec,
          typename TScoreValue, typename TScoreSpec,
          typename TAlgoTag>
TScoreValue globalAlignment(Align<TSequence, TAlignSpec> & align,
                            Score<TScoreValue, TScoreSpec> const & scoringScheme,
                            TAlgoTag const & algoTag)
{
    AlignConfig<> alignConfig;
    return globalAlignment(align, scoringScheme, alignConfig, algoTag);
}

// Interface without algorithm tag.

template <typename TSequence, typename TAlignSpec,
          typename TScoreValue, typename TScoreSpec,
          bool TOP, bool LEFT, bool RIGHT, bool BOTTOM, typename TACSpec>
TScoreValue globalAlignment(Align<TSequence, TAlignSpec> & align,
                            Score<TScoreValue, TScoreSpec> const & scoringScheme,
                            AlignConfig<TOP, LEFT, RIGHT, BOTTOM, TACSpec> const & alignConfig)
{
    if (scoreGapOpen(scoringScheme) == scoreGapExtend(scoringScheme))
        return globalAlignment(align, scoringScheme, alignConfig, NeedlemanWunsch());
    else
        return globalAlignment(align, scoringScheme, alignConfig, Gotoh());
}

// Interface without AlignConfig<> and algorithm tag.

template <typename TSequence, typename TAlignSpec,
          typename TScoreValue, typename TScoreSpec>
TScoreValue globalAlignment(Align<TSequence, TAlignSpec> & align,
                            Score<TScoreValue, TScoreSpec> const & scoringScheme)
{
    AlignConfig<> alignConfig;
    return globalAlignment(align, scoringScheme, alignConfig);
}

// ----------------------------------------------------------------------------
// Function globalAlignment()                                  [unbanded, Gaps]
// ----------------------------------------------------------------------------

template <typename TSequenceH, typename TGapsSpecH,
          typename TSequenceV, typename TGapsSpecV,
          typename TScoreValue, typename TScoreSpec,
          bool TOP, bool LEFT, bool RIGHT, bool BOTTOM, typename TACSpec,
          typename TAlgoTag>
TScoreValue globalAlignment(Gaps<TSequenceH, TGapsSpecH> & gapsH,
                            Gaps<TSequenceV, TGapsSpecV> & gapsV,
                            Score<TScoreValue, TScoreSpec> const & scoringScheme,
                            AlignConfig<TOP, LEFT, RIGHT, BOTTOM, TACSpec> const & alignConfig,
                            TAlgoTag const & algoTag)
{
    typedef typename Size<TSequenceH>::Type TSize;
    typedef typename Position<TSequenceH>::Type TPosition;
    typedef TraceSegment_<TPosition, TSize> TTraceSegment;

    String<TTraceSegment> traceSegments;

    // We do not need string ids for this variant and set them to 0u.  They are
    // only required for the Fragment String and the Alignment Graph variant.
    TScoreValue res = _setUpAndRunAlignment(traceSegments, source(gapsH), source(gapsV), scoringScheme, alignConfig,
                                            algoTag);
    _adaptTraceSegmentsTo(gapsH, gapsV, traceSegments);
    return res;
}

// Interface without AlignConfig<>.

template <typename TSequenceH, typename TGapsSpecH,
          typename TSequenceV, typename TGapsSpecV,
          typename TScoreValue, typename TScoreSpec,
          typename TAlgoTag>
TScoreValue globalAlignment(Gaps<TSequenceH, TGapsSpecH> & gapsH,
                            Gaps<TSequenceV, TGapsSpecV> & gapsV,
                            Score<TScoreValue, TScoreSpec> const & scoringScheme,
                            TAlgoTag const & algoTag)
{
    AlignConfig<> alignConfig;
    return globalAlignment(gapsH, gapsV, scoringScheme, alignConfig, algoTag);
}

// Interface without algorithm tag.

template <typename TSequenceH, typename TGapsSpecH,
          typename TSequenceV, typename TGapsSpecV,
          typename TScoreValue, typename TScoreSpec,
          bool TOP, bool LEFT, bool RIGHT, bool BOTTOM, typename TACSpec>
TScoreValue globalAlignment(Gaps<TSequenceH, TGapsSpecH> & gapsH,
                            Gaps<TSequenceV, TGapsSpecV> & gapsV,
                            Score<TScoreValue, TScoreSpec> const & scoringScheme,
                            AlignConfig<TOP, LEFT, RIGHT, BOTTOM, TACSpec> const & alignConfig)
{
    if (scoreGapOpen(scoringScheme) == scoreGapExtend(scoringScheme))
        return globalAlignment(gapsH, gapsV, scoringScheme, alignConfig, NeedlemanWunsch());
    else
        return globalAlignment(gapsH, gapsV, scoringScheme, alignConfig, Gotoh());
}

// Interface without AlignConfig<> and algorithm tag.

template <typename TSequenceH, typename TGapsSpecH,
          typename TSequenceV, typename TGapsSpecV,
          typename TScoreValue, typename TScoreSpec>
TScoreValue globalAlignment(Gaps<TSequenceH, TGapsSpecH> & gapsH,
                            Gaps<TSequenceV, TGapsSpecV> & gapsV,
                            Score<TScoreValue, TScoreSpec> const & scoringScheme)
{
    AlignConfig<> alignConfig;
    return globalAlignment(gapsH, gapsV, scoringScheme, alignConfig);
}

// ----------------------------------------------------------------------------
// Function globalAlignment()                   [unbanded, Graph<Alignment<> >]
// ----------------------------------------------------------------------------

// Full interface.

template <typename TStringSet, typename TCargo, typename TGraphSpec,
          typename TScoreValue, typename TScoreSpec,
          bool TOP, bool LEFT, bool RIGHT, bool BOTTOM, typename TACSpec,
          typename TAlgoTag>
TScoreValue globalAlignment(Graph<Alignment<TStringSet, TCargo, TGraphSpec> > & alignmentGraph,
                            Score<TScoreValue, TScoreSpec> const & scoringScheme,
                            AlignConfig<TOP, LEFT, RIGHT, BOTTOM, TACSpec> const & alignConfig,
                            TAlgoTag const & algoTag)
{
    typedef Graph<Alignment<TStringSet, TCargo, TGraphSpec> > TGraph;
    typedef typename Position<TGraph>::Type TPosition;
    typedef typename Size<TGraph>::Type TSize;
    typedef TraceSegment_<TPosition, TSize> TTraceSegment;

    String<TTraceSegment> traceSegments;

    TScoreValue res = _setUpAndRunAlignment(traceSegments, value(stringSet(alignmentGraph), 0),
                                            value(stringSet(alignmentGraph), 1), scoringScheme, alignConfig, algoTag);
    _adaptTraceSegmentsTo(alignmentGraph, positionToId(stringSet(alignmentGraph), 0),
                          positionToId(stringSet(alignmentGraph), 1), traceSegments);
    return res;
}

// Interface without AlignConfig<>.

template <typename TStringSet, typename TCargo, typename TGraphSpec,
          typename TScoreValue, typename TScoreSpec,
          typename TAlgoTag>
TScoreValue globalAlignment(Graph<Alignment<TStringSet, TCargo, TGraphSpec> > & alignmentGraph,
                            Score<TScoreValue, TScoreSpec> const & scoringScheme,
                            TAlgoTag const & algoTag)
{
    AlignConfig<> alignConfig;
    return globalAlignment(alignmentGraph, scoringScheme, alignConfig, algoTag);
}

// Interface without algorithm tag.

template <typename TStringSet, typename TCargo, typename TGraphSpec,
          typename TScoreValue, typename TScoreSpec,
          bool TOP, bool LEFT, bool RIGHT, bool BOTTOM, typename TACSpec>
TScoreValue globalAlignment(Graph<Alignment<TStringSet, TCargo, TGraphSpec> > & alignmentGraph,
                            Score<TScoreValue, TScoreSpec> const & scoringScheme,
                            AlignConfig<TOP, LEFT, RIGHT, BOTTOM, TACSpec> const & alignConfig)
{
    if (scoreGapOpen(scoringScheme) == scoreGapExtend(scoringScheme))
        return globalAlignment(alignmentGraph, scoringScheme, alignConfig, NeedlemanWunsch());
    else
        return globalAlignment(alignmentGraph, scoringScheme, alignConfig, Gotoh());
}

// Interface without AlignConfig<> and algorithm tag.

template <typename TStringSet, typename TCargo, typename TGraphSpec,
          typename TScoreValue, typename TScoreSpec>
TScoreValue globalAlignment(Graph<Alignment<TStringSet, TCargo, TGraphSpec> > & alignmentGraph,
                            Score<TScoreValue, TScoreSpec> const & scoringScheme)
{
    AlignConfig<> alignConfig;
    return globalAlignment(alignmentGraph, scoringScheme, alignConfig);
}

// ----------------------------------------------------------------------------
// Function globalAlignment()                   [unbanded, String<Fragment<> >]
// ----------------------------------------------------------------------------

// Full interface.

template <typename TSize, typename TFragmentSpec, typename TStringSpec,
          typename TSequence, typename TStringSetSpec,
          typename TScoreValue, typename TScoreSpec,
          bool TOP, bool LEFT, bool RIGHT, bool BOTTOM, typename TACSpec,
          typename TAlgoTag>
TScoreValue globalAlignment(String<Fragment<TSize, TFragmentSpec>, TStringSpec> & fragmentString,
                            StringSet<TSequence, TStringSetSpec> const & strings,
                            Score<TScoreValue, TScoreSpec> const & scoringScheme,
                            AlignConfig<TOP, LEFT, RIGHT, BOTTOM, TACSpec> const & alignConfig,
                            TAlgoTag const & algoTag)
{
    typedef String<Fragment<TSize, TFragmentSpec>, TStringSpec> TFragments;
    typedef typename Position<TFragments>::Type TPosition;
    typedef TraceSegment_<TPosition, TSize> TTraceSegment;

    String<TTraceSegment> traceSegments;

    TScoreValue res = _setUpAndRunAlignment(traceSegments, value(strings, 0), value(strings, 1), scoringScheme,
                                            alignConfig, algoTag);
    _adaptTraceSegmentsTo(fragmentString, positionToId(strings, 0), positionToId(strings, 1), traceSegments);
    return res;
}

// Interface without AlignConfig<>.

template <typename TSize, typename TFragmentSpec, typename TStringSpec,
          typename TSequence, typename TStringSetSpec,
          typename TScoreValue, typename TScoreSpec,
          typename TAlgoTag>
TScoreValue globalAlignment(String<Fragment<TSize, TFragmentSpec>, TStringSpec> & fragmentString,
                            StringSet<TSequence, TStringSetSpec> const & strings,
                            Score<TScoreValue, TScoreSpec> const & scoringScheme,
                            TAlgoTag const & algoTag)
{
    AlignConfig<> alignConfig;
    return globalAlignment(fragmentString, strings, scoringScheme, alignConfig, algoTag);
}

// Interface without algorithm tag.

template <typename TSize, typename TFragmentSpec, typename TStringSpec,
          typename TSequence, typename TStringSetSpec,
          typename TScoreValue, typename TScoreSpec,
          bool TOP, bool LEFT, bool RIGHT, bool BOTTOM, typename TACSpec>
TScoreValue globalAlignment(String<Fragment<TSize, TFragmentSpec>, TStringSpec> & fragmentString,
                            StringSet<TSequence, TStringSetSpec> const & strings,
                            Score<TScoreValue, TScoreSpec> const & scoringScheme,
                            AlignConfig<TOP, LEFT, RIGHT, BOTTOM, TACSpec> const & alignConfig)
{
    if (scoreGapOpen(scoringScheme) == scoreGapExtend(scoringScheme))
        return globalAlignment(fragmentString, strings, scoringScheme, alignConfig, NeedlemanWunsch());
    else
        return globalAlignment(fragmentString, strings, scoringScheme, alignConfig, Gotoh());
}

// Interface without AlignConfig<> and algorithm tag.

template <typename TSize, typename TFragmentSpec, typename TStringSpec,
          typename TSequence, typename TStringSetSpec,
          typename TScoreValue, typename TScoreSpec>
TScoreValue globalAlignment(String<Fragment<TSize, TFragmentSpec>, TStringSpec> & fragmentString,
                            StringSet<TSequence, TStringSetSpec> const & strings,
                            Score<TScoreValue, TScoreSpec> const & scoringScheme)
{
    AlignConfig<> alignConfig;
    return globalAlignment(fragmentString, strings, scoringScheme, alignConfig);
}

// ----------------------------------------------------------------------------
// Function globalAlignmentScore()
// ----------------------------------------------------------------------------

/**
.Function.globalAlignmentScore
..summary:Computes the best global pairwise alignment score.
..cat:Alignments
..signature:globalAlignmentScore(seqH, seqV, scoringScheme, [alignConfig,] [lowerDiag, upperDiag,] [algorithmTag])
..signature:globalAlignmentScore(strings,    scoringScheme, [alignConfig,] [lowerDiag, upperDiag,] [algorithmTag])
..signature:globalAlignmentScore(seqH, seqV, {MyersBitVector | MyersHirschberg})
..signature:globalAlignmentScore(strings,    {MyersBitVector | MyersHirschberg})
..param.seqH:Horizontal gapped sequence in alignment matrix.
...type:Class.String
..param.seqV:Vertical gapped sequence in alignment matrix.
...type:Class.String
..param.strings:A @Class.StringSet@ containing two sequences.
...type:Class.StringSet
..param.scoringScheme:
The scoring scheme to use for the alignment.
Note that the user is responsible for ensuring that the scoring scheme is compatible with $algorithmTag$.
...type:Class.Score
..param.alignConfig:The @Class.AlignConfig@ to use for the alignment.
...type:Class.AlignConfig
..param.lowerDiag:Optional lower diagonal.
...type:nolink:$int$
..param.upperDiag:Optional upper diagonal.
...type:nolink:$int$
..param.algorithmTag:The Tag for picking the alignment algorithm.
...type:Tag.Pairwise Global Alignment Algorithms.tag.Gotoh
...type:Tag.Pairwise Global Alignment Algorithms.tag.NeedlemanWunsch
...type:Tag.Pairwise Global Alignment Algorithms.tag.Hirschberg
...type:Tag.Pairwise Global Alignment Algorithms.tag.MyersHirschberg
...type:Tag.Pairwise Global Alignment Algorithms.tag.MyersBitVector
..returns:An integer with the alignment score, as given by the @Metafunction.Value@ metafunction of the @Class.Score@ type.
..remarks:
This function does not perform the (linear time) traceback step after the (mostly quadratic time) dynamic programming step.
Note that Myers' bit-vector algorithm does not compute an alignment (only in the Myers-Hirschberg variant) but scores can be computed using $globalAlignmentScore$.
..remarks:
The same limitations to algorithms as in @Function.globalAlignment@ apply.
Furthermore, the $MyersBitVector$ and $MyersHirschberg$ variants can only be used without any other parameter.
..see:Function.globalAlignment
..wiki:Tutorial/PairwiseSequenceAlignment
*/

// ----------------------------------------------------------------------------
// Function globalAlignmentScore()                        [unbanded, 2 Strings]
// ----------------------------------------------------------------------------

template <typename TAlphabetH, typename TSpecH,
          typename TAlphabetV, typename TSpecV,
          typename TScoreValue, typename TScoreSpec,
          bool TOP, bool LEFT, bool RIGHT, bool BOTTOM, typename TACSpec,
          typename TAlgoTag>
TScoreValue globalAlignmentScore(String<TAlphabetH, TSpecH> const & seqH,
                                 String<TAlphabetV, TSpecV> const & seqV,
                                 Score<TScoreValue, TScoreSpec> const & scoringScheme,
                                 AlignConfig<TOP, LEFT, RIGHT, BOTTOM, TACSpec> const & alignConfig,
                                 TAlgoTag const & algoTag)
{
    return _setUpAndRunAlignment(seqH, seqV, scoringScheme, alignConfig, algoTag);
}

// Interface without AlignConfig<>.

template <typename TAlphabetH, typename TSpecH,
          typename TAlphabetV, typename TSpecV,
          typename TScoreValue, typename TScoreSpec,
          typename TAlgoTag>
TScoreValue globalAlignmentScore(String<TAlphabetH, TSpecH> const & seqH,
                                 String<TAlphabetV, TSpecV> const & seqV,
                                 Score<TScoreValue, TScoreSpec> const & scoringScheme,
                                 TAlgoTag const & algoTag)
{
    AlignConfig<> alignConfig;
    return globalAlignmentScore(seqH, seqV, scoringScheme, alignConfig, algoTag);
}

// Interface without algorithm tag.

template <typename TAlphabetH, typename TSpecH,
          typename TAlphabetV, typename TSpecV,
          typename TScoreValue, typename TScoreSpec,
          bool TOP, bool LEFT, bool RIGHT, bool BOTTOM, typename TACSpec>
TScoreValue globalAlignmentScore(String<TAlphabetH, TSpecH> const & seqH,
                                 String<TAlphabetV, TSpecV> const & seqV,
                                 Score<TScoreValue, TScoreSpec> const & scoringScheme,
                                 AlignConfig<TOP, LEFT, RIGHT, BOTTOM, TACSpec> const & alignConfig)
{
    if (scoreGapOpen(scoringScheme) == scoreGapExtend(scoringScheme))
        return globalAlignmentScore(seqH, seqV, scoringScheme, alignConfig, NeedlemanWunsch());
    else
        return globalAlignmentScore(seqH, seqV, scoringScheme, alignConfig, Gotoh());
}

// Interface without AlignConfig<> and algorithm tag.

template <typename TAlphabetH, typename TSpecH,
          typename TAlphabetV, typename TSpecV,
          typename TScoreValue, typename TScoreSpec>
TScoreValue globalAlignmentScore(String<TAlphabetH, TSpecH> const & seqH,
                                 String<TAlphabetV, TSpecV> const & seqV,
                                 Score<TScoreValue, TScoreSpec> const & scoringScheme)
{
    AlignConfig<> alignConfig;
    return globalAlignmentScore(seqH, seqV, scoringScheme, alignConfig);
}

// ----------------------------------------------------------------------------
// Function globalAlignmentScore()                        [unbanded, StringSet]
// ----------------------------------------------------------------------------

template <typename TString, typename TSpec,
          typename TScoreValue, typename TScoreSpec,
          bool TOP, bool LEFT, bool RIGHT, bool BOTTOM, typename TACSpec,
          typename TAlgoTag>
TScoreValue globalAlignmentScore(StringSet<TString, TSpec> const & strings,
                                 Score<TScoreValue, TScoreSpec> const & scoringScheme,
                                 AlignConfig<TOP, LEFT, RIGHT, BOTTOM, TACSpec> const & alignConfig,
                                 TAlgoTag const & algoTag)
{
    SEQAN_ASSERT_EQ(length(strings), 2u);
    return _setUpAndRunAlignment(strings[0], strings[1], scoringScheme, alignConfig, algoTag);
}

// Interface without AlignConfig<>.

template <typename TString, typename TSpec,
          typename TScoreValue, typename TScoreSpec,
          typename TAlgoTag>
TScoreValue globalAlignmentScore(StringSet<TString, TSpec> const & strings,
                                 Score<TScoreValue, TScoreSpec> const & scoringScheme,
                                 TAlgoTag const & algoTag)
{
    SEQAN_ASSERT_EQ(length(strings), 2u);

    AlignConfig<> alignConfig;
    return globalAlignmentScore(strings[0], strings[1], scoringScheme, alignConfig, algoTag);
}

// Interface without algorithm tag.

template <typename TString, typename TSpec,
          typename TScoreValue, typename TScoreSpec,
          bool TOP, bool LEFT, bool RIGHT, bool BOTTOM, typename TACSpec>
TScoreValue globalAlignmentScore(StringSet<TString, TSpec> const & strings,
                                 Score<TScoreValue, TScoreSpec> const & scoringScheme,
                                 AlignConfig<TOP, LEFT, RIGHT, BOTTOM, TACSpec> const & alignConfig)
{
    SEQAN_ASSERT_EQ(length(strings), 2u);

    if (scoreGapOpen(scoringScheme) == scoreGapExtend(scoringScheme))
        return globalAlignmentScore(strings[0], strings[1], scoringScheme, alignConfig, NeedlemanWunsch());
    else
        return globalAlignmentScore(strings[0], strings[1], scoringScheme, alignConfig, Gotoh());
}

// Interface without AlignConfig<> and algorithm tag.

template <typename TString, typename TSpec,
          typename TScoreValue, typename TScoreSpec>
TScoreValue globalAlignmentScore(StringSet<TString, TSpec> const & strings,
                                 Score<TScoreValue, TScoreSpec> const & scoringScheme)
{
    SEQAN_ASSERT_EQ(length(strings), 2u);

    AlignConfig<> alignConfig;
    return globalAlignmentScore(strings[0], strings[1], scoringScheme, alignConfig);
}

}  // namespace seqan

#endif  // #ifndef SEQAN_CORE_INCLUDE_SEQAN_ALIGN_GLOBAL_ALIGNMENT_UNBANDED_H_
