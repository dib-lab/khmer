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
// Author: Rene Rahn <rene.rahn@fu-berlin.de>
// ==========================================================================
// The global interfaces to call to compute the banded chain alignment
// algorithm.
// ==========================================================================

#ifndef INCLUDE_SEQAN_SEEDS_BANDED_CHAIN_ALIGNMENT_H_
#define INCLUDE_SEQAN_SEEDS_BANDED_CHAIN_ALIGNMENT_H_

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

/*!
 * @fn bandedChainAlignment
 * @headerfile <seqan/seeds.h>
 * @brief Computes the best global pairwise alignment between two sequences given a non-empty seed chain.
 *
 * @signature TValue bandedChainAlignment(align,          seedChain, scoringScheme1[, scoringScheme2] [, alignConfig] [, k]);
 * @signature TValue bandedChainAlignment(gapsH, gapsV,   seedChain, scoringScheme1[, scoringScheme2] [, alignConfig] [, k]);
 * @signature TValue bandedChainAlignment(frags, strings, seedChain, scoringScheme1[, scoringScheme2] [, alignConfig] [, k]);
 * @signature TValue bandedChainAlignment(alignmentGraph, seedChain, scoringScheme1[, scoringScheme2] [, alignConfig] [, k]);
 *
 * @param[in,out] align   An @link Align @endlink object that stores the alignment. The number of rows must be 2 and the
 *                        sequences must have already been set.  <tt>row(align, 0)</tt> is the horizontal sequence in the
 *                        alignment matrix, <tt>row(align, 1)</tt> is the vertical sequence.
 * @param[in,out] gapsH   Horizontal gapped sequence in alignment matrix. Type: @link Gaps @endlink.
 * @param[in,out] gapsV   Vertical gapped sequence in alignment matrix. Type: @link Gaps @endlink.
 * @param[out]    frags   @link String @endlink of @link Fragment @endlink objects.  The sequence with id <tt>0</tt>
 *                        is the horizontal one, the sequence with id <tt>1</tt> is the vertical one.
 * @param[in]     strings A @link StringSet @endlink containing two sequences.
 * @param[in,out] alignmentGraph
 *                       @link AlignmentGraph @endlink object to store the alignment in.  The underlying @link
 *                       StringSet @endlink must be an @link DependentStringSet @endlink. Types: AlignmentGraph.
 * @param[in] seedChain  The container holding the @link Seed seeds @endlink.  Note that the @link Seed seeds @endlink
 *                       have to be in montonic non-decreasing order and the container has to implement a forward-iterator.
 *                       Type: @link SeedSet @endlink.
 * @param[in] scoringScheme1
 *                       The scoring scheme used for the alignment. If <tt>scoringScheme2</tt> is specified, then
 *                       <tt>scoringScheme1</tt> is used for the regions around the seeds and <tt>scoringScheme2</tt>
 *                       for the gap regions between two consecutive seeds.  Types: Score
 * @param[in] scoringScheme2
 *                       The optional scoring scheme for the gap regions between two anchors. Types: Score
 * @param[in] k          Optional extension of the band around the seeds.  At the moment only band extensions greater
 *                       or equal <tt>1</tt> are allowed. Type: nolink:<tt>int</tt>. Default: 15.
 * @param[in] alignConfig
 *                       The @link AlignConfig @endlink to use for the alignment.
 *
 * @return TValue An integer with the alignment score, as given by the @link Score#Value @endlink metafunction of the
 *                @link Score @endlink type.  If the seed chain is empty then the @link OrderedAlphabetConcept#MinValue
 *                smallest value of the score type @endlink is used to return the minimal value of the selected score
 *                type and no alignment is computed.
 *
 * There exist multiple overloads for this function with four configuration dimensions.
 *
 * First, you can select whether begin and end gaps are free in either sequence using <tt>alignConfig</tt>.
 *
 * Second, you can select the type of the target storing the alignment. This can be either an @link Align @endlink
 * object, two @link Gaps @endlink objects, a @link AlignmentGraph @endlink, or a string of @link Fragment @endlink
 * objects. @link Align @endlink objects provide an interface to tabular alignments with the restriction of all rows
 * having the same type. Using two @link Gaps @endlink objects has the advantage that you can align sequences with
 * different types, for example @link DnaString @endlink and @link Dna5String @endlink. @link AlignmentGraph Alignment
 * Graphs @endlink provide a graph-based representation of segment-based colinear alignments. Using @link Fragment
 * @endlink strings is useful for collecting many pairwise alignments, for example in the construction of @link
 * AlignmentGraph Alignment Graphs @endlink for multiple-sequence alignments (MSA).
 *
 * Third, you can optionally give a second scoring scheme to fill the gaps between two consecutive seeds. Note that
 * based on the specified scores either an affine or linear gap cost function is used. This only depends on whether for
 * one of the scoring schemes the scores for gap opening and gap extension differ or not. If only one scoring scheme is
 * defined the complete region is computed with the same scoring scheme.
 *
 * Fourth, you can optinally select a proper band extension for the bands around the seeds.  At the moment only band
 * extensions of at least <tt>1</tt> are allowed.  The default value is <tt>15</tt> and is based on the default values
 * for the LAGAN-algorithm described by Brudno et al., 2003.
 *
 * The examples below show some common use cases.
 *
 * @section Examples
 *
 * Banded chain alignment of two sequences using an @link Align @endlink object and using only one scoring scheme and no
 * free end-gaps.
 *
 * @code{.cpp}
 * Dna5String seqH = "CGAATCCATCCCACACA";
 * Dna5String seqV = "GGCGATNNNCATGGCACA";
 *
 * String<Seed<Simple> > seedChain;
 * appendValue(seedChain, Seed<Simple>(2, 0, 6, 5));
 * appendValue(seedChain, Seed<Simple>(9, 6, 12, 9));
 * appendValue(seedChain, Seed<Simple>(14, 11, 16, 17));
 *
 * Align<Dna5String, ArrayGaps> alignment;
 * resize(rows(alignment), 2);
 * assignSource(row(alignment, 0), seqH);
 * assignSource(row(alignment, 1), seqV);
 *
 * Score<int, Simple> scoringScheme(2, -1, -2);
 *
 * int result = bandedChainAlignment(alignment, seedChain, scoringScheme, 2);
 * @endcode
 *
 * Banded chain alignment of two sequences using two @link Gaps @endlink objects, an unordered seed set to hold the
 * seeds, two different scoring schemes for the gaps between the seeds and the seeds and free end-gaps.
 *
 * @code{.cpp}
 * DnaString seqH = "CGAATCCATCCCACACA";
 * Dna5String seqV = "GGCGATNNNCATGGCACA";
 *
 * SeedSet<Simple, Unordered> seedChain;
 * addSeed(seedChain, Seed<Simple>(2, 0, 6, 5), Single());
 * addSeed(seedChain, Seed<Simple>(9, 6, 12, 9), Single());
 * addSeed(seedChain, Seed<Simple>(14, 11, 16, 17), Single());
 *
 * Gaps<DnaString, ArrayGaps> gapsH(seqH);
 * Gaps<Dna5String, AnchorGaps<> > gapsV(seqV);
 *
 * Score<int, Simple> scoringSchemeSeed(2, -1, -2);
 * Score<int, Simple> scoringSchemeGap(5, -3, -1, -5);
 * AlignConfig<true, true, true, true> alignConfig;
 *
 * int result = globalAlignment(gapsH, gapsV, scoringSchemeSeed, scoringSchemeGap, alignConfig, 2);
 * @endcode
 *
 * @section Tutorial
 *
 * Also see the <a href="http://seqan.readthedocs.org/en/develop/Tutorial/SeedAndExtend.html">Seed-and-Extend Tutorial</a>
 *
 * @section Reference
 *
 * <ul>
 *   <li>Brudno M, Do CB, Cooper GM, et al.: LAGAN and Multi-LAGAN: Efficient Tools for Large-Scale Multiple Alignment
 *       of Genomic DNA. Genome Research 2003, 13: 721-731.</li>
 * </ul>
 */

// ----------------------------------------------------------------------------
// Function bandedChainAlignment()                                      [Align]
// ----------------------------------------------------------------------------

template <typename TSequence, typename TAlignSpec, typename TSeeds, typename TScoreValue, typename TScoreSpecAnchor,
          typename TScoreSpecGap, bool TFirstRow, bool TFirstColumn, bool TLastColumn, bool TLastRow, typename TACSpec>
inline TScoreValue
bandedChainAlignment(Align<TSequence, TAlignSpec> & align,
                     TSeeds const & seedSet,
                     Score<TScoreValue, TScoreSpecAnchor> const & scoreSchemeAnchor,
                     Score<TScoreValue, TScoreSpecGap> const & scoreSchemeGap,
                     AlignConfig<TFirstRow, TFirstColumn, TLastColumn, TLastRow, TACSpec> const & alignConfig,
                     unsigned bandExtension = 15)
{
    typedef typename Position<TSequence>::Type TPosition;
    typedef typename Size<TSequence>::Type TSize;
    typedef StringSet<String<TraceSegment_<TPosition, TSize> > > TTraceSegmentSet;

    TTraceSegmentSet traceSet;
    TScoreValue score =
        _setupAndRunBandedChainAlignment(traceSet, seedSet, source(row(align, 0)), source(row(align, 1)),
                                         scoreSchemeAnchor, scoreSchemeGap, alignConfig, bandExtension, GapsLeft());

    if (empty(traceSet))
        return score;

    _adaptTraceSegmentsTo(row(align,0), row(align,1), value(traceSet, 0));
    return score;
}

// With only one scoring scheme.
template <typename TSequence, typename TAlignSpec, typename TSeeds, typename TScoreValue, typename TScoreSpec,
          bool TFirstRow, bool TFirstColumn, bool TLastColumn, bool TLastRow, typename TACSpec>
inline TScoreValue
bandedChainAlignment(Align<TSequence, TAlignSpec> & align,
                     TSeeds const & seedSet,
                     Score<TScoreValue, TScoreSpec> const & scoreScheme,
                     AlignConfig<TFirstRow, TFirstColumn, TLastColumn, TLastRow, TACSpec> const & alignConfig,
                     unsigned bandExtension = 15)
{
    return bandedChainAlignment(align, seedSet, scoreScheme, scoreScheme, alignConfig, bandExtension);
}

// Without AlignConfig.
template <typename TSequence, typename TAlignSpec, typename TSeeds, typename TScoreValue, typename TScoreSpecAnchor,
          typename TScoreSpecGap>
inline TScoreValue
bandedChainAlignment(Align<TSequence, TAlignSpec> & align,
                     TSeeds const & seedSet,
                     Score<TScoreValue, TScoreSpecAnchor> const & scoreSchemeAnchor,
                     Score<TScoreValue, TScoreSpecGap> const & scoreSchemeGap,
                     unsigned bandExtension = 15)
{
    return bandedChainAlignment(align, seedSet, scoreSchemeAnchor, scoreSchemeGap, AlignConfig<>(), bandExtension);
}

// Without AlignConfig and with only one scoring scheme.
template <typename TSequence, typename TAlignSpec, typename TSeeds, typename TScoreValue, typename TScoreSpec>
inline TScoreValue
bandedChainAlignment(Align<TSequence, TAlignSpec> & align,
                     TSeeds const & seedSet,
                     Score<TScoreValue, TScoreSpec> const & scoreScheme,
                     unsigned bandExtension = 15)
{
    return bandedChainAlignment(align, seedSet, scoreScheme, scoreScheme, AlignConfig<>(), bandExtension);
}

// ----------------------------------------------------------------------------
// Function bandedChainAlignment()                                       [Gaps]
// ----------------------------------------------------------------------------

template <typename TSequenceH, typename TGapSpecH, typename TSequenceV, typename TGapSpecV, typename TSeeds,
          typename TScoreValue, typename TScoreSpecAnchor, typename TScoreSpecGap, bool TFirstRow, bool TFirstColumn,
          bool TLastColumn, bool TLastRow, typename TACSpec>
inline TScoreValue
bandedChainAlignment(Gaps<TSequenceH, TGapSpecH> & gapsHorizontal,
                     Gaps<TSequenceV, TGapSpecV> & gapsVertical,
                     TSeeds const & seedSet,
                     Score<TScoreValue, TScoreSpecAnchor> const & scoreSchemeAnchors,
                     Score<TScoreValue, TScoreSpecGap> const & scoreSchemeGaps,
                     AlignConfig<TFirstRow, TFirstColumn, TLastColumn, TLastRow, TACSpec> const & alignConfig,
                     unsigned bandExtension = 15)
{
    typedef typename Position<TSequenceH>::Type TPosition;
    typedef typename Size<TSequenceH>::Type TSize;
    typedef StringSet<String<TraceSegment_<TPosition, TSize> > > TTraceSegmentSet;

    TTraceSegmentSet traceSet;
    TScoreValue score =
        _setupAndRunBandedChainAlignment(traceSet, seedSet, source(gapsHorizontal), source(gapsVertical),
                                         scoreSchemeAnchors, scoreSchemeGaps, alignConfig, bandExtension, GapsLeft());

    if (empty(traceSet))
        return score;

    _adaptTraceSegmentsTo(gapsHorizontal, gapsVertical, value(traceSet, 0));
    return score;
}

// With only one scoring scheme.
template <typename TSequenceH, typename TGapSpecH, typename TSequenceV, typename TGapSpecV, typename TSeeds,
          typename TScoreValue, typename TScoreSpec, bool TFirstRow, bool TFirstColumn, bool TLastColumn, bool TLastRow,
          typename TACSpec>
inline TScoreValue
bandedChainAlignment(Gaps<TSequenceH, TGapSpecH> & gapsHorizontal,
                     Gaps<TSequenceV, TGapSpecV> & gapsVertical,
                     TSeeds const & seedSet,
                     Score<TScoreValue, TScoreSpec> const & scoreScheme,
                     AlignConfig<TFirstRow, TFirstColumn, TLastColumn, TLastRow, TACSpec> const & alignConfig,
                     unsigned bandExtension = 15)
{
    return bandedChainAlignment(gapsHorizontal, gapsVertical, seedSet, scoreScheme, scoreScheme, alignConfig,
                                bandExtension);
}

// Without AlignConfig.
template <typename TSequenceH, typename TGapSpecH, typename TSequenceV, typename TGapSpecV, typename TSeeds,
          typename TScoreValue, typename TScoreSpecAnchor, typename TScoreSpecGap>
inline TScoreValue
bandedChainAlignment(Gaps<TSequenceH, TGapSpecH> & gapsHorizontal,
                     Gaps<TSequenceV, TGapSpecV> & gapsVertical,
                     TSeeds const & seedSet,
                     Score<TScoreValue, TScoreSpecAnchor> const & scoreSchemeAnchor,
                     Score<TScoreValue, TScoreSpecGap> const & scoreSchemeGap,
                     unsigned bandExtension = 15)
{
    return bandedChainAlignment(gapsHorizontal, gapsVertical, seedSet, scoreSchemeAnchor, scoreSchemeGap,
                                AlignConfig<>(), bandExtension);
}

// Without AlignConfig and with only one scoring scheme.
template <typename TSequenceH, typename TGapSpecH, typename TSequenceV, typename TGapSpecV, typename TSeeds,
          typename TScoreValue, typename TScoreSpec>
inline TScoreValue
bandedChainAlignment(Gaps<TSequenceH, TGapSpecH> & gapsHorizontal,
                     Gaps<TSequenceV, TGapSpecV> & gapsVertical,
                     TSeeds const & seedSet,
                     Score<TScoreValue, TScoreSpec> const & scoreScheme,
                     unsigned bandExtension = 15)
{
    return bandedChainAlignment(gapsHorizontal, gapsVertical, seedSet, scoreScheme, scoreScheme, AlignConfig<>(),
                                bandExtension);
}

// ----------------------------------------------------------------------------
// Function bandedChainAlignment()                        [Graph<Alignment<> >]
// ----------------------------------------------------------------------------

template <typename TStringSet, typename TCargo, typename TGraphSpec, typename TSeeds, typename TScoreValue,
          typename TScoreSpecAnchor, typename TScoreSpecGap, bool TFirstRow, bool TFirstColumn, bool TLastColumn,
          bool TLastRow, typename TACSpec>
inline TScoreValue
bandedChainAlignment(Graph<Alignment<TStringSet, TCargo, TGraphSpec> > & alignmentGraph,
                     TSeeds const & seedSet,
                     Score<TScoreValue, TScoreSpecAnchor> const & scoreSchemeAnchors,
                     Score<TScoreValue, TScoreSpecGap> const & scoreSchemeGaps,
                     AlignConfig<TFirstRow, TFirstColumn, TLastColumn, TLastRow, TACSpec> const & alignConfig,
                     unsigned bandExtension = 15)
{
    typedef typename Position<TStringSet>::Type TPosition;
    typedef typename Size<TStringSet>::Type TSize;
    typedef StringSet<String<TraceSegment_<TPosition, TSize> > > TTraceSegmentSet;

    TTraceSegmentSet traceSet;
    TScoreValue score =
        _setupAndRunBandedChainAlignment(traceSet, seedSet, value(stringSet(alignmentGraph), 0),
                                         value(stringSet(alignmentGraph), 1), scoreSchemeAnchors, scoreSchemeGaps,
                                         alignConfig, bandExtension, GapsLeft());

    if (empty(traceSet))
        return score;

    _adaptTraceSegmentsTo(alignmentGraph, positionToId(stringSet(alignmentGraph), 0),
                          positionToId(stringSet(alignmentGraph), 1), value(traceSet, 0));
    return score;
}

// With only one scoring scheme.
template <typename TStringSet, typename TCargo, typename TGraphSpec, typename TSeeds, typename TScoreValue,
          typename TScoreSpec, bool TFirstRow, bool TFirstColumn, bool TLastColumn,
          bool TLastRow, typename TACSpec>
inline TScoreValue
bandedChainAlignment(Graph<Alignment<TStringSet, TCargo, TGraphSpec> > & alignmentGraph,
                     TSeeds const & seedSet,
                     Score<TScoreValue, TScoreSpec> const & scoreScheme,
                     AlignConfig<TFirstRow, TFirstColumn, TLastColumn, TLastRow, TACSpec> const & alignConfig,
                     unsigned bandExtension = 15)
{
    return bandedChainAlignment(alignmentGraph, seedSet, scoreScheme, scoreScheme, alignConfig, bandExtension);
}

// Without AlignConfig.
template <typename TStringSet, typename TCargo, typename TGraphSpec, typename TSeeds, typename TScoreValue,
          typename TScoreSpecAnchor, typename TScoreSpecGap>
inline TScoreValue
bandedChainAlignment(Graph<Alignment<TStringSet, TCargo, TGraphSpec> > & alignmentGraph,
                     TSeeds const & seedSet,
                     Score<TScoreValue, TScoreSpecAnchor> const & scoreSchemeAnchor,
                     Score<TScoreValue, TScoreSpecGap> const & scoreSchemeGap,
                     unsigned bandExtension = 15)
{
    return bandedChainAlignment(alignmentGraph, seedSet, scoreSchemeAnchor, scoreSchemeGap, AlignConfig<>(),
                                bandExtension);
}

// Without AlignConfig and with only one scoring scheme.
template <typename TStringSet, typename TCargo, typename TGraphSpec, typename TSeeds, typename TScoreValue,
          typename TScoreSpec>
inline TScoreValue
bandedChainAlignment(Graph<Alignment<TStringSet, TCargo, TGraphSpec> > & alignmentGraph,
                     TSeeds const & seedSet,
                     Score<TScoreValue, TScoreSpec> const & scoreScheme,
                     unsigned bandExtension = 15)
{
    return bandedChainAlignment(alignmentGraph, seedSet, scoreScheme, scoreScheme, AlignConfig<>(), bandExtension);
}

// ----------------------------------------------------------------------------
// Function bandedChainAlignment()                          [String<Fragments>]
// ----------------------------------------------------------------------------

template <typename TSize, typename TFragmentSpec, typename TStringSpec, typename TSequence, typename TStringSetSpec,
          typename TSeeds, typename TScoreValue, typename TScoreSpecAnchor, typename TScoreSpecGap,
          bool TFirstRow, bool TFirstColumn, bool TLastColumn, bool TLastRow, typename TACSpec>
inline TScoreValue
bandedChainAlignment(String<Fragment<TSize, TFragmentSpec>, TStringSpec> & fragmentString,
                     StringSet<TSequence, TStringSetSpec> const & strings,
                     TSeeds const & seedSet,
                     Score<TScoreValue, TScoreSpecAnchor> const & scoreSchemeAnchors,
                     Score<TScoreValue, TScoreSpecGap> const & scoreSchemeGaps,
                     AlignConfig<TFirstRow, TFirstColumn, TLastColumn, TLastRow, TACSpec> const & alignConfig,
                     unsigned bandExtension = 15)
{
    typedef typename Position<TSequence>::Type TPosition;
    typedef StringSet<String<TraceSegment_<TPosition, TSize> > > TTraceSegmentSet;

    TTraceSegmentSet traceSet;
    TScoreValue score =
        _setupAndRunBandedChainAlignment(traceSet, seedSet, value(strings, 0), value(strings, 1), scoreSchemeAnchors,
                                         scoreSchemeGaps, alignConfig, bandExtension, GapsLeft());
    if (empty(traceSet))
        return score;

    _adaptTraceSegmentsTo(fragmentString, positionToId(strings, 0), positionToId(strings, 1), value(traceSet, 0));
    return score;
}

// With only one scoring scheme.
template <typename TSize, typename TFragmentSpec, typename TStringSpec, typename TSequence, typename TStringSetSpec,
          typename TSeeds, typename TScoreValue, typename TScoreSpec, bool TFirstRow, bool TFirstColumn,
          bool TLastColumn, bool TLastRow, typename TACSpec>
inline TScoreValue
bandedChainAlignment(String<Fragment<TSize, TFragmentSpec>, TStringSpec> & fragmentString,
                     StringSet<TSequence, TStringSetSpec> const & strings,
                     TSeeds const & seedSet,
                     Score<TScoreValue, TScoreSpec> const & scoreScheme,
                     AlignConfig<TFirstRow, TFirstColumn, TLastColumn, TLastRow, TACSpec> const & alignConfig,
                     unsigned bandExtension = 15)
{
    return bandedChainAlignment(fragmentString, strings, seedSet, scoreScheme, scoreScheme, alignConfig, bandExtension);
}

// Without AlignConfig.
template <typename TSize, typename TFragmentSpec, typename TStringSpec, typename TSequence, typename TStringSetSpec,
          typename TSeeds, typename TScoreValue, typename TScoreSpecAnchor, typename TScoreSpecGap>
inline TScoreValue
bandedChainAlignment(String<Fragment<TSize, TFragmentSpec>, TStringSpec> & fragmentString,
                     StringSet<TSequence, TStringSetSpec> const & strings,
                     TSeeds const & seedSet,
                     Score<TScoreValue, TScoreSpecAnchor> const & scoreSchemeAnchor,
                     Score<TScoreValue, TScoreSpecGap> const & scoreSchemeGap,
                     unsigned bandExtension = 15)
{
    return bandedChainAlignment(fragmentString, strings, seedSet, scoreSchemeAnchor, scoreSchemeGap, AlignConfig<>(),
                                bandExtension);
}

// Without AlignConfig and with only one scoring scheme.
template <typename TSize, typename TFragmentSpec, typename TStringSpec, typename TSequence, typename TStringSetSpec,
          typename TSeeds, typename TScoreValue, typename TScoreSpec>
inline TScoreValue
bandedChainAlignment(String<Fragment<TSize, TFragmentSpec>, TStringSpec> & fragmentString,
                     StringSet<TSequence, TStringSetSpec> const & strings,
                     TSeeds const & seedSet,
                     Score<TScoreValue, TScoreSpec> const & scoreScheme,
                     unsigned bandExtension = 15)
{
    return bandedChainAlignment(fragmentString, strings, seedSet, scoreScheme, scoreScheme, AlignConfig<>(),
                                bandExtension);
}

}  // namespace seqan

#endif  // #ifndef INCLUDE_SEQAN_SEEDS_BANDED_CHAIN_ALIGNMENT_H_
