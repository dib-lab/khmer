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

#ifndef SEQAN_CORE_INCLUDE_SEQAN_ALIGN_LOCAL_ALIGNMENT_ENUMERATION_H_
#define SEQAN_CORE_INCLUDE_SEQAN_ALIGN_LOCAL_ALIGNMENT_ENUMERATION_H_

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// ----------------------------------------------------------------------------
// Class LocalAlignmentEnumerator
// ----------------------------------------------------------------------------

template <typename TScore, typename TSpec>
class LocalAlignmentEnumerator;

/*!
 * @class LocalAlignmentEnumerator
 * @headerfile <seqan/align.h>
 *
 * @signature template <typename TScore, typename TSpec>
 *            class LocalAlignmentEnumerator;
 *
 * @tparam TScore The type of the @link Score @endlink to use for the local alignment.
 * @tparam TSpec  The tag to use for specializing the enumerator.
 *
 * See the documentation of the specializations for examples.
 *
 * @section References
 *
 * <ul>
 *   <li>Waterman MS, Eggert M: A new algorithm for best subsequence alignments with application to tRNA-rRNA
 *       comparisons. J Mol Biol 1987, 197(4):723-728.</lI.
 * </ul>
 */

/*!
 * @class UnbandedLocalAlignmentEnumerator
 * @extends LocalAlignmentEnumerator
 * @headerfile <seqan/align.h>
 * @brief Unbanded enumeration of local alignments using the Waterman-Eggert algorithm.
 *
 * @signature template <typename TScore>
 *            class LocalAlignmentEnumerator<TScore, Unbanded>;
 *
 * @tparam TScore The @link Score @endlink type.
 *
 * @section Example
 *
 * Enumerate all alignments into a @link Align @endlink object.
 *
 * @code{.cpp}
 * SimpleScore scoringScheme(2, -1, -1, -2);
 * LocalAlignmentEnumerator<SimpleScore, Unbanded> enumerator(scoringScheme, 5);
 * 
 * Dna5String seqH = "CGAGAGAGACCGAGA";
 * Dna5String seqV = "TTCTGAGATCCGTTTTT";
 * 
 * Align<Dna5String> align;
 * resize(rows(align), 2);
 * assignSource(row(align), 0, seqH);
 * assignSource(row(align), 1, seqV);
 * 
 * int i = 0;
 * while (nextLocalAlignment(align, enumerator))
 * {
 *     std::cout << i << "-th alignment:\n";
 *     std::cout << align << "\n\n";
 *     std::cout << "score == " << getScore(enumerator) << "\n";
 * }
 * @endcode
 */

/*!
 * @fn UnbandedLocalAlignmentEnumerator::LocalAlignmentEnumerator
 * @brief Constructor.
 *
 * @signature LocalAlignmentEnumerator::LocalAlignmentEnumerator(scheme[, cutoff]);
 *
 * @param scheme    The @link Score @endlink object to use for the alignment score.
 * @param cutoff    Alignments with scores <tt>&lt; cutoff</tt> will be discarded (<tt>int</tt>, default 0).
 */

/*!
 * @class BandedLocalAlignmentEnumerator
 * @extends LocalAlignmentEnumerator
 * @headerfile <seqan/align.h>
 * @brief Banded enumeration of local alignments using the Waterman-Eggert algorithm.
 *
 * @signature template <typename TScore>
 *            class LocalAlignmentEnumerator<TScore, Banded>;
 *
 * @tparam TScore The @link Score @endlink type.
 *
 * @section Example
 *
 * Enumerate all alignments in the band between -3 and 0 into an @link Align @endlink object.
 *
 * @code{.cpp}
 * SimpleScore scoringScheme(2, -1, -1, -2);
 * LocalAlignmentEnumerator<SimpleScore, Banded> enumerator(scoringScheme, -3, 0, 5);
 * 
 * Dna5String seqH = "CGAGAGAGACCGAGA";
 * Dna5String seqV = "TTCTGAGATCCGTTTTT";
 * 
 * Align<Dna5String> align;
 * resize(rows(align), 2);
 * assignSource(row(align), 0, seqH);
 * assignSource(row(align), 1, seqV);
 * 
 * int i = 0;
 * while (nextLocalAlignment(align, enumerator))
 * {
 *     std::cout << i << "-th alignment:\n";
 *     std::cout << align << "\n\n";
 *     std::cout << "score == " << getScore(enumerator) << "\n";
 * }
 * @endcode
 */

/*!
 * @fn BandedLocalAlignmentEnumerator::LocalAlignmentEnumerator
 * @brief Constructor.
 *
 * @signature LocalAlignmentEnumerator::LocalAlignmentEnumerator(scheme, upperDiag, lowerDiag[, cutoff]);
 *
 * @param scheme    The @link Score @endlink object to use for the alignment score.
 * @param upperDiag An <tt>int</tt> with the upper diagonal.
 * @param lowerDiag An <tt>int</tt> with the lower diagonal.
 * @param cutoff    Alignments with scores <tt>&lt; cutoff</tt> will be discarded (<tt>int</tt>, default 0).
 */

/**
.Class.LocalAlignmentEnumerator
..cat:Alignments
..summary:Enumerate local alignments using the Waterman-Eggert algorithm.
..description:This is an abstract base class for the alignment enumeration; the specializations provide the actual implementaiton of banded and unbanded local aligment search.
..signature:LocalAlignmentEnumerator<TScore, TSpec>
..param.TScore:The @Class.Score@ type to use.
...type:Class.Score
..param.TSpec:Specialization tag.
..example.text:See the specializations for usage examples.
..cite:Waterman MS, Eggert M: A new algorithm for best subsequence alignments with application to tRNA-rRNA comparisons. J Mol Biol 1987, 197(4):723-728.
..include:seqan/align.h

.Spec.Unbanded LocalAlignmentEnumerator
..cat:Alignments
..general:Class.LocalAlignmentEnumerator
..summary:Unbanded enumeration of local alignments using the Waterman-Eggert algorithm.
..signature:LocalAlignmentEnumerator<TScore, Unbanded>
..example.text:Enumerate all alignments into an @Class.Align@ object.
..example.code:
SimpleScore scoringScheme(2, -1, -1, -2);
LocalAlignmentEnumerator<SimpleScore, Unbanded> enumerator(scoringScheme, 5);

Dna5String seqH = "CGAGAGAGACCGAGA";
Dna5String seqV = "TTCTGAGATCCGTTTTT";

Align<Dna5String> align;
resize(rows(align), 2);
assignSource(row(align), 0, seqH);
assignSource(row(align), 1, seqV);

int i = 0;
while (nextLocalAlignment(align, enumerator))
{
    std::cout << i << "-th alignment:\n";
    std::cout << align << "\n\n";
    std::cout << "score == " << getScore(enumerator) << "\n";
}
..include:seqan/align.h

.Memfunc.Unbanded LocalAlignmentEnumerator#LocalAlignmentEnumerator
..class:Spec.Unbanded LocalAlignmentEnumerator
..summary:Constructor
..signature:LocalAlignmentEnumerator(score, [cutoff])
..param.score:The scoring scheme to use for the alignments.
...type:Class.Score
..param.cutoff:Alignments with scores < $cutoff$ will be discarded.
...default:0
...type:nolink:$int$

.Spec.Banded LocalAlignmentEnumerator
..cat:Alignments
..general:Class.LocalAlignmentEnumerator
..signature:LocalAlignmentEnumerator<TScore, Banded>
..summary:Banded enumeration of local alignments using the Waterman-Eggert algorithm.
..example.text:Enumerate all alignments in the band between $-3$ and $0$ into an @Class.Align@ object.
..example.code:
SimpleScore scoringScheme(2, -1, -1, -2);
LocalAlignmentEnumerator<SimpleScore, Banded> enumerator(scoringScheme, 5, -3, 0);

Dna5String seqH = "CGAGAGAGACCGAGA";
Dna5String seqV = "TTCTGAGATCCGTTTTT";

Align<Dna5String> align;
resize(rows(align), 2);
assignSource(row(align), 0, seqH);
assignSource(row(align), 1, seqV);

int i = 0;
while (nextLocalAlignment(align, enumerator))
{
    std::cout << i << "-th alignment:\n";
    std::cout << align << "\n\n";
    std::cout << "score == " << getScore(enumerator) << "\n";
}
..include:seqan/align.h

.Memfunc.Banded LocalAlignmentEnumerator#LocalAlignmentEnumerator
..class:Spec.Banded LocalAlignmentEnumerator
..summary:Constructor
..signature:LocalAlignmentEnumerator(score, upperDiag, lowerDiag, [cutoff])
..param.score:The scoring scheme to use for the alignments.
...type:Class.Score
..param.upperDiag:Upper diagonal of the band.
...type:nolink:$int$
..param.lowerDiag:Lower diagonal of the band.
...type:nolink:$int$
..param.cutoff:Alignments with scores < $cutoff$ will be discarded.
...type:nolink:$int$
...default:0
*/

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function getScore()
// ----------------------------------------------------------------------------

/*!
 * @fn LocalAlignmentEnumerator#getScore
 * @headerfile <seqan/align.h>
 * @brief Get current alignment score.
 *
 * @signature TScoreVal getScore(enumerator);
 *
 * @param enumerator The LocalAlignmentEnumerator to query.
 * 
 * @return TScoreVal The current alignment score.
 */

/**
.Function.LocalAlignmentEnumerator#getScore
..cat:Alignments
..summary:Compute next suboptimal local alignment.
..signature:getScore(enumerator)
..param.enumerator:The local alignment enumerator to use.
...type:Class.LocalAlignmentEnumerator
..returns:
The score of the previously computed alignment.
(Type: @Metafunction.Value@ of $enumerator$'s class.)
..include:seqan/align.h
 */

// ----------------------------------------------------------------------------
// Function nextLocalAlignment()
// ----------------------------------------------------------------------------

/*!
 * @fn LocalAlignmentEnumerator#nextLocalAlignment
 * @headerfile <seqan/align.h>
 * @brief Compute next suboptimal local alignment.
 *
 * @signature bool nextLocalAlignment(align,        enumerator);
 * @signature bool nextLocalAlignment(gapsH, gapsV, enumerator);
 *
 * @param align      @link Align @endlink object to use for the alignment representation.
 * @param gapsH      @link Gaps @endlink object to use for the first/horizontal sequence in the alignment matrix.
 * @param gapsV      @link Gaps @endlink object to use for the second/vertical sequence in the alignment matrix.
 * @param enumerator The LocalAlignmentEnumerator to advance.
 * 
 * @return bool <tt>true</tt> if another suboptimal alignment above the given threshold was found and <tt> false
 *              otherwise.
 */

/**
.Function.nextLocalAlignment
..cat:Alignments
..summary:Compute next suboptimal local alignment.
..signature:nextLocalAlignment(align,        enumerator)
..signature:nextLocalAlignment(gapsH, gapsV, enumerator)
..param.align:The @Class.Align@ object to use for the alignment representation.
...type:Class.Align
..param.gapsH:The @Class.Gaps@ object to use for the horizontal sequence in the alignment matrix.
...type:Class.Gaps
..param.gapsV:The @Class.Gaps@ object to use for the vertical sequence in the alignment matrix.
...type:Class.Gaps
..param.enumerator:The @Class.LocalAlignmentEnumerator@ object to use.
...type:Class.LocalAlignmentEnumerator
..returns:$true$ if another suboptimal alignment above the given threshold was found, $false$ otherwise.
...type:nolink:$bool$
..include:seqan/align.h
*/

}  // namespace seqan

#endif  // #ifndef SEQAN_CORE_INCLUDE_SEQAN_ALIGN_LOCAL_ALIGNMENT_ENUMERATION_H_
