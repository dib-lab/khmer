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

#ifndef SEQAN_INCLUDE_SEQAN_ALIGN_LOCAL_ALIGNMENT_ENUMERATION_H_
#define SEQAN_INCLUDE_SEQAN_ALIGN_LOCAL_ALIGNMENT_ENUMERATION_H_

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
 * @brief Enumeration of local alignments.
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
 *       comparisons. J Mol Biol 1987, 197(4):723-728.</li>
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
 * @section Examples
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
 * resize(rows(align), 2);@s
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
 * @param[in] scheme    The @link Score @endlink object to use for the alignment score.
 * @param[in] cutoff    Alignments with scores <tt>&lt; cutoff</tt> will be discarded (<tt>int</tt>, default 0).
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
 * @section Examples
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
 * @param[in] scheme    The @link Score @endlink object to use for the alignment score.
 * @param[in] upperDiag An <tt>int</tt> with the upper diagonal.
 * @param[in] lowerDiag An <tt>int</tt> with the lower diagonal.
 * @param[in] cutoff    Alignments with scores <tt>&lt; cutoff</tt> will be discarded (<tt>int</tt>, default 0).
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
 * @param[in] enumerator The LocalAlignmentEnumerator to query.
 *
 * @return TScoreVal The current alignment score (@link Score#Value @endlink of <tt>TScore</tt>).
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
 * @param[in] align      @link Align @endlink object to use for the alignment representation.
 * @param[in] gapsH      @link Gaps @endlink object to use for the first/horizontal sequence in the alignment matrix.
 * @param[in] gapsV      @link Gaps @endlink object to use for the second/vertical sequence in the alignment matrix.
 * @param[in] enumerator The LocalAlignmentEnumerator to advance.
 *
 * @return bool <tt>true</tt> if another suboptimal alignment above the given threshold was found and <tt> false
 *              otherwise.
 */

}  // namespace seqan

#endif  // #ifndef SEQAN_INCLUDE_SEQAN_ALIGN_LOCAL_ALIGNMENT_ENUMERATION_H_
