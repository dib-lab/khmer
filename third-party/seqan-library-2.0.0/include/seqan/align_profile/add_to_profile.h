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
// Algorithm for profile-to-sequence alignment that then adds the resulting
// sequence alignment to the profile.
// ==========================================================================

#ifndef SEQAN_INCLUDE_SEQAN_ALIGN_PROFILE_ADD_TO_PROFILE_H_
#define SEQAN_INCLUDE_SEQAN_ALIGN_PROFILE_ADD_TO_PROFILE_H_

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
// Function addToProfile()
// ----------------------------------------------------------------------------

/*!
 * @fn addToProfile
 * @headerfile <seqan/align_profile.h>
 * @brief Sequence-to-profile alignment with sequence integration.
 *
 * <tt>seq</tt> will be aligned to <tt>profile</tt> using @link globalAlignment @endlink. To compute a banded alignment
 * <tt>lDiag</tt> and <tt>uDiag</tt> have to be given.
 * The resulting alignment will then be integrated into <tt>profile</tt> by incrementing the <tt>count</tt> values of
 * the <tt>profile</tt> entries that correspond to the <tt>seq</tt> row of the alignment.
 *
 * @signature void addToProfile(profile, seq[, lDiag, uDiag]);
 *
 * @param[in,out] profile The profile to add <tt>seq</tt> to. A @link String @endlink over @link ProfileChar @endlink.
 * @param[in]     seq     @link String @endlink to align to <tt>profile</tt>.
 * @param[in]     lDiag   Lower diagonal for alignment (@link SignedIntegerConcept @endlink).
 * @param[in]     uDiag   Upper diagonal for alignment (@link SignedIntegerConcept @endlink).
 */

// TODO(holtgrew): Add example?

template <typename TChar, typename TValue, typename TSpec, typename TSpec2>
void addToProfile(String<ProfileChar<TChar, TValue, TSpec> > & profile,
                  String<TChar, TSpec2> /*const*/ & seq,
                  int lDiag = minValue<int>(),
                  int uDiag = maxValue<int>())  // non-const because of holder issues
{
    typedef ProfileChar<TChar, TValue, TSpec> TProfileChar;

    typedef String<TProfileChar> TProfileString;
    typedef String<TChar, TSpec2> TSequence;

    // Define gaps and scoring scheme.
    Gaps<TProfileString> gapsH(profile);
    Gaps<TSequence> gapsV(seq);
    seqan::Score<int, seqan::ProfileSeqScore> sScheme(profile);

    // Perform the global alignment.
    if (lDiag == minValue<int>() || uDiag == maxValue<int>())
        globalAlignment(gapsH, gapsV, sScheme, Gotoh());
    else
        globalAlignment(gapsH, gapsV, sScheme, lDiag, uDiag, Gotoh());

    // Construct a new profile from the alignment into buffer and finally swap out the new profile.
    TProfileString buffer;
    reserve(buffer, length(gapsH));
    typename Iterator<Gaps<TProfileString>, Standard>::Type
            itH = begin(gapsH, Standard()),
            itHEnd = end(gapsH, Standard());
    typename Iterator<Gaps<TSequence>, Standard>::Type itV = begin(gapsV, Standard());

    // std::cout << "--- BEFORE ----\n";
    // std::cout << "HORIZONTAL\n";
    // for (; itH != itHEnd; ++itH)
    // {
    //     if (isGap(itH))
    //     {
    //         std::cerr << "-\n";
    //     }
    //     else
    //     {
    //         TProfileChar c = (TProfileChar)*itH;
    //         std::cout << "(" << c.count[0] << ", " << c.count[1] << ", " << c.count[2] << ", " << c.count[3] << ", " << c.count[4] << ")\n";
    //     }
    //}
    //std::cout << "VERTICAL\n";
    //std::cout << gapsV << "\n";

    itH = begin(gapsH, Standard());
    SEQAN_ASSERT_EQ(length(gapsH), length(gapsV));
    for (; itH != itHEnd; ++itH, ++itV)
    {
        if (isGap(itH))
        {
            SEQAN_ASSERT_NOT_MSG(isGap(itV), "Must not generate gaps columns!");
            appendValue(buffer, TProfileChar());
            back(buffer).count[valueSize<TChar>()] += 1;
            back(buffer).count[ordValue(TChar(*itV))] += 1;
        }
        else
        {
            TProfileChar c = *itH;
            appendValue(buffer, c);
            if (isGap(itV))
                back(buffer).count[valueSize<TChar>()] += 1;
            else
                back(buffer).count[ordValue(TChar(*itV))] += 1;
        }
    }

    //std::cout << "--- AFTER ----\n";
    //std::cout << "HORIZONTAL\n";
    //for (typename Iterator<TProfileString>::Type it = begin(buffer, Standard()); it != end(buffer, Standard()); ++it)
    //{
    //    std::cout << "(" << it->count[0] << ", " << it->count[1] << ", " << it->count[2] << ", "
    //              << it->count[3] << ", " << it->count[4] << ")\n";
    //}

    swap(buffer, profile);
}

}  // namespace seqan

#endif  // #ifndef SEQAN_INCLUDE_SEQAN_ALIGN_PROFILE_ADD_TO_PROFILE_H_
