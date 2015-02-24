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
// Author: Hannes Hauswedell <hannes.hauswedell@fu-berlin.de>
// ==========================================================================
// Reduced Versions of the 24-letter amino acid alphabet
// ==========================================================================

#ifndef SEQAN_REDUCED_AMINOACID_BASE_LATE_H_
#define SEQAN_REDUCED_AMINOACID_BASE_LATE_H_

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

// -----------------------------------------------------------------------
// Function unknownValueImpl()
// -----------------------------------------------------------------------

template <typename TRedSpec>
inline SimpleType<unsigned char, ReducedAminoAcid_<TRedSpec> >
unknownValueImpl(SimpleType<unsigned char, ReducedAminoAcid_<TRedSpec> > *)
{
    static const SimpleType<unsigned char, ReducedAminoAcid_<TRedSpec> >
      _result = SimpleType<unsigned char, ReducedAminoAcid_<TRedSpec> >('X');
    return _result;
}

// -----------------------------------------------------------------------
// Function assign()
// -----------------------------------------------------------------------

template <typename TRedSpec>
inline void assign(char & c_target, SimpleType<unsigned char,
                   ReducedAminoAcid_<TRedSpec> > const & source)
{
    c_target = TranslateTableRedAAToChar_<TRedSpec>::VALUE[source.value];
}

template <typename TRedSpec>
inline void assign(SimpleType<unsigned char,
                   ReducedAminoAcid_<TRedSpec> > & target,
                   __uint8 c_source)
{
    target.value = TranslateTableByteToRedAA_<TRedSpec>::VALUE[c_source];
}

template <typename TRedSpec>
inline void assign(SimpleType<unsigned char,
                   ReducedAminoAcid_<TRedSpec> > & target,
                   char c_source)
{
    target.value = TranslateTableCharToRedAA_<TRedSpec>::VALUE[(unsigned char) c_source];
}

template <typename TRedSpec>
inline void assign(SimpleType<unsigned char,
                   ReducedAminoAcid_<TRedSpec> > & target,
                   AminoAcid c_source)
{
    target.value = TranslateTableAAToRedAA_<TRedSpec>::VALUE[c_source.value];
}

}
#endif // def SEQAN_REDUCED_AMINOACID_BASE_LATE_H_
