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
// Reduced Versions of the 22-letter amino acid alphabet
// ==========================================================================

#ifndef SEQAN_REDUCED_AMINOACID_CLUSTER_RED_TABLES_22_to_n_B62_H_
#define SEQAN_REDUCED_AMINOACID_CLUSTER_RED_TABLES_22_to_n_B62_H_

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

// ========================== DECLARATIONS FIRST ===========================

// ---------------------------------- N = 21 ------------------------------
template <typename TSpec>
struct TranslateTableRedAAToChar_<ClusterReduction<21, 22, Blosum62>, TSpec>
{
    static const char VALUE[ValueSize<SimpleType<unsigned char, ReducedAminoAcid_<ClusterReduction<21, 22, Blosum62> > > >::VALUE];
};

template <typename TSpec>
struct TranslateTableCharToRedAA_<ClusterReduction<21, 22, Blosum62>, TSpec>
{
    static const char VALUE[256];
};

template <typename TSpec>
struct TranslateTableAAToRedAA_<ClusterReduction<21, 22, Blosum62>, TSpec>
{
    static const char VALUE[24];
};

template <typename TSpec>
struct TranslateTableByteToRedAA_<ClusterReduction<21, 22, Blosum62>, TSpec>
{
    static const char VALUE[256];
};

// ---------------------------------- N = 20 ------------------------------
template <typename TSpec>
struct TranslateTableRedAAToChar_<ClusterReduction<20, 22, Blosum62>, TSpec>
{
    static const char VALUE[ValueSize<SimpleType<unsigned char, ReducedAminoAcid_<ClusterReduction<20, 22, Blosum62> > > >::VALUE];
};

template <typename TSpec>
struct TranslateTableCharToRedAA_<ClusterReduction<20, 22, Blosum62>, TSpec>
{
    static const char VALUE[256];
};

template <typename TSpec>
struct TranslateTableAAToRedAA_<ClusterReduction<20, 22, Blosum62>, TSpec>
{
    static const char VALUE[24];
};

template <typename TSpec>
struct TranslateTableByteToRedAA_<ClusterReduction<20, 22, Blosum62>, TSpec>
{
    static const char VALUE[256];
};

// ---------------------------------- N = 19 ------------------------------
template <typename TSpec>
struct TranslateTableRedAAToChar_<ClusterReduction<19, 22, Blosum62>, TSpec>
{
    static const char VALUE[ValueSize<SimpleType<unsigned char, ReducedAminoAcid_<ClusterReduction<19, 22, Blosum62> > > >::VALUE];
};

template <typename TSpec>
struct TranslateTableCharToRedAA_<ClusterReduction<19, 22, Blosum62>, TSpec>
{
    static const char VALUE[256];
};

template <typename TSpec>
struct TranslateTableAAToRedAA_<ClusterReduction<19, 22, Blosum62>, TSpec>
{
    static const char VALUE[24];
};

template <typename TSpec>
struct TranslateTableByteToRedAA_<ClusterReduction<19, 22, Blosum62>, TSpec>
{
    static const char VALUE[256];
};

// ---------------------------------- N = 18 ------------------------------
template <typename TSpec>
struct TranslateTableRedAAToChar_<ClusterReduction<18, 22, Blosum62>, TSpec>
{
    static const char VALUE[ValueSize<SimpleType<unsigned char, ReducedAminoAcid_<ClusterReduction<18, 22, Blosum62> > > >::VALUE];
};

template <typename TSpec>
struct TranslateTableCharToRedAA_<ClusterReduction<18, 22, Blosum62>, TSpec>
{
    static const char VALUE[256];
};

template <typename TSpec>
struct TranslateTableAAToRedAA_<ClusterReduction<18, 22, Blosum62>, TSpec>
{
    static const char VALUE[24];
};

template <typename TSpec>
struct TranslateTableByteToRedAA_<ClusterReduction<18, 22, Blosum62>, TSpec>
{
    static const char VALUE[256];
};

// ---------------------------------- N = 17 ------------------------------
template <typename TSpec>
struct TranslateTableRedAAToChar_<ClusterReduction<17, 22, Blosum62>, TSpec>
{
    static const char VALUE[ValueSize<SimpleType<unsigned char, ReducedAminoAcid_<ClusterReduction<17, 22, Blosum62> > > >::VALUE];
};

template <typename TSpec>
struct TranslateTableCharToRedAA_<ClusterReduction<17, 22, Blosum62>, TSpec>
{
    static const char VALUE[256];
};

template <typename TSpec>
struct TranslateTableAAToRedAA_<ClusterReduction<17, 22, Blosum62>, TSpec>
{
    static const char VALUE[24];
};

template <typename TSpec>
struct TranslateTableByteToRedAA_<ClusterReduction<17, 22, Blosum62>, TSpec>
{
    static const char VALUE[256];
};

// ---------------------------------- N = 16 ------------------------------
template <typename TSpec>
struct TranslateTableRedAAToChar_<ClusterReduction<16, 22, Blosum62>, TSpec>
{
    static const char VALUE[ValueSize<SimpleType<unsigned char, ReducedAminoAcid_<ClusterReduction<16, 22, Blosum62> > > >::VALUE];
};

template <typename TSpec>
struct TranslateTableCharToRedAA_<ClusterReduction<16, 22, Blosum62>, TSpec>
{
    static const char VALUE[256];
};

template <typename TSpec>
struct TranslateTableAAToRedAA_<ClusterReduction<16, 22, Blosum62>, TSpec>
{
    static const char VALUE[24];
};

template <typename TSpec>
struct TranslateTableByteToRedAA_<ClusterReduction<16, 22, Blosum62>, TSpec>
{
    static const char VALUE[256];
};

// ---------------------------------- N = 15 ------------------------------
template <typename TSpec>
struct TranslateTableRedAAToChar_<ClusterReduction<15, 22, Blosum62>, TSpec>
{
    static const char VALUE[ValueSize<SimpleType<unsigned char, ReducedAminoAcid_<ClusterReduction<15, 22, Blosum62> > > >::VALUE];
};

template <typename TSpec>
struct TranslateTableCharToRedAA_<ClusterReduction<15, 22, Blosum62>, TSpec>
{
    static const char VALUE[256];
};

template <typename TSpec>
struct TranslateTableAAToRedAA_<ClusterReduction<15, 22, Blosum62>, TSpec>
{
    static const char VALUE[24];
};

template <typename TSpec>
struct TranslateTableByteToRedAA_<ClusterReduction<15, 22, Blosum62>, TSpec>
{
    static const char VALUE[256];
};

// ---------------------------------- N = 14 ------------------------------
template <typename TSpec>
struct TranslateTableRedAAToChar_<ClusterReduction<14, 22, Blosum62>, TSpec>
{
    static const char VALUE[ValueSize<SimpleType<unsigned char, ReducedAminoAcid_<ClusterReduction<14, 22, Blosum62> > > >::VALUE];
};

template <typename TSpec>
struct TranslateTableCharToRedAA_<ClusterReduction<14, 22, Blosum62>, TSpec>
{
    static const char VALUE[256];
};

template <typename TSpec>
struct TranslateTableAAToRedAA_<ClusterReduction<14, 22, Blosum62>, TSpec>
{
    static const char VALUE[24];
};

template <typename TSpec>
struct TranslateTableByteToRedAA_<ClusterReduction<14, 22, Blosum62>, TSpec>
{
    static const char VALUE[256];
};

// ---------------------------------- N = 13 ------------------------------
template <typename TSpec>
struct TranslateTableRedAAToChar_<ClusterReduction<13, 22, Blosum62>, TSpec>
{
    static const char VALUE[ValueSize<SimpleType<unsigned char, ReducedAminoAcid_<ClusterReduction<13, 22, Blosum62> > > >::VALUE];
};

template <typename TSpec>
struct TranslateTableCharToRedAA_<ClusterReduction<13, 22, Blosum62>, TSpec>
{
    static const char VALUE[256];
};

template <typename TSpec>
struct TranslateTableAAToRedAA_<ClusterReduction<13, 22, Blosum62>, TSpec>
{
    static const char VALUE[24];
};

template <typename TSpec>
struct TranslateTableByteToRedAA_<ClusterReduction<13, 22, Blosum62>, TSpec>
{
    static const char VALUE[256];
};

// ---------------------------------- N = 12 ------------------------------
template <typename TSpec>
struct TranslateTableRedAAToChar_<ClusterReduction<12, 22, Blosum62>, TSpec>
{
    static const char VALUE[ValueSize<SimpleType<unsigned char, ReducedAminoAcid_<ClusterReduction<12, 22, Blosum62> > > >::VALUE];
};

template <typename TSpec>
struct TranslateTableCharToRedAA_<ClusterReduction<12, 22, Blosum62>, TSpec>
{
    static const char VALUE[256];
};

template <typename TSpec>
struct TranslateTableAAToRedAA_<ClusterReduction<12, 22, Blosum62>, TSpec>
{
    static const char VALUE[24];
};

template <typename TSpec>
struct TranslateTableByteToRedAA_<ClusterReduction<12, 22, Blosum62>, TSpec>
{
    static const char VALUE[256];
};

// ---------------------------------- N = 11 ------------------------------
template <typename TSpec>
struct TranslateTableRedAAToChar_<ClusterReduction<11, 22, Blosum62>, TSpec>
{
    static const char VALUE[ValueSize<SimpleType<unsigned char, ReducedAminoAcid_<ClusterReduction<11, 22, Blosum62> > > >::VALUE];
};

template <typename TSpec>
struct TranslateTableCharToRedAA_<ClusterReduction<11, 22, Blosum62>, TSpec>
{
    static const char VALUE[256];
};

template <typename TSpec>
struct TranslateTableAAToRedAA_<ClusterReduction<11, 22, Blosum62>, TSpec>
{
    static const char VALUE[24];
};

template <typename TSpec>
struct TranslateTableByteToRedAA_<ClusterReduction<11, 22, Blosum62>, TSpec>
{
    static const char VALUE[256];
};

// ---------------------------------- N = 10 ------------------------------
template <typename TSpec>
struct TranslateTableRedAAToChar_<ClusterReduction<10, 22, Blosum62>, TSpec>
{
    static const char VALUE[ValueSize<SimpleType<unsigned char, ReducedAminoAcid_<ClusterReduction<10, 22, Blosum62> > > >::VALUE];
};

template <typename TSpec>
struct TranslateTableCharToRedAA_<ClusterReduction<10, 22, Blosum62>, TSpec>
{
    static const char VALUE[256];
};

template <typename TSpec>
struct TranslateTableAAToRedAA_<ClusterReduction<10, 22, Blosum62>, TSpec>
{
    static const char VALUE[24];
};

template <typename TSpec>
struct TranslateTableByteToRedAA_<ClusterReduction<10, 22, Blosum62>, TSpec>
{
    static const char VALUE[256];
};

// ---------------------------------- N = 9 ------------------------------
template <typename TSpec>
struct TranslateTableRedAAToChar_<ClusterReduction<9, 22, Blosum62>, TSpec>
{
    static const char VALUE[ValueSize<SimpleType<unsigned char, ReducedAminoAcid_<ClusterReduction<9, 22, Blosum62> > > >::VALUE];
};

template <typename TSpec>
struct TranslateTableCharToRedAA_<ClusterReduction<9, 22, Blosum62>, TSpec>
{
    static const char VALUE[256];
};

template <typename TSpec>
struct TranslateTableAAToRedAA_<ClusterReduction<9, 22, Blosum62>, TSpec>
{
    static const char VALUE[24];
};

template <typename TSpec>
struct TranslateTableByteToRedAA_<ClusterReduction<9, 22, Blosum62>, TSpec>
{
    static const char VALUE[256];
};

// ---------------------------------- N = 8 ------------------------------
template <typename TSpec>
struct TranslateTableRedAAToChar_<ClusterReduction<8, 22, Blosum62>, TSpec>
{
    static const char VALUE[ValueSize<SimpleType<unsigned char, ReducedAminoAcid_<ClusterReduction<8, 22, Blosum62> > > >::VALUE];
};

template <typename TSpec>
struct TranslateTableCharToRedAA_<ClusterReduction<8, 22, Blosum62>, TSpec>
{
    static const char VALUE[256];
};

template <typename TSpec>
struct TranslateTableAAToRedAA_<ClusterReduction<8, 22, Blosum62>, TSpec>
{
    static const char VALUE[24];
};

template <typename TSpec>
struct TranslateTableByteToRedAA_<ClusterReduction<8, 22, Blosum62>, TSpec>
{
    static const char VALUE[256];
};

// ---------------------------------- N = 7 ------------------------------
template <typename TSpec>
struct TranslateTableRedAAToChar_<ClusterReduction<7, 22, Blosum62>, TSpec>
{
    static const char VALUE[ValueSize<SimpleType<unsigned char, ReducedAminoAcid_<ClusterReduction<7, 22, Blosum62> > > >::VALUE];
};

template <typename TSpec>
struct TranslateTableCharToRedAA_<ClusterReduction<7, 22, Blosum62>, TSpec>
{
    static const char VALUE[256];
};

template <typename TSpec>
struct TranslateTableAAToRedAA_<ClusterReduction<7, 22, Blosum62>, TSpec>
{
    static const char VALUE[24];
};

template <typename TSpec>
struct TranslateTableByteToRedAA_<ClusterReduction<7, 22, Blosum62>, TSpec>
{
    static const char VALUE[256];
};

// ---------------------------------- N = 6 ------------------------------
template <typename TSpec>
struct TranslateTableRedAAToChar_<ClusterReduction<6, 22, Blosum62>, TSpec>
{
    static const char VALUE[ValueSize<SimpleType<unsigned char, ReducedAminoAcid_<ClusterReduction<6, 22, Blosum62> > > >::VALUE];
};

template <typename TSpec>
struct TranslateTableCharToRedAA_<ClusterReduction<6, 22, Blosum62>, TSpec>
{
    static const char VALUE[256];
};

template <typename TSpec>
struct TranslateTableAAToRedAA_<ClusterReduction<6, 22, Blosum62>, TSpec>
{
    static const char VALUE[24];
};

template <typename TSpec>
struct TranslateTableByteToRedAA_<ClusterReduction<6, 22, Blosum62>, TSpec>
{
    static const char VALUE[256];
};

// ---------------------------------- N = 5 ------------------------------
template <typename TSpec>
struct TranslateTableRedAAToChar_<ClusterReduction<5, 22, Blosum62>, TSpec>
{
    static const char VALUE[ValueSize<SimpleType<unsigned char, ReducedAminoAcid_<ClusterReduction<5, 22, Blosum62> > > >::VALUE];
};

template <typename TSpec>
struct TranslateTableCharToRedAA_<ClusterReduction<5, 22, Blosum62>, TSpec>
{
    static const char VALUE[256];
};

template <typename TSpec>
struct TranslateTableAAToRedAA_<ClusterReduction<5, 22, Blosum62>, TSpec>
{
    static const char VALUE[24];
};

template <typename TSpec>
struct TranslateTableByteToRedAA_<ClusterReduction<5, 22, Blosum62>, TSpec>
{
    static const char VALUE[256];
};

// ---------------------------------- N = 4 ------------------------------
template <typename TSpec>
struct TranslateTableRedAAToChar_<ClusterReduction<4, 22, Blosum62>, TSpec>
{
    static const char VALUE[ValueSize<SimpleType<unsigned char, ReducedAminoAcid_<ClusterReduction<4, 22, Blosum62> > > >::VALUE];
};

template <typename TSpec>
struct TranslateTableCharToRedAA_<ClusterReduction<4, 22, Blosum62>, TSpec>
{
    static const char VALUE[256];
};

template <typename TSpec>
struct TranslateTableAAToRedAA_<ClusterReduction<4, 22, Blosum62>, TSpec>
{
    static const char VALUE[24];
};

template <typename TSpec>
struct TranslateTableByteToRedAA_<ClusterReduction<4, 22, Blosum62>, TSpec>
{
    static const char VALUE[256];
};

// ---------------------------------- N = 3 ------------------------------
template <typename TSpec>
struct TranslateTableRedAAToChar_<ClusterReduction<3, 22, Blosum62>, TSpec>
{
    static const char VALUE[ValueSize<SimpleType<unsigned char, ReducedAminoAcid_<ClusterReduction<3, 22, Blosum62> > > >::VALUE];
};

template <typename TSpec>
struct TranslateTableCharToRedAA_<ClusterReduction<3, 22, Blosum62>, TSpec>
{
    static const char VALUE[256];
};

template <typename TSpec>
struct TranslateTableAAToRedAA_<ClusterReduction<3, 22, Blosum62>, TSpec>
{
    static const char VALUE[24];
};

template <typename TSpec>
struct TranslateTableByteToRedAA_<ClusterReduction<3, 22, Blosum62>, TSpec>
{
    static const char VALUE[256];
};

// ---------------------------------- N = 2 ------------------------------
template <typename TSpec>
struct TranslateTableRedAAToChar_<ClusterReduction<2, 22, Blosum62>, TSpec>
{
    static const char VALUE[ValueSize<SimpleType<unsigned char, ReducedAminoAcid_<ClusterReduction<2, 22, Blosum62> > > >::VALUE];
};

template <typename TSpec>
struct TranslateTableCharToRedAA_<ClusterReduction<2, 22, Blosum62>, TSpec>
{
    static const char VALUE[256];
};

template <typename TSpec>
struct TranslateTableAAToRedAA_<ClusterReduction<2, 22, Blosum62>, TSpec>
{
    static const char VALUE[24];
};

template <typename TSpec>
struct TranslateTableByteToRedAA_<ClusterReduction<2, 22, Blosum62>, TSpec>
{
    static const char VALUE[256];
};

// ========================== DEFINITIONS SECOND ===========================

// ---------------------------------- N = 21 ------------------------------

template <typename TVoidSpec>
char const TranslateTableRedAAToChar_<
                ClusterReduction<21, 22, Blosum62>, TVoidSpec>::VALUE[21] =
{
    'A', // A
    'R', // R
    'N', // N
    'D', // D
    'C', // C
    'Q', // Q
    'E', // E Z
    'G', // G
    'H', // H
    'I', // I
    'L', // L
    'K', // K
    'M', // M
    'F', // F
    'P', // P
    'S', // S
    'T', // T
    'W', // W
    'Y', // Y
    'V', // V
    'B'  // B
};

template <typename TVoidSpec>
char const TranslateTableCharToRedAA_<
                ClusterReduction<21, 22, Blosum62>, TVoidSpec>::VALUE[256] =
{
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0, 20,  4,  3,  6, 13,  7,  8,  9,  0,
    11, 10, 12,  2,  0, 14,  5,  1, 15, 16,  0, 19, 17,  0, 18,
     6,  0,  0,  0,  0,  0,  0,  0, 20,  4,  3,  6, 13,  7,  8,
     9,  0, 11, 10, 12,  2,  0, 14,  5,  1, 15, 16,  0, 19, 17,
     0, 18,  6,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0
};

template <typename TVoidSpec>
char const TranslateTableAAToRedAA_<
                ClusterReduction<21, 22, Blosum62>, TVoidSpec>::VALUE[24] =
{
     0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11,
    12, 13, 14, 15, 16, 17, 18, 19, 20,  6,  0,  0
};

template <typename TVoidSpec>
char const TranslateTableByteToRedAA_<
                ClusterReduction<21, 22, Blosum62>, TVoidSpec>::VALUE[256] =
{
     0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14,
    15, 16, 17, 18, 19, 20,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0
};

// ---------------------------------- N = 20 ------------------------------

template <typename TVoidSpec>
char const TranslateTableRedAAToChar_<
                ClusterReduction<20, 22, Blosum62>, TVoidSpec>::VALUE[20] =
{
    'A', // A
    'R', // R
    'N', // N
    'D', // D B
    'C', // C
    'Q', // Q
    'E', // E Z
    'G', // G
    'H', // H
    'I', // I
    'L', // L
    'K', // K
    'M', // M
    'F', // F
    'P', // P
    'S', // S
    'T', // T
    'W', // W
    'Y', // Y
    'V'  // V
};

template <typename TVoidSpec>
char const TranslateTableCharToRedAA_<
                ClusterReduction<20, 22, Blosum62>, TVoidSpec>::VALUE[256] =
{
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  3,  4,  3,  6, 13,  7,  8,  9,  0,
    11, 10, 12,  2,  0, 14,  5,  1, 15, 16,  0, 19, 17,  0, 18,
     6,  0,  0,  0,  0,  0,  0,  0,  3,  4,  3,  6, 13,  7,  8,
     9,  0, 11, 10, 12,  2,  0, 14,  5,  1, 15, 16,  0, 19, 17,
     0, 18,  6,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0
};

template <typename TVoidSpec>
char const TranslateTableAAToRedAA_<
                ClusterReduction<20, 22, Blosum62>, TVoidSpec>::VALUE[24] =
{
     0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11,
    12, 13, 14, 15, 16, 17, 18, 19,  3,  6,  0,  0
};

template <typename TVoidSpec>
char const TranslateTableByteToRedAA_<
                ClusterReduction<20, 22, Blosum62>, TVoidSpec>::VALUE[256] =
{
     0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14,
    15, 16, 17, 18, 19,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0
};

// ---------------------------------- N = 19 ------------------------------

template <typename TVoidSpec>
char const TranslateTableRedAAToChar_<
                ClusterReduction<19, 22, Blosum62>, TVoidSpec>::VALUE[19] =
{
    'A', // A
    'R', // R
    'N', // N
    'D', // D B
    'C', // C
    'Q', // Q
    'E', // E Z
    'G', // G
    'H', // H
    'I', // I V
    'L', // L
    'K', // K
    'M', // M
    'F', // F
    'P', // P
    'S', // S
    'T', // T
    'W', // W
    'Y'  // Y
};

template <typename TVoidSpec>
char const TranslateTableCharToRedAA_<
                ClusterReduction<19, 22, Blosum62>, TVoidSpec>::VALUE[256] =
{
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  3,  4,  3,  6, 13,  7,  8,  9,  0,
    11, 10, 12,  2,  0, 14,  5,  1, 15, 16,  0,  9, 17,  0, 18,
     6,  0,  0,  0,  0,  0,  0,  0,  3,  4,  3,  6, 13,  7,  8,
     9,  0, 11, 10, 12,  2,  0, 14,  5,  1, 15, 16,  0,  9, 17,
     0, 18,  6,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0
};

template <typename TVoidSpec>
char const TranslateTableAAToRedAA_<
                ClusterReduction<19, 22, Blosum62>, TVoidSpec>::VALUE[24] =
{
     0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11,
    12, 13, 14, 15, 16, 17, 18,  9,  3,  6,  0,  0
};

template <typename TVoidSpec>
char const TranslateTableByteToRedAA_<
                ClusterReduction<19, 22, Blosum62>, TVoidSpec>::VALUE[256] =
{
     0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14,
    15, 16, 17, 18,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0
};

// ---------------------------------- N = 18 ------------------------------

template <typename TVoidSpec>
char const TranslateTableRedAAToChar_<
                ClusterReduction<18, 22, Blosum62>, TVoidSpec>::VALUE[18] =
{
    'A', // A
    'R', // R
    'N', // N
    'D', // D B
    'C', // C
    'Q', // Q E Z
    'G', // G
    'H', // H
    'I', // I V
    'L', // L
    'K', // K
    'M', // M
    'F', // F
    'P', // P
    'S', // S
    'T', // T
    'W', // W
    'Y'  // Y
};

template <typename TVoidSpec>
char const TranslateTableCharToRedAA_<
                ClusterReduction<18, 22, Blosum62>, TVoidSpec>::VALUE[256] =
{
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  3,  4,  3,  5, 12,  6,  7,  8,  0,
    10,  9, 11,  2,  0, 13,  5,  1, 14, 15,  0,  8, 16,  0, 17,
     5,  0,  0,  0,  0,  0,  0,  0,  3,  4,  3,  5, 12,  6,  7,
     8,  0, 10,  9, 11,  2,  0, 13,  5,  1, 14, 15,  0,  8, 16,
     0, 17,  5,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0
};

template <typename TVoidSpec>
char const TranslateTableAAToRedAA_<
                ClusterReduction<18, 22, Blosum62>, TVoidSpec>::VALUE[24] =
{
     0,  1,  2,  3,  4,  5,  5,  6,  7,  8,  9, 10,
    11, 12, 13, 14, 15, 16, 17,  8,  3,  5,  0,  0
};

template <typename TVoidSpec>
char const TranslateTableByteToRedAA_<
                ClusterReduction<18, 22, Blosum62>, TVoidSpec>::VALUE[256] =
{
     0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14,
    15, 16, 17,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0
};

// ---------------------------------- N = 17 ------------------------------

template <typename TVoidSpec>
char const TranslateTableRedAAToChar_<
                ClusterReduction<17, 22, Blosum62>, TVoidSpec>::VALUE[17] =
{
    'A', // A
    'R', // R
    'N', // N
    'D', // D B
    'C', // C
    'Q', // Q E Z
    'G', // G
    'H', // H
    'I', // I V
    'L', // L M
    'K', // K
    'F', // F
    'P', // P
    'S', // S
    'T', // T
    'W', // W
    'Y'  // Y
};

template <typename TVoidSpec>
char const TranslateTableCharToRedAA_<
                ClusterReduction<17, 22, Blosum62>, TVoidSpec>::VALUE[256] =
{
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  3,  4,  3,  5, 11,  6,  7,  8,  0,
    10,  9,  9,  2,  0, 12,  5,  1, 13, 14,  0,  8, 15,  0, 16,
     5,  0,  0,  0,  0,  0,  0,  0,  3,  4,  3,  5, 11,  6,  7,
     8,  0, 10,  9,  9,  2,  0, 12,  5,  1, 13, 14,  0,  8, 15,
     0, 16,  5,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0
};

template <typename TVoidSpec>
char const TranslateTableAAToRedAA_<
                ClusterReduction<17, 22, Blosum62>, TVoidSpec>::VALUE[24] =
{
     0,  1,  2,  3,  4,  5,  5,  6,  7,  8,  9, 10,
     9, 11, 12, 13, 14, 15, 16,  8,  3,  5,  0,  0
};

template <typename TVoidSpec>
char const TranslateTableByteToRedAA_<
                ClusterReduction<17, 22, Blosum62>, TVoidSpec>::VALUE[256] =
{
     0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14,
    15, 16,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0
};

// ---------------------------------- N = 16 ------------------------------

template <typename TVoidSpec>
char const TranslateTableRedAAToChar_<
                ClusterReduction<16, 22, Blosum62>, TVoidSpec>::VALUE[16] =
{
    'A', // A
    'R', // R K
    'N', // N
    'D', // D B
    'C', // C
    'Q', // Q E Z
    'G', // G
    'H', // H
    'I', // I V
    'L', // L M
    'F', // F
    'P', // P
    'S', // S
    'T', // T
    'W', // W
    'Y'  // Y
};

template <typename TVoidSpec>
char const TranslateTableCharToRedAA_<
                ClusterReduction<16, 22, Blosum62>, TVoidSpec>::VALUE[256] =
{
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  3,  4,  3,  5, 10,  6,  7,  8,  0,
     1,  9,  9,  2,  0, 11,  5,  1, 12, 13,  0,  8, 14,  0, 15,
     5,  0,  0,  0,  0,  0,  0,  0,  3,  4,  3,  5, 10,  6,  7,
     8,  0,  1,  9,  9,  2,  0, 11,  5,  1, 12, 13,  0,  8, 14,
     0, 15,  5,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0
};

template <typename TVoidSpec>
char const TranslateTableAAToRedAA_<
                ClusterReduction<16, 22, Blosum62>, TVoidSpec>::VALUE[24] =
{
     0,  1,  2,  3,  4,  5,  5,  6,  7,  8,  9,  1,
     9, 10, 11, 12, 13, 14, 15,  8,  3,  5,  0,  0
};

template <typename TVoidSpec>
char const TranslateTableByteToRedAA_<
                ClusterReduction<16, 22, Blosum62>, TVoidSpec>::VALUE[256] =
{
     0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14,
    15,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0
};

// ---------------------------------- N = 15 ------------------------------

template <typename TVoidSpec>
char const TranslateTableRedAAToChar_<
                ClusterReduction<15, 22, Blosum62>, TVoidSpec>::VALUE[15] =
{
    'A', // A
    'R', // R K
    'N', // N
    'D', // D B
    'C', // C
    'Q', // Q E Z
    'G', // G
    'H', // H
    'I', // I V
    'L', // L M
    'F', // F Y
    'P', // P
    'S', // S
    'T', // T
    'W'  // W
};

template <typename TVoidSpec>
char const TranslateTableCharToRedAA_<
                ClusterReduction<15, 22, Blosum62>, TVoidSpec>::VALUE[256] =
{
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  3,  4,  3,  5, 10,  6,  7,  8,  0,
     1,  9,  9,  2,  0, 11,  5,  1, 12, 13,  0,  8, 14,  0, 10,
     5,  0,  0,  0,  0,  0,  0,  0,  3,  4,  3,  5, 10,  6,  7,
     8,  0,  1,  9,  9,  2,  0, 11,  5,  1, 12, 13,  0,  8, 14,
     0, 10,  5,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0
};

template <typename TVoidSpec>
char const TranslateTableAAToRedAA_<
                ClusterReduction<15, 22, Blosum62>, TVoidSpec>::VALUE[24] =
{
     0,  1,  2,  3,  4,  5,  5,  6,  7,  8,  9,  1,
     9, 10, 11, 12, 13, 14, 10,  8,  3,  5,  0,  0
};

template <typename TVoidSpec>
char const TranslateTableByteToRedAA_<
                ClusterReduction<15, 22, Blosum62>, TVoidSpec>::VALUE[256] =
{
     0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13, 14,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0
};

// ---------------------------------- N = 14 ------------------------------

template <typename TVoidSpec>
char const TranslateTableRedAAToChar_<
                ClusterReduction<14, 22, Blosum62>, TVoidSpec>::VALUE[14] =
{
    'A', // A
    'R', // R K
    'N', // N D B
    'C', // C
    'Q', // Q E Z
    'G', // G
    'H', // H
    'I', // I V
    'L', // L M
    'F', // F Y
    'P', // P
    'S', // S
    'T', // T
    'W'  // W
};

template <typename TVoidSpec>
char const TranslateTableCharToRedAA_<
                ClusterReduction<14, 22, Blosum62>, TVoidSpec>::VALUE[256] =
{
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  2,  3,  2,  4,  9,  5,  6,  7,  0,
     1,  8,  8,  2,  0, 10,  4,  1, 11, 12,  0,  7, 13,  0,  9,
     4,  0,  0,  0,  0,  0,  0,  0,  2,  3,  2,  4,  9,  5,  6,
     7,  0,  1,  8,  8,  2,  0, 10,  4,  1, 11, 12,  0,  7, 13,
     0,  9,  4,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0
};

template <typename TVoidSpec>
char const TranslateTableAAToRedAA_<
                ClusterReduction<14, 22, Blosum62>, TVoidSpec>::VALUE[24] =
{
     0,  1,  2,  2,  3,  4,  4,  5,  6,  7,  8,  1,
     8,  9, 10, 11, 12, 13,  9,  7,  2,  4,  0,  0
};

template <typename TVoidSpec>
char const TranslateTableByteToRedAA_<
                ClusterReduction<14, 22, Blosum62>, TVoidSpec>::VALUE[256] =
{
     0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12, 13,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0
};

// ---------------------------------- N = 13 ------------------------------

template <typename TVoidSpec>
char const TranslateTableRedAAToChar_<
                ClusterReduction<13, 22, Blosum62>, TVoidSpec>::VALUE[13] =
{
    'A', // A
    'R', // R K
    'N', // N D B
    'C', // C
    'Q', // Q E Z
    'G', // G
    'H', // H
    'I', // I L M V
    'F', // F Y
    'P', // P
    'S', // S
    'T', // T
    'W'  // W
};

template <typename TVoidSpec>
char const TranslateTableCharToRedAA_<
                ClusterReduction<13, 22, Blosum62>, TVoidSpec>::VALUE[256] =
{
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  2,  3,  2,  4,  8,  5,  6,  7,  0,
     1,  7,  7,  2,  0,  9,  4,  1, 10, 11,  0,  7, 12,  0,  8,
     4,  0,  0,  0,  0,  0,  0,  0,  2,  3,  2,  4,  8,  5,  6,
     7,  0,  1,  7,  7,  2,  0,  9,  4,  1, 10, 11,  0,  7, 12,
     0,  8,  4,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0
};

template <typename TVoidSpec>
char const TranslateTableAAToRedAA_<
                ClusterReduction<13, 22, Blosum62>, TVoidSpec>::VALUE[24] =
{
     0,  1,  2,  2,  3,  4,  4,  5,  6,  7,  7,  1,
     7,  8,  9, 10, 11, 12,  8,  7,  2,  4,  0,  0
};

template <typename TVoidSpec>
char const TranslateTableByteToRedAA_<
                ClusterReduction<13, 22, Blosum62>, TVoidSpec>::VALUE[256] =
{
     0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11, 12,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0
};

// ---------------------------------- N = 12 ------------------------------

template <typename TVoidSpec>
char const TranslateTableRedAAToChar_<
                ClusterReduction<12, 22, Blosum62>, TVoidSpec>::VALUE[12] =
{
    'A', // A
    'R', // R Q E K Z
    'N', // N D B
    'C', // C
    'G', // G
    'H', // H
    'I', // I L M V
    'F', // F Y
    'P', // P
    'S', // S
    'T', // T
    'W'  // W
};

template <typename TVoidSpec>
char const TranslateTableCharToRedAA_<
                ClusterReduction<12, 22, Blosum62>, TVoidSpec>::VALUE[256] =
{
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  2,  3,  2,  1,  7,  4,  5,  6,  0,
     1,  6,  6,  2,  0,  8,  1,  1,  9, 10,  0,  6, 11,  0,  7,
     1,  0,  0,  0,  0,  0,  0,  0,  2,  3,  2,  1,  7,  4,  5,
     6,  0,  1,  6,  6,  2,  0,  8,  1,  1,  9, 10,  0,  6, 11,
     0,  7,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0
};

template <typename TVoidSpec>
char const TranslateTableAAToRedAA_<
                ClusterReduction<12, 22, Blosum62>, TVoidSpec>::VALUE[24] =
{
     0,  1,  2,  2,  3,  1,  1,  4,  5,  6,  6,  1,
     6,  7,  8,  9, 10, 11,  7,  6,  2,  1,  0,  0
};

template <typename TVoidSpec>
char const TranslateTableByteToRedAA_<
                ClusterReduction<12, 22, Blosum62>, TVoidSpec>::VALUE[256] =
{
     0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10, 11,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0
};

// ---------------------------------- N = 11 ------------------------------

template <typename TVoidSpec>
char const TranslateTableRedAAToChar_<
                ClusterReduction<11, 22, Blosum62>, TVoidSpec>::VALUE[11] =
{
    'A', // A
    'R', // R Q E K Z
    'N', // N D B
    'C', // C
    'G', // G
    'H', // H
    'I', // I L M V
    'F', // F W Y
    'P', // P
    'S', // S
    'T'  // T
};

template <typename TVoidSpec>
char const TranslateTableCharToRedAA_<
                ClusterReduction<11, 22, Blosum62>, TVoidSpec>::VALUE[256] =
{
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  2,  3,  2,  1,  7,  4,  5,  6,  0,
     1,  6,  6,  2,  0,  8,  1,  1,  9, 10,  0,  6,  7,  0,  7,
     1,  0,  0,  0,  0,  0,  0,  0,  2,  3,  2,  1,  7,  4,  5,
     6,  0,  1,  6,  6,  2,  0,  8,  1,  1,  9, 10,  0,  6,  7,
     0,  7,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0
};

template <typename TVoidSpec>
char const TranslateTableAAToRedAA_<
                ClusterReduction<11, 22, Blosum62>, TVoidSpec>::VALUE[24] =
{
     0,  1,  2,  2,  3,  1,  1,  4,  5,  6,  6,  1,
     6,  7,  8,  9, 10,  7,  7,  6,  2,  1,  0,  0
};

template <typename TVoidSpec>
char const TranslateTableByteToRedAA_<
                ClusterReduction<11, 22, Blosum62>, TVoidSpec>::VALUE[256] =
{
     0,  1,  2,  3,  4,  5,  6,  7,  8,  9, 10,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0
};

// ---------------------------------- N = 10 ------------------------------

template <typename TVoidSpec>
char const TranslateTableRedAAToChar_<
                ClusterReduction<10, 22, Blosum62>, TVoidSpec>::VALUE[10] =
{
    'A', // A
    'R', // R Q E K Z
    'N', // N D S B
    'C', // C
    'G', // G
    'H', // H
    'I', // I L M V
    'F', // F W Y
    'P', // P
    'T'  // T
};

template <typename TVoidSpec>
char const TranslateTableCharToRedAA_<
                ClusterReduction<10, 22, Blosum62>, TVoidSpec>::VALUE[256] =
{
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  2,  3,  2,  1,  7,  4,  5,  6,  0,
     1,  6,  6,  2,  0,  8,  1,  1,  2,  9,  0,  6,  7,  0,  7,
     1,  0,  0,  0,  0,  0,  0,  0,  2,  3,  2,  1,  7,  4,  5,
     6,  0,  1,  6,  6,  2,  0,  8,  1,  1,  2,  9,  0,  6,  7,
     0,  7,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0
};

template <typename TVoidSpec>
char const TranslateTableAAToRedAA_<
                ClusterReduction<10, 22, Blosum62>, TVoidSpec>::VALUE[24] =
{
     0,  1,  2,  2,  3,  1,  1,  4,  5,  6,  6,  1,
     6,  7,  8,  2,  9,  7,  7,  6,  2,  1,  0,  0
};

template <typename TVoidSpec>
char const TranslateTableByteToRedAA_<
                ClusterReduction<10, 22, Blosum62>, TVoidSpec>::VALUE[256] =
{
     0,  1,  2,  3,  4,  5,  6,  7,  8,  9,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0
};

// ---------------------------------- N = 9 ------------------------------

template <typename TVoidSpec>
char const TranslateTableRedAAToChar_<
                ClusterReduction<9, 22, Blosum62>, TVoidSpec>::VALUE[9] =
{
    'A', // A
    'R', // R N D Q E K S B Z
    'C', // C
    'G', // G
    'H', // H
    'I', // I L M V
    'F', // F W Y
    'P', // P
    'T'  // T
};

template <typename TVoidSpec>
char const TranslateTableCharToRedAA_<
                ClusterReduction<9, 22, Blosum62>, TVoidSpec>::VALUE[256] =
{
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  1,  2,  1,  1,  6,  3,  4,  5,  0,
     1,  5,  5,  1,  0,  7,  1,  1,  1,  8,  0,  5,  6,  0,  6,
     1,  0,  0,  0,  0,  0,  0,  0,  1,  2,  1,  1,  6,  3,  4,
     5,  0,  1,  5,  5,  1,  0,  7,  1,  1,  1,  8,  0,  5,  6,
     0,  6,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0
};

template <typename TVoidSpec>
char const TranslateTableAAToRedAA_<
                ClusterReduction<9, 22, Blosum62>, TVoidSpec>::VALUE[24] =
{
     0,  1,  1,  1,  2,  1,  1,  3,  4,  5,  5,  1,
     5,  6,  7,  1,  8,  6,  6,  5,  1,  1,  0,  0
};

template <typename TVoidSpec>
char const TranslateTableByteToRedAA_<
                ClusterReduction<9, 22, Blosum62>, TVoidSpec>::VALUE[256] =
{
     0,  1,  2,  3,  4,  5,  6,  7,  8,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0
};

// ---------------------------------- N = 8 ------------------------------

template <typename TVoidSpec>
char const TranslateTableRedAAToChar_<
                ClusterReduction<8, 22, Blosum62>, TVoidSpec>::VALUE[8] =
{
    'A', // A T
    'R', // R N D Q E K S B Z
    'C', // C
    'G', // G
    'H', // H
    'I', // I L M V
    'F', // F W Y
    'P'  // P
};

template <typename TVoidSpec>
char const TranslateTableCharToRedAA_<
                ClusterReduction<8, 22, Blosum62>, TVoidSpec>::VALUE[256] =
{
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  1,  2,  1,  1,  6,  3,  4,  5,  0,
     1,  5,  5,  1,  0,  7,  1,  1,  1,  0,  0,  5,  6,  0,  6,
     1,  0,  0,  0,  0,  0,  0,  0,  1,  2,  1,  1,  6,  3,  4,
     5,  0,  1,  5,  5,  1,  0,  7,  1,  1,  1,  0,  0,  5,  6,
     0,  6,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0
};

template <typename TVoidSpec>
char const TranslateTableAAToRedAA_<
                ClusterReduction<8, 22, Blosum62>, TVoidSpec>::VALUE[24] =
{
     0,  1,  1,  1,  2,  1,  1,  3,  4,  5,  5,  1,
     5,  6,  7,  1,  0,  6,  6,  5,  1,  1,  0,  0
};

template <typename TVoidSpec>
char const TranslateTableByteToRedAA_<
                ClusterReduction<8, 22, Blosum62>, TVoidSpec>::VALUE[256] =
{
     0,  1,  2,  3,  4,  5,  6,  7,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0
};

// ---------------------------------- N = 7 ------------------------------

template <typename TVoidSpec>
char const TranslateTableRedAAToChar_<
                ClusterReduction<7, 22, Blosum62>, TVoidSpec>::VALUE[7] =
{
    'A', // A T
    'R', // R N D Q E H K S B Z
    'C', // C
    'G', // G
    'I', // I L M V
    'F', // F W Y
    'P'  // P
};

template <typename TVoidSpec>
char const TranslateTableCharToRedAA_<
                ClusterReduction<7, 22, Blosum62>, TVoidSpec>::VALUE[256] =
{
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  1,  2,  1,  1,  5,  3,  1,  4,  0,
     1,  4,  4,  1,  0,  6,  1,  1,  1,  0,  0,  4,  5,  0,  5,
     1,  0,  0,  0,  0,  0,  0,  0,  1,  2,  1,  1,  5,  3,  1,
     4,  0,  1,  4,  4,  1,  0,  6,  1,  1,  1,  0,  0,  4,  5,
     0,  5,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0
};

template <typename TVoidSpec>
char const TranslateTableAAToRedAA_<
                ClusterReduction<7, 22, Blosum62>, TVoidSpec>::VALUE[24] =
{
     0,  1,  1,  1,  2,  1,  1,  3,  1,  4,  4,  1,
     4,  5,  6,  1,  0,  5,  5,  4,  1,  1,  0,  0
};

template <typename TVoidSpec>
char const TranslateTableByteToRedAA_<
                ClusterReduction<7, 22, Blosum62>, TVoidSpec>::VALUE[256] =
{
     0,  1,  2,  3,  4,  5,  6,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0
};

// ---------------------------------- N = 6 ------------------------------

template <typename TVoidSpec>
char const TranslateTableRedAAToChar_<
                ClusterReduction<6, 22, Blosum62>, TVoidSpec>::VALUE[6] =
{
    'A', // A T
    'R', // R N D Q E H K S B Z
    'C', // C I L M V
    'G', // G
    'F', // F W Y
    'P'  // P
};

template <typename TVoidSpec>
char const TranslateTableCharToRedAA_<
                ClusterReduction<6, 22, Blosum62>, TVoidSpec>::VALUE[256] =
{
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  1,  2,  1,  1,  4,  3,  1,  2,  0,
     1,  2,  2,  1,  0,  5,  1,  1,  1,  0,  0,  2,  4,  0,  4,
     1,  0,  0,  0,  0,  0,  0,  0,  1,  2,  1,  1,  4,  3,  1,
     2,  0,  1,  2,  2,  1,  0,  5,  1,  1,  1,  0,  0,  2,  4,
     0,  4,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0
};

template <typename TVoidSpec>
char const TranslateTableAAToRedAA_<
                ClusterReduction<6, 22, Blosum62>, TVoidSpec>::VALUE[24] =
{
     0,  1,  1,  1,  2,  1,  1,  3,  1,  2,  2,  1,
     2,  4,  5,  1,  0,  4,  4,  2,  1,  1,  0,  0
};

template <typename TVoidSpec>
char const TranslateTableByteToRedAA_<
                ClusterReduction<6, 22, Blosum62>, TVoidSpec>::VALUE[256] =
{
     0,  1,  2,  3,  4,  5,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0
};

// ---------------------------------- N = 5 ------------------------------

template <typename TVoidSpec>
char const TranslateTableRedAAToChar_<
                ClusterReduction<5, 22, Blosum62>, TVoidSpec>::VALUE[5] =
{
    'A', // A T
    'R', // R N D Q E G H K S B Z
    'C', // C I L M V
    'F', // F W Y
    'P'  // P
};

template <typename TVoidSpec>
char const TranslateTableCharToRedAA_<
                ClusterReduction<5, 22, Blosum62>, TVoidSpec>::VALUE[256] =
{
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  1,  2,  1,  1,  3,  1,  1,  2,  0,
     1,  2,  2,  1,  0,  4,  1,  1,  1,  0,  0,  2,  3,  0,  3,
     1,  0,  0,  0,  0,  0,  0,  0,  1,  2,  1,  1,  3,  1,  1,
     2,  0,  1,  2,  2,  1,  0,  4,  1,  1,  1,  0,  0,  2,  3,
     0,  3,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0
};

template <typename TVoidSpec>
char const TranslateTableAAToRedAA_<
                ClusterReduction<5, 22, Blosum62>, TVoidSpec>::VALUE[24] =
{
     0,  1,  1,  1,  2,  1,  1,  1,  1,  2,  2,  1,
     2,  3,  4,  1,  0,  3,  3,  2,  1,  1,  0,  0
};

template <typename TVoidSpec>
char const TranslateTableByteToRedAA_<
                ClusterReduction<5, 22, Blosum62>, TVoidSpec>::VALUE[256] =
{
     0,  1,  2,  3,  4,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0
};

// ---------------------------------- N = 4 ------------------------------

template <typename TVoidSpec>
char const TranslateTableRedAAToChar_<
                ClusterReduction<4, 22, Blosum62>, TVoidSpec>::VALUE[4] =
{
    'A', // A T
    'R', // R N D Q E G H K P S B Z
    'C', // C I L M V
    'F'  // F W Y
};

template <typename TVoidSpec>
char const TranslateTableCharToRedAA_<
                ClusterReduction<4, 22, Blosum62>, TVoidSpec>::VALUE[256] =
{
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  1,  2,  1,  1,  3,  1,  1,  2,  0,
     1,  2,  2,  1,  0,  1,  1,  1,  1,  0,  0,  2,  3,  0,  3,
     1,  0,  0,  0,  0,  0,  0,  0,  1,  2,  1,  1,  3,  1,  1,
     2,  0,  1,  2,  2,  1,  0,  1,  1,  1,  1,  0,  0,  2,  3,
     0,  3,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0
};

template <typename TVoidSpec>
char const TranslateTableAAToRedAA_<
                ClusterReduction<4, 22, Blosum62>, TVoidSpec>::VALUE[24] =
{
     0,  1,  1,  1,  2,  1,  1,  1,  1,  2,  2,  1,
     2,  3,  1,  1,  0,  3,  3,  2,  1,  1,  0,  0
};

template <typename TVoidSpec>
char const TranslateTableByteToRedAA_<
                ClusterReduction<4, 22, Blosum62>, TVoidSpec>::VALUE[256] =
{
     0,  1,  2,  3,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0
};

// ---------------------------------- N = 3 ------------------------------

template <typename TVoidSpec>
char const TranslateTableRedAAToChar_<
                ClusterReduction<3, 22, Blosum62>, TVoidSpec>::VALUE[3] =
{
    'A', // A T
    'R', // R N D Q E G H K P S B Z
    'C'  // C I L M F W Y V
};

template <typename TVoidSpec>
char const TranslateTableCharToRedAA_<
                ClusterReduction<3, 22, Blosum62>, TVoidSpec>::VALUE[256] =
{
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  1,  2,  1,  1,  2,  1,  1,  2,  0,
     1,  2,  2,  1,  0,  1,  1,  1,  1,  0,  0,  2,  2,  0,  2,
     1,  0,  0,  0,  0,  0,  0,  0,  1,  2,  1,  1,  2,  1,  1,
     2,  0,  1,  2,  2,  1,  0,  1,  1,  1,  1,  0,  0,  2,  2,
     0,  2,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0
};

template <typename TVoidSpec>
char const TranslateTableAAToRedAA_<
                ClusterReduction<3, 22, Blosum62>, TVoidSpec>::VALUE[24] =
{
     0,  1,  1,  1,  2,  1,  1,  1,  1,  2,  2,  1,
     2,  2,  1,  1,  0,  2,  2,  2,  1,  1,  0,  0
};

template <typename TVoidSpec>
char const TranslateTableByteToRedAA_<
                ClusterReduction<3, 22, Blosum62>, TVoidSpec>::VALUE[256] =
{
     0,  1,  2,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0
};

// ---------------------------------- N = 2 ------------------------------

template <typename TVoidSpec>
char const TranslateTableRedAAToChar_<
                ClusterReduction<2, 22, Blosum62>, TVoidSpec>::VALUE[2] =
{
    'A', // A R N D Q E G H K P S T B Z
    'C'  // C I L M F W Y V
};

template <typename TVoidSpec>
char const TranslateTableCharToRedAA_<
                ClusterReduction<2, 22, Blosum62>, TVoidSpec>::VALUE[256] =
{
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  1,  0,  0,  1,  0,
     0,  1,  1,  0,  0,  0,  0,  0,  0,  0,  0,  1,  1,  0,  1,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  1,  0,  0,  1,  0,  0,
     1,  0,  0,  1,  1,  0,  0,  0,  0,  0,  0,  0,  0,  1,  1,
     0,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0
};

template <typename TVoidSpec>
char const TranslateTableAAToRedAA_<
                ClusterReduction<2, 22, Blosum62>, TVoidSpec>::VALUE[24] =
{
     0,  0,  0,  0,  1,  0,  0,  0,  0,  1,  1,  0,
     1,  1,  0,  0,  0,  1,  1,  1,  0,  0,  0,  0
};

template <typename TVoidSpec>
char const TranslateTableByteToRedAA_<
                ClusterReduction<2, 22, Blosum62>, TVoidSpec>::VALUE[256] =
{
     0,  1,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
     0
};

// ============================================================================
// Functions
// ============================================================================

} // namespace

#endif // SEQAN_REDUCED_AMINOACID_CLUSTER_RED_TABLES_22_to_n_B62_H_
