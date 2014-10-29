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
// Author: Andreas Gogol-Doering <andreas.doering@mdc-berlin.de>
// Author: Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>
// ==========================================================================
// Static data definitions of the following matrices:
//
//  - BLOSUM30, BLOSUM45, BLOSUM62 and BLOSUM80
//  - PAM40, PAM120, PAM200, PAM250; choice following [Altschul, 1991].
//  - VTML200; choice following [Edgar, 2009].
//
// [Altschul, 1991]  Altschul SF.  Amino acid substitution matrices from an
// information theoretic perspective.  Journal of molecular biology.
// 1991;219(3):555-65.
//
// [Edgar, 2009]  Edgar RC.  Optimizing substitution matrix choice and gap
// parameters for sequence alignment.  BMC bioinformatics.  2009;10:396.
//
// Note that there is a script mat2cpp.py that allows the easy conversion of
// scoring matrices into C++ fragments.
// ==========================================================================

// TODO(holtgrew): Maybe also set gap penalties when setting matrices?

#ifndef SEQAN_SCORE_SCORE_MATRIX_DATA_H_
#define SEQAN_SCORE_SCORE_MATRIX_DATA_H_

namespace SEQAN_NAMESPACE_MAIN {

/*
.Tag.Blosum30_:
..cat:Scoring
..summary:Tag for Retrieving a BLOSUM30 matrix.
..include:seqan/score.h
 */
struct Blosum30_ {};


/**
.Shortcut.Blosum30:
..cat:Scoring
..summary:Blosum30 scoring matrix.
..signature:Blosum30
..shortcutfor:Spec.Score Matrix
...signature:Score<int, ScoreMatrix<AminoAcid, Blosum30_> >
..include:seqan/score.h
*/
typedef Score<int, ScoreMatrix<AminoAcid, Blosum30_> > Blosum30;


template <>
struct ScoringMatrixData_<int, AminoAcid, Blosum30_> {
    enum {
        VALUE_SIZE = ValueSize<AminoAcid>::VALUE,
        TAB_SIZE = VALUE_SIZE * VALUE_SIZE
    };

    static inline int const * getData() {
        SEQAN_CHECKPOINT;
        // The matrix data, ordered by amino acid alphabet.
        // Matrix made by matblas from blosum30.iij
        // * column uses minimum score
        // BLOSUM Clustered Scoring Matrix in 1/5 Bit Units
        // Blocks Database = /data/blocks_5.0/blocks.dat
        // Cluster Percentage: >= 30
        // Entropy =   0.1424, Expected =  -0.1074
        static int const _data[TAB_SIZE] = {
             4, -1,  0,  0, -3,  1,  0,  0, -2,  0, -1,  0,  1, -2, -1,  1,  1, -5, -4,  1,  0,  0,  0, -7,
            -1,  8, -2, -1, -2,  3, -1, -2, -1, -3, -2,  1,  0, -1, -1, -1, -3,  0,  0, -1, -2,  0, -1, -7,
             0, -2,  8,  1, -1, -1, -1,  0, -1,  0, -2,  0,  0, -1, -3,  0,  1, -7, -4, -2,  4, -1,  0, -7,
             0, -1,  1,  9, -3, -1,  1, -1, -2, -4, -1,  0, -3, -5, -1,  0, -1, -4, -1, -2,  5,  0, -1, -7,
            -3, -2, -1, -3, 17, -2,  1, -4, -5, -2,  0, -3, -2, -3, -3, -2, -2, -2, -6, -2, -2,  0, -2, -7,
             1,  3, -1, -1, -2,  8,  2, -2,  0, -2, -2,  0, -1, -3,  0, -1,  0, -1, -1, -3, -1,  4,  0, -7,
             0, -1, -1,  1,  1,  2,  6, -2,  0, -3, -1,  2, -1, -4,  1,  0, -2, -1, -2, -3,  0,  5, -1, -7,
             0, -2,  0, -1, -4, -2, -2,  8, -3, -1, -2, -1, -2, -3, -1,  0, -2,  1, -3, -3,  0, -2, -1, -7,
            -2, -1, -1, -2, -5,  0,  0, -3, 14, -2, -1, -2,  2, -3,  1, -1, -2, -5,  0, -3, -2,  0, -1, -7,
             0, -3,  0, -4, -2, -2, -3, -1, -2,  6,  2, -2,  1,  0, -3, -1,  0, -3, -1,  4, -2, -3,  0, -7,
            -1, -2, -2, -1,  0, -2, -1, -2, -1,  2,  4, -2,  2,  2, -3, -2,  0, -2,  3,  1, -1, -1,  0, -7,
             0,  1,  0,  0, -3,  0,  2, -1, -2, -2, -2,  4,  2, -1,  1,  0, -1, -2, -1, -2,  0,  1,  0, -7,
             1,  0,  0, -3, -2, -1, -1, -2,  2,  1,  2,  2,  6, -2, -4, -2,  0, -3, -1,  0, -2, -1,  0, -7,
            -2, -1, -1, -5, -3, -3, -4, -3, -3,  0,  2, -1, -2, 10, -4, -1, -2,  1,  3,  1, -3, -4, -1, -7,
            -1, -1, -3, -1, -3,  0,  1, -1,  1, -3, -3,  1, -4, -4, 11, -1,  0, -3, -2, -4, -2,  0, -1, -7,
             1, -1,  0,  0, -2, -1,  0,  0, -1, -1, -2,  0, -2, -1, -1,  4,  2, -3, -2, -1,  0, -1,  0, -7,
             1, -3,  1, -1, -2,  0, -2, -2, -2,  0,  0, -1,  0, -2,  0,  2,  5, -5, -1,  1,  0, -1,  0, -7,
            -5,  0, -7, -4, -2, -1, -1,  1, -5, -3, -2, -2, -3,  1, -3, -3, -5, 20,  5, -3, -5, -1, -2, -7,
            -4,  0, -4, -1, -6, -1, -2, -3,  0, -1,  3, -1, -1,  3, -2, -2, -1,  5,  9,  1, -3, -2, -1, -7,
             1, -1, -2, -2, -2, -3, -3, -3, -3,  4,  1, -2,  0,  1, -4, -1,  1, -3,  1,  5, -2, -3,  0, -7,
             0, -2,  4,  5, -2, -1,  0,  0, -2, -2, -1,  0, -2, -3, -2,  0,  0, -5, -3, -2,  5,  0, -1, -7,
             0,  0, -1,  0,  0,  4,  5, -2,  0, -3, -1,  1, -1, -4,  0, -1, -1, -1, -2, -3,  0,  4,  0, -7,
             0, -1,  0, -1, -2,  0, -1, -1, -1,  0,  0,  0,  0, -1, -1,  0,  0, -2, -1,  0, -1,  0, -1, -7,
            -7, -7, -7, -7, -7, -7, -7, -7, -7, -7, -7, -7, -7, -7, -7, -7, -7, -7, -7, -7, -7, -7, -7,  1,
        };
        return _data;
    }
};


/*
.Tag.Blosum45_:
..cat:Scoring
..summary:Tag for Retrieving a BLOSUM45 matrix.
..include:seqan/score.h
 */
struct Blosum45_ {};


/**
.Shortcut.Blosum45:
..cat:Scoring
..summary:Blosum45 scoring matrix.
..signature:Blosum45
..shortcutfor:Spec.Score Matrix
...signature:Score<int, ScoreMatrix<AminoAcid, Blosum45_> >
..include:seqan/score.h
*/
typedef Score<int, ScoreMatrix<AminoAcid, Blosum45_> > Blosum45;


template <>
struct ScoringMatrixData_<int, AminoAcid, Blosum45_> {
    enum {
        VALUE_SIZE = ValueSize<AminoAcid>::VALUE,
        TAB_SIZE = VALUE_SIZE * VALUE_SIZE
    };

    static inline int const * getData() {
        SEQAN_CHECKPOINT;
        // Matrix made by matblas from blosum45.iij
        // * column uses minimum score
        // BLOSUM Clustered Scoring Matrix in 1/3 Bit Units
        // Blocks Database = /data/blocks_5.0/blocks.dat
        // Cluster Percentage: >= 45
        // Entropy =   0.3795, Expected =  -0.2789
        static int const _data[TAB_SIZE] = {
             5, -2, -1, -2, -1, -1, -1,  0, -2, -1, -1, -1, -1, -2, -1,  1,  0, -2, -2,  0, -1, -1,  0, -5,
            -2,  7,  0, -1, -3,  1,  0, -2,  0, -3, -2,  3, -1, -2, -2, -1, -1, -2, -1, -2, -1,  0, -1, -5,
            -1,  0,  6,  2, -2,  0,  0,  0,  1, -2, -3,  0, -2, -2, -2,  1,  0, -4, -2, -3,  4,  0, -1, -5,
            -2, -1,  2,  7, -3,  0,  2, -1,  0, -4, -3,  0, -3, -4, -1,  0, -1, -4, -2, -3,  5,  1, -1, -5,
            -1, -3, -2, -3, 12, -3, -3, -3, -3, -3, -2, -3, -2, -2, -4, -1, -1, -5, -3, -1, -2, -3, -2, -5,
            -1,  1,  0,  0, -3,  6,  2, -2,  1, -2, -2,  1,  0, -4, -1,  0, -1, -2, -1, -3,  0,  4, -1, -5,
            -1,  0,  0,  2, -3,  2,  6, -2,  0, -3, -2,  1, -2, -3,  0,  0, -1, -3, -2, -3,  1,  4, -1, -5,
             0, -2,  0, -1, -3, -2, -2,  7, -2, -4, -3, -2, -2, -3, -2,  0, -2, -2, -3, -3, -1, -2, -1, -5,
            -2,  0,  1,  0, -3,  1,  0, -2, 10, -3, -2, -1,  0, -2, -2, -1, -2, -3,  2, -3,  0,  0, -1, -5,
            -1, -3, -2, -4, -3, -2, -3, -4, -3,  5,  2, -3,  2,  0, -2, -2, -1, -2,  0,  3, -3, -3, -1, -5,
            -1, -2, -3, -3, -2, -2, -2, -3, -2,  2,  5, -3,  2,  1, -3, -3, -1, -2,  0,  1, -3, -2, -1, -5,
            -1,  3,  0,  0, -3,  1,  1, -2, -1, -3, -3,  5, -1, -3, -1, -1, -1, -2, -1, -2,  0,  1, -1, -5,
            -1, -1, -2, -3, -2,  0, -2, -2,  0,  2,  2, -1,  6,  0, -2, -2, -1, -2,  0,  1, -2, -1, -1, -5,
            -2, -2, -2, -4, -2, -4, -3, -3, -2,  0,  1, -3,  0,  8, -3, -2, -1,  1,  3,  0, -3, -3, -1, -5,
            -1, -2, -2, -1, -4, -1,  0, -2, -2, -2, -3, -1, -2, -3,  9, -1, -1, -3, -3, -3, -2, -1, -1, -5,
             1, -1,  1,  0, -1,  0,  0,  0, -1, -2, -3, -1, -2, -2, -1,  4,  2, -4, -2, -1,  0,  0,  0, -5,
             0, -1,  0, -1, -1, -1, -1, -2, -2, -1, -1, -1, -1, -1, -1,  2,  5, -3, -1,  0,  0, -1,  0, -5,
            -2, -2, -4, -4, -5, -2, -3, -2, -3, -2, -2, -2, -2,  1, -3, -4, -3, 15,  3, -3, -4, -2, -2, -5,
            -2, -1, -2, -2, -3, -1, -2, -3,  2,  0,  0, -1,  0,  3, -3, -2, -1,  3,  8, -1, -2, -2, -1, -5,
             0, -2, -3, -3, -1, -3, -3, -3, -3,  3,  1, -2,  1,  0, -3, -1,  0, -3, -1,  5, -3, -3, -1, -5,
            -1, -1,  4,  5, -2,  0,  1, -1,  0, -3, -3,  0, -2, -3, -2,  0,  0, -4, -2, -3,  4,  2, -1, -5,
            -1,  0,  0,  1, -3,  4,  4, -2,  0, -3, -2,  1, -1, -3, -1,  0, -1, -2, -2, -3,  2,  4, -1, -5,
             0, -1, -1, -1, -2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1,  0,  0, -2, -1, -1, -1, -1, -1, -5,
            -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5, -5,  1,
        };
        return _data;
    }
};


/*
.Tag.Blosum62_:
..cat:Scoring
..summary:Tag for Retrieving a BLOSUM62 matrix.
..include:seqan/score.h
 */
struct Blosum62_ {};


/**
.Shortcut.Blosum62:
..cat:Scoring
..summary:Blosum62 scoring matrix.
..signature:Blosum62
..shortcutfor:Spec.Score Matrix
...signature:Score<int, ScoreMatrix<AminoAcid, Blosum62_> >
..include:seqan/score.h
 */
typedef Score<int, ScoreMatrix<AminoAcid, Blosum62_> > Blosum62;


template <>
struct ScoringMatrixData_<int, AminoAcid, Blosum62_> {
    enum {
        VALUE_SIZE = ValueSize<AminoAcid>::VALUE,
        TAB_SIZE = VALUE_SIZE * VALUE_SIZE
    };

    static inline int const * getData() {
        SEQAN_CHECKPOINT;
        // Matrix made by matblas from blosum62.iij
        // * column uses minimum score
        // BLOSUM Clustered Scoring Matrix in 1/2 Bit Units
        // Blocks Database = /data/blocks_5.0/blocks.dat
        // Cluster Percentage: >= 62
        // Entropy =   0.6979, Expected =  -0.5209
        static int const _data[TAB_SIZE] = {
             4, -1, -2, -2,  0, -1, -1,  0, -2, -1, -1, -1, -1, -2, -1,  1,  0, -3, -2,  0, -2, -1,  0, -4,
            -1,  5,  0, -2, -3,  1,  0, -2,  0, -3, -2,  2, -1, -3, -2, -1, -1, -3, -2, -3, -1,  0, -1, -4,
            -2,  0,  6,  1, -3,  0,  0,  0,  1, -3, -3,  0, -2, -3, -2,  1,  0, -4, -2, -3,  3,  0, -1, -4,
            -2, -2,  1,  6, -3,  0,  2, -1, -1, -3, -4, -1, -3, -3, -1,  0, -1, -4, -3, -3,  4,  1, -1, -4,
             0, -3, -3, -3,  9, -3, -4, -3, -3, -1, -1, -3, -1, -2, -3, -1, -1, -2, -2, -1, -3, -3, -2, -4,
            -1,  1,  0,  0, -3,  5,  2, -2,  0, -3, -2,  1,  0, -3, -1,  0, -1, -2, -1, -2,  0,  3, -1, -4,
            -1,  0,  0,  2, -4,  2,  5, -2,  0, -3, -3,  1, -2, -3, -1,  0, -1, -3, -2, -2,  1,  4, -1, -4,
             0, -2,  0, -1, -3, -2, -2,  6, -2, -4, -4, -2, -3, -3, -2,  0, -2, -2, -3, -3, -1, -2, -1, -4,
            -2,  0,  1, -1, -3,  0,  0, -2,  8, -3, -3, -1, -2, -1, -2, -1, -2, -2,  2, -3,  0,  0, -1, -4,
            -1, -3, -3, -3, -1, -3, -3, -4, -3,  4,  2, -3,  1,  0, -3, -2, -1, -3, -1,  3, -3, -3, -1, -4,
            -1, -2, -3, -4, -1, -2, -3, -4, -3,  2,  4, -2,  2,  0, -3, -2, -1, -2, -1,  1, -4, -3, -1, -4,
            -1,  2,  0, -1, -3,  1,  1, -2, -1, -3, -2,  5, -1, -3, -1,  0, -1, -3, -2, -2,  0,  1, -1, -4,
            -1, -1, -2, -3, -1,  0, -2, -3, -2,  1,  2, -1,  5,  0, -2, -1, -1, -1, -1,  1, -3, -1, -1, -4,
            -2, -3, -3, -3, -2, -3, -3, -3, -1,  0,  0, -3,  0,  6, -4, -2, -2,  1,  3, -1, -3, -3, -1, -4,
            -1, -2, -2, -1, -3, -1, -1, -2, -2, -3, -3, -1, -2, -4,  7, -1, -1, -4, -3, -2, -2, -1, -2, -4,
             1, -1,  1,  0, -1,  0,  0,  0, -1, -2, -2,  0, -1, -2, -1,  4,  1, -3, -2, -2,  0,  0,  0, -4,
             0, -1,  0, -1, -1, -1, -1, -2, -2, -1, -1, -1, -1, -2, -1,  1,  5, -2, -2,  0, -1, -1,  0, -4,
            -3, -3, -4, -4, -2, -2, -3, -2, -2, -3, -2, -3, -1,  1, -4, -3, -2, 11,  2, -3, -4, -3, -2, -4,
            -2, -2, -2, -3, -2, -1, -2, -3,  2, -1, -1, -2, -1,  3, -3, -2, -2,  2,  7, -1, -3, -2, -1, -4,
             0, -3, -3, -3, -1, -2, -2, -3, -3,  3,  1, -2,  1, -1, -2, -2,  0, -3, -1,  4, -3, -2, -1, -4,
            -2, -1,  3,  4, -3,  0,  1, -1,  0, -3, -4,  0, -3, -3, -2,  0, -1, -4, -3, -3,  4,  1, -1, -4,
            -1,  0,  0,  1, -3,  3,  4, -2,  0, -3, -3,  1, -1, -3, -1,  0, -1, -3, -2, -2,  1,  4, -1, -4,
             0, -1, -1, -1, -2, -1, -1, -1, -1, -1, -1, -1, -1, -1, -2,  0,  0, -2, -1, -1, -1, -1, -1, -4,
            -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4, -4,  1,
        };
        return _data;
    }
};


/*
.Tag.Blosum80_:
..cat:Scoring
..summary:Tag for Retrieving a BLOSUM80 matrix.
..include:seqan/score.h
 */
struct Blosum80_ {};


/**
.Shortcut.Blosum80:
..cat:Scoring
..summary:Blosum80 scoring matrix.
..signature:Blosum80
..shortcutfor:Spec.Score Matrix
...signature:Score<int, ScoreMatrix<AminoAcid, Blosum80_> >
..include:seqan/score.h
*/
typedef Score<int, ScoreMatrix<AminoAcid, Blosum80_> > Blosum80;


template <>
struct ScoringMatrixData_<int, AminoAcid, Blosum80_> {
    enum {
        VALUE_SIZE = ValueSize<AminoAcid>::VALUE,
        TAB_SIZE = VALUE_SIZE * VALUE_SIZE
    };

    static inline int const * getData() {
        SEQAN_CHECKPOINT;
        // Matrix made by matblas from blosum80_3.iij
        // * column uses minimum score
        // BLOSUM Clustered Scoring Matrix in 1/3 Bit Units
        // Blocks Database = /data/blocks_5.0/blocks.dat
        // Cluster Percentage: >= 80
        // Entropy =   0.9868, Expected =  -0.7442
        static int const _data[TAB_SIZE] = {
             7, -3, -3, -3, -1, -2, -2,  0, -3, -3, -3, -1, -2, -4, -1,  2,  0, -5, -4, -1, -3, -2, -1, -8,
            -3,  9, -1, -3, -6,  1, -1, -4,  0, -5, -4,  3, -3, -5, -3, -2, -2, -5, -4, -4, -2,  0, -2, -8,
            -3, -1,  9,  2, -5,  0, -1, -1,  1, -6, -6,  0, -4, -6, -4,  1,  0, -7, -4, -5,  5, -1, -2, -8,
            -3, -3,  2, 10, -7, -1,  2, -3, -2, -7, -7, -2, -6, -6, -3, -1, -2, -8, -6, -6,  6,  1, -3, -8,
            -1, -6, -5, -7, 13, -5, -7, -6, -7, -2, -3, -6, -3, -4, -6, -2, -2, -5, -5, -2, -6, -7, -4, -8,
            -2,  1,  0, -1, -5,  9,  3, -4,  1, -5, -4,  2, -1, -5, -3, -1, -1, -4, -3, -4, -1,  5, -2, -8,
            -2, -1, -1,  2, -7,  3,  8, -4,  0, -6, -6,  1, -4, -6, -2, -1, -2, -6, -5, -4,  1,  6, -2, -8,
             0, -4, -1, -3, -6, -4, -4,  9, -4, -7, -7, -3, -5, -6, -5, -1, -3, -6, -6, -6, -2, -4, -3, -8,
            -3,  0,  1, -2, -7,  1,  0, -4, 12, -6, -5, -1, -4, -2, -4, -2, -3, -4,  3, -5, -1,  0, -2, -8,
            -3, -5, -6, -7, -2, -5, -6, -7, -6,  7,  2, -5,  2, -1, -5, -4, -2, -5, -3,  4, -6, -6, -2, -8,
            -3, -4, -6, -7, -3, -4, -6, -7, -5,  2,  6, -4,  3,  0, -5, -4, -3, -4, -2,  1, -7, -5, -2, -8,
            -1,  3,  0, -2, -6,  2,  1, -3, -1, -5, -4,  8, -3, -5, -2, -1, -1, -6, -4, -4, -1,  1, -2, -8,
            -2, -3, -4, -6, -3, -1, -4, -5, -4,  2,  3, -3,  9,  0, -4, -3, -1, -3, -3,  1, -5, -3, -2, -8,
            -4, -5, -6, -6, -4, -5, -6, -6, -2, -1,  0, -5,  0, 10, -6, -4, -4,  0,  4, -2, -6, -6, -3, -8,
            -1, -3, -4, -3, -6, -3, -2, -5, -4, -5, -5, -2, -4, -6, 12, -2, -3, -7, -6, -4, -4, -2, -3, -8,
             2, -2,  1, -1, -2, -1, -1, -1, -2, -4, -4, -1, -3, -4, -2,  7,  2, -6, -3, -3,  0, -1, -1, -8,
             0, -2,  0, -2, -2, -1, -2, -3, -3, -2, -3, -1, -1, -4, -3,  2,  8, -5, -3,  0, -1, -2, -1, -8,
            -5, -5, -7, -8, -5, -4, -6, -6, -4, -5, -4, -6, -3,  0, -7, -6, -5, 16,  3, -5, -8, -5, -5, -8,
            -4, -4, -4, -6, -5, -3, -5, -6,  3, -3, -2, -4, -3,  4, -6, -3, -3,  3, 11, -3, -5, -4, -3, -8,
            -1, -4, -5, -6, -2, -4, -4, -6, -5,  4,  1, -4,  1, -2, -4, -3,  0, -5, -3,  7, -6, -4, -2, -8,
            -3, -2,  5,  6, -6, -1,  1, -2, -1, -6, -7, -1, -5, -6, -4,  0, -1, -8, -5, -6,  6,  0, -3, -8,
            -2,  0, -1,  1, -7,  5,  6, -4,  0, -6, -5,  1, -3, -6, -2, -1, -2, -5, -4, -4,  0,  6, -1, -8,
            -1, -2, -2, -3, -4, -2, -2, -3, -2, -2, -2, -2, -2, -3, -3, -1, -1, -5, -3, -2, -3, -1, -2, -8,
            -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8, -8,  1,
        };
        return _data;
    }
};


/*
.Tag.Pam40_:
..cat:Scoring
..summary:Tag for Retrieving a PAM40 matrix.
..include:seqan/score.h
 */
struct Pam40_ {};


/**
.Shortcut.Pam40:
..cat:Scoring
..summary:Pam40 scoring matrix.
..signature:Pam40
..shortcutfor:Spec.Score Matrix
...signature:Score<int, ScoreMatrix<AminoAcid, Pam40_> >
..include:seqan/score.h
*/
typedef Score<int, ScoreMatrix<AminoAcid, Pam40_> > Pam40;


template <>
struct ScoringMatrixData_<int, AminoAcid, Pam40_> {
    enum {
        VALUE_SIZE = ValueSize<AminoAcid>::VALUE,
        TAB_SIZE = VALUE_SIZE * VALUE_SIZE
    };

    static inline int const * getData() {
        SEQAN_CHECKPOINT;
        // 
        // This matrix was produced by "pam" Version 1.0.6 [28-Jul-93]
        // 
        // PAM 40 substitution matrix, scale = ln(2)/2 = 0.346574
        // 
        // Expected score = -4.27, Entropy = 2.26 bits
        // 
        // Lowest score = -15, Highest score = 13
        // 
        static int const _data[TAB_SIZE] = {
              6,  -6,  -3,  -3,  -6,  -3,  -2,  -1,  -6,  -4,  -5,  -6,  -4,  -7,  -1,   0,   0, -12,  -7,  -2,  -3,  -2,  -3, -15,
             -6,   8,  -5,  -9,  -7,  -1,  -8,  -8,  -1,  -5,  -8,   1,  -3,  -8,  -3,  -2,  -5,  -1,  -9,  -7,  -6,  -3,  -5, -15,
             -3,  -5,   7,   2,  -9,  -3,  -1,  -2,   1,  -4,  -6,   0,  -7,  -8,  -5,   0,  -1,  -7,  -4,  -7,   6,  -2,  -3, -15,
             -3,  -9,   2,   7, -12,  -2,   3,  -3,  -3,  -6, -11,  -4,  -9, -13,  -7,  -3,  -4, -13, -10,  -7,   6,   2,  -5, -15,
             -6,  -7,  -9, -12,   9, -12, -12,  -8,  -7,  -5, -13, -12, -12, -11,  -7,  -2,  -7, -14,  -3,  -5, -11, -12,  -8, -15,
             -3,  -1,  -3,  -2, -12,   8,   2,  -6,   1,  -7,  -4,  -2,  -3, -11,  -2,  -4,  -5, -11, -10,  -6,  -2,   6,  -4, -15,
             -2,  -8,  -1,   3, -12,   2,   7,  -3,  -4,  -5,  -8,  -4,  -6, -12,  -5,  -4,  -5, -15,  -8,  -6,   2,   6,  -4, -15,
             -1,  -8,  -2,  -3,  -8,  -6,  -3,   6,  -8,  -9,  -9,  -6,  -7,  -8,  -5,  -1,  -5, -13, -12,  -5,  -2,  -4,  -4, -15,
             -6,  -1,   1,  -3,  -7,   1,  -4,  -8,   9,  -8,  -5,  -5,  -9,  -5,  -3,  -5,  -6,  -6,  -3,  -6,  -1,   0,  -4, -15,
             -4,  -5,  -4,  -6,  -5,  -7,  -5,  -9,  -8,   8,  -1,  -5,   0,  -2,  -7,  -6,  -2, -12,  -5,   2,  -5,  -5,  -4, -15,
             -5,  -8,  -6, -11, -13,  -4,  -8,  -9,  -5,  -1,   7,  -7,   1,  -2,  -6,  -7,  -6,  -5,  -6,  -2,  -8,  -6,  -5, -15,
             -6,   1,   0,  -4, -12,  -2,  -4,  -6,  -5,  -5,  -7,   6,  -1, -12,  -6,  -3,  -2, -10,  -8,  -8,  -2,  -3,  -4, -15,
             -4,  -3,  -7,  -9, -12,  -3,  -6,  -7,  -9,   0,   1,  -1,  11,  -3,  -7,  -5,  -3, -11, -10,  -1,  -8,  -4,  -4, -15,
             -7,  -8,  -8, -13, -11, -11, -12,  -8,  -5,  -2,  -2, -12,  -3,   9,  -9,  -6,  -8,  -4,   2,  -7,  -9, -12,  -7, -15,
             -1,  -3,  -5,  -7,  -7,  -2,  -5,  -5,  -3,  -7,  -6,  -6,  -7,  -9,   8,  -1,  -3, -12, -12,  -5,  -6,  -3,  -4, -15,
              0,  -2,   0,  -3,  -2,  -4,  -4,  -1,  -5,  -6,  -7,  -3,  -5,  -6,  -1,   6,   1,  -4,  -6,  -5,  -1,  -4,  -2, -15,
              0,  -5,  -1,  -4,  -7,  -5,  -5,  -5,  -6,  -2,  -6,  -2,  -3,  -8,  -3,   1,   7, -11,  -6,  -2,  -2,  -5,  -3, -15,
            -12,  -1,  -7, -13, -14, -11, -15, -13,  -6, -12,  -5, -10, -11,  -4, -12,  -4, -11,  13,  -4, -14,  -9, -13,  -9, -15,
             -7,  -9,  -4, -10,  -3, -10,  -8, -12,  -3,  -5,  -6,  -8, -10,   2, -12,  -6,  -6,  -4,  10,  -6,  -6,  -8,  -7, -15,
             -2,  -7,  -7,  -7,  -5,  -6,  -6,  -5,  -6,   2,  -2,  -8,  -1,  -7,  -5,  -5,  -2, -14,  -6,   7,  -7,  -6,  -4, -15,
             -3,  -6,   6,   6, -11,  -2,   2,  -2,  -1,  -5,  -8,  -2,  -8,  -9,  -6,  -1,  -2,  -9,  -6,  -7,   6,   1,  -4, -15,
             -2,  -3,  -2,   2, -12,   6,   6,  -4,   0,  -5,  -6,  -3,  -4, -12,  -3,  -4,  -5, -13,  -8,  -6,   1,   6,  -4, -15,
             -3,  -5,  -3,  -5,  -8,  -4,  -4,  -4,  -4,  -4,  -5,  -4,  -4,  -7,  -4,  -2,  -3,  -9,  -7,  -4,  -4,  -4,  -4, -15,
            -15, -15, -15, -15, -15, -15, -15, -15, -15, -15, -15, -15, -15, -15, -15, -15, -15, -15, -15, -15, -15, -15, -15,   1,
        };
        return _data;
    }
};


/*
.Tag.Pam120_:
..cat:Scoring
..summary:Tag for Retrieving a PAM120 matrix.
..include:seqan/score.h
 */
struct Pam120_ {};


/**
.Shortcut.Pam120:
..cat:Scoring
..summary:Pam120 scoring matrix.
..signature:Pam120
..shortcutfor:Spec.Score Matrix
...signature:Score<int, ScoreMatrix<AminoAcid, Pam120_> >
..include:seqan/score.h
*/
typedef Score<int, ScoreMatrix<AminoAcid, Pam120_> > Pam120;


template <>
struct ScoringMatrixData_<int, AminoAcid, Pam120_> {
    enum {
        VALUE_SIZE = ValueSize<AminoAcid>::VALUE,
        TAB_SIZE = VALUE_SIZE * VALUE_SIZE
    };

    static inline int const * getData() {
        SEQAN_CHECKPOINT;
        // This matrix was produced by "pam" Version 1.0.6 [28-Jul-93]
        // 
        // PAM 120 substitution matrix, scale = ln(2)/2 = 0.346574
        // 
        // Expected score = -1.64, Entropy = 0.979 bits
        // 
        // Lowest score = -8, Highest score = 12
        // 
        static int const _data[TAB_SIZE] = {
              3,  -3,  -1,   0,  -3,  -1,   0,   1,  -3,  -1,  -3,  -2,  -2,  -4,   1,   1,   1,  -7,  -4,   0,   0,  -1,  -1,  -8,
             -3,   6,  -1,  -3,  -4,   1,  -3,  -4,   1,  -2,  -4,   2,  -1,  -5,  -1,  -1,  -2,   1,  -5,  -3,  -2,  -1,  -2,  -8,
             -1,  -1,   4,   2,  -5,   0,   1,   0,   2,  -2,  -4,   1,  -3,  -4,  -2,   1,   0,  -4,  -2,  -3,   3,   0,  -1,  -8,
              0,  -3,   2,   5,  -7,   1,   3,   0,   0,  -3,  -5,  -1,  -4,  -7,  -3,   0,  -1,  -8,  -5,  -3,   4,   3,  -2,  -8,
             -3,  -4,  -5,  -7,   9,  -7,  -7,  -4,  -4,  -3,  -7,  -7,  -6,  -6,  -4,   0,  -3,  -8,  -1,  -3,  -6,  -7,  -4,  -8,
             -1,   1,   0,   1,  -7,   6,   2,  -3,   3,  -3,  -2,   0,  -1,  -6,   0,  -2,  -2,  -6,  -5,  -3,   0,   4,  -1,  -8,
              0,  -3,   1,   3,  -7,   2,   5,  -1,  -1,  -3,  -4,  -1,  -3,  -7,  -2,  -1,  -2,  -8,  -5,  -3,   3,   4,  -1,  -8,
              1,  -4,   0,   0,  -4,  -3,  -1,   5,  -4,  -4,  -5,  -3,  -4,  -5,  -2,   1,  -1,  -8,  -6,  -2,   0,  -2,  -2,  -8,
             -3,   1,   2,   0,  -4,   3,  -1,  -4,   7,  -4,  -3,  -2,  -4,  -3,  -1,  -2,  -3,  -3,  -1,  -3,   1,   1,  -2,  -8,
             -1,  -2,  -2,  -3,  -3,  -3,  -3,  -4,  -4,   6,   1,  -3,   1,   0,  -3,  -2,   0,  -6,  -2,   3,  -3,  -3,  -1,  -8,
             -3,  -4,  -4,  -5,  -7,  -2,  -4,  -5,  -3,   1,   5,  -4,   3,   0,  -3,  -4,  -3,  -3,  -2,   1,  -4,  -3,  -2,  -8,
             -2,   2,   1,  -1,  -7,   0,  -1,  -3,  -2,  -3,  -4,   5,   0,  -7,  -2,  -1,  -1,  -5,  -5,  -4,   0,  -1,  -2,  -8,
             -2,  -1,  -3,  -4,  -6,  -1,  -3,  -4,  -4,   1,   3,   0,   8,  -1,  -3,  -2,  -1,  -6,  -4,   1,  -4,  -2,  -2,  -8,
             -4,  -5,  -4,  -7,  -6,  -6,  -7,  -5,  -3,   0,   0,  -7,  -1,   8,  -5,  -3,  -4,  -1,   4,  -3,  -5,  -6,  -3,  -8,
              1,  -1,  -2,  -3,  -4,   0,  -2,  -2,  -1,  -3,  -3,  -2,  -3,  -5,   6,   1,  -1,  -7,  -6,  -2,  -2,  -1,  -2,  -8,
              1,  -1,   1,   0,   0,  -2,  -1,   1,  -2,  -2,  -4,  -1,  -2,  -3,   1,   3,   2,  -2,  -3,  -2,   0,  -1,  -1,  -8,
              1,  -2,   0,  -1,  -3,  -2,  -2,  -1,  -3,   0,  -3,  -1,  -1,  -4,  -1,   2,   4,  -6,  -3,   0,   0,  -2,  -1,  -8,
             -7,   1,  -4,  -8,  -8,  -6,  -8,  -8,  -3,  -6,  -3,  -5,  -6,  -1,  -7,  -2,  -6,  12,  -2,  -8,  -6,  -7,  -5,  -8,
             -4,  -5,  -2,  -5,  -1,  -5,  -5,  -6,  -1,  -2,  -2,  -5,  -4,   4,  -6,  -3,  -3,  -2,   8,  -3,  -3,  -5,  -3,  -8,
              0,  -3,  -3,  -3,  -3,  -3,  -3,  -2,  -3,   3,   1,  -4,   1,  -3,  -2,  -2,   0,  -8,  -3,   5,  -3,  -3,  -1,  -8,
              0,  -2,   3,   4,  -6,   0,   3,   0,   1,  -3,  -4,   0,  -4,  -5,  -2,   0,   0,  -6,  -3,  -3,   4,   2,  -1,  -8,
             -1,  -1,   0,   3,  -7,   4,   4,  -2,   1,  -3,  -3,  -1,  -2,  -6,  -1,  -1,  -2,  -7,  -5,  -3,   2,   4,  -1,  -8,
             -1,  -2,  -1,  -2,  -4,  -1,  -1,  -2,  -2,  -1,  -2,  -2,  -2,  -3,  -2,  -1,  -1,  -5,  -3,  -1,  -1,  -1,  -2,  -8,
             -8,  -8,  -8,  -8,  -8,  -8,  -8,  -8,  -8,  -8,  -8,  -8,  -8,  -8,  -8,  -8,  -8,  -8,  -8,  -8,  -8,  -8,  -8,   1,
        };
        return _data;
    }
};


/*
.Tag.Pam200_:
..cat:Scoring
..summary:Tag for Retrieving a PAM200 matrix.
..include:seqan/score.h
 */
struct Pam200_ {};


/**
.Shortcut.Pam200:
..cat:Scoring
..summary:Pam200 scoring matrix.
..signature:Pam200
..shortcutfor:Spec.Score Matrix
...signature:Score<int, ScoreMatrix<AminoAcid, Pam200_> >
..include:seqan/score.h
*/
typedef Score<int, ScoreMatrix<AminoAcid, Pam200_> > Pam200;


template <>
struct ScoringMatrixData_<int, AminoAcid, Pam200_> {
    enum {
        VALUE_SIZE = ValueSize<AminoAcid>::VALUE,
        TAB_SIZE = VALUE_SIZE * VALUE_SIZE
    };

    static inline int const * getData() {
        SEQAN_CHECKPOINT;
        // This matrix was produced by "pam" Version 1.0.6 [28-Jul-93]
        // 
        // PAM 200 substitution matrix, scale = ln(2)/3 = 0.231049
        // 
        // Expected score = -1.23, Entropy = 0.507 bits
        // 
        // Lowest score = -9, Highest score = 18
        // 
        static int const _data[TAB_SIZE] = {
              3,  -2,   0,   0,  -3,  -1,   0,   1,  -2,  -1,  -2,  -2,  -2,  -4,   1,   1,   1,  -7,  -4,   0,   0,   0,   0,  -9,
             -2,   7,   0,  -2,  -4,   1,  -2,  -4,   2,  -2,  -4,   4,  -1,  -5,   0,  -1,  -1,   2,  -5,  -3,  -1,   0,  -1,  -9,
              0,   0,   3,   3,  -5,   1,   2,   0,   2,  -2,  -4,   1,  -2,  -4,  -1,   1,   0,  -5,  -2,  -2,   3,   1,   0,  -9,
              0,  -2,   3,   5,  -6,   2,   4,   0,   0,  -3,  -5,   0,  -4,  -7,  -2,   0,   0,  -8,  -5,  -3,   4,   3,  -1,  -9,
             -3,  -4,  -5,  -6,  12,  -7,  -7,  -4,  -4,  -3,  -7,  -7,  -6,  -6,  -4,   0,  -3,  -9,   0,  -2,  -5,  -7,  -4,  -9,
             -1,   1,   1,   2,  -7,   5,   3,  -2,   3,  -3,  -2,   1,  -1,  -6,   0,  -1,  -1,  -6,  -5,  -3,   1,   4,  -1,  -9,
              0,  -2,   2,   4,  -7,   3,   5,   0,   0,  -3,  -4,   0,  -3,  -7,  -1,   0,  -1,  -9,  -5,  -2,   3,   4,  -1,  -9,
              1,  -4,   0,   0,  -4,  -2,   0,   6,  -3,  -3,  -5,  -2,  -4,  -6,  -1,   1,   0,  -8,  -6,  -2,   0,  -1,  -1,  -9,
             -2,   2,   2,   0,  -4,   3,   0,  -3,   8,  -3,  -3,  -1,  -3,  -2,  -1,  -1,  -2,  -3,   0,  -3,   1,   2,  -1,  -9,
             -1,  -2,  -2,  -3,  -3,  -3,  -3,  -3,  -3,   6,   2,  -2,   2,   1,  -3,  -2,   0,  -6,  -2,   4,  -3,  -3,  -1,  -9,
             -2,  -4,  -4,  -5,  -7,  -2,  -4,  -5,  -3,   2,   7,  -4,   4,   2,  -3,  -4,  -2,  -2,  -2,   2,  -4,  -3,  -2,  -9,
             -2,   4,   1,   0,  -7,   1,   0,  -2,  -1,  -2,  -4,   6,   1,  -7,  -2,   0,   0,  -4,  -5,  -3,   0,   0,  -1,  -9,
             -2,  -1,  -2,  -4,  -6,  -1,  -3,  -4,  -3,   2,   4,   1,   8,   0,  -3,  -2,  -1,  -5,  -3,   2,  -3,  -2,  -1,  -9,
             -4,  -5,  -4,  -7,  -6,  -6,  -7,  -6,  -2,   1,   2,  -7,   0,  10,  -6,  -4,  -4,   0,   7,  -2,  -6,  -6,  -3,  -9,
              1,   0,  -1,  -2,  -4,   0,  -1,  -1,  -1,  -3,  -3,  -2,  -3,  -6,   7,   1,   0,  -7,  -6,  -2,  -1,  -1,  -1,  -9,
              1,  -1,   1,   0,   0,  -1,   0,   1,  -1,  -2,  -4,   0,  -2,  -4,   1,   2,   2,  -3,  -3,  -1,   1,  -1,   0,  -9,
              1,  -1,   0,   0,  -3,  -1,  -1,   0,  -2,   0,  -2,   0,  -1,  -4,   0,   2,   4,  -6,  -3,   0,   0,  -1,   0,  -9,
             -7,   2,  -5,  -8,  -9,  -6,  -9,  -8,  -3,  -6,  -2,  -4,  -5,   0,  -7,  -3,  -6,  18,  -1,  -8,  -6,  -7,  -5,  -9,
             -4,  -5,  -2,  -5,   0,  -5,  -5,  -6,   0,  -2,  -2,  -5,  -3,   7,  -6,  -3,  -3,  -1,  11,  -3,  -4,  -5,  -3,  -9,
              0,  -3,  -2,  -3,  -2,  -3,  -2,  -2,  -3,   4,   2,  -3,   2,  -2,  -2,  -1,   0,  -8,  -3,   5,  -3,  -2,  -1,  -9,
              0,  -1,   3,   4,  -5,   1,   3,   0,   1,  -3,  -4,   0,  -3,  -6,  -1,   1,   0,  -6,  -4,  -3,   3,   2,  -1,  -9,
              0,   0,   1,   3,  -7,   4,   4,  -1,   2,  -3,  -3,   0,  -2,  -6,  -1,  -1,  -1,  -7,  -5,  -2,   2,   4,  -1,  -9,
              0,  -1,   0,  -1,  -4,  -1,  -1,  -1,  -1,  -1,  -2,  -1,  -1,  -3,  -1,   0,   0,  -5,  -3,  -1,  -1,  -1,  -1,  -9,
             -9,  -9,  -9,  -9,  -9,  -9,  -9,  -9,  -9,  -9,  -9,  -9,  -9,  -9,  -9,  -9,  -9,  -9,  -9,  -9,  -9,  -9,  -9,   1,
        };
        return _data;
    }
};


/*
.Tag.Pam250_:
..cat:Scoring
..summary:Tag for Retrieving a PAM250 matrix.
..include:seqan/score.h
 */
struct Pam250_ {};


/**
.Shortcut.Pam250:
..cat:Scoring
..summary:Pam250 scoring matrix.
..signature:Pam250
..shortcutfor:Spec.Score Matrix
...signature:Score<int, ScoreMatrix<AminoAcid, Pam250_> >
..include:seqan/score.h
*/
typedef Score<int, ScoreMatrix<AminoAcid, Pam250_> > Pam250;


template <>
struct ScoringMatrixData_<int, AminoAcid, Pam250_> {
    enum {
        VALUE_SIZE = ValueSize<AminoAcid>::VALUE,
        TAB_SIZE = VALUE_SIZE * VALUE_SIZE
    };

    static inline int const * getData() {
        SEQAN_CHECKPOINT;
        // This matrix was produced by "pam" Version 1.0.6 [28-Jul-93]
        // 
        // PAM 250 substitution matrix, scale = ln(2)/3 = 0.231049
        // 
        // Expected score = -0.844, Entropy = 0.354 bits
        // 
        // Lowest score = -8, Highest score = 17
        // 
        static int const _data[TAB_SIZE] = {
              2,  -2,   0,   0,  -2,   0,   0,   1,  -1,  -1,  -2,  -1,  -1,  -3,   1,   1,   1,  -6,  -3,   0,   0,   0,   0,  -8,
             -2,   6,   0,  -1,  -4,   1,  -1,  -3,   2,  -2,  -3,   3,   0,  -4,   0,   0,  -1,   2,  -4,  -2,  -1,   0,  -1,  -8,
              0,   0,   2,   2,  -4,   1,   1,   0,   2,  -2,  -3,   1,  -2,  -3,   0,   1,   0,  -4,  -2,  -2,   2,   1,   0,  -8,
              0,  -1,   2,   4,  -5,   2,   3,   1,   1,  -2,  -4,   0,  -3,  -6,  -1,   0,   0,  -7,  -4,  -2,   3,   3,  -1,  -8,
             -2,  -4,  -4,  -5,  12,  -5,  -5,  -3,  -3,  -2,  -6,  -5,  -5,  -4,  -3,   0,  -2,  -8,   0,  -2,  -4,  -5,  -3,  -8,
              0,   1,   1,   2,  -5,   4,   2,  -1,   3,  -2,  -2,   1,  -1,  -5,   0,  -1,  -1,  -5,  -4,  -2,   1,   3,  -1,  -8,
              0,  -1,   1,   3,  -5,   2,   4,   0,   1,  -2,  -3,   0,  -2,  -5,  -1,   0,   0,  -7,  -4,  -2,   3,   3,  -1,  -8,
              1,  -3,   0,   1,  -3,  -1,   0,   5,  -2,  -3,  -4,  -2,  -3,  -5,   0,   1,   0,  -7,  -5,  -1,   0,   0,  -1,  -8,
             -1,   2,   2,   1,  -3,   3,   1,  -2,   6,  -2,  -2,   0,  -2,  -2,   0,  -1,  -1,  -3,   0,  -2,   1,   2,  -1,  -8,
             -1,  -2,  -2,  -2,  -2,  -2,  -2,  -3,  -2,   5,   2,  -2,   2,   1,  -2,  -1,   0,  -5,  -1,   4,  -2,  -2,  -1,  -8,
             -2,  -3,  -3,  -4,  -6,  -2,  -3,  -4,  -2,   2,   6,  -3,   4,   2,  -3,  -3,  -2,  -2,  -1,   2,  -3,  -3,  -1,  -8,
             -1,   3,   1,   0,  -5,   1,   0,  -2,   0,  -2,  -3,   5,   0,  -5,  -1,   0,   0,  -3,  -4,  -2,   1,   0,  -1,  -8,
             -1,   0,  -2,  -3,  -5,  -1,  -2,  -3,  -2,   2,   4,   0,   6,   0,  -2,  -2,  -1,  -4,  -2,   2,  -2,  -2,  -1,  -8,
             -3,  -4,  -3,  -6,  -4,  -5,  -5,  -5,  -2,   1,   2,  -5,   0,   9,  -5,  -3,  -3,   0,   7,  -1,  -4,  -5,  -2,  -8,
              1,   0,   0,  -1,  -3,   0,  -1,   0,   0,  -2,  -3,  -1,  -2,  -5,   6,   1,   0,  -6,  -5,  -1,  -1,   0,  -1,  -8,
              1,   0,   1,   0,   0,  -1,   0,   1,  -1,  -1,  -3,   0,  -2,  -3,   1,   2,   1,  -2,  -3,  -1,   0,   0,   0,  -8,
              1,  -1,   0,   0,  -2,  -1,   0,   0,  -1,   0,  -2,   0,  -1,  -3,   0,   1,   3,  -5,  -3,   0,   0,  -1,   0,  -8,
             -6,   2,  -4,  -7,  -8,  -5,  -7,  -7,  -3,  -5,  -2,  -3,  -4,   0,  -6,  -2,  -5,  17,   0,  -6,  -5,  -6,  -4,  -8,
             -3,  -4,  -2,  -4,   0,  -4,  -4,  -5,   0,  -1,  -1,  -4,  -2,   7,  -5,  -3,  -3,   0,  10,  -2,  -3,  -4,  -2,  -8,
              0,  -2,  -2,  -2,  -2,  -2,  -2,  -1,  -2,   4,   2,  -2,   2,  -1,  -1,  -1,   0,  -6,  -2,   4,  -2,  -2,  -1,  -8,
              0,  -1,   2,   3,  -4,   1,   3,   0,   1,  -2,  -3,   1,  -2,  -4,  -1,   0,   0,  -5,  -3,  -2,   3,   2,  -1,  -8,
              0,   0,   1,   3,  -5,   3,   3,   0,   2,  -2,  -3,   0,  -2,  -5,   0,   0,  -1,  -6,  -4,  -2,   2,   3,  -1,  -8,
              0,  -1,   0,  -1,  -3,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -1,  -2,  -1,   0,   0,  -4,  -2,  -1,  -1,  -1,  -1,  -8,
             -8,  -8,  -8,  -8,  -8,  -8,  -8,  -8,  -8,  -8,  -8,  -8,  -8,  -8,  -8,  -8,  -8,  -8,  -8,  -8,  -8,  -8,  -8,   1,
        };
        return _data;
    }
};


/*
.Tag.Vtml200_:
..cat:Scoring
..summary:Tag for Retrieving a PAM200 matrix.
..include:seqan/score.h
 */
struct Vtml200_ {};


/**
.Shortcut.Vtml200:
..cat:Scoring
..summary:Vtml200 scoring matrix.
..signature:Vtml200
..shortcutfor:Spec.Score Matrix
...signature:Score<int, ScoreMatrix<AminoAcid, Vtml200_> >
..include:seqan/score.h
*/
typedef Score<int, ScoreMatrix<AminoAcid, Vtml200_> > Vtml200;


template <>
struct ScoringMatrixData_<int, AminoAcid, Vtml200_> {
    enum {
        VALUE_SIZE = ValueSize<AminoAcid>::VALUE,
        TAB_SIZE = VALUE_SIZE * VALUE_SIZE
    };

    static inline int const * getData() {
        SEQAN_CHECKPOINT;
        // VTML200
        // 
        // This matrix was produced by scripts written by
        // Kai Kneutgen and Tobias Mueller [Mai-2002]
        // 
        // VTML200  substitution matrix, Units = Third-Bits
        // 
        // Expected Score = -1.038 Third-Bits
        // Lowest Score   = -6.432
        // Highest Score  = 15.127
        // Entropy H      = 0.420 Bits
        // 
        // 
        // For further information on the VTML substitution model, see
        // 
        // Estimating Amino Acid Substitution Models:
        // A Comparison of Dayhoff's Estimator, the Resolvent Approach and a Maximum Likelihood Method.
        // T. Mueller, R. Spang and M. Vingron
        // Mol Biol Evol 19(1): 8-13. 2002.
        // 
        // or mail to
        // Tobias.Mueller@molgen.mpg.de
        // 
        // The latest version of this perl script can be downloaded at
        // http://www.molgen.mpg.de/~muelle_t
        static int const _data[TAB_SIZE] = {
              4,  -2,  -1,  -1,   1,  -1,  -1,   0,  -2,  -1,  -2,  -1,  -1,  -3,   0,   1,   1,  -4,  -3,   0,  -1,  -1,  -1,  -6,
             -2,   7,   0,  -2,  -3,   2,  -1,  -2,   1,  -3,  -3,   4,  -2,  -4,  -1,  -1,  -1,  -3,  -2,  -3,  -1,   0,  -1,  -6,
             -1,   0,   6,   3,  -2,   1,   1,   0,   1,  -4,  -4,   1,  -3,  -4,  -2,   1,   0,  -5,  -2,  -3,   4,   1,  -1,  -6,
             -1,  -2,   3,   6,  -4,   1,   3,  -1,   0,  -5,  -5,   0,  -4,  -6,  -1,   0,  -1,  -6,  -4,  -4,   5,   2,  -2,  -6,
              1,  -3,  -2,  -4,  12,  -3,  -4,  -2,  -2,   0,  -3,  -4,  -1,  -3,  -3,   1,   0,  -6,   0,   1,  -3,  -4,  -1,  -6,
             -1,   2,   1,   1,  -3,   5,   2,  -2,   2,  -3,  -2,   2,  -1,  -3,  -1,   0,   0,  -6,  -3,  -2,   1,   4,  -1,  -6,
             -1,  -1,   1,   3,  -4,   2,   5,  -1,   0,  -4,  -4,   1,  -3,  -5,  -1,   0,  -1,  -6,  -3,  -3,   2,   4,  -1,  -6,
              0,  -2,   0,  -1,  -2,  -2,  -1,   8,  -2,  -6,  -5,  -2,  -4,  -5,  -2,   0,  -2,  -5,  -5,  -4,   0,  -2,  -2,  -6,
             -2,   1,   1,   0,  -2,   2,   0,  -2,   8,  -3,  -2,   0,  -3,   0,  -2,   0,  -1,  -1,   3,  -3,   1,   1,   0,  -6,
             -1,  -3,  -4,  -5,   0,  -3,  -4,  -6,  -3,   5,   3,  -3,   2,   0,  -4,  -3,  -1,  -2,  -2,   4,  -4,  -3,  -1,  -6,
             -2,  -3,  -4,  -5,  -3,  -2,  -4,  -5,  -2,   3,   5,  -3,   3,   2,  -3,  -3,  -2,  -1,  -1,   2,  -5,  -3,  -1,  -6,
             -1,   4,   1,   0,  -4,   2,   1,  -2,   0,  -3,  -3,   5,  -2,  -5,  -1,   0,   0,  -4,  -3,  -3,   0,   2,  -1,  -6,
             -1,  -2,  -3,  -4,  -1,  -1,  -3,  -4,  -3,   2,   3,  -2,   6,   1,  -3,  -2,  -1,  -3,  -2,   2,  -3,  -2,  -1,  -6,
             -3,  -4,  -4,  -6,  -3,  -3,  -5,  -5,   0,   0,   2,  -5,   1,   8,  -4,  -3,  -3,   3,   5,  -1,  -5,  -4,  -1,  -6,
              0,  -1,  -2,  -1,  -3,  -1,  -1,  -2,  -2,  -4,  -3,  -1,  -3,  -4,   9,   0,  -1,  -4,  -5,  -3,  -1,  -1,  -2,  -6,
              1,  -1,   1,   0,   1,   0,   0,   0,   0,  -3,  -3,   0,  -2,  -3,   0,   4,   2,  -4,  -2,  -2,   1,   0,  -1,  -6,
              1,  -1,   0,  -1,   0,   0,  -1,  -2,  -1,  -1,  -2,   0,  -1,  -3,  -1,   2,   4,  -5,  -3,   0,   0,  -1,  -1,  -6,
             -4,  -3,  -5,  -6,  -6,  -6,  -6,  -5,  -1,  -2,  -1,  -4,  -3,   3,  -4,  -4,  -5,  15,   4,  -4,  -5,  -6,  -2,  -6,
             -3,  -2,  -2,  -4,   0,  -3,  -3,  -5,   3,  -2,  -1,  -3,  -2,   5,  -5,  -2,  -3,   4,   9,  -2,  -3,  -3,  -1,  -6,
              0,  -3,  -3,  -4,   1,  -2,  -3,  -4,  -3,   4,   2,  -3,   2,  -1,  -3,  -2,   0,  -4,  -2,   4,  -3,  -3,  -1,  -6,
             -1,  -1,   4,   5,  -3,   1,   2,   0,   1,  -4,  -5,   0,  -3,  -5,  -1,   1,   0,  -5,  -3,  -3,   4,   1,  -1,  -6,
             -1,   0,   1,   2,  -4,   4,   4,  -2,   1,  -3,  -3,   2,  -2,  -4,  -1,   0,  -1,  -6,  -3,  -3,   1,   4,  -1,  -6,
             -1,  -1,  -1,  -2,  -1,  -1,  -1,  -2,   0,  -1,  -1,  -1,  -1,  -1,  -2,  -1,  -1,  -2,  -1,  -1,  -1,  -1,  -1,  -6,
             -6,  -6,  -6,  -6,  -6,  -6,  -6,  -6,  -6,  -6,  -6,  -6,  -6,  -6,  -6,  -6,  -6,  -6,  -6,  -6,  -6,  -6,  -6,   1,
        };
        return _data;
    }
};

}  // namespace SEQAN_NAMESPACE_MAIN

#endif  // SEQAN_SCORE_SCORE_MATRIX_DATA_H_
