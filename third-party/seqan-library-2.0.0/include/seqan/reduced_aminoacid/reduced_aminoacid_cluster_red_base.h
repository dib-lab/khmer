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

#ifndef SEQAN_REDUCED_AMINOACID_CLUSTER_RED_BASE_H_
#define SEQAN_REDUCED_AMINOACID_CLUSTER_RED_BASE_H_

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// -----------------------------------------------------------------------
// Tag ClusterReduction
// -----------------------------------------------------------------------

/*
 * @class ClusterReduction
 * @brief Specialization for @link ReducedAminoAcid @endlink
 * @headerfile seqan/reduced_aminoacid.h
 *
 * @signature template <unsigned char n, unsigned char m = 24, typename TMatrix = Blosum62>
 * struct ClusterReduction;
 *
 * @tparam n the size of the reduced alphabet (between 2 and m-1)
 * @tparam m size to truncate alphabet to, <b>before</b> clustering
 * (one of 20, 22, 24; default 24)
 * @tparam TMatrix Matrix used for clustering (default @link Blosum62 @endlink,
 * none other supported right now)
 *
 * @section WhenToUse
 *
 * Use m = 24 when you expect 'X' and '*' in the dataset you reduce
 * from. This is especially the case on translated genomic reads.
 *
 * If you have validated protein sequences, you can use can use m = 20 or
 * m = 22, which will not include special characters (see
 * @link AminoAcid @endlink for details).
 *
 * @section Background
 *
 * The method employed for reducing the alphabet is similar to Murphy et al,
 * 2000, <a href="http://www.ncbi.nlm.nih.gov/pubmed/10775656">http://www.ncbi.nlm.nih.gov/pubmed/10775656</a>
 *
 * Correlation coefficients for the Blosum62 scores of all pairs of amino
 * acids in the alphabet were computed and clustered with WPGMA (using
 * UPGMA as second criterium when WPGMA yields the same
 * distance between two clusters).
 *
 * <img src="ClusterReduction.png">
 *
 * The exact clustering for m = 24.
 *
 */

template <unsigned char n, unsigned char m = 24, typename TMatrix = Blosum62>
struct ClusterReduction;

// ============================================================================
// Metafunctions
// ============================================================================

// -----------------------------------------------------------------------
// Metafunction ValueSize
// -----------------------------------------------------------------------

template <unsigned char n, unsigned char m, typename TMatrix>
struct ValueSize<
        SimpleType<unsigned char,
                   ReducedAminoAcid_<ClusterReduction<n, m, TMatrix> > > >
{
    typedef uint8_t Type;
    static const Type VALUE = n;
};

// -----------------------------------------------------------------------
// Metafunction BitPerValue
// -----------------------------------------------------------------------

template <unsigned char n, unsigned char m, typename TMatrix>
struct BitsPerValue<
         SimpleType<unsigned char,
                    ReducedAminoAcid_<ClusterReduction<n, m, TMatrix> > > >
{
    typedef uint8_t Type;
    static const Type VALUE = Log2<n>::VALUE;
};

// -----------------------------------------------------------------------
// Translation Tables (implementations see extra files)
// -----------------------------------------------------------------------

// ============================================================================
// Functions
// ============================================================================

}
#endif // def SEQAN_REDUCED_AMINOACID_CLUSTER_RED_BASE_H_
