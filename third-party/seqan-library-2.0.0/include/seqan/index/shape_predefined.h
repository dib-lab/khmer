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
// Author: David Weese <david.weese@fu-berlin.de>
// ==========================================================================

#ifndef SEQAN_HEADER_SHAPE_PREDEFINED_H
#define SEQAN_HEADER_SHAPE_PREDEFINED_H

//SEQAN_NO_GENERATED_FORWARDS: no forwards are generated for this file

namespace SEQAN_NAMESPACE_MAIN
{

    //////////////////////////////////////////////////////////////////////////////
    // some predefined gapped shapes


    //////////////////////////////////////////////////////////////////////////////
    // Single seed of
    // B.Ma and J.Tromp and M.Li,
    // "PatternHunter: faster and more sensitive homology search"
    // Bioinformatics 18, 2002
    //
    // weight:11
    // length:18
    //
    // shape:
    // 111010010100110111

    typedef GappedShape<
        HardwiredShape< 1, 1, 2, 3, 2, 3, 1, 2, 1, 1 >
    > ShapePatternHunter;



    //////////////////////////////////////////////////////////////////////////////
    // Multiple seeds of
    // L.Ilie and S.Ilie, "Fast Computation of Good Multiple Spaced Seeds"
    // WABI, 2007
    //
    // weight:9
    // length:15
    //
    // shapes:
    // 111010100100111
    // 110100110011101
    // 111010001011011
    //
    // sensitivity:
    // 65% 0.747975        70% 0.897741
    // 75% 0.973134        80% 0.996226

    typedef GappedShape<
        HardwiredShape< 1, 1, 2, 2, 3, 3, 1, 1 >
    > ShapeIlieA1;

    typedef GappedShape<
        HardwiredShape< 1, 2, 3, 1, 3, 1, 1, 2 >
    > ShapeIlieA2;

    typedef GappedShape<
        HardwiredShape< 1, 1, 2, 4, 2, 1, 2, 1 >
    > ShapeIlieA3;



    //////////////////////////////////////////////////////////////////////////////
    // Multiple seeds of
    // L.Ilie and S.Ilie, "Fast Computation of Good Multiple Spaced Seeds"
    // WABI 2007
    //
    // weight:9
    // length:13..23
    //
    // shapes:
    // 1110110100111
    // 11010000110010111
    // 11100010010000101011
    //
    // sensitivity:
    // 65% 0.767413        70% 0.910949
    // 75% 0.978558        80% 0.997357

    typedef GappedShape<
        HardwiredShape< 1, 1, 2, 1, 2, 3, 1, 1 >
    > ShapeIlieB1;

    typedef GappedShape<
        HardwiredShape< 1, 2, 5, 1, 3, 2, 1, 1 >
    > ShapeIlieB2;

    typedef GappedShape<
        HardwiredShape< 1, 1, 4, 3, 5, 2, 2, 1 >
    > ShapeIlieB3;


}    // namespace seqan

#endif
