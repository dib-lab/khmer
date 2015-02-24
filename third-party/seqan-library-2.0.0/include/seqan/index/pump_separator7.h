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

#ifndef SEQAN_HEADER_INDEX_PUMP_SEPARATOR_H
#define SEQAN_HEADER_INDEX_PUMP_SEPARATOR_H

namespace SEQAN_NAMESPACE_MAIN
{


    template <
        typename TInput, typename TFunctor,
        typename TOut1, typename TOut2, typename TOut4
    >
    static void _skew7SeparateSlices(
        TInput &in, TFunctor const &funcSlice,
        TOut1 &out1, TOut2 &out2, TOut4 &out4)
    {
        beginRead(in);

        resize(out1, funcSlice.n1);
        resize(out2, funcSlice.n2);
        resize(out4, funcSlice.n4);

        beginWrite(out1);
        beginWrite(out2);
        beginWrite(out4);

        typename Value<TInput>::Type i;
        while (!eof(in)) {
            pop(in, i);
            if (i.i1 < funcSlice.n4) {
                push(out4, i);
            } else
                if (i.i1 < funcSlice.n24) {
                    i.i1 -= funcSlice.n4;
                    push(out2, i);
                } else {
                    i.i1 -= funcSlice.n24;
                    push(out1, i);
                }
        }

        endWrite(out4);
        endWrite(out2);
        endWrite(out1);
        endRead(in);
    }

//}

}

#endif
