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

#ifndef SEQAN_HEADER_PUMP_EXTENDER3_H
#define SEQAN_HEADER_PUMP_EXTENDER3_H

namespace SEQAN_NAMESPACE_MAIN
{

//namespace SEQAN_NAMESPACE_PIPELINING
//{

    struct Extender3;

    template < typename TTextInput, typename TNameInput >
    struct Pipe< Bundle2< TTextInput, TNameInput >, Extender3 >
    {
        enum { maxShift = 2 };
        typedef typename Size<Pipe>::Type                   TSize;
        typedef typename Value<TTextInput>::Type            TValue;
        typedef typename Value<TNameInput>::Type            TNameInputValue;
        typedef typename Value<TNameInputValue, 2>::Type    TName;

        typedef Tuple<TValue, maxShift>             XTuple;
        typedef Tuple<TName, maxShift>              NTuple;
        typedef Triple<TSize, NTuple, XTuple, Pack> OutType0;
        typedef Triple<TSize, NTuple, XTuple, Pack> OutType12;

        // pipeline interfaces to ease specialization
        typedef Pipe< void, AbstractSource< OutType0, TSize > > Out0;
        typedef Pipe< void, AbstractSource< OutType12, TSize > > Out12;
    };

    template < typename TTextInput, typename TNameInput, typename TOut0, typename TOut12 >
    static bool _skew3Extend(TTextInput &textIn, TNameInput &nameIn, TOut0 &out0, TOut12 &out12)
    {
        resize(out0, length(textIn) / 3);
        resize(out12, length(nameIn));
        if (!(
            beginRead(textIn) &&
            beginRead(nameIn) &&
            beginWrite(out0) &&
            beginWrite(out12))) return false;

        typename Value<TOut0>::Type  o0 = typename Value<TOut0>::Type();
        typename Value<TOut12>::Type o1 = typename Value<TOut12>::Type();
        typename Value<TOut12>::Type o2 = typename Value<TOut12>::Type();

        unsigned r = (unsigned)(length(textIn) % 3);
        bool filled = (r != 0);

        if (r == 2)
        {
            o2.i1 = (*nameIn).i1;
            o2.i2[0] = (*nameIn).i2; ++nameIn;
            o2.i3[0] = *textIn; ++textIn;
        }

        if (r >= 1)
        {
            o1.i1 = (*nameIn).i1;
            o1.i2[0] = (*nameIn).i2; ++nameIn;
            o1.i3[0] = *textIn; ++textIn;
        }

        if (r == 2)
        {
            o2.i2[1] = o1.i2[0];
            o2.i3[1] = o1.i3[0];
            push(out12, o2);
        }

        while (!eof(nameIn))
        {
            o1.i3[1] = o0.i3[0] = *textIn; ++textIn;
            o2.i3[0] = o0.i3[1] = *textIn; ++textIn;

            o2.i1 = (*nameIn).i1;
            o0.i2[0] = o2.i2[0] = o1.i2[1] = (*nameIn).i2; ++nameIn;
            o2.i3[1] = *textIn; ++textIn;

            if (filled)
                push(out12, o1);
            else
                filled = true;

            o1.i1 = (*nameIn).i1;
            o0.i2[1] = o2.i2[1] = o1.i2[0] = (*nameIn).i2; ++nameIn;
            o1.i3[0] = o2.i3[1];
            push(out12, o2);
            o0.i1 = o2.i1 + 1;
            push(out0, o0);
        }

        o1.i2[1] = 0;
        o1.i3[1] = 0;
        if (filled) push(out12, o1);

        endWrite(out12);
        endWrite(out0);
        endRead(nameIn);
        endRead(textIn);
        return true;
    }

//}

}

#endif
