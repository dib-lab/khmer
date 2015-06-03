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

#ifndef SEQAN_HEADER_PIPE_COUNTER_H
#define SEQAN_HEADER_PIPE_COUNTER_H

namespace SEQAN_NAMESPACE_MAIN
{

//namespace SEQAN_NAMESPACE_PIPELINING
//{

    struct Counter;

    template < typename TInput >
    struct Value< Pipe< TInput, Counter > > {
        typedef Pair<
            typename Value<TInput>::Type,
            typename Size<TInput>::Type,
            Pack
        > Type;
    };

/*!
 * @class Counter
 * @extends Pipe
 * @headerfile <seqan/pipe.h>
 * @brief Extends the input stream by a second field which enumerates the elements.
 *
 * @signature template <typename TInput>
 *            class Pipe<TInput, Counter>;
 *
 * @tparam TInput The type of the pipeline module this module reads from.
 *
 * The output type is a Pair of input type and size type (i.e.  <tt>Pair&lt;Value&lt;TInput&gt;::Type, Size&lt;TInput&gt;::Type</tt>).
 *
 * The first output field is the original input stream.
 *
 * The second output field begins with 0 and increases by 1 per element.
 */

    //////////////////////////////////////////////////////////////////////////////
    // counter class
    template < typename TInput >
    struct Pipe< TInput, Counter >
    {
        TInput                      &in;
        typename Value<Pipe>::Type    tmp;

        Pipe(TInput& _in):
            in(_in) {}

        inline typename Value<Pipe>::Type const & operator*() {
            tmp.i1 = *in;
            return tmp;
        }

        inline Pipe& operator++() {
            ++in;
            ++tmp.i2;
            return *this;
        }
    };


    //////////////////////////////////////////////////////////////////////////////
    // global pipe functions
    template < typename TInput >
    inline bool control(Pipe< TInput, Counter > &me, ControlBeginRead const &command) {
        if (!control(me.in, command)) return false;
        me.tmp.i2 = 0;
        return true;
    }

//}

}

#endif
