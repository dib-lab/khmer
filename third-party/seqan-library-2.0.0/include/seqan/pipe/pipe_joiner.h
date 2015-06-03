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

#ifndef SEQAN_HEADER_PIPE_JOINER_H
#define SEQAN_HEADER_PIPE_JOINER_H

namespace SEQAN_NAMESPACE_MAIN
{

//namespace SEQAN_NAMESPACE_PIPELINING
//{

    struct Joiner;

    template < typename TInput1, typename TInput2 >
    struct Value< Pipe< Bundle2< TInput1, TInput2 >, Joiner > > {
        typedef Pair<
            typename Value<TInput1>::Type,
            typename Value<TInput2>::Type
        > Type;
    };

    template < typename TInput1, typename TInput2, typename TInput3 >
    struct Value< Pipe< Bundle3< TInput1, TInput2, TInput3 >, Joiner > > {
        typedef Triple<
            typename Value<TInput1>::Type,
            typename Value<TInput2>::Type,
            typename Value<TInput3>::Type
        > Type;
    };

/*!
 * @class Joiner
 * @extends Pipe
 * @headerfile <seqan/pipe.h>
 * @brief Joins two or three input streams.
 *
 * @signature template <typename TInput1, typename TInput2>
 *            class Pipe<Bundle2<TInput1, TInput2>, Joiner>;
 * @signature template <typename TInput1, typename TInput2, typename TInput3>
 *            class Pipe<Bundle3<TInput1, TInput2, TInput3>, Joiner>;
 *
 * @tparam TInput1 The type of the first pipeline module this module reads from.
 * @tparam TInput2 The type of the second pipeline module this module reads from.
 * @tparam TInput3 The type of the third pipeline module this module reads from.
 *
 * The output type is a packed @link  Pair  @endlink or @link Triple @endlink of the input types <tt>Value&lt;TInputX&gt;::Type</tt>.
 */

    //////////////////////////////////////////////////////////////////////////////
    // joiner class
    template < typename TInput1, typename TInput2 >
    struct Pipe< Bundle2< TInput1, TInput2 >, Joiner >
    {
        Bundle2< TInput1, TInput2 >    in;
        typename Value<Pipe>::Type    tmp;

        Pipe(Bundle2< TInput1, TInput2 > _in):
            in(_in) {}

        inline typename Value<Pipe>::Type const & operator*() {
            tmp.i1 = *in.in1;
            tmp.i2 = *in.in2;
            return tmp;
        }

        inline Pipe& operator++() {
            ++in.in1;
            ++in.in2;
            return *this;
        }
    };

    template < typename TInput1, typename TInput2, typename TInput3 >
    struct Pipe< Bundle3< TInput1, TInput2, TInput3 >, Joiner >
    {
        Bundle3< TInput1, TInput2, TInput3 >    in;
        typename Value<Pipe>::Type                tmp;

        Pipe(Bundle3< TInput1, TInput2, TInput3 > _in):
            in(_in) {}

        inline typename Value<Pipe>::Type const & operator*() {
            tmp.i1 = *in.in1;
            tmp.i2 = *in.in2;
            tmp.i3 = *in.in3;
            return tmp;
        }

        inline Pipe& operator++() {
            ++in.in1;
            ++in.in2;
            ++in.in3;
            return *this;
        }
    };

//}

}

#endif
