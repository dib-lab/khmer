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

#ifndef SEQAN_HEADER_PIPE_ECHOER_H
#define SEQAN_HEADER_PIPE_ECHOER_H

namespace SEQAN_NAMESPACE_MAIN
{

//namespace SEQAN_NAMESPACE_PIPELINING
//{

    //////////////////////////////////////////////////////////////////////////////
    // some metaprogramming to unrool fixed-size loops
    struct EchoerFillWorker_ {
        template <typename Arg>
        static inline void body(Arg &arg, unsigned I) {
            arg.tmp.i2[I-1] = *(arg.in); ++(arg.in);
        }
    };

    struct EchoerClearWorker_ {
        template <typename Arg>
        static inline void body(Arg &arg, unsigned I) {
            arg.i2[I] = typename Value< typename Value<Arg, 2>::Type >::Type ();
        }
    };

    struct EchoerShiftWorker_ {
        template <typename Arg>
        static inline void body(Arg &arg, unsigned I) {
            arg.i2[I] = arg.i2[I-1];
        }
    };


    template < unsigned echoRepeats, bool omitFirst >
    struct Echoer;

    template < typename TInput, unsigned echoRepeats, bool omitFirst >
    struct Value< Pipe< TInput, Echoer< echoRepeats, omitFirst > > > {
        typedef Tuple<typename Value<TInput>::Type, echoRepeats>    EchoType;
        typedef Pair<typename Size<TInput>::Type, EchoType>            Type;
    };

/*!
 * @class Echoer
 * @extends Pipe
 * @headerfile <seqan/pipe.h>
 *
 * @brief Outputs tuples of the <tt>echoRepeats</tt> last elements of the input stream.
 *
 * @signature template <typename Input, unsigned ECHO_REPEATS, bool OMIT_FIRST>
 *            class Pipe;
 *
 * @tparam TInput       The type of the pipeline module this module reads from.
 * @tparam ECHO_REPEATS The tuple length.The tuples contain elements <tt>in[i]in[i-1]...in[i-(echoRepeats-1)]</tt>.
 * @tparam OMIT_FIRST   Omit half filled tuples.  If <tt>true</tt>, the output stream is <tt>echoRepeats-1</tt>
 *                      elements shorter than the input stream.  If <tt>false</tt>, the lengths are identical and the
 *                      tuple is filled with blanks (default constructed elements) for undefined entries.
 *
 * The output type is a @link Tuple @endlink of input elements and length <tt>echoRepeats</tt> (i.e.
 * <tt>Tuple&lt;Value&lt;TInput&lt;::Type, echoRepeats&gt;</tt>).
 *
 * The tuples are sequences of the form <tt>in[i]in[i-1]in[i-2]..in[i-echoRepeats+1]</tt>. For <tt>omitFirst=false</tt>
 * <tt>i</tt> begins with 0 and for <tt>omitFirst=true</tt> <tt>i</tt> begins with <tt>echoRepeats-1</tt>.
 */

    //////////////////////////////////////////////////////////////////////////////
    // echoer class
    template < typename TInput, unsigned echoRepeats, bool omitFirst >
    struct Pipe< TInput, Echoer<echoRepeats, omitFirst> >
    {
        typedef typename Value<Pipe>::Type TValue;

        TInput    &in;
        TValue    tmp;

        Pipe(TInput& _in):
            in(_in),
            tmp(0, typename Value<TValue, 2>::Type()) {}

        inline typename Value<Pipe>::Type const & operator*() const {
            return tmp;
        }

        inline Pipe& operator++() {
            ++in;
            if (eof(in)) return *this;
            LoopReverse<EchoerShiftWorker_, echoRepeats - 1>::run(this->tmp);
            ++tmp.i1;
            tmp.i2[0] = *in;
            return *this;
        }
    };


    //////////////////////////////////////////////////////////////////////////////
    // global pipe functions
    template < typename TInput, unsigned echoRepeats, bool omitFirst >
    inline bool control(Pipe< TInput, Echoer< echoRepeats, omitFirst > > &me, ControlBeginRead const &command) {
        if (!control(me.in, command)) return false;
        me.tmp.i1 = 0;
        Loop<EchoerClearWorker_, echoRepeats - 1>::run(me.tmp);
        if (!eof(me.in)) me.tmp.i2[0] = *me.in;
        return true;
    }

    template < typename TInput, unsigned echoRepeats >
    inline bool control(Pipe< TInput, Echoer< echoRepeats, true > > &me, ControlBeginRead const &command) {
        if (!control(me.in, command) || length(me.in) < echoRepeats - 1) return false;
        me.tmp.i1 = 0;
        LoopReverse<EchoerFillWorker_, echoRepeats - 1>::run(me);
        if (!eof(me.in)) me.tmp.i2[0] = *me.in;
        return true;
    }

    template < typename TInput, unsigned echoRepeats >
    inline Size< Pipe< TInput, Echoer< echoRepeats, true > > >
    length(Pipe< TInput, Echoer< echoRepeats, true > > const &me) {
        return length(me.in) - (echoRepeats - 1);
    }

//}

}

#endif
