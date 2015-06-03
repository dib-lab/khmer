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

#ifndef SEQAN_HEADER_PIPE_SHIFTER_H
#define SEQAN_HEADER_PIPE_SHIFTER_H

namespace SEQAN_NAMESPACE_MAIN
{

//namespace SEQAN_NAMESPACE_PIPELINING
//{

//    template < int delta, bool omitBlank = false, bool _echoing = (delta < 0) >
//    struct Shifter;
    template < int delta, bool omitBlank = false, bool _echoing = true >
    struct Shifter;

/*!
 * @class Shifter
 *
 * @extends Pipe
 *
 * @headerfile <seqan/pipe.h>
 *
 * @brief Shifts the input stream by <tt>delta</tt> elements.
 *
 * @signature template <typename TInput, int DELTA[, bool OMIT_BLANK]>
 *            class Pipe<TInput, Shifter<DELTA, OMIT_BLANK> >;
 *
 * @tparam TInput     The type of the pipeline module this module reads from.
 * @tparam OMIT_BLANK Omit undefined entries.  If <tt>true</tt>, the output stream is <tt>|delta|</tt> elements shorter
 *                    than the input stream.  If <tt>false</tt>, the lengths are equal and blanks (default constructed
 *                    elements) are inserted on the cut-off-opposite side, defaults to false.
 * @tparam DELTA      The shift size. For the output stream holds <tt>out[i]=in[i+delta]</tt>.  For <tt>delta>0</tt>
 *                    the input stream is cut of at the beginning and for <tt>delta<0</tt> at the end.
 *
 * The output type equals the input type.
 */

    //////////////////////////////////////////////////////////////////////////////
    // echoer class
    template < typename TInput, int delta, bool omitBlank >
    struct Pipe< TInput, Shifter<delta, omitBlank, true> >
    {
        TInput                      &in;
        typename Size<Pipe>::Type    blankCounter, charCounter;
        typename Value<Pipe>::Type    blank;

        Pipe(TInput& _in):
            in(_in),
            blank()    {}

        inline typename Value<Pipe>::Type const & operator*() {
            if (blankCounter)    return blank;
            else                return *in;
        }

        inline Pipe& operator++() {
            if (blankCounter)
                --blankCounter;
            else {
                ++in;
                --charCounter;
            }
            return *this;
        }
    };


    template < typename TInput, int delta, bool omitBlank >
    struct Pipe< TInput, Shifter<delta, omitBlank, false> >
    {
        TInput                      &in;
        typename Size<Pipe>::Type    blankCounter, charCounter;
        typename Value<Pipe>::Type    blank;

        Pipe(TInput& _in):
            in(_in),
            blank()    {}

        inline typename Value<Pipe>::Type const & operator*() {
            if (charCounter)    return *in;
            else                return blank;
        }

        inline Pipe& operator++() {
            if (charCounter) {
                ++in;
                --charCounter;
            } else
                --blankCounter;
            return *this;
        }
    };


    //////////////////////////////////////////////////////////////////////////////
    // global pipe functions
    template < typename TInput, int delta, bool omitBlank >
    inline bool control(Pipe< TInput, Shifter< delta, omitBlank, false > > &me, ControlBeginRead const &command) {
        if (!control(me.in, command)) return false;
        for(typename Size<TInput>::Type i = 0; i < delta && !eof(me.in); ++i)
            ++me.in;
        me.blankCounter = (omitBlank)? 0: delta;
        me.charCounter = length(me.in) - delta;
        return !eof(me.in);
    }

    template < typename TInput, int delta, bool omitBlank >
    inline bool control(Pipe< TInput, Shifter< delta, omitBlank, true > > &me, ControlBeginRead const &command) {
        if (!control(me.in, command)) return false;
        me.blankCounter = (omitBlank)? 0: -delta;
        me.charCounter = length(me.in) + delta;
        return true;
    }




    template < typename TInput, int delta, bool _echoing >
    inline Size< Pipe< TInput, Shifter< delta, true, _echoing > > >
    length(Pipe< TInput, Shifter< delta, true, _echoing > > const &me) {
        return length(me.in) - abs(delta);
    }

    template < typename TInput, int delta, bool omitBlank, bool _echoing >
    inline bool control(Pipe< TInput, Shifter< delta, omitBlank, _echoing > > &me, ControlEof const &/*command*/) {
        return me.charCounter == 0 && me.blankCounter == 0;
    }

    template < typename TInput, int delta, bool omitBlank, bool _echoing >
    inline bool control(Pipe< TInput, Shifter< delta, omitBlank, _echoing > > &me, ControlEos const &/*command*/) {
        return control(me, ControlEof());
    }

//}

}

#endif
