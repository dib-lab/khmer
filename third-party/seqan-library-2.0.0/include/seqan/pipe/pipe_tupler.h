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

#ifndef SEQAN_HEADER_PIPE_TUPLER_H
#define SEQAN_HEADER_PIPE_TUPLER_H

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////

    template < unsigned SIZE, bool omitLast = false, typename TPack = void >
    struct Tupler;

    template < typename TInput, unsigned SIZE, bool omitLast, typename TPack >
    struct Value< Pipe< TInput, Tupler<SIZE, omitLast, TPack> > >
    {
        typedef Tuple<typename Value<TInput>::Type, SIZE, TPack>    TTuple;
        typedef Pair<typename Size<TInput>::Type, TTuple, Pack>     Type;
    };

//////////////////////////////////////////////////////////////////////////////

    template <
        typename TInput,
        unsigned SIZE,
        bool omitLast,
        typename TPack,
        typename TPair,
        typename TLimitsString >
    struct Value< Pipe< TInput, Multi< Tupler<SIZE, omitLast, TPack>, TPair, TLimitsString > > >
    {
        typedef Tuple<typename Value<TInput>::Type, SIZE, TPack>    TTuple;
        typedef Pair<TPair, TTuple, Pack>                           Type;
    };

//////////////////////////////////////////////////////////////////////////////


    // output only fully filled tuples
    template <unsigned SIZE, bool omitLast>
    struct TuplerNumberOfLastTuples_
    {
        enum { VALUE = 1 };
    };

    // output SIZE-1 half filled tuples at the end
    template <unsigned SIZE>
    struct TuplerNumberOfLastTuples_<SIZE, false>
    {
        enum { VALUE = SIZE };
    };

    struct ShiftLeftWorker_ {
        template <typename Arg>
        static inline void body(Arg &arg, unsigned I) {
            arg[I-1] = arg[I];
        }
    };


    // shift left by 1 character
    template <typename TValue, unsigned SIZE, typename TSpec>
    inline void
    _tuplerShiftLeft(Tuple<TValue, SIZE, TSpec> & tuple)
    {
        Loop<ShiftLeftWorker_, SIZE - 1>::run(tuple);
    }
    template <typename TValue, unsigned SIZE>
    inline void
    _tuplerShiftLeft(Tuple<TValue, SIZE, BitPacked<> > & tuple)
    {
        tuple <<= 1;
    }

    // assign last character
    template <typename TValue, unsigned SIZE, typename TSpec>
    inline void
    _tuplerAssignLast(Tuple<TValue, SIZE, TSpec> & tuple, TValue const & val)
    {
        tuple[SIZE - 1] = val;
    }
    template <typename TValue, unsigned SIZE>
    inline void
    _tuplerAssignLast(Tuple<TValue, SIZE, BitPacked<> > & tuple, TValue const & val)
    {
        tuple |= val;
    }



/*!
 * @class Tupler
 * @extends Pipe
 * @headerfile <seqan/pipe.h>
 *
 * @brief Outputs tuples of the <tt>SIZE</tt> consecutive elements of the input stream.
 *
 * @signature template <typename TInput, unsigned TUPLE_LEN, bool OMIT_LAST>
 *            class Pipe<TInput, Tupler<TUPLE_LEN, OMIT_LAST> >;
 *
 * @tparam TInput    The type of the pipeline module this module reads from.
 * @tparam TUPLE_LEN The tuple length.The tuples contain elements <tt>in[i]in[i+1]...in[i+(SIZE-1)]</tt>.
 * @tparam OMIT_LAST Omit half filled tuples.  If <tt>true</tt>, the output stream is <tt>SIZE-1</tt> elements
 *                   shorter than the input stream.  If <tt>false</tt>, the lengths are identical and the last tuples
 *                   are filled with blanks (default constructed elements) for undefined entries.
 *
 * The output type is a @link Tuple @endlink of input elements and length
 * <tt>SIZE</tt> (i.e. <tt>Tuple&lt;Value&lt;TInput&gt;::Type, TUPLE_LEN&gt;</tt>).
 *
 * The tuples are sequences of the form
 * <tt>in[i]in[i-1]in[i-2]..in[i-SIZE+1]</tt>. For <tt>omitLast=false</tt>
 * <tt>i</tt> begins with 0 and for <tt>omitLast=true</tt> <tt>i</tt> begins
 * with <tt>SIZE-1</tt>.
 */

    template <typename TInput, unsigned SIZE, bool omitLast, typename TPack>
    struct Pipe<TInput, Tupler<SIZE, omitLast, TPack> >
    {
        typedef typename Value<TInput>::Type TValue;

        TInput                      &in;
        typename Value<Pipe>::Type    tmp;
        typename Size<TInput>::Type    lastTuples;

        Pipe(TInput& _in):
            in(_in) {}

        inline typename Value<Pipe>::Type const & operator*() const
        {
            return tmp;
        }

        inline Pipe& operator++()
        {
            if (eof(in)) --lastTuples;
            _tuplerShiftLeft(tmp.i2);
            ++tmp.i1;
            if (lastTuples < TuplerNumberOfLastTuples_<SIZE, omitLast>::VALUE)
                _tuplerAssignLast(tmp.i2, TValue());
            else {
                _tuplerAssignLast(tmp.i2, *in);
                ++in;
            }
            return *this;
        }


        template <typename TPack_>
        inline unsigned _tryFill(TPack_ *)
        {
            unsigned i;
            for (i = 0; i < SIZE && !eof(in); ++i, ++in)
                tmp.i2.i[i] = *in;
            for (unsigned j = i; j < SIZE; ++j)
                tmp.i2.i[j] = TValue();
            return i;
        }

        inline unsigned _tryFill(BitPacked<> *)
        {
            unsigned i;
            clear(tmp.i2);
            for (i = 0; i < SIZE && !eof(in); ++i, ++in)
            {
                tmp.i2 <<= 1;
                tmp.i2 |= *in;
            }
            tmp.i2 <<= (SIZE - i);
            return i;
        }

        inline void fill()
        {
            unsigned charsRead = _tryFill(static_cast<TPack*>(NULL));
            if (TuplerNumberOfLastTuples_<SIZE, omitLast>::VALUE > SIZE - charsRead)
                lastTuples = TuplerNumberOfLastTuples_<SIZE, omitLast>::VALUE - (SIZE - charsRead);
            else
                lastTuples = 0;
            tmp.i1 = 0;
        }
    };

//____________________________________________________________________________


    template <
        typename TInput,
        unsigned SIZE,
        bool omitLast,
        typename TPack,
        typename TPair,
        typename TLimitsString >
    struct Pipe<TInput, Multi<Tupler<SIZE, omitLast, TPack>, TPair, TLimitsString> >
    {
        typedef typename Value<typename Value<Pipe>::Type, 2>::Type TTuple;
        typedef typename Value<TTuple>::Type                        TValue;

        typedef PairIncrementer_<TPair, TLimitsString>              Incrementer;

        TInput                      &in;
        Incrementer                    localPos;
        typename Value<Pipe>::Type    tmp;
        typename Size<TInput>::Type    seqLength, lastTuples;

        TLimitsString const &limits;

        template <typename TLimitsString_>
        Pipe(TInput& _in, TLimitsString_ &_limits):  // const &_limits is intentionally omitted to suppress implicit casts (if types mismatch) and taking refs of them
            in(_in),
            tmp(),
            seqLength(),
            lastTuples(),
            limits(_limits)
        {}

        inline typename Value<Pipe>::Type const & operator*() const
        {
            return tmp;
        }

        inline Pipe& operator++()
        {
            // process next sequence
            if (eos())
                if (--lastTuples == 0)
                {
                    assignValueI1(tmp.i1, getValueI1(tmp.i1) + 1);
                    fill();
                    return *this;
                }

            // shift left by 1 character
            _tuplerShiftLeft(tmp.i2);
            assignValueI2(tmp.i1, getValueI2(tmp.i1) + 1);
            if (lastTuples < TuplerNumberOfLastTuples_<SIZE, omitLast>::VALUE)
            {
                _tuplerAssignLast(tmp.i2, TValue());
            }
            else
            {
                _tuplerAssignLast(tmp.i2, *in);
                ++in;
                ++localPos;
            }
            return *this;
        }

        template <typename TPack_>
        inline unsigned _tryFill(TPack_ *)
        {
            unsigned i = 0;
            if (!eof(in))
            {
                do {
                    tmp.i2.i[i] = *in;
                    ++in;
                    ++i;
                    ++localPos;
                } while ((i < SIZE) && !eos());
            }

            // fill up with null chars
            for (unsigned j = i; j < SIZE; ++j)
                tmp.i2.i[j] = TValue();
            return i;
        }

        inline unsigned _tryFill(BitPacked<> *)
        {
            unsigned i = 0;
            if (!eof(in))
            {
                do
                {
                    tmp.i2 <<= 1;
                    tmp.i2 |= *in;
                    ++in;
                    ++i;
                    ++localPos;
                } while ((i < SIZE) && !eos());
            }

            // fill up with null chars
            tmp.i2 <<= (SIZE - i);
            return i;
        }

        inline void fill()
        {
            assignValueI2(tmp.i1, 0);
            do
            {
                assignValueI1(tmp.i1, getValueI1(value(localPos)));
                unsigned charsRead = _tryFill(static_cast<TPack*>(NULL));
                lastTuples = TuplerNumberOfLastTuples_<SIZE, omitLast>::VALUE;

                // eventually, reduce the number of half-filled tuples
                if (lastTuples <= SIZE - charsRead)
                    lastTuples = 0;
                else
                    lastTuples -= SIZE - charsRead;

            } while ((lastTuples == 0) && !eof(in));
        }

        inline bool eos()
        {
            return (getValueI1(value(localPos)) > 0) && (getValueI2(value(localPos)) == 0);
        }
    };


    //////////////////////////////////////////////////////////////////////////////
    // global pipe functions
    template < typename TInput, unsigned SIZE, bool omitLast, typename TPack >
    inline bool
    control(
        Pipe< TInput, Tupler<SIZE, omitLast, TPack> > &me,
        ControlBeginRead const &command)
    {
        if (!control(me.in, command)) return false;
        me.fill();
        return true;
    }

    template <
        typename TInput,
        unsigned SIZE,
        bool omitLast,
        typename TPack,
        typename TPair,
        typename TLimitsString >
    inline bool
    control(
        Pipe< TInput, Multi<Tupler<SIZE, omitLast, TPack>, TPair, TLimitsString> > &me,
        ControlBeginRead const &command)
    {
        if (!control(me.in, command)) return false;
        setHost(me.localPos, me.limits);
        assignValueI1(me.tmp.i1, 0);
        me.fill();
        return true;
    }

    template < typename TInput, unsigned SIZE, bool omitLast, typename TPack >
    inline bool
    control(
        Pipe< TInput, Tupler<SIZE, omitLast, TPack> > &me,
        ControlEof const &)
    {
        return me.lastTuples == 0;
    }

    template <
        typename TInput,
        unsigned SIZE,
        bool omitLast,
        typename TPack,
        typename TPair,
        typename TLimitsString >
    inline bool
    control(
        Pipe< TInput, Multi<Tupler<SIZE, omitLast, TPack>, TPair, TLimitsString> > &me,
        ControlEof const &)
    {
        return me.lastTuples == 0;
    }

    template < typename TInput, unsigned SIZE, bool omitLast, typename TPack >
    inline bool
    control(
        Pipe< TInput, Tupler<SIZE, omitLast, TPack> > &me,
        ControlEos const &)
    {
        return control(me, ControlEof());
    }

    template <
        typename TInput,
        unsigned SIZE,
        bool omitLast,
        typename TPack,
        typename TPair,
        typename TLimitsString >
    inline bool
    control(
        Pipe< TInput, Multi<Tupler<SIZE, omitLast, TPack>, TPair, TLimitsString> > &me,
        ControlEos const &)
    {
        return (getValueI1(me.tmp.i1) > 0) && (getValueI2(me.tmp.i1) == 0);
    }

    template < typename TInput, unsigned SIZE, bool omitLast, typename TPack >
    inline typename Size< Pipe< TInput, Tupler<SIZE, omitLast, TPack> > >::Type
    length(Pipe< TInput, Tupler<SIZE, omitLast, TPack> > const &me)
    {
        if (length(me.in) > (SIZE - TuplerNumberOfLastTuples_<SIZE, omitLast>::VALUE))
            return length(me.in) - (SIZE - TuplerNumberOfLastTuples_<SIZE, omitLast>::VALUE);
        else
            return 0;
    }

    template <
        typename TInput,
        unsigned SIZE,
        bool omitLast,
        typename TPack,
        typename TPair,
        typename TLimitsString >
    inline typename Size<Pipe<TInput, Multi<Tupler<SIZE, omitLast, TPack>, TPair, TLimitsString> > >::Type
    length(Pipe<TInput, Multi<Tupler<SIZE, omitLast, TPack>, TPair, TLimitsString> > const &me)
    {
        typedef Pipe<TInput, Multi<Tupler<SIZE, omitLast, TPack>, TPair, TLimitsString> >   TPipe;
        typedef typename Size<TPipe>::Type                                                  TSize;

        // count all overlapping q-grams
        TSize count = 0;
        TSize seqs = countSequences(me);
        for (TSize i = 0; i < seqs; ++i)
        {
            TSize seqLength = me.limits[i + 1] - me.limits[i];
            if (seqLength > (SIZE - TuplerNumberOfLastTuples_<SIZE, omitLast>::VALUE))
                count += seqLength - (SIZE - TuplerNumberOfLastTuples_<SIZE, omitLast>::VALUE);
        }
        return count;
    }

    template < typename TInput, unsigned SIZE, bool omitLast, typename TPack >
    inline unsigned
    countSequences(Pipe< TInput, Tupler<SIZE, omitLast, TPack> > const &)
    {
        return 1;
    }

    template <
        typename TInput,
        unsigned SIZE,
        bool omitLast,
        typename TPack,
        typename TPair,
        typename TLimitsString >
    inline unsigned
    countSequences(Pipe< TInput, Multi<Tupler<SIZE, omitLast, TPack>, TPair, TLimitsString> > const &me)
    {
        return length(me.limits) - 1;
    }

}

#endif
