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

#ifndef SEQAN_HEADER_PIPE_SAMPLER_H
#define SEQAN_HEADER_PIPE_SAMPLER_H

namespace SEQAN_NAMESPACE_MAIN
{

//namespace SEQAN_NAMESPACE_PIPELINING
//{

    template <int I, typename T = void>
    struct SkewDC_;

//////////////////////////////////////////////////////////////////////////////

    template < unsigned m, typename TPack = Pack >
    struct Sampler;

    template < typename TInput, unsigned m, typename TPack >
    struct Value< Pipe< TInput, Sampler<m, TPack> > >
    {
        typedef Tuple<typename Value<TInput>::Type, m, TPack>    mTuple;
        typedef Pair<typename Size<TInput>::Type, mTuple, Pack> Type;
    };

//////////////////////////////////////////////////////////////////////////////

    template < typename TInput, unsigned m, typename TPack, typename TPair, typename TLimitsString >
    struct Value< Pipe< TInput, Multi<Sampler<m, TPack>, TPair, TLimitsString> > >
    {
        typedef Tuple<typename Value<TInput>::Type, m, TPack>    mTuple;
        typedef Pair<TPair, mTuple, Pack> Type;
    };

//////////////////////////////////////////////////////////////////////////////

// TODO(holtgrew): Documentation bug with DC!

/*!
 * @class Sampler
 * @extends Pipe
 * @headerfile <seqan/pipe.h>
 * @brief Outputs m-tuples beginning at a position of difference cover DC.
 *
 * @signature template <typename TInput, unsigned M[, typename TPack]>
 *            class Pipe<TInput, Sampler<M, TPack> >;
 *
 * @tparam TInput The type of the pipeline module this module reads from.
 * @tparam m      The tuple size.
 * @tparam TPack  Specifies the packing method of the tuples (<tt>void</tt> = no packing), default is <tt>Pack</tt>.
 *
 * The output type is a Pair of size type and Tuple of input elements and length m (i.e. <tt>Pair&lt;Size&lt;TInput&gt;::Type,
 * Tuple&lt;Value&lt;TInput&gt;::Type, m, TPack&gt; &gt;</tt>).
 *
 * The first output field contains the number of remaining pipe elements. The m-tuple in the second field contains the
 * first m elements of them. The m-tuples are substrings of the input stream beginning at positions <tt>i</tt>, with
 * <tt>(n-i) mod m</tt> is element of the set DC (n is the input stream length).
 *
 * @section Examples
 *
 * The set <tt>{1,2,4}</tt> is represented by <tt>int DC[] = { 3, 1, 2, 4 }</tt>.
 *
 * @see BitPackedTuple
 */

    //////////////////////////////////////////////////////////////////////////////
    // sampler class
    template < typename TInput, unsigned m, typename TPack >
    struct Pipe< TInput, Sampler<m, TPack> >
    {
        typedef typename Value<Pipe>::Type  TOutValue;
        typedef typename Size<Pipe>::Type   TSize;

        TInput        &in;
        bool        filter[m];
        TSize       idx, _size, _rest;
        unsigned    idxMod;
        TOutValue   tmp1, tmp2;
        TOutValue   *outRef, *tmpRef;
        bool        last;

        Pipe(TInput& _in):
            in(_in),
            outRef(&tmp1),
            tmpRef(&tmp2) {}

        inline void prepare()
        {
            for (unsigned i = 0; i < m; ++i)
                filter[i] = 0;
            for(unsigned i = 1; i <= SkewDC_<m>::VALUE[0]; i++)
                filter[SkewDC_<m>::VALUE[i]] = true;

            idx = length(in);
            idxMod = idx % m;

            while (!filter[idxMod] && !eof(in))
            {
                ++in;
                if (idxMod == 0) idxMod = m;
                --idxMod; --idx;
            }
            _rest = length(*this);
            fill();
            swap();
        }

        inline void fill()
        {
            unsigned i;
            for(i = 0; i < m && !eof(in); ++i, ++in)
                tmpRef->i2.i[i] = *in;
            last = eof(in);
            for(; i < m; ++i)
                tmpRef->i2.i[i] = 0;
            tmpRef->i1 = idx;
        }

        inline void rotate(unsigned r)
        {
            for(unsigned i = 0; i < m; ++i, ++r)
            {
                if (r == m) r = 0;
                tmpRef->i2.i[i] = outRef->i2.i[r];
            }
        }

        inline void swap()
        {
            TOutValue *newOutRef = tmpRef;
            tmpRef = outRef;
            outRef = newOutRef;
        }

        inline TOutValue const& operator*()
        {
            return *outRef;
        }

        Pipe& operator++()
        {
            unsigned skipped = 0;
            if (--_rest)
            {
                if (!last)
                {
                    do
                    {
                        outRef->i2.i[skipped++] = *in;
                        ++in;
                        if (idxMod == 0) idxMod = m;
                        --idxMod; --idx;
                        if (eof(in))
                        {
                            last = true;
                            while (!filter[idxMod])
                            {
                                outRef->i2.i[skipped++] = 0;
                                if (idxMod == 0) idxMod = m;
                                --idxMod; --idx;
                            }
                            break;
                        }
                    } while (!filter[idxMod]);
                }
                else
                {
                    do
                    {
                        outRef->i2.i[skipped++] = 0;
                        if (idxMod == 0) idxMod = m;
                        --idxMod; --idx;
                    }
                    while (!filter[idxMod]);
                }
            }
            rotate(skipped);
            tmpRef->i1 = idx;
            swap();
            return *this;
        }
    };


    //////////////////////////////////////////////////////////////////////////////
    // sampler class (uses packing)
    template < typename TInput, unsigned m >
    struct Pipe< TInput, Sampler<m, BitPacked<> > >
    {
        typedef typename Value<Pipe>::Type          TOutValue;
        typedef typename Size<Pipe>::Type           TSize;
        typedef typename Value<TOutValue, 2>::Type  TTuple;

        TInput        &in;
        bool        filter[m];
        TSize       _size, _rest;
        unsigned    idxMod;
        TOutValue   tmp;
        bool        last;

        Pipe(TInput & _in) : in(_in), _size(0), _rest(0), idxMod(0), last(false)
        {}

        inline void prepare()
        {
            for (unsigned i = 0; i < m; ++i)
                filter[i] = 0;
            for (unsigned i = 1; i <= SkewDC_<m>::VALUE[0]; i++)
                filter[SkewDC_<m>::VALUE[i]] = true;

            tmp.i1 = length(in);
            idxMod = tmp.i1 % m;

            while (!filter[idxMod] && !eof(in))
            {
                ++in;
                if (idxMod == 0) idxMod = m;
                --idxMod; --tmp.i1;
            }
            _rest = length(*this);
            fill();
        }

        inline void fill()
        {
            unsigned i;
            clear(tmp.i2);
            for(i = 0; i < m && !eof(in); ++i, ++in)
            {
                tmp.i2 <<= 1;
                tmp.i2 |= *in;
            }
            last = eof(in);
            tmp.i2 <<= (m - i);
        }

        inline TOutValue const& operator*()
        {
            return tmp;
        }

        Pipe& operator++()
        {
            if (--_rest)
            {
                if (!last)
                {
                    do
                    {
                        tmp.i2 <<= 1;
                        tmp.i2 |= *in;
                        ++in;
                        if (idxMod == 0) idxMod = m;
                        --idxMod; --tmp.i1;
                        if (eof(in))
                        {
                            last = true;
                            while (!filter[idxMod])
                            {
                                tmp.i2 <<= 1;
                                if (idxMod == 0) idxMod = m;
                                --idxMod; --tmp.i1;
                            }
                            break;
                        }
                    } while (!filter[idxMod]);
                }
                else
                {
                    do
                    {
                        tmp.i2 <<= 1;
                        if (idxMod == 0) idxMod = m;
                        --idxMod; --tmp.i1;
                    }
                    while (!filter[idxMod]);
                }
            }
            return *this;
        }
    };



    //////////////////////////////////////////////////////////////////////////////
    // global pipe functions
    template < typename TInput, unsigned m, typename TPack >
    inline bool control(Pipe< TInput, Sampler<m, TPack> > &me, ControlBeginRead const &command)
    {
        if (!control(me.in, command)) return false;
        me.prepare();
        return true;
    }

    template < typename TInput, unsigned m, typename TPack >
    inline bool control(Pipe< TInput, Sampler<m, TPack> > &me, ControlEof const &/*command*/)
    {
        return me._rest == 0;
    }

    template < typename TInput, unsigned m, typename TPack >
    inline typename Size< Pipe< TInput, Sampler<m, TPack> > >::Type
    length(Pipe< TInput, Sampler<m, TPack> > const &me)
    {
        typedef typename Size< Pipe< TInput, Sampler<m> > >::Type TSize;

        TSize sum = 0;
        TSize size = length(me.in);

        // sum up the number of tuples in each residue class
        // for a string of length n there are 1+(n-x)/m suffixes
        // whose lengths are in residue class x
        for (unsigned i = 1; i <= SkewDC_<m>::VALUE[0]; ++i)
        {
            sum += (size + m - SkewDC_<m>::VALUE[i]) / m;
            if (SkewDC_<m>::VALUE[i] == 0)
                --sum;  // we don't consider empty suffixes
        }
        return sum;
    }



//////////////////////////////////////////////////////////////////////////////

    template < typename TInput, unsigned m, typename TPack, typename TPair, typename TLimitsString >
    struct Size< Pipe<TInput, Multi<Sampler<m, TPack>, TPair, TLimitsString> > > :
        public Value<TLimitsString> {};

    //////////////////////////////////////////////////////////////////////////////
    // sampler class
    template < typename TInput, unsigned m, typename TPack, typename TPair, typename TLimitsString >
    struct Pipe< TInput, Multi<Sampler<m, TPack>, TPair, TLimitsString> >
    {
        typedef typename Value<Pipe>::Type  TOutValue;
        typedef typename Size<Pipe>::Type   TSize;

        typedef PairDecrementer_<TPair, TLimitsString, m> Decrementer;

        TInput        &in;
        bool        filter[m];
        Decrementer    localPos;
        TSize       _size, _rest;
        TOutValue   tmp1, tmp2;
        TOutValue   *outRef, *tmpRef;
        bool        last;

        TLimitsString const &limits;

        template <typename TLimitsString_>
        Pipe(TInput& _in, TLimitsString_ &_limits):  // const &_limits is intentionally omitted to suppress implicit casts (if types mismatch) and taking refs of them
            in(_in),
            outRef(&tmp1),
            tmpRef(&tmp2),
            limits(_limits) {}

        inline void prepare()
        {
            for (unsigned i = 0; i < m; ++i)
                filter[i] = 0;
            for(unsigned i = 1; i <= SkewDC_<m>::VALUE[0]; i++)
                filter[SkewDC_<m>::VALUE[i]] = true;

            setHost(localPos, limits);

            while (!filter[localPos.residue] && !eof(in))
            {
                ++in;
                --localPos;
            }
            _rest = length(*this);
            fill();
            swap();
        }

        inline void fill()
        {
            unsigned i;
            for(i = 0; i < m && !eof(in); ++i, ++in)
                tmpRef->i2.i[i] = *in;
            last = eof(in);
            for(; i < m; ++i)
                tmpRef->i2.i[i] = 0;
            tmpRef->i1 = localPos;
        }

        inline void rotate(unsigned r)
        {
            for(unsigned i = 0; i < m; ++i, ++r)
            {
                if (r == m) r = 0;
                tmpRef->i2.i[i] = outRef->i2.i[r];
            }
        }

        inline void swap()
        {
            TOutValue *newOutRef = tmpRef;
            tmpRef = outRef;
            outRef = newOutRef;
        }

        inline TOutValue const& operator*()
        {
            return *outRef;
        }

        Pipe& operator++()
        {
            unsigned skipped = 0;
            if (--_rest)
            {
                if (!last)
                {
                    do
                    {
                        outRef->i2.i[skipped++] = *in;
                        ++in;
                        --localPos;
                        if (eof(in))
                        {
                            last = true;
                            while (!filter[localPos.residue])
                            {
                                outRef->i2.i[skipped++] = 0;
                                --localPos;
                            }
                            break;
                        }
                    }
                    while (!filter[localPos.residue]);
                }
                else
                {
                    do
                    {
                        outRef->i2.i[skipped++] = 0;
                        --localPos;
                    }
                    while (!filter[localPos.residue]);
                }
                rotate(skipped);
                tmpRef->i1 = localPos;
                swap();
            }
            return *this;
        }
    };


    //////////////////////////////////////////////////////////////////////////////
    // sampler class (uses bit-packing)
    template < typename TInput, unsigned m, typename TPair, typename TLimitsString >
    struct Pipe< TInput, Multi<Sampler<m, BitPacked<> >, TPair, TLimitsString> >
    {
        typedef typename Value<Pipe>::Type          TOutValue;
        typedef typename Size<Pipe>::Type           TSize;
        typedef typename Value<TOutValue, 2>::Type  TTuple;

        typedef PairDecrementer_<TPair, TLimitsString, m> Decrementer;

        TInput        &in;
        bool        filter[m];
        TSize       _size, _rest;
        Decrementer localPos;
        TOutValue   tmp;
        bool        last;

        TLimitsString const &limits;

        template <typename TLimitsString_>
        Pipe(TInput& _in, TLimitsString_ &_limits):  // const &_limits is intentionally omitted to suppress implicit casts (if types mismatch) and taking refs of them
            in(_in),
            limits(_limits) {}

        inline void prepare()
        {
            for (unsigned i = 0; i < m; ++i)
                filter[i] = 0;
            for(unsigned i = 1; i <= SkewDC_<m>::VALUE[0]; i++)
                filter[SkewDC_<m>::VALUE[i]] = true;

            setHost(localPos, limits);

            while (!filter[localPos.residue] && !eof(in)) {
                ++in;
                --localPos;
            }
            _rest = length(*this);
            fill();
        }

        inline void fill()
        {
            unsigned i;
            clear(tmp.i2);
            for(i = 0; i < m && !eof(in); ++i, ++in) {
                tmp.i2 <<= 1;
                tmp.i2 |= *in;
            }
            last = eof(in);
            tmp.i2 <<= (m - i);
            tmp.i1 = localPos;
        }

        inline TOutValue const& operator*()
        {
            return tmp;
        }

        Pipe& operator++()
        {
            if (--_rest)
            {
                if (!last)
                {
                    do
                    {
                        tmp.i2 <<= 1;
                        tmp.i2 |= *in;
                        ++in;
                        --localPos;
                        if (eof(in))
                        {
                            last = true;
                            while (!filter[localPos.residue])
                            {
                                tmp.i2 <<= 1;
                                --localPos;
                            }
                            break;
                        }
                    }
                    while (!filter[localPos.residue]);
                }
                else
                {
                    do
                    {
                        tmp.i2 <<= 1;
                        --localPos;
                    } while (!filter[localPos.residue]);
                }
            }
            tmp.i1 = localPos;
            return *this;
        }
    };



    //////////////////////////////////////////////////////////////////////////////
    // global pipe functions
    template < typename TInput, unsigned m, typename TPack, typename TPair, typename TLimitsString >
    inline bool control(Pipe< TInput, Multi<Sampler<m, TPack>, TPair, TLimitsString> > &me, ControlBeginRead const &command)
    {
        if (!control(me.in, command))
            return false;
        me.prepare();
        return true;
    }

    template < typename TInput, unsigned m, typename TPack, typename TPair, typename TLimitsString >
    inline bool control(Pipe< TInput, Multi<Sampler<m, TPack>, TPair, TLimitsString> > &me, ControlEof const &/*command*/)
    {
        return me._rest == 0;
    }

    template < typename TInput, unsigned m, typename TPack, typename TPair, typename TLimitsString >
    inline bool control(Pipe< TInput, Multi<Sampler<m, TPack>, TPair, TLimitsString> > &me, ControlEos const &/*command*/)
    {
        return control(me, ControlEof());
    }

    template < typename TInput, unsigned m, typename TPack, typename TPair, typename TLimitsString >
    inline typename Size< Pipe< TInput, Multi<Sampler<m, TPack>, TPair, TLimitsString> > >::Type
    length(Pipe< TInput, Multi<Sampler<m, TPack>, TPair, TLimitsString> > const &me)
    {
        typedef typename Size< Pipe< TInput, Multi<Sampler<m, TPack>, TPair, TLimitsString> > >::Type TSize;

        if (empty(me.limits))
            return 0;

        TSize sum = 0;
        TLimitsString const &limits = me.limits;
        __int64 seqCountPlusOne = length(me.limits);

        SEQAN_OMP_PRAGMA(parallel for reduction(+:sum))
        for (__int64 i = 1; i < seqCountPlusOne; ++i)
        {
            TSize prev = limits[i - 1];
            TSize cur = limits[i];

            SEQAN_ASSERT_LEQ(prev, cur);
            TSize size = cur - prev;

            // sum up the number of tuples in each residue class
            // for a string of length n there are 1+(n-x)/m suffixes
            // whose lengths are in residue class x
            for (unsigned i = 1; i <= SkewDC_<m>::VALUE[0]; ++i)
            {
                sum += (size + m - SkewDC_<m>::VALUE[i]) / m;
                if (SkewDC_<m>::VALUE[i] == 0)
                    --sum;  // we don't consider empty suffixes
            }
        }

        return sum;
    }
//}

}

#endif
