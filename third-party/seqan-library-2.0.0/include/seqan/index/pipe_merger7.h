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

#ifndef SEQAN_HEADER_INDEX_MERGER7_H
#define SEQAN_HEADER_INDEX_MERGER7_H

namespace SEQAN_NAMESPACE_MAIN
{

//namespace SEQAN_NAMESPACE_PIPELINING
//{

    template <int I, typename T = void>
    struct SkewShift_;

    template <int I, typename T = void>
    struct SkewNIndex_;


    template <typename T>
    struct SkewShift_<7, T>
    {
        static const unsigned VALUE[7][7];
    };

    template <typename T>
    struct SkewNIndex_<7, T>
    {
        static const unsigned VALUE[7][7];
    };


    template <typename T>
    const unsigned SkewShift_<7, T>::VALUE[7][7] =
                                  {{0,6,5,6,3,3,5},
                                   {6,0,0,6,0,4,4},
                                   {5,0,0,1,0,1,5},
                                   {6,6,1,0,2,1,2},
                                   {3,0,0,2,0,3,2},
                                   {3,4,1,1,3,0,4},
                                   {5,4,5,2,2,4,0}};


    template <typename T>
    const unsigned SkewNIndex_<7, T>::VALUE[7][7] =
                                  {{0,0,0,0,0,1,2},
                                   {0,0,0,0,1,1,2},
                                   {0,1,1,1,1,2,2},
                                   {0,0,1,1,1,1,2},
                                   {0,0,1,2,2,2,2},
                                   {0,0,0,1,2,2,2},
                                   {0,0,0,0,1,2,2}};


    template <typename TValue>
    struct SkewDCStream
    {
        TValue      i;
        unsigned    stream;
    };

    template <typename TValue>
    std::ostream& operator<<(std::ostream &out, const SkewDCStream<TValue> &s)
    {
        out << "< " << s.i.i1 << " , [ ";
        for(unsigned i = 0; i < 3; ++i)
            out << s.i.i2[i] << " ";
        out << "] , [ ";
        for(unsigned i = 0; i < 6; ++i)
            out << s.i.i3[i] << " ";
        out << "] , " << s.stream << " >";
        return out;
    }


    // greater-operator for two SkewDCStreams
    template <typename TValue>
    struct CompareSkewDCStream :
        public std::binary_function < SkewDCStream<TValue>,
                                        SkewDCStream<TValue>,
                                        bool >
    {
        inline bool operator()(const SkewDCStream<TValue> &a,
                               const SkewDCStream<TValue> &b) const
        {
            typedef typename Value<TValue, 2>::Type     TNameTuple;
            typedef typename Value<TNameTuple>::Type    TName;

            int shft = SkewShift_<7>::VALUE[a.stream][b.stream];

            for (int i = 0; i < shft; ++i)
            {
                if (a.i.i3[i] < b.i.i3[i]) return false;
                if (a.i.i3[i] > b.i.i3[i]) return true;
            }

            TName na = a.i.i2[SkewNIndex_<7>::VALUE[a.stream][shft]];
            TName nb = b.i.i2[SkewNIndex_<7>::VALUE[b.stream][shft]];
            if (na < nb) return false;
            if (na > nb) return true;

            // we get here, only if a septet crosses the border of
            // the single text or a sequence in a multiple sequence text
            return a.stream - 1 > b.stream - 1; // this is NOT the same as a.stream > b.stream, 0 should be bigger than 6
        }
    };


    // greater-operator for two SkewDCStreams (optimized for bit-packed character tuples)
    template <typename T1, typename T2, typename T, unsigned SIZE>
    struct CompareSkewDCStream< Triple<T1, T2, Tuple<T, SIZE, BitPacked<> >, Pack> > :
        public std::binary_function < SkewDCStream<Triple<T1, T2, Tuple<T, SIZE, BitPacked<> >, Pack> >,
                                        SkewDCStream<Triple<T1, T2, Tuple<T, SIZE, BitPacked<> >, Pack> >,
                                        bool >
    {
        typedef Tuple<T,SIZE, BitPacked<> > T3;
        inline bool operator()(const SkewDCStream<Triple<T1, T2, T3, Pack> > &a,
                               const SkewDCStream<Triple<T1, T2, T3, Pack> > &b) const
        {
            typedef typename Value<T2>::Type TName;
            typedef typename Value<T3>::Type TValue3;

            int shft = SkewShift_<7>::VALUE[a.stream][b.stream];
            typename T3::TBitVector mask = ~((1 << ((SIZE - shft) * BitsPerValue<TValue3>::VALUE)) - 1);

            if ((a.i.i3.i & mask) < (b.i.i3.i & mask)) return false;
            if ((a.i.i3.i & mask) > (b.i.i3.i & mask)) return true;

            TName na = a.i.i2[SkewNIndex_<7>::VALUE[a.stream][shft]];
            TName nb = b.i.i2[SkewNIndex_<7>::VALUE[b.stream][shft]];
            if (na < nb) return false;
            if (na > nb) return true;

            // we get here, only if a septet crosses the border of
            // the single text or a sequence in a multiple sequence text
            return a.stream - 1 > b.stream - 1; // this is NOT the same as a.stream > b.stream, 0 should be bigger than 6
        }
    };




    struct Merger7;

    template < typename TInput0, typename TInput3, typename TInput5, typename TInput6, typename TInput124 >
    struct Value< Pipe< Bundle5< TInput0, TInput3, TInput5, TInput6, TInput124 >, Merger7 > >
    {
        typedef typename Value<TInput0>::Type       TInValue0;
        typedef typename Value<TInValue0, 1>::Type  Type;
    };


    template <typename TLimitsString>
    struct Merger7Multi;

    template < typename TInput0, typename TInput3, typename TInput5, typename TInput6, typename TInput124, typename TLimitsString >
    struct Value< Pipe< Bundle5< TInput0, TInput3, TInput5, TInput6, TInput124 >, Merger7Multi<TLimitsString> > >
    {
        typedef typename Value<TInput0>::Type       TInValue0;
        typedef typename Value<TInValue0, 1>::Type  Type;
    };


    //////////////////////////////////////////////////////////////////////////////
    // merger7 class
    template < typename TInput0, typename TInput3, typename TInput5, typename TInput6, typename TInput124 >
    struct Pipe< Bundle5< TInput0, TInput3, TInput5, TInput6, TInput124 >, Merger7 >
    {
        typedef typename Value<TInput0>::Type       TInValue0;
        typedef typename Value<TInput3>::Type       TInValue3;
        typedef typename Value<TInput5>::Type       TInValue5;
        typedef typename Value<TInput6>::Type       TInValue6;
        typedef typename Value<TInput124>::Type     TInValue124;
        typedef typename Size<Pipe>::Type           TSize;

        typedef typename Value<TInValue0, 3>::Type  TInTuple0;
//        typedef typename Value<TInTuple0>::Type     Type;

        typedef SkewDCStream<TInValue0>             TSkewDCStream;

        TSkewDCStream                                inValue[7];
        unsigned                                    rank[5];
        unsigned                                    first;
        CompareSkewDCStream<TInValue0>                streamGreater;

        Bundle5 <
            TInput0,
            TInput3,
            TInput5,
            TInput6,
            TInput124 > in;
        TSize N;

        Pipe(Bundle5< TInput0, TInput3, TInput5, TInput6, TInput124 > _in):
            first(0), in(_in), N(0)
        {}

        template <typename T1, typename T2, typename T3>
        inline static void _copy(TSkewDCStream &dst, Triple<T1,T2,T3, Pack> const &src)
        {
            memcpy(&dst.i.i2, &src.i2, sizeof(T2));
            memcpy(&dst.i.i3, &src.i3, sizeof(T3));
        }

        template <typename T1, typename T2, typename T, unsigned SIZE>
        inline static void _copy(TSkewDCStream &dst, Triple<T1,T2,Tuple<T,SIZE, BitPacked<> >, Pack> const &src)
        {
            typedef typename Value<TInValue0, 3>::Type TCharTuple;
            memcpy(&dst.i.i2, &src.i2, sizeof(T2));
            dst.i.i3 = src.i3;
            dst.i.i3 <<= LENGTH<TCharTuple>::VALUE - SIZE;
        }

        inline typename Value<Pipe>::Type const operator*()
        {
            return getValueI1(inValue[rank[first]].i);
        }

        inline void insertStream(unsigned stream)
        {
            switch (stream) {
                case 0:
                    if (eof(in.in1)) {
                        ++first;
                        return;
                    }
                    inValue[0].i.i1 = N - (*in.in1).i1;
                    _copy(inValue[0], *in.in1);
                    ++in.in1;
                    break;

                case 3:
                    if (eof(in.in2)) {
                        ++first;
                        return;
                    }
                    inValue[3].i.i1 = N - (*in.in2).i1;
                    _copy(inValue[3], *in.in2);
                    ++in.in2;
                    break;

                case 5:
                    if (eof(in.in3)) {
                        ++first;
                        return;
                    }
                    inValue[5].i.i1 = N - (*in.in3).i1;
                    _copy(inValue[5], *in.in3);
                    ++in.in3;
                    break;

                case 6:
                    if (eof(in.in4)) {
                        ++first;
                        return;
                    }
                    inValue[6].i.i1 = N - (*in.in4).i1;
                    _copy(inValue[6], *in.in4);
                    ++in.in4;
                    break;

                default:    // case 1, 2, or 4
                    if (eof(in.in5)) {
                        ++first;
                        return;
                    }
                    // calculate residue class from the suffix length
                    inValue[1].stream = (*in.in5).i1 % 7;
                    inValue[1].i.i1 = N - (*in.in5).i1;
                    _copy(inValue[1], *in.in5);
                    ++in.in5;
            }

            // linear search
            int right;
            for(right = first + 1;  right < 5;  ++right)
                if (!streamGreater(inValue[stream],inValue[rank[right]])) break;

            // remove the least suffix ...
            for(int i = first + 1;  i < right;  ++i)
                rank[i-1] = rank[i];

            // ... and insert the new one
            rank[right-1] = stream;
        }

        Pipe& operator++()
        {
            insertStream(rank[first]);
            return *this;
        }

        void fill()
        {
            inValue[0].stream = 0;
            inValue[3].stream = 3;
            inValue[5].stream = 5;
            inValue[6].stream = 6;

            first = 4;    insertStream(0);
            --first;    insertStream(3);
            --first;    insertStream(5);
            --first;    insertStream(6);
            --first;    insertStream(1);
        }
    };


    //////////////////////////////////////////////////////////////////////////////
    // merger7 class for multiple sequences
    template < typename TInput0, typename TInput3, typename TInput5, typename TInput6, typename TInput124, typename TLimitsString >
    struct Pipe< Bundle5< TInput0, TInput3, TInput5, TInput6, TInput124 >, Merger7Multi<TLimitsString> >
    {
        typedef typename Value<TInput0>::Type       TInValue0;
        typedef typename Value<TInput3>::Type       TInValue3;
        typedef typename Value<TInput5>::Type       TInValue5;
        typedef typename Value<TInput6>::Type       TInValue6;
        typedef typename Value<TInput124>::Type     TInValue124;
        typedef typename Size<Pipe>::Type           TSize;

        typedef typename Value<TInValue0, 3>::Type  TInTuple0;
//        typedef typename Value<TInTuple0>::Type     Type;

        typedef SkewDCStream<TInValue0>             TSkewDCStream;

        TSkewDCStream                                inValue[7];
        unsigned                                    rank[5];
        unsigned                                    first;
        CompareSkewDCStream<TInValue0>                streamGreater;

        Bundle5 <
            TInput0,
            TInput3,
            TInput5,
            TInput6,
            TInput124 > in;
        TLimitsString const &limits;

        Pipe(Bundle5< TInput0, TInput3, TInput5, TInput6, TInput124 > _in, TLimitsString const &_limits):
            in(_in), limits(_limits)
        {}

        template <typename T1, typename T2, typename T3>
        inline static void _copy(TSkewDCStream &dst, Triple<T1,T2,T3, Pack> const &src)
        {
            memcpy(&dst.i.i1, &src.i1, sizeof(T1));
            memcpy(&dst.i.i2, &src.i2, sizeof(T2));
            memcpy(&dst.i.i3, &src.i3, sizeof(T3));
        }

        template <typename T1, typename T2, typename T, unsigned SIZE>
        inline static void _copy(TSkewDCStream &dst, Triple<T1,T2,Tuple<T,SIZE, BitPacked<> >, Pack> const &src)
        {
            typedef typename Value<TInValue0, 3>::Type TCharTuple;
            memcpy(&dst.i.i1, &src.i1, sizeof(T1));
            memcpy(&dst.i.i2, &src.i2, sizeof(T2));
            dst.i.i3 = src.i3;
            dst.i.i3 <<= LENGTH<TCharTuple>::VALUE - SIZE;
        }

        inline typename Value<Pipe>::Type const operator*()
        {
            return getValueI1(inValue[rank[first]].i);
        }

        inline void insertStream(unsigned stream)
        {
            switch (stream) {
                case 0:
                    if (eof(in.in1)) {
                        ++first;
                        return;
                    }
                    _copy(inValue[0], *in.in1);
                    ++in.in1;
                    break;

                case 3:
                    if (eof(in.in2)) {
                        ++first;
                        return;
                    }
                    _copy(inValue[3], *in.in2);
                    ++in.in2;
                    break;

                case 5:
                    if (eof(in.in3)) {
                        ++first;
                        return;
                    }
                    _copy(inValue[5], *in.in3);
                    ++in.in3;
                    break;

                case 6:
                    if (eof(in.in4)) {
                        ++first;
                        return;
                    }
                    _copy(inValue[6], *in.in4);
                    ++in.in4;
                    break;

                default:    // case 1, 2, or 4
                    if (eof(in.in5)) {
                        ++first;
                        return;
                    }
                    // calculate residue class from the suffix length
                    typename Value<TInValue124,1>::Type in_in5_i1 = getValueI1(*in.in5);
                    unsigned seqNo = getValueI1(in_in5_i1);
                    inValue[1].stream = ((limits[seqNo + 1] - limits[seqNo]) - getValueI2(in_in5_i1)) % 7;

                    _copy(inValue[1], *in.in5);
                    ++in.in5;
            }

            // linear search
            int right;
            for(right = first + 1;  right < 5;  ++right)
                if (!streamGreater(inValue[stream],inValue[rank[right]])) break;

            // remove the least suffix ...
            for(int i = first + 1;  i < right;  ++i)
                rank[i-1] = rank[i];

            // ... and insert the new one
            rank[right-1] = stream;
        }

        Pipe& operator++()
        {
            insertStream(rank[first]);
            return *this;
        }

        void fill()
        {
            inValue[0].stream = 0;
            inValue[3].stream = 3;
            inValue[5].stream = 5;
            inValue[6].stream = 6;

            first = 4;    insertStream(0);
            --first;    insertStream(3);
            --first;    insertStream(5);
            --first;    insertStream(6);
            --first;    insertStream(1);
        }
    };


    //////////////////////////////////////////////////////////////////////////////
    // global pipe functions
    template < typename TInput >
    inline bool control(Pipe< TInput, Merger7 > &me, ControlBeginRead const &command)
    {
        if (!control(me.in, command)) return false;
        me.N = length(me);
        me.fill();
        return true;
    }

    template < typename TInput >
    inline bool control(Pipe< TInput, Merger7 > &me, ControlEof const &)
    {
        return me.first == 5;
    }

    template < typename TInput >
    inline bool control(Pipe< TInput, Merger7 > &me, ControlEos const &)
    {
        return control(me, ControlEof());
    }

    template < typename TInput >
    inline typename Size< Pipe< TInput, Merger7 > >::Type
    length(Pipe< TInput, Merger7 > const &me)
    {
        return length(me.in.in1) +
               length(me.in.in2) +
               length(me.in.in3) +
               length(me.in.in4) +
               length(me.in.in5);
    }




    template < typename TInput, typename TLimitsString >
    inline bool control(Pipe< TInput, Merger7Multi<TLimitsString> > &me, ControlBeginRead const &command)
    {
        if (!control(me.in, command)) return false;
        me.fill();
        return true;
    }

    template < typename TInput, typename TLimitsString >
    inline bool control(Pipe< TInput, Merger7Multi<TLimitsString> > &me, ControlEof const &)
    {
        return me.first == 5;
    }

    template < typename TInput, typename TLimitsString >
    inline bool control(Pipe< TInput, Merger7Multi<TLimitsString> > &me, ControlEos const &)
    {
        return control(me, ControlEof());
    }

    template < typename TInput, typename TLimitsString >
    inline typename Size< Pipe< TInput, Merger7Multi<TLimitsString> > >::Type
    length(Pipe< TInput, Merger7Multi<TLimitsString> > const &me) {
        return length(me.in.in1) +
               length(me.in.in2) +
               length(me.in.in3) +
               length(me.in.in4) +
               length(me.in.in5);
    }

//}

}

#endif
