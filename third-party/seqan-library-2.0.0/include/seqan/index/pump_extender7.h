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

#ifndef SEQAN_HEADER_PUMP_EXTENDER7_H
#define SEQAN_HEADER_PUMP_EXTENDER7_H

namespace SEQAN_NAMESPACE_MAIN
{

//namespace SEQAN_NAMESPACE_PIPELINING
//{

//////////////////////////////////////////////////////////////////////////////

    template <typename TPack = void>
    struct Extender7;

    template <typename TTextInput, typename TNameInput, typename TPack >
    struct Pipe< Bundle2< TTextInput, TNameInput >, Extender7<TPack> >
    {
        typedef typename Size<Pipe>::Type                   TSize;
        typedef typename Value<TTextInput>::Type            TValue;
        typedef typename Value<TNameInput>::Type            TNameInputValue;
        typedef typename Value<TNameInputValue, 2>::Type    TName;

        typedef Tuple<TValue, 4, TPack>                 X4Tuple;
        typedef Tuple<TValue, 5, TPack>                 X5Tuple;
        typedef Tuple<TValue, 6, TPack>                 X6Tuple;
        typedef Tuple<TName, 3>                         NTuple;
        typedef Triple<TSize, NTuple, X6Tuple, Pack>    OutType0;
        typedef Triple<TSize, NTuple, X6Tuple, Pack>    OutType3;
        typedef Triple<TSize, NTuple, X4Tuple, Pack>    OutType5;
        typedef Triple<TSize, NTuple, X5Tuple, Pack>    OutType6;
        typedef Triple<TSize, NTuple, X6Tuple, Pack>    OutType124;

        // pipeline interfaces to ease specialization
        typedef Pipe< void, AbstractSource< OutType0, TSize > > Out0;
        typedef Pipe< void, AbstractSource< OutType3, TSize > > Out3;
        typedef Pipe< void, AbstractSource< OutType5, TSize > > Out5;
        typedef Pipe< void, AbstractSource< OutType6, TSize > > Out6;
        typedef Pipe< void, AbstractSource< OutType124, TSize > > Out124;
    };


//////////////////////////////////////////////////////////////////////////////


    // little copy helper
    // which is needed for bitpacked structs (no =sign possible)
    template <typename Dest, typename Ofs, typename Src>
    static finline Src const cp___(Dest &dst, Ofs const ofs, Src const src) {
        return dst.i3.assignValue(ofs, src);
    }

    template < typename TTextInput, typename TNameInput,
               typename TOut0, typename TOut3, typename TOut5, typename TOut6, typename TOut124 >
    static bool _skew7Extend(TTextInput &textIn, TNameInput &nameIn,
                             TOut0 &out0, TOut3 &out3, TOut5 &out5, TOut6 &out6, TOut124 &out124)
    {
        resize(out0, length(textIn) / 7);
        resize(out3, (length(textIn) + 4) / 7);
        resize(out5, (length(textIn) + 2) / 7);
        resize(out6, (length(textIn) + 1) / 7);
        resize(out124, length(nameIn));
        if (!(
            beginRead(textIn) &&
            beginRead(nameIn) &&
            beginWrite(out0) &&
            beginWrite(out3) &&
            beginWrite(out5) &&
            beginWrite(out6) &&
            beginWrite(out124))) return false;

        // not necessary, but this hides an 'uninitialized' warning...
        typename Value<TOut0>::Type   o0;
        typename Value<TOut124>::Type o1;
        o1.i1 = typename Value<typename Value<TOut124>::Type, 1>::Type();
        o1.i2 = typename Value<typename Value<TOut124>::Type, 2>::Type();
        o1.i3 = typename Value<typename Value<TOut124>::Type, 3>::Type();

        typename Value<TOut124>::Type o2 = o1;
        typename Value<TOut3>::Type   o3;
        typename Value<TOut124>::Type o4 = o1;
        typename Value<TOut5>::Type   o5;
        typename Value<TOut6>::Type   o6;

        o0.i1 = typename Value<typename Value<TOut0>::Type, 1>::Type();
        o0.i2 = typename Value<typename Value<TOut0>::Type, 2>::Type();
        o0.i3 = typename Value<typename Value<TOut0>::Type, 3>::Type();

        o3.i1 = typename Value<typename Value<TOut3>::Type, 1>::Type();
        o3.i2 = typename Value<typename Value<TOut3>::Type, 2>::Type();
        o3.i3 = typename Value<typename Value<TOut3>::Type, 3>::Type();

        o5.i1 = typename Value<typename Value<TOut5>::Type, 1>::Type();
        o5.i2 = typename Value<typename Value<TOut5>::Type, 2>::Type();
        o5.i3 = typename Value<typename Value<TOut5>::Type, 3>::Type();

        o6.i1 = typename Value<typename Value<TOut6>::Type, 1>::Type();
        o6.i2 = typename Value<typename Value<TOut6>::Type, 2>::Type();
        o6.i3 = typename Value<typename Value<TOut6>::Type, 3>::Type();

        typename Size<TTextInput>::Type p = length(textIn);
        unsigned r = (unsigned)(p % 7);


        // BEGIN I: PREFILL

        switch (r) {
        case 6:
/* 6 */                                                                    cp___(o6,0,    *textIn); ++textIn; o6.i1 = p--;

        case 5:
/* 5 */                                                         cp___(o5,0,cp___(o6,1,    *textIn)); ++textIn; o5.i1 = p--;

        case 4:
/* 4 */                                                 o4.i2[0] = o5.i2[0] = o6.i2[0] = (*nameIn).i2; ++nameIn; o4.i1 = p--;
                                                     cp___(o4,0,cp___(o5,1,cp___(o6,2,    *textIn))); ++textIn;

        case 3:
/* 3 */                                   cp___(o3,0,cp___(o4,1,cp___(o5,2,cp___(o6,3,   *textIn)))); ++textIn; o3.i1 = p--;

        case 2:
/* 2 */                           o2.i2[0] = o3.i2[0] = o4.i2[1] = o5.i2[1] = o6.i2[1] = (*nameIn).i2; ++nameIn; o2.i1 = p--;
                               cp___(o2,0,cp___(o3,1,cp___(o4,2,cp___(o5,3,cp___(o6,4,    *textIn))))); ++textIn;

        case 1:
/* 1 */                o1.i2[0] = o2.i2[1] = o3.i2[1] = o4.i2[2] = o5.i2[2] = o6.i2[2] = (*nameIn).i2; ++nameIn; o1.i1 = p--;
                    cp___(o1,0,cp___(o2,1,cp___(o3,2,                                     *textIn))); ++textIn;
            if (r >= 6) push(out6, o6);
            if (r >= 5) push(out5, o5);
            if (r >= 4) push(out124, o4);

        case 0:;
        }

            // BEGIN II: PREFILL And PUSH FULLY FILLED TRIPLES

        if (!eof(nameIn)) {
/* 0 */  cp___(o0,0,cp___(o1,1,cp___(o2,2,cp___(o3,3,                                     *textIn)))); ++textIn; o0.i1 = p--;

/* 6 */  cp___(o0,1,cp___(o1,2,cp___(o2,3,cp___(o3,4,                      cp___(o6,0,    *textIn))))); ++textIn; o6.i1 = p--;

/* 5 */  cp___(o0,2,cp___(o1,3,cp___(o2,4,cp___(o3,5,           cp___(o5,0,cp___(o6,1,    *textIn)))))); ++textIn; o5.i1 = p--;

/* 4 */     o0.i2[0] = o1.i2[1] = o2.i2[2] = o3.i2[2] = o4.i2[0] = o5.i2[0] = o6.i2[0] = (*nameIn).i2; ++nameIn; o4.i1 = p--;
         cp___(o0,3,cp___(o1,4,                      cp___(o4,0,cp___(o5,1,cp___(o6,2,    *textIn))))); ++textIn;
            if (r >= 3) push(out3, o3);
            if (r >= 2) push(out124, o2);

/* 3 */  cp___(o0,4,cp___(o1,5,           cp___(o3,0,cp___(o4,1,cp___(o5,2,cp___(o6,3,    *textIn)))))); ++textIn; o3.i1 = p--;

/* 2 */     o0.i2[1] = o1.i2[2] = o2.i2[0] = o3.i2[0] = o4.i2[1] = o5.i2[1] = o6.i2[1] = (*nameIn).i2; ++nameIn; o2.i1 = p--;
         cp___(o0,5,           cp___(o2,0,cp___(o3,1,cp___(o4,2,cp___(o5,3,cp___(o6,4,    *textIn)))))); ++textIn;
            if (r >= 1) push(out124, o1);

/* 1 */     o0.i2[2] = o1.i2[0] = o2.i2[1] = o3.i2[1] = o4.i2[2] = o5.i2[2] = o6.i2[2] = (*nameIn).i2; ++nameIn; o1.i1 = p--;
                    cp___(o1,0,cp___(o2,1,cp___(o3,2,                                     *textIn))); ++textIn;
            push(out0, o0);
            push(out6, o6);
            push(out5, o5);
            push(out124, o4);
            r = 7;
        }

        // MAIN Loop: PUSH FULLY FILLED TRIPLES

        while (!eof(nameIn)) {
/* 0 */  cp___(o0,0,cp___(o1,1,cp___(o2,2,cp___(o3,3,                                     *textIn)))); ++textIn; o0.i1 -= 7;

/* 6 */  cp___(o0,1,cp___(o1,2,cp___(o2,3,cp___(o3,4,                      cp___(o6,0,    *textIn))))); ++textIn; o6.i1 -= 7;

/* 5 */  cp___(o0,2,cp___(o1,3,cp___(o2,4,cp___(o3,5,           cp___(o5,0,cp___(o6,1,    *textIn)))))); ++textIn; o5.i1 -= 7;

/* 4 */     o0.i2[0] = o1.i2[1] = o2.i2[2] = o3.i2[2] = o4.i2[0] = o5.i2[0] = o6.i2[0] = (*nameIn).i2; ++nameIn; o4.i1 -= 7;
         cp___(o0,3,cp___(o1,4,                      cp___(o4,0,cp___(o5,1,cp___(o6,2,    *textIn))))); ++textIn;
            push(out3, o3);
            push(out124, o2);

/* 3 */  cp___(o0,4,cp___(o1,5,           cp___(o3,0,cp___(o4,1,cp___(o5,2,cp___(o6,3,    *textIn)))))); ++textIn; o3.i1 -= 7;

/* 2 */     o0.i2[1] = o1.i2[2] = o2.i2[0] = o3.i2[0] = o4.i2[1] = o5.i2[1] = o6.i2[1] = (*nameIn).i2; ++nameIn; o2.i1 -= 7;
         cp___(o0,5,           cp___(o2,0,cp___(o3,1,cp___(o4,2,cp___(o5,3,cp___(o6,4,    *textIn)))))); ++textIn;
            push(out124, o1);

/* 1 */     o0.i2[2] = o1.i2[0] = o2.i2[1] = o3.i2[1] = o4.i2[2] = o5.i2[2] = o6.i2[2] = (*nameIn).i2; ++nameIn; o1.i1 -= 7;
                    cp___(o1,0,cp___(o2,1,cp___(o3,2,                                     *textIn))); ++textIn;
            push(out0, o0);
            push(out6, o6);
            push(out5, o5);
            push(out124, o4);
        }

        // END: FLUSH PARTIALLY FILLED TRIPLES

        {
/* 0 */             cp___(o1,1,cp___(o2,2,cp___(o3,3,   0)));

/* 6 */             cp___(o1,2,cp___(o2,3,cp___(o3,4,   0)));

/* 5 */             cp___(o1,3,cp___(o2,4,cp___(o3,5,   0)));

/* 4 */                o1.i2[1] = o2.i2[2] = o3.i2[2] = 0;
                    cp___(o1,4,                         0);
            if (r >= 3) push(out3, o3);
            if (r >= 2) push(out124, o2);

/* 3 */             cp___(o1,5,                         0);

/* 2 */                o1.i2[2] =                        0;
            if (r >= 1) push(out124, o1);

        }

        endWrite(out124);
        endWrite(out6);
        endWrite(out5);
        endWrite(out3);
        endWrite(out0);
        endRead(nameIn);
        endRead(textIn);
        return false;
    }


//////////////////////////////////////////////////////////////////////////////

    template <typename TPair, typename TPack = void>
    struct Extender7Multi;

    template <typename TTextInput, typename TNameInput, typename TPair, typename TPack >
    struct Pipe< Bundle2< TTextInput, TNameInput >, Extender7Multi<TPair, TPack> >
    {
        typedef typename Size<Pipe>::Type                   TSize;
        typedef typename Value<TTextInput>::Type            TValue;
        typedef typename Value<TNameInput>::Type            TNameInputValue;
        typedef typename Value<TNameInputValue, 2>::Type    TName;

        typedef Tuple<TValue, 4, TPack>                 X4Tuple;
        typedef Tuple<TValue, 5, TPack>                 X5Tuple;
        typedef Tuple<TValue, 6, TPack>                 X6Tuple;
        typedef Tuple<TName, 3>                         NTuple;
        typedef Triple<TPair, NTuple, X6Tuple, Pack>    OutType0;
        typedef Triple<TPair, NTuple, X6Tuple, Pack>    OutType3;
        typedef Triple<TPair, NTuple, X4Tuple, Pack>    OutType5;
        typedef Triple<TPair, NTuple, X5Tuple, Pack>    OutType6;
        typedef Triple<TPair, NTuple, X6Tuple, Pack>    OutType124;

        // pipeline interfaces to ease specialization
        typedef Pipe< void, AbstractSource< OutType0, TSize > > Out0;
        typedef Pipe< void, AbstractSource< OutType3, TSize > > Out3;
        typedef Pipe< void, AbstractSource< OutType5, TSize > > Out5;
        typedef Pipe< void, AbstractSource< OutType6, TSize > > Out6;
        typedef Pipe< void, AbstractSource< OutType124, TSize > > Out124;
    };


//////////////////////////////////////////////////////////////////////////////


    template < typename TTextInput, typename TLimitsString, typename TNameInput,
               typename TOut0, typename TOut3, typename TOut5, typename TOut6, typename TOut124 >
    static bool _skew7ExtendMulti(
        TTextInput &textIn, TLimitsString const &limits,
        TNameInput &nameIn1, TNameInput &nameIn2, TNameInput &nameIn4,
        TOut0 &out0, TOut3 &out3, TOut5 &out5, TOut6 &out6, TOut124 &out124)
    {
        typedef typename Value<TLimitsString>::Type TSize;
        typename Iterator<TLimitsString const>::Type it = begin(limits), itEnd = end(limits);

        if (it == itEnd) return true;

        {
            TSize n0 = 0, n3 = 0, n5 = 0, n6 = 0, n124 = 0;

            // count the numbers of septets in residue class 1, 2, and 4

            TSize size;
            TSize old = *it; ++it;

            while (it != itEnd) {
                size = *it - old;
                old = *it;

                n0   +=  size      / 7;
                n3   += (size + 4) / 7;
                n5   += (size + 2) / 7;
                n6   += (size + 1) / 7;
                n124 += (size + 6) / 7 + (size + 5) / 7 + (size + 3) / 7;

                ++it;
            }

            resize(out0, n0);
            resize(out3, n3);
            resize(out5, n5);
            resize(out6, n6);
            resize(out124, n124);
        }


        if (!(
            beginRead(textIn) &&
            beginRead(nameIn1) &&
            beginRead(nameIn2) &&
            beginRead(nameIn4) &&
            beginWrite(out0) &&
            beginWrite(out3) &&
            beginWrite(out5) &&
            beginWrite(out6) &&
            beginWrite(out124))) return false;

        // not necessary, but this hides an 'uninitialized' warning...
        typename Value<TOut0>::Type   o0;
        typename Value<TOut124>::Type o1;
        o1.i1 = typename Value<typename Value<TOut124>::Type, 1>::Type();
        o1.i2 = typename Value<typename Value<TOut124>::Type, 2>::Type();
        o1.i3 = typename Value<typename Value<TOut124>::Type, 3>::Type();

        typename Value<TOut124>::Type o2 = o1;
        typename Value<TOut3>::Type   o3;
        typename Value<TOut124>::Type o4 = o1;
        typename Value<TOut5>::Type   o5;
        typename Value<TOut6>::Type   o6;

        o0.i1 = typename Value<typename Value<TOut0>::Type, 1>::Type();
        o0.i2 = typename Value<typename Value<TOut0>::Type, 2>::Type();
        o0.i3 = typename Value<typename Value<TOut0>::Type, 3>::Type();

        o3.i1 = typename Value<typename Value<TOut3>::Type, 1>::Type();
        o3.i2 = typename Value<typename Value<TOut3>::Type, 2>::Type();
        o3.i3 = typename Value<typename Value<TOut3>::Type, 3>::Type();

        o5.i1 = typename Value<typename Value<TOut5>::Type, 1>::Type();
        o5.i2 = typename Value<typename Value<TOut5>::Type, 2>::Type();
        o5.i3 = typename Value<typename Value<TOut5>::Type, 3>::Type();

        o6.i1 = typename Value<typename Value<TOut6>::Type, 1>::Type();
        o6.i2 = typename Value<typename Value<TOut6>::Type, 2>::Type();
        o6.i3 = typename Value<typename Value<TOut6>::Type, 3>::Type();

        it = begin(limits);

        TSize rounds, cur;
        TSize old = *it; ++it;

        typedef typename Value<TOut0>::Type         TOutValue0;
        typedef typename Value<TOutValue0, 1>::Type TPair;

        PairIncrementer_<TPair, TLimitsString> p;
        setHost(p, limits);

        while (it != itEnd) {
            cur = *it; ++it;
            unsigned r = (unsigned)((cur - old) % 7);
            rounds = (cur - old) / 7;
            old = cur;


            // BEGIN I: PREFILL

            switch (r) {
            case 6:
    /* 6 */                                                                    cp___(o6,0,    *textIn); ++textIn; o6.i1 = p; ++p;

            case 5:
    /* 5 */                                                         cp___(o5,0,cp___(o6,1,    *textIn)); ++textIn; o5.i1 = p; ++p;

            case 4:
    /* 4 */                                                 o4.i2[0] = o5.i2[0] = o6.i2[0] = (*nameIn4).i2; ++nameIn4; o4.i1 = p; ++p;
                                                        cp___(o4,0,cp___(o5,1,cp___(o6,2,    *textIn))); ++textIn;

            case 3:
    /* 3 */                                   cp___(o3,0,cp___(o4,1,cp___(o5,2,cp___(o6,3,   *textIn)))); ++textIn; o3.i1 = p; ++p;

            case 2:
    /* 2 */                           o2.i2[0] = o3.i2[0] = o4.i2[1] = o5.i2[1] = o6.i2[1] = (*nameIn2).i2; ++nameIn2; o2.i1 = p; ++p;
                                cp___(o2,0,cp___(o3,1,cp___(o4,2,cp___(o5,3,cp___(o6,4,    *textIn))))); ++textIn;

            case 1:
    /* 1 */                o1.i2[0] = o2.i2[1] = o3.i2[1] = o4.i2[2] = o5.i2[2] = o6.i2[2] = (*nameIn1).i2; ++nameIn1; o1.i1 = p; ++p;
                        cp___(o1,0,cp___(o2,1,cp___(o3,2,                                     *textIn))); ++textIn;
                if (r >= 6) push(out6, o6);
                if (r >= 5) push(out5, o5);
                if (r >= 4) push(out124, o4);

            case 0:;
            }

                // BEGIN II: PREFILL And PUSH FULLY FILLED TRIPLES

            if (rounds != 0) {
    /* 0 */  cp___(o0,0,cp___(o1,1,cp___(o2,2,cp___(o3,3,                                     *textIn)))); ++textIn; o0.i1 = p; ++p;

    /* 6 */  cp___(o0,1,cp___(o1,2,cp___(o2,3,cp___(o3,4,                      cp___(o6,0,    *textIn))))); ++textIn; o6.i1 = p; ++p;

    /* 5 */  cp___(o0,2,cp___(o1,3,cp___(o2,4,cp___(o3,5,           cp___(o5,0,cp___(o6,1,    *textIn)))))); ++textIn; o5.i1 = p; ++p;

    /* 4 */     o0.i2[0] = o1.i2[1] = o2.i2[2] = o3.i2[2] = o4.i2[0] = o5.i2[0] = o6.i2[0] = (*nameIn4).i2; ++nameIn4; o4.i1 = p; ++p;
            cp___(o0,3,cp___(o1,4,                      cp___(o4,0,cp___(o5,1,cp___(o6,2,    *textIn))))); ++textIn;
                if (r >= 3) push(out3, o3);
                if (r >= 2) push(out124, o2);

    /* 3 */  cp___(o0,4,cp___(o1,5,           cp___(o3,0,cp___(o4,1,cp___(o5,2,cp___(o6,3,    *textIn)))))); ++textIn; o3.i1 = p; ++p;

    /* 2 */     o0.i2[1] = o1.i2[2] = o2.i2[0] = o3.i2[0] = o4.i2[1] = o5.i2[1] = o6.i2[1] = (*nameIn2).i2; ++nameIn2; o2.i1 = p; ++p;
            cp___(o0,5,           cp___(o2,0,cp___(o3,1,cp___(o4,2,cp___(o5,3,cp___(o6,4,    *textIn)))))); ++textIn;
                if (r >= 1) push(out124, o1);

    /* 1 */     o0.i2[2] = o1.i2[0] = o2.i2[1] = o3.i2[1] = o4.i2[2] = o5.i2[2] = o6.i2[2] = (*nameIn1).i2; ++nameIn1; o1.i1 = p; ++p;
                        cp___(o1,0,cp___(o2,1,cp___(o3,2,                                     *textIn))); ++textIn;
                push(out0, o0);
                push(out6, o6);
                push(out5, o5);
                push(out124, o4);
                r = 7;
                --rounds;
            }

            // MAIN Loop: PUSH FULLY FILLED TRIPLES

            while (rounds != 0) {
    /* 0 */  cp___(o0,0,cp___(o1,1,cp___(o2,2,cp___(o3,3,                                     *textIn)))); ++textIn; o0.i1 = p; ++p;

    /* 6 */  cp___(o0,1,cp___(o1,2,cp___(o2,3,cp___(o3,4,                      cp___(o6,0,    *textIn))))); ++textIn; o6.i1 = p; ++p;

    /* 5 */  cp___(o0,2,cp___(o1,3,cp___(o2,4,cp___(o3,5,           cp___(o5,0,cp___(o6,1,    *textIn)))))); ++textIn; o5.i1 = p; ++p;

    /* 4 */     o0.i2[0] = o1.i2[1] = o2.i2[2] = o3.i2[2] = o4.i2[0] = o5.i2[0] = o6.i2[0] = (*nameIn4).i2; ++nameIn4; o4.i1 = p; ++p;
            cp___(o0,3,cp___(o1,4,                      cp___(o4,0,cp___(o5,1,cp___(o6,2,    *textIn))))); ++textIn;
                push(out3, o3);
                push(out124, o2);

    /* 3 */  cp___(o0,4,cp___(o1,5,           cp___(o3,0,cp___(o4,1,cp___(o5,2,cp___(o6,3,    *textIn)))))); ++textIn; o3.i1 = p; ++p;

    /* 2 */     o0.i2[1] = o1.i2[2] = o2.i2[0] = o3.i2[0] = o4.i2[1] = o5.i2[1] = o6.i2[1] = (*nameIn2).i2; ++nameIn2; o2.i1 = p; ++p;
            cp___(o0,5,           cp___(o2,0,cp___(o3,1,cp___(o4,2,cp___(o5,3,cp___(o6,4,    *textIn)))))); ++textIn;
                push(out124, o1);

    /* 1 */     o0.i2[2] = o1.i2[0] = o2.i2[1] = o3.i2[1] = o4.i2[2] = o5.i2[2] = o6.i2[2] = (*nameIn1).i2; ++nameIn1; o1.i1 = p; ++p;
                        cp___(o1,0,cp___(o2,1,cp___(o3,2,                                     *textIn))); ++textIn;
                push(out0, o0);
                push(out6, o6);
                push(out5, o5);
                push(out124, o4);
                --rounds;
            }

            // END: FLUSH PARTIALLY FILLED TRIPLES

            {
    /* 0 */             cp___(o1,1,cp___(o2,2,cp___(o3,3,   0)));

    /* 6 */             cp___(o1,2,cp___(o2,3,cp___(o3,4,   0)));

    /* 5 */             cp___(o1,3,cp___(o2,4,cp___(o3,5,   0)));

    /* 4 */                o1.i2[1] = o2.i2[2] = o3.i2[2] = 0;
                        cp___(o1,4,                         0);
                if (r >= 3) push(out3, o3);
                if (r >= 2) push(out124, o2);

    /* 3 */             cp___(o1,5,                         0);

    /* 2 */                o1.i2[2] =                        0;
                if (r >= 1) push(out124, o1);

            }
        }

        endWrite(out124);
        endWrite(out6);
        endWrite(out5);
        endWrite(out3);
        endWrite(out0);
        endRead(nameIn4);
        endRead(nameIn2);
        endRead(nameIn1);
        endRead(textIn);
        return false;
    }

//}

}

#endif
