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

#ifndef SEQAN_HEADER_INDEX_BWT_H
#define SEQAN_HEADER_INDEX_BWT_H

namespace SEQAN_NAMESPACE_MAIN
{

//namespace SEQAN_NAMESPACE_PIPELINING
//{

    struct Bwt {};


    //////////////////////////////////////////////////////////////////////////////
    // external Bwt algorithm
    //////////////////////////////////////////////////////////////////////////////


    template < typename TTextInput, typename TSuffixArrayInput >
    struct Value< Pipe< Bundle2< TTextInput, TSuffixArrayInput >, Bwt > > {
        typedef typename Value<TTextInput>::Type Type;
    };

    //////////////////////////////////////////////////////////////////////////////
    // bwt class
    template < typename TTextInput, typename TSuffixArrayInput >
    struct Pipe< Bundle2< TTextInput, TSuffixArrayInput >, Bwt >
    {
        // *** SPECIALIZATION ***

        typedef Pipe< TSuffixArrayInput, Counter > TSA;
                                        typedef typename Size<TTextInput>::Type    TSize;
        typedef Pool< TypeOf_(TSA), MapperSpec< MapperConfigSize< filterI1<TypeOf_(TSA)>, TSize> > > TInverter;
        typedef Pipe< TInverter, Filter< filterI2<TypeOf_(TInverter)> > > TCounterFilter;
        typedef Pipe< TTextInput, Shifter< -1, false > > TShiftText;

        typedef Pipe< Bundle2< TCounterFilter, TShiftText >, Joiner > TJoiner;
        typedef Pool< TypeOf_(TJoiner), MapperSpec< MapperConfigSize< filterI1<TypeOf_(TJoiner)>, TSize> > > TLinearMapper;
        typedef Pipe< TLinearMapper, Filter< filterI2<TypeOf_(TLinearMapper)> > > TFilter;

        TLinearMapper        mapper;
        TFilter                in;

        Pipe():
            in(mapper) {}

        Pipe(Bundle2< TTextInput, TSuffixArrayInput > const &_bundleIn):
            in(mapper)
        {
            process(_bundleIn.in1, _bundleIn.in2);
        }

        template < typename TTextInput_, typename TSuffixArrayInput_ >
        bool process(TTextInput_ &textIn, TSuffixArrayInput_ &suffixArrayIn) {

            // *** INSTANTIATION ***

            TSA                            sa(suffixArrayIn);
            TInverter                    inverter;
            TCounterFilter                filter(inverter);

            #ifdef SEQAN_DEBUG_INDEX
                std::cerr << "  invert suffix array" << std::endl;
            #endif
            inverter << sa;
            SEQAN_PROMARK("Suffix-Array invertiert");

            TShiftText                    shifter(textIn);
            TJoiner                        joiner(bundle2(filter, shifter));

            #ifdef SEQAN_DEBUG_INDEX
                std::cerr << "  de-invert suffix array" << std::endl;
            #endif
            mapper << joiner;
            SEQAN_PROMARK("Suffix-Array linearisiert");

            return true;
        }

        inline typename Value<Pipe>::Type const operator*() const {
            return *in;
        }

        inline Pipe& operator++() {
            ++in;
            return *this;
        }
    };

    // not sure which interface is more intuitive, we support both
    // you can call "skew << pipe" or "skew_t skew(pipe); skew.process()"
    // for the first we would need no _in member
    template < typename TInput, typename TTextInput_, typename TSuffixArrayInput_ >
    inline bool operator<<(Pipe< TInput, Bwt > &me, Bundle2< TTextInput_, TSuffixArrayInput_ > const &bundleIn) {
         return me.process(bundleIn.in1, bundleIn.in2);
    }



    //////////////////////////////////////////////////////////////////////////////
    // external Bwt algorithm (optimized for multiple sequences)
    //////////////////////////////////////////////////////////////////////////////


    template < typename TTextInput, typename TSuffixArrayInput, typename TPair, typename TLimitsString >
    struct Value< Pipe< Bundle2< TTextInput, TSuffixArrayInput >, Multi<Bwt, TPair, TLimitsString> > > {
        typedef typename Value<TTextInput>::Type Type;
    };

    template <typename InType, typename TLimitsString, typename Result = typename Value<TLimitsString>::Type>
    struct _filterGlobalizer : public std::unary_function<InType,Result> {
        TLimitsString const &limits;
        _filterGlobalizer(TLimitsString const &_limits) : limits(_limits) {}
        inline Result operator()(const InType& x) const
        {
            return posGlobalize(x, limits);
        }
    };


    //////////////////////////////////////////////////////////////////////////////
    // bwt class
    template < typename TTextInput, typename TSuffixArrayInput, typename TPair, typename TLimitsString >
    struct Pipe< Bundle2< TTextInput, TSuffixArrayInput >, Multi<Bwt, TPair, TLimitsString> >
    {
        // *** SPECIALIZATION ***

                                        typedef _filterGlobalizer<TypeOf_(TSuffixArrayInput), TLimitsString, TSizeOf_(TSuffixArrayInput)> filter_globalizer_t;
        typedef Pipe< TSuffixArrayInput, Filter<filter_globalizer_t> > TGlobalizer;
        typedef Pipe< TGlobalizer, Counter > TSA;
                                        typedef typename Size<TTextInput>::Type    TSize;
        typedef Pool< TypeOf_(TSA), MapperSpec< MapperConfigSize< filterI1<TypeOf_(TSA)>, TSize> > > TInverter;
        typedef Pipe< TInverter, Filter< filterI2<TypeOf_(TInverter)> > > TCounterFilter;
        typedef Pipe< TTextInput, Shifter< -1, false > > TShiftText;

        typedef Pipe< Bundle2< TCounterFilter, TShiftText >, Joiner > TJoiner;
        typedef Pool< TypeOf_(TJoiner), MapperSpec< MapperConfigSize< filterI1<TypeOf_(TJoiner)>, TSize> > > TLinearMapper;
        typedef Pipe< TLinearMapper, Filter< filterI2<TypeOf_(TLinearMapper)> > > TFilter;

        TTextInput            *textIn;
        TSuffixArrayInput    *suffixArrayIn;
        TLinearMapper        mapper;
        TFilter                in;

        TLimitsString const    &limits;

        Pipe(TLimitsString const &_limits):
            in(mapper),
            limits(_limits)    {}

        Pipe(Bundle2< TTextInput, TSuffixArrayInput > const &_bundleIn, TLimitsString const &_limits):
            textIn(&_bundleIn.in1),
            suffixArrayIn(&_bundleIn.in2),
            in(mapper),
            limits(_limits)
        {
            process();
        }

        inline void process() {
            process(*textIn, *suffixArrayIn);
        }

        template < typename TTextInput_, typename TSuffixArrayInput_ >
        bool process(TTextInput_ &textIn, TSuffixArrayInput_ &suffixArrayIn) {

            // *** INSTANTIATION ***

            for(int i=0;i<length(limits);++i)
                std::cout << limits[i]<<"  ";
            std::cout<<std::endl;

            TGlobalizer                    globalizer(suffixArrayIn, limits);
            TSA                            sa(globalizer);
            TInverter                    inverter;
            TCounterFilter                filter(inverter);

            #ifdef SEQAN_DEBUG_INDEX
                std::cerr << "  invert suffix array" << std::endl;
            #endif
            inverter << sa;
            SEQAN_PROMARK("Suffix-Array invertiert");

            TShiftText                    shifter(textIn);
            TJoiner                        joiner(bundle2(filter, shifter));

            #ifdef SEQAN_DEBUG_INDEX
                std::cerr << "  de-invert suffix array" << std::endl;
            #endif
            mapper << joiner;
            SEQAN_PROMARK("Suffix-Array linearisiert");

            return true;
        }

        inline typename Value<Pipe>::Type const operator*() const {
            return *in;
        }

        inline Pipe& operator++() {
            ++in;
            return *this;
        }
    };

    // not sure which interface is more intuitive, we support both
    // you can call "bwt << pipe" or "bwt_t bwt(pipe); bwt.process()"
    // for the first we would need no _in member
    template < typename TInput, typename TTextInput_, typename TSuffixArrayInput_, typename TPair, typename TLimitsString >
    inline bool operator<<(Pipe< TInput, Multi<Bwt, TPair, TLimitsString> > &me, Bundle2< TTextInput_, TSuffixArrayInput_ > const &bundleIn) {
         return me.process(bundleIn.in1, bundleIn.in2);
    }



    //////////////////////////////////////////////////////////////////////////////
    // internal Bwt algorithm
    //////////////////////////////////////////////////////////////////////////////


    template < typename TBWT,
               typename TText,
               typename TSA >
    void createBWTableInt(
        TBWT &bwt,
        TText const &s,
        TSA const &SA)
    {
        typedef typename Value<TBWT>::Type        TValue;
        typedef typename GetValue<TSA>::Type    TSAValue;
        typedef typename Size<TSA>::Type        TSize;

        TSize n = length(s);

        for (TSize i = 0; i < n; ++i)
        {
            TSAValue sa = getValue(SA, i);
            if (sa != 0)
                bwt[i] = getValue(s, sa - 1);
            else
                bwt[i] = TValue();
//                bwt[i] = getValue(s, n - 1);
        }
    }

    template < typename TBWT,
               typename TString,
               typename TSpec,
               typename TSA >
    void createBWTableInt(
        TBWT &bwt,
        StringSet<TString, TSpec> const &s,
        TSA const &SA)
    {
        typedef typename Value<TBWT>::Type        TValue;
        typedef typename Size<TSA>::Type        TSize;

        TSize n = lengthSum(s);
        Pair<unsigned, typename Size<TString>::Type> loc;

        for (TSize i = 0; i < n; ++i)
        {
            posLocalize(loc, getValue(SA, i), stringSetLimits(s));
            if (loc.i2 != 0)
                bwt[i] = s[loc.i1][loc.i2 - 1];
            else
                if (loc.i1 != 0)
                    bwt[i] = s[loc.i1 - 1][length(s[loc.i1 - 1]) - 1];
                else
                    bwt[i] = TValue();
//                    bwt[i] = getValue(s, n - 1);
        }
    }

//}

}

#endif
