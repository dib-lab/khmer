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

#ifndef SEQAN_HEADER_INDEX_SKEW7_MULTI_H
#define SEQAN_HEADER_INDEX_SKEW7_MULTI_H

namespace SEQAN_NAMESPACE_MAIN
{

    //////////////////////////////////////////////////////////////////////////////
    // external Skew7 algorithm (optimized for multiple sequences)
    //  - creates a suffix array of pairs (seqNo,seqPos)
    //////////////////////////////////////////////////////////////////////////////


    template < typename TInput, typename TPair, typename TLimitsString >
    struct Value< Pipe< TInput, Multi<Skew7, TPair, TLimitsString> > >
    {
        typedef TPair Type;
    };


    // *** COMPARATORS & MAPS ***

    template <typename TValue, typename TResult = int>
    struct _skew7NCompMulti : public std::binary_function<TValue, TValue, TResult> {
        inline TResult operator() (const TValue &a, const TValue &b) const
        {
            typedef typename Value<TValue, 1>::Type                 TPos;
            typedef typename Value<TPos, 2>::Type                   TSize;
            typedef typename Value<TValue, 2>::Type                 TSeptet;
            typedef typename Value<TSeptet>::Type                   TSeptetValue;
            typedef typename StoredTupleValue_<TSeptetValue>::Type  TStoredValue;

            const TStoredValue *sa = a.i2.i;
            const TStoredValue *sb = b.i2.i;

            TSize n = LENGTH<TSeptet>::VALUE;
            TSize aLeft = getValueI2(getValueI1(a));
            TSize bLeft = getValueI2(getValueI1(b));
            if (aLeft < n) n = aLeft;
            if (bLeft < n) n = bLeft;

            // compare the overlap of septets (the first n bases)
            for(TSize i = 0; i < n; i++, ++sa, ++sb)
            {
                if (*sa == *sb) continue;
                return (*sa < *sb)? -1 : 1;
            }

            if (n < LENGTH<TSeptet>::VALUE)
            {
                // overlap is equal, now suffix lengths decide
                if (aLeft < bLeft) return -1;
                if (aLeft > bLeft) return 1;
                return (getValueI1(getValueI1(a)) > getValueI1(getValueI1(b)))? -1 : 1;
            }

            return 0;
        }
    };

    // optimized for bitvectors
    template <typename T1, typename TValue, typename TResult>
    struct _skew7NCompMulti< Pair<T1, Tuple<TValue, 7, BitPacked<> >, Pack >, TResult > :
        public std::binary_function<
            Pair<T1, Tuple<TValue, 7, BitPacked<> >, Pack >,
            Pair<T1, Tuple<TValue, 7, BitPacked<> >, Pack >,
            TResult> {
        inline TResult operator() (
            const Pair<T1, Tuple<TValue, 7, BitPacked<> >, Pack > &a,
            const Pair<T1, Tuple<TValue, 7, BitPacked<> >, Pack > &b) const
        {
            typedef typename Value<T1, 2>::Type TSize;
            TSize aLeft = getValueI2(getValueI1(a));
            TSize bLeft = getValueI2(getValueI1(b));

            if (aLeft >= 7 && bLeft >= 7)
            {
                if (a.i2 < b.i2) return -1;
                if (a.i2 > b.i2) return 1;
                return 0;
            }

            TSize n = 7;
            if (aLeft < n) n = aLeft;
            if (bLeft < n) n = bLeft;

            // compare the overlap of septets (the first n bases)
            for(TSize i = 0; i < n; i++)
            {
                if (a.i2[i] == b.i2[i]) continue;
                return (a.i2[i] < b.i2[i])? -1 : 1;
            }

            // overlap is equal, now suffix lengths decide
            if (aLeft < bLeft) return -1;
            if (aLeft > bLeft) return 1;
            return (getValueI1(getValueI1(a)) > getValueI1(getValueI1(b)))? -1 : 1;
        }
    };

    template <
        typename TValue,
        typename TLimitsString,
        typename TResultSize = typename Value<TLimitsString>::Type,
        typename TResult = Pair<TResultSize, typename Value<TValue,2>::Type, typename Spec<TValue>::Type>
    >
    struct _skew7GlobalSlicedMulti :
        public std::unary_function<TValue, TResult>
    {

        typedef TResultSize TSize;

        TSize            n4, n2, n1;
        TSize            n24;
        TLimitsString    off[5];

        _skew7GlobalSlicedMulti(TLimitsString const &limits)
        {
            typename Iterator<TLimitsString const, Standard>::Type it = begin(limits, Standard()), itEnd = end(limits, Standard());

            n4 = n2 = n1 = n24 = 0;

            if (it == itEnd) return;

            // count the numbers of septets in residue class 1, 2, and 4

            TSize size, cur;
            TSize old = *it; ++it;

            while (it != itEnd) {
                cur = *it;
                size = cur - old;
                old = cur;

                n4 += (size + 3) / 7;
                n2 += (size + 5) / 7;
                n1 += (size + 6) / 7;

                ++it;
            }

            n24 = n2 + n4;

            // precompute the begin positions (off) of septet names
            // in the sliced string for every sequence and every residue class

            resize(off[1], length(limits) - 1);
            resize(off[2], length(limits) - 1);
            resize(off[4], length(limits) - 1);

            it = begin(limits, Standard());
            old = *it; ++it;

            TSize seqNo = 0;
            TSize off4 = 0, off2 = n4, off1 = n24;
            while (it != itEnd) {
                cur = *it;
                size = cur - old;
                old = cur;

                off4 += (size + 3) / 7;
                off2 += (size + 5) / 7;
                off1 += (size + 6) / 7;

                off[4][seqNo] = (off4)? off4-1 : 0;
                off[2][seqNo] = (off2)? off2-1 : 0;
                off[1][seqNo] = (off1)? off1-1 : 0;

                ++it; ++seqNo;
            }
        }

        inline TResult operator() (const TValue &x) const {
            typename Value<TValue,1>::Type x_i1 = x.i1;
            TSize seqOfs = getValueI2(x_i1);
            return TResult(off[seqOfs % 7][getValueI1(x_i1)] - seqOfs/7, x.i2);
        }
    };



    //////////////////////////////////////////////////////////////////////////////
    // skew7 class
    template < typename TInput, typename TPair_, typename TLimitsString >
    struct Pipe< TInput, Multi<Skew7, TPair_, TLimitsString> >
    {

        // *** SPECIALIZATION ***

        // use packing if lessorequal 16 different values per char
        typedef typename MakePacked<TPair_>::Type TPair;
        typedef typename IfC<
            (BitsPerValue<TypeOf_(TInput)>::VALUE > 0) &&
            (BitsPerValue<TypeOf_(TInput)>::VALUE <= 4),
            BitPacked<>,
            Pack>::Type TPack;

        // step 1
        typedef Pipe< TInput, Multi<Sampler<7, TPack>, TPair, TLimitsString> >  TSamplerDC7;
                                        typedef _skew7NCompMulti<TypeOf_(TSamplerDC7)> ncomp_t;
        typedef Pool< TypeOf_(TSamplerDC7), SorterSpec< SorterConfigSize<ncomp_t, TSizeOf_(TSamplerDC7) > > > TSortTuples;
        typedef Pipe< TSortTuples, Namer<ncomp_t> > TNamer;
                                        typedef _skew7GlobalSlicedMulti<TypeOf_(TNamer), TLimitsString, TSizeOf_(TNamer)> func_slice_t;

        typedef Pipe< TNamer, Filter<func_slice_t> > TSlicedPos;
        typedef Pool< TypeOf_(TSlicedPos), MapperSpec< MapperConfigSize< filterI1<TypeOf_(TSlicedPos)>, TSizeOf_(TSlicedPos) > > > TNames_Sliced;

        // non-unique names
        typedef Pipe< TNames_Sliced, Filter< filterI2<TypeOf_(TNames_Sliced)> > > TFilter;

            // recursion
            typedef Pipe< TFilter, Skew7 > TRecurse;
            typedef Pipe< TRecurse, Counter > TRenamer;

            // no recursion inMemory shortcut
            typedef Pipe< TFilter, LarssonSadakane > TInMem;
            typedef Pipe< TInMem, Counter > TRenamerInMem;

        typedef Pool< TypeOf_(TRenamer), MapperSpec< MapperConfigSize< filterI1<TypeOf_(TRenamer)>, TSizeOf_(TRenamer) > > > TNames_Linear;

        // step 2
        typedef Pipe< Bundle2< TInput, TNames_Linear >, Extender7Multi<TPair, TPack> > TExtender;
                                        typedef _skew7ExtendComp<TypeOf_(typename TExtender::Out0),3> extend0_comp_t;
                                        typedef _skew7ExtendComp<TypeOf_(typename TExtender::Out6),2> extend6_comp_t;
                                        typedef _skew7ExtendComp<TypeOf_(typename TExtender::Out5),1> extend5_comp_t;
                                        typedef _skew7ExtendComp<TypeOf_(typename TExtender::Out3),1> extend3_comp_t;
        typedef Pool< TypeOf_(typename TExtender::Out0), SorterSpec< SorterConfigSize< extend0_comp_t, TSizeOf_(typename TExtender::Out0) > > > TSorterS0;
        typedef Pool< TypeOf_(typename TExtender::Out6), SorterSpec< SorterConfigSize< extend6_comp_t, TSizeOf_(typename TExtender::Out6) > > > TSorterS6;
        typedef Pool< TypeOf_(typename TExtender::Out5), SorterSpec< SorterConfigSize< extend5_comp_t, TSizeOf_(typename TExtender::Out5) > > > TSorterS5;
        typedef Pool< TypeOf_(typename TExtender::Out3), SorterSpec< SorterConfigSize< extend3_comp_t, TSizeOf_(typename TExtender::Out3) > > > TSorterS3;

        // step 3
                                        typedef _skew7NMapExtended<TypeOf_(typename TExtender::Out124)> nmap_extended_t;
        typedef Pool< TypeOf_(typename TExtender::Out124), MapperSpec< MapperConfigSize< nmap_extended_t, TSizeOf_(typename TExtender::Out124) > > > TSorterS124;
        typedef Pipe< Bundle5< TSorterS0, TSorterS3, TSorterS5, TSorterS6, TSorterS124 >, Merger7Multi<TLimitsString> > TMerger;

        TSorterS0            sortedS0;
        TSorterS3            sortedS3;
        TSorterS5            sortedS5;
        TSorterS6            sortedS6;
        TSorterS124            sortedS124;
        TMerger                in;
        TLimitsString       const &limits;

        // (weese): The SEQAN_CTOR_ENABLE_IF is necessary to avoid implicit casts and references to temporaries

        template <typename TLimitsString_>
        Pipe(TLimitsString_ const &_limits, SEQAN_CTOR_ENABLE_IF(IsSameType<TLimitsString, TLimitsString_>)) :
            in(bundle5(sortedS0, sortedS3, sortedS5, sortedS6, sortedS124), _limits),
            limits(_limits)
        {
            ignoreUnusedVariableWarning(dummy);
        }

        template <typename TLimitsString_>
        Pipe(TInput& _textIn, TLimitsString_ const &_limits, SEQAN_CTOR_ENABLE_IF(IsSameType<TLimitsString, TLimitsString_>)) :
            in(bundle5(sortedS0, sortedS3, sortedS5, sortedS6, sortedS124), _limits),
            limits(_limits)
        {
            ignoreUnusedVariableWarning(dummy);
            process(_textIn);
        }

        template < typename TInput_ >
        bool process(TInput_ &textIn) {

            SEQAN_PROADD(SEQAN_PRODEPTH, 1);
            SEQAN_PROMARK("Enter recursion");
            #ifdef SEQAN_DEBUG_INDEX
                std::cerr << "enter level " << SEQAN_PROVAL(SEQAN_PRODEPTH) << " bit-packing: ";
                std::cerr << IsSameType<TPack, BitPacked<> >::VALUE << " "<<ValueSize<TypeOf_(TInput)>::VALUE<<"" << std::endl;
            #endif
            {


            // *** INSTANTIATION ***

            // step 1
            TSamplerDC7                 sampler(textIn, limits);
            TSortTuples                 sorter;
            #ifdef SEQAN_DEBUG_INDEX
                std::cerr << "  sort names (" << length(sampler)<< ")" << std::endl;
            #endif
            sorter << sampler;
            SEQAN_PROMARK("Sorter (2) - sort 7-mers");

            TNamer                      namer(sorter);
            func_slice_t                func_slice(limits);

            TSlicedPos                    slicedPos(namer, func_slice);
            TNames_Sliced               names_sliced;
            #ifdef SEQAN_DEBUG_INDEX
                std::cerr << "  slice names" << std::endl;
            #endif
            names_sliced << slicedPos;

            if (namer.unique() || empty(names_sliced)) {
                // unique names

                clear(sorter);
                SEQAN_PROMARK("Mapper (4) - construct s124");

                TNames_Linear               names_S1, names_S2, names_S4;

                #ifdef SEQAN_DEBUG_INDEX
                    std::cerr << "  make names linear" << std::endl;
                #endif

                _skew7SeparateSlices(
                    names_sliced, func_slice,
                    names_S1, names_S2, names_S4);

                clear(names_sliced);
                SEQAN_PROMARK("Mapper (10) - construct ISA124");

                // step 2
                _skew7ExtendMulti(
                    textIn, limits,
                    names_S1, names_S2, names_S4,
                    sortedS0, sortedS3, sortedS5, sortedS6, sortedS124);

            } else {
                // non-unique names

                clear(sorter);
                SEQAN_PROMARK("Mapper (4) - construct s124");

                TFilter                     filter(names_sliced);
                TNames_Linear               names_S1, names_S2, names_S4;

//                if (length(filter) > 128*1024*1024)
                {
                    // recursion
                    TRecurse                    recurse(filter);

                    #ifdef SEQAN_TEST_SKEW7
                    {
                        String<typename Value<TFilter>::Type, Alloc<> > _text;
                        _text << filter;
                        SEQAN_ASSERT(isSuffixArray(recurse, _text));
                    }
                    #endif

                    clear(filter);
                    TRenamer                    renamer(recurse);

                    // partition SA by residue classes

                    #ifdef SEQAN_DEBUG_INDEX
                        std::cerr << "  rename names" << std::endl;
                    #endif

                    _skew7SeparateSlices(
                        renamer, func_slice,
                        names_S1, names_S2, names_S4);
                }
/*                else
                {
                    TInMem                        inMem(filter);
                    clear(filter);
                    TRenamerInMem               renamer(inMem);
                    _skew7SeparateSlices(
                        renamer, func_slice,
                        names_S1, names_S2, names_S4);
                }
*/
                SEQAN_PROMARK("Mapper (10) - construct ISA124");

                // step 2
                #ifdef SEQAN_DEBUG_INDEX
                    std::cerr << "  prepare merge" << std::endl;
                #endif
                _skew7ExtendMulti(
                    textIn, limits,
                    names_S1, names_S2, names_S4,
                    sortedS0, sortedS3, sortedS5, sortedS6, sortedS124);

                SEQAN_PROMARK("Mapper (12), Sorter (13-16) - merge SA124, SA3, SA5, SA6, SA0");
            }

            // step 3
            // ... is done on-demand by merger
            }
            #ifdef SEQAN_DEBUG_INDEX
                std::cerr << "left level " << SEQAN_PROVAL(SEQAN_PRODEPTH) << std::endl;
            #endif
            SEQAN_PROMARK("Left recursion");
            SEQAN_PROSUB(SEQAN_PRODEPTH, 1);

            return true;
        }

        inline typename Value<Pipe>::Type const operator*() {
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
    template < typename TInput, typename TObject, typename TPair, typename TLimitsString >
    inline bool operator<<(Pipe< TInput, Multi<Skew7, TPair, TLimitsString> > &me, TObject &textIn) {
        return me.process(textIn);
    }

}

#endif
