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

#ifndef SEQAN_HEADER_INDEX_SKEW7_H
#define SEQAN_HEADER_INDEX_SKEW7_H

namespace SEQAN_NAMESPACE_MAIN
{

//namespace SEQAN_NAMESPACE_PIPELINING
//{

    struct Skew7 {};


    //////////////////////////////////////////////////////////////////////////////
    // external Skew7 algorithm
    //////////////////////////////////////////////////////////////////////////////

    template <typename T>
    struct SkewDC_<7, T>
    {
        static const unsigned VALUE[];
    };

    template <typename T>
    const unsigned SkewDC_<7, T>::VALUE[] = { 3,   1, 2, 4 };


    // *** COMPARATORS & MAPS ***

    template <typename TValue, typename TResult = int>
    struct _skew7NComp :
        public std::binary_function<TValue, TValue, TResult>
    {
        inline TResult operator() (const TValue &a, const TValue &b) const
        {
            typedef typename Value<TValue, 1>::Type                 TSize;
            typedef typename Value<TValue, 2>::Type                 TSeptet;
            typedef typename Value<TSeptet>::Type                   TSeptetValue;
            typedef typename StoredTupleValue_<TSeptetValue>::Type  TStoredValue;

            const TStoredValue *sa = a.i2.i;
            const TStoredValue *sb = b.i2.i;

            TSize n = LENGTH<TSeptet>::VALUE;
            if (a.i1 < n) n = a.i1;
            if (b.i1 < n) n = b.i1;

            // compare the overlap of septets (the first n bases)
            for (TSize i = 0; i < n; i++, ++sa, ++sb)
            {
                if (*sa == *sb) continue;
                return (*sa < *sb)? -1 : 1;
            }

            // overlap is equal, now suffix lengths decide
            if (n < LENGTH<TSeptet>::VALUE)
                return (a.i1 < b.i1)? -1 : 1;
            else
                return 0;
        }
    };

    // optimized for bitvectors
    template <typename T1, typename TValue, typename TResult>
    struct _skew7NComp< Pair<T1, Tuple<TValue, 7, BitPacked<> >, Pack >, TResult > :
        public std::binary_function<
            Pair<T1, Tuple<TValue, 7, BitPacked<> >, Pack >,
            Pair<T1, Tuple<TValue, 7, BitPacked<> >, Pack >,
            TResult>
    {
        inline TResult operator() (
            const Pair<T1, Tuple<TValue, 7, BitPacked<> >, Pack > &a,
            const Pair<T1, Tuple<TValue, 7, BitPacked<> >, Pack > &b) const
        {
            if (a.i2 < b.i2) return -1;
            if (a.i2 > b.i2) return 1;
            if (a.i1 < 7 || b.i1 < 7)
                return (a.i1 < b.i1)? -1 : 1;
            return 0;
        }
    };

    template <typename TValue, typename TResult = typename Value<TValue, 1>::Type>
    struct _skew7NMapLinear :
        public std::unary_function<TValue, TResult>
    {
        TResult BN4, BN;

        _skew7NMapLinear(TResult BN_): BN4(BN_+1), BN(BN_) {}

        inline TResult operator() (const TValue& x) const
        {
            TResult i = x.i1; return (i%7 == 4)? BN4-(i-(i/7)*4): BN-(i-(i/7)*4);
        }
    };

    template <typename TValue, typename TResult = typename Value<TValue, 1>::Type>
    struct _skew7NMapSliced :
        public std::unary_function<TValue, TResult>
    {
        TResult off[5];

        _skew7NMapSliced(TResult BN_)
        {
            off[0] = 0;
            off[1] = BN_ - 1;
            off[2] = (2*BN_)/3 - 1;
            off[3] = 0;
            off[4] = BN_/3 - 1;
        }

        inline TResult operator() (const TValue& x) const
        {
            return off[x.i1 % 7] - x.i1/7;
        }
    };


    template <typename TValue, typename TResult = TValue>
    struct _skew7UnslicerFunc :
        public std::unary_function<TValue, TResult>
    {
        TResult o1, o2, o4, n4, n24;

        _skew7UnslicerFunc(TResult N) :
            o1(N - (N + 6) % 7),
            o2(N - (N + 5) % 7),
            o4(N - (N + 3) % 7),
            n4((N + 3) / 7),
            n24((N + 5) / 7 + n4) {}

        inline TResult operator() (const TValue& x) const
        {
            return (x < n4)  ? o4 -  x        * 7:
                   (x < n24) ? o2 - (x - n4)  * 7:
                               o1 - (x - n24) * 7;
        }
    };

    template <typename TValue, typename TResult = typename Value<typename Value<TValue, 2>::Type>::Type>
    struct _skew7NMapExtended :
        public std::unary_function<TValue, TResult>
    {
        inline TResult operator() (const TValue& x) const
        {
            return x.i2[0];
        }
    };

    template <typename TValue, unsigned EXT_LENGTH, typename TResult = int>
    struct _skew7ExtendComp :
        public std::binary_function<TValue, TValue, TResult>
    {
        inline TResult operator() (const TValue &a, const TValue &b) const
        {
            for (unsigned i = 0; i < EXT_LENGTH; i++)
            {
                if (a.i3[i] == b.i3[i]) continue;
                return (a.i3[i] <  b.i3[i])? -1 : 1;
            }
            return (a.i2[0] < b.i2[0])? -1 : 1;
        }
    };

    // optimized for bitvectors
    /*
    template <typename T1, typename T2, typename T, const int _size, const int EXT_LENGTH, typename TResult>
    struct _skew7ExtendComp< Triple<T1,T2,Tuple<T,_size, BitPacked<> >, Pack>, EXT_LENGTH, TResult> :
        public std::binary_function<
            Triple<T1,T2,Tuple<T,_size, BitPacked<> >, Pack>,
            Triple<T1,T2,Tuple<T,_size, BitPacked<> >, Pack>,
            TResult>
    {
        inline TResult operator() (
            const Triple<T1,T2,Tuple<T,_size, BitPacked<> >, Pack> &a,
            const Triple<T1,T2,Tuple<T,_size, BitPacked<> >, Pack> &b) const
        {
            if (a.i3 < b.i3) return -1;
            if (a.i3 > b.i3) return 1;
            return (a.i2[0] < b.i2[0])? -1 : 1;
        }
    };
    */

    template < typename TInput >
    struct Value< Pipe< TInput, Skew7 > > :
        public Size<TInput> {};

    template <typename T>
    struct Skew7StringSpec_:
        public Spec<T> {};

    template <typename T, typename TStringSpec>
    struct Skew7StringSpec_<String<T, TStringSpec> >
    {
        typedef TStringSpec Type;
    };

    template <typename T, typename TSegmentSpec>
    struct Skew7StringSpec_<Segment<T, TSegmentSpec> >:
        public Skew7StringSpec_<T> {};

    template <typename T>
    struct Skew7StringSpec_<T const>:
        public Skew7StringSpec_<T> {};

    //////////////////////////////////////////////////////////////////////////////
    // skew7 class
    template < typename TInput >
    struct Pipe< TInput, Skew7 >
    {

        // *** SPECIALIZATION ***

        // use packing if lessorequal 16 different values per char
        typedef typename IfC<
            (BitsPerValue<TypeOf_(TInput)>::VALUE > 0) &&
            (BitsPerValue<TypeOf_(TInput)>::VALUE <= 4),
            BitPacked<>,
            Pack>::Type TPack;

        // step 1
        typedef Pipe< TInput, Sampler<7, TPack> >  TSamplerDC7;
                                        typedef _skew7NComp<TypeOf_(TSamplerDC7)> ncomp_t;
        typedef Pool< TypeOf_(TSamplerDC7), SorterSpec< SorterConfigSize<ncomp_t, TSizeOf_(TSamplerDC7) > > > TSortTuples;
        typedef Pipe< TSortTuples, Namer<ncomp_t> > TNamer;
                                        typedef _skew7NMapSliced<TypeOf_(TNamer)> nmap_sliced_t;
                                        typedef _skew7NMapLinear<TypeOf_(TNamer)> nmap_linear_t;
        typedef Pool< TypeOf_(TNamer), MapperSpec< MapperConfigSize< nmap_sliced_t, TSizeOf_(TNamer) > > > TNames_Sliced;

        // unique names - shortcut
        typedef Pool< TypeOf_(TNames_Sliced), MapperSpec< MapperConfigSize< nmap_linear_t, TSizeOf_(TNames_Sliced) > > > TNames_Linear_Unique;

        // non-unique names
        typedef Pipe< TNames_Sliced, Filter< filterI2<TypeOf_(TNames_Sliced)> > > TFilter;

            // recursion
            typedef Pipe< TFilter, Skew3 > TRecurse;
                                        typedef _skew7UnslicerFunc<TypeOf_(TRecurse)> unslicer_func_t;
            typedef Pipe< TRecurse, Filter<unslicer_func_t> > TUnslicer;
            typedef Pipe< TUnslicer, Counter > TRenamer;

            // no recursion inMemory shortcut
            typedef Pipe< TFilter, LarssonSadakane > TInMem;
            typedef Pipe< TInMem, Filter<unslicer_func_t> > TUnslicerInMem;
            typedef Pipe< TUnslicerInMem, Counter > TRenamerInMem;

        typedef Pool< TypeOf_(TRenamer), MapperSpec< MapperConfigSize< nmap_linear_t, TSizeOf_(TRenamer) > > > TNames_Linear;

        // step 2
        typedef Pipe< Bundle2< TInput, TNames_Linear >, Extender7<TPack> > TExtender;
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
        typedef Pipe< Bundle5< TSorterS0, TSorterS3, TSorterS5, TSorterS6, TSorterS124 >, Merger7 > TMerger;

        TSorterS0   sortedS0;
        TSorterS3   sortedS3;
        TSorterS5   sortedS5;
        TSorterS6   sortedS6;
        TSorterS124 sortedS124;
        TMerger     in;

        Pipe():
            in(bundle5(sortedS0, sortedS3, sortedS5, sortedS6, sortedS124)) {}

        Pipe(TInput& _textIn):
            in(bundle5(sortedS0, sortedS3, sortedS5, sortedS6, sortedS124))
        {
            process(_textIn);
        }

        template < typename TInput_ >
        bool process(TInput_ &textIn) {

            SEQAN_PROADD(SEQAN_PRODEPTH, 1);
            SEQAN_PROMARK("Enter recursion");
            #ifdef SEQAN_DEBUG_INDEX
                std::cerr << "enter level " << SEQAN_PROVAL(SEQAN_PRODEPTH) << " bit-packing: ";
                std::cerr << IsSameType<TPack, BitPacked<> >::VALUE << " " << BitsPerValue<TypeOf_(TInput)>::VALUE << std::endl;
            #endif
            {


            // *** INSTANTIATION ***

            // step 1
            TSamplerDC7                 sampler(textIn);
            TSortTuples                 sorter;
            #ifdef SEQAN_DEBUG_INDEX
                std::cerr << "  sort names (" << length(sampler)<< ")" << std::endl;
            #endif
            sorter << sampler;
            SEQAN_PROMARK("Sorter (2) - sort 7-mers");

            TNamer                      namer(sorter);
            nmap_sliced_t               map_sliced(length(namer));
            nmap_linear_t               map_linear(length(namer));
            TNames_Sliced               names_sliced(map_sliced);
            #ifdef SEQAN_DEBUG_INDEX
                std::cerr << "  slice names" << std::endl;
            #endif
            names_sliced << namer;

            if (namer.unique() || empty(names_sliced)) {
                // unique names

                clear(sorter);
                SEQAN_PROMARK("Mapper (4) - construct s124");
                TNames_Linear_Unique        names_linear(map_linear);

                #ifdef SEQAN_DEBUG_INDEX
                    std::cerr << "  make names linear" << std::endl;
                #endif
                names_linear << names_sliced;
                clear(names_sliced);
                SEQAN_PROMARK("Mapper (10) - construct ISA124");

                // step 2
                _skew7Extend(textIn, names_linear, sortedS0, sortedS3, sortedS5, sortedS6, sortedS124);

            } else {
                // non-unique names

                clear(sorter);
                SEQAN_PROMARK("Mapper (4) - construct s124");

                TFilter                     filter(names_sliced);
                TNames_Linear               names_linear(map_linear);

                if (length(filter) > 128*1024*1024)
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
                    unslicer_func_t             func(length(textIn));
                    TUnslicer                   unslicer(recurse, func);
                    TRenamer                    renamer(unslicer);

                    #ifdef SEQAN_DEBUG_INDEX
                        std::cerr << "  rename names" << std::endl;
                    #endif

                    names_linear << renamer;
                }
                else
                {
                    TInMem                        inMem(filter);

                    clear(filter);
                    unslicer_func_t                func(length(textIn));
                    TUnslicerInMem              unslicer(inMem, func);
                    TRenamerInMem               renamer(unslicer);

                    #ifdef SEQAN_DEBUG_INDEX
                        std::cerr << "  rename names" << std::endl;
                    #endif

                    names_linear << renamer;
                }

                SEQAN_PROMARK("Mapper (10) - ISA124 konstruieren");

                // step 2
                #ifdef SEQAN_DEBUG_INDEX
                    std::cerr << "  prepare merge" << std::endl;
                #endif
                _skew7Extend(textIn, names_linear, sortedS0, sortedS3, sortedS5, sortedS6, sortedS124);
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
    template < typename TInput, typename TObject >
    inline bool operator<<(Pipe< TInput, Skew7 > &me, TObject &textIn) {
        return me.process(textIn);
    }


    //////////////////////////////////////////////////////////////////////////////
    // internal Skew7 algorithm
    //////////////////////////////////////////////////////////////////////////////

    //////////////////////////////////////////////////////////////////////////////
    // typedefs and helpers

    // compares n characters and in case of equality the names a2 and b2 (no clipping)
    template <typename TTextIter, typename TSize> inline
    bool _leqSkew7(TTextIter a1, TSize a2,   TTextIter b1, TSize b2,   TSize n)
    {
        // lexic. order for n-tupels
        for (; n != 0; --n, ++a1, ++b1) {
            if (ordLess(*a1, *b1)) return true;
            if (ordLess(*b1, *a1)) return false;
        }
        return (a2 <= b2);
    }

    // compares at most the last n characters (a) with b (clipping)
    template <typename TTextIter, typename TSize> inline
    bool _leqSkew7(TTextIter a,   TTextIter b,   TSize n)
    { // lexic. order for n-tupels
        for (; n != 0; --n, ++a, ++b) {
            if (ordLess(*a, *b)) return true;
            if (ordLess(*b, *a)) return false;
        }
        return true;    // a is shorter than b
    }

    // compares two suffixes of residue classes a and b
    template <typename TTextIter, typename TSize, typename TString> inline
    bool _leqSkew7(unsigned a, unsigned b,   TTextIter spos[], const TSize tpos[], const bool islast[], const TString &s124, const TSize adjust[7][7])
    {
        TTextIter sa = spos[a];
        TTextIter sb = spos[b];
        TSize shft = SkewShift_<7>::VALUE[a][b];
        if (sa > sb) {
            if ((a != 0) && (a < shft) && islast[a]) // do we need to clip?
                return _leqSkew7 (sa,   sb,   a);
        } else {
            if ((b != 0) && (b < shft) && islast[b]) // do we need to clip?
                return !_leqSkew7 (sb,   sa,   b);
        }
        return _leqSkew7 (sa, s124[tpos[a] + adjust[a][b]],   sb, s124[tpos[b] + adjust[b][a]],   shft);
    }

    template <typename TSAOut, typename TText, typename TSAIn, typename TEnable>
    inline bool _fastTupleSortSkew7(TSAOut &, TText const &, TSAIn const &, TEnable)
    {
        return false;
    }

    template <typename TSAOut, typename TText, typename TSAIn>
    inline bool _fastTupleSortSkew7(TSAOut &SA124, TText const &s, TSAIn const &s124, True)
    {
        typedef typename Value<TText>::Type TValue;
        typedef typename Value<TSAOut>::Type TSize;
        typedef typename Iterator<TText const, Standard>::Type TValueIter;
        typedef typename Iterator<TSAIn const, Standard>::Type TSAIter;

        // optimized tuple sort for short alphabets
        Shape<TValue, UngappedShape<7> > shape;
        String<TSize, Alloc<> > cnt;
        resize(cnt, _fullDir2Length(shape), 0, Exact());

        TValueIter textBegin = begin(s, Standard());
        TSAIter it = begin(s124, Standard());
        TSAIter itEnd = end(s124, Standard());

        TSize n = length(s);
        for(; it != itEnd; ++it)
            ++cnt[hash2(shape, textBegin + *it, n - *it)];

        _qgramCummulativeSum(cnt, False());

        it = begin(s124, Standard());
        for(; it != itEnd; ++it)
            SA124[cnt[hash2(shape, textBegin + *it, n - *it) + 1]++] = *it;

        return true;
    }


    //////////////////////////////////////////////////////////////////////////////
    // Skew algorithm with difference cover of Z_7
    //
    // Construct the suffix array SA of s[0..n-1], where s is a string over
    // the alphabet {0..K}.
    //
    // The following algorithm divides the suffixes in seven residue classes
    // according to their lengths and uses a difference cover {1,2,4}.
    // That approach results in an algorithm that is more space and time efficient
    // than the original skew algorithm with a difference cover {1,2} of Z_3.
    //
    // * no trailing 0's required
    // * no dummy triples in special cases

    template < typename TSA,
               typename TText >
    void createSuffixArray(
        TSA &SA,
        TText &s,
        Skew7 const &,
        unsigned K,
        unsigned maxdepth,
        unsigned depth)
    {
        typedef typename Value<TSA>::Type TSize;
        typedef typename Value<TText>::Type TValue;
        typedef String<TSize, Alloc<> > TBuffer;
        typedef typename Iterator<TSA, Standard>::Type TSAIter;
        //typedef typename Iterator<TBuffer, Standard>::Type TBufferIter;
        typedef typename Iterator<TText const, Standard>::Type TValueIter;

        SEQAN_ASSERT(AllowsFastRandomAccess<TText>::VALUE == true);
        SEQAN_ASSERT(AllowsFastRandomAccess<TSA>::VALUE == true);

        #ifdef SEQAN_DEBUG_INDEX
            std::cerr << "--- CREATE SUFFIX ARRAY ---" << std::endl;
            std::cerr << "Skew7 [random access]" << std::endl;
        #endif

        TSize n = length(s);
        if (n < 1) return;

        TSize _n[7];
        TSize _o[7];

        _n[0] = n/7;
        _o[0] = n%7;
        TSize j = n + 6;
        for(int i = 1; i < 7; ++i, --j) {
            _n[i] = j/7;
            _o[i] = j%7;
        }

        TSize _n24  = _n[2]+_n[4];
        TSize _n124 = _n[1]+_n24;

        SEQAN_PROSET(SEQAN_PRODEPTH, depth);
        SEQAN_PROSET(SEQAN_PROEXTRA1, K);
        SEQAN_PROMARK("Rekursionsabstieg");
        #ifdef SEQAN_DEBUG_INDEX
            std::cerr << "enter level " << depth << " (" << n << ")" << std::endl;
        #endif

        TBuffer s124;
        resize(s124, _n124, Exact());
        // we use SA[n-n124..n-1] as a temporary buffer instead of allocating one
        typename Suffix<TSA>::Type SA124 = suffix(SA, n - _n124);


        // generate positions of mod 3, mod 5 and mod 6 suffixes
        {
            TSize j = 0;
            if (_n[2] > _n[4]) s124[j++] = _o[2];
            if (_n[1] > _n[4]) s124[j++] = _o[1];

            for (TSize i=_o[4];  i < n;  i+=7) {
                s124[j++] = i;
                s124[j++] = i + 2;
                s124[j++] = i + 3;
            }
        }


        // lsb radix sort the mod 3, mod 5 and mod 6 7-tupels
        if (!_fastTupleSortSkew7(SA124, s, s124, typename Eval<BitsPerValue<TValue>::VALUE < 4>::Type()))
        {
            String<TSize, Alloc<> > cnt;
            resize(cnt, K, Exact());    // counter array

            radixPass(SA124, s124,  s, cnt, K, 6);
            radixPass(s124,  SA124, s, cnt, K, 5);
            radixPass(SA124, s124,  s, cnt, K, 4);
            radixPass(s124,  SA124, s, cnt, K, 3);
            radixPass(SA124, s124,  s, cnt, K, 2);
            radixPass(s124,  SA124, s, cnt, K, 1);
            radixPass(SA124, s124,  s, cnt, K);
        }
        SEQAN_PROMARK("7-lets sortiert");

        // find lexicographic names of 7-tupel
        TSize name = 0;
        {
            TSize ofs[7] = {0, _n24, _n[4], 0, 0, 0, 0};
            bool differ = true;
            TValue c0 = TValue(0), c1 = TValue(0), c2 = TValue(0), c3 = TValue(0), c4 = TValue(0), c5 = TValue(0), c6 = TValue(0);
            for (TSize i = 0, clip = _max(n, (TSize)6) - 6, l;  i < _n124;  i++) {
                if ((l = SA124[i]) < clip) {
                    if (differ || s[l] != c0 || s[l+1] != c1 || s[l+2] != c2 || s[l+3] != c3 ||
                                                s[l+4] != c4 || s[l+5] != c5 || s[l+6] != c6) {
                        name++;  c0 = s[l];  c1 = s[l+1];  c2 = s[l+2];  c3 = s[l+3];
                                             c4 = s[l+4];  c5 = s[l+5];  c6 = s[l+6];
                        differ = false;
                    }
                } else {
                    name++;
                    differ = true;  // the last 6 7-tupels always differ from the rest
                }
                s124[(l/7) + ofs[(n-l) % 7]] = name - 1;   // select a third
            }
        }
        SEQAN_PROMARK("s12 konstruiert");

        // recurse if names are not yet unique
        if (name < _n124) {
            if (depth != maxdepth)
            {
                createSuffixArray(SA124, s124, Skew7(), name, maxdepth, depth + 1);
                #ifdef SEQAN_TEST_SKEW7
                    SEQAN_ASSERT(isSuffixArray(SA124, s124));
                #endif
            }
            // store unique names in s124 using the suffix array
            for (TSize i = 0;  i < _n124;  i++) s124[SA124[i]] = i;
        } else // generate the suffix array of s124 directly
            for (TSize i = 0;  i < _n124;  i++) SA124[s124[i]] = i;


        // use SA[0...n3-1] and SA[n3...n3+n5-1] as a temporary buffers instead of allocating some
        // and allocate SA0, SA3, SA5 and SA6

        {
            typename Infix<TSA>::Type s3 = infix(SA, 0, _n[3]), s5 = infix(SA, _n[3], _n[3] + _n[5]);
            String<TSize, typename Skew7StringSpec_<TSA>::Type> SA0, SA3, SA5, SA6;

            resize(SA0, _n[0], Exact());
            resize(SA3, _n[3], Exact());
            resize(SA5, _n[5], Exact());
            resize(SA6, _n[6], Exact());

            // stably sort the mod 5 and mod 3 suffixes from SA124 by their first character
            {
                for (TSize i=0, j3=0, j5=0, l;  i < _n124;  i++) {
                    l = SA124[i];
                    if (l < _n[4]) {
                        if ((l = _o[4] + (7 * l)) > 0)
                            s5[j5++] = l - 1;
                    } else if (l < _n24) {
                        if ((l = _o[2] + (7 * (l - _n[4]))) > 0)
                            s3[j3++] = l - 1;
                    }
                }

                {
                    String<TSize, Alloc<> > cnt;
                    resize(cnt, K, Exact());    // counter array

                    radixPass(SA3, s3, s, cnt, K);
                    SEQAN_PROMARK("SA3 konstruiert");

                    radixPass(SA5, s5, s, cnt, K);
                    SEQAN_PROMARK("SA5 konstruiert");

                    // stably sort the mod 6 suffixes from SA5 by their first character

                    if (_n[5] == _n[6]) radixExtend    (SA6, SA5, s, cnt, K);
                    else                radixExtendClip(SA6, SA5, s, cnt, K);
                    SEQAN_PROMARK("SA6 konstruiert");

                    // stably sort the mod 0 suffixes from SA6 by their first character

                    if (_n[6] == _n[0]) radixExtend    (SA0, SA6, s, cnt, K);
                    else                radixExtendClip(SA0, SA6, s, cnt, K);
                    SEQAN_PROMARK("SA0 konstruiert");
                }
            }

            // MULTIWAY MERGE all SA_ streams
            {
                // a helper matrix to lex-name-compare every combination of suffixes
                TSize adjust[7][7] =
                    //      0               1              2             3             4              5               6
                   {{0             , _n124-_n[0]   , _n24-_n[0]  , _n124-_n[0] , _n[4]-_n[0] , _n[4]-_n[0]   , _n24-_n[0]    },  // 0
                    {1-_n[1]       , 0             , 0           , 1-_n[1]     , 0           , 1-_n[1]-_n[2] , 1-_n[1]-_n[2] },  // 1*
                    {1-_n[2]       , 0             , 0           , _n[1]       , 0           , _n[1]         , 1-_n[2]       },  // 2*
                    {1+_n[4]-_n[3] , 1+_n[4]-_n[3] , _n24-_n[3]  , 0           , _n124-_n[3] , _n24-_n[3]    , _n124-_n[3]   },  // 3
                    {_n[1]+_n[2]   , 0             , 0           ,_n[2]        , 0           , _n[1]+_n[2]   , _n[2]         },  // 4*
                    {_n24-_n[5]    , _n124-_n[5]   , _n[4]-_n[5] , _n[4]-_n[5] , _n24-_n[5]  , 0             , _n124-_n[5]   },  // 5
                    {_n124-_n[6]   , _n24-_n[6]    , _n124-_n[6] , _n[4]-_n[6] , _n[4]-_n[6] , _n24-_n[6]    , 0             }}; // 6

                TSAIter pos[7] = {begin(SA0, Standard())      , begin(SA124, Standard())      , begin(SA124, Standard())      ,
                                  begin(SA3, Standard())      , begin(SA124, Standard())      , begin(SA5, Standard())        , begin(SA6, Standard())      };
                TSAIter max[7] = {begin(SA0, Standard())+_n[0], begin(SA124, Standard())+_n124, begin(SA124, Standard())+_n124,
                                  begin(SA3, Standard())+_n[3], begin(SA124, Standard())+_n124, begin(SA5, Standard())+_n[5]  , begin(SA6, Standard())+_n[6]};
                TValueIter spos[7];
                TSize tpos[7];
                bool islast[7];

                int a, b, rank[5];
                int fill = 0;
                TSize k = 0;

                #define SEQAN_GET_ISKEW7(ii) (ii < _n[4] ? (ii * 7) + _o[4] : (ii < _n24 ? ((ii - _n[4]) * 7) + _o[2] : ((ii - _n24) * 7) + _o[1]))
                #define SEQAN_GET_ASKEW7(ii) (ii < _n[4] ? 4 : (ii < _n24 ? 2 : 1))

                // fill the stream ranking list
                for(int i = 0; i < 7; ++i) {
                    if (!_n[i]) continue;
                    if (i == 2 || i == 4) continue; // insert only the least suffix of SA124

                    if (i == 1) {

                        TSize ii  = *(pos[1]);
                        a         = SEQAN_GET_ASKEW7(ii);
                        TSize j      = SEQAN_GET_ISKEW7(ii);

                        tpos[a]   = ii;
                        spos[a]   = begin(s, Standard()) + j;
                        islast[a] = (j + 7 >= n);

                    } else {

                        a = i;
                        TSize j   = *(pos[a]);

                        tpos[a]   = j / 7;
                        spos[a]   = begin(s, Standard()) + j;
                        islast[a] = (j + 7 >= n);

                    }

                    // get the rank of stream a's suffix
                    int j;
                    for(j = 0;  j < fill;  ++j) {
                        b = rank[j];
                        if (_leqSkew7 (a,  b,   spos, tpos, islast, s124, adjust)) break;
                    }

                    // insert the suffix
                    for(int i = fill; i > j; --i)
                        rank[i] = rank[i-1];
                    rank[j] = a;
                    fill++;
                }

                // main merge loop
                // in order to find the least suffix in every step we use a stream ranking list
                // so we only need to keep up the ordering and thus rank[0] is always the least suffix
                while (fill > 1) {

                    // add the least suffix to SA and get the next of the corresponding stream
                    a = rank[0];
                    SA[k++] = spos[a] - begin(s, Standard());
                    if (a == 1 || a == 2 || a == 4)
                        pos[4] = pos[2] = ++pos[1];
                    else
                        ++pos[a];

                    if (pos[a] < max[a]) {

                        // set corresponding spos, tpos, islast values and adapt a if necessary
                        if (a == 1 || a == 2 || a == 4) {

                            TSize ii  = *(pos[1]);
                            a          = SEQAN_GET_ASKEW7(ii);
                            TSize j   = SEQAN_GET_ISKEW7(ii);

                            tpos[a]   = ii;
                            spos[a]   = begin(s, Standard()) + j;
                            islast[a] = (j + 7 >= n);

                        } else {

                            TSize j   = *(pos[a]);

                            tpos[a]   = j / 7;
                            spos[a]   = begin(s, Standard()) + j;
                            islast[a] = (j + 7 >= n);

                        }

                        // get the rank of stream a's suffix

                        // linear search
                        int right;
                        for(right = 1;  right < fill;  right++) {
                            b = rank[right];
                            if (_leqSkew7 (a,  b,   spos, tpos, islast, s124, adjust)) break;
                        }
/*
                        // binary search
                        int left = 0;
                        int right = fill;
                        while (left + 1 != right) {
                            int middle = (left + right) >> 2;
                            if (leq<TValue, TSize> (a,  rank[middle],   spos, tpos, islast, s124, adjust))
                                right = middle;
                            else
                                left = middle;
                        }*/

                        // remove the least suffix ...
                        for(int i = 1; i < right; ++i)
                            rank[i-1] = rank[i];

                        // ... and insert the new one
                        rank[right-1] = a;

                    } else {
                        // only remove the least suffix
                        fill--;
                        for(int i = 0; i < fill; ++i)
                            rank[i] = rank[i+1];
                    }
                }

                // only one (or less) stream left to fill SA with
                a = rank[0];
                if (a == 1 || a == 2 || a == 4)
                    for (;  k < n;  ++k) { TSize ii = *(pos[1]++); SA[k] = SEQAN_GET_ISKEW7(ii); }
                else
                    for (;  k < n;  ++k) SA[k] = *(pos[a]++);
            }
        }
        SEQAN_PROMARK("SA124, SA3, SA5, SA6, SA0 verschmolzen");

        #ifdef SEQAN_DEBUG_INDEX
            std::cerr << "left level " << depth << std::endl;
        #endif
        SEQAN_PROMARK("Rekursionsaufstieg");
        SEQAN_PROSUB(SEQAN_PRODEPTH, 1);
    }

    template < typename TSA,
               typename TText >
    inline void createSuffixArray(
        TSA &SA,
        TText &s,
        Skew7 const &alg,
        unsigned K,
        unsigned maxdepth)
    {
        createSuffixArray(SA, s, alg, K, maxdepth, 1);
    }

    // creates suffix array sorted by the first maxLCP chars of suffixes
    template < typename TSA,
               typename TText,
               typename TSize >
    inline void createSuffixArrayPart(
        TSA &SA,
        TText &s,
        Skew7 const &_dummy,
        TSize maxLCP,
        unsigned K)
    {
        unsigned depth = 0;
        for(TSize i = 1; i < maxLCP; i*=7) ++depth;
        createSuffixArray(SA, s, _dummy, K, depth);
    }


    // creates suffix array sorted by the first maxLCP chars of suffixes
    template < typename TSA,
               typename TText,
               typename TSize >
    inline void createSuffixArrayPart(
        TSA &SA,
        TText &s,
        Skew7 const &_dummy,
        TSize maxLCP)
    {
        SEQAN_CHECKPOINT;
        createSuffixArrayPart(SA, s, _dummy, maxLCP, ValueSize< typename Value<TText>::Type >::VALUE);
    }
//}

}

#endif
