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

#ifndef SEQAN_HEADER_INDEX_SKEW3_H
#define SEQAN_HEADER_INDEX_SKEW3_H

namespace SEQAN_NAMESPACE_MAIN
{

//namespace SEQAN_NAMESPACE_PIPELINING
//{

    struct Skew3 {};


    //////////////////////////////////////////////////////////////////////////////
    // external Skew3 algorithm
    //////////////////////////////////////////////////////////////////////////////

    template <typename T>
    struct SkewDC_<3, T>
    {
        static const unsigned VALUE[];
    };

    template <typename T>
    const unsigned SkewDC_<3, T>::VALUE[] = { 2,   1, 2 };


    // *** COMPARATORS & MAPS ***

    template <typename TValue, typename TResult = int>
    struct _skew3NComp :
        public std::binary_function<TValue, TValue, TResult>
    {
        inline TResult operator() (const TValue &a, const TValue &b) const
        {
            typedef typename Value<TValue, 1>::Type                 TSize;
            typedef typename Value<TValue, 2>::Type                 TTriplet;
            typedef typename Value<TTriplet>::Type                  TTripletValue;
            typedef typename StoredTupleValue_<TTripletValue>::Type TStoredValue;

            const TStoredValue *sa = a.i2.i;
            const TStoredValue *sb = b.i2.i;

            TSize n = LENGTH<TTriplet>::VALUE;
            if (a.i1 < n) n = a.i1;
            if (b.i1 < n) n = b.i1;

            // compare the overlap of septets (the first n bases)
            for (TSize i = 0; i < n; i++, ++sa, ++sb)
            {
                if (*sa == *sb) continue;
                return (*sa < *sb)? -1 : 1;
            }

            // overlap is equal, now suffix lengths decide
            if (n < LENGTH<TTriplet>::VALUE)
                return (a.i1 < b.i1)? -1 : 1;
            else
                return 0;
        }
    };

    template <typename TValue, typename TResult = typename Value<TValue, 1>::Type>
    struct _skew3NMapLinear :
        public std::unary_function<TValue, TResult>
    {
        TResult BN;

        _skew3NMapLinear(TResult BN_) : BN(BN_) {}

        inline TResult operator() (const TValue & x) const
        {
            TResult i = x.i1; return BN - (i - i / 3);
        }
    };

    template <typename TValue, typename TResult = typename Value<TValue, 1>::Type>
    struct _skew3NMapSliced :
        public std::unary_function<TValue, TResult>
    {
        TResult BN, BN2;

        _skew3NMapSliced(TResult BN_) : BN(BN_-1), BN2(BN_/2-1) {}

        inline TResult operator() (const TValue & x) const
        {
            return (x.i1 % 3 == 1)? BN - x.i1/3 : BN2 - x.i1/3;
        }
    };


    template <typename TValue, typename TResult = TValue>
    struct _skew3UnslicerFunc :
        public std::unary_function<TValue, TResult>
    {
        TResult o1, o2, n2;

        _skew3UnslicerFunc(TResult N) :
            o1(N - (N + 2) % 3),
            o2(N - (N + 1) % 3),
            n2((N + 1) / 3) { }

        inline TResult operator() (const TValue & x) const
        {
            return (x < n2) ? o2 - x * 3 : o1 - (x - n2) * 3;
        }
    };

    template <typename TValue, typename TResult = typename Value<typename Value<TValue, 2>::Type>::Type>
    struct _skew3NMapExtended :
        public std::unary_function<TValue, TResult>
    {
        inline TResult operator() (const TValue & x) const
        {
            return x.i2[0];
        }
    };

    template <typename TValue, typename TResult = int>
    struct _skew3ExtendComp :
        public std::binary_function<TValue, TValue, TResult>
    {
        inline TResult operator() (const TValue &a, const TValue &b) const
        {
            return (a.i3[0] <  b.i3[0] ||
                   (a.i3[0] == b.i3[0] && a.i2[0] < b.i2[0])) ? -1 : 1;
        }
    };


    template < typename TInput >
    struct Value< Pipe< TInput, Skew3 > > :
        public Size<TInput> {};


    //////////////////////////////////////////////////////////////////////////////
    // Skew3 pipeline module

    template < typename TInput >
    struct Pipe< TInput, Skew3 >
    {

        // *** SPECIALIZATION ***

        // step 1
        typedef Pipe< TInput, Sampler<3> >  TSamplerDC3;
                                        typedef _skew3NComp<TypeOf_(TSamplerDC3)> ncomp_t;
        typedef Pool< TypeOf_(TSamplerDC3), SorterSpec< SorterConfigSize<ncomp_t, TSizeOf_(TSamplerDC3) > > > TSortTuples;
        typedef Pipe< TSortTuples, Namer<ncomp_t> > TNamer;
                                        typedef _skew3NMapSliced<TypeOf_(TNamer)> nmap_sliced_t;
                                        typedef _skew3NMapLinear<TypeOf_(TNamer)> nmap_linear_t;
        typedef Pool< TypeOf_(TNamer), MapperSpec< MapperConfigSize< nmap_sliced_t, TSizeOf_(TNamer) > > > TNames_Sliced;

        // unique names - shortcut
        typedef Pool< TypeOf_(TNames_Sliced), MapperSpec< MapperConfigSize< nmap_linear_t, TSizeOf_(TNames_Sliced) > > > TNames_Linear_Unique;

        // non-unique names - recursion
        typedef Pipe< TNames_Sliced, Filter< filterI2<TypeOf_(TNames_Sliced)> > > TFilter;
        typedef Pipe< TFilter, Skew3 > TRecurse;
                                        typedef _skew3UnslicerFunc<TypeOf_(TRecurse)> unslicer_func_t;
        typedef Pipe< TRecurse, Filter<unslicer_func_t> > TUnslicer;
        typedef Pipe< TUnslicer, Counter > TRenamer;
        typedef Pool< TypeOf_(TRenamer), MapperSpec< MapperConfigSize< nmap_linear_t, TSizeOf_(TRenamer) > > > TNames_Linear;

        // step 2
        typedef Pipe< Bundle2< TInput, TNames_Linear >, Extender3 > TExtender;
                                        typedef _skew3ExtendComp<TypeOf_(typename TExtender::Out0)> extend_comp_t;
        typedef Pool< TypeOf_(typename TExtender::Out0), SorterSpec< SorterConfigSize< extend_comp_t, TSizeOf_(typename TExtender::Out0) > > > TSorterS0;

        // step 3
                                        typedef _skew3NMapExtended<TypeOf_(typename TExtender::Out12)> nmap_extended_t;
        typedef Pool< TypeOf_(typename TExtender::Out12), MapperSpec< MapperConfigSize< nmap_extended_t, TSizeOf_(typename TExtender::Out12) > > > TSorterS12;
        typedef Pipe< Bundle2< TSorterS0, TSorterS12 >, Merger3 > TMerger;

        TSorterS0   sortedS0;
        TSorterS12  sortedS12;
        TMerger     in;

        Pipe() :
            in(bundle2(sortedS0, sortedS12)) {}

        Pipe(TInput& _textIn) :
            in(bundle2(sortedS0, sortedS12))
        {
            process(_textIn);
        }

        template < typename TInput_ >
        bool process(TInput_ &textIn)
        {

            SEQAN_PROADD(SEQAN_PRODEPTH, 1);
            SEQAN_PROMARK("Rekursionsabstieg");
            #ifdef SEQAN_DEBUG_INDEX
                std::cerr << "enter level " << SEQAN_PROVAL(SEQAN_PRODEPTH) << std::endl;
            #endif
            {


            // *** INSTANTIATION ***

            // step 1
            TSamplerDC3                 sampler(textIn);
            TSortTuples                 sorter;
            #ifdef SEQAN_DEBUG_INDEX
                std::cerr << "  sort names (" << length(sampler)<< ")" << std::endl;
            #endif
            sorter << sampler;
            SEQAN_PROMARK("Sorter (2) - Triplets sortieren");

            TNamer                      namer(sorter);
            nmap_sliced_t               map_sliced(length(namer));
            nmap_linear_t               map_linear(length(namer));
            TNames_Sliced               names_sliced(map_sliced);
            #ifdef SEQAN_DEBUG_INDEX
                std::cerr << "  slice names" << std::endl;
            #endif
            names_sliced << namer;

            if (namer.unique() || empty(names_sliced))
            {
                // unique names

                clear(sorter);
                SEQAN_PROMARK("Mapper (4) - s12 konstruieren");
                TNames_Linear_Unique        names_linear(map_linear);

                #ifdef SEQAN_DEBUG_INDEX
                    std::cerr << "  make names linear" << std::endl;
                #endif
                names_linear << names_sliced;
                clear(names_sliced);
                SEQAN_PROMARK("Mapper (10) - ISA12 konstruieren");

                // step 2
                _skew3Extend(textIn, names_linear, sortedS0, sortedS12);

            }
            else
            {
                // non-unique names

                clear(sorter);
                SEQAN_PROMARK("Mapper (4) - s12 konstruieren");

                TFilter                     filter(names_sliced);
                TRecurse                    recurse(filter);

                #ifdef SEQAN_TEST_SKEW3
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

                TNames_Linear               names_linear(map_linear);
                #ifdef SEQAN_DEBUG_INDEX
                    std::cerr << "  rename names" << std::endl;
                #endif
                names_linear << renamer;
                clear(renamer);
                SEQAN_PROMARK("Mapper (10) - ISA12 konstruieren");

                // step 2
                #ifdef SEQAN_DEBUG_INDEX
                    std::cerr << "  prepare merge" << std::endl;
                #endif
                _skew3Extend(textIn, names_linear, sortedS0, sortedS12);
                SEQAN_PROMARK("Mapper (12), Sorter (13) - SA12 und SA0 verschmelzen");
            }

            // step 3
            // ... is done on-demand by merger
            }
            #ifdef SEQAN_DEBUG_INDEX
                std::cerr << "left level " << SEQAN_PROVAL(SEQAN_PRODEPTH) << std::endl;
            #endif
            SEQAN_PROMARK("Rekursionsaufstieg");
            SEQAN_PROSUB(SEQAN_PRODEPTH, 1);
            return true;
        }

        inline typename Value<Pipe>::Type const & operator*()
        {
            return *in;
        }

        inline Pipe& operator++()
        {
            ++in;
            return *this;
        }
    };

    // not sure which interface is more intuitive, we support both
    // you can call "skew << pipe" or "skew_t skew(pipe); skew.process()"
    // for the first we would need no _in member
    template < typename TInput, typename TObject >
    inline bool operator<<(Pipe< TInput, Skew3 > &me, TObject &textIn)
    {
        return me.process(textIn);
    }

    template <
        typename TSA,
        typename TValue,
        typename TConfig >
    inline void createSuffixArray(
        TSA &SA,
        String< TValue, External<TConfig> > &s,
        Skew3 const &spec,
        unsigned /*K*/,
        unsigned /*maxdepth*/)
    {
        createSuffixArrayExt(SA, s, spec);
    }


    //////////////////////////////////////////////////////////////////////////////
    // internal Skew3 algorithm
    //////////////////////////////////////////////////////////////////////////////

    //////////////////////////////////////////////////////////////////////////////
    // typedefs and helpers

    template <typename T, typename ST> inline
    bool _leqSkew3(T a1, ST a2,   T b1, ST b2)
    { // lexic. order for pairs
        return (ordLess(a1, b1) || (a1 == b1 && a2 <= b2));
    }

    template <typename T, typename ST> inline
    bool _leqSkew3(T a1, T a2, ST a3,   T b1, T b2, ST b3)
    { // and triples
        return (ordLess(a1, b1) || (a1 == b1 && _leqSkew3(a2,a3, b2,b3)));
    }


    //////////////////////////////////////////////////////////////////////////////
    // Skew Implementation
    // find the suffix array SA of s[0..n-1] in {0..K}^n
    //
    // * no trailing 0's required
    // * no dummy triples in special cases

    // creates suffix array SA of s
    // chars have to be in the range [0,K)
    template < typename TSA,
               typename TText >
    void createSuffixArray(
        TSA &SA,
        TText &s,
        Skew3 const &,
        unsigned K,
        unsigned maxdepth,
        unsigned depth)
    {
        typedef typename Value<TSA>::Type TSize;
        typedef typename Value<TText>::Type TValue;

        #ifdef SEQAN_DEBUG_INDEX
            std::cerr << "--- CREATE SUFFIX ARRAY ---" << std::endl;
            std::cerr << "Skew3 [random access]" << std::endl;
        #endif

        TSize n = length(s);
        if (n < 1) return;

        TSize n0=n/3, n1=(n+2)/3, n2=(n+1)/3, n12=n1+n2;
        TSize         o1=(n+2)%3, o2=(n+1)%3;

        SEQAN_PROSET(SEQAN_PRODEPTH, depth);
        SEQAN_PROSET(SEQAN_PROEXTRA1, K);
        SEQAN_PROMARK("Rekursionsabstieg");
        #ifdef SEQAN_DEBUG_INDEX
            std::cerr << "enter level " << depth << " (" << n << ")" << std::endl;
        #endif

        String<TSize, Alloc<> > s12;
        resize(s12, n12, Exact());
        // we use SA[n0..n-1] as a temporary buffer instead of allocating one;
        typename Suffix<TSA>::Type SA12 = suffix(SA, n0);


        // generate positions of mod 1 and mod 2 suffixes
        {
            //for (TSize i=0, j=0;  i < n;  i++) if ((n-i)%3) s12[j++] = i;
            s12[0] = o1;
            for (TSize i=o2, j=n1-n2;  i < n;  i++)
            {
                s12[j++] = i++;
                s12[j++] = i++;
            }
        }


        // lsb radix sort the mod 1 and mod 2 triples
        {
            String<TSize, Alloc<> > cnt;
            resize(cnt, K, Exact());    // counter array

            radixPass(SA12, s12, s, cnt, K, 2);
            radixPass(s12, SA12, s, cnt, K, 1);
            radixPass(SA12, s12, s, cnt, K);
        }
        SEQAN_PROMARK("Triplets sortiert");

        // find lexicographic names of triples
        TSize name = 0;

        bool differ = true;
        TValue c0 = TValue(), c1 = TValue(), c2 = TValue();
        for (TSize i = 0, clip = _max(n, (TSize)2) - 2, l;  i < n12;  i++)
        {
            if ((l = SA12[i]) < clip)
            {
                if (differ || s[l] != c0 || s[l+1] != c1 || s[l+2] != c2)
                {
                    name++;  c0 = s[l];  c1 = s[l+1];  c2 = s[l+2];
                    differ = false;
                }
            }
            else
            {
                name++;
                differ = true;  // the last 2 triples always differ from the rest
            }
            s12[(n-l) % 3 == 2 ? l/3 : l/3 + n2] = name - 1;    // left or right half
        }
        SEQAN_PROMARK("s12 konstruiert");

        // recurse if names are not yet unique
        if (name < n12)
        {
            if (depth != maxdepth)
            {
                createSuffixArray(SA12, s12, Skew3(), name, maxdepth, depth + 1);
                #ifdef SEQAN_TEST_SKEW3
                    SEQAN_ASSERT(isSuffixArray(SA12, s12));
                #endif
            }
            // store unique names in s12 using the suffix array
            for (TSize i = 0;  i < n12;  i++) s12[SA12[i]] = i;
        }
        else // generate the suffix array of s12 directly
            for (TSize i = 0;  i < n12;  i++) SA12[s12[i]] = i;



        // use SA[0...n0-1] as a temporary buffer instead of s0
        // and allocate SA0

        {
            String<TSize, Alloc<> > SA0;
            resize(SA0, n0, Exact());
            typename Infix<TSA>::Type s0 = infix(SA, 0, n0);

            // stably sort the mod 0 suffixes from SA12 by their first character
            {
                String<TSize, Alloc<> > cnt;
                resize(cnt, K, Exact());    // counter array

                for (TSize i=0, j=0, l;  i < n12;  i++)
                    if ((l = SA12[i]) < n2)
                        if ((l = o2 + (3 * l) ) > 0)
                            s0[j++] = l - 1;
                radixPass(SA0, s0, s, cnt, K);
            }
            SEQAN_PROMARK("SA0 konstruiert");

            // merge sorted SA0 suffixes and sorted SA12 suffixes
            #define SEQAN_GET_ISKEW3(ii) (ii < n2 ? ii * 3 + o2 : (ii - n2) * 3 + o1)
            if (n0)
            {
                for (TSize p=0,  t=0,  k=0,  clip = n - 1;  k < n;  k++)
                {
                    TSize ii = SA12[t];                    // pos of current interleave offset 12 suffix
                    TSize i  = SEQAN_GET_ISKEW3(ii);    // pos of current offset 12 suffix
                    TSize j  = SA0[p];                    // pos of current offset 0  suffix
                    if (ii < n2 ?
                        _leqSkew3<TValue, TSize> (s[i],       s12[ii + n1],     s[j],        s12[j/3 + n2  - n0]) :
                    (i < clip) ?     // clip if 12 suffix is the last
                        _leqSkew3<TValue, TSize> (s[i],s[i+1],s12[ii + 1 - n1], s[j],s[j+1], s12[j/3 + n12 - n0]) :
                        s[i] <= s[j])

                    { // suffix from SA12 is smaller
                        SA[k] = i;
                        if (++t == n12) // done --- only SA0 suffixes left
                            for (;  p < n0;  p++) SA[++k] = SA0[p];
                    }
                    else
                    {
                        SA[k] = j;
                        if (++p == n0)  // done --- only SA12 suffixes left
                            for (;  t < n12;  t++) { ii = SA12[t]; SA[++k] = SEQAN_GET_ISKEW3(ii); }
                    }
                }
            } else
                for (TSize t = 0;  t < n12;  t++) { TSize ii = SA12[t]; SA[t] = SEQAN_GET_ISKEW3(ii); }
        }
        SEQAN_PROMARK("SA12 und SA0 verschmolzen");

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
        Skew3 const &alg,
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
        Skew3 const &_dummy,
        TSize maxLCP,
        unsigned K)
    {
        unsigned depth = 0;
        for(TSize i = 1; i < maxLCP; i*=3) ++depth;
        createSuffixArray(SA, s, _dummy, K, depth);
    }


    template < typename TSA,
               typename TText,
               typename TSize >
    inline void createSuffixArrayPart(
        TSA &SA,
        TText &s,
        Skew3 const &_dummy,
        TSize maxLCP)
    {
        SEQAN_CHECKPOINT;
        createSuffixArrayPart(SA, s, _dummy, maxLCP, ValueSize< typename Value<TText>::Type >::VALUE);
    }

//}

}

#endif
