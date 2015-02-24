// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2015, Knut Reinert, FU Berlin
// Copyright (c) 2013 NVIDIA Corporation
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

#ifndef SEQAN_HEADER_INDEX_SHIMS_H
#define SEQAN_HEADER_INDEX_SHIMS_H

namespace SEQAN_NAMESPACE_MAIN
{

    //////////////////////////////////////////////////////////////////////////////
    // Suffix Array creation wrappers
    //////////////////////////////////////////////////////////////////////////////

    // build suffix array with an external pipeling algorithm (skew3, skew7, ...)
    template <
        typename TSA,
        typename TObject,
        typename TAlgSpec >
    void _createSuffixArrayPipelining(
        TSA &suffixArray,
        TObject const &text,
        TAlgSpec const)
    {
    SEQAN_CHECKPOINT
        // signed characters behave different than unsigned when compared
        // to get the same index with signed or unsigned chars we simply cast them to unsigned
        // before feeding them into the pipeline
        typedef typename MakeUnsigned_< typename Value<TObject>::Type >::Type TUValue;

        // specialization
        typedef Pipe< TObject, Source<> >                src_t;
        typedef Pipe< src_t, Caster<TUValue> >          unsigner_t;
        typedef Pipe< unsigner_t, TAlgSpec >            creator_t;

        // instantiation and processing
        src_t        src(text);
        unsigner_t  unsigner(src);
        creator_t    creator(unsigner);

        suffixArray << creator;
        #ifdef SEQAN_TEST_INDEX
            isSuffixArray(suffixArray, text);
        #endif
    }


    // build suffix array (external) for mutliple sequences
    template <
        typename TSA,
        typename TString,
        typename TSpec,
        typename TAlgSpec >
    void _createSuffixArrayPipelining(
        TSA &suffixArray,
        StringSet<TString, TSpec> const &stringSet,
        TAlgSpec const)
    {
    SEQAN_CHECKPOINT
        // signed characters behave different than unsigned when compared
        // to get the same index with signed or unsigned chars we simply cast them to unsigned
        // before feeding them into the pipeline
        typedef typename Concatenator<StringSet<TString, TSpec> >::Type            TConcat;
        typedef typename MakeUnsigned_< typename Value<TConcat>::Type >::Type    TUValue;
        typedef Multi<
            TAlgSpec,
            typename Value<TSA>::Type,
            typename StringSetLimits<StringSet<TString, TSpec> >::Type >        MultiConstrSpec;

        // specialization
        typedef Pipe< TConcat, Source<> >                src_t;
        typedef Pipe< src_t, Caster<TUValue> >          unsigner_t;
        typedef Pipe< unsigner_t, MultiConstrSpec >        creator_t;

        // instantiation and processing
        src_t        src(concat(stringSet));
        unsigner_t  unsigner(src);
        creator_t    creator(unsigner, stringSetLimits(stringSet));

        suffixArray << creator;
        #ifdef SEQAN_TEST_INDEX
            //isSuffixArray(suffixArray, stringSet);
        #endif
    }


/*!
 * @fn createSuffixArray
 * @headerfile <seqan/index.h>
 * @brief Creates a suffix array from a given text.
 *
 * @signature void createSuffixArray(suffixArray, text[, algoTag]);
 *
 * @param[out] suffix  Array The resulting suffix array.
 * @param[in]  text    A given text. Types: @link ContainerConcept @endlink
 * @param[in]  algoTag A tag that identifies the algorithm which is used for creation.
 *
 * This function should not be called directly.  Please use @link Index#indexCreate
 * @endlink or @link Index#indexRequire @endlink.  The size of <tt>suffixArray</tt>
 * must be at least <tt>length(text)</tt> before calling this function.
 *
 * @see DemoSuffixArray
 */
    template < typename TSA,
               typename TText,
               typename TAlgSpec >
    inline void _createSuffixArrayRandomAccess(
        TSA &sa,
        TText const &s,
        TAlgSpec const &alg)
    {
    SEQAN_CHECKPOINT
        // -> call internal memory algorithm with an extended interface (+ alphabet size, max_depth)
        if (BitsPerValue< typename Value<TText>::Type >::VALUE > 16)
            createSuffixArray(sa, s, alg, length(s), 0);
        else
            createSuffixArray(sa, s, alg, ValueSize< typename Value<TText>::Type >::VALUE, 0);
    }

    template <
        typename TSA,
        typename TText,
        typename TAlgSpec >
    inline void _createSuffixArrayWrapper(
        TSA &sa,
        TText const &s,
        TAlgSpec const &alg,
        True)
    {
    SEQAN_CHECKPOINT
        _createSuffixArrayRandomAccess(sa, s, alg);
    }

    // always use external Skew7 for multiple strings
    template <
        typename TSA,
        typename TSequence,
        typename TSetSpec,
        typename TAlgSpec >
    inline void _createSuffixArrayWrapper(
        TSA &sa,
        StringSet<TSequence, TSetSpec> const &s,
        TAlgSpec const &,
        True)
    {
    SEQAN_CHECKPOINT
        _createSuffixArrayPipelining(sa, s, Skew7());
    }

    template <
        typename TSA,
        typename TText,
        typename TAlgSpec >
    inline void _createSuffixArrayWrapper(
        TSA &sa,
        TText const &s,
        TAlgSpec const &alg,
        False)
    {
    SEQAN_CHECKPOINT
        _createSuffixArrayPipelining(sa, s, alg);
    }

    template <
        typename TSA,
        typename TText,
        typename TAlgSpec >
    inline void createSuffixArray(
        TSA &sa,
        TText const &s,
        TAlgSpec const &alg)
    {
    SEQAN_CHECKPOINT
        _createSuffixArrayWrapper(sa, s, alg, typename SACreatorRandomAccess_<TSA, TText, TAlgSpec>::Type());
    }

//____________________________________________________________________________

/*!
 * @fn createInvSuffixArray
 * @headerfile <seqan/index.h>
 * @brief Creates the inverse suffix array from a given suffix array.
 *
 * @signature void createInvSuffixArray(invSuffixArray, suffixArray);
 *
 * @param[out] invSuffixArray  The resulting inverse suffix array.
 * @param[in]  suffixArray     The precomputed suffix array for some text.
 *
 * This function should not be called directly. Please use @link Index#indexCreate
 * @endlink or @link Index#indexRequire @endlink. The size of <tt>invSuffixArray</tt> must be at
 * least <tt>length(suffixArray)</tt> before calling this function.
 *
 * The complexity is linear in size of the suffix array.
 */

    template <typename TIsa, typename TSa, typename TParallel>
    inline void
    createInvSuffixArray(TIsa &isa,
                         TSa const &sa,
                         FromSortedSa<TParallel> const &/*alg*/)
    {
        typedef typename Size<TSa>::Type                     TSize;
        typedef typename MakeSigned<TSize>::Type             TSignedSize;
        typedef typename StringSetLimits<TSa>::Type          TLimits;
        typedef typename Iterator<TSa const, Standard>::Type TIter;

        TLimits const & limits = stringSetLimits(sa);
        Splitter<TSize> splitter(0, length(sa), TParallel());

        SEQAN_OMP_PRAGMA(parallel for)
        for (TSignedSize job = 0; job < static_cast<TSignedSize>(length(splitter)); ++job)
        {
            TIter saIt = begin(sa, Standard()) + splitter[job];
            TIter saItEnd = begin(sa, Standard()) + splitter[job + 1];
            TSize pos = splitter[job];

            for (; saIt != saItEnd; ++saIt)
                isa[posGlobalize(*saIt, limits)] = pos++;
        }
    }

//____________________________________________________________________________



    //////////////////////////////////////////////////////////////////////////////
    // LCP Table creation wrappers
    //////////////////////////////////////////////////////////////////////////////

    template <
        typename TLCPTable,
        typename TObject,
        typename TSA,
        typename TAlgSpec >
    void _createLCPTablePipelining(
        TLCPTable &LCP,
        TObject const &text,
        TSA const &suffixArray,
        TAlgSpec const)
    {
    SEQAN_CHECKPOINT
        // specialization
        typedef Pipe< TObject, Source<> >                            srcText_t;
        typedef Pipe< TSA, Source<> >                               srcSA_t;
        typedef Pipe< Bundle2< srcText_t, srcSA_t >, TAlgSpec >        creator_t;

        // instantiation and processing
        srcText_t    srcText(text);
        srcSA_t        srcSA(suffixArray);
        creator_t    creator(bundle2(srcText, srcSA));

        LCP << creator;
        #ifdef SEQAN_TEST_INDEX
            isLCPTable(LCP, suffixArray, text);
        #endif
    }


    // build lcp table (external) for mutliple sequences
    template <
        typename TLCPTable,
        typename TString,
        typename TSpec,
        typename TSA,
        typename TAlgSpec >
    void _createLCPTablePipelining(
        TLCPTable &LCP,
        StringSet<TString, TSpec> const &stringSet,
        TSA const &suffixArray,
        TAlgSpec const)
    {
    SEQAN_CHECKPOINT
        typedef typename Concatenator<StringSet<TString, TSpec> >::Type TConcat;
        typedef Multi<
            TAlgSpec,
            typename Value<TSA>::Type,
            typename StringSetLimits<StringSet<TString, TSpec> >::Type > MultiConstrSpec;

        // specialization
        typedef Pipe< TConcat, Source<> >                                srcText_t;
        typedef Pipe< TSA, Source<> >                                   srcSA_t;
        typedef Pipe< Bundle2< srcText_t, srcSA_t >, MultiConstrSpec >    creator_t;

        // instantiation and processing
        srcText_t    srcText(concat(stringSet));
        srcSA_t        srcSA(suffixArray);
        creator_t    creator(bundle2(srcText, srcSA), stringSetLimits(stringSet));

        LCP << creator;
        #ifdef SEQAN_TEST_INDEX
            isLCPTable(LCP, suffixArray, text);
        #endif
    }


/*!
 * @fn createLcpTable
 * @headerfile <seqan/index.h>
 * @brief Creates a LCP table from a given text and suffix array.
 *
 * @signature void createLcpTable(lcp, text, suffixArray[, algoTag]);
 *
 * @param[out] lcp         The resulting lcp table.
 * @param[in]  text        A given text. Types: @link ContainerConcept @endlink
 * @param[in]  suffixArray The suffix array of <tt>text</tt>.
 * @param[in]  algoTag     A tag that identifies the algorithm which is used for creation.
 *
 * This function should not be called directly.  Please use @link Index#indexCreate
 * @endlink or @link Index#indexRequire @endlink.  The size of <tt>lcp</tt> must be at
 * least <tt>length(text)</tt> before calling this function.
 */

    template <
        typename TLCP,
        typename TText,
        typename TSA,
        typename TAlgSpec >
    inline void _createLCPTableWrapper(
        TLCP &lcp,
        TText const &s,
        TSA const &sa,
        TAlgSpec const &alg,
        True)
    {
    SEQAN_CHECKPOINT
        _createLCPTableRandomAccess(lcp, s, sa, alg);
    }

    template <
        typename TLCP,
        typename TText,
        typename TSA,
        typename TAlgSpec >
    inline void _createLCPTableWrapper(
        TLCP &lcp,
        TText const &s,
        TSA const &sa,
        TAlgSpec const &alg,
        False)
    {
    SEQAN_CHECKPOINT
        _createLCPTablePipelining(lcp, s, sa, alg);
    }

    template <
        typename TLCP,
        typename TText,
        typename TSA,
        typename TAlgSpec >
    inline void createLcpTable(
        TLCP &lcp,
        TText const &s,
        TSA const &sa,
        TAlgSpec const &alg)
    {
    SEQAN_CHECKPOINT
        _createLCPTableWrapper(lcp, s, sa, alg, typename LcpCreatorRandomAccess_<TLCP, TText, TSA, TAlgSpec>::Type());
    }


//____________________________________________________________________________



    //////////////////////////////////////////////////////////////////////////////
    // Enhanced LCP Table creation wrappers
    //////////////////////////////////////////////////////////////////////////////

    // build enhanced LCP table with an external pipelining algorithm (ext kasai, ...)
    // and a dynamic programming tree construction alg
    // (in contrast to the LCP table the enhanced LCP table contains the lcp-values
    // of suffix intervals used in the binary search)
    template <
        typename TValue,
        typename TSpec,
        typename TObject,
        typename TSA,
        typename TLCP,
        typename TAlgSpec >
    void createLcpeTableExt(
        String< TValue, TSpec > &LCPE,
        TObject const &text,
        TSA const &suffixArray,
        TLCP const & /*LCP*/,
        TAlgSpec const)
    {
    SEQAN_CHECKPOINT
        typedef typename Concatenator<TObject>::Type                TConcat;

        // specialization
        typedef Pipe< TConcat, Source<> >                            srcText_t;
        typedef Pipe< TSA, Source<> >                                srcSA_t;
        typedef Pipe< Bundle2< srcText_t, srcSA_t >, TAlgSpec >        creator_t;

        // instantiation and processing
        srcText_t    srcText(concat(text));
        srcSA_t        srcSA(suffixArray);
        creator_t    creator(bundle2(srcText, srcSA));

        #ifdef SEQAN_TEST_INDEX
            isLCPTable(creator, suffixArray, text);
        #endif
        createLcpBinTree(LCPE, creator);
    }

    // build enhanced LCP table with an lcp algorithm
    // and a dynamic programming tree construction alg
    template <
        typename TValue,
        typename TSpec,
        typename TText,
        typename TSA,
        typename TLCP,
        typename TAlgSpec >
    void createLcpeTable(
        String< TValue, TSpec > &LCPE,
        TText const &s,
        TSA const &,
        TLCP const &LCP,
        TAlgSpec const)
    {
    SEQAN_CHECKPOINT
        //TSA LCP;
        //resize(LCP, length(s), Exact());
        // we use LCPE[n-lcpSize..n-1] as a temporary buffer instead of allocating one
        typename Size<TText>::Type lcpSize = length(s) > 1? length(s) - 1: 0;
        typename Suffix<String< TValue, TSpec > >::Type LCPcopy = suffix(LCPE, length(LCPE) - lcpSize);
        LCPcopy = prefix(LCP, lcpSize);
        createLcpBinTree(LCPE, LCP);
    }

    template <
        typename TValue,
        typename TConfig,
        typename TText,
        typename TSA,
        typename TLCP,
        typename TAlgSpec >
    void createLcpeTable(
        String< TValue, External<TConfig> > &LCPE,
        TText const &s,
        TSA const &SA,
        TLCP const &LCP,
        TAlgSpec const alg)
    {
    SEQAN_CHECKPOINT
        createLcpeTableExt(LCPE, s, SA, LCP, alg);
    }

    template <
        typename TValue,
        typename TSpec,
        typename TText,
        typename TSA,
        typename TLCP>
    inline void createLcpeTable(
        String< TValue, TSpec > &LCPE,
        TText &s,
        TSA &SA,
        TLCP &LCP)
    {
    SEQAN_CHECKPOINT
        createLcpeTable(LCPE, s, SA, LCP, Kasai());
    }


//____________________________________________________________________________



    //////////////////////////////////////////////////////////////////////////////
    // Burrows-Wheeler-Table creation wrappers
    //////////////////////////////////////////////////////////////////////////////

    template < typename TBWT,
               typename TText,
               typename TSA >
    void createBWTableExt(
        TBWT &bwt,
        TText const &s,
        TSA const &SA)
    {
    SEQAN_CHECKPOINT
        // specialization
        typedef Pipe< TText, Source<> >                        srcText_t;
        typedef Pipe< TSA, Source<> >                       srcSA_t;
        typedef Pipe< Bundle2< srcText_t, srcSA_t >, Bwt >    creator_t;

        // instantiation and processing
        srcText_t    srcText(s);
        srcSA_t        srcSA(SA);
        creator_t    creator(bundle2(srcText, srcSA));

        bwt << creator;
    }

/*!
 * @fn createBWTable
 * @headerfile <seqan/index.h>
 * @brief Creates a Burrows-Wheeler table from a given text and suffix array.
 *
 * @signature void createBWTable(bwt, text, suffixArray[, algoTag]);
 *
 * @param[out] bwt     The resulting Burrows-Wheeler table.
 * @param[in]  suffix  Array The suffix array of <tt>text</tt>.
 * @param[in]  text    A given text. Types: @link ContainerConcept @endlink
 * @param[in]  algoTag A tag that identifies the algorithm which is used for creation.
 *
 * This function should not be called directly.  Please use @link Index#indexCreate
 * @endlink or @link Index#indexRequire @endlink.  The size of <tt>bwt</tt> must be at
 * least <tt>length(text)</tt> before calling this function.
 */
    // default
    template < typename TBWT, typename TText, typename TSA, typename TTextRandom_ >
    inline void _createBWTableWrapper(TBWT &bwt, TText const &s, TSA const &sa,        TTextRandom_ const)
    {
    SEQAN_CHECKPOINT
        createBWTableExt(bwt, concat(s), sa);
    }

    // text supports fast random access
    template < typename TBWT, typename TText, typename TSA >
    inline void _createBWTableWrapper(TBWT &bwt, TText const &s, TSA const &sa,        True const)
    {
    SEQAN_CHECKPOINT
        createBWTableInt(bwt, concat(s), sa);
    }

    template < typename TBWT, typename TText, typename TSA >
    inline void createBWTable(TBWT &bwt, TText const &s, TSA const &sa)
    {
    SEQAN_CHECKPOINT
        _createBWTableWrapper(bwt, s, sa, typename AllowsFastRandomAccess<TText>::Type());
    }


//////////////////////////////////////////////////////////////////////////////

    template <typename TOccValue>
    struct SAValueLess_:
        public std::less<TOccValue> {};

    template <typename T1, typename T2, typename TPack>
    struct SAValueLess_< Pair<T1,T2,TPack> >:
        public std::binary_function< Pair<T1,T2,TPack>, Pair<T1,T2,TPack>, bool>
    {
        inline bool operator()(Pair<T1,T2,TPack> const &a, Pair<T1,T2,TPack> const &b) const {
            return    (getValueI1(a) < getValueI1(b)) ||
                    ((getValueI1(a) == getValueI1(b)) && (getValueI2(a) < getValueI2(b)));
        }
    };

//TODO(singer): Make this internal
/*!
 * @fn orderOccurrences
 * @headerfile <seqan/index.h>
 * @brief Sorts a string of occurrences.
 *
 * @signature void orderOccurrences(occString);
 *
 * @param[in,out] occString String of occurrences.
 *
 * The occurrences are sorted by increasing positions.
 *
 * @see DemoMummy
 * @see DemoSupermaximalRepeats
 * @see DemoMaximalUniqueMatches
 *
 * @see VSTreeIterator#getOccurrences
 * @see IndexQGram#getOccurrences
 * @see SAValue
 */
    template <typename TValue, typename TSpec>
    inline void orderOccurrences(String<TValue, TSpec> &occString)
    {
    SEQAN_CHECKPOINT
        std::sort(begin(occString, Standard()), end(occString, Standard()), SAValueLess_<TValue>());
    }


//////////////////////////////////////////////////////////////////////////////
// fibre creators

/*!
 * @fn Index#indexCreate
 * @headerfile <seqan/index.h>
 * @brief Creates a specific @link Fibre @endlink.
 *
 * @signature bool indexCreate(index, fibreTag[, algoTag]);
 *
 * @param[in]     fibreTag A tag that identifies the @link Fibre @endlink
 * @param[in]     algoTag  A tag that identifies the algorithm which is used to create the fibre.  Default: The
 *                         result of @link Index#DefaultIndexCreator @endlink.
 * @param[in,out] index    The @link Index @endlink object holding the fibre.
 *
 * @return bool <tt>true</tt> on a success and false <tt>otherwise</tt>
 *
 * <tt>indexCreate</tt> calls the fibre corresponding <tt>createXXX(..)</tt> function (e.g. @link createSuffixArray
 * @endlink).
 */
    template <typename TText, typename TSpec, typename TSpecAlg>
    inline bool indexCreate(Index<TText, TSpec> &index, FibreSA, TSpecAlg const alg) {
    SEQAN_CHECKPOINT
        resize(indexSA(index), length(indexRawText(index)), Exact());
        createSuffixArray(indexSA(index), indexText(index), alg);
        return true;
    }

    template <typename TText, typename TSpec, typename TParallel>
    inline bool indexCreate(Index<TText, TSpec> &index, FibreIsa, FromSortedSa<TParallel> const alg)
    {
        resize(indexIsa(index), length(indexRawText(index)), Exact());
        createInvSuffixArray(indexIsa(index), indexSA(index), alg);
        return true;
    }

    template <typename TText, typename TSpec, typename TSpecAlg>
    inline bool indexCreate(Index<TText, TSpec> &index, FibreLcp, TSpecAlg const alg) {
    SEQAN_CHECKPOINT
        resize(indexLcp(index), length(indexRawText(index)), Exact());
        createLcpTable(indexLcp(index), indexText(index), indexSA(index), alg);
        return true;
    }

    template <typename TText, typename TSpec, typename TSpecAlg>
    inline bool indexCreate(Index<TText, TSpec> &index, FibreLcpe, TSpecAlg const alg) {
    SEQAN_CHECKPOINT
    //TODO: separate LCP from LCPE (for now LCPE = LCP + extra)
        resize(indexLcpe(index), sizeofLcpe(lengthSum(index)), Exact());
        createLcpeTable(indexLcpe(index), indexRawText(index), indexSA(index), indexLcp(index), alg);
        return true;
//        return false;
    }

    template <typename TText, typename TSpec>
    inline bool indexCreate(Index<TText, TSpec> &index, FibreBwt, Bwt const) {
    SEQAN_CHECKPOINT
        resize(indexBwt(index), length(indexRawText(index)), Exact());
        createBWTable(indexBwt(index), indexText(index), indexRawSA(index));
        return true;
    }

    template <typename TText, typename TSpec>
    inline bool indexCreate(Index<TText, TSpec> &index, FibreChildtab, Childtab const) {
    SEQAN_CHECKPOINT
        resize(indexChildtab(index), length(indexRawText(index)), Exact());
        createChildtab(indexChildtab(index), indexLcp(index));
        return true;
    }

    template <typename TText, typename TSpec, typename TFibre>
    inline bool indexCreate(Index<TText, TSpec> &index, Tag<TFibre> const fibre) {
    SEQAN_CHECKPOINT
        return indexCreate(index, fibre, typename DefaultIndexCreator<Index<TText, TSpec>, Tag<TFibre> const>::Type());
    }


// ----------------------------------------------------------------------------
// Function indexCreate()
// ----------------------------------------------------------------------------

template <typename TText, typename TSpec>
inline bool indexCreate(Index<TText, TSpec> & index, FibreSA, Trie)
{
    typedef Index<TText, TSpec>                     TIndex;
    typedef typename Fibre<TIndex, FibreSA>::Type   TSA;
    typedef typename Value<TSA>::Type               TSAValue;
    typedef typename Size<TText>::Type              TSize;
    typedef QGramLess_<TSAValue, TText const>       TLess;

    TText const & text = indexText(index);
    TSA & sa = indexSA(index);
    TSize textLen = length(text);

    resize(sa, textLen, Exact());

    // Fill the suffix array with (i, 0).
    for (TSize i = 0; i < textLen; i++)
        sa[i] = TSAValue(i, 0);

    // Sort the suffix array using quicksort.
    sort(sa, TLess(text, maxLength(text)));

    return true;
}

//////////////////////////////////////////////////////////////////////////////
// automatic fibre creation

/*!
 * @fn Index#indexSupplied
 * @headerfile <seqan/index.h>
 * @brief Returns whether a specific @link Fibre @endlink is present.
 *
 * @signature bool indexSupplied(index, fibreTag);
 *
 * @param[in] index    The @link Index @endlink object holding the fibre.
 * @param[in] fibreTag A tag that identifies the @link Fibre @endlink Index Fibres.
 *
 * @return bool <tt>true</tt>, iff the fibre is present.
 */
    template <typename TText, typename TSpec, typename TFibre>
    SEQAN_HOST_DEVICE inline bool indexSupplied(Index<TText, TSpec> &index, Tag<TFibre> const fibre) {
    SEQAN_CHECKPOINT
        return !empty(getFibre(index, fibre));
    }


//////////////////////////////////////////////////////////////////////////////
/*!
 * @fn Index#indexRequire
 * @headerfile <seqan/index.h>
 * @brief On-demand creation of a specific @link Fibre @endlink.
 *
 * @signature bool indexRequire(index, fibreTag);
 *
 * @param[in,out] index    The @link Index @endlink object holding the fibre.
 * @param[in]     fibreTag A tag that identifies the @link Fibre @endlink
 *
 * @return bool <tt>true</tt> on a successful creation.
 *
 * If the fibre already exists (@link Index#indexSupplied @endlink is true) then <tt>indexRequire</tt> does nothing. If
 * the fibre doesn't exist then @link Index#indexCreate @endlink is called to create it.
 *
 * @section Examples
 *
 * The following code shows how the BWT of an text can be computed.
 *
 * @include demos/index/index_textAt_indexText_saAt_indexRequire.cpp
 *
 * The output is as follows:
 *
 * @include demos/index/index_textAt_indexText_saAt_indexRequire.cpp.stdout
 */

    template <typename TText, typename TSpec, typename TFibre>
    inline bool indexRequire(Index<TText, TSpec> &index, Tag<TFibre> const fibre) {
    SEQAN_CHECKPOINT
        if (indexSupplied(index, fibre)) return true;                // if the table doesn't exist,
        if (!indexSolveDependencies(index, fibre)) return false;    // fulfill requirements
        return indexCreate(index, fibre);                            // and create table
    }


//////////////////////////////////////////////////////////////////////////////
// index cargo interface

    template <typename TText, typename TSpec>
    inline typename Reference< typename Cargo<Index<TText, TSpec> >::Type >::Type
    cargo(Index<TText, TSpec> & me)
    {
    SEQAN_CHECKPOINT
        return me.cargo;
    }

    template <typename TText, typename TSpec>
    inline typename Reference< typename Cargo<Index<TText, TSpec> const>::Type >::Type
    cargo(Index<TText, TSpec> const & me)
    {
    SEQAN_CHECKPOINT
        return me.cargo;
    }

//////////////////////////////////////////////////////////////////////////////
// solve dependencies

    template <typename TText, typename TSpec, typename TFibre>
    inline bool indexSolveDependencies(Index<TText, TSpec> &, Tag<TFibre> const) {
    SEQAN_CHECKPOINT
        return true;
    }

    template <typename TText, typename TSpec>
    inline bool indexSolveDependencies(Index<TText, TSpec> &index, FibreIsa) {
    SEQAN_CHECKPOINT
        return indexRequire(index, FibreSA());
    }

    template <typename TText, typename TSpec>
    inline bool indexSolveDependencies(Index<TText, TSpec> &index, FibreLcp) {
    SEQAN_CHECKPOINT
        return indexRequire(index, FibreSA());
    }

    template <typename TText, typename TSpec>
    inline bool indexSolveDependencies(Index<TText, TSpec> &index, FibreLcpe) {
    SEQAN_CHECKPOINT
        return indexRequire(index, FibreLcp());
    }

    template <typename TText, typename TSpec>
    inline bool indexSolveDependencies(Index<TText, TSpec> &index, FibreChildtab) {
    SEQAN_CHECKPOINT
        return indexRequire(index, FibreLcp());
    }

    template <typename TText, typename TSpec>
    inline bool indexSolveDependencies(Index<TText, TSpec> &index, FibreBwt) {
    SEQAN_CHECKPOINT
        return indexRequire(index, FibreSA());
    }


//////////////////////////////////////////////////////////////////////////////
// open

    // TODO(esiragusa): Move open() / save() to respective String/StringSet classes.

    template <typename TValue>
    inline bool open(TValue & value, const char *fileName, int openMode)
    {
        String<TValue, External< ExternalConfigLarge<> > > extString;
        if (!open(extString, fileName, openMode & ~OPEN_CREATE)) return false;
        if (!empty(extString)) assign(value, front(extString));
        return true;
    }

    template <typename TValue>
    inline bool open(TValue & value, const char *fileName)
    {
        return open(value, fileName, OPEN_RDONLY);
    }

    template < typename TValue, typename TSpec >
    inline bool open(String<TValue, TSpec> &string, const char *fileName, int openMode) {
    SEQAN_CHECKPOINT
        String<TValue, External< ExternalConfigLarge<> > > extString;
        if (!open(extString, fileName, openMode & ~OPEN_CREATE)) return false;
        assign(string, extString, Exact());
        return true;
    }
    template < typename TValue, typename TSpec >
    inline bool open(String<TValue, TSpec> &string, const char *fileName) {
    SEQAN_CHECKPOINT
        return open(string, fileName, OPEN_RDONLY);
    }

    template < typename THost, typename TSpec >
    inline bool open(Segment<THost, TSpec> &string, const char *fileName, int openMode) {
    SEQAN_CHECKPOINT
        String<typename Value<THost>::Type, External< ExternalConfigLarge<> > > extString;
        if (!open(extString, fileName, openMode & ~OPEN_CREATE)) return false;
        assign(string, extString, Exact());
        return true;
    }
    template < typename THost, typename TSpec >
    inline bool open(Segment<THost, TSpec> &string, const char *fileName) {
    SEQAN_CHECKPOINT
        return open(string, fileName, OPEN_RDONLY);
    }

#ifdef PLATFORM_CUDA
    template <typename TChar, typename TAlloc>
    inline bool open(thrust::device_vector<TChar, TAlloc> & me, const char *fileName, int openMode)
    {
        String<TChar> str;
        if (!open(str, fileName, openMode)) return false;
        assign(me, str);
        return true;
    }

    template <typename TChar, typename TAlloc>
    inline bool open(thrust::device_vector<TChar, TAlloc> & me, const char *fileName)
    {
        return open(me, fileName, OPEN_RDONLY);
    }
#endif

    // ATTENTION:
    // This implementation of open doesn't work with external memory StringSets (External<>, MMap<>)
    // If you need a persistent external StringSet you have to use a Owner<ConcatDirect<> > StringSet.
    template < typename TString, typename TSSSpec >
    inline bool open(StringSet<TString, TSSSpec> &multi, const char *fileName, int openMode) {
    SEQAN_CHECKPOINT
        char id[12]; // 2^32 has 10 decimal digits + 1 (0x00)
        unsigned i = 0;
        clear(multi);
        CharString name;
        while (true)
        {
            sprintf(id, ".%u", i);
            name = fileName;
            append(name, id);
            {
                resize(multi, i + 1);
                if (!open(multi[i], toCString(name), (openMode & ~OPEN_CREATE) | OPEN_QUIET))
                {
                    resize(multi, i);
                    break;
                }
            }
            ++i;
        }
        return i > 0;
    }

    template < typename TValue, typename TSpec, typename TSSSpec >
        inline bool open(StringSet<String<TValue, TSpec>, Dependent<TSSSpec> > &, const char *, int) {
        SEQAN_CHECKPOINT
        // Do nothing for dependent string sets
        return true;
    }

    template < typename TString, typename TSSSpec >
    inline bool open(StringSet<TString, Owner<ConcatDirect<TSSSpec> > > &multi, const char *fileName, int openMode) {
    SEQAN_CHECKPOINT
        CharString name;
        name = fileName;
        append(name, ".concat");
        if (!open(multi.concat, toCString(name), openMode | OPEN_QUIET)) return false;
        name = fileName;
        append(name, ".limits");
        if (!open(multi.limits, toCString(name), openMode | OPEN_QUIET) && !empty(multi.concat))
        {
            clear(multi);
            return false;
        }
        // limits file was just created
        if (empty(multi.limits))
            appendValue(multi.limits, 0);
        return true;
    }

    template < typename TValue, typename TSpec, typename TSSSpec>
    inline bool open(StringSet<String<TValue, TSpec>, TSSSpec> &multi, const char *fileName) {
    SEQAN_CHECKPOINT
        return open(multi, fileName, OPEN_RDONLY);
    }


//////////////////////////////////////////////////////////////////////////////
// save

    template <typename TValue>
    inline bool save(TValue const &val, const char *fileName, int openMode)
    {
        String<TValue, External< ExternalConfigLarge<> > > extString;
        if (!open(extString, fileName, openMode)) return false;
        clear(extString);
        appendValue(extString, val);
        return true;
    }

    template <typename TValue>
    inline bool save(TValue const &val, const char *fileName)
    {
        return save(val, fileName, OPEN_WRONLY | OPEN_CREATE);
    }

    template < typename TValue, typename TSpec >
    inline bool save(String<TValue, TSpec> const &string, const char *fileName, int openMode) {
    SEQAN_CHECKPOINT
//
//        if (length(string) == 0) return true;
        String<TValue, External< ExternalConfigLarge<> > > extString;
        if (!open(extString, fileName, openMode)) return false;
        assign(extString, string, Exact());
        return true;
    }
    template < typename TValue, typename TSpec >
    inline bool save(String<TValue, TSpec> const &string, const char *fileName) {
    SEQAN_CHECKPOINT
        return save(string, fileName, OPEN_WRONLY | OPEN_CREATE);
    }

    template < typename THost, typename TSpec >
    inline bool save(Segment<THost, TSpec> const &string, const char *fileName, int openMode) {
    SEQAN_CHECKPOINT
        if (length(string) == 0) return true;
        String<typename Value<THost>::Type, External< ExternalConfigLarge<> > > extString;
        if (!open(extString, fileName, openMode)) return false;
        assign(extString, string, Exact());
        return true;
    }
    template < typename THost, typename TSpec >
    inline bool save(Segment<THost, TSpec> const &string, const char *fileName) {
    SEQAN_CHECKPOINT
        return save(string, fileName, OPEN_WRONLY | OPEN_CREATE);
    }

    template < typename TString, typename TSSSpec>
    inline bool save(StringSet<TString, TSSSpec> const &multi, const char *fileName, int openMode) {
    SEQAN_CHECKPOINT
        if (length(multi) == 0) return true;
        char id[12]; // 2^32 has 10 decimal digits + 2 ('.' and 0x00)
        CharString name;
        for(unsigned i = 0; i < length(multi); ++i)
        {
            sprintf(id, ".%u", i);
            name = fileName;
            append(name, &(id[0]));
            if (!save(multi[i], toCString(name), openMode))
                return false;
        }
        return true;
    }

    template < typename TValue, typename TSpec, typename TSSSpec >
    inline bool save(StringSet<String<TValue, TSpec>, Dependent<TSSSpec> > const &, const char *, int) {
        SEQAN_CHECKPOINT
        // Do nothing for dependent string sets
        return true;
    }

    template < typename TString, typename TSSSpec >
    inline bool save(StringSet<TString, Owner<ConcatDirect<TSSSpec> > > const &multi, const char *fileName, int openMode) {
    SEQAN_CHECKPOINT
        CharString name;
        name = fileName;
        append(name, ".concat");
        if (!save(multi.concat, toCString(name), openMode)) return false;
        name = fileName;
        append(name, ".limits");
        if (!save(multi.limits, toCString(name), openMode)) return false;
        return true;
    }
    template < typename TValue, typename TSpec, typename TSSSpec>
    inline bool save(StringSet<String<TValue, TSpec>, TSSSpec> const &multi, const char *fileName) {
    SEQAN_CHECKPOINT
        return save(multi, fileName, OPEN_WRONLY | OPEN_CREATE);
    }

}

#endif
