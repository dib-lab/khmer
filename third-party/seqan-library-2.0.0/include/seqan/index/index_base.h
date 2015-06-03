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

#ifndef SEQAN_HEADER_INDEX_BASE_H
#define SEQAN_HEADER_INDEX_BASE_H

//#define SEQAN_TEST_INDEX

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////
// needful forward declarations

    // suffix array construction specs
    struct Skew3;
    struct Skew7;
    struct LarssonSadakane;
    struct ManberMyers;
    struct SAQSort;
    struct QGramAlg;

    // inverse suffix array construction specs
    template <typename TParallel>
    struct FromSortedSa{};

    // lcp table construction algorithms
    struct Kasai;
    struct KasaiOriginal;    // original, but more space-consuming algorithm

    // enhanced suffix array construction algorithms
    struct Childtab;
    struct Bwt;

/*!
 * @defgroup IndexFindAlgorithm Index Find Algorithm
 * @brief Tag to specify the index search algorithm.
 *
 * These tags can be used to specify the @link Finder#find @endlink algorithm for @link
 * Index @endlink based substring searches.
 *
 * @see Finder
 *
 * @tag IndexFindAlgorithm#FinderSTree
 * @brief Suffix tree search.
 *
 * Exact string matching using a suffix tree.
 *
 * @tag IndexFindAlgorithm#PizzaChiliFinder
 * @brief Finds an occurrence in a @link PizzaChiliIndex @endlink index.
 *
 * The actual algorithm used for searching depends on the @link PizzaChiliIndexTags @endlink used.
 *
 * @tag IndexFindAlgorithm#QGramFindLookup
 * @brief q-gram search. Finds q-grams in a @link IndexQGram @endlink index using the hash table.
 *
 * @tag IndexFindAlgorithm#FinderLcpe
 * @brief Binary search using lcp values.
 *
 *
 * Exact string matching using a suffix array binary search and a lcp-interval tree.
 *
 * @tag IndexFindAlgorithm#FinderMlr
 * @brief Binary search with mlr-heuristic.
 *
 * Exact string matching using a suffix array binary search with the mlr-heuristic.
 */
    // finder tags
    struct FinderMlr_;     // simple Suffix Array finder with mlr-heuristic
    struct FinderLcpe_;    // Suffix Array finder using an enhanced LCP-Table
    struct FinderSTree_;    // Suffix Array finder using an enhanced LCP-Table

    typedef Tag<FinderMlr_> const FinderMlr;
    typedef Tag<FinderLcpe_> const FinderLcpe;
    typedef Tag<FinderSTree_> const FinderSTree;

    template <typename TSpec = void>
    struct IndexEsa {};


// ----------------------------------------------------------------------------
// Metafunction DefaultIndexSpec
// ----------------------------------------------------------------------------
/*!
 * @mfn Index#DefaultIndexSpec
 * @headerfile <seqan/index.h>
 * @brief Default @link Index @endlink specialization type.
 *
 * @signature DefaultIndexSpec<TText>::Type;
 *
 * @tparam TText The given text type.
 *
 * @return TReturn Currently the return type is @link IndexEsa @endlink.
 *
 * @section Examples
 *
 * The following will define <tt>TIndex</tt> to be of type @link IndexEsa @endlink, because the
 * default index type for @link String @endlink or @link StringSet @endlink is @link IndexEsa @endlink.
 *
 * @code{.cpp}
 * typedef DefaultIndexSpec<String<Dna> >::Type TIndex;
 * @endcode
 */
    template < typename TObject >
    struct DefaultIndexSpec {
        typedef IndexEsa<> Type;
    };

// ----------------------------------------------------------------------------
// Metafunction DefaultIndexStringSpec
// ----------------------------------------------------------------------------
//NOTE(esiragusa): Deprecated in favor of StringSpec.

/*!
 * @mfn Index#DefaultIndexStringSpec
 * @headerfile <seqan/index.h>
 * @brief Default @link String @endlink specialization type of the @link Fibre
 *        @endlink of an @link Index @endlink.
 *
 * @signature DefaultIndexStringSpec<TIndex>::Type;
 *
 * @tparam TIndex An @link Index @endlink Type.
 *
 * @return TReturn If the underlying text is a @link String @endlink or a set of
 *                 Strings (see @link StringSet @endlink) the String's spec.
 *                 type is returned.
 *
 * @section Remarks
 *
 * Most of the @link Index @endlink fibres are strings. The @link String
 * @endlink specialization type is chosen by this meta-function.
 */

    template <typename TObject>
    struct DefaultIndexStringSpec : StringSpec<TObject> {};


//////////////////////////////////////////////////////////////////////////////
template <
        typename TObject,
        typename TSpec = typename DefaultIndexSpec<TObject>::Type >
    class Index;

/*!
 * @class Index
 * @headerfile <seqan/index.h>
 * @brief Indices are data structures which contain preprocessing data of a fixed text (or set of texts). In combination
 *        with a @link Finder @endlink or an @link VSTreeIterator @endlink it allows fast dictionary look-up and advanced computations.
 *
 * @signature class Index<TText[, TSpec]>;
 *
 * @tparam TSpec The index type.
 * @tparam TText The text type. Types: @link String @endlink, @link StringSet @endlink
 *
 *
 * An index contains various arrays or objects, also called fibres (see @link Fibre @endlink).
 * These fibres are created on demand depending on the requirements of an algorithm. To force the fibre creation you can
 * use the function @link Index#indexCreate @endlink.
 * The list of fibres is available in @link Fibre @endlink
 *
 * @section Examples
 *
 * The following code shows how to search for exact matches between the reference "tobeornottobe" and the
 * pattern "to" with the means of a Finder.
 *
 * @include demos/index/index_finder.cpp
 *
 * The result is as follows
 *
 * @include demos/index/index_finder.cpp.stdout
 *
 * This code shows how an index can be used with iterators to achieve a pre-order tree like traversal
 * in DFS of the text "tobeornottobe". In order to do so a Top-Down History iterator is used.
 *
 * @include demos/index/index_iterator.cpp
 *
 * The result is as follows
 *
 * @include demos/index/index_iterator.cpp.stdout
 *
 * Note that you can also use specialized iterators such as:
 *
 * @code{.cpp}
 * Iterator<TIndex, TopDown<ParentLinks<PreOrder> > >::Type
 * @endcode
 *
 * or
 *
 * @code{.cpp}
 * Iterator<TIndex, TopDown<ParentLinks<PostOrder> > >::Type
 * @endcode
 *
 * You can achieve a post-order traversal like this:
 *
 * @snippet demos/index/index_iterator_short.cpp iteration
 *
 * The result is as follows
 *
 * @include demos/index/index_iterator_short.cpp.stdout
 */

    template <typename TObject, typename TSpec>
    struct Host< Index<TObject, TSpec> > {
        typedef TObject Type;
    };

    template <typename TObject, typename TSpec>
    struct Spec< Index<TObject, TSpec> > {
        typedef TSpec Type;
    };
/*!
 * @mfn Fibre
 * @headerfile <seqan/index.h>
 * @brief Type of a specific container member (fibre).
 *
 * @signature Fibre<TObject, TSpec>::Type;
 *
 * @tparam TSpec Tag to specify the fibre. Types: @link IndexEsaFibres @endlink, @link FMIndexFibres @endlink, @link
 *               QGramIndexFibres @endlink, @link WOTDIndexFibres @endlink, @link WaveletTreeFibres @endlink, @link
 *               RightArrayBinaryTreeFibres @endlink, @link RankDictionaryFibres @endlink, @link
 *               SentinelRankDictionaryFibres @endlink
 * @tparam TObject The container type. Types: @link Index @endlink, @link RankDictionary @endlink, @link SparseString
 *                 @endlink, @link CompressedSA @endlink
 *
 * @return Type Fibre type.
 *
 * @section Naming
 *
 * Some containers, such as the @link Index @endlink or the @link RankDictionary @endlink, can be seen as a collection
 * of tables.  However, each table alone is just a collection of information.  They only become powerful if used
 * together.  Therefore, a more appropriate label for the tables is fibre, like the fibres of a rope.
 *
 * In addition, sometimes a fibre can be a single value and calling it a table would be misleading.
 */

    // meta function to get the type of a bundle fibre
    template < typename TIndex, typename TSpec >
    struct Fibre {
        typedef String< typename Size<TIndex>::Type > Type;
    };

    template < typename TIndex, typename TSpec >
    struct Fibre<TIndex const, TSpec> {
        typedef typename Fibre<TIndex, TSpec>::Type const Type;
    };

    struct FibreRecord {
        unsigned    id;
        void*        ptr;
        bool        owner;
    };

    // less function to search in sorted list for fibre id
    struct FibreLess: public std::binary_function<FibreRecord, unsigned, bool>
    {    // functor for operator>
        inline bool operator()(FibreRecord const & _Left, unsigned const Right_) const
        {    // apply operator> to operands
            return (_Left.id < Right_);
        }
    };

//////////////////////////////////////////////////////////////////////////////
/*!
 * @mfn Index#DefaultIndexCreator
 * @headerfile <seqan/index.h>
 * @note Advanced functionality, not commonly used.
 * @brief Default algorithm to create a demanded and not yet existing @link Fibre @endlink.
 *
 * @signature DefaultIndexCreator<TIndex, TFibre>::Type;
 *
 * @tparam TIndex An @link Index @endlink Type.
 * @tparam TFibre A tag specifying the fibre (e.g. @link IndexEsaFibres#EsaSA
 *                @endlink).
 *
 * @return Type A tag specifying the default algorithm to create the fibre with.
 */

    template < typename TIndex, typename TFibre >
    struct DefaultIndexCreator {
        typedef Default Type;
    };


//////////////////////////////////////////////////////////////////////////////

    template <
        typename TSA,
        typename TText,
        typename TAlgSpec >
    struct SACreatorRandomAccess_
    {
        typedef typename AllowsFastRandomAccess<TSA>::Type   TRandomSA;
        typedef typename AllowsFastRandomAccess<TText>::Type TRandomText;
        typedef typename And<TRandomText,TRandomSA>::Type Type;
    };

    template <
        typename TLCP,
        typename TText,
        typename TSA,
        typename TAlgSpec >
    struct LcpCreatorRandomAccess_
    {
        typedef typename AllowsFastRandomAccess<TText>::Type TRandomText;
        typedef typename AllowsFastRandomAccess<TLCP>::Type  TRandomLCP;
        typedef typename AllowsFastRandomAccess<TSA>::Type   TRandomSA;
        typedef typename And<TRandomLCP, typename And<TRandomText,TRandomSA>::Type>::Type Type;
    };


//////////////////////////////////////////////////////////////////////////////

/*
    template < typename TSpec = void >
    truct Bundle {
        typedef std::vector<FibreRecord>    TFibreRecords;
        TFibreRecords                        fibres;
    };

    template < typename TBundleSpec, typename TFibreSpec >
    inline FibreRecord& getRecord(Bundle<TBundleSpec> &bundle, TFibreSpec const) {
        unsigned id = (unsigned)ClassIdentifier_<TFibreSpec>::getID();

        typename Bundle<TBundleSpec>::TFibreRecords::iterator first = lower_bound(bundle.fibres.begin(), bundle.fibres.end(), id, FibreLess());
        if (!first->id != id) {
            FibreRecord rec;
            rec.id = id;
            rec.ptr = NULL;
            rec.owner = true;
            bundle.fibres.insert(first, rec);
        } else
            return *first;
    }

    template < typename TBundleSpec, typename TFibreSpec >
    inline typename Fibre<Bundle<TBundleSpec>, TFibreSpec>::Type & getFibre(Bundle<TBundleSpec> &bundle, TFibreSpec const) {
        typedef typename Fibre<Bundle<TBundleSpec>, TFibreSpec>::Type Type;
        unsigned id = (unsigned)ClassIdentifier_<TFibreSpec>::getID();

        FibreRecord &rec = getRecord(bundle, TFibreSpec());
        if (!rec.ptr)
            rec.ptr = new Type();
        return *reinterpret_cast<Type*>(rec.ptr);
    }

    template < typename TBundleSpec, typename TFibreSpec >
    inline typename Fibre<Bundle<TBundleSpec>, TFibreSpec>::Type const & getFibre(Bundle<TBundleSpec> const &bundle, TFibreSpec const) {
        typedef typename Fibre<Bundle<TBundleSpec>, TFibreSpec>::Type Type;
        unsigned id = (unsigned)ClassIdentifier_<TFibreSpec>::getID();

        FibreRecord &rec = getRecord(bundle, TFibreSpec());
        return *reinterpret_cast<Type*>(rec.ptr);
    }
*/

//////////////////////////////////////////////////////////////////////////////
// various fibre specs for enhanced suffix arrays

    struct FibreText_;      // Original text. Can be a String or a StringSet
    struct FibreRawText_;   // Concatenation of the strings above
    struct FibreSA_;        // suffix array (of raw text with virtual $-delimiters) with Pair entries
    struct FibreRawSA_;     // suffix array with integer entries
    struct FibreIsa_;       // inverse suffix array.
    struct FibreSae_;       // suffix array reordered in a b-tree
    struct FibreLcp_;       // lcp table of raw text
    struct FibreLcpe_;      // lcp interval tree
    struct FibreChildtab_;  // childtab (Kurtz et al.) of raw text
    struct FibreBwt_;       // burrows wheeler table of raw text

    typedef Tag<FibreText_> const      FibreText;
    typedef Tag<FibreRawText_> const   FibreRawText;
    typedef Tag<FibreSA_> const        FibreSA;
    typedef Tag<FibreRawSA_> const     FibreRawSA;
    typedef Tag<FibreIsa_> const       FibreIsa;
    typedef Tag<FibreSae_> const       FibreSae;
    typedef Tag<FibreLcp_> const       FibreLcp;
    typedef Tag<FibreLcpe_> const      FibreLcpe;
    typedef Tag<FibreChildtab_> const  FibreChildtab;
    typedef Tag<FibreBwt_> const       FibreBwt;

//////////////////////////////////////////////////////////////////////////////
/*!
 * @mfn SAValue
 * @headerfile <seqan/index.h>
 * @brief The default alphabet type of a suffix array, i.e. the type to store a
 *        position of a string or string set.
 *
 * @signature SAValue<TObject>::Type;
 *
 * @tparam TObject A string, string set, or index type. Types: @link String @endlink, @link StringSet @endlink,
 *                 @link Index @endlink
 * @return TReturn A type to store a position.If <tt>TObject</tt> is a @link String @endlink, it is a single integer
 *                 value. By default this is the @link Size @endlink type of <tt>TObject</tt>.If <tt>TObject</tt> is a
 *                 @link StringSet @endlink, it could be a single integer too (called global position, see @link
 *                 ConcatDirectStringSet @endlink) or a @link Pair @endlink (called local position, see
 *                 @link OwnerStringSet @endlink).  Currently SeqAn defaults to a local position for @link StringSet
 *                 @endlink classes (index_base.h).
 *
 * SAValue is the return type of various functions, e.g. @link Finder#position @endlink for the Index @link Finder @endlink
 * class, @link VSTreeIterator#getOccurrence @endlink, @link VSTreeIterator#getOccurrences @endlink, @link IndexQGram#getOccurrences
 * @endlink etc. You should always use the type of this meta-function to store the return values. If you want to write
 * algorithms for both variants (local and global positions) you should use the functions @link TextConcept#posLocalize @endlink,
 * @link TextConcept#posGlobalize @endlink, @link TextConcept#getSeqNo @endlink and @link TextConcept#getSeqOffset @endlink.
 *
 * @note If <tt>TObject</tt> is an @link Index @endlink, @link Position @endlink returns the same value as <tt>SAValue</tt>.
 *       You can change the position type of an index by overloading <tt>SAValue</tt>, not @link Position @endlink.
 *
 * @section Examples
 *
 * The following code snippet demonstrates the usage of @link SAValue @endlink.
 * @code{.cpp}
 * template < typename TString, typename TSpec >
 * struct SAValue< StringSet<TString, TSpec> > {
 *     typedef Pair<
 *         typename Size< StringSet<TString, TSpec> >::Type,
 *         typename SAValue<TString>::Type,
 *         Pack
 *     > Type;
 * };
 * @endcode
 * @see orderOccurrences
 */

    template <typename TObject>
    struct SAValue:
        Position<TObject> {};

    template <typename TObject>
    struct SAValue<TObject const>:
        SAValue<TObject> {};

    // to speed up sequence number computation
    // we use a pair of seqNo and localPosition
    template <typename TString, typename TSpec>
    struct SAValue<StringSet<TString, TSpec> > :
        StringSetPosition<StringSet<TString, TSpec> > {};

/*
    template < typename TString, typename TSpec >
    struct SAValue< StringSet<TString, TSpec> > {
        typedef Pair<
            typename Size< StringSet<TString, TSpec> >::Type,
            typename SAValue<TString>::Type,
            BitPacked<2,30>                            // max. 4 sequences
        > Type;                                        // max. 2^30 characters each
    };
*/
    template < typename TText, typename TSpec >
    struct SAValue< Index<TText, TSpec> >:
        SAValue<TText> {};

    template <typename TText, typename TSpec>
    struct StringSpec<Index<TText, TSpec> > : StringSpec<TText> {};

//////////////////////////////////////////////////////////////////////////////
// value and size type of an index

    template < typename TText, typename TSpec >
    struct Value< Index<TText, TSpec> > {
        typedef typename Value<
            typename Fibre< Index<TText, TSpec>, FibreRawText>::Type
        >::Type Type;
    };

/*!
 * @mfn Index#Size
 * @headerfile <seqan/index.h>
 * @brief Returns the size type of an Index.
 *
 * @signature Size<TIndex>::Type;
 *
 * @tparam TIndex The Index specialization.
 *
 * @return Type The resulting size type of the index.
 */

    template < typename TText, typename TSpec >
    struct Size< Index<TText, TSpec> > {
        typedef typename LengthSum<TText>::Type Type;
    };

    template < typename TText, typename TSpec >
    struct Position< Index<TText, TSpec> >:
        SAValue< Index<TText, TSpec> > {};

//////////////////////////////////////////////////////////////////////////////
// infix of an index

    template < typename TText, typename TSpec >
    struct Infix< Index<TText, TSpec> >:
        public Infix<TText> {};

    template < typename TText, typename TSpec >
    struct Infix< Index<TText, TSpec> const >:
        public Infix<TText> {};

//////////////////////////////////////////////////////////////////////////////
// default table type

    template < typename TObject, typename TSpec, typename TFibre >
    struct Fibre< Index<TObject, TSpec>, Tag<TFibre> const > {
        typedef String<
            typename Size< Index<TObject, TSpec> >::Type,
            typename StringSpec< Index<TObject, TSpec> >::Type
        > Type;
    };

//////////////////////////////////////////////////////////////////////////////
// original text

    template < typename TText, typename TSpec >
    struct Fibre< Index<TText, TSpec>, FibreText> {
        typedef TText Type;
    };

//////////////////////////////////////////////////////////////////////////////
// type of the text member

    template <typename TText, typename TSpec>
    struct Member<Index<TText, TSpec>, FibreText>
    {
        typedef Holder<typename Fibre<Index<TText, TSpec>, FibreText>::Type>    Type;
    };

//////////////////////////////////////////////////////////////////////////////
// concatenated text

    template < typename TText, typename TSpec >
    struct Fibre< Index<TText, TSpec>, FibreRawText> {
        typedef typename Concatenator<TText>::Type Type;
    };

//////////////////////////////////////////////////////////////////////////////
// suffix array type

    template < typename TText, typename TSpec >
    struct Fibre< Index<TText, TSpec>, FibreSA> {
        typedef String<
            typename SAValue< Index<TText, TSpec> >::Type,
            typename StringSpec< Index<TText, TSpec> >::Type
        > Type;
    };

//////////////////////////////////////////////////////////////////////////////
// lcp type

    template <typename TText, typename TSSetSpec, typename TSpec>
    struct Fibre<Index<StringSet<TText, TSSetSpec>, TSpec>, FibreLcp>
    {
        typedef String<
            typename Size<TText>::Type,
            typename StringSpec<Index<StringSet<TText, TSSetSpec>, TSpec> >::Type
        > Type;
    };

//////////////////////////////////////////////////////////////////////////////
// globalize functor

    template <typename InType, typename TLimitsString, typename Result = typename Value<TLimitsString>::Type>
    struct FunctorGlobalize : public std::unary_function<InType,Result>
    {
        TLimitsString const * limits;

        FunctorGlobalize() : limits() {}
        FunctorGlobalize(TLimitsString const &_limits) : limits(&_limits) {}

        inline Result operator()(InType const &x) const
        {
            return posGlobalize(x, *limits);
        }
    };

    template <typename InType, typename Result>
    struct FunctorGlobalize<InType, Nothing, Result> : public std::unary_function<InType,InType>
    {
        FunctorGlobalize() {}
        FunctorGlobalize(Nothing const &) {}

        inline InType operator()(InType const &x) const
        {
            return x;
        }
    };

//////////////////////////////////////////////////////////////////////////////
// raw suffix array contains integer offsets relative to raw text
/*
    template < typename TString, typename TSpec >
    struct Fibre< Index<TString, TSpec>, FibreRawSA>:
        public Fibre< Index<TString, TSpec> const, FibreSA> {};

    template < typename TString, typename TSSetSpec, typename TSpec >
    struct Fibre< Index<StringSet<TString, TSSetSpec>, TSpec>, FibreRawSA>
    {
        typedef Index< StringSet<TString, TSSetSpec>, TSpec> TIndex;
        typedef ModifiedString<
            typename Fibre<TIndex, FibreSA>::Type,
            ModView< FunctorGlobalize<
                typename Value< typename Fibre<TIndex, FibreSA>::Type >::Type,
                typename StringSetLimits<TString>::Type >
            >
        > Type;
    };
*/
    template < typename TText, typename TSpec >
    struct Fibre< Index<TText, TSpec>, FibreRawSA>
    {
        typedef Index<TText, TSpec> TIndex;
        typedef ModifiedString<
            typename Fibre<TIndex, FibreSA>::Type,
            ModView< FunctorGlobalize<
                typename Value< typename Fibre<TIndex, FibreSA>::Type >::Type,
                typename StringSetLimits<TText>::Type >
            >
        > Type;
    };

//////////////////////////////////////////////////////////////////////////////
// default burrows-wheeler table

    template < typename TText, typename TSpec >
    struct Fibre< Index<TText, TSpec>, FibreBwt> {
        typedef String <
            typename Value< Index<TText, TSpec> >::Type,
            typename StringSpec< Index<TText, TSpec> >::Type
        > Type;
    };


//////////////////////////////////////////////////////////////////////////////
// default fibre creators

    template < typename TText, typename TSpec >
    struct DefaultIndexCreator<Index<TText, TSpec>, FibreSA> {
        typedef Skew7 Type;                            // standard suffix array creator is skew7
    };

    template < typename TText, typename TSpec >
    struct DefaultIndexCreator<Index<TText, TSpec>, FibreIsa> {
        typedef FromSortedSa<Serial> Type;
    };

    template < typename TText, typename TSpec >
    struct DefaultIndexCreator<Index<TText, TSpec>, FibreLcp> {
        typedef Kasai Type;
    };

    template < typename TText, typename TSpec >
    struct DefaultIndexCreator<Index<TText, TSpec>, FibreBwt> {
        typedef Bwt Type;
    };

    template < typename TText, typename TSpec >
    struct DefaultIndexCreator<Index<TText, TSpec>, FibreChildtab> {
        typedef Childtab Type;
    };

    template <typename TText, typename TSpec>
    inline typename Member<Index<TText, TSpec>, FibreText>::Type &
    _dataHost(Index<TText, TSpec> & index) {
        return index.text;
    }
    template <typename TText, typename TSpec>
    inline typename Member<Index<TText, TSpec> const, FibreText>::Type &
    _dataHost(Index<TText, TSpec> const &index) {
        return index.text;
    }

//////////////////////////////////////////////////////////////////////////////
/*!
 * @fn Index#getFibre
 * @headerfile <seqan/index.h>
 * @brief Returns a specific fibre of a container.
 *
 * @signature TFibre getFibre(index, fibreTag);
 *
 * @param[in] index    The container holding the fibre.
 * @param[in] fibreTag A tag that identifies the @link Fibre @endlink. Types: @link IndexEsaFibres Index Esa Fibres
 *                     @endlink, @link FMIndexFibres FM Index Fibres @endlink, @link WOTDIndexFibres Index Wotd Fibres @endlink,
 *                     and @link QGramIndexFibres Index QGram Fibres @endlink.
 *
 * @return TFibre A reference to the @link Fibre @endlink object.
 *
 * @section Examples
 *
 * The following example shows how to search for a pattern in a string.
 *
 * @include demos/index/index_getOccurrences_getFrequency_range_getFibre.cpp
 *
 * The output is as follows:
 *
 * @include demos/index/index_getOccurrences_getFrequency_range_getFibre.cpp.stdout
 *
 * @see Fibre
 */

    template <typename TText, typename TSpec>
    inline typename Fibre<Index<TText, TSpec>, FibreText>::Type &
    getFibre(Index<TText, TSpec> &index, FibreText) {
        return value(index.text);
    }
    template <typename TText, typename TSpec>
    inline typename Fibre<Index<TText, TSpec> const, FibreText>::Type &
    getFibre(Index<TText, TSpec> const &index, FibreText) {
        return value(index.text);
    }

//////////////////////////////////////////////////////////////////////////////

    template <typename TText, typename TSpec>
    inline typename Fibre<Index<TText, TSpec>, FibreRawText>::Type &
    getFibre(Index<TText, TSpec> &index, FibreRawText) {
        return concat(getFibre(index, FibreText()));
    }
    template <typename TText, typename TSpec>
    inline typename Fibre<Index<TText, TSpec> const, FibreRawText>::Type &
    getFibre(Index<TText, TSpec> const &index, FibreRawText) {
        return concat(getFibre(index, FibreText()));
    }

//////////////////////////////////////////////////////////////////////////////

    template <typename TText, typename TSpec>
    SEQAN_HOST_DEVICE inline typename Fibre<Index<TText, TSpec>, FibreSA>::Type &
    getFibre(Index<TText, TSpec> &index, FibreSA) {
        return index.sa;
    }
    template <typename TText, typename TSpec>
    SEQAN_HOST_DEVICE inline typename Fibre<Index<TText, TSpec> const, FibreSA>::Type &
    getFibre(Index<TText, TSpec> const &index, FibreSA) {
        return index.sa;
    }

//////////////////////////////////////////////////////////////////////////////

    template <typename TText, typename TSpec>
    SEQAN_HOST_DEVICE inline typename Fibre<Index<TText, TSpec>, FibreIsa>::Type &
    getFibre(Index<TText, TSpec> &index, FibreIsa) {
        return index.isa;
    }
    template <typename TText, typename TSpec>
    SEQAN_HOST_DEVICE inline typename Fibre<Index<TText, TSpec> const, FibreIsa>::Type &
    getFibre(Index<TText, TSpec> const &index, FibreIsa) {
        return index.isa;
    }

//////////////////////////////////////////////////////////////////////////////
/*
    template <typename TText, typename TSpec>
    inline typename Fibre<Index<TText, TSpec>, FibreRawSA>::Type const &
    getFibre(Index<TText, TSpec> &index, FibreRawSA) {
        return indexSA(index);
    }

    template <typename TString, typename TSSetSpec, typename TSpec>
    inline typename Fibre<Index<StringSet<TString, TSSetSpec>, TSpec>, FibreRawSA>::Type
    getFibre(Index<StringSet<TString, TSSetSpec>, TSpec> &index, FibreRawSA)
    {
        typedef Index< StringSet<TString, TSSetSpec>, TSpec> TIndex;

        typedef FunctorGlobalize<
            typename Value< typename Fibre<TIndex, FibreSA>::Type >::Type,
            typename StringSetLimits<StringSet<TString, TSSetSpec> >::Type
        > TFunctor;

        typedef ModifiedString<
            typename Fibre<Index<StringSet<TString, TSSetSpec>, TSpec>, FibreSA>::Type,
            ModView< TFunctor >
        > ModString;

        return ModString(indexSA(index), TFunctor(stringSetLimits(indexText(index))));
    }
*/

    template <typename TText, typename TSpec>
    inline typename Fibre<Index<TText, TSpec>, FibreRawSA>::Type
    getFibre(Index<TText, TSpec> &index, FibreRawSA)
    {
        typedef Index<TText, TSpec> TIndex;

        typedef FunctorGlobalize<
            typename Value< typename Fibre<TIndex, FibreSA>::Type >::Type,
            typename StringSetLimits<TText>::Type
        > TFunctor;

        typedef ModifiedString<
            typename Fibre<Index<TText, TSpec>, FibreSA>::Type,
            ModView< TFunctor >
        > ModString;

        return ModString(indexSA(index), TFunctor(stringSetLimits(indexText(index))));
    }
//////////////////////////////////////////////////////////////////////////////

    template <typename TText, typename TSpec>
    SEQAN_HOST_DEVICE inline typename Fibre<Index<TText, TSpec>, FibreLcp>::Type &
    getFibre(Index<TText, TSpec> &index, FibreLcp) {
        return index.lcp;
    }
    template <typename TText, typename TSpec>
    SEQAN_HOST_DEVICE inline typename Fibre<Index<TText, TSpec> const, FibreLcp>::Type &
    getFibre(Index<TText, TSpec> const &index, FibreLcp) {
        return index.lcp;
    }

//////////////////////////////////////////////////////////////////////////////

    template <typename TText, typename TSpec>
    inline typename Fibre<Index<TText, TSpec>, FibreLcpe>::Type &
    getFibre(Index<TText, TSpec> &index, FibreLcpe) {
        return index.lcpe;
    }
    template <typename TText, typename TSpec>
    inline typename Fibre<Index<TText, TSpec> const, FibreLcpe>::Type &
    getFibre(Index<TText, TSpec> const &index, FibreLcpe) {
        return index.lcpe;
    }

//////////////////////////////////////////////////////////////////////////////

    template <typename TText, typename TSpec>
    SEQAN_HOST_DEVICE inline typename Fibre<Index<TText, TSpec>, FibreChildtab>::Type &
    getFibre(Index<TText, TSpec> &index, FibreChildtab) {
        return index.childtab;
    }
    template <typename TText, typename TSpec>
    SEQAN_HOST_DEVICE inline typename Fibre<Index<TText, TSpec> const, FibreChildtab>::Type &
    getFibre(Index<TText, TSpec> const &index, FibreChildtab) {
        return index.childtab;
    }

//////////////////////////////////////////////////////////////////////////////

    template <typename TText, typename TSpec>
    SEQAN_HOST_DEVICE inline typename Fibre<Index<TText, TSpec>, FibreBwt>::Type &
    getFibre(Index<TText, TSpec> &index, FibreBwt) {
        return index.bwt;
    }
    template <typename TText, typename TSpec>
    SEQAN_HOST_DEVICE inline typename Fibre<Index<TText, TSpec> const, FibreBwt>::Type &
    getFibre(Index<TText, TSpec> const &index, FibreBwt) {
        return index.bwt;
    }

//////////////////////////////////////////////////////////////////////////////
/*!
 * @fn Index#length
 * @headerfile <seqan/index.h>
 * @brief Returns the number of characters in the underlying text of the index.
 *
 * @signature TSize length(index);
 *
 * @param[in] index An index of a text.
 *
 * @return TSize Returns the number of characters in the raw underlying text of the
 *               index with TSize being the result of the @link Size @endlink meta-function
 *               of @link Index @endlink.
 *
 * @section Examples
 *
 * The following example shows how to count characters of an index, determine the number of sequences involved and how
 * to search for a pattern.
 *
 * @include demos/index/index_length_countSequences.cpp
 *
 * The output is as follows:
 *
 * @include demos/index/index_length_countSequences.cpp.stdout
 */

    template <typename TText, typename TSpec>
    SEQAN_HOST_DEVICE inline typename Size<Index<TText, TSpec> >::Type
    length(Index<TText, TSpec> const &index) {
        return length(indexRawText(index));
    }

//////////////////////////////////////////////////////////////////////////////
/*!
 * @fn Index#countSequences
 * @headerfile <seqan/index.h>
 * @brief Return the number of sequences in an index' underlying text.
 *
 * @signature TSize countSequences(index);
 *
 * @param[in] index The index to return the number of sequences of.
 *
 * @return TSize The number of sequences in the index' underlying text with TSize being the result @link Size @endlink.
 *
 * @section Examples
 *
 * The following example shows how to count characters of an index, determine the number of sequences involved and how
 * to search for a pattern.
 *
 * @include demos/index/index_length_countSequences.cpp
 *
 * The output is as follows:
 *
 * @code{.output}
 * The text has 25 characters and consists of 2 sequences.
 * Hit at position: < 1 , 2 >
 * Hit at position: < 0 , 0 >
 * @endcode
 */

    template <typename TText, typename TSpec>
    inline typename Size<TText>::Type
    countSequences(Index<TText, TSpec> const &index) {
        return countSequences(indexText(index));
    }

//////////////////////////////////////////////////////////////////////////////
// TODO(singer): Since this is a public function it should be documented
    template <typename TText, typename TSpec>
    struct GetSequenceByNo< Index<TText, TSpec> >
    {
        typedef typename GetSequenceByNo<TText>::Type Type;
    };

    template <typename TText, typename TSpec>
    struct GetSequenceByNo< Index<TText, TSpec> const>
    {
        typedef typename GetSequenceByNo<TText const>::Type Type;
    };

//////////////////////////////////////////////////////////////////////////////
// TODO(singer): Since this is a public function it should be documented
    template <typename TSeqNo, typename TText, typename TSpec>
    SEQAN_HOST_DEVICE inline typename GetSequenceByNo< Index<TText, TSpec> >::Type
    getSequenceByNo(TSeqNo seqNo, Index<TText, TSpec> &index)
    {
        return getSequenceByNo(seqNo, indexText(index));
    }

    template <typename TSeqNo, typename TText, typename TSpec>
    SEQAN_HOST_DEVICE inline typename GetSequenceByNo< Index<TText, TSpec> const>::Type
    getSequenceByNo(TSeqNo seqNo, Index<TText, TSpec> const &index)
    {
        return getSequenceByNo(seqNo, indexText(index));
    }

//////////////////////////////////////////////////////////////////////////////
// TODO(singer): Since this is a public function it should be documented
    template <typename TSeqNo, typename TText, typename TSpec>
    SEQAN_HOST_DEVICE inline typename Size<Index<TText, TSpec> >::Type
    sequenceLength(TSeqNo seqNo, Index<TText, TSpec> const &index) {
        return sequenceLength(seqNo, indexText(index));
    }

//////////////////////////////////////////////////////////////////////////////
// TODO(singer): Since this is a public function it should be documented
    template <typename TPos, typename TText, typename TSpec>
    SEQAN_HOST_DEVICE inline typename Size<Index<TText, TSpec> >::Type
    suffixLength(TPos pos, Index<TText, TSpec> const &index)
    {
        return length(indexText(index)) - pos;
    }

    template <typename TPos, typename TString, typename TSSetSpec, typename TSpec>
    SEQAN_HOST_DEVICE inline typename Size<Index<StringSet<TString, TSSetSpec>, TSpec> >::Type
    suffixLength(TPos pos, Index<StringSet<TString, TSSetSpec>, TSpec> const &index)
    {
        typename StringSetLimits<StringSet<TString, TSSetSpec> >::Type const &limits = stringSetLimits(index);
        return sequenceLength(getSeqNo(pos, limits), index) - getSeqOffset(pos, limits);
    }


//////////////////////////////////////////////////////////////////////////////
/*!
 * @fn IndexEsa#textAt
 * @headerfile <seqan/index.h>
 * @brief Shortcut for <tt>value(indexText(..), ..)</tt>.
 *
 * @signature TValue textAt(position, index);
 *
 * @param[in] index The @link Index @endlink object holding the fibre.
 * @param[in] position A position in the fibre on which the value should be accessed.
 *
 * @return TValue A reference or proxy to the value in the text.
 */

    template <typename TPos, typename TIndex>
    inline typename Reference<typename Fibre<TIndex, FibreRawText>::Type>::Type
    textAt(TPos i, TIndex &index) {
        return value(getFibre(index, FibreRawText()), i);
    }
    template <typename TPos, typename TString, typename TSSetSpec, typename TSpec>
    inline typename Reference<typename Fibre< Index< StringSet<TString, TSSetSpec>, TSpec>, FibreRawText>::Type>::Type
    textAt(TPos i, Index< StringSet<TString, TSSetSpec>, TSpec> &index) {
        return value(getFibre(index, FibreRawText()), posGlobalize(i, stringSetLimits(index)));
    }
    template <typename TPos, typename TString, typename TSSetSpec, typename TSpec>
    inline typename Reference<typename Fibre< Index< StringSet<TString, TSSetSpec>, TSpec> const, FibreRawText>::Type>::Type
    textAt(TPos i, Index< StringSet<TString, TSSetSpec>, TSpec> const &index) {
        return value(getFibre(index, FibreRawText()), posGlobalize(i, stringSetLimits(index)));
    }
    template <typename TPos, typename TString, typename TSpec>
    inline typename Reference<typename Fibre< Index< StringSet<TString, Owner<Default> >, TSpec>, FibreRawText>::Type>::Type
    textAt(TPos i, Index< StringSet<TString, Owner<Default> >, TSpec> &index) {
        Pair <
            typename Size< StringSet<TString, Owner<Default> > >::Type,
            typename Size< TString >::Type > locPos;
        posLocalize(locPos, i, stringSetLimits(index));
        return value(value(getFibre(index, FibreText()), getValueI1(locPos)), getValueI2(locPos));
    }
    template <typename TPos, typename TString, typename TSpec>
    inline typename Reference<typename Fibre< Index< StringSet<TString, Owner<Default> >, TSpec> const, FibreRawText>::Type>::Type
    textAt(TPos i, Index< StringSet<TString, Owner<Default> >, TSpec> const &index) {
        Pair <
        typename Size< StringSet<TString, Owner<Default> > >::Type,
        typename Size< TString >::Type > locPos;
        posLocalize(locPos, i, stringSetLimits(index));
        return value(value(getFibre(index, FibreText()), getValueI1(locPos)), getValueI2(locPos));
    }

//////////////////////////////////////////////////////////////////////////////
// infix

    template <typename TText, typename TSpec, typename TPosBegin, typename TPosEnd>
    inline typename Infix<TText>::Type
    infix(Index<TText, TSpec> &index, TPosBegin pos_begin, TPosEnd pos_end)
    {
        return infix(indexText(index), pos_begin, pos_end);
    }

    template <typename TText, typename TSpec, typename TPosBegin, typename TPosEnd>
    inline typename Infix<TText>::Type
    infix(Index<TText, TSpec> const &index, TPosBegin pos_begin, TPosEnd pos_end)
    {
        return infix(indexText(index), pos_begin, pos_end);
    }

//////////////////////////////////////////////////////////////////////////////
/*!
 * @fn Index#rawtextAt
 * @headerfile <seqan/index.h>
 * @brief Shortcut for <tt>value(indexRawText(..), ..)</tt>.
 *
 * @signature TValue rawtextAt(position, index);
 *
 * @param[in] index    The @link Index @endlink object. Types: @link Index @endlink
 * @param[in] position A position in the array on which the value should be accessed.
 *
 * @return TValue A reference or proxy to the value. The type is the result if the meta-function @link Reference
 *                @endlink.
 *
 * @note The result of this function when used on an Index&lt;TText, FMIndex&lt;TOccSpec, Compress&gt; &gt; is not
 *       defined.
 */

    template <typename TPos, typename TIndex>
    inline typename Reference<typename Fibre<TIndex, FibreRawText>::Type>::Type rawtextAt(TPos i, TIndex &index) {
        return value(getFibre(index, FibreRawText()), i);
    }
    template <typename TPos, typename TIndex>
    inline typename Reference<typename Fibre<TIndex const, FibreRawText>::Type>::Type rawtextAt(TPos i, TIndex const &index) {
        return value(getFibre(index, FibreRawText()), i);
    }

//////////////////////////////////////////////////////////////////////////////
/*!
 * @fn IndexEsa#saAt
 * @headerfile <seqan/index.h>
 * @brief Shortcut for <tt>value(indexSA(..), ..)</tt>.
 *
 * @signature TValue saAt(position, index);
 *
 * @param[in] index The @link Index @endlink object holding the fibre.
 * @param[in] position A position in the array on which the value should be accessed.
 *
 * @return TValue A reference or proxy to the value in the suffix array.
 */

    template <typename TPos, typename TIndex>
    SEQAN_HOST_DEVICE inline typename Reference<typename Fibre<TIndex, FibreSA>::Type>::Type saAt(TPos i, TIndex &index) {
        return value(getFibre(index, FibreSA()), i);
    }
    template <typename TPos, typename TIndex>
    SEQAN_HOST_DEVICE inline typename Reference<typename Fibre<TIndex const, FibreSA>::Type>::Type saAt(TPos i, TIndex const &index) {
        return value(getFibre(index, FibreSA()), i);
    }

//////////////////////////////////////////////////////////////////////////////
// TODO(weese): I disabled the doc, as we don't want to encourage users to use it

/*
 * @fn IndexEsa#rawsaAt
 * @headerfile <seqan/index.h>
 * @note Advanced functionality, not commonly used.
 * @brief Shortcut for <tt>value(indexRawSA(..), ..)</tt>.
 *
 * @signature TValue rawsaAt(position, index);
 *
 * @param[in] index    The @link Index @endlink object holding the fibre.
 * @param[in] position A position in the array on which the value should be accessed.
 *
 * @return TValue A reference or proxy to the value.  To be more precise, a reference to a position
 *                containing a value of type @link SAValue @endlink is returned (or a proxy).
 */

    template <typename TPos, typename TIndex>
    inline typename Value<typename Fibre<TIndex const, FibreRawSA>::Type>::Type rawsaAt(TPos i, TIndex const &index) {
        return posGlobalize(saAt(i, index), stringSetLimits(indexText(index)));
    }

//////////////////////////////////////////////////////////////////////////////
/*!
 * @fn IndexEsa#isaAt
 * @headerfile <seqan/index.h>
 * @brief Shortcut for <tt>value(indexIsa(..), ..)</tt>.
 *
 * @signature TValue isaAt(position, index);
 *
 * @param[in] index The @link Index @endlink object holding the fibre.
 * @param[in] position A position in the array on which the value should be accessed.
 *
 * @return TValue A reference or proxy to the value in the inverse suffix array.
 */

    template <typename TPos, typename TIndex>
    SEQAN_HOST_DEVICE inline typename Reference<typename Fibre<TIndex, FibreIsa>::Type>::Type isaAt(TPos i, TIndex &index) {
        return value(getFibre(index, FibreIsa()), i);
    }
    template <typename TPos, typename TIndex>
    SEQAN_HOST_DEVICE inline typename Reference<typename Fibre<TIndex const, FibreIsa>::Type>::Type isaAt(TPos i, TIndex const &index) {
        return value(getFibre(index, FibreIsa()), i);
    }
//////////////////////////////////////////////////////////////////////////////
/*!
 * @fn IndexEsa#lcpAt
 * @headerfile <seqan/index.h>
 * @brief Shortcut for <tt>value(indexLcp(..), ..)</tt>.
 *
 * @signature TValue lcpAt(position, index);
 *
 * @param[in] index    The @link Index @endlink object holding the fibre.
 * @param[in] position A position in the array on which the value should be accessed.
 *
 * @return TValue A reference or proxy to the value.
 */

    template <typename TPos, typename TIndex>
    SEQAN_HOST_DEVICE inline typename Reference<typename Fibre<TIndex, FibreLcp>::Type>::Type lcpAt(TPos i, TIndex &index) {
        return value(getFibre(index, FibreLcp()), i);
    }
    template <typename TPos, typename TIndex>
    SEQAN_HOST_DEVICE inline typename Reference<typename Fibre<TIndex const, FibreLcp>::Type>::Type lcpAt(TPos i, TIndex const &index) {
        return value(getFibre(index, FibreLcp()), i);
    }

//////////////////////////////////////////////////////////////////////////////
/*!
 * @fn IndexEsa#lcpeAt
 * @headerfile <seqan/index.h>
 * @brief Shortcut for <tt>value(indexLcpe(..), ..)</tt>.
 *
 * @signature TValue lcpeAt(position, index);
 *
 * @param[in] index The @link Index @endlink object holding the fibre.
 * @param[in] position A position in the array on which the value should be accessed.
 *
 * @return TValue A reference or proxy to the value.
 */

    template <typename TPos, typename TIndex>
    inline typename Reference<typename Fibre<TIndex, FibreLcpe>::Type>::Type lcpeAt(TPos i, TIndex &index) {
        return value(getFibre(index, FibreLcpe()), i);
    }
    template <typename TPos, typename TIndex>
    inline typename Reference<typename Fibre<TIndex const, FibreLcpe>::Type>::Type lcpeAt(TPos i, TIndex const &index) {
        return value(getFibre(index, FibreLcpe()), i);
    }

//////////////////////////////////////////////////////////////////////////////
/*!
 * @fn IndexEsa#childAt
 * @headerfile <seqan/index.h>
 * @brief Shortcut for <tt>value(indexChildtab(..), ..)</tt>.
 *
 * @signature TValue childAt(position, index);
 *
 * @param[in] position A position in the array on which the value should be accessed.
 * @param[in] index    The @link IndexEsa @endlink object holding the fibre.
 *
 * @return TValue A reference or proxy to the value.
 */

    template <typename TPos, typename TIndex>
    SEQAN_HOST_DEVICE inline typename Reference<typename Fibre<TIndex, FibreChildtab>::Type>::Type childAt(TPos i, TIndex &index) {
        return value(getFibre(index, FibreChildtab()), i);
    }
    template <typename TPos, typename TIndex>
    SEQAN_HOST_DEVICE inline typename Reference<typename Fibre<TIndex const, FibreChildtab>::Type>::Type childAt(TPos i, TIndex const &index) {
        return value(getFibre(index, FibreChildtab()), i);
    }

//////////////////////////////////////////////////////////////////////////////
/*!
 * @fn IndexEsa#bwtAt
 * @headerfile <seqan/index.h>
 * @brief Shortcut for <tt>value(indexBwt(..), ..)</tt>.
 *
 * @signature TValue bwtAt(position, index);
 *
 * @param[in] index    The @link IndexEsa @endlink object holding the fibre.
 * @param[in] position A position in the array on which the value should be accessed.
 *
 * @return TValue A reference or proxy to the value.
 */

    template <typename TPos, typename TIndex>
    inline typename Reference<typename Fibre<TIndex, FibreBwt>::Type>::Type bwtAt(TPos i, TIndex &index) {
        return value(getFibre(index, FibreBwt()), i);
    }
    template <typename TPos, typename TIndex>
    inline typename Reference<typename Fibre<TIndex const, FibreBwt>::Type>::Type bwtAt(TPos i, TIndex const &index) {
        return value(getFibre(index, FibreBwt()), i);
    }

//////////////////////////////////////////////////////////////////////////////

    template <typename TIndex, typename TPos, typename TSize>
    inline typename SAValue<TIndex>::Type toSuffixPosition(TIndex &, TPos i, TSize) {
        return i;
    }
    template <typename TIndex, typename TPos, typename TSize>
    inline typename SAValue<TIndex const>::Type toSuffixPosition(TIndex const &, TPos i, TSize) {
        return i;
    }

//////////////////////////////////////////////////////////////////////////////
// interface for infinity/invalid values

    template <typename TValue>
    SEQAN_HOST_DEVICE inline void _setSizeInval(TValue &v) {
        v = MaxValue<TValue>::VALUE;
    }

    template <typename TValue>
    SEQAN_HOST_DEVICE inline bool _isSizeInval(TValue const &v) {
//IOREV _notio_
        return v == MaxValue<TValue>::VALUE;
    }

//////////////////////////////////////////////////////////////////////////////
/*!
 * @fn Index#indexText
 * @headerfile <seqan/index.h>
 * @brief Shortcut for <tt>getFibre(.., EsaText)</tt>.
 *
 * @signature TFibre indexText(index);
 *
 * @param[in] index The @link Index @endlink object holding the fibre.
 *
 * @return TFibre A reference to the text of the index.
 *
 * @note The result of this function when used on an <tt>Index&lt;TText, FMIndex&lt;TOccSpec, CompressText&gt; &gt;
 *       is not defined.
 *
 * @section Examples
 *
 * The following code shows how the BWT of a text can be computed.
 *
 * @include demos/index/index_textAt_indexText_saAt_indexRequire.cpp
 *
 * The output is as follows:
 *
 * @include demos/index/index_textAt_indexText_saAt_indexRequire.cpp.stdout
 */

    template <typename TText, typename TSpec>
    inline typename Fibre<Index<TText, TSpec>, FibreText>::Type & indexText(Index<TText, TSpec> &index) { return getFibre(index, FibreText()); }
    template <typename TText, typename TSpec>
    inline typename Fibre<Index<TText, TSpec> const, FibreText>::Type & indexText(Index<TText, TSpec> const &index) { return getFibre(index, FibreText()); }

//////////////////////////////////////////////////////////////////////////////

    template <typename TText, typename TSpec>
    SEQAN_HOST_DEVICE inline typename StringSetLimits<TText const>::Type
    stringSetLimits(Index<TText, TSpec> &) {
        return Nothing();
    }

    template <typename TText, typename TSpec>
    SEQAN_HOST_DEVICE inline typename StringSetLimits<TText const>::Type
    stringSetLimits(Index<TText, TSpec> const &) {
        return Nothing();
    }

    template <typename TString, typename TSSetSpec, typename TSpec>
    SEQAN_HOST_DEVICE inline typename StringSetLimits< StringSet<TString, TSSetSpec> const >::Type &
    stringSetLimits(Index<StringSet<TString, TSSetSpec>, TSpec> &index) {
        return stringSetLimits(indexText(index));
    }

    template <typename TString, typename TSSetSpec, typename TSpec>
    SEQAN_HOST_DEVICE inline typename StringSetLimits< StringSet<TString, TSSetSpec> const >::Type &
    stringSetLimits(Index<StringSet<TString, TSSetSpec>, TSpec> const &index) {
        return stringSetLimits(indexText(index));
    }

    template <typename TString, typename TSSetSpec, typename TSpec>
    inline typename StringSetLimits< StringSet<TString, TSSetSpec> const >::Type &
    stringSetLimits(Index<StringSet<TString, TSSetSpec> const, TSpec> &index) {
        return stringSetLimits(indexText(index));
    }

    template <typename TString, typename TSSetSpec, typename TSpec>
    inline typename StringSetLimits< StringSet<TString, TSSetSpec> const >::Type &
    stringSetLimits(Index<StringSet<TString, TSSetSpec> const, TSpec> const &index) {
        return stringSetLimits(indexText(index));
    }

//////////////////////////////////////////////////////////////////////////////
/*!
 * @fn Index#indexRawText
 * @headerfile <seqan/index.h>
 * @brief Shortcut for <tt>$getFibre(.., FibreRawText)</tt>.
 *
 * @signature TFibre indexRawText(position, index);
 *
 * @param[in] position A position in the array on which the value should be accessed.
 * @param[in] index The @link Index @endlink object.
 *
 * @return TFibre A reference or proxy to the value.
 */

    template <typename TText, typename TSpec>
    inline SEQAN_HOST_DEVICE
    typename Fibre<Index<TText, TSpec>, FibreRawText>::Type &
    indexRawText(Index<TText, TSpec> &index)
    {
        return getFibre(index, FibreRawText());
    }

    template <typename TText, typename TSpec>
    inline SEQAN_HOST_DEVICE
    typename Fibre<Index<TText, TSpec> const, FibreRawText>::Type &
    indexRawText(Index<TText, TSpec> const &index)
    {
        return getFibre(index, FibreRawText());
    }

//////////////////////////////////////////////////////////////////////////////

/*!
 * @fn IndexEsa#indexSA
 * @headerfile <seqan/index.h>
 * @brief Shortcut for <tt>getFibre(.., EsaSA)</tt>.
 *
 * @signature TSA indexSA(index);
 *
 * @param[in] index The @link Index @endlink object holding the fibre.
 *
 * @return TSA A reference to the @link IndexEsaFibres#EsaSA @endlink fibre (suffix array).
 */

    template <typename TText, typename TSpec>
    SEQAN_HOST_DEVICE inline typename Fibre<Index<TText, TSpec>, FibreSA>::Type & indexSA(Index<TText, TSpec> &index) { return getFibre(index, FibreSA()); }
    template <typename TText, typename TSpec>
    SEQAN_HOST_DEVICE inline typename Fibre<Index<TText, TSpec> const, FibreSA>::Type & indexSA(Index<TText, TSpec> const &index) { return getFibre(index, FibreSA()); }

//////////////////////////////////////////////////////////////////////////////

// NOTE(esiragusa): indexRawSA() is internal and must stay invisible in the doc.

/*
 * @fn IndexEsa#indexRawSA
 * @headerfile <seqan/index.h>
 * @brief Shortcut for <tt>getFibre(.., EsaRawSA)</tt>.
 *
 * @signature TSA indexRawSA(index);
 *
 * @param[in] index The @link Index @endlink object holding the fibre.
 *
 * @return TSA A reference to the @link IndexEsaFibres#EsaRawSA @endlink fibre (suffix array).
 */

    template <typename TText, typename TSpec>
    inline typename Fibre<Index<TText, TSpec>, FibreRawSA>::Type indexRawSA(Index<TText, TSpec> &index) { return getFibre(index, FibreRawSA()); }
    template <typename TText, typename TSpec>
    inline typename Fibre<Index<TText, TSpec> const, FibreRawSA>::Type indexRawSA(Index<TText, TSpec> const &index) { return getFibre(index, FibreRawSA()); }

//////////////////////////////////////////////////////////////////////////////
/*!
 * @fn IndexEsa#indexIsa
 * @headerfile <seqan/index.h>
 * @brief Shortcut for <tt>getFibre(.., EsaIsa)</tt>.
 *
 * @signature TIsa indexIsa(index);
 *
 * @param[in] index The @link Index @endlink object holding the fibre.
 *
 * @return TIsa A reference to the inverse suffix array fibre.
 */

    template <typename TText, typename TSpec>
    SEQAN_HOST_DEVICE inline typename Fibre<Index<TText, TSpec>, FibreIsa>::Type & indexIsa(Index<TText, TSpec> &index) { return getFibre(index, FibreIsa()); }
    template <typename TText, typename TSpec>
    SEQAN_HOST_DEVICE inline typename Fibre<Index<TText, TSpec> const, FibreIsa>::Type & indexIsa(Index<TText, TSpec> const &index) { return getFibre(index, FibreIsa()); }

//////////////////////////////////////////////////////////////////////////////
/*!
 * @fn IndexEsa#indexLcp
 * @headerfile <seqan/index.h>
 * @brief Shortcut for <tt>getFibre(.., EsaLcp)</tt>.
 *
 * @signature TLcp indexLcp(index);
 *
 * @param[in] index The @link Index @endlink object holding the fibre.
 *
 * @return TLcp A reference to the @link IndexEsaFibres#EsaLcp @endlink fibre (lcp table).
 */

    template <typename TText, typename TSpec>
    SEQAN_HOST_DEVICE inline typename Fibre<Index<TText, TSpec>, FibreLcp>::Type & indexLcp(Index<TText, TSpec> &index) { return getFibre(index, FibreLcp()); }
    template <typename TText, typename TSpec>
    SEQAN_HOST_DEVICE inline typename Fibre<Index<TText, TSpec> const, FibreLcp>::Type & indexLcp(Index<TText, TSpec> const &index) { return getFibre(index, FibreLcp()); }

//////////////////////////////////////////////////////////////////////////////
/*!
 * @fn IndexEsa#indexLcpe
 * @headerfile <seqan/index.h>
 * @brief Shortcut for <tt>getFibre(.., EsaLcpe)</tt>.
 *
 * @signature TLcpe indexLcpe(index);
 *
 * @param[in] index The @link Index @endlink object holding the fibre.
 *
 * @return TLcpe A reference to the @link IndexEsaFibres#EsaLcpe @endlink fibre (enhanced lcp table).
 */

    template <typename TText, typename TSpec>
    inline typename Fibre<Index<TText, TSpec>, FibreLcpe>::Type & indexLcpe(Index<TText, TSpec> &index) { return getFibre(index, FibreLcpe()); }
    template <typename TText, typename TSpec>
    inline typename Fibre<Index<TText, TSpec> const, FibreLcpe>::Type & indexLcpe(Index<TText, TSpec> const &index) { return getFibre(index, FibreLcpe()); }

//////////////////////////////////////////////////////////////////////////////
/*!
 * @fn IndexEsa#indexBwt
 * @headerfile <seqan/index.h>
 * @brief Shortcut for <tt>getFibre(.., EsaBwt)</tt>.
 *
 * @signature TBwt indexBwt(index);
 *
 * @param[in] index The @link IndexEsa @endlink object holding the fibre.
 *
 * @return TBwt A reference to the @link IndexEsaFibres#EsaBwt @endlink fibre (Burrows-Wheeler table).
 */

    template <typename TText, typename TSpec>
    SEQAN_HOST_DEVICE inline typename Fibre<Index<TText, TSpec>, FibreBwt>::Type & indexBwt(Index<TText, TSpec> &index) { return getFibre(index, FibreBwt()); }
    template <typename TText, typename TSpec>
    SEQAN_HOST_DEVICE inline typename Fibre<Index<TText, TSpec> const, FibreBwt>::Type & indexBwt(Index<TText, TSpec> const &index) { return getFibre(index, FibreBwt()); }

//////////////////////////////////////////////////////////////////////////////
/*!
 * @fn IndexEsa#indexChildtab
 * @headerfile <seqan/index.h>
 * @brief Shortcut for <tt>getFibre(.., EsaChildtab)</tt>.
 *
 * @signature TChildTab indexChildtab(index);
 *
 * @param[in] index The @link Index @endlink object holding the fibre.
 *
 * @return TChildTab A reference to the @link IndexEsaFibres#EsaChildtab @endlink fibre (child table).
 */

    template <typename TText, typename TSpec>
    SEQAN_HOST_DEVICE inline typename Fibre<Index<TText, TSpec>, FibreChildtab>::Type & indexChildtab(Index<TText, TSpec> &index) { return getFibre(index, FibreChildtab()); }
    template <typename TText, typename TSpec>
    SEQAN_HOST_DEVICE inline typename Fibre<Index<TText, TSpec> const, FibreChildtab>::Type & indexChildtab(Index<TText, TSpec> const &index) { return getFibre(index, FibreChildtab()); }


// ----------------------------------------------------------------------------
// Function open
// ----------------------------------------------------------------------------
/*!
 * @fn Index#open
 * @headerfile <seqan/index.h>
 * @brief This functions opens an index from disk.
 *
 * @signature bool open(index, fileName[, mode]);
 *
 * @param[in,out]  index The index to be opened.
 * @param[in]      fileName C-style character string containing the file name.
 * @param[in] mode The combination of flags defining how the file should be opened.To open a file read-only, write-only or
 *                 to read and write use <tt>OPEN_RDONLY</tt>, <tt>OPEN_WRONLY</tt>, or <tt>OPEN_RDWR</tt>.To create or
 *                 overwrite a file add <tt>OPEN_CREATE</tt>.To append a file if existing add <tt>OPEN_APPEND</tt>.To
 *                 circumvent problems, files are always opened in binary mode.
 *                 Default: <tt>OPEN_RDWR | OPEN_CREATE | OPEN_APPEND</tt>.
 *
 * @return bool <tt>true</tt> on success and <tt>false</tt> otherwise.
 *
 * @section Examples
 *
 * The following code shows how the function @link Index#open @endlink is used with indices.
 *
 * @include demos/index/index_open_save.cpp
 *
 * The output is as follows:
 *
 * @include demos/index/index_open_save.cpp.stdout
 */

// ----------------------------------------------------------------------------
// Function save
// ----------------------------------------------------------------------------

/*!
 * @fn Index#save
 * @headerfile <seqan/index.h>
 * @brief This functions saves an index to disk.
 *
 * @signature bool save(index, fileName[, mode]);
 *
 * @param[in,out]  index The index to be opened.
 * @param[in]      fileName C-style character string containing the file name.
 * @param[in] mode The combination of flags defining how the file should be opened.To open a file read-only, write-only or
 *                 to read and write use <tt>OPEN_RDONLY</tt>, <tt>OPEN_WRONLY</tt>, or <tt>OPEN_RDWR</tt>.To create or
 *                 overwrite a file add <tt>OPEN_CREATE</tt>.To append a file if existing add <tt>OPEN_APPEND</tt>.To
 *                 circumvent problems, files are always opened in binary mode.
 *                 Default: <tt>OPEN_RDWR | OPEN_CREATE | OPEN_APPEND</tt>.
 *
 * @return bool <tt>true</tt> on success and <tt>false</tt> otherwise.
 *
 * @section Examples
 *
 * The following code shows how the function @link Index#save @endlink is used with indices.
 *
 * @include demos/index/index_open_save.cpp
 *
 * The output is as follows:
 *
 * @include demos/index/index_open_save.cpp.stdout
 */
}

#endif

