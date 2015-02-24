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

#ifndef SEQAN_HEADER_INDEX_ESA_BASE_H
#define SEQAN_HEADER_INDEX_ESA_BASE_H

namespace SEQAN_NAMESPACE_MAIN
{

    // dfs order
    struct Preorder_;
    struct Postorder_;

    template <typename TDfsOrder = Postorder_, typename THideEmptyEdges = True>
    struct VSTreeIteratorTraits {
        typedef TDfsOrder DfsOrder;
        typedef THideEmptyEdges HideEmptyEdges;
    };

/*!
 * @defgroup DfsOrder DFS Order
 * @brief Pre/postorder selection for depth-first search.
 *
 * These tags are given to @link InputIteratorConcept#goNext @endlink and trigger post-order or pre-
 * order traversal of a suffix tree. In case of <tt>PreorderEmptyEdges</tt> and
 * <tt>PostorderEmptyEdges</tt>, the empty edges are also traversed.
 *
 * @tag DfsOrder#Preorder
 * @brief Visit the node before its children.
 *
 * @tag DfsOrder#PostorderEmptyEdges
 * @brief Visit the node after its children, visit empty edges.
 *
 * @tag DfsOrder#PreorderEmptyEdges
 * @brief Visit the node before its children, visit empty edges.
 *
 * @tag DfsOrder#Postorder
 * @brief Visit the node after its children.
 */
    // predefined iterator traits
    struct Preorder:            VSTreeIteratorTraits<Preorder_,  True> {};
    struct Postorder:            VSTreeIteratorTraits<Postorder_, True> {};
    struct PreorderEmptyEdges:    VSTreeIteratorTraits<Preorder_,  False> {};    // also iterate over
    struct PostorderEmptyEdges:    VSTreeIteratorTraits<Postorder_, False> {};    // empty edges (with $-label)

    // traits for TopDown iterators (w/o ParentLinks) for which postorder/preorder is ignored
    struct HideEmptyEdges:        VSTreeIteratorTraits<Postorder_, True> {};
    struct EmptyEdges:            VSTreeIteratorTraits<Postorder_, False> {};    // empty edges (with $-label)

    // MultiMems are more specialized MaxRepeats
    template <typename TSpec = void>
    struct MaxRepeats_;    // base class
    struct MultiMems_;    // subclass tag



    // virtual suffix tree iterators
    template <typename TSpec = void>
    struct VSTree;

/*!
 * @defgroup TopDown Top-Down Iteration
 * @brief Tag that specifies a @link VSTreeIterator @endlink to traverse the virtual string tree from the root towards
 *        the leafs.
 *
 * @section Examples
 *
 * The following example shows how a the @link TopDown @endlink tag is used.
 *
 * @include demos/index/index_begin_atEnd_representative.cpp
 *
 * @code{.output}
 * A
 * AA
 * ATAA
 * TA
 * TAA
 * TATAA
 * --------------------------------
 * AA
 * ATAA
 * A
 * TAA
 * TATAA
 * TA
 * @endcode
 *
 * @tag TopDown#ParentLinks
 * @brief A top down iterator with the possibility to go back up again.
 *
 * @tag TopDown#Preorder
 * @brief Pre-order traversal of the virtual string tree.
 *
 * @tag TopDown#Postorder
 * @brief Post-order traversal of the virtual string tree.
 */

        // top down traversal iterators
        template <typename TSpec = Preorder>
        struct TopDown {};                    // starts in the suffix tree root and can go down and go right

            // allows an top-down iterator to go up
            template < typename TSpec = Preorder >
            struct ParentLinks {};            // .. can also go up

/*!
 * @defgroup BottomUp Bottom-Up Iteration
 * @brief Tag that specifies a @link VSTreeIterator @endlink to traverse the
 *        virtual string tree from the root towards the leafs.
 *
 * @section Examples
 *
 * The following example shows how the @link BottomUp @endlink tag is used.
 *
 * @include demos/index/index_begin_atEnd_representative_bottomUp.cpp
 *
 * @code{.txt}
 * AA
 * ATAA
 * A
 * TAA
 * TATAA
 * TA
 * @endcode
 *
 * @tag BottomUp#Postorder
 * @brief Post-order traversal of the virtual string tree.
 */
        // bottom up traversal iterators
        template <typename TSpec = Postorder>
        struct BottomUp {};                    // starts in the first node of a depth-first-search and can go next

            struct    SuperMaxRepeats;                    // maximal repeat and not part of a longer repeat
            struct    SuperMaxRepeatsFast;
            struct    Mums;                                // Maximal Unique Match (unique in every sequence)

            typedef MaxRepeats_<void>        MaxRepeats;    // maximal repeat
            struct    MaxRepeatOccurrences;
            typedef MaxRepeats_<MultiMems_> MultiMems;    // Multiple Maximal Exact Match
            struct    MultiMemOccurences;                    // i.e. maximal match over different sequences


/*!
 * @mfn Index#GetVSTreeIteratorTraits
 *
 * @headerfile <seqan/index.h>
 *
 * @brief Default behaviour of @link InputIteratorConcept#goNext @endlink when no second parameter is given.
 *
 * @signature GetVSTreeIteratorTraits<TIterator>::Type
 *
 * @tparam TIterator A @link VSTreeIterator @endlink.
 *
 * @return TReturn @link DfsOrder#Postorder @endlink by default and @link DfsOrder#Preorder @endlink
 *                 if <tt>TIterator</tt> is <tt>VSTree&lt;TopDown&lt;ParentLinks&lt;&gt; &gt; &gt;</tt>
 *                 or <tt>VSTree&lt;TopDown&lt;ParentLinks&lt;Preorder&gt; &gt; &gt;</tt>.
 */

    template <typename TIterator>
    struct GetVSTreeIteratorTraits:
        DeepestSpec<TIterator> {};

//////////////////////////////////////////////////////////////////////////////

    template <typename TSize>
    struct VertexEsa {
        Pair<TSize> range;            // current SA interval of hits (unique node identifier)
        TSize        parentRight;    // right boundary of parent node's range (allows to go right)

        SEQAN_HOST_DEVICE
        VertexEsa() : range(0, 0), parentRight(0) {}

        SEQAN_HOST_DEVICE
        VertexEsa(MinimalCtor):
            range(0,0),
            parentRight(0) {}

        SEQAN_HOST_DEVICE
        VertexEsa(TSize otherRangeLeft, TSize otherRangeRight, TSize otherParentRight):
            range(Pair<TSize>(otherRangeLeft, otherRangeRight)),
            parentRight(otherParentRight) {}

        SEQAN_HOST_DEVICE
        VertexEsa(Pair<TSize> const &otherRange, TSize otherParentRight):
            range(otherRange),
            parentRight(otherParentRight) {}

        SEQAN_HOST_DEVICE
        VertexEsa(VertexEsa const &other):
            range(other.range),
            parentRight(other.parentRight) {}
    };

    template <typename TSize>
    SEQAN_HOST_DEVICE inline bool operator==(VertexEsa<TSize> const &a, VertexEsa<TSize> const &b)
    {
        return a.range == b.range;
    }

    template <typename TSize>
    SEQAN_HOST_DEVICE inline bool operator!=(VertexEsa<TSize> const &a, VertexEsa<TSize> const &b)
    {
        return a.range != b.range;
    }

//////////////////////////////////////////////////////////////////////////////
/*!
 * @mfn StringTreeConcept#VertexDescriptor
 * @headerfile <seqan/index.h>
 * @brief Returns the type of an object that represents a string tree node.
 *
 * @signature VertexDescriptor<TIndex>::Type;
 *
 * @tparam TIndex The index type.
 *
 * @return Type The resulting vertex descriptor type.
 */

    template < typename TText, typename TSpec >
    struct VertexDescriptor< Index<TText, IndexEsa<TSpec> > > {
        typedef typename Size< Index<TText, IndexEsa<TSpec> > >::Type TSize;
        typedef VertexEsa<TSize> Type;
    };


//////////////////////////////////////////////////////////////////////////////
// needful forward declarations

    struct ArrayGaps_;
    typedef Tag<ArrayGaps_> ArrayGaps;

    template <typename TSource, typename TSpec>
    class Align;


//////////////////////////////////////////////////////////////////////////////
// ESA fibres

/*!
 * @defgroup IndexEsaFibres Index Esa Fibres
 * @brief Tag to select a specific fibre (e.g. table, object, ...) of an @link
 *        IndexEsa @endlink index.
 *
 * These tags can be used to get @link Fibre Fibres @endlink of an Enhanced
 * Suffix Array based @link IndexEsa @endlink.
 *
 * @see Fibre
 * @see Index#getFibre
 * @see IndexEsa
 *
 * @tag IndexEsaFibres#EsaSA
 * @headerfile <seqan/index.h>
 * @brief The suffix array.
 *
 * The suffix array contains the indices of all suffices of <tt>EsaRawText</tt>
 * in lexicographical order.
 *
 * @link Fibre @endlink returns a @link String @endlink over the alphabet of the
 * @link SAValue @endlink of <tt>TIndex</tt>.
 *
 * @tag IndexEsaFibres#EsaIsa
 * @headerfile <seqan/index.h>
 * @brief The inverse suffix array.
 *
 * The inverse suffix array stores the lexicographic rank of each suffix of <tt>EsaRawText</tt>.
 *
 * @link Fibre @endlink returns a @link String @endlink over the alphabet of a
 * size type.
 *
 * @tag IndexEsaFibres#EsaChildtab
 * @headerfile <seqan/index.h>
 * @brief The child table.
 *
 * The child table contains structural information of the suffix tree (see
 * Abhouelda et al.).
 *
 * @link Fibre @endlink returns a @link String @endlink over the alphabet of a
 * size type.
 *
 * @tag IndexEsaFibres#EsaRawText
 * @headerfile <seqan/index.h>
 * @brief The raw text the index is really based on.
 *
 * <tt>EsaText</tt> and <tt>EsaRawText</tt> fibres are equal by default. They
 * differ if the index text is a set of strings. Then, raw text is the
 * concatenation of all strings in this set.
 *
 * @tag IndexEsaFibres#EsaText
 * @headerfile <seqan/index.h>
 * @brief The original text the index should be based on.
 *
 * @tag IndexEsaFibres#EsaBwt
 * @headerfile <seqan/index.h>
 * @brief The Burrows-Wheeler table.
 *
 * The Burrows-Wheeler table contains the Burrows-Wheeler transformation of
 * <tt>EsaRawText</tt>. The entries are the characters left of the corresponding
 * suffix in the suffix array <tt>EsaSA</tt>.
 *
 * @link Fibre @endlink returns the same type for <tt>EsaRawText</tt> and for
 * <tt>EsaBwt</tt>.
 *
 * @tag IndexEsaFibres#EsaLcp
 * @headerfile <seqan/index.h>
 * @brief The lcp table.
 *
 * The lcp table contains the lcp-value of two adjacent suffices in the suffix
 * array <tt>EsaSA</tt>.
 *
 * @link Fibre @endlink returns a @link String @endlink over the alphabet of a
 * size type.
 *
 * @tag IndexEsaFibres#EsaLcpe
 * @headerfile <seqan/index.h>
 * @brief The lcpe table.
 */

    typedef FibreText       EsaText;
    typedef FibreRawText    EsaRawText;
    typedef FibreSA         EsaSA;
    typedef FibreIsa        EsaIsa;
    typedef FibreRawSA      EsaRawSA;
    typedef FibreSae        EsaSae;
    typedef FibreLcp        EsaLcp;
    typedef FibreLcpe       EsaLcpe;
    typedef FibreChildtab   EsaChildtab;
    typedef FibreBwt        EsaBwt;

//////////////////////////////////////////////////////////////////////////////
// ESA index

/*!
 * @class IndexEsa
 * @extends Index
 * @implements StringTreeConcept
 * @headerfile <seqan/index.h>
 * @brief An index based on an enhanced suffix array.
 *
 * @signature template <typename TText, typename TSpec>
 *            class Index<TText, IndexEsa<TSpec> >;
 *
 * @tparam TText The @link TextConcept text type @endlink.
 * @tparam TSpec The specialization, defaults to <tt>void</tt>.
 *
 * The fibres (see @link Index @endlink and @link Fibre @endlink) of this index are a suffix array (see @link
 * IndexEsaFibres#EsaSA @endlink), a lcp table (see @link IndexEsaFibres#EsaLcp @endlink), etc.
 *
 * This index can be accessed as a Suffix Tree using the @link VSTreeIterator @endlink classes.
 *
 * @see IndexEsaFibres
 */

/*
    already defined in index_base.h

    template <typename TSpec = void>
    struct IndexEsa;s
*/

    template < typename TText, typename TSpec >
    class Index<TText, IndexEsa<TSpec> > {
    public:
        typename Member<Index, EsaText>::Type       text;
        typename Fibre<Index, EsaSA>::Type          sa;            // suffix array
        typename Fibre<Index, EsaIsa>::Type         isa;        // inverse suffix array
        typename Fibre<Index, EsaLcp>::Type         lcp;        // longest-common-prefix table
        typename Fibre<Index, EsaLcpe>::Type        lcpe;        // extended lcp table
        typename Fibre<Index, EsaChildtab>::Type    childtab;    // child table (tree topology)
        typename Fibre<Index, EsaBwt>::Type         bwt;        // burrows-wheeler table
        typename Cargo<Index>::Type                 cargo;        // user-defined cargo

        Index() {}

        Index(Index &other):
            text(other.text),
            sa(other.sa),
            isa(other.isa),
            lcp(other.lcp),
            lcpe(other.lcpe),
            childtab(other.childtab),
            bwt(other.bwt),
            cargo(other.cargo) {}

        Index(Index const &other):
            text(other.text),
            sa(other.sa),
            isa(other.isa),
            lcp(other.lcp),
            lcpe(other.lcpe),
            childtab(other.childtab),
            bwt(other.bwt),
            cargo(other.cargo) {}

        template <typename TText_>
        Index(TText_ &_text):
            text(_text) {}

        template <typename TText_>
        Index(TText_ const &_text):
            text(_text) {}
    };

//////////////////////////////////////////////////////////////////////////////

    template < typename TText, typename TSpec >
    SEQAN_HOST_DEVICE inline void _indexRequireTopDownIteration(Index<TText, IndexEsa<TSpec> > &index)
    {
        indexRequire(index, EsaSA());
        indexRequire(index, EsaLcp());
        indexRequire(index, EsaChildtab());
    }

    template < typename TText, typename TSpec >
    void _indexRequireBottomUpIteration(Index<TText, IndexEsa<TSpec> > &index)
    {
        indexRequire(index, EsaSA());
        indexRequire(index, EsaLcp());
    }

//////////////////////////////////////////////////////////////////////////////
/*!
 * @fn Index#clear
 * @brief Resets all fibres of an index.
 *
 * @signature void clear(index);
 *
 * @param[in,out] index The index to be cleared.
 */

    template <typename TText, typename TSpec>
    inline void clear(Index<TText, IndexEsa<TSpec> > &index) {
        clear(getFibre(index, EsaSA()));
        clear(getFibre(index, EsaIsa()));
        clear(getFibre(index, EsaLcp()));
        clear(getFibre(index, EsaLcpe()));
        clear(getFibre(index, EsaChildtab()));
        clear(getFibre(index, EsaBwt()));
    }

// ----------------------------------------------------------------------------
// Function open
// ----------------------------------------------------------------------------

    template < typename TObject, typename TSpec >
    inline bool open(
        Index< TObject, IndexEsa<TSpec> > &index,
        const char *fileName,
        int openMode)
    {
        String<char> name;

        name = fileName;    append(name, ".txt");
        if ((!open(getFibre(index, EsaText()), toCString(name), openMode)) &&
            (!open(getFibre(index, EsaText()), fileName, openMode))) return false;

        name = fileName;    append(name, ".sa");
        if (!open(getFibre(index, EsaSA()), toCString(name), openMode)) return false;

        name = fileName;    append(name, ".isa");
        if (!open(getFibre(index, EsaIsa()), toCString(name), openMode)) return false;

        name = fileName;    append(name, ".lcp");
        if (!open(getFibre(index, EsaLcp()), toCString(name), openMode)) return false;

        name = fileName;    append(name, ".child");
        if (!open(getFibre(index, EsaChildtab()), toCString(name), openMode)) return false;

        name = fileName;    append(name, ".bwt");
        if (!open(getFibre(index, EsaBwt()), toCString(name), openMode)) return false;

        return true;
    }
    template < typename TObject, typename TSpec >
    inline bool open(
        Index< TObject, IndexEsa<TSpec> > &index,
        const char *fileName)
    {
        return open(index, fileName, DefaultOpenMode<Index< TObject, IndexEsa<TSpec> > >::VALUE);
    }

// ----------------------------------------------------------------------------
// Function save
// ----------------------------------------------------------------------------

    template < typename TObject, typename TSpec >
    inline bool save(
        Index< TObject, IndexEsa<TSpec> > const &index,
        const char *fileName,
        int openMode)
    {
        String<char> name;

        name = fileName;    append(name, ".txt");
        if ((!save(getFibre(index, EsaText()), toCString(name), openMode)) &&
            (!save(getFibre(index, EsaText()), fileName, openMode))) return false;

        name = fileName;    append(name, ".sa");
        if (!save(getFibre(index, EsaSA()), toCString(name), openMode)) return false;

        name = fileName;    append(name, ".isa");
        if (!save(getFibre(index, EsaIsa()), toCString(name), openMode)) return false;

        name = fileName;    append(name, ".lcp");
        if (!save(getFibre(index, EsaLcp()), toCString(name), openMode)) return false;

        name = fileName;    append(name, ".child");
        if (!save(getFibre(index, EsaChildtab()), toCString(name), openMode)) return false;

        name = fileName;    append(name, ".bwt");
        if (!save(getFibre(index, EsaBwt()), toCString(name), openMode)) return false;

        return true;
    }
    template < typename TObject, typename TSpec >
    inline bool save(
        Index< TObject, IndexEsa<TSpec> > const &index,
        const char *fileName)
    {
        return save(index, fileName, DefaultOpenMode<Index< TObject, IndexEsa<TSpec> > >::VALUE);
    }

}

#endif
