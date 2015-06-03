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

#ifndef SEQAN_HEADER_INDEX_QGRAM_H
#define SEQAN_HEADER_INDEX_QGRAM_H

namespace SEQAN_NAMESPACE_MAIN
{


//////////////////////////////////////////////////////////////////////////////
// q-gram index fibres

/*!
 * @defgroup QGramIndexFibres Index QGram Fibres
 * @brief Tag to select a specific fibre (e.g. table, object, ...) of a @link
 *        IndexQGram @endlink.
 *
 * @see Fibre
 * @see Index#getFibre
 * @see IndexQGram
 *
 * @tag QGramIndexFibres#QGramDir
 * @brief The directory/hash table.
 *
 * The directory contains for every possible q-gram hash value the start index
 * of the q-gram bucket. A q-gram bucket is a contiguous interval in the suffix
 * array (<tt>QGramSA</tt>, see above). Each suffix in this interval begins with
 * the same q-gram. The end index is the start index of the next bucket.
 *
 * @link Fibre @endlink returns a @link String @endlink over the alphabet of a
 * size type.
 *
 * @tag QGramIndexFibres#QGramBucketMap
 * @brief Maps q-gram hashes to buckets. This fibre is used by the @link
 *        OpenAddressingQGramIndex @endlink index and stores all parameters of the open
 *        addressing hash function and hash value occupancy in the QGramDir
 *        fibre. In contrast to @link OpenAddressingQGramIndex @endlink, @link IndexQGram
 *        @endlink uses a trivial 1-to-1 mapping from q-gram hash values to
 *        buckets. For that index the fibre is of type @link Nothing @endlink.
 *
 * @tag QGramIndexFibres#QGramCountsDir
 * @brief The counts directory.
 *
 * The counts directory contains for every possible q-gram hash value the start
 * index of the q-gram count bucket. A q-gram count bucket is a contiguous
 * interval in the counts array (<tt>QGramCounts</tt>, see above). The end index
 * is the start index of the next bucket.
 *
 * @link Fibre @endlink returns a @link String @endlink over the alphabet of a
 * size type.
 *
 * @tag QGramIndexFibres#QGramText
 * @brief The original text the index should be based on.
 *
 * @tag QGramIndexFibres#QGramShape
 * @brief The shape the index is based on.
 *
 * The q-gram index needs an underlying @link Shape @endlink. This shape can be
 * gapped or ungapped. The number of '1's (relevant positions) in the shape
 * determines <tt>q</tt> and the size of the directory table.
 *
 * Dynamic shapes (@link SimpleShape @endlink, @link GenericShape @endlink, ...)
 * must be initialized before the index can be used.
 *
 * @tag QGramIndexFibres#QGramSADir
 * @brief The union of suffix array and directory.
 *
 * In most applications a q-gram index consisting of both of these table is
 * required. To efficiently create them at once use this tag for @link
 * Index#indexRequire @endlink or @link Index#indexCreate @endlink.
 *
 * @tag QGramIndexFibres#QGram_RawText
 * @brief The concatenation of all text sequences.
 *
 * <tt>QGramText</tt> and <tt>QGram_RawText</tt> fibres are equal by default.
 * They differ if the index text is a set of strings. Then, raw text is the
 * concatenation of all strings in this set.
 *
 * @tag QGramIndexFibres#QGramCounts
 * @brief The counts array.
 *
 * Contains the numbers of occurrences per sequence of each q-gram, s.t. the
 * numbers of the same q-gram are stored in a contiguous block (q-gram count
 * bucket). A bucket contains entries (seqNo,count) of sequences with at least
 * one q-gram occurrence. q-grams exceeding the end of the text are ignored. The
 * beginning of each count bucket can be determined by the q-gram counts
 * directory (<tt>QGramCountsDir</tt>, see below).
 *
 * @link Fibre @endlink returns a @link String @endlink over the alphabet of the
 * @link SAValue @endlink of <tt>TIndex</tt>.
 *
 * @tag QGramIndexFibres#QGramSA
 * @brief The suffix array.
 *
 * Contains all occurrences of q-grams, s.t. the occurrences of a single q-gram
 * are stored in a contiguous block (q-gram bucket). q-grams exceeding the end
 * of the text are ignored. The beginning of each bucket can be determined by
 * the q-gram directory (<tt>QGramDir</tt>, see below).
 *
 * It corresponds to a suffix array which is sorted by the first q-gram.
 *
 * @link Fibre @endlink returns a @link String @endlink over the alphabet of the
 * @link SAValue @endlink of <tt>TIndex</tt>.
 */

    struct FibreDir_;            // directory/hash table, contains start indices of buckets
    struct FibreSADir_;        // identifies algorithm to construct both SA and directory at once
    struct Fibre_Shape_;        // underlying shape
    struct FibreCounts_;        // counts each q-gram
    struct FibreCountsDir_;    // directory for counts buckets
    struct FibreBucketMap_;    // stores a q-gram hash value for each directory entry (-1 if empty)

    typedef Tag<FibreDir_> const        FibreDir;
    typedef Tag<FibreSADir_> const        FibreSADir;
    typedef Tag<Fibre_Shape_> const        FibreShape;
    typedef Tag<FibreCounts_> const    FibreCounts;
    typedef Tag<FibreCountsDir_> const    FibreCountsDir;
    typedef Tag<FibreBucketMap_> const    FibreBucketMap;

//////////////////////////////////////////////////////////////////////////////

    typedef FibreText        QGramText;
    typedef FibreRawText    QGram_RawText;
    typedef FibreSA         QGramSA;
    typedef FibreRawSA        QGramRawSA;
    typedef FibreDir        QGramDir;
    typedef FibreSADir        QGramSADir;
    typedef FibreShape        QGramShape;
    typedef FibreCounts     QGramCounts;
    typedef FibreCountsDir    QGramCountsDir;
    typedef FibreBucketMap    QGramBucketMap;

//////////////////////////////////////////////////////////////////////////////
// q-gram index

/*!
 * @defgroup OpenAdressingTags Open Adressing Tags
 * @brief Tags for selecting open adressing.
 *
 * @tag OpenAdressingTags#OpenAdressing
 * @brief Tag for enabling open adressing.
 */

/*!
 * @class IndexQGram
 * @extends Index
 * @headerfile <seqan/index.h>
 *
 * @brief An index based on an array of sorted q-grams.
 *
 * @signature template <typename TIndex, typename TShapeSpec, typename TSpec>
 *            class Index<TText, IndexQGram<TShapeSpec[, TSpec]> >;
 *
 * @tparam TSpec The specializing type. Types: @link OpenAdressingTags#OpenAdressing @endlink, Default: void
 * @tparam TText The text type. Types: @link String @endlink
 * @tparam TShapeSpec The @link Shape @endlink specialization type.
 *
 * The fibres (see @link Index @endlink and @link Fibre @endlink) of this index are a suffix array sorted by the first
 * <i>q</i> characters (see @link QGramIndexFibres#QGramSA @endlink) and a <i>q</i>-gram directory (see @link
 * QGramIndexFibres#QGramDir @endlink).
 *
 * The size of the q-gram directory is <i>&Sigma;<sup>q</sup></i>.  On a 32 bit system the <i>q</i>-gram length is
 * limited to 3 for <tt>char</tt> alphabets or 13-14 for @link Dna @endlink alphabets.
 *
 * Consider to use the @link OpenAddressingQGramIndex @endlink for longer <i>q</i>-grams if you don't need
 * <i>q</i>-grams to be sorted.
 *
 * @see QGramIndexFibres
 */

    template < typename TShapeSpec, typename TSpec = Default >
    struct IndexQGram {};

    // use the index value type as shape value type
    template < typename TText, typename TShapeSpec, typename TSpec >
    struct Fibre< Index<TText, IndexQGram<TShapeSpec, TSpec> >, FibreShape>
    {
        typedef Index< TText, IndexQGram<TShapeSpec, TSpec> >    TIndex;
        typedef Shape< typename Value<TIndex>::Type, TShapeSpec >    Type;
    };

    // allow different value types for the shape
    template < typename TText, typename TShapeValue, typename TShapeSpec, typename TSpec >
    struct Fibre< Index<TText, IndexQGram<Shape<TShapeValue, TShapeSpec>, TSpec> >, FibreShape>
    {
        typedef Shape<TShapeValue, TShapeSpec>    Type;
    };

    template < typename TText, typename TShapeSpec, typename TSpec >
    struct Fibre< Index<TText, IndexQGram<TShapeSpec, TSpec> >, FibreBucketMap>
    {
        typedef Nothing Type;
    };

#ifdef PLATFORM_WINDOWS_VS
#pragma warning( push )
// Disable warning C4521 locally (multiple copy constructors).
#pragma warning( disable: 4521 )
// Disable warning C4522 locally (multiple assignment operators).
#pragma warning( disable: 4522 )
#endif  // PLATFORM_WINDOWS_VS

    template < typename TText_, typename TShapeSpec, typename TSpec >
    class Index<TText_, IndexQGram<TShapeSpec, TSpec> > {
    public:
        typedef typename Member<Index, QGramText>::Type         TTextMember;
        typedef typename Fibre<Index, QGramText>::Type            TText;
        typedef typename Fibre<Index, QGramSA>::Type            TSA;
        typedef typename Fibre<Index, QGramDir>::Type            TDir;
        typedef typename Fibre<Index, QGramCounts>::Type        TCounts;
        typedef typename Fibre<Index, QGramCountsDir>::Type     TCountsDir;
        typedef typename Fibre<Index, QGramShape>::Type         TShape;
        typedef typename Fibre<Index, QGramBucketMap>::Type     TBucketMap;
        typedef typename Cargo<Index>::Type                        TCargo;
        typedef typename Size<Index>::Type                        TSize;

        TTextMember     text;        // underlying text
        TSA                sa;            // suffix array sorted by the first q chars
        TDir            dir;        // bucket directory
        TCounts            counts;        // counts each q-gram per sequence
        TCountsDir        countsDir;    // directory for count buckets
        TShape            shape;        // underlying shape
        TCargo            cargo;        // user-defined cargo
        TBucketMap        bucketMap;    // bucketMap table (used by open-addressing index)
        TSize            stepSize;    // store every <stepSize>'th q-gram in the index

        Index():
            stepSize(1) {}

        Index(Index &other):
            text(other.text),
            sa(other.sa),
            dir(other.dir),
            counts(other.counts),
            countsDir(other.countsDir),
            shape(other.shape),
            cargo(other.cargo),
            stepSize(1) {}

        Index(Index const &other):
            text(other.text),
            sa(other.sa),
            dir(other.dir),
            counts(other.counts),
            countsDir(other.countsDir),
            shape(other.shape),
            cargo(other.cargo),
            stepSize(1) {}

        template <typename TText__>
        Index(TText__ &_text):
            text(_text),
            stepSize(1) {}

        template <typename TText__>
        Index(TText__ const &_text):
            text(_text),
            stepSize(1) {}

        template <typename TText__, typename TShape_>
        Index(TText__ &_text, TShape_ const &_shape):
            text(_text),
            shape(_shape),
            stepSize(1) {}

        template <typename TText__, typename TShape_>
        Index(TText__ const &_text, TShape_ const &_shape):
            text(_text),
            shape(_shape),
            stepSize(1) {}
    };

#ifdef PLATFORM_WINDOWS_VS
// Reset warning state to previous values for C4521, C4522.
#pragma warning( pop )
#endif  // PLATFORM_WINDOWS_VS

    template < typename TText, typename TShapeSpec, typename TSpec >
    struct Value< Index<TText, IndexQGram<TShapeSpec, TSpec> > > {
        typedef typename Value< typename Fibre< Index<TText, IndexQGram<TShapeSpec, TSpec> >, QGram_RawText >::Type >::Type Type;
    };

    template < typename TText, typename TShapeSpec, typename TSpec >
    struct Size< Index<TText, IndexQGram<TShapeSpec, TSpec> > > {
        typedef typename Size< typename Fibre< Index<TText, IndexQGram<TShapeSpec, TSpec> >, QGram_RawText >::Type >::Type Type;
    };


//////////////////////////////////////////////////////////////////////////////
// default fibre creators

    template < typename TText, typename TShapeSpec, typename TSpec >
    struct DefaultIndexCreator<Index<TText, IndexQGram<TShapeSpec, TSpec> >, FibreSA> {
        typedef Default Type;
    };

//////////////////////////////////////////////////////////////////////////////
// counts array type

    template < typename TText, typename TSpec >
    struct Fibre< Index<TText, TSpec>, FibreCounts> {
        typedef String<
                Pair<
                    typename Size< TText >::Type,
                    typename Size< Index<TText, TSpec> >::Type
                >,
                typename DefaultIndexStringSpec< Index<TText, TSpec> >::Type
        > Type;
    };


//////////////////////////////////////////////////////////////////////////////

    template <typename TText, typename TSpec>
    inline typename Fibre<Index<TText, TSpec>, FibreDir>::Type &
    getFibre(Index<TText, TSpec> &index, FibreDir) {
        return index.dir;
    }
    template <typename TText, typename TSpec>
    inline typename Fibre<Index<TText, TSpec> const, FibreDir>::Type &
    getFibre(Index<TText, TSpec> const &index, FibreDir) {
        return index.dir;
    }

    template <typename TText, typename TSpec>
    inline typename Fibre<Index<TText, TSpec>, FibreCounts>::Type &
    getFibre(Index<TText, TSpec> &index, FibreCounts) {
        return index.counts;
    }
    template <typename TText, typename TSpec>
    inline typename Fibre<Index<TText, TSpec> const, FibreCounts>::Type &
    getFibre(Index<TText, TSpec> const &index, FibreCounts) {
        return index.counts;
    }

    template <typename TText, typename TSpec>
    inline typename Fibre<Index<TText, TSpec>, FibreCountsDir>::Type &
    getFibre(Index<TText, TSpec> &index, FibreCountsDir) {
        return index.countsDir;
    }
    template <typename TText, typename TSpec>
    inline typename Fibre<Index<TText, TSpec> const, FibreCountsDir>::Type &
    getFibre(Index<TText, TSpec> const &index, FibreCountsDir) {
        return index.countsDir;
    }

    template <typename TText, typename TSpec>
    inline typename Fibre<Index<TText, TSpec>, FibreBucketMap>::Type &
    getFibre(Index<TText, TSpec> &index, FibreBucketMap) {
        return index.bucketMap;
    }
    template <typename TText, typename TSpec>
    inline typename Fibre<Index<TText, TSpec> const, FibreBucketMap>::Type &
    getFibre(Index<TText, TSpec> const &index, FibreBucketMap) {
        return index.bucketMap;
    }

    template <typename TText, typename TSpec>
    inline typename Fibre<Index<TText, TSpec>, FibreShape>::Type &
    getFibre(Index<TText, TSpec> &index, FibreShape) {
        return index.shape;
    }
    template <typename TText, typename TSpec>
    inline typename Fibre<Index<TText, TSpec> const, FibreShape>::Type &
    getFibre(Index<TText, TSpec> const &index, FibreShape) {
        return index.shape;
    }

//////////////////////////////////////////////////////////////////////////////
/*!
 * @fn IndexQGram#indexSA
 * @headerfile <seqan/index.h>
 * @brief Shortcut for <tt>getFibre(.., QGramSA)</tt>.
 *
 * @signature TSa indexSA(index);
 *
 * @param[in] index The @link IndexQGram @endlink object holding the fibre.
 *
 * @return TSa A reference to the @link QGramIndexFibres#QGramSA @endlink fibre (q-gram positions).
 */

/*!
 * @fn IndexQGram#saAt
 * @headerfile <seqan/index.h>
 * @note Advanced functionality, not commonly used.
 * @brief Shortcut for <tt>value(indexSA(..), ..)</tt>.
 *
 * @signature TValue saAt(position, index);
 *
 * @param[in] index The @link IndexQGram @endlink object holding the fibre.
 * @param[in] position A position in the array on which the value should be accessed.
 *
 * @return TValue A reference or proxy to the value in the @link QGramIndexFibres#QGramSA @endlink fibre.
 *                To be more precise, a reference to a position containing a value of type
 *                @link SAValue @endlink is returned (or a proxy).
 */

//////////////////////////////////////////////////////////////////////////////
/*!
 * @fn IndexQGram#indexDir
 * @headerfile <seqan/index.h>
 * @brief Shortcut for <tt>getFibre(.., QGramDir())</tt>.
 * @signature TFibre indexDir(index);
 *
 * @param[in] index The @link IndexQGram @endlink object holding the fibre.
 *
 * @return TFibre A reference to the @link QGramIndexFibres#QGramDir @endlink fibre (q-gram directory).
 */

    template <typename TText, typename TSpec>
    inline typename Fibre<Index<TText, TSpec>, FibreDir>::Type &
    indexDir(Index<TText, TSpec> &index) {
        return getFibre(index, FibreDir());
    }
    template <typename TText, typename TSpec>
    inline typename Fibre<Index<TText, TSpec> const, FibreDir>::Type &
    indexDir(Index<TText, TSpec> const &index) {
        return getFibre(index, FibreDir());
    }

//////////////////////////////////////////////////////////////////////////////
/*!
 * @fn IndexQGram#dirAt
 * @headerfile <seqan/index.h>
 * @brief Shortcut for <tt>value(indexDir(index), position)</tt>.
 *
 * @signature TFibre dirAt(position, index);
 *
 * @param[in] index    The IndexQGram object holding the fibre.
 * @param[in] position A position in the array on which the value should be accessed.
 *
 * @return TFibre A reference to the @link QGramIndexFibres#QGramDir @endlink fibre.
 */

    template <typename TPos, typename TIndex>
    inline typename Reference<typename Fibre<TIndex, FibreDir>::Type>::Type
    dirAt(TPos i, TIndex &index) {
        return value(getFibre(index, FibreDir()), i);
    }
    template <typename TPos, typename TIndex>
    inline typename Reference<typename Fibre<TIndex const, FibreDir>::Type>::Type
    dirAt(TPos i, TIndex const &index) {
        return value(getFibre(index, FibreDir()), i);
    }


//////////////////////////////////////////////////////////////////////////////
/*!
 * @fn IndexQGram#indexCounts
 * @headerfile <seqan/index.h>
 * @brief Shortcut for <tt>getFibre(index, QGramCounts())</tt>.
 *
 * @signature TFibre indexCounts(index);
 *
 * @param[in] index The @link IndexQGram @endlink object holding the fibre.
 *
 * @return TFibre A reference to the @link QGramIndexFibres#QGramCounts @endlink fibre (counts array).
 */
    template <typename TText, typename TSpec>
    inline typename Fibre<Index<TText, TSpec>, FibreCounts>::Type &
    indexCounts(Index<TText, TSpec> &index) {
        return getFibre(index, FibreCounts());
    }
    template <typename TText, typename TSpec>
    inline typename Fibre<Index<TText, TSpec> const, FibreCounts>::Type &
    indexCounts(Index<TText, TSpec> const &index) {
        return getFibre(index, FibreCounts());
    }

//////////////////////////////////////////////////////////////////////////////
/*!
 * @fn IndexQGram#indexCountsDir
 * @headerfile <seqan/index.h>
 * @brief Shortcut for <tt>getFibre(index, QGramCountsDir())</tt>.
 *
 * @signature TFibre indexCountsDir(index);
 *
 * @param[in] index The IndexQGram object holding the fibre.
 *
 * @return TFibre A reference to the @link QGramIndexFibres#QGramCountsDir @endlink fibre (counts directory).
 */

    template <typename TText, typename TSpec>
    inline typename Fibre<Index<TText, TSpec>, FibreCountsDir>::Type &
    indexCountsDir(Index<TText, TSpec> &index) {
        return getFibre(index, FibreCountsDir());
    }
    template <typename TText, typename TSpec>
    inline typename Fibre<Index<TText, TSpec> const, FibreCountsDir>::Type &
    indexCountsDir(Index<TText, TSpec> const &index) {
        return getFibre(index, FibreCountsDir());
    }

//////////////////////////////////////////////////////////////////////////////
/*!
 * @fn IndexQGram#indexBucketMap
 * @headerfile <seqan/index.h>
 * @brief Shortcut for <tt>getFibre(index, QGramBucketMap())</tt>.
 *
 * @signature TFibre indexBucketMap(index);
 *
 * @param[in] index The IndexQGram object holding the fibre.
 *
 * @return TFibre A reference to the @link QGramIndexFibres#QGramBucketMap @endlink fibre (maps q-gram hashes to
 *                buckets).
 */

    template <typename TText, typename TSpec>
    inline typename Fibre<Index<TText, TSpec>, FibreBucketMap>::Type &
    indexBucketMap(Index<TText, TSpec> &index) {
        return getFibre(index, FibreBucketMap());
    }
    template <typename TText, typename TSpec>
    inline typename Fibre<Index<TText, TSpec> const, FibreBucketMap>::Type &
    indexBucketMap(Index<TText, TSpec> const &index) {
        return getFibre(index, FibreBucketMap());
    }

//////////////////////////////////////////////////////////////////////////////
/*!
 * @fn IndexQGram#indexShape
 * @headerfile <seqan/index.h>
 * @brief Shortcut for <tt>getFibre(index, QGramShape())</tt>.
 * @signature TFibre indexShape(index);
 *
 * @param[in] index The @link Index @endlink object holding the fibre. Types: @link IndexQGram @endlink
 *
 * @return TFibre Returns a reference to the @link Shape @endlink object of a q-gram index.  Formally, this is a
 *                reference to the @link QGramIndexFibres#QGramShape @endlink fibre.
 */
    template <typename TText, typename TSpec>
    inline typename Fibre<Index<TText, TSpec>, FibreShape>::Type &
    indexShape(Index<TText, TSpec> &index) {
        return getFibre(index, FibreShape());
    }
    template <typename TText, typename TSpec>
    inline typename Fibre<Index<TText, TSpec> const, FibreShape>::Type &
    indexShape(Index<TText, TSpec> const &index) {
        return getFibre(index, FibreShape());
    }

//////////////////////////////////////////////////////////////////////////////
/*!
 * @fn IndexQGram#getStepSize
 * @headerfile <seqan/index.h>
 * @brief Return the <i>q</i>-gram step size used for index creation.
 * @signature TSize getStepSize(index);
 *
 * @param[in] index A IndexQGram object.
 *
 * @return TSize The step size of type @link Index#Size @endlink.  If <i>x</i> is returned every <i>x</i><sup>th</sup>
 *               <i>q</i>-gram is stored in the index.
 *
 * @see IndexQGram#setStepSize
 */
    template <typename TText, typename TShapeSpec, typename TSpec>
    inline typename Size<Index<TText, IndexQGram<TShapeSpec, TSpec> > const>::Type
    getStepSize(Index<TText, IndexQGram<TShapeSpec, TSpec> > const &index)
    {
        return (index.stepSize != 0)? index.stepSize: length(indexShape(index));
    }

//////////////////////////////////////////////////////////////////////////////
/*!
 * @fn IndexQGram#setStepSize
 * @headerfile <seqan/index.h>
 * @brief Change the <i>q</i>-gram step size used for index creation.
 *
 * @signature void setStepSize(index, stepSize);
 *
 * @param[in,out] index    The IndexQGram to modify.
 * @param[in]     stepSize The step size, @link IntegerConcept @endlink;  Each <tt>stepSize</tt><sup>th</sup>
 *                         <i>q</i>-gram will be stored.
 *
 * The default step size of a <i>q</i>-gram index is 1, which corresponds to all overlapping <i>q</i>-grams.  To take
 * effect of changing the <tt>stepSize</tt> the <i>q</i>-gram index should be empty or recreated.
 *
 * A <tt>stepSize</tt> of 0 corresponds to <tt>stepSize = length(indexShape(index))</tt>, i.e. all non-overlapping
 * q-grams.
 *
 * @see IndexQGram#getStepSize
 */

    template <typename TText, typename TShapeSpec, typename TSpec, typename TSize>
    inline void
    setStepSize(Index<TText, IndexQGram<TShapeSpec, TSpec> > &index, TSize stepSize)
    {
        index.stepSize = stepSize;
    }

//////////////////////////////////////////////////////////////////////////////

    template <typename TValue, typename TSpec>
    inline __int64 _fullDirLength(Shape<TValue, TSpec> const &shape)
    {
        return _intPow((__int64)ValueSize<TValue>::VALUE, weight(shape)) + 1;
    }

    template <typename TValue, typename TSpec>
    inline __int64 _fullDir2Length(Shape<TValue, TSpec> const &shape)
    {
        return (_intPow(
                    (__int64)ValueSize<TValue>::VALUE,
                    weight(shape) + 1) - 1)
                / ((unsigned)ValueSize<TValue>::VALUE - 1) + 1;
    }

    template <typename TIndex>
    inline __int64 _fullDirLength(TIndex const &index)
    {
        return _fullDirLength(indexShape(index));
    }

    template <typename TIndex>
    inline __int64 _fullDir2Length(TIndex const &index)
    {
        return _fullDir2Length(indexShape(index));
    }


    //////////////////////////////////////////////////////////////////////////////
    // QGramLess
    //
    // compare two q-grams of a given text (q-grams can be smaller than q)
    template < typename TSAValue, typename TText >
    struct QGramLess_ :
        public std::binary_function < TSAValue, TSAValue, bool >
    {
        typedef typename Iterator<TText, Standard>::Type TIter;
        TIter _begin, _end;
        typename Size<TText>::Type _q;

        template <typename TSize>
        QGramLess_(TText &text, TSize q):
            _begin(begin(text, Standard())),
            _end(end(text, Standard())),
            _q(q) {}

        // skip the first <offset> characters
        template <typename TSize1, typename TSize2>
        QGramLess_(TText &text, TSize1 q, TSize2 offset):
            _begin(begin(text, Standard()) + offset),
            _end(end(text, Standard())),
            _q(q) {}

        inline bool operator() (TSAValue const a, TSAValue const b) const
        {
            if (a == b) return false;
            TIter itA = _begin + a;
            TIter itB = _begin + b;
            if (a <= b) {
                TIter itEnd = itB + _q;
                if (itEnd > _end)
                    itEnd = _end;
                for(; itB != itEnd; ++itB, ++itA) {
                    if (ordLess(*itA, *itB)) return true;
                    if (ordLess(*itB, *itA)) return false;
                }
                return false;
            } else {
                TIter itEnd = itA + _q;
                if (itEnd > _end)
                    itEnd = _end;
                for(; itA != itEnd; ++itA, ++itB) {
                    if (ordLess(*itA, *itB)) return true;
                    if (ordLess(*itB, *itA)) return false;
                }
                return true;
            }
        }
    };

    template < typename TSAValue, typename TString, typename TSpec >
    struct QGramLess_<TSAValue, StringSet<TString, TSpec> const> :
        public std::binary_function < TSAValue, TSAValue, bool >
    {
        typedef typename Iterator<TString const, Standard>::Type TIter;
        typedef typename Size<TString>::Type                     TSize;
        typedef StringSet<TString, TSpec>                        TStringSet;
        TStringSet const &_stringSet;
        typename Size<TString>::Type _q;

        template <typename TSize>
        QGramLess_(StringSet<TString, TSpec> const &text, TSize q):
            _stringSet(text),
            _q(q) {}

        inline bool operator() (TSAValue const &a, TSAValue const &b) const
        {
            if (a == b) return false;

            typename Size<TStringSet>::Type seqNoA = getSeqNo(a, stringSetLimits(_stringSet));
            typename Size<TStringSet>::Type seqNoB = getSeqNo(b, stringSetLimits(_stringSet));
            TIter itA = begin(_stringSet[seqNoA], Standard()) + getSeqOffset(a, stringSetLimits(_stringSet));
            TIter itB = begin(_stringSet[seqNoB], Standard()) + getSeqOffset(b, stringSetLimits(_stringSet));
            TIter itAEnd = end(_stringSet[seqNoA], Standard());
            TIter itBEnd = end(_stringSet[seqNoB], Standard());

            if (itAEnd - itA < itBEnd - itB)
            {
                TIter itQEnd = itA + _q;
                TIter itEnd = _min(itQEnd, itAEnd);
                for(; itA < itEnd; ++itA, ++itB) {
                    if (ordLess(*itA, *itB)) return true;
                    if (ordLess(*itB, *itA)) return false;
                }
                // if qgram a is shorter than b => a < b
                if (itA != itQEnd) return true;
            }
            else
            {
                TIter itQEnd = itB + _q;
                TIter itEnd = _min(itQEnd, itBEnd);
                for(; itB < itEnd; ++itA, ++itB) {
                    if (ordLess(*itA, *itB)) return true;
                    if (ordLess(*itB, *itA)) return false;
                }
                // if qgram b is shorter or equal than a => a >= b
                if (itB != itQEnd) return false;
            }
            if (seqNoA < seqNoB) return true;
            if (seqNoA > seqNoB) return false;
            return suffixLength(a, _stringSet) > suffixLength(b, _stringSet);
        }
    };


    //////////////////////////////////////////////////////////////////////////////
    // QGramLessOffset
    //
    // compare two q-grams of a given text and skip the first <offset> characters
    template < typename TSAValue, typename TText >
    struct QGramLessOffset_ :
        QGramLess_<TSAValue, TText>
    {
        template <typename TSize1, typename TSize2>
        QGramLessOffset_(TText &text, TSize1 q, TSize2 offset):
            QGramLess_<TSAValue, TText> (text, q, offset) {}
    };

    template < typename TSAValue, typename TString, typename TSpec >
    struct QGramLessOffset_<TSAValue, StringSet<TString, TSpec> const> :
        public std::binary_function < TSAValue, TSAValue, bool >
    {
        typedef typename Iterator<TString const, Standard>::Type TIter;
        typedef typename Size<TString>::Type                     TSize;
        typedef StringSet<TString, TSpec>                        TStringSet;
        typedef typename StringSetLimits<TStringSet>::Type       TLimits;

        TStringSet const &_stringSet;
        TLimits const &limits;
        TSize _q, _offset;

        template <typename TSize1, typename TSize2>
        QGramLessOffset_(StringSet<TString, TSpec> const &text, TSize1 q, TSize2 offset):
            _stringSet(text),
            limits(stringSetLimits(text)),
            _q(q),
            _offset(offset) {}

        inline bool operator() (TSAValue const &a, TSAValue const &b) const
        {
            if (a == b) return false;

            typename Size<TStringSet const>::Type seqNoA = getSeqNo(a, limits);
            typename Size<TStringSet const>::Type seqNoB = getSeqNo(b, limits);

            TIter itA    = begin(_stringSet[seqNoA], Standard()) + getSeqOffset(a, limits);
            TIter itAEnd =   end(_stringSet[seqNoA], Standard());
            TIter itB    = begin(_stringSet[seqNoB], Standard()) + getSeqOffset(b, limits);
            TIter itBEnd =   end(_stringSet[seqNoB], Standard());

            if (itAEnd - itA < itBEnd - itB)
            {
                itA += _offset;
                itB += _offset;
                TIter itQEnd = itA + _q;
                TIter itEnd = _min(itQEnd, itAEnd);
                for(; itA < itEnd; ++itA, ++itB) {
                    if (ordLess(*itA, *itB)) return true;
                    if (ordLess(*itB, *itA)) return false;
                }
                // if qgram a is shorter than b => a < b
                if (itA != itQEnd) return true;
            }
            else
            {
                itA += _offset;
                itB += _offset;
                TIter itQEnd = itB + _q;
                TIter itEnd = _min(itQEnd, itBEnd);
                for(; itB < itEnd; ++itA, ++itB) {
                    if (ordLess(*itA, *itB)) return true;
                    if (ordLess(*itB, *itA)) return false;
                }
                // if qgram b is shorter or equal than a => a >= b
                if (itB != itQEnd) return false;
            }
            if (seqNoA < seqNoB) return true;
            if (seqNoA > seqNoB) return false;
            return suffixLength(a, _stringSet) > suffixLength(b, _stringSet);
        }
    };


    //////////////////////////////////////////////////////////////////////////////
    // QGramLessNoCheck
    //
    // compare two q-grams of a given text (no check for q-grams smaller than q)
    template < typename TSAValue, typename TText >
    struct QGramLessNoCheck_ :
        public std::binary_function < TSAValue, TSAValue, bool >
    {
        typedef typename Iterator<TText, Standard>::Type TIter;
        TIter _begin;
        typename Size<TText>::Type _q;

        template <typename TSize>
        QGramLessNoCheck_(TText &text, TSize q):
            _begin(begin(text, Standard())),
            _q(q) {}

        // skip the first <offset> characters
        template <typename TSize>
        QGramLessNoCheck_(TText &text, TSize q, TSize offset):
            _begin(begin(text, Standard()) + offset),
            _q(q) {}

        inline bool operator() (TSAValue const a, TSAValue const b) const
        {
            if (a == b) return false;
            TIter itA = _begin + a;
            TIter itB = _begin + b;

            TIter itEnd = itA + _q;
            for(; itA != itEnd; ++itA, ++itB) {
                if (ordLess(*itA, *itB)) return true;
                if (ordLess(*itB, *itA)) return false;
            }
            return a < b;
        }
    };

    template < typename TSAValue, typename TString, typename TSpec >
    struct QGramLessNoCheck_<TSAValue, StringSet<TString, TSpec> const> :
        public std::binary_function < TSAValue, TSAValue, bool >
    {
        typedef typename Iterator<TString, Standard>::Type  TIter;
        typedef StringSet<TString, TSpec>                   TStringSet;
        TStringSet const &_stringSet;
        typename Size<TString>::Type _q;

        template <typename TSize>
        QGramLessNoCheck_(StringSet<TString, TSpec> const &text, TSize q):
            _stringSet(text),
            _q(q) {}

        inline bool operator() (TSAValue const &a, TSAValue const &b) const
        {
            if (a == b) return false;

            typename Size<TStringSet>::Type seqNoA = getSeqNo(a, stringSetLimits(_stringSet));
            typename Size<TStringSet>::Type seqNoB = getSeqNo(b, stringSetLimits(_stringSet));

            TIter itA = begin(_stringSet[seqNoA], Standard()) + getSeqOffset(a, stringSetLimits(_stringSet));
            TIter itB = begin(_stringSet[seqNoB], Standard()) + getSeqOffset(b, stringSetLimits(_stringSet));

            TIter itEnd = itA + _q;
            for(; itA != itEnd; ++itA, ++itB) {
                if (ordLess(*itA, *itB)) return true;
                if (ordLess(*itB, *itA)) return false;
            }
            if (seqNoA < seqNoB) return true;
            if (seqNoA > seqNoB) return false;
            return suffixLength(a, _stringSet) > suffixLength(b, _stringSet);
        }
    };


    //////////////////////////////////////////////////////////////////////////////
    // QGramLessNoCheckOffset_
    //
    // compare two q-grams of a given text and skip the first <offset> characters
    template < typename TSAValue, typename TText >
    struct QGramLessNoCheckOffset_: QGramLessNoCheck_<TSAValue, TText>
    {
        template <typename TSize1, typename TSize2>
        QGramLessNoCheckOffset_(TText &text, TSize1 q, TSize2 offset):
            QGramLessNoCheck_<TSAValue, TText> (text, q, offset) {}
    };

    template < typename TSAValue, typename TString, typename TSpec >
    struct QGramLessNoCheckOffset_<TSAValue, StringSet<TString, TSpec> const> :
        public std::binary_function < TSAValue, TSAValue, bool >
    {
        typedef typename Iterator<TString, Standard>::Type TIter;
        StringSet<TString, TSpec> const &_stringSet;
        typename Size<TString>::Type _q, _offset;

        template <typename TSize1, typename TSize2>
        QGramLessNoCheckOffset_(StringSet<TString, TSpec> const &text, TSize1 q, TSize2 offset):
            _stringSet(text),
            _q(q),
            _offset(offset) {}

        inline bool operator() (TSAValue const &a, TSAValue const &b) const
        {
            if (a == b) return false;
            TIter itA = begin(_stringSet[getValueI1(a)], Standard()) + getValueI2(a) + _offset;
            TIter itB = begin(_stringSet[getValueI1(b)], Standard()) + getValueI2(b) + _offset;
            if (a <= b) {
                TIter itEnd = itB + _q;
                for(; itB != itEnd; ++itB, ++itA) {
                    if (ordLess(*itA, *itB)) return true;
                    if (ordLess(*itB, *itA)) return false;
                }
                return false;
            } else {
                TIter itEnd = itA + _q;
                for(; itA != itEnd; ++itA, ++itB) {
                    if (ordLess(*itA, *itB)) return true;
                    if (ordLess(*itB, *itA)) return false;
                }
                return true;
            }
        }
    };

    //////////////////////////////////////////////////////////////////////////////
    // Counting sort - little helpers
    //

    // map hashes 1:1 to directory
    template < typename THashValue, typename TParallelTag >
    inline THashValue
    requestBucket(Nothing &, THashValue hash, Tag<TParallelTag>)
    {
        return hash;
    }

    template < typename THashValue, typename TParallelTag >
    inline THashValue
    getBucket(Nothing const &, THashValue hash, Tag<TParallelTag>)
    {
        return hash;
    }

    // backward compatibility wrappers
    template < typename TBucketMap, typename THashValue >
    inline THashValue
    requestBucket(TBucketMap &bucketMap, THashValue hash)
    {
        return requestBucket(bucketMap, hash, Serial());
    }

    template < typename TBucketMap, typename THashValue >
    inline THashValue
    getBucket(TBucketMap const &bucketMap, THashValue hash)
    {
        return getBucket(bucketMap, hash, Serial());
    }

    //////////////////////////////////////////////////////////////////////////////
    // Counting sort - Step 1: Clear directory
    template < typename TDir, typename TParallelTag >
    inline void
    _qgramClearDir(TDir &dir, Nothing &, Tag<TParallelTag> parallelTag)
    {
        arrayFill(begin(dir, Standard()), end(dir, Standard()), 0, parallelTag);
    }

    template < typename TDir, typename TBucketMap >
    inline void _qgramClearDir(TDir &dir, TBucketMap &bucketMap)
    {
        _qgramClearDir(dir, bucketMap, Serial());
    }


    //////////////////////////////////////////////////////////////////////////////
    // Counting sort - Step 2: Count q-grams
    template < typename TDir, typename TBucketMap, typename TText, typename TShape, typename TStepSize >
    inline void
    _qgramCountQGrams(TDir &dir, TBucketMap &bucketMap, TText const &text, TShape shape, TStepSize stepSize)
    {
    SEQAN_CHECKPOINT
        typedef typename Iterator<TText const, Standard>::Type    TIterator;
        typedef typename Value<TDir>::Type                        TSize;

        if (length(text) < length(shape) || empty(shape)) return;
        TSize num_qgrams = (length(text) - length(shape)) / stepSize + 1;

        TIterator itText = begin(text, Standard());
        ++dir[requestBucket(bucketMap, hash(shape, itText))];
        if (stepSize == 1)
            for(TSize i = 1; i < num_qgrams; ++i)
            {
                ++itText;
                ++dir[requestBucket(bucketMap, hashNext(shape, itText))];
            }
        else
            for(TSize i = 1; i < num_qgrams; ++i)
            {
                itText += stepSize;
                ++dir[requestBucket(bucketMap, hash(shape, itText))];
            }
    }

    template < typename TDir, typename TBucketMap, typename TString, typename TSpec, typename TShape, typename TStepSize >
    inline void
    _qgramCountQGrams(TDir &dir, TBucketMap &bucketMap, StringSet<TString, TSpec> const &stringSet, TShape shape, TStepSize stepSize)
    {
    SEQAN_CHECKPOINT
        typedef typename Iterator<TString const, Standard>::Type    TIterator;
        typedef typename Value<TDir>::Type                            TSize;

        if (empty(shape)) return;

        if (stepSize == 1)
            for(unsigned seqNo = 0; seqNo < length(stringSet); ++seqNo)
            {
                TString const &sequence = value(stringSet, seqNo);
                if (length(sequence) < length(shape)) continue;
                TSize num_qgrams = length(sequence) - length(shape) + 1;

                TIterator itText = begin(sequence, Standard());
                ++dir[requestBucket(bucketMap, hash(shape, itText))];
                for(TSize i = 1; i < num_qgrams; ++i)
                {
                    ++itText;
                    ++dir[requestBucket(bucketMap, hashNext(shape, itText))];
                }
            }
        else
            for(unsigned seqNo = 0; seqNo < length(stringSet); ++seqNo)
            {
                TString const &sequence = value(stringSet, seqNo);
                if (length(sequence) < length(shape)) continue;
                TSize num_qgrams = (length(sequence) - length(shape)) / stepSize + 1;

                TIterator itText = begin(sequence, Standard());
                for(TSize i = 0; i < num_qgrams; ++i)
                {
                    ++dir[requestBucket(bucketMap, hash(shape, itText))];
                    itText += stepSize;
                }
            }
    }

    //////////////////////////////////////////////////////////////////////////////
    // Counting sort - Step 3: Cumulative sum
    //
    // a disabled bucket 3,2,x,4   (x = disabled)               y_i=sum_{j<=i} x_j
    // results in        0,0,3,x,5 (_qgramCummulativeSum)       z_i=y_{i-2}
    // or in             0,3,5,5,9 (_qgramCummulativeSumAlt)    z_i=y_{i-1}

    //    3 2 1 x x 5
    // 0 0 3 5 x x 6 11
    //   0 3 5 6 6 6 11

    // First two entries are 0.
    // Step 4 increments the entries hash(qgram)+1 on-the-fly while filling the SA table.
    // After step 4 each entry (0..n-1) is the beginning of a qgram bucket.
    template < typename TDir, typename TWithConstraints >
    inline typename Value<TDir>::Type
    _qgramCummulativeSum(TDir &dir, TWithConstraints)
    {
    SEQAN_CHECKPOINT
        typedef typename Iterator<TDir, Standard>::Type TDirIterator;
        typedef typename Value<TDir>::Type              TSize;

        TDirIterator it = begin(dir, Standard());
        TDirIterator itEnd = end(dir, Standard());
        TSize prevDiff = 0, prev2Diff = 0, sum = 0;
        while (it != itEnd)
        {
            if (TWithConstraints::VALUE && prevDiff == (TSize)-1)
            {
                sum += prev2Diff;
                prev2Diff = 0;
                prevDiff = *it;
                *it = (TSize)-1;                                // disable bucket
            } else {
                sum += prev2Diff;
                prev2Diff = prevDiff;
                prevDiff = *it;
                *it = sum;
            }
            ++it;
        }
        return sum + prev2Diff;
    }

    // The first entry is 0.
    // This function is used when Steps 4 and 5 (fill SA, correct disabled buckets) are ommited.
    template < typename TDir, typename TWithConstraints >
    inline typename Value<TDir>::Type
    _qgramCummulativeSumAlt(TDir &dir, TWithConstraints const)
    {
    SEQAN_CHECKPOINT
        typedef typename Iterator<TDir, Standard>::Type TDirIterator;
        typedef typename Value<TDir>::Type              TSize;

        TDirIterator it = begin(dir, Standard());
        TDirIterator itEnd = end(dir, Standard());
        TSize prevDiff = 0, sum = 0;
        for (; it != itEnd; ++it)
        {
            if (!TWithConstraints::VALUE || prevDiff != (TSize)-1)
            {
                sum += prevDiff;
                prevDiff = *it;
                *it = sum;
            }
            else
            {
                prevDiff = *it;
                *it = (TSize)-1;    // disable bucket
            }
        }
        return sum + prevDiff;
    }

    //////////////////////////////////////////////////////////////////////////////
    // Counting sort - Step 4: Fill suffix array
    // w/o constraints
    template <
        typename TSA,
        typename TText,
        typename TShape,
        typename TDir,
        typename TBucketMap,
        typename TWithConstraints,
        typename TStepSize >
    inline void
    _qgramFillSuffixArray(
        TSA &sa,
        TText const &text,
        TShape shape,
        TDir &dir,
        TBucketMap &bucketMap,
        TStepSize stepSize,
        TWithConstraints const)
    {
    SEQAN_CHECKPOINT
        typedef typename Iterator<TText const, Standard>::Type    TIterator;
        typedef typename Value<TDir>::Type                        TSize;

        if (empty(shape) || length(text) < length(shape)) return;

        TSize num_qgrams = length(text) - length(shape) + 1;
        TIterator itText = begin(text, Standard());

        if (TWithConstraints::VALUE) {
            TSize bktNo = getBucket(bucketMap, hash(shape, itText)) + 1;            // first hash
            if (dir[bktNo] != (TSize)-1) sa[dir[bktNo]++] = 0;                        // if bucket is enabled
        } else
            sa[dir[getBucket(bucketMap, hash(shape, itText)) + 1]++] = 0;            // first hash

        if (stepSize == 1)
            for(TSize i = 1; i < num_qgrams; ++i)
            {
                ++itText;
                if (TWithConstraints::VALUE) {
                    TSize bktNo = getBucket(bucketMap, hashNext(shape, itText)) + 1;    // next hash
                    if (dir[bktNo] != (TSize)-1) sa[dir[bktNo]++] = i;                    // if bucket is enabled
                } else
                    sa[dir[getBucket(bucketMap, hashNext(shape, itText)) + 1]++] = i;    // next hash
            }
        else
            for(TSize i = stepSize; i < num_qgrams; i += stepSize)
            {
                itText += stepSize;
                if (TWithConstraints::VALUE) {
                    TSize bktNo = getBucket(bucketMap, hash(shape, itText)) + 1;        // next hash (we mustn't use hashNext here)
                    if (dir[bktNo] != (TSize)-1) sa[dir[bktNo]++] = i;                    // if bucket is enabled
                } else
                    sa[dir[getBucket(bucketMap, hash(shape, itText)) + 1]++] = i;        // next hash
            }
    }

    // multiple sequences
    template <
        typename TSA,
        typename TString,
        typename TSpec,
        typename TShape,
        typename TDir,
        typename TBucketMap,
        typename TStepSize,
        typename TWithConstraints >
    inline void
    _qgramFillSuffixArray(
        TSA &sa,
        StringSet<TString, TSpec> const &stringSet,
        TShape shape,
        TDir &dir,
        TBucketMap &bucketMap,
        TStepSize stepSize,
        TWithConstraints const)
    {
    SEQAN_CHECKPOINT
        typedef typename Iterator<TString const, Standard>::Type    TIterator;
        typedef typename Value<TDir>::Type                            TSize;

        if (empty(shape)) return;

        if (stepSize == 1)
            for(unsigned seqNo = 0; seqNo < length(stringSet); ++seqNo)
            {
                TString const &sequence = value(stringSet, seqNo);
                if (length(sequence) < length(shape)) continue;
                TSize num_qgrams = length(sequence) - length(shape) + 1;

                typename Value<TSA>::Type localPos;
                assignValueI1(localPos, seqNo);
                assignValueI2(localPos, 0);

                TIterator itText = begin(sequence, Standard());
                if (TWithConstraints::VALUE) {
                    TSize bktNo = getBucket(bucketMap, hash(shape, itText)) + 1;                    // first hash
                    if (dir[bktNo] != (TSize)-1) sa[dir[bktNo]++] = localPos;                        // if bucket is enabled
                } else
                    sa[dir[getBucket(bucketMap, hash(shape, itText)) + 1]++] = localPos;            // first hash

                for(TSize i = 1; i < num_qgrams; ++i)
                {
                    ++itText;
                    assignValueI2(localPos, i);
                    if (TWithConstraints::VALUE) {
                        TSize bktNo = getBucket(bucketMap, hashNext(shape, itText)) + 1;            // next hash
                        if (dir[bktNo] != (TSize)-1) sa[dir[bktNo]++] = localPos;                    // if bucket is enabled
                    } else
                        sa[dir[getBucket(bucketMap, hashNext(shape, itText)) + 1]++] = localPos;    // next hash
                }
            }
        else
            for(unsigned seqNo = 0; seqNo < length(stringSet); ++seqNo)
            {
                TString const &sequence = value(stringSet, seqNo);
                if (length(sequence) < length(shape)) continue;
                TSize num_qgrams = length(sequence) - length(shape) + 1;

                typename Value<TSA>::Type localPos;
                assignValueI1(localPos, seqNo);
                assignValueI2(localPos, 0);

                TIterator itText = begin(sequence, Standard());
                for(TSize i = 0; i < num_qgrams; i += stepSize)
                {
                    assignValueI2(localPos, i);
                    if (TWithConstraints::VALUE) {
                        TSize bktNo = getBucket(bucketMap, hash(shape, itText)) + 1;                // hash
                        if (dir[bktNo] != (TSize)-1) sa[dir[bktNo]++] = localPos;                    // if bucket is enabled
                    } else
                        sa[dir[getBucket(bucketMap, hash(shape, itText)) + 1]++] = localPos;        // hash
                    itText += stepSize;
                }
            }
    }

    //////////////////////////////////////////////////////////////////////////////
    // Step 5: Correct disabled buckets
    template < typename TDir >
    inline void
    _qgramPostprocessBuckets(TDir &dir)
    {
    SEQAN_CHECKPOINT
        typedef typename Iterator<TDir, Standard>::Type            TDirIterator;
        typedef typename Value<TDir>::Type                        TSize;

        TDirIterator it = begin(dir, Standard());
        TDirIterator itEnd = end(dir, Standard());
        TSize prev = 0;
        for (; it != itEnd; ++it)
            if (*it == (TSize)-1)   // end positions
                *it = prev;
            else
                prev = *it;
    }


//////////////////////////////////////////////////////////////////////////////
/*!
 * @fn IndexQGram#createQGramIndex
 * @headerfile <seqan/index.h>
 * @brief Builds a <i>q</i>-gram index on a sequence.
 *
 * @signature void createQGramIndex(index);
 * @signature void createQGramIndex(sa, dir, bucketMap, text, shape, stepSize); [DEPRECATED]
 *
 * @param[out] index     The IndexQGram to create.
 * @param[out] sa        The resulting list in which all <i>q</i>-grams are sorted alphabetically.
 * @param[out] dir       The resulting array that indicates at which position in index the corresponding <i>q</i>-grams
 * @param[in]  bucketMap Stores the <i>q</i>-gram hashes for the openaddressing hash maps, see
 *                       @link IndexQGram#indexBucketMap @endlink. If bucketMap is of the type @link Nothing @endlink
 *                       the <i>q</i>-gram hash determines the bucket address in the index.
 * @param[in]  text      The @link TextConcept @endlink object to build the index for.
 * @param[in]  shape     The shape to be used. Types: @link Shape @endlink
 *                       can be found.
 * @param[in]  stepSize  Store every <tt>stepSize</tt>'th <i>q</i>-gram in the index, @link IntegerConcept @endlink.
 *
 * The resulting <i>q</i>-gram <tt>index</tt> contains the sorted list of qgrams. For each <i>q</i>-gram <tt>dir</tt> contains the
 * first position in index that corresponds to this <i>q</i>-gram.
 *
 * @warning This function should not be called directly. Please use @link Index#indexCreate @endlink or @link
 *          Index#indexRequire @endlink.  The resulting tables must have appropriate size before calling this function.
 */
    template < typename TIndex >
    inline bool _qgramDisableBuckets(TIndex &)
    {
        return false;    // we disable no buckets by default
    }

    template < typename TIndex >
    void createQGramIndex(TIndex &index)
    {
    SEQAN_CHECKPOINT
        typename Fibre<TIndex, QGramText>::Type const &text      = indexText(index);
        typename Fibre<TIndex, QGramSA>::Type         &sa        = indexSA(index);
        typename Fibre<TIndex, QGramDir>::Type        &dir       = indexDir(index);
        typename Fibre<TIndex, QGramShape>::Type      &shape     = indexShape(index);
        typename Fibre<TIndex, QGramBucketMap>::Type  &bucketMap = index.bucketMap;

        // 1. clear counters
        _qgramClearDir(dir, bucketMap);

        // 2. count q-grams
        _qgramCountQGrams(dir, bucketMap, text, shape, getStepSize(index));

        if (_qgramDisableBuckets(index))
        {
            // 3. cumulative sum
            _qgramCummulativeSum(dir, True());

            // 4. fill suffix array
            _qgramFillSuffixArray(sa, text, shape, dir, bucketMap, getStepSize(index), True());

            // 5. correct disabled buckets
            _qgramPostprocessBuckets(dir);
        }
        else
        {
            // 3. cumulative sum
            _qgramCummulativeSum(dir, False());

            // 4. fill suffix array
            _qgramFillSuffixArray(sa, text, shape, dir, bucketMap, getStepSize(index), False());
        }
    }

    // DEPRECATED
    // better use createQGramIndex(index) (above)
    template <
        typename TSA,
        typename TDir,
        typename TBucketMap,
        typename TText,
        typename TShape,
        typename TStepSize >
    void createQGramIndex(
        TSA &sa,
        TDir &dir,
        TBucketMap &bucketMap,
        TText const &text,
        TShape &shape,
        TStepSize stepSize)
    {
    SEQAN_CHECKPOINT

        // 1. clear counters
        _qgramClearDir(dir, bucketMap);

        // 2. count q-grams
        _qgramCountQGrams(dir, bucketMap, text, shape, stepSize);

        // 3. cumulative sum
        __int64 res = _qgramCummulativeSum(dir, False());
        SEQAN_ASSERT_EQ(res, static_cast<__int64>(length(sa)));
        ignoreUnusedVariableWarning(res);

        // 4. fill suffix array
        _qgramFillSuffixArray(sa, text, shape, dir, bucketMap, stepSize, False());
    }

    // DEPRECATED
    // better use createQGramIndex(index) (above)
    template <
        typename TSA,
        typename TDir,
        typename TBucketMap,
        typename TText,
        typename TShape,
        typename TStepSize >
    void createQGramIndex(
        TSA &sa,
        TDir &dir,
        TBucketMap &bucketMap,
        TText const &text,
        TShape &shape)
    {
        createQGramIndex(sa, dir, bucketMap, text, shape, 1);
    }

//////////////////////////////////////////////////////////////////////////////
/*!
 * @fn IndexQGram#createQGramIndexSAOnly
 * @headerfile <seqan/index.h>
 * @brief Builds the suffix array of a q-gram index on a sequence.
 *
 * @signature void createQGramIndexSAOnly(sa, text, shape, stepSize)
 *
 * @param[out] sa       The resulting list in which all q-grams are sorted alphabetically, @link String @endlink
 *                      object of @link SAValue @endlink of the underlying @link TextConcept text @endlink.
 * @param[in]  text     The sequence, a @link TextConcept @endlink object.
 * @param[in]  shape    The @link Shape @endlink to be used. q is the length of this shape.
 * @param[in]  stepSize Store every <tt>stepSize</tt>'th q-gram in the index, @link IntegerConcept @endlink.
 *
 * @warning This function should not be called directly. Please use @link Index#indexCreate @endlink or @link
 *          Index#indexRequire @endlink. The resulting tables must have appropriate size before calling this function.
 */
    template <
        typename TSA,
        typename TText,
        typename TShape,
        typename TStepSize >
    void createQGramIndexSAOnly(
        TSA &sa,
        TText const &text,
        TShape &shape,
        TStepSize stepSize)
    {
    SEQAN_CHECKPOINT
        typedef typename Size<TSA>::Type TSize;
        typedef typename Iterator<TSA, Standard>::Type TIter;

        // 1. Fill suffix array with a permutation (the identity)
        TIter it = begin(sa, Standard());
        TIter itEnd = end(sa, Standard());
        TSize i = 0;
        for(; it != itEnd; ++it, i += stepSize)
            *it = i;

        // 2. Sort suffix array with quicksort
        TSize span = length(shape);
        if (i + span > length(text) + 1)
            std::sort(
                begin(sa, Standard()),
                end(sa, Standard()),
                QGramLess_<typename Value<TSA>::Type, TText const>(text, span));
        else
            std::sort(
                begin(sa, Standard()),
                end(sa, Standard()),
                QGramLessNoCheck_<typename Value<TSA>::Type, TText const>(text, span));
    }

    template <
        typename TSA,
        typename TString,
        typename TSpec,
        typename TShape,
        typename TStepSize >
    void createQGramIndexSAOnly(
        TSA &sa,
        StringSet<TString, TSpec> const &stringSet,
        TShape &shape,
        TStepSize stepSize)
    {
    SEQAN_CHECKPOINT
        typedef typename Iterator<TSA, Standard>::Type    TIter;
        typedef typename Value<TSA>::Type                TValue;
        typedef typename Size<TString>::Type            TSize;
        typedef StringSet<TString, TSpec>                TStringSet;

        // 1. Fill suffix array with a permutation (the identity)
        TIter it = begin(sa, Standard());
        TValue pair;
        unsigned int const q1 = length(shape) - 1;
        for(unsigned seqNo = 0; seqNo < length(stringSet); ++seqNo)
        {
            unsigned int const strlen = length(stringSet[seqNo]);
            if (strlen > q1) {
                assignValueI1(pair, seqNo);
                for (TSize i = 0; i < strlen - q1; ++it, i += stepSize) {
                    assignValueI2(pair, i);
                    *it = pair;
                }
            }
        }
        SEQAN_ASSERT_EQ(it, end(sa, Standard()));

        // 2. Sort suffix array with quicksort
        TSize q = length(shape);
        if (lengthSum(stringSet) == length(sa))
            std::sort(
                begin(sa, Standard()),
                end(sa, Standard()),
                QGramLess_<typename Value<TSA>::Type, TStringSet const>(stringSet, q));
        else
            std::sort(
                begin(sa, Standard()),
                end(sa, Standard()),
                QGramLessNoCheck_<typename Value<TSA>::Type, TStringSet const>(stringSet, q));
    }

    template <
        typename TSA,
        typename TText,
        typename TSize1,
        typename TSize2 >
    inline void
    _refineQGramIndexBucket(
        TSA &sa,
        TText const &text,
        TSize1 oldQ,
        TSize2 newQ)
    {
    SEQAN_CHECKPOINT
        if (length(sa) <= 1)
            return;

        std::sort(
            begin(sa, Standard()),
            end(sa, Standard()),
            QGramLessOffset_<typename Value<TSA>::Type, TText const>(text, newQ - oldQ, oldQ));
    }

    template <
        typename TSA,
        typename TDir,
        typename TText,
        typename TSize1,
        typename TSize2 >
    inline void _refineQGramIndex(
        TSA &sa,
        TDir const &dir,
        TText const &text,
        TSize1 oldQ,
        TSize2 newQ)
    {
        //typedef typename Size<TSA>::Type                        TSize;
        typedef typename Iterator<TSA, Standard>::Type          TIter;
        typedef typename Iterator<TDir const, Standard>::Type   TDirIter;

        if (newQ <= (TSize2)oldQ) return;

        if (length(dir) < 2) {
            std::sort(
                begin(sa, Standard()),
                end(sa, Standard()),
                QGramLess_<typename Value<TSA>::Type, TText const>(text, newQ));
            return;
        }

        // 1. Sort each bucket with quicksort and compare substrings s[i+oldQ..i+newQ)
        TDirIter dirIt = begin(dir, Standard());
        TDirIter dirItEnd = end(dir, Standard());
        TIter itBegin = begin(sa, Standard());
        TIter itBktBegin = itBegin + *dirIt;
        TIter itBktEnd;
        ++dirIt;
        for(; dirIt != dirItEnd; ++dirIt, itBktBegin = itBktEnd) {
            itBktEnd = itBegin + *dirIt;
            if (itBktEnd - itBktBegin < 2) continue;
            std::sort(
                itBktBegin,
                itBktEnd,
                QGramLessOffset_<typename Value<TSA>::Type, TText const>(text, newQ - oldQ, oldQ));
        }
    }


//////////////////////////////////////////////////////////////////////////////
/*!
 * @fn IndexQGram#createQGramIndexDirOnly
 * @headerfile <seqan/index.h>
 * @brief Builds the directory of a <i>q</i>-gram index on a sequence.
 *
 * @signature void createQGramIndexDirOnly(dir, bucketMap, text, shape, stepSize);
 *
 * @param[out] dir       The resulting array that indicates at which position in index the corresponding q-grams can be
 *                       found.
 * @param[out] bucketMap Stores the q-gram hashes for the openaddressing hash maps, see @link IndexQGram#indexBucketMap
 *                       @endlink.  If bucketMap is of the type @link Nothing @endlink the q-gram hash determines the
 *                       bucket address in the index.
 * @param[in]  text      The sequence, @link TextConcept @endlink.
 * @param[in]  stepSize  Store every <tt>stepSize</tt>'th q-gram in the index, @link IntegerConcept @endlink.
 * @param[in]  shape     The @link Shape @endlink to be used.
 *
 * The resulting <tt>index</tt> contains the sorted list of qgrams. For each possible q-gram pos contains the first
 * position in index that corresponds to this q-gram.
 *
 * @warning This function should not be called directly. Please use @link Index#indexCreate @endlink or @link
 *          Index#indexRequire @endlink. The resulting tables must have appropriate size before calling this function.
 */

    template <
        typename TDir,
        typename TBucketMap,
        typename TText,
        typename TShape,
        typename TStepSize >
    void createQGramIndexDirOnly(
        TDir &dir,
        TBucketMap &bucketMap,
        TText const &text,
        TShape &shape,
        TStepSize stepSize)
    {
    SEQAN_CHECKPOINT

        // 1. clear counters
        _qgramClearDir(dir, bucketMap);

        // 2. count q-grams
        _qgramCountQGrams(dir, bucketMap, text, shape, stepSize);

        // 3. cumulative sum (Step 4 is ommited)
        _qgramCummulativeSumAlt(dir, False());
    }


//////////////////////////////////////////////////////////////////////////////
/*!
 * @fn IndexQGram#createCountArray
 * @headerfile <seqan/index.h>
 * @brief Builds an index on a StringSet storing how often a <i>q</i>-gram occurs in each sequence.
 *
 * @signature void createCountsArray(counts, dir, bucketMap, stringSet, shape, stepSize);
 *
 * @param[out] counts    The resulting @link String @endlink of pairs <i>(seqNo, count)</i>.
 * @param[out] dir       The resulting array that indicates at which position in the count
 *                       table the corresponding a certain q-gram can be found.
 * @param[out] bucketMap Stores the q-gram hashes for the openaddressing hash maps,
 *                       see @link IndexQGram#indexBucketMap @endlink. If bucketMap is of the
 *                       type @link Nothing @endlink the q-gram hash determines the
 *                       bucket address in the index.
 * @param[in] stringSet  The @link StringSet @endlink.
 * @param[in] shape      The @link Shape @endlink to be used.
 * @param[in] stepSize   Store every <tt>stepSize</tt><sup>th</sup> <i>q</i>-gram in the index,
 *                       @link IntegerConcept @endlink.
 *
 * @warning This function should not be called directly.  Please use @link Index#indexCreate @endlink or @link
 *          Index#indexRequire @endlink.  The resulting tables must have appropriate size before calling this function.
 */
    template <
        typename TCounts,
        typename TDir,
        typename TBucketMap,
        typename TString,
        typename TSpec,
        typename TShape,
        typename TStepSize >
    void createCountsArray(
        TCounts &counts,
        TDir &dir,
        TBucketMap &bucketMap,
        StringSet<TString, TSpec> const &stringSet,
        TShape shape,
        TStepSize stepSize)
    {
        typedef typename Iterator<TString const, Standard>::Type    TIterator;
        //typedef typename Iterator<TDir, Standard>::Type                TDirIterator;
        typedef typename Value<TShape>::Type                        TValue;
        typedef typename Size<TString>::Type                        TSize;

        TDir lastSeqSeen;
        resize(lastSeqSeen, length(dir), Exact());

        // 1. clear counters
        _qgramClearDir(dir, bucketMap);
        if (empty(shape))
        {
            clear(counts);
            return;
        }
        arrayFill(begin(lastSeqSeen, Standard()), end(lastSeqSeen, Standard()), -1);

        // 2. count distinct sequences for each q-gram
        for(unsigned seqNo = 0; seqNo < length(stringSet); ++seqNo)
        {
            TString const &sequence = value(stringSet, seqNo);
            if (length(sequence) < length(shape)) continue;
            TSize num_qgrams = (length(sequence) - length(shape)) / stepSize + 1;

            TIterator itText = begin(sequence, Standard());
            TValue bktNo = requestBucket(bucketMap, hash(shape, itText));
            lastSeqSeen[bktNo] = seqNo;
            ++dir[bktNo];
            if (stepSize == 1)
                for(TSize i = 1; i < num_qgrams; ++i)
                {
                    ++itText;
                    bktNo = requestBucket(bucketMap, hashNext(shape, itText));
                    if (seqNo != lastSeqSeen[bktNo]) {
                        lastSeqSeen[bktNo] = seqNo;
                        ++dir[bktNo];
                    }
                }
            else
                for(TSize i = 1; i < num_qgrams; ++i)
                {
                    itText += stepSize;
                    bktNo = requestBucket(bucketMap, hash(shape, itText));
                    if (seqNo != lastSeqSeen[bktNo]) {
                        lastSeqSeen[bktNo] = seqNo;
                        ++dir[bktNo];
                    }
                }
        }

        // 3. cumulative sum
        resize(counts, _qgramCummulativeSum(dir, False()), Exact());

        // 4. fill count array
        arrayFill(begin(lastSeqSeen, Standard()), end(lastSeqSeen, Standard()), -1);
        for(unsigned seqNo = 0; seqNo < length(stringSet); ++seqNo)
        {
            TString const &sequence = value(stringSet, seqNo);
            if (length(sequence) < length(shape)) continue;
            TSize num_qgrams = (length(sequence) - length(shape)) / stepSize + 1;

            TIterator itText = begin(sequence, Standard());
            TValue bktNo = getBucket(bucketMap, hash(shape, itText));
            lastSeqSeen[bktNo] = seqNo;
            TSize pos = dir[bktNo + 1]++;
            counts[pos].i1 = seqNo;                // first hash
            counts[pos].i2 = 1;

            if (stepSize == 1)
                for(TSize i = 1; i < num_qgrams; ++i)
                {
                    ++itText;
                    bktNo = getBucket(bucketMap, hashNext(shape, itText));
                    if (seqNo == lastSeqSeen[bktNo])
                        ++(counts[dir[bktNo + 1] - 1].i2);
                    else {
                        lastSeqSeen[bktNo] = seqNo;
                        pos = dir[bktNo + 1]++;
                        counts[pos].i1 = seqNo;                // next hash
                        counts[pos].i2 = 1;
                    }
                }
            else
                for(TSize i = 1; i < num_qgrams; ++i)
                {
                    itText += stepSize;
                    bktNo = getBucket(bucketMap, hash(shape, itText));
                    if (seqNo == lastSeqSeen[bktNo])
                        ++(counts[dir[bktNo + 1] - 1].i2);
                    else {
                        lastSeqSeen[bktNo] = seqNo;
                        pos = dir[bktNo + 1]++;
                        counts[pos].i1 = seqNo;                // next hash
                        counts[pos].i2 = 1;
                    }
                }
        }
    }


//////////////////////////////////////////////////////////////////////////////
// create q-gram index of *one* sequence in external memory

    // *** COMPARATORS & MAPS ***

    template <typename InType, typename Result = int>
    struct _qgramComp : public std::binary_function<InType,InType,Result> {
        inline Result operator()(InType const &a, InType const &b) const
        {
            typedef typename Value<InType, 2>::Type TQGram;
            typename Value<TQGram>::Type const *sa = a.i2.i;
            typename Value<TQGram>::Type const *sb = b.i2.i;
            typename Value<TQGram>::Type const *saEnd = sa + length(a.i2);

            for(; sa != saEnd; ++sa, ++sb) {
                if (*sa == *sb) continue;
                return (*sa < *sb)? -1 : 1;
            }
            return posCompare(a.i1, b.i1);
        }
    };

    // optimized for bitvectors
    template <typename T1, typename TValue, unsigned _size, typename Result>
    struct _qgramComp< Pair<T1, Tuple<TValue, _size, BitPacked<> >, Pack >, Result > :
        public std::binary_function<
            Pair<T1, Tuple<TValue, _size, BitPacked<> >, Pack >,
            Pair<T1, Tuple<TValue, _size, BitPacked<> >, Pack >,
            Result> {
        inline Result operator()(
            const Pair<T1, Tuple<TValue, _size, BitPacked<> >, Pack > &a,
            const Pair<T1, Tuple<TValue, _size, BitPacked<> >, Pack > &b) const
        {
            if (a.i2 < b.i2) return -1;
            if (a.i2 > b.i2) return 1;
            return posCompare(a.i1, b.i1);
        }
    };

    template <
        typename TText_,
        typename TShapeSpec,
        typename TSpec >
    void createQGramIndexExt(Index<TText_, IndexQGram<TShapeSpec, TSpec> > &index)
    {
    SEQAN_CHECKPOINT
        typedef Index<TText_, IndexQGram<TShapeSpec, TSpec> >       TIndex;
        typedef typename Fibre<TIndex, QGramText>::Type             TText;
        typedef typename Fibre<TIndex, QGramSA>::Type               TSA;
        typedef typename Fibre<TIndex, QGramDir>::Type              TDir;
        typedef typename Fibre<TIndex, QGramShape>::Type            TShape;
        typedef typename Fibre<TIndex, QGramBucketMap>::Type        TBucketMap;

        TText const &text = indexText(index);
        TSA &sa = indexSA(index);
        TDir &dir = indexDir(index);
        TShape &shape = indexShape(index);
        TBucketMap &bucketMap = index.bucketMap;

        // signed characters behave different than unsigned when compared
        // to get the same index with signed or unsigned chars we simply cast them to unsigned
        // before feeding them into the pipeline
        typedef typename Value<TIndex>::Type                        TValue;
        typedef typename MakeUnsigned_<TValue>::Type                TUValue;

        // *** SPECIALIZATION ***

        typedef Pipe< TText, Source<> >                                TSource;
        typedef Pipe< TSource, Caster<TUValue, CasterConvert> >        TUnsigner;
        typedef Pipe< TUnsigner, Tupler<LENGTH<TShape>::VALUE> >    TTupler;
                                        typedef _qgramComp<TypeOf_(TTupler)> qcomp_t;
        typedef Pool<
                    TypeOf_(TTupler),
                    SorterSpec< SorterConfigSize<qcomp_t, TSizeOf_(TTupler) > >
                > TSortTuples;

        // *** INSTANTIATION ***

        TSource            src(text);
        TUnsigner        unsigner(src);
        TTupler            tupler(unsigner);
        TSortTuples        sorter;

        // sort q-grams
        sorter << tupler;

        // fill sa and dir
        if (!beginRead(sorter)) return;

        typename Iterator<TSA>::Type itSA = begin(sa);
        typename Iterator<TDir>::Type itDir = begin(dir);

        qcomp_t    qcomp;

        typename Size<TSortTuples>::Type leftToRead = length(sorter);
        typename Value<TShape>::Type code, old_code = 0;

        if (leftToRead > 0)
        {
            typename Value<TSortTuples>::Type old_qgram = *sorter;

            old_code = hash(shape, getValueI2(old_qgram));
            *itSA = (*sorter).i1;
            *itDir = 0;
            --leftToRead, ++sorter, ++itSA, ++itDir;

            for (; leftToRead > 0; --leftToRead, ++sorter, ++itSA)
            {
                // copy occurrence position
                *itSA = (*sorter).i1;
                if (qcomp(old_qgram, *sorter) != 0)
                {
                    old_qgram = *sorter;
                    code = hash(shape, getValueI2(old_qgram));

                    SEQAN_ASSERT_LT(old_code, code);

                    // copy bucket begin
                    typename Size<TSortTuples>::Type i = length(sorter) - leftToRead;
                    for(; old_code < code; ++old_code, ++itDir)
                        *itDir = i;
                }
            }
        }

        // fill bucket table
        typename Size<TSortTuples>::Type i = length(sorter);
        code = length(dir);
        for(; old_code < code; ++old_code, ++itDir)
            *itDir = i;

        endRead(sorter);
    }


//////////////////////////////////////////////////////////////////////////////
// create q-gram index of *multiple* sequences in external memory

    // uses a Sorter to (1) sort q-grams and (2) afterwards create the bucket directory
    template <
        typename TString,
        typename TSSSpec,
        typename TShapeSpec,
        typename TSpec >
    void createQGramIndexExt(Index<StringSet<TString, TSSSpec>, IndexQGram<TShapeSpec, TSpec> > &index)
    {
    SEQAN_CHECKPOINT
        typedef StringSet<TString, TSSSpec>                         TText_;
        typedef Index<TText_, IndexQGram<TShapeSpec, TSpec> >       TIndex;
        typedef typename Fibre<TIndex, QGramText>::Type             TText;
        typedef typename Fibre<TIndex, QGramSA>::Type               TSA;
        typedef typename Fibre<TIndex, QGramDir>::Type              TDir;
        typedef typename Fibre<TIndex, QGramShape>::Type            TShape;
        typedef typename Fibre<TIndex, QGramBucketMap>::Type        TBucketMap;

        TText const &stringSet = indexText(index);
        TSA &sa = indexSA(index);
        TDir &dir = indexDir(index);
        TShape &shape = indexShape(index);
        TBucketMap &bucketMap = index.bucketMap;

        typedef typename Value<TDir>::Type TSize;
        const TSize DISABLED_BUCKET = (TSize)-1;

        // signed characters behave different than unsigned when compared
        // to get the same index with signed or unsigned chars we simply cast them to unsigned
        // before feeding them into the pipeline
        typedef typename Concatenator<TText>::Type                  TConcat;
        typedef typename Value<TIndex>::Type                        TValue;
        typedef typename MakeUnsigned_<TValue>::Type                TUValue;
        typedef Multi<
            Tupler<LENGTH<TShape>::VALUE, true, BitPacked<> >,
            typename Value<TSA>::Type,
            typename StringSetLimits<TText>::Type >                 TTuplerSpec;

        // *** SPECIALIZATION ***

        typedef Pipe< TConcat, Source<> >                            TSource;
        typedef Pipe< TSource, Caster<TUValue, CasterConvert> >        TUnsigner;
        typedef Pipe< TUnsigner, TTuplerSpec >                        TTupler;
                                        typedef _qgramComp<TypeOf_(TTupler)> qcomp_t;
        typedef Pool<
                    TypeOf_(TTupler),
                    SorterSpec< SorterConfigSize<qcomp_t, TSizeOf_(TTupler) > >
                > TSortTuples;

        // *** INSTANTIATION ***

        TSource            src(concat(stringSet));
        TUnsigner        unsigner(src);
        TTupler            tupler(unsigner, stringSetLimits(stringSet));
        TSortTuples        sorter;

        if (empty(shape)) return;

#ifdef _OPENMP
std::cout << "start\n";
double start = omp_get_wtime();
#endif

        // 1. clear counters
        _qgramClearDir(dir, bucketMap);

        // 2. count q-grams
        if (!beginRead(tupler))
            return;

        for (; !eof(tupler); ++tupler)
            ++dir[requestBucket(bucketMap, hash(shape, getValueI2(*tupler)))];

        endRead(tupler);

#ifdef _OPENMP
std::cout << "scan took " << omp_get_wtime() - start << " seconds.\n";
#endif

        // 3a. optionally disable some buckets
        bool bucketsDisabled = _qgramDisableBuckets(index);

        // 3b. cumulative sum
        TSize qgramCount;
        if (bucketsDisabled)
            qgramCount = _qgramCummulativeSumAlt(dir, True());
        else
            qgramCount = _qgramCummulativeSumAlt(dir, False());

        // 4. fill suffix array
        resize(sorter, qgramCount);
        if (!beginRead(tupler) || !beginWrite(sorter))
            return;

#ifdef _OPENMP
std::cout << "start write\n";
start = omp_get_wtime();
#endif

        for (; !eof(tupler); ++tupler)
            if (dir[getBucket(bucketMap, hash(shape, getValueI2(*tupler))) + 1] != DISABLED_BUCKET)
                push(sorter, *tupler);

        endRead(tupler);
        endWrite(sorter);

#ifdef _OPENMP
std::cout << "writing took " << omp_get_wtime() - start << " seconds.\n";
start = omp_get_wtime();
#endif

        beginRead(sorter);
        typename Iterator<TSA, Standard>::Type saIt = begin(sa, Standard());
        for (; !eof(sorter); ++sorter, ++saIt)
            *saIt = getValueI1(*sorter);
        endRead(sorter);

#ifdef _OPENMP
std::cout << "reading took " << omp_get_wtime() - start << " seconds.\n";
#endif

        // 5. correct disabled buckets
        if (bucketsDisabled)
            _qgramPostprocessBuckets(dir);


/*
        // sort q-grams
        sorter << tupler;

        // fill sa and dir
        if (!beginRead(sorter)) return;

        typename Iterator<TSA>::Type itSA = begin(sa);
        typename Iterator<TDir>::Type itDir = begin(dir);

        qcomp_t    qcomp;

        typename Size<TSortTuples>::Type leftToRead = length(sorter);
        typename Value<TShape>::Type code, old_code = 0;

        if (leftToRead > 0)
        {
            typename Value<TSortTuples>::Type old_qgram = *sorter;

            old_code = hash(shape, getValueI2(old_qgram));
            *itSA = (*sorter).i1;
            *itDir = 0;
            --leftToRead, ++sorter, ++itSA, ++itDir;

            for (; leftToRead > 0; --leftToRead, ++sorter, ++itSA)
            {
                // copy occurrence position
                *itSA = (*sorter).i1;

                if (qcomp(old_qgram, *sorter) != 0)
                {
                    old_qgram = *sorter;
                    code = hash(shape, getValueI2(old_qgram));

                    // copy bucket begin
                    typename Size<TSortTuples>::Type i = length(sorter) - leftToRead;
                    for(; old_code < code; ++old_code, ++itDir)
                        *itDir = i;
                }
            }
        }

        // fill bucket table
        typename Size<TSortTuples>::Type i = length(sorter);
        code = length(dir);
        for(; old_code < code; ++old_code, ++itDir)
            *itDir = i;

        endRead(sorter);
*/
    }

    // uses a Mapper to (1) create the bucket directory and (2) afterwards map q-grams to their final position
    // dir must allow random accesses
    template < typename TIndex >
    void createQGramIndexExtSA(TIndex &index)
    {
    SEQAN_CHECKPOINT
        typedef typename Fibre<TIndex, QGramText>::Type             TText;
        typedef typename Fibre<TIndex, QGramSA>::Type               TSA;
        typedef typename Fibre<TIndex, QGramDir>::Type              TDir;
        typedef typename Fibre<TIndex, QGramShape>::Type            TShape;
        typedef typename Fibre<TIndex, QGramBucketMap>::Type        TBucketMap;

        TText const &stringSet = indexText(index);
        TSA &sa = indexSA(index);
        TDir &dir = indexDir(index);
        TShape &shape = indexShape(index);
        TBucketMap &bucketMap = index.bucketMap;

        typedef typename Value<TDir>::Type TSize;
        const TSize DISABLED_BUCKET = (TSize)-1;

        // signed characters behave different than unsigned when compared
        // to get the same index with signed or unsigned chars we simply cast them to unsigned
        // before feeding them into the pipeline
        typedef typename Concatenator<TText>::Type                  TConcat;
        typedef typename Value<TIndex>::Type                        TValue;
        typedef typename MakeUnsigned_<TValue>::Type                TUValue;
        typedef Multi<
            Tupler<LENGTH<TShape>::VALUE, true, BitPacked<> >,
            typename Value<TSA>::Type,
            typename StringSetLimits<TText>::Type >                 TTuplerSpec;

        // *** SPECIALIZATION ***

        typedef Pipe< TConcat, Source<> >                            TSource;
        typedef Pipe< TSource, Caster<TUValue, CasterConvert> >        TUnsigner;
        typedef Pipe< TUnsigner, TTuplerSpec >                        TTupler;
        typedef Pair<
            typename Value<TSA>::Type,
            typename Size<TSA>::Type,
            Pack >                                                  TPosWithRank;       // (pos_pair,rank)

        typedef Pool<
                    TPosWithRank,
                    MapperSpec< MapperConfig<filterI2<TPosWithRank> > >
                > TMapTuples;


        if (empty(shape)) return;

        // *** INSTANTIATION ***

        TSource            src(concat(stringSet));
        TUnsigner        unsigner(src);
        TTupler            tupler(unsigner, stringSetLimits(stringSet));
        TMapTuples      mapper;

#ifdef _OPENMP
std::cout << "start\n";
double start = omp_get_wtime();
#endif

        // 1. clear counters
        _qgramClearDir(dir, bucketMap);

        // 2. count q-grams
        if (!beginRead(tupler))
            return;

        for (; !eof(tupler); ++tupler)
            ++dir[requestBucket(bucketMap, hash(shape, getValueI2(*tupler)))];

        endRead(tupler);

#ifdef _OPENMP
std::cout << "scan took " << omp_get_wtime() - start << " seconds.\n";
#endif

        // 3a. optionally disable some buckets
        bool bucketsDisabled = _qgramDisableBuckets(index);

        // 3b. cumulative sum
        TSize qgramCount;
        if (bucketsDisabled)
            qgramCount = _qgramCummulativeSum(dir, True());
        else
            qgramCount = _qgramCummulativeSum(dir, False());


        // 4. fill suffix array
        resize(mapper, qgramCount);
        if (!beginRead(tupler) || !beginWrite(mapper))
            return;

#ifdef _OPENMP
std::cout << "start write\n";
start = omp_get_wtime();
#endif

        for (; !eof(tupler); ++tupler)
        {
            TSize bktNo = getBucket(bucketMap, hash(shape, getValueI2(*tupler))) + 1;
            if (dir[bktNo] != DISABLED_BUCKET)
                push(mapper, TPosWithRank(getValueI1(*tupler), dir[bktNo]++));
        }

        endRead(tupler);
        endWrite(mapper);

#ifdef _OPENMP
std::cout << "writing took " << omp_get_wtime() - start << " seconds.\n";
start = omp_get_wtime();
#endif

        beginRead(mapper);
        typename Iterator<TSA, Standard>::Type saIt = begin(sa, Standard());
        for (; !eof(mapper); ++mapper, ++saIt)
            *saIt = getValueI1(*mapper);
        endRead(mapper);

#ifdef _OPENMP
std::cout << "reading took " << omp_get_wtime() - start << " seconds.\n";
#endif

        // 5. correct disabled buckets
        if (bucketsDisabled)
            _qgramPostprocessBuckets(dir);
    }



//////////////////////////////////////////////////////////////////////////////
// interface for automatic index creation

    template <typename TText, typename TShape, typename TStepSize>
    inline typename LengthSum<TText const>::Type
    _qgramQGramCount(TText const &text, TShape const &shape, TStepSize stepSize)
    {
        typedef typename LengthSum<TText const>::Type   TSize;
        typedef typename Size<TText const>::Type        TSetSize;

        if (empty(shape)) return 0;

        // count all overlapping q-grams
        TSize qgramCount = 0;
        for(TSetSize i = 0; i < countSequences(text); ++i)
            if (sequenceLength(i, text) >= length(shape))
                qgramCount += (sequenceLength(i, text) - length(shape)) / stepSize + 1;
        return qgramCount;
    }

    template <typename TText, typename TShapeSpec, typename TSpec>
    inline typename Size<Index<TText, IndexQGram<TShapeSpec, TSpec> > >::Type
    _qgramQGramCount(Index<TText, IndexQGram<TShapeSpec, TSpec> > const &index)
    {
        return _qgramQGramCount(indexText(index), indexShape(index), getStepSize(index));
    }

    template <typename TText, typename TShapeSpec, typename TSpec>
    inline bool indexCreate(
        Index<TText, IndexQGram<TShapeSpec, TSpec> > &index,
        FibreSADir,
        Default const)
    {
        resize(indexSA(index), _qgramQGramCount(index), Exact());
        resize(indexDir(index), _fullDirLength(index), Exact());
        createQGramIndex(index);
        resize(indexSA(index), back(indexDir(index)), Exact());     // shrink if some buckets were disabled
        return true;
    }

    template <typename TText, typename TSpec>
    inline bool indexSupplied(Index<TText, TSpec> &index, FibreSADir) {
        return !(empty(getFibre(index, FibreSA())) || empty(getFibre(index, FibreDir())));
    }

    template <typename TText, typename TShapeSpec, typename TSpec>
    inline bool indexCreate(
        Index<TText, IndexQGram<TShapeSpec, TSpec> > &index,
        FibreSA,
        Default const)
    {
        resize(indexSA(index), _qgramQGramCount(index), Exact());
        createQGramIndexSAOnly(indexSA(index), indexText(index), indexShape(index), getStepSize(index));
        return true;
    }

    template <typename TText, typename TShapeSpec, typename TSpec>
    inline bool indexCreate(
        Index<TText, IndexQGram<TShapeSpec, TSpec> > &index,
        FibreCounts,
        Default const)
    {
        resize(indexCountsDir(index), _fullDirLength(index), Exact());
        createCountsArray(indexCounts(index), indexCountsDir(index), indexBucketMap(index), indexText(index), indexShape(index), getStepSize(index));
        return true;
    }

    template <typename TText, typename TShapeSpec, typename TSpec>
    inline bool indexCreate(
        Index<TText, IndexQGram<TShapeSpec, TSpec> > &index,
        FibreDir,
        Default const)
    {
        resize(indexDir(index), _fullDirLength(index), Exact());
        createQGramIndexDirOnly(indexDir(index), indexBucketMap(index), indexText(index), indexShape(index), getStepSize(index));
        return true;
    }


//////////////////////////////////////////////////////////////////////////////
/*!
 * @fn IndexQGram#getKmerSimilarityMatrix
 * @headerfile <seqan/index.h>
 * @brief Creates a matrix storing the number of common q-grams between all pairs of sequences.
 *
 * @signature void getKmerSimilarityMatrix(index, distMat[, seqSet]);
 *
 * @param[out] distMat The resulting q-gram similarity matrix.  Types: @link ContainerConcept @endlink
 * @param[in]  seqSet  Contains sequence numbers if only a subset of sequences should be compared.  Types: @link
 *                     ContainerConcept @endlink
 * @param[in]  index   The @link IndexQGram @endlink to use.
 *
 * <tt>distMat</tt> need to be a container of a floating point type and will be resized to <tt>seqCount * seqCount</tt>,
 * where <tt>seqCount</tt> is the number of sequences in the index/in <tt>seqSet</tt>.  The fraction of common q-grams
 * between sequence <tt>i</tt> and <tt>j</tt> is stored at position <tt>i*seqCount + j</tt>.  It sums up the minimum
 * number of q-gram occurrences between both sequences for each q-gram and normalizes it.
 */
    template < typename TText, typename TShapeSpec, typename TSpec, typename TDistMatrix >
    inline void getKmerSimilarityMatrix(
        Index< TText, IndexQGram<TShapeSpec, TSpec> > &index,
        TDistMatrix &distMat)
    {
        typedef Index< TText, IndexQGram<TShapeSpec,TSpec> >    TIndex;
        typedef typename Size<TIndex>::Type                        TSize;
        typedef typename Size<TDistMatrix>::Type                TSizeMat;
        typedef typename Value<TDistMatrix>::Type                TValueMat;

        typedef typename Fibre<TIndex, QGramCountsDir>::Type    TCountsDir;
        typedef typename Iterator<TCountsDir, Standard>::Type    TIterCountsDir;
        typedef typename Fibre<TIndex, QGramCounts>::Type        TCounts;
        typedef typename Iterator<TCounts, Standard>::Type        TIterCounts;

        // declare requirements
        indexRequire(index, QGramCounts());

        // initialize distance matrix
        TSizeMat seqNoLength = countSequences(index);
        clear(distMat);
        resize(distMat, seqNoLength * seqNoLength);
        arrayFill(begin(distMat, Standard()), end(distMat, Standard()), 0);

        TIterCountsDir itCountsDir = begin(indexCountsDir(index), Standard());
        TIterCountsDir itCountsDirEnd = end(indexCountsDir(index), Standard());
        TIterCounts itCountsBegin = begin(indexCounts(index), Standard());

        // for each bucket count common q-grams for each sequence pair
        TSize bucketBegin = *itCountsDir;
        for(++itCountsDir; itCountsDir != itCountsDirEnd; ++itCountsDir)
        {
            TSize bucketEnd = *itCountsDir;

            // q-gram must occur in at least 2 different sequences
            if (bucketBegin != bucketEnd)
            {
                TIterCounts itA = itCountsBegin + bucketBegin;
                TIterCounts itEnd = itCountsBegin + bucketEnd;
                for(; itA != itEnd; ++itA)
                {
                    TSizeMat ofs = (*itA).i1 * seqNoLength;
                    TSize countA = (*itA).i2;
                    TIterCounts itB = itA;

                    for(; itB != itEnd; ++itB)
                    {
                        TSize countB = (*itB).i2;
                        if (countA < countB)
                            distMat[ofs + (*itB).i1] += countA;
                        else
                            distMat[ofs + (*itB).i1] += countB;
                    }
                }
            }
            bucketBegin = bucketEnd;
        }

        // copy upper triangle to lower triangle and scale
        for(TSizeMat row = 0; row < seqNoLength; ++row)
        {
            TValueMat maxValRow = distMat[row * (seqNoLength + 1)];
            for(TSizeMat col = row + 1; col < seqNoLength; ++col)
            {
                // fractional common kmer count
                TValueMat maxValCol = distMat[col * (seqNoLength + 1)];
                TValueMat val = distMat[row * seqNoLength + col];

                // number of common q-grams / Number of possible common q-grams
                if (maxValRow < maxValCol) {
                    if (maxValRow != 0)
                        val /= maxValRow;
                    distMat[col * seqNoLength + row] = val;
                    distMat[row * seqNoLength + col] = val;
                } else {
                    if (maxValCol != 0)
                        val /= maxValCol;
                    distMat[col * seqNoLength + row] = val;
                    distMat[row * seqNoLength + col] = val;
                }
            }
        }

        // set diagonal to 1
        for(TSizeMat i = 0; i < seqNoLength; ++i)
            distMat[i * (seqNoLength + 1)] = 1;
    }


//////////////////////////////////////////////////////////////////////////////
// getKmerSimilarityMatrix
//     for a subset of the StringSet given as a sorted string of sequence numbers

    template <
        typename TText,
        typename TShapeSpec,
        typename TSpec,
        typename TDistMatrix,
        typename TSeqNoString
    >
    inline void getKmerSimilarityMatrix(
        Index< TText, IndexQGram<TShapeSpec, TSpec> > &index,
        TDistMatrix &distMat,
        TSeqNoString const &seqNo)
    {
        typedef Index< TText, IndexQGram<TShapeSpec,TSpec> >    TIndex;
        typedef typename Size<TIndex>::Type                        TSize;
        typedef typename Size<TDistMatrix>::Type                TSizeMat;
        typedef typename Value<TDistMatrix>::Type                TValueMat;

        typedef typename Fibre<TIndex, QGramCountsDir>::Type    TCountsDir;
        typedef typename Iterator<TCountsDir, Standard>::Type    TIterCountsDir;
        typedef typename Fibre<TIndex, QGramCounts>::Type        TCounts;
        typedef typename Iterator<TCounts, Standard>::Type        TIterCounts;
        typedef typename Iterator<TSeqNoString, Standard>::Type    TIterSeqNo;

        // declare requirements
        indexRequire(index, QGramCounts());

        // initialize distance matrix
        TSizeMat seqNoLength = length(seqNo);
        clear(distMat);
        resize(distMat, seqNoLength * seqNoLength);
        arrayFill(begin(distMat, Standard()), end(distMat, Standard()), 0);

        TIterCountsDir itCountsDir = begin(indexCountsDir(index), Standard());
        TIterCountsDir itCountsDirEnd = end(indexCountsDir(index), Standard());
        TIterCounts itCountsBegin = begin(indexCounts(index), Standard());
        TIterSeqNo itSetBegin = begin(seqNo, Standard());
        TIterSeqNo itSetEnd = end(seqNo, Standard());

        // for each bucket count common q-grams for each sequence pair
        TSize bucketBegin = *itCountsDir;
        for(++itCountsDir; itCountsDir != itCountsDirEnd; ++itCountsDir)
        {
            TSize bucketEnd = *itCountsDir;

            // q-gram must occur in at least 2 different sequences
            if (bucketBegin != bucketEnd)
            {
                TIterCounts itA = itCountsBegin + bucketBegin;
                TIterCounts itEnd = itCountsBegin + bucketEnd;
                TIterSeqNo itSetA = itSetBegin;

                while (itA != itEnd && itSetA != itSetEnd)
                {
                    if ((*itA).i1 < *itSetA)
                        ++itA;
                    else if ((*itA).i1 > *itSetA)
                        ++itSetA;
                    else
                    {
                        TSizeMat ofs = (itSetA - itSetBegin) * seqNoLength;
                        TSize countA = (*itA).i2;
                        TIterCounts itB = itA;
                        TIterSeqNo itSetB = itSetA;

                        while (itB != itEnd && itSetB != itSetEnd)
                        {
                            if ((*itB).i1 < *itSetB)
                                ++itB;
                            else if ((*itB).i1 > *itSetB)
                                ++itSetB;
                            else
                            {
                                TSize countB = (*itB).i2;
                                if (countA < countB)
                                    distMat[ofs + (itSetB - itSetBegin)] += countA;
                                else
                                    distMat[ofs + (itSetB - itSetBegin)] += countB;
                                ++itB;
                                ++itSetB;
                            }
                        }
                        ++itA;
                        ++itSetA;
                    }
                }
            }
            bucketBegin = bucketEnd;
        }

        // copy upper triangle to lower triangle and scale
        for(TSizeMat row = 0; row < seqNoLength; ++row)
        {
            TValueMat maxValRow = distMat[row * (seqNoLength + 1)];
            for(TSizeMat col = row + 1; col < seqNoLength; ++col)
            {
                // fractional common kmer count
                TValueMat maxValCol = distMat[col * (seqNoLength + 1)];
                TValueMat val = distMat[row * seqNoLength + col];

                // number of common q-grams / Number of possible common q-grams
                if (maxValRow < maxValCol) {
                    if (maxValRow != 0)
                        val /= maxValRow;
                    distMat[col * seqNoLength + row] = val;
                    distMat[row * seqNoLength + col] = val;
                } else {
                    if (maxValCol != 0)
                        val /= maxValCol;
                    distMat[col * seqNoLength + row] = val;
                    distMat[row * seqNoLength + col] = val;
                }
            }
        }

        // set diagonal to 1
        for(TSizeMat i = 0; i < seqNoLength; ++i)
            distMat[i * (seqNoLength + 1)] = 1;
    }


//////////////////////////////////////////////////////////////////////////////
/*!
 * @fn IndexQGram#range
 * @headerfile <seqan/index.h>
 * @brief Returns the suffix array interval borders of a q-gram in the index text.
 *
 * @signature TPair range(index, shape);
 *
 * @param[in] index  The IndexQGram object to query.
 * @param[in] shape  The @link Shape @endlink object to use for the query.
 *
 * @return TPair All positions where the q-gram stored in <tt>shape</tt> occurs in the text (see @link
 *               QGramIndexFibres#QGramText @endlink) are stored in a contiguous range of the suffix
 *               array.  <tt>range</tt> returns begin and end position of this range. If the type of <tt>index</tt> is
 *               <tt>TIndex</tt> the return type is <tt>Pair&lt;Size&lt;TIndex&gt;::Type&gt;</tt>.
 *
 * The necessary index tables are built on-demand via @link Index#indexRequire @endlink if index is not <tt>const</tt>.
 */


    template < typename TText, typename TShapeSpec, typename TSpec, typename TShapeSpec2, typename TValue >
    inline Pair<typename Size< Index< TText, IndexQGram<TShapeSpec, TSpec> > >::Type>
    range(
        Index< TText, IndexQGram<TShapeSpec, TSpec> > const &index,
        Shape< TValue, TShapeSpec2 > const &shape)
    {
        typedef typename Size< Index< TText, IndexQGram<TShapeSpec, TSpec> > >::Type TSize;
        typedef typename Size< typename Fibre< Index< TText, IndexQGram<TShapeSpec, TSpec> >, FibreDir>::Type>::Type TDirSize;
        TDirSize bucket = getBucket(indexBucketMap(index), value(shape));
        return Pair<TSize>(indexDir(index)[bucket], indexDir(index)[bucket + 1]);
    }

    template < typename TText, typename TShapeSpec, typename TSpec, typename TShapeSpec2, typename TValue >
    inline Pair<typename Size< Index< TText, IndexQGram<TShapeSpec, TSpec> > >::Type>
    range(
        Index< TText, IndexQGram<TShapeSpec, TSpec> > &index,
        Shape< TValue, TShapeSpec2 > const &shape)
    {
        indexRequire(index, QGramDir());
        return range(const_cast<Index< TText, IndexQGram<TShapeSpec, TSpec> > const &>(index), shape);
    }

//////////////////////////////////////////////////////////////////////////////
/*!
 * @fn IndexQGram#getOccurrence
 * @headerfile <seqan/index.h>
 * @brief Returns an occurrence a q-gram in the index text.
 *
 * @signature TSAValue getOccurrence(index, shape);
 *
 * @param[in] index The IndexQGram object to query.
 * @param[in] shape The @link Shape @endlink object to use for the query.  The shape stores the <i>q</i>-gram of the
 *                  last call ot @link Shape#hash @endlink or @link Shape#hashNext @endlink.
 *
 * @return TSAValue A position where the <i>q</i>-gram stored in <tt>shape</tt> occurs in the text (see @link
 *                  QGramIndexFibres#QGramText @endlink).  Type: @link SAValue @endlink of the index type.
 */


    template < typename TText, typename TShapeSpec, typename TSpec, typename TShapeSpec2, typename TValue >
    inline typename SAValue< Index< TText, IndexQGram<TShapeSpec, TSpec> > >::Type
    getOccurrence(
        Index< TText, IndexQGram<TShapeSpec, TSpec> > const &index,
        Shape< TValue, TShapeSpec2 > const &shape)
    {
        return saAt(indexDir(index)[getBucket(indexBucketMap(index), value(shape))], index);
    }

    template < typename TText, typename TShapeSpec, typename TSpec, typename TShapeSpec2, typename TValue >
    inline typename SAValue< Index< TText, IndexQGram<TShapeSpec, TSpec> > >::Type
    getOccurrence(
        Index< TText, IndexQGram<TShapeSpec, TSpec> > &index,
        Shape< TValue, TShapeSpec2 > const &shape)
    {
        indexRequire(index, QGramSADir());
        return getOccurrence(const_cast<Index< TText, IndexQGram<TShapeSpec, TSpec> > const &>(index), shape);
    }

//////////////////////////////////////////////////////////////////////////////
/*!
 * @fn IndexQGram#getOccurrences
 * @headerfile <seqan/index.h>
 * @brief Returns all occurrences of a q-gram in the index text.
 *
 * @signature TSAInfix getOccurrences(index, shape);
 *
 * @param[in] index A IndexQGram to query.
 * @param[in] shape A @link Shape @endlink object.
 *
 * @return TSAInfix All positions where the q-gram stored in <tt>shape</tt> occurs in the text (see @link
 *                 QGramIndexFibres#QGramText @endlink).  Tupes: @link SegmentableConcept#Infix @endlink<@link Fibre
 *                 @endlink<TIndex, QGramSA>::Type>::Type>.
 *
 * The necessary index tables are built on-demand via @link Index#indexRequire @endlink if index is not <tt>const</tt>.
 *
 * @see DemoMummy
 * @see DemoSupermaximalRepeats
 * @see DemoMaximalUniqueMatches
 *
 * @see VSTreeIterator#isUnique
 * @see VSTreeIterator#getFrequency
 * @see VSTreeIterator#isPartiallyLeftExtensible
 * @see VSTreeIterator#isLeftMaximal
 */

    template < typename TText, typename TShapeSpec, typename TSpec, typename TShapeSpec2, typename TValue >
    inline typename Infix< typename Fibre< Index< TText, IndexQGram<TShapeSpec, TSpec> >, FibreSA>::Type const >::Type
    getOccurrences(
        Index< TText, IndexQGram<TShapeSpec, TSpec> > const &index,
        Shape< TValue, TShapeSpec2 > const &shape)
    {
        typedef typename Size<typename Fibre< Index< TText, IndexQGram<TShapeSpec, TSpec> >, FibreDir>::Type>::Type TDirSize;
        TDirSize bucket = getBucket(indexBucketMap(index), value(shape));
        return infix(indexSA(index), indexDir(index)[bucket], indexDir(index)[bucket + 1]);
    }

    template < typename TText, typename TShapeSpec, typename TSpec, typename TShapeSpec2, typename TValue >
    inline typename Infix< typename Fibre< Index< TText, IndexQGram<TShapeSpec, TSpec> >, FibreSA>::Type const >::Type
    getOccurrences(
        Index< TText, IndexQGram<TShapeSpec, TSpec> > &index,
        Shape< TValue, TShapeSpec2 > const &shape)
    {
        indexRequire(index, QGramSADir());
        return getOccurrences(const_cast<Index< TText, IndexQGram<TShapeSpec, TSpec> > const &>(index), shape);
    }

//////////////////////////////////////////////////////////////////////////////
/*!
 * @fn IndexQGram#countOccurrences
 * @headerfile <seqan/index.h>
 * @brief Returns the number of occurrences of a q-gram in the index text.
 *
 * @signature TSize countOccurrences(index, shape);
 *
 * @param[in] index IndexQGram object to query.
 * @param[in] shape A @link Shape @endlink object to use for the query.
 *
 * @return TSize The number of positions where the q-gram stored in <tt>shape</tt> occurs in the text (see @link
 *               QGramIndexFibres#QGramText @endlink).  Metafunction: @link Index#Size @endlink.
 *
 * The necessary index tables are built on-demand via @link Index#indexRequire @endlink if index is not <tt>const</tt>.
 *
 * Demo: Demo.Supermaximal Repeats
 *
 * Demo: Demo.Index countChildren
 */

    //TODO(singer): Why Size of FibreSA and not Size<Index>?
    template < typename TText, typename TShapeSpec, typename TSpec, typename TShapeSpec2, typename TValue >
    inline typename Size< typename Fibre< Index< TText, IndexQGram<TShapeSpec, TSpec> >, FibreSA>::Type const >::Type
    countOccurrences(
        Index< TText, IndexQGram<TShapeSpec, TSpec> > const &index,
        Shape< TValue, TShapeSpec2 > const &shape)
    {
        typedef typename Size<typename Fibre< Index< TText, IndexQGram<TShapeSpec, TSpec> >, FibreDir>::Type>::Type TDirSize;
        TDirSize bucket = getBucket(indexBucketMap(index), value(shape));
        return indexDir(index)[bucket + 1] - indexDir(index)[bucket];
    }

    template < typename TText, typename TShapeSpec, typename TSpec, typename TShapeSpec2, typename TValue >
    inline typename Size< typename Fibre< Index< TText, IndexQGram<TShapeSpec, TSpec> >, FibreSA>::Type const >::Type
    countOccurrences(
        Index< TText, IndexQGram<TShapeSpec, TSpec> > &index,
        Shape< TValue, TShapeSpec2 > const &shape)
    {
        indexRequire(index, QGramDir());
        return countOccurrences(const_cast<Index< TText, IndexQGram<TShapeSpec, TSpec> > const &>(index), shape);
    }

//////////////////////////////////////////////////////////////////////////////
/*!
 * @fn IndexQGram#countOccurrencesMultiple
 * @headerfile <seqan/index.h>
 * @brief Returns the number of occurrences of a q-gram for every sequence of a @link StringSet @endlink .
 *
 * @signature TCountInfix countOccurrencesMultiple(index, shape);
 *
 * @param[in] index A IndexQGram of a @link StringSet @endlink.
 * @param[in] shape A @link Shape @endlink.
 *
 * @return TCountInfix A sequence of @link Pair pairs @endlink (seqNo,count), count &gt; 0.  For every @link StringSet
 *                     @endlink sequence the q-gram occurs in, seqNo is the sequence number and count the number of
 *                     occurrences.  If the type of <tt>index</tt> is <tt>TIndex</tt> the return type is
 *                     <tt>Infix&lt;Fibre&lt;TIndex, QGramCounts&gt;::Type const&gt;::Type</tt>.
 *
 * The necessary index tables are built on-demand via @link Index#indexRequire @endlink if index is not <tt>const</tt>.
 *
 * @see VSTreeIterator#countOccurrences
 */

    template < typename TText, typename TShapeSpec, typename TSpec, typename TShapeSpec2, typename TValue >
    inline typename Infix< typename Fibre< Index< TText, IndexQGram<TShapeSpec, TSpec> >, FibreCounts>::Type const >::Type
    countOccurrencesMultiple(
        Index< TText, IndexQGram<TShapeSpec, TSpec> > const &index,
        Shape< TValue, TShapeSpec2 > const &shape)
    {
        typedef typename Size<typename Fibre< Index< TText, IndexQGram<TShapeSpec, TSpec> >, FibreCountsDir>::Type>::Type TDirSize;
        TDirSize bucket = getBucket(indexBucketMap(index), value(shape));
        return infix(indexCounts(index), indexCountsDir(index)[bucket], indexCountsDir(index)[bucket + 1]);
    }

    template < typename TText, typename TShapeSpec, typename TSpec, typename TShapeSpec2, typename TValue >
    inline typename Infix< typename Fibre< Index< TText, IndexQGram<TShapeSpec, TSpec> >, FibreCounts>::Type const >::Type
    countOccurrencesMultiple(
        Index< TText, IndexQGram<TShapeSpec, TSpec> > &index,
        Shape< TValue, TShapeSpec2 > const &shape)
    {
        indexRequire(index, QGramCounts());
        return countOccurrencesMultiple(const_cast<Index< TText, IndexQGram<TShapeSpec, TSpec> > const &>(index), shape);
    }


//////////////////////////////////////////////////////////////////////////////
// clear

    template < typename TText, typename TShapeSpec, typename TSpec >
    inline void clear(Index<TText, IndexQGram<TShapeSpec, TSpec> > &index)
    {
        clear(getFibre(index, QGramSA()));
        clear(getFibre(index, QGramDir()));
        clear(getFibre(index, QGramCounts()));
        clear(getFibre(index, QGramCountsDir()));
    }

//////////////////////////////////////////////////////////////////////////////
// open

    template < typename TText, typename TShapeSpec, typename TSpec >
    inline bool open(
        Index< TText, IndexQGram<TShapeSpec, TSpec> > &index,
        const char *fileName,
        int openMode)
    {
        String<char> name;

        name = fileName;    append(name, ".txt");
        if ((!open(getFibre(index, QGramText()), toCString(name), openMode)) &&
            (!open(getFibre(index, QGramText()), fileName, openMode))) return false;

        name = fileName;    append(name, ".sa");
        if (!open(getFibre(index, QGramSA()), toCString(name), openMode)) return false;

        name = fileName;    append(name, ".dir");
        if (!open(getFibre(index, QGramDir()), toCString(name), openMode)) return false;

        return true;
    }
    template < typename TText, typename TShapeSpec, typename TSpec >
    inline bool open(
        Index< TText, IndexQGram<TShapeSpec, TSpec> > &index,
        const char *fileName)
    {
        return open(index, fileName, OPEN_RDONLY);
    }


//////////////////////////////////////////////////////////////////////////////
// save

    template < typename TText, typename TShapeSpec, typename TSpec >
    inline bool save(
        Index< TText, IndexQGram<TShapeSpec, TSpec> > &index,
        const char *fileName,
        int openMode)
    {
        String<char> name;

        name = fileName;    append(name, ".txt");
        if ((!save(getFibre(index, QGramText()), toCString(name), openMode)) &&
            (!save(getFibre(index, QGramText()), fileName, openMode))) return false;

        name = fileName;    append(name, ".sa");
        if (!save(getFibre(index, QGramSA()), toCString(name), openMode)) return false;

        name = fileName;    append(name, ".dir");
        if (!save(getFibre(index, QGramDir()), toCString(name), openMode)) return false;

        return true;
    }
    template < typename TText, typename TShapeSpec, typename TSpec >
    inline bool save(
        Index< TText, IndexQGram<TShapeSpec, TSpec> > &index,
        const char *fileName)
    {
        return save(index, fileName, OPEN_WRONLY | OPEN_CREATE);
    }

}

#endif //#ifndef SEQAN_HEADER_...
