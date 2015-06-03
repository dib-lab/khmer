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
// Author: Jochen Singer <jochen.singer@fu-berlin.de>
// ==========================================================================

//SEQAN_NO_DDDOC:do not generate documentation for this file


#ifndef INDEX_FM_COMPRESSED_SA_H_
#define INDEX_FM_COMPRESSED_SA_H_

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

template <typename TText, typename TSpec, typename TConfig>
struct CompressedSA;

template <typename TText, typename TSpec, typename TConfig>
struct LF;

struct FibreLF_;
typedef Tag<FibreLF_> const     FibreLF;

// ============================================================================
// Tags
// ============================================================================

/*!
 * @defgroup CompressedSAFibres CompressedSA Fibres
 * @brief Tag to select a specific fibre of a @link CompressedSA @endlink.
 *
 * @tag CompressedSAFibres#FibreSparseString
 * @brief The sparse string.
 *
 * @tag CompressedSAFibres#FibreLF
 * @brief A @link LF @endlink to recompute the missing values of the compressed suffix array.
 *
 * @see Fibre
 * @see CompressedSA#getFibre
 */
// TODO(esiragusa): Rename FibreSparseString as FibreSparseValues.
struct FibreSparseString_;
typedef Tag<FibreSparseString_> const FibreSparseString;

// ============================================================================
// Metafunctions
// ============================================================================

// ----------------------------------------------------------------------------
// Metafunction DefaultIndexStringSpec
// ----------------------------------------------------------------------------

template <typename TText, typename TSpec, typename TConfig>
struct DefaultIndexStringSpec<CompressedSA<TText, TSpec, TConfig> > :
    DefaultIndexStringSpec<TText> {};

// ----------------------------------------------------------------------------
// Metafunction Fibre
// ----------------------------------------------------------------------------

template <typename TText, typename TSpec, typename TConfig>
struct Fibre<CompressedSA<TText, TSpec, TConfig>, FibreSparseString>
{
    // TODO(esiragusa): Change SparseString spec to be SparseString<TValue, TSpec, TConfig>.
    typedef CompressedSA<TText, TSpec, TConfig>         TCSA;
    typedef typename SAValue<TText>::Type               TSAValue_;
    typedef typename DefaultIndexStringSpec<TCSA>::Type TSASpec_;
    typedef String<TSAValue_, TSASpec_>                 TSA_;
    typedef SparseString<TSA_, TConfig>                 Type;
};

template <typename TText, typename TSpec, typename TConfig>
struct Fibre<CompressedSA<TText, TSpec, TConfig>, FibreLF>
{
    typedef LF<TText, TSpec, TConfig>   Type;
};

// ----------------------------------------------------------------------------
// Metafunction Member
// ----------------------------------------------------------------------------

template <typename TText, typename TSpec, typename TConfig>
struct Member<CompressedSA<TText, TSpec, TConfig>, FibreLF>
{
    typedef typename Fibre<CompressedSA<TText, TSpec, TConfig>, FibreLF>::Type *    Type;
};

// ----------------------------------------------------------------------------
// Metafunction Reference
// ----------------------------------------------------------------------------

template <typename TText, typename TSpec, typename TConfig>
struct Reference<CompressedSA<TText, TSpec, TConfig> >
{
    // TODO(singer): We actually need a proxy here.
    typedef typename Value<CompressedSA<TText, TSpec, TConfig> >::Type Type;
};

template <typename TText, typename TSpec, typename TConfig>
struct Reference<CompressedSA<TText, TSpec, TConfig> const>
{
    // TODO(singer): We actually need a proxy here.
    typedef typename Value<CompressedSA<TText, TSpec, TConfig> >::Type /*const*/ Type;
};

// ----------------------------------------------------------------------------
// Metafunction Value
// ----------------------------------------------------------------------------

template <typename TText, typename TSpec, typename TConfig>
struct Value<CompressedSA<TText, TSpec, TConfig> >
{
    typedef typename Value<typename Fibre<CompressedSA<TText, TSpec, TConfig>, FibreSparseString>::Type>::Type   Type;
};

template <typename TText, typename TSpec, typename TConfig>
struct Value<CompressedSA<TText, TSpec, TConfig> const>
{
    typedef typename Value<typename Fibre<CompressedSA<TText, TSpec, TConfig> const, FibreSparseString>::Type>::Type    Type;
};

// ============================================================================
// Classes
// ============================================================================

// ----------------------------------------------------------------------------
// Class CompressedSA
// ----------------------------------------------------------------------------

/*!
 * @class CompressedSA
 * @implements ContainerConcept
 * @headerfile <seqan/index.h>
 *
 * @brief A suffix array storing only a few suffix array entries and computing the remaining on demand.
 *
 * @signature template <typename TText, typename TSpec, typename TConfig>
 *            class CompressedSA;
 *
 * @tparam TSpec Possibility to specialize a compressed suffix array. Default: void.
 * @tparam TText The type of the text the compressed suffix array is created from. Types: @link String @endlink, @link
 *               StringSet @endlink
 * @tparam TConfig A configuration object that can be used to change the types of the fibres easily. This possibility is
 *                 provided for convenience.
 *
 * @section Remarks
 *
 * The compressed suffix array can only be used together with a @link LF @endlink.
 */
template <typename TText, typename TSpec, typename TConfig>
struct CompressedSA
{
    typename Fibre<CompressedSA, FibreSparseString>::Type   sparseString;
    typename Member<CompressedSA, FibreLF>::Type            lf;

    CompressedSA() :
        lf()
    {}

    template <typename TLF>
    CompressedSA(TLF & lf)
    {
        setFibre(*this, lf, FibreLF());
    }

    template <typename TPos>
    SEQAN_HOST_DEVICE inline typename Value<CompressedSA>::Type const
    operator[](TPos pos)
    {
        return value(*this, pos);
    }

    template <typename TPos>
    SEQAN_HOST_DEVICE inline typename Value<CompressedSA>::Type
    operator[](TPos pos) const
    {
        return value(*this, pos);
    }
};

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function clear()
// ----------------------------------------------------------------------------

/*!
 * @fn CompressedSA#clear
 * @headerfile <seqan/index.h>
 * @brief Resets the compressed suffix array.
 *
 * @signature void clear(compressedSA);
 *
 * @param[in,out] compressesSA The compressed suffix array to be cleared.
 */

template <typename TText, typename TSpec, typename TConfig>
inline void clear(CompressedSA<TText, TSpec, TConfig> & compressedSA)
{
    clear(getFibre(compressedSA, FibreSparseString()));
//    clear(getFibre(compressedSA, FibreLF()));
}

// ----------------------------------------------------------------------------
// Function empty()
// ----------------------------------------------------------------------------

/*!
 * @fn CompressedSA#empty
 * @headerfile <seqan/index.h>
 * @brief Checks whether or not a compressed suffix array contains any elements.
 *
 * @signature bool empty(compressedSA);
 *
 * @param[in] compressesSA The compressed suffix array to be cleared.
 *
 * @return bool Returns true if the compressed suffix array is empty and false otherwise.
 */
template <typename TText, typename TSpec, typename TConfig>
SEQAN_HOST_DEVICE inline bool empty(CompressedSA<TText, TSpec, TConfig> & compressedSA)
{
    return empty(getFibre(compressedSA, FibreSparseString()));
    //    && empty(getFibre(compressedSA, FibreLF()));
}

// ----------------------------------------------------------------------------
// Function createCompressedSa()
// ----------------------------------------------------------------------------
// This function creates a compressed suffix array using a normal one.

// TODO(holtgrew): Rename to createCompressedSA
/*!
 * @fn CompressedSA#createCompressedSa
 * @headerfile <seqan/index.h>
 * @brief This function creates a compressed suffix array with a specified compression factor.
 *
 * @signature void createCompressedSa(compressedSA, completeSA, compressionFactor[, offset]);
 *
 * @param[out] compressedSA      The compressed suffix array.
 * @param[in]  completeSA        A complete suffix array containing all values. Types: @link String @endlink
 * @param[in]  compressionFactor The compression factor. A compression factor of x means that the compressed suffix array
 *                               specifically stores a value for every x values in the complete suffix array. Types: @link
 *                               UnsignedIntegerConcept @endlink
 * @param[in] offset             The offset determines how many empty values should be inserted into the compressed suffix array at the
 *                               beginning. This possibility accounts for the sentinel positions of the @link FMIndex @endlink.
 */

template <typename TText, typename TSpec, typename TConfig, typename TSA, typename TSize>
void createCompressedSa(CompressedSA<TText, TSpec, TConfig> & compressedSA, TSA const & sa, TSize offset)
{
    typedef CompressedSA<TText, TSpec, TConfig>                     TCompressedSA;
    typedef typename Size<TSA>::Type                                TSASize;
    typedef typename Fibre<TCompressedSA, FibreSparseString>::Type  TSparseSA;
    typedef typename Fibre<TSparseSA, FibreIndicators>::Type        TIndicators;
    typedef typename Fibre<TSparseSA, FibreValues>::Type            TValues;
    typedef typename Iterator<TSA const, Standard>::Type            TSAIter;

    TSparseSA & sparseString = getFibre(compressedSA, FibreSparseString());
    TIndicators & indicators = getFibre(sparseString, FibreIndicators());
    TValues & values = getFibre(sparseString, FibreValues());

    TSASize saLen = length(sa);
    resize(compressedSA, saLen + offset, Exact());

    TSAIter saIt = begin(sa, Standard());
    TSAIter saItEnd = end(sa, Standard());

    for (TSASize pos = 0; pos < offset; ++pos)
        setValue(indicators, pos, false);

    for (TSASize pos = offset; saIt != saItEnd; ++saIt, ++pos)
    {
        if (getSeqOffset(getValue(saIt)) % TConfig::SAMPLING == 0)
            setValue(indicators, pos, true);
        else
            setValue(indicators, pos, false);
    }
    updateRanks(indicators);

    resize(values, getRank(indicators, length(sparseString) - 1), Exact());

    saIt = begin(sa, Standard());
    for (TSASize pos = offset, counter = 0; saIt != saItEnd; ++saIt, ++pos)
    {
        if (getValue(indicators, pos))
        {
            assignValue(values, counter, getValue(saIt));
            ++counter;
        }
    }
}

template <typename TText, typename TSpec, typename TConfig, typename TSA>
void createCompressedSa(CompressedSA<TText, TSpec, TConfig> & compressedSA, TSA const & sa)
{
    createCompressedSa(compressedSA, sa, 0);
}

// ----------------------------------------------------------------------------
// Function getFibre()
// ----------------------------------------------------------------------------

/*!
 * @fn CompressedSA#getFibre
 *
 * @headerfile <seqan/index.h>
 *
 * @brief Returns a specific fibre of a compressed suffix array.
 *
 * @signature TFibre getFibre(compressedSA, fibreTag);
 *
 * @param[in]  fibreTag A tag that identifies the @link Fibre @endlink.  Types: @link CompressedSAFibres @endlink.
 * @param[out] compressedSA The container holding the fibre.
 *
 * @return TFibre A reference to the specified fibre of type @link Fibre @endlink&lt;CompressedSA&lt;TText, TSpec, TConfig&gt;, FibreSparseString&gt;::Type.
 */

template <typename TText, typename TSpec, typename TConfig>
SEQAN_HOST_DEVICE inline typename Fibre<CompressedSA<TText, TSpec, TConfig>, FibreSparseString>::Type const &
getFibre(CompressedSA<TText, TSpec, TConfig> const & compressedSA, FibreSparseString)
{
    return compressedSA.sparseString;
}

template <typename TText, typename TSpec, typename TConfig>
SEQAN_HOST_DEVICE inline typename Fibre<CompressedSA<TText, TSpec, TConfig>, FibreSparseString>::Type &
getFibre(CompressedSA<TText, TSpec, TConfig> & compressedSA, FibreSparseString)
{
    return compressedSA.sparseString;
}

template <typename TText, typename TSpec, typename TConfig>
inline typename Fibre<CompressedSA<TText, TSpec, TConfig>, FibreLF>::Type const &
getFibre(CompressedSA<TText, TSpec, TConfig> const & compressedSA, FibreLF)
{
    return value(compressedSA.lf);
}

template <typename TText, typename TSpec, typename TConfig>
inline typename Fibre<CompressedSA<TText, TSpec, TConfig>, FibreLF>::Type &
getFibre(CompressedSA<TText, TSpec, TConfig> & compressedSA, FibreLF)
{
    return value(compressedSA.lf);
}

// ----------------------------------------------------------------------------
// Function setFibre()
// ----------------------------------------------------------------------------

/*!
 * @fn CompressedSA#setFibre
 *
 * @headerfile <seqan/index.h>
 *
 * @brief Set the LF of the compressed suffix array. *
 * @signature void setFibre(compressedSa, lf, fibreLF);
 *
 * @param[in] compressesSa The compressed suffix array.
 * @param[in] lf The LF table to be used by the compressed suffix array.
 * @param[in] fibreLF A tag to specify the LF table Fibre
 */

template <typename TText, typename TSpec, typename TConfig, typename TLF>
void setFibre(CompressedSA<TText, TSpec, TConfig> & compressedSA, TLF & lf, FibreLF)
{
    setValue(compressedSA.lf, lf);
}

// ----------------------------------------------------------------------------
// Function length()
// ----------------------------------------------------------------------------

/*!
 * @fn CompressedSA#length
 *
 * @headerfile <seqan/index.h>
 *
 * @brief Returns the number of elements in the compressed suffix array.
 *
 * @signature TSize length(compressedSA);
 *
 * @param[in] compressesSA The compressed suffix array.
 *
 * @return TSize The number of elements in the compressed suffix array.  Types: The result of @link Size @endlink
 *               of the compressed suffix array.
 */

template <typename TText, typename TSpec, typename TConfig>
SEQAN_HOST_DEVICE inline typename Size<typename Fibre<CompressedSA<TText, TSpec, TConfig>, FibreSparseString>::Type>::Type
length(CompressedSA<TText, TSpec, TConfig> const & compressedSA)
{
    return length(getFibre(compressedSA, FibreSparseString()));
}

// ----------------------------------------------------------------------------
// Function resize()
// ----------------------------------------------------------------------------

/*!
 * @fn CompressedSA#resize
 *
 * @headerfile <seqan/index.h>
 *
 * @brief Resets the number of elements in the compressed suffix array.
 *
 * @signature TSize resize(compressedSA, newLength);
 *
 * @param[in,out] compressesSA The compressed suffix array.
 * @param[in]     newLength    The number of elements which should be stored in the compressed suffix array.
 *                             Types: @link UnsignedIntegerConcept @endlink.
 *
 * @return TSize The number of elements in the compressed suffix array.  Types: The result of @link Size @endlink
 *               of the compressed suffix array.
 *
 * If the new length is smaller than the actual one then the last <tt>x<tt> items of the compressed suffix array
 * are deleted with <tt>x = oldLength - newLength</tt>.
 */

template <typename TText, typename TSpec, typename TConfig, typename TSize, typename TExpand>
inline typename Size<typename Fibre<CompressedSA<TText, TSpec, TConfig>, FibreSparseString>::Type>::Type
resize(CompressedSA<TText, TSpec, TConfig> & compressedSA, TSize size, Tag<TExpand> tag)
{
    return resize(getFibre(compressedSA, FibreSparseString()), size, tag);
}

// ----------------------------------------------------------------------------
// Function value()
// ----------------------------------------------------------------------------

/*!
 * @fn CompressedSA#value
 *
 * @brief Returns the value stored at a specified position in the compressed suffix-array.
 *
 * @signature TValue value(compressedSA, pos);
 *
 * @param[in] compressedSA The compressed suffix array to access.
 * @param[in] pos          Position at which to access the suffix array. Types: @link UnsignedIntegerConcept @endlink.
 *
 * Note that the compressed suffix array is read only. Therefore a const reference is return by this function.
 */

template <typename TText, typename TSpec, typename TConfig, typename TPos>
SEQAN_HOST_DEVICE inline typename Value<CompressedSA<TText, TSpec, TConfig> >::Type
value(CompressedSA<TText, TSpec, TConfig> & compressedSA, TPos pos)
{
    typedef typename Fibre<CompressedSA<TText, TSpec, TConfig>, FibreSparseString>::Type     TSparseString;
    typedef typename Fibre<TSparseString, FibreIndicators>::Type    TIndicators;
    typedef typename Fibre<TSparseString, FibreValues>::Type        TValues;

    TIndicators const & indicators = getFibre(compressedSA.sparseString, FibreIndicators());
    TValues const & values = getFibre(compressedSA.sparseString, FibreValues());

    TPos counter = 0;
    for (; !getValue(indicators, pos); ++counter)
        pos = getFibre(compressedSA, FibreLF())(pos);

    return posAdd(getValue(values, getRank(indicators, pos) - 1), counter);
}

template <typename TText, typename TSpec, typename TConfig, typename TPos>
SEQAN_HOST_DEVICE inline typename Value<CompressedSA<TText, TSpec, TConfig> >::Type const
value(CompressedSA<TText, TSpec, TConfig> const & compressedSA, TPos pos)
{
    typedef typename Fibre<CompressedSA<TText, TSpec, TConfig>, FibreSparseString>::Type     TSparseString;
    typedef typename Fibre<TSparseString, FibreIndicators>::Type    TIndicators;
    typedef typename Fibre<TSparseString, FibreValues>::Type        TValues;

    TIndicators const & indicators = getFibre(compressedSA.sparseString, FibreIndicators());
    TValues const & values = getFibre(compressedSA.sparseString, FibreValues());

    TPos counter = 0;
    for (; !getValue(indicators, pos); ++counter)
        pos = getFibre(compressedSA, FibreLF())(pos);

    return posAdd(getValue(values, getRank(indicators, pos) - 1), counter);
}

// ----------------------------------------------------------------------------
// Function open()
// ----------------------------------------------------------------------------

/*!
 * @fn CompressedSA#open
 * @headerfile <seqan/index.h>
 * @brief This functions opens a compressed suffix array from disk.
 *
 * @signature bool open(compressedSA, fileName[, mode]);
 *
 * @param[in,out]  compressedSA The compressed suffix array to be opened.
 * @param[in]      fileName     <tt>char const *</tt> containing the file name.
 * @param[in]      mode         The combination of flags defining how the file should be opened.  To open a file
 *                              read-only, write-only or to read and write use <tt>OPEN_RDONLY</tt>,
 *                              <tt>OPEN_WRONLY</tt>, or <tt>OPEN_RDWR</tt>.To create or overwrite a file add
 *                              <tt>OPEN_CREATE</tt>.  To append a file if existing add <tt>OPEN_APPEND</tt>.  To
 *                              circumvent problems, files are always opened in binary mode.
 *                              Default: <tt>OPEN_RDWR | OPEN_CREATE | OPEN_APPEND</tt>
 *
 * @return bool <tt>true</tt> on success.
 */

template <typename TText, typename TSpec, typename TConfig>
inline bool open(CompressedSA<TText, TSpec, TConfig> & compressedSA, const char * fileName, int openMode)
{
    return open(getFibre(compressedSA, FibreSparseString()), fileName, openMode);
}

template <typename TText, typename TSpec, typename TConfig>
inline bool open(CompressedSA<TText, TSpec, TConfig> & compressedSA, const char * fileName)
{
    return open(compressedSA, fileName, DefaultOpenMode<CompressedSA<TText, TSpec, TConfig> >::VALUE);
}

// ----------------------------------------------------------------------------
// Function save()
// ----------------------------------------------------------------------------

/*!
 * @fn CompressedSA#save
 *
 * @headerfile <seqan/index.h>
 *
 * @brief This functions saves a compressed suffix array to disk.
 *
 * @signature bool save(compressedSA, fileName[, mode]);
 *
 * @param[in,out]  compressedSA The compressed suffix array to be opened.
 * @param[in]      fileName     <tt>char const *</tt> containing the file name.
 * @param[in]      mode         The combination of flags defining how the file should be opened.  To open a file
 *                              read-only, write-only or to read and write use <tt>OPEN_RDONLY</tt>,
 *                              <tt>OPEN_WRONLY</tt>, or <tt>OPEN_RDWR</tt>.To create or overwrite a file add
 *                              <tt>OPEN_CREATE</tt>.  To append a file if existing add <tt>OPEN_APPEND</tt>.  To
 *                              circumvent problems, files are always opened in binary mode.
 *                              Default: <tt>OPEN_RDWR | OPEN_CREATE | OPEN_APPEND</tt>
 *
 * @return bool <tt>true</tt> on success.
 */

template <typename TText, typename TSpec, typename TConfig>
inline bool save(CompressedSA<TText, TSpec, TConfig> const & compressedSA, const char * fileName, int openMode)
{
    return save(getFibre(compressedSA, FibreSparseString()), fileName, openMode);
}

template <typename TText, typename TSpec, typename TConfig>
inline bool save(CompressedSA<TText, TSpec, TConfig> const & compressedSA, const char * fileName)
{
    return save(compressedSA, fileName, DefaultOpenMode<CompressedSA<TText, TSpec, TConfig> >::VALUE);
}

}

#endif // INDEX_FM_COMPRESSED_SA_H_
