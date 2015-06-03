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
// Author: Enrico Siragusa <enrico.siragusa@fu-berlin.de>
// ==========================================================================
// LF stores all necessary information for the LF-mapping.
// ============================================================================

//SEQAN_NO_DDDOC:do not generate documentation for this file

#ifndef INDEX_FM_LF_TABLE_H_
#define INDEX_FM_LF_TABLE_H_

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

template <typename TText, typename TSpec, typename TConfig>
struct LF;

// ============================================================================
// Tags
// ============================================================================

/*!
 * @defgroup LFTableFibres LF Table Fibres
 * @brief Tag to select a specific fibre of a @link LF @endlink.
 *
 * These tags can be used to get @link Fibre Fibres @endlink of a @link LF @endlink.
 *
 * @see Fibre
 * @see Index#getFibre
 *
 * @tag LFTableFibres#FibrePrefixSums
 * @brief The prefix sum table of the lf table.
 *
 * @tag LFTableFibres#FibreBwt
 * @brief The occurrence table of the lf table.
 *
 * @tag LFTableFibres#FibreSentinels
 * @brief The type of the senitnels.
 */

struct FibrePrefixSums_;
struct FibreSentinels_;
struct FibreTempBwt_;

typedef Tag<FibrePrefixSums_>   const FibrePrefixSums;
typedef Tag<FibreSentinels_>    const FibreSentinels;
typedef Tag<FibreTempBwt_>      const FibreTempBwt;

// ============================================================================
// Metafunctions
// ============================================================================

// ----------------------------------------------------------------------------
// Metafunction Size
// ----------------------------------------------------------------------------

template <typename TText, typename TSpec, typename TConfig>
struct Size<LF<TText, TSpec, TConfig> >
{
    typedef typename TConfig::LengthSum    Type;
};

// ----------------------------------------------------------------------------
// Metafunction Value
// ----------------------------------------------------------------------------
// TODO(esiragusa): define Value<LF> == Size<LF> because LF[pos-i] = pos-j

template <typename TText, typename TSpec, typename TConfig>
struct Value<LF<TText, TSpec, TConfig> > : Value<TText> {};

template <typename TText, typename TSpec, typename TConfig>
struct Value<LF<TText, TSpec, TConfig> const> :
    public Value<LF<TText, TSpec, TConfig> > {};

template <typename TText, typename TSSetSpec, typename TSpec, typename TConfig>
struct Value<LF<StringSet<TText, TSSetSpec>, TSpec, TConfig> > : Value<TText> {};

template <typename TText, typename TSSetSpec, typename TSpec, typename TConfig>
struct Value<LF<StringSet<TText, TSSetSpec>, TSpec, TConfig> const> :
    public Value<LF<StringSet<TText, TSSetSpec>, TSpec, TConfig> > {};

// ----------------------------------------------------------------------------
// Metafunction Fibre<FibrePrefixSums>
// ----------------------------------------------------------------------------

template <typename TText, typename TSpec, typename TConfig>
struct Fibre<LF<TText, TSpec, TConfig>, FibrePrefixSums>
{
//    typedef typename Value<LF<TText, TSpec, TConfig> >::Type  TValue_;
//    typedef Tuple<TSize_, ValueSize<TValue_>::VALUE>          Type;

    typedef typename Size<LF<TText, TSpec, TConfig> >::Type TSize_;
    typedef typename DefaultIndexStringSpec<TText>::Type    TSpec_;
    typedef String<TSize_,  TSpec_>                         Type;
};

// ----------------------------------------------------------------------------
// Metafunction Fibre<FibreBwt>
// ----------------------------------------------------------------------------

template <typename TText, typename TSpec, typename TConfig>
struct Fibre<LF<TText, TSpec, TConfig>, FibreBwt>
{
    typedef typename Value<LF<TText, TSpec, TConfig> >::Type    TValue_;
    typedef RankDictionary<TValue_, typename TConfig::Bwt>      Type;
};

// ----------------------------------------------------------------------------
// Metafunction Fibre<FibreSentinels>
// ----------------------------------------------------------------------------

template <typename TText, typename TSpec, typename TConfig>
struct Fibre<LF<TText, TSpec, TConfig>, FibreSentinels>
{
    typedef typename Size<TText>::Type   Type;
};

template <typename TText, typename TSSetSpec, typename TSpec, typename TConfig>
struct Fibre<LF<StringSet<TText, TSSetSpec>, TSpec, TConfig>, FibreSentinels>
{
    typedef RankDictionary<bool, typename TConfig::Sentinels>   Type;
};

// ----------------------------------------------------------------------------
// Metafunction Fibre<FibreTempBwt>
// ----------------------------------------------------------------------------

template <typename TText, typename TSpec, typename TConfig>
struct Fibre<LF<TText, TSpec, TConfig>, FibreTempBwt>
{
    typedef typename Value<LF<TText, TSpec, TConfig> >::Type        TValue_;
    typedef String<TValue_, External<ExternalConfigLarge<> > >      Type;
};

// ============================================================================
// Classes
// ============================================================================

// ----------------------------------------------------------------------------
// Class LF
// ----------------------------------------------------------------------------

/*!
 * @class LF
 *
 * @headerfile <seqan/Index.h>
 *
 * @signature template <typename TText, typename TSpec, typename TConfig>
 *            struct LF;
 *
 * @brief LF is an object storing all necessary information for the LF-mapping.
 *
 * @tparam TText The type of the text the LF table is constructed from.
 * @tparam TSpec A possibility to specialize the LF table. Default: <tt>void</tt>
 * @tparam TConfig A configuration object for easily defining the LF table fibres.
 */

template <typename TText, typename TSpec, typename TConfig>
struct LF
{
    typename Fibre<LF, FibrePrefixSums>::Type   sums;
    typename Fibre<LF, FibreBwt>::Type          bwt;
    typename Fibre<LF, FibreSentinels>::Type    sentinels;
    typename Value<LF>::Type                    sentinelSubstitute;

    LF() : sentinels(), sentinelSubstitute(0)
    {}

    LF(TText const & text) : sentinels(), sentinelSubstitute(0)
    {
        createLF(text);
    }

    template <typename TPos>
    SEQAN_HOST_DEVICE typename Size<LF const>::Type
    operator[] (TPos pos) const
    {
        return _getBwtRank(*this, pos);
    }

    template <typename TPos>
    SEQAN_HOST_DEVICE typename Size<LF const>::Type
    operator() (TPos pos) const
    {
        return _getBwtRank(*this, pos);
    }

    template <typename TPos, typename TValue>
    SEQAN_HOST_DEVICE typename Size<LF const>::Type
    operator() (TPos pos, TValue val) const
    {
        return _getBwtRank(*this, pos, val);
    }
};

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function bwtLength()
// ----------------------------------------------------------------------------
// This function returns the length of the BWT of a text or a text collection.

template <typename TText>
inline typename LengthSum<TText>::Type
bwtLength(TText const & text)
{
    return lengthSum(text) + countSequences(text);
}

// ----------------------------------------------------------------------------
// Function getFibre()
// ----------------------------------------------------------------------------

/*!
 * @fn LF#getFibre
 *
 * @headerfile <seqan/index.h>
 *
 * @brief Returns a specific fibre of a LF table.
 *
 * @signature TFibre getFibre(lfTable, fibreTag);
 *
 * @param[in] fibreTag A tag that identifies the @link Fibre @endlink. Types: @link LFTableFibres @endlink
 * @param[in] lfTable  The LF table.
 *
 * @return TFibre A reference to the @link Fibre @endlink object of type @link Fibre @endlink&lt;@link LF @endlink&lt;TText, TSpec, TConfig&gt;, FibrePrefixSums&gt;::Type
 */
template <typename TText, typename TSpec, typename TConfig>
SEQAN_HOST_DEVICE inline typename Fibre<LF<TText, TSpec, TConfig>, FibrePrefixSums>::Type &
getFibre(LF<TText, TSpec, TConfig> & lf, FibrePrefixSums)
{
    return lf.sums;
}

template <typename TText, typename TSpec, typename TConfig>
SEQAN_HOST_DEVICE inline typename Fibre<LF<TText, TSpec, TConfig>, FibrePrefixSums>::Type const &
getFibre(LF<TText, TSpec, TConfig> const & lf, FibrePrefixSums)
{
    return lf.sums;
}

template <typename TText, typename TSpec, typename TConfig>
SEQAN_HOST_DEVICE inline typename Fibre<LF<TText, TSpec, TConfig>, FibreBwt>::Type &
getFibre(LF<TText, TSpec, TConfig> & lf, FibreBwt)
{
    return lf.bwt;
}

template <typename TText, typename TSpec, typename TConfig>
SEQAN_HOST_DEVICE inline typename Fibre<LF<TText, TSpec, TConfig>, FibreBwt>::Type const &
getFibre(LF<TText, TSpec, TConfig> const & lf, FibreBwt)
{
    return lf.bwt;
}

template <typename TText, typename TSpec, typename TConfig>
SEQAN_HOST_DEVICE inline typename Fibre<LF<TText, TSpec, TConfig>, FibreSentinels>::Type &
getFibre(LF<TText, TSpec, TConfig> & lf, FibreSentinels)
{
    return lf.sentinels;
}

template <typename TText, typename TSpec, typename TConfig>
SEQAN_HOST_DEVICE inline typename Fibre<LF<TText, TSpec, TConfig>, FibreSentinels>::Type const &
getFibre(LF<TText, TSpec, TConfig> const & lf, FibreSentinels)
{
    return lf.sentinels;
}

// ----------------------------------------------------------------------------
// Function empty()
// ----------------------------------------------------------------------------

/*!
 * @fn LF#empty
 *
 * @headerfile <seqan/index.h>
 *
 * @brief Clears the LF table.
 *
 * @signature bool empty(lfTable);
 *
 * @param[in] lfTable The LF table to be checked.
 *
 * @return bool <tt>true</tt> if the LF table is empty, <tt>false</tt> otherwise.
 */


template <typename TText, typename TSpec, typename TConfig>
SEQAN_HOST_DEVICE inline bool empty(LF<TText, TSpec, TConfig> const & lf)
{
    return empty(lf.bwt) &&
           empty(lf.sums);
}

template <typename TText, typename TSSetSpec, typename TSpec, typename TConfig>
SEQAN_HOST_DEVICE inline bool empty(LF<StringSet<TText, TSSetSpec>, TSpec, TConfig> const & lf)
{
    return empty(lf.bwt) &&
           empty(lf.sentinels) &&
           empty(lf.sums);
}

// ----------------------------------------------------------------------------
// Function clear()
// ----------------------------------------------------------------------------

/*!
 * @fn LF#clear
 *
 * @headerfile <seqan/index.h>
 *
 * @brief Resets the LF table.
 *
 * @signature void clear(lfTable);
 *
 * @param[in,out] lfTable The LF table to be cleared.
 */

template <typename TText, typename TSpec, typename TConfig>
inline void clear(LF<TText, TSpec, TConfig> & lf)
{
    clear(lf.bwt);
    _clearSentinels(lf);
    clear(lf.sums);
    lf.sentinelSubstitute = 0;
}

// ----------------------------------------------------------------------------
// Function _clearSentinels()
// ----------------------------------------------------------------------------

template <typename TText, typename TSpec, typename TConfig>
inline void _clearSentinels(LF<TText, TSpec, TConfig> & lf)
{
    lf.sentinels = 0;
}

template <typename TText, typename TSSetSpec, typename TSpec, typename TConfig>
inline void _clearSentinels(LF<StringSet<TText, TSSetSpec>, TSpec, TConfig> & lf)
{
    clear(lf.sentinels);
}

// ----------------------------------------------------------------------------
// Function _getSentinelsRank()
// ----------------------------------------------------------------------------

template <typename TText, typename TSpec, typename TConfig, typename TPos>
SEQAN_HOST_DEVICE inline typename Size<LF<TText, TSpec, TConfig> const>::Type
_getSentinelsRank(LF<TText, TSpec, TConfig> const & lf, TPos pos)
{
    return pos >= lf.sentinels;
}

template <typename TText, typename TSSetSpec, typename TSpec, typename TConfig, typename TPos>
SEQAN_HOST_DEVICE inline typename Size<LF<StringSet<TText, TSSetSpec>, TSpec, TConfig> const>::Type
_getSentinelsRank(LF<StringSet<TText, TSSetSpec>, TSpec, TConfig> const & lf, TPos pos)
{
    return getRank(lf.sentinels, pos);
}

// ----------------------------------------------------------------------------
// Function isSentinel()
// ----------------------------------------------------------------------------

template <typename TText, typename TSpec, typename TConfig, typename TPos>
SEQAN_HOST_DEVICE inline bool isSentinel(LF<TText, TSpec, TConfig> const & lf, TPos pos)
{
    return lf.sentinels == pos;
}

template <typename TText, typename TSSetSpec, typename TSpec, typename TConfig, typename TPos>
SEQAN_HOST_DEVICE inline bool isSentinel(LF<StringSet<TText, TSSetSpec>, TSpec, TConfig> const & lf, TPos pos)
{
    return getValue(lf.sentinels, pos);
}

// ----------------------------------------------------------------------------
// Function _getPrefixSum()
// ----------------------------------------------------------------------------

template <typename TText, typename TSpec, typename TConfig, typename TValue>
SEQAN_HOST_DEVICE inline
typename Size<LF<TText, TSpec, TConfig> const>::Type
_getPrefixSum(LF<TText, TSpec, TConfig> const & lf, TValue val)
{
    typedef LF<TText, TSpec, TConfig> const             TLF;
    typedef typename Value<TLF>::Type                   TTextValue;

    return getValue(lf.sums, ordValue(static_cast<TTextValue>(val)));
}

// ----------------------------------------------------------------------------
// Function _getBwtRank(pos, val)
// ----------------------------------------------------------------------------

template <typename TText, typename TSpec, typename TConfig, typename TPos, typename TValue>
SEQAN_HOST_DEVICE inline
typename Size<LF<TText, TSpec, TConfig> >::Type
_getBwtRank(LF<TText, TSpec, TConfig> const & lf, TPos pos, TValue val)
{
    typedef LF<TText, TSpec, TConfig> const                TLF;
    typedef typename Size<TLF>::Type                       TSize;

    TSize rank = _getPrefixSum(lf, val);

    if (pos > 0)
    {
        rank += getRank(lf.bwt, pos - 1, val);

        if (ordEqual(lf.sentinelSubstitute, val))
            rank -= _getSentinelsRank(lf, pos - 1);
    }

    return rank;
}

// ----------------------------------------------------------------------------
// Function _getBwtRank(pos)
// ----------------------------------------------------------------------------

template <typename TText, typename TSpec, typename TConfig, typename TPos>
SEQAN_HOST_DEVICE inline
typename Size<LF<TText, TSpec, TConfig> const>::Type
_getBwtRank(LF<TText, TSpec, TConfig> const & lf, TPos pos)
{
    return _getBwtRank(lf, pos, getValue(lf.bwt, pos));
}

// ----------------------------------------------------------------------------
// Function _setSentinelSubstitute()
// ----------------------------------------------------------------------------
// This function determines the '$' substitute.
// The character with the smallest number of occurrences greater 0 is chosen.

template <typename TText, typename TSpec, typename TConfig>
inline void _setSentinelSubstitute(LF<TText, TSpec, TConfig> & lf)
{
    typedef LF<TText, TSpec, TConfig>                   TLF;
    typedef typename Fibre<TLF, FibrePrefixSums>::Type  TPrefixSums;
    typedef typename Value<TPrefixSums>::Type           TValue;
    typedef typename Size<TPrefixSums>::Type            TSize;

    TValue minOcc = MaxValue<TValue>::VALUE;
    TSize ordVal = 0;

    for (TSize i = 0; i < length(lf.sums) - 1; ++i)
    {
        TSize occ = getValue(lf.sums, i + 1) - getValue(lf.sums, i);
        if (occ > 0 && occ < minOcc)
        {
            minOcc = occ;
            ordVal = i;
        }
    }

    lf.sentinelSubstitute = ordVal;
}

// ----------------------------------------------------------------------------
// Function _createBwt()
// ----------------------------------------------------------------------------
// This function computes the BWT of a text. Note that the sentinel sign is substituted and its position stored.

template <typename TText, typename TSpec, typename TConfig, typename TBwt, typename TOtherText, typename TSA>
inline void
_createBwt(LF<TText, TSpec, TConfig> & lf, TBwt & bwt, TOtherText const & text, TSA const & sa)
{
    typedef typename GetValue<TSA>::Type                    TSAValue;
    typedef typename Iterator<TSA const, Standard>::Type    TSAIter;
    typedef typename Iterator<TBwt, Standard>::Type         TBwtIter;

    TSAIter saIt = begin(sa, Standard());
    TSAIter saItEnd = end(sa, Standard());
    TBwtIter bwtIt = begin(bwt, Standard());

    assignValue(bwtIt, back(text));
    ++bwtIt;

    for (; saIt != saItEnd; ++saIt, ++bwtIt)
    {
        TSAValue pos = getValue(saIt);

        if (pos != 0)
        {
            assignValue(bwtIt, getValue(text, pos - 1));
        }
        else
        {
            assignValue(bwtIt, lf.sentinelSubstitute);
            lf.sentinels = bwtIt - begin(bwt, Standard());
        }
    }
}

// ----------------------------------------------------------------------------
// Function _createBwt()
// ----------------------------------------------------------------------------
// This function computes the BWT of a text collection. Note that the sentinel sign is substituted and its position stored.

template <typename TText, typename TSSetSpec, typename TSpec, typename TConfig, typename TBwt, typename TOtherText, typename TSA>
inline void
_createBwt(LF<StringSet<TText, TSSetSpec>, TSpec, TConfig> & lf, TBwt & bwt, TOtherText const & text, TSA const & sa)
{
    typedef typename Value<TSA>::Type                       TSAValue;
    typedef typename Size<TSA>::Type                        TSize;
    typedef typename Iterator<TSA const, Standard>::Type    TSAIter;
    typedef typename Iterator<TBwt, Standard>::Type         TBwtIter;

    TSize seqNum = countSequences(text);
    TSize totalLen = lengthSum(text);

    resize(lf.sentinels, seqNum + totalLen, Exact());

    TSAIter saIt = begin(sa, Standard());
    TSAIter saItEnd = end(sa, Standard());
    TBwtIter bwtItBeg = begin(bwt, Standard());
    TBwtIter bwtIt = bwtItBeg;

    // Fill the sentinel positions (they are all at the beginning of the bwt).
    for (TSize i = 1; i <= seqNum; ++i, ++bwtIt)
    {
        assignValue(bwtIt, back(text[seqNum - i]));
        setValue(lf.sentinels, bwtIt - bwtItBeg, false);
    }

    // Compute the rest of the bwt.
    for (; saIt != saItEnd; ++saIt, ++bwtIt)
    {
        TSAValue pos;    // = SA[i];
        posLocalize(pos, getValue(saIt), stringSetLimits(text));

        if (getSeqOffset(pos) != 0)
        {
            assignValue(bwtIt, getValue(getValue(text, getSeqNo(pos)), getSeqOffset(pos) - 1));
            setValue(lf.sentinels, bwtIt - bwtItBeg, false);
        }
        else
        {
            assignValue(bwtIt, lf.sentinelSubstitute);
            setValue(lf.sentinels, bwtIt - bwtItBeg, true);
        }
    }

    // Update the auxiliary RankDictionary of sentinel positions.
    updateRanks(lf.sentinels);
}

// ----------------------------------------------------------------------------
// Function createLF()
// ----------------------------------------------------------------------------
/*!
 * @fn LF#createLF
 *
 * @headerfile <seqan/index.h>
 *
 * @brief Creates the LF table
 *
 * @signature void createLF(lfTable, text, sa);
 *
 * @param[out] lfTable The LF table to be constructed.
 * @param[in]  text    The underlying text Types: @link String @endlink.
 * @param[in]  sa      The suffix array of the LF table underlying text. Types: @link String @endlink,
 *                     @link StringSet @endlink.
 *
 * @return TReturn Returns a <tt>bool</tt> which is <tt>true</tt> on successes and <tt>false</tt> otherwise.
 */
// This function creates all table of the lf table given a text and a suffix array.
template <typename TText, typename TSpec, typename TConfig, typename TOtherText, typename TSA>
inline void createLF(LF<TText, TSpec, TConfig> & lf, TOtherText const & text, TSA const & sa)
{
    typedef LF<TText, TSpec, TConfig>                          TLF;
    typedef typename Fibre<TLF, FibreTempBwt>::Type            TBwt;
    typedef typename Value<TLF>::Type                          TValue;
    typedef typename Size<TLF>::Type                           TSize;

    // Clear assuming undefined state.
    clear(lf);

    // Compute prefix sum.
    prefixSums<TValue>(lf.sums, text);

    // Choose the sentinel substitute.
    _setSentinelSubstitute(lf);

    // Create BWT and mark sentinels.
    TBwt bwt;
    resize(bwt, bwtLength(text), Exact());
    _createBwt(lf, bwt, text, sa);

    // Index BWT bwt for rank queries.
    createRankDictionary(lf.bwt, bwt);

    // Add sentinels to prefix sum.
    TSize sentinelsCount = countSequences(text);
    for (TSize i = 0; i < length(lf.sums); ++i)
        lf.sums[i] += sentinelsCount;
}

// ----------------------------------------------------------------------------
// Function open()
// ----------------------------------------------------------------------------

/*!
 * @fn LF#open
 * @headerfile <seqan/index.h>
 * @brief This functions loads a LF table from disk.
 *
 * @signature bool open(lfTable, fileName[, openMode]);
 *
 * @param[in,out] lfTable  The LF object.
 * @param[in]     fileName C-style character string containing the file name.
 * @param[in]      openMode
 *                     The combination of flags defining how the file should be opened.  To open a file
 *                     read-only, write-only or to read and write use <tt>OPEN_RDONLY</tt>, <tt>OPEN_WRONLY</tt>,
 *                     or <tt>OPEN_RDWR</tt>.  To create or overwrite a file add <tt>OPEN_CREATE</tt>.  To append
 *                     a file if existing add <tt>OPEN_APPEND</tt>.  To circumvent problems, files are always
 *                     opened in binary mode.  Default: <tt>OPEN_RDWR | OPEN_CREATE | OPEN_APPEND</tt>.
 *
 * @return bool <tt>true</tt> on success.
 */


template <typename TText, typename TSpec, typename TConfig>
inline bool open(LF<TText, TSpec, TConfig> & lf, const char * fileName, int openMode)
{
    String<char> name;

    name = fileName;
    append(name, ".pst");
    if (!open(lf.sums, toCString(name), openMode)) return false;

    name = fileName;
    append(name, ".drv");
    if (!open(lf.bwt, toCString(name), openMode)) return false;

    name = fileName;
    append(name, ".drp");
    if (!open(lf.sentinels, toCString(name), openMode)) return false;

    name = fileName;
    append(name, ".drs");
    if (!open(lf.sentinelSubstitute, toCString(name), openMode)) return false;

    return true;
}

template <typename TText, typename TSpec, typename TConfig>
inline bool open(LF<TText, TSpec, TConfig> & lf, const char * fileName)
{
    return open(lf, fileName, DefaultOpenMode<LF<TText, TSpec, TConfig> >::VALUE);
}

// ----------------------------------------------------------------------------
// Function save()
// ----------------------------------------------------------------------------

/*!
 * @fn LF#save
 *
 * @headerfile <seqan/index.h>
 *
 * @brief This functions saves a LF table to disk.
 *
 * @signature bool save(lfTable, fileName[, openMode]);
 *
 * @param[in] lfTable  The LF object to save.
 * @param[in] fileName C-style character string containing the file name.
 * @param[in] openMode
 *                     The combination of flags defining how the file should be opened.  To open a file
 *                     read-only, write-only or to read and write use <tt>OPEN_RDONLY</tt>, <tt>OPEN_WRONLY</tt>,
 *                     or <tt>OPEN_RDWR</tt>.  To create or overwrite a file add <tt>OPEN_CREATE</tt>.  To append
 *                     a file if existing add <tt>OPEN_APPEND</tt>.  To circumvent problems, files are always
 *                     opened in binary mode.  Default: <tt>OPEN_RDWR | OPEN_CREATE | OPEN_APPEND</tt>.
 *
 * @return bool <tt>true</tt> on success.
 */

template <typename TText, typename TSpec, typename TConfig>
inline bool save(LF<TText, TSpec, TConfig> const & lf, const char * fileName, int openMode)
{
    String<char> name;

    name = fileName;
    append(name, ".pst");
    if (!save(lf.sums, toCString(name), openMode)) return false;

    name = fileName;
    append(name, ".drv");
    if (!save(lf.bwt, toCString(name), openMode)) return false;

    name = fileName;
    append(name, ".drp");
    if (!save(lf.sentinels, toCString(name), openMode)) return false;

    name = fileName;
    append(name, ".drs");
    if (!save(lf.sentinelSubstitute, toCString(name), openMode)) return false;

    return true;
}

template <typename TText, typename TSpec, typename TConfig>
inline bool save(LF<TText, TSpec, TConfig> const & lf, const char * fileName)
{
    return save(lf, fileName, DefaultOpenMode<LF<TText, TSpec, TConfig> >::VALUE);
}

}
#endif // INDEX_FM_LF_TABLE_H_
