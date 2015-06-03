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

#ifndef INDEX_FM_RANKDICTIONARY_WT
#define INDEX_FM_RANKDICTIONARY_WT

namespace seqan {

// ==========================================================================
// Tags
// ==========================================================================

// --------------------------------------------------------------------------
// Tag WaveletTreeConfig
// --------------------------------------------------------------------------

template <typename TSize = size_t, typename TFibre = Alloc<>, unsigned LEVELS = 1, unsigned ARITY_ = 2>
struct WTRDConfig : LevelsRDConfig<TSize, TFibre, LEVELS>
{
    static const unsigned ARITY = ARITY_;
};

// --------------------------------------------------------------------------
// Tag WaveletTree
// --------------------------------------------------------------------------

template <typename TSpec = void, typename TConfig = WTRDConfig<> >
struct WaveletTree {};

struct FibreTreeStructure_;
typedef Tag<FibreTreeStructure_>    const FibreTreeStructure;

// ==========================================================================
// Metafunctions
// ==========================================================================
/*!
 * @defgroup WaveletTreeFibres WaveletTree Fibres
 * @brief Tag to select a specific fibre (e.g. table, object, ...) of a @link WaveletTree @endlink.
 *
 * These tags can be used to get @link Fibre Fibres @endlink of a @link WaveletTree @endlink.
 *
 * @see Fibre
 * @see Index#getFibre
 *
 * @tag WaveletTreeFibres#FibreTreeStructure
 * @brief The wavelet tree structure of the wavelet tree.
 *
 * @tag WaveletTreeFibres#FibreRanks
 * @brief A string set containing a rank support bit string for each node in the tree.
 *
 */

// ----------------------------------------------------------------------------
// Metafunction Fibre
// ----------------------------------------------------------------------------

/*!
 * @class WaveletTree
 * @extends RankDictionary
 * @headerfile <seqan/index.h>
 *
 * @brief A WaveletTree is a hierarchical @link RankDictionary @endlink.
 *
 * @signature template <typename TValue, typename TSpec, typename TConfig>
 *            class RankDictionary<TValue, WaveletTree<TSpec, TConfig> >;
 *
 * @tparam TValue The alphabet type of the wavelet tree.
 * @tparam TSpec A tag for specialization purposes. Default: <tt>void</tt>
 *
 * The nodes of a wavelet tree consist of a bit string as well as a character c.
 * In each level of the tree, characters smaller than c are represented as a 0
 * while character greater or equal to c are represented with a 1. The
 * characters represented by a 0 form the string to be represented by the left
 * subtree while characters represented by a 1 form the string of the right
 * subtree. Therefore, only the bit string of the root node represents all
 * characters while all other nodes represent subsets.
 */
// TODO(esiragusa): update doc


template <typename TValue, typename TSpec, typename TConfig>
struct Fibre<RankDictionary<TValue, WaveletTree<TSpec, TConfig> >, FibreRanks>
{
    typedef String<RankDictionary<bool, Levels<TSpec, TConfig> > > Type;
};

template <typename TValue, typename TSpec, typename TConfig>
struct Fibre<RankDictionary<TValue, WaveletTree<TSpec, TConfig> >, FibreTreeStructure>
{
    typedef typename MakeUnsigned<TValue>::Type TUChar_;
    typedef RightArrayBinaryTree<TUChar_, void>  Type;
};

// ----------------------------------------------------------------------------
// Metafunction Value
// ----------------------------------------------------------------------------
// TODO(esiragusa): remove this - const version should be Value const (as by default)

template <typename TValue, typename TSpec, typename TConfig>
struct Value<RankDictionary<TValue, WaveletTree<TSpec, TConfig> > >
{
    typedef TValue Type;
};

template <typename TValue, typename TSpec, typename TConfig>
struct Value<RankDictionary<TValue, WaveletTree<TSpec, TConfig> > const> :
    public Value<RankDictionary<TValue, WaveletTree<TSpec, TConfig> > > {};

// ----------------------------------------------------------------------------
// Metafunction RankDictionaryBlock_
// ----------------------------------------------------------------------------

template <typename TValue, typename TSpec, typename TConfig>
struct RankDictionaryBlock_<TValue, WaveletTree<TSpec, TConfig> >
{
    typedef RankDictionary<TValue, WaveletTree<TSpec, TConfig> >    TRankDictionary_;
    typedef typename Size<TRankDictionary_>::Type                   TSize_;
//    typedef Tuple<TSize_, ValueSize<TValue>::VALUE>                 Type;
    typedef String<TSize_>                                          Type;
};

// ==========================================================================
// Classes
// ==========================================================================

// ----------------------------------------------------------------------------
// Spec WaveletTree
// ----------------------------------------------------------------------------

template <typename TValue, typename TSpec, typename TConfig>
struct RankDictionary<TValue, WaveletTree<TSpec, TConfig> >
{
    typename Fibre<RankDictionary, FibreRanks>::Type            ranks;
    typename Fibre<RankDictionary, FibreTreeStructure>::Type    waveletTreeStructure;

    RankDictionary() {}

    template <typename TText>
    RankDictionary(TText const & text)
    {
        createRankDictionary(*this, text);
    }
};

// ==========================================================================
// Functions
// ==========================================================================

// ----------------------------------------------------------------------------
// Function getFibre()
// ----------------------------------------------------------------------------

template <typename TValue, typename TSpec, typename TConfig>
inline typename Fibre<RankDictionary<TValue, WaveletTree<TSpec, TConfig> >, FibreTreeStructure>::Type &
getFibre(RankDictionary<TValue, WaveletTree<TSpec, TConfig> > & dict, FibreTreeStructure)
{
    return dict.waveletTreeStructure;
}

template <typename TValue, typename TSpec, typename TConfig>
inline typename Fibre<RankDictionary<TValue, WaveletTree<TSpec, TConfig> >, FibreTreeStructure>::Type const &
getFibre(RankDictionary<TValue, WaveletTree<TSpec, TConfig> > const & dict, FibreTreeStructure)
{
    return dict.waveletTreeStructure;
}

// ----------------------------------------------------------------------------
// Function clear()
// ----------------------------------------------------------------------------

template <typename TValue, typename TSpec, typename TConfig>
inline void clear(RankDictionary<TValue, WaveletTree<TSpec, TConfig> > & dict)
{
    clear(getFibre(dict, FibreRanks()));
    clear(getFibre(dict, FibreTreeStructure()));
}

// ----------------------------------------------------------------------------
// Function empty()
// ----------------------------------------------------------------------------

template <typename TValue, typename TSpec, typename TConfig>
inline bool empty(RankDictionary<TValue, WaveletTree<TSpec, TConfig> > const & dict)
{
    return empty(getFibre(dict, FibreRanks())) && empty(getFibre(dict, FibreTreeStructure()));
}

// ----------------------------------------------------------------------------
// Function getValue()
// ----------------------------------------------------------------------------

template <typename TValue, typename TSpec, typename TConfig, typename TPos>
inline TValue getValue(RankDictionary<TValue, WaveletTree<TSpec, TConfig> > & dict, TPos pos)
{
    typedef typename Fibre<RankDictionary<TValue, WaveletTree<TSpec, TConfig> >, FibreTreeStructure>::Type const    TWaveletTreeStructure;
    typedef typename Fibre<TWaveletTreeStructure, FibreTreeStructureEncoding>::Type                 TWaveletTreeStructureString;
    typedef typename Value<TWaveletTreeStructureString>::Type                                       TWaveletTreeStructureEntry;
    typedef typename Value<TWaveletTreeStructureEntry, 1>::Type                                     TChar;

    unsigned treePos = 0;
    typename Iterator<TWaveletTreeStructure, TopDown<> >::Type iter(dict.waveletTreeStructure, treePos);

    // initialize the return value with the smallest possible value
    TChar character = dict.waveletTreeStructure.minCharValue;

    // while the current node is not a leaf, go right if the bit at the current position is 1
    // go left otherwise
    // note that when going right the return value changes
    while (true)
    {
        TPos rank1 = getRank(dict.ranks[treePos], pos);
        if (getValue(dict.ranks[treePos], pos))
        {
            character = getCharacter(iter);
            pos = rank1 - 1;  // -1 because strings start at 0
            if (!goRightChild(iter))
                break;
        }
        else
        {
            pos -= rank1;
            if (!goLeftChild(iter))
                break;
        }
        treePos = getPosition(iter);
    }

    return character;
}

template <typename TValue, typename TSpec, typename TConfig, typename TPos>
inline TValue getValue(RankDictionary<TValue, WaveletTree<TSpec, TConfig> > const & dict, TPos pos)
{
    return getValue(const_cast<RankDictionary<TValue, WaveletTree<TSpec, TConfig> > &>(dict), pos);
}

// ----------------------------------------------------------------------------
// Function getRank()
// ----------------------------------------------------------------------------

template <typename TValue, typename TSpec, typename TConfig, typename TPos, typename TChar>
inline typename Size<RankDictionary<TValue, WaveletTree<TSpec, TConfig> > >::Type
getRank(RankDictionary<TValue, WaveletTree<TSpec, TConfig> > const & dict, TPos pos, TChar character)
{
    typedef typename Fibre<RankDictionary<TValue, WaveletTree<TSpec, TConfig> >, FibreTreeStructure>::Type  TWaveletTreeStructure;
    typedef typename Fibre<TWaveletTreeStructure, FibreTreeStructureEncoding>::Type         TWaveletTreeStructureString;
    typedef typename Value<TWaveletTreeStructureString>::Type                               TWaveletTreeStructureEntry;
    typedef typename Value<TWaveletTreeStructureEntry, 1>::Type                             TChar_;

    TPos sum = pos;
    TPos treePos = 0;

    // determine the leaf containing the character
    // count the number of 1 or 0 up to the computed position
    typename Iterator<TWaveletTreeStructure const, TopDown<> >::Type it(dict.waveletTreeStructure, treePos);
    TChar_ charInTree = dict.waveletTreeStructure.minCharValue;

    while (true)
    {
        TPos addValue = getRank(dict.ranks[treePos], sum);
        if (ordGreater(getCharacter(it), character))
        {
            if (addValue > sum) return 0;

            sum -= addValue;
            if (!goLeftChild(it))
                break;
        }
        else
        {
            if (addValue == 0) return 0;

            charInTree = getCharacter(it);
            sum = addValue - 1;
            if (!goRightChild(it))
                break;
        }
        treePos = getPosition(it);
    }

    if (ordEqual(charInTree, character))
        return sum + 1;

    return 0;
}

// ----------------------------------------------------------------------------
// Function _fillStructure()
// ----------------------------------------------------------------------------

// This function is used to fill the bit strings of the wavelet tree.
template <typename TValue, typename TSpec, typename TConfig, typename TText>
inline void _fillStructure(RankDictionary<TValue, WaveletTree<TSpec, TConfig> > & dict, TText const & text)
{
    typedef RankDictionary<TValue, WaveletTree<TSpec, TConfig> >        TRankDictionary;
    typedef typename Fibre<TRankDictionary, FibreTreeStructure>::Type   TWaveletTreeStructure;
    typedef typename Iterator<TWaveletTreeStructure, TopDown<> >::Type  TWaveletTreeIterator;
    typedef typename Size<TRankDictionary>::Type                        TSize;
    typedef typename Iterator<TText const, Standard>::Type              TTextIterator;

    resize(dict.ranks, length(dict.waveletTreeStructure), Exact());

    for (TSize i = 0; i < length(dict.ranks); ++i)
        resize(dict.ranks[i], 0);

    TTextIterator textBegin = begin(text, Standard());
    TTextIterator textEnd = end(text, Standard());

    for (TTextIterator textIt = textBegin; textIt != textEnd; ++textIt)
    {
        TWaveletTreeIterator it(dict.waveletTreeStructure, 0);

        while (true)
        {
            // decide whether the character is smaller then the pivot element of the current node
            if (ordGreater(getCharacter(it), value(textIt)))
            {
                // TODO(esiragusa): use resize() & setValue() instead of appendValue().
                appendValue(dict.ranks[getPosition(it)], false);
                if (!goLeftChild(it))
                    break;
            }
            else
            {
                appendValue(dict.ranks[getPosition(it)], true);
                if (!goRightChild(it))
                    break;
            }
        }
    }
}

// ----------------------------------------------------------------------------
// Function updateRanks()                                      [RankDictionary]
// ----------------------------------------------------------------------------

template <typename TValue, typename TSpec, typename TConfig>
inline void updateRanks(RankDictionary<TValue, WaveletTree<TSpec, TConfig> > & dict)
{
    typedef RankDictionary<TValue, WaveletTree<TSpec, TConfig> >        TRankDictionary;
    typedef typename Size<TRankDictionary>::Type                        TSize;

    for (TSize i = 0; i < length(getFibre(dict, FibreRanks())); ++i)
        updateRanks(getFibre(dict, FibreRanks())[i]);
}

// ----------------------------------------------------------------------------
// Function createRankDictionary()
// ----------------------------------------------------------------------------
template <typename TValue, typename TSpec, typename TConfig, typename TText, typename TPrefixSums>
inline void
createRankDictionary(RankDictionary<TValue, WaveletTree<TSpec, TConfig> > & dict, TText const & text, TPrefixSums const & sums)
{
    createRightArrayBinaryTree(getFibre(dict, FibreTreeStructure()), sums);
//    _resizeStructure(dict, text);
    _fillStructure(dict, text);
    updateRanks(dict);
}

template <typename TValue, typename TSpec, typename TConfig, typename TText>
inline void
createRankDictionary(RankDictionary<TValue, WaveletTree<TSpec, TConfig> > & dict, TText const & text)
{
    typename RankDictionaryBlock_<TValue, WaveletTree<TSpec, TConfig> >::Type sums;
    prefixSums<TValue>(sums, text);
    createRankDictionary(dict, text, sums);
}

// ----------------------------------------------------------------------------
// Function open()
// ----------------------------------------------------------------------------

template <typename TValue, typename TSpec, typename TConfig>
inline bool open(RankDictionary<TValue, WaveletTree<TSpec, TConfig> > & dict, const char * fileName, int openMode)
{
    String<char> name;

    name = fileName;
    append(name, ".wtc");
    if (!open(getFibre(dict, FibreRanks()), toCString(name), openMode)) return false;

    name = fileName;
    append(name, ".wts");
    if (!open(getFibre(dict, FibreTreeStructure()), toCString(name), openMode)) return false;

    return true;
}

// ----------------------------------------------------------------------------
// Function save()
// ----------------------------------------------------------------------------

template <typename TValue, typename TSpec, typename TConfig>
inline bool save(RankDictionary<TValue, WaveletTree<TSpec, TConfig> > const & dict, const char * fileName, int openMode)
{
    String<char> name;

    name = fileName;
    append(name, ".wtc");
    if (!save(getFibre(dict, FibreRanks()), toCString(name), openMode)) return false;

    name = fileName;
    append(name, ".wts");
    if (!save(getFibre(dict, FibreTreeStructure()), toCString(name), openMode)) return false;

    return true;
}

}
#endif  // INDEX_FM_RANKDICTIONARY_WT
