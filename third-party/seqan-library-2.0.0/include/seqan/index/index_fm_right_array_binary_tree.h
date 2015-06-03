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

#ifndef INDEX_FM_RIGHT_ARRAY_BINARY_TREE_H
#define INDEX_FM_RIGHT_ARRAY_BINARY_TREE_H

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

template <typename TChar, typename TSpec = void>
class RightArrayBinaryTree;

// ============================================================================
// Tags
// ============================================================================
//
/*!
 * @defgroup RightArrayBinaryTreeFibres RightArrayBinaryTree Fibres
 * @brief Tag to select a specific fibre (e.g. table, object, ...) of a @link
 *        RightArrayBinaryTree @endlink.
 *
 * These tags can be used to get @link Fibre Fibres @endlink of a @link RightArrayBinaryTree @endlink.
 *
 * @see Fibre
 * @see Index#getFibre
 *
 * @tag RightArrayBinaryTreeFibres#FibreTreeStructureEncoding
 *
 * @brief The string encoding the wavelet tree structure.
 */

struct FibreTreeStructureEncoding_;
typedef Tag<FibreTreeStructureEncoding_> const FibreTreeStructureEncoding;

struct FibreTreeStructure_;
typedef Tag<FibreTreeStructure_> const FibreTreeStructure;

// ============================================================================
// Metafunctions
// ============================================================================

// ----------------------------------------------------------------------------
// Metafunction Fibre
// ----------------------------------------------------------------------------

template <typename TChar, typename TSpec>
struct Fibre<RightArrayBinaryTree<TChar, TSpec>, FibreTreeStructureEncoding>
{
    // TODO (singer): Description why we need + 2
    typedef typename BitVector_<Log2<ValueSize<TChar>::VALUE + 2>::VALUE>::Type TPos;
    typedef String<Pair<TChar, TPos> > Type;
};

// ----------------------------------------------------------------------------
// Metafunction Value
// ----------------------------------------------------------------------------

template <typename TChar, typename TSpec>
struct Value<RightArrayBinaryTree<TChar, TSpec> >
{
    typedef typename Fibre<RightArrayBinaryTree<TChar, TSpec>, FibreTreeStructureEncoding>::Type TWaveletTreeVertices_;
    typedef typename Value<TWaveletTreeVertices_>::Type TWaveletTreeVertex_;
    typedef typename Value<TWaveletTreeVertex_, 2>::Type TPos;

    typedef Pair<TChar, TPos> Type;
};

template <typename TChar, typename TSpec>
struct Value<RightArrayBinaryTree<TChar, TSpec> const> :
    Value<RightArrayBinaryTree<TChar, TSpec> > {};

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// ----------------------------------------------------------------------------
// Class RightArrayBinaryTree
// ----------------------------------------------------------------------------
/*!
 * @class RightArrayBinaryTree
 * @headerfile <seqan/index.h>
 * @brief A special format to encode the structure of a wavelet tree.  The structure is very space efficient because
 *        only one position is stored which encodes where the left and right subtree of a given node exist.
 *
 * @signature template <typename TChar, typename TSpec>
 *            class RightArrayBinaryTree;
 *
 * @tparam TValue The type of the stored characters.
 * @tparam TSpec  The wavelet tree structure specialisation. Default: void.
 */

template <typename TChar, typename TSpec>
class RightArrayBinaryTree
{
public:
    typename Fibre<RightArrayBinaryTree, FibreTreeStructureEncoding>::Type treeVertices;
    TChar minCharValue;

    RightArrayBinaryTree() :
        treeVertices(),
        minCharValue()
    {}

    template <typename TText>
    explicit RightArrayBinaryTree(TText const & text) :
        treeVertices(),
        minCharValue()
    {
        createRightArrayBinaryTree(*this, text);
    }
};

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function clear()
// ----------------------------------------------------------------------------
/*!
 * @fn RightArrayBinaryTree#clear
 * @headerfile <seqan/index.h>
 * @brief Resets a right-array-binary tree.
 *
 * @signature void clear(rightArrayBinaryTree);
 *
 * @param[in,out] rightArrayBinaryTree The RightArrayBinaryTree to be cleared.
 */

template <typename TChar, typename TSpec>
inline void clear(RightArrayBinaryTree<TChar, TSpec> & treeStructure)
{
    clear(treeStructure.treeVertices);
}

// ----------------------------------------------------------------------------
// Function createRightArrayBinaryTree()
// ----------------------------------------------------------------------------
/*!
 * @fn RightArrayBinaryTree#createRightArrayBinaryTree
 * @headerfile <seqan/index.h>
 * @brief Computes the right-array-binary tree of a text.
 *
 * @signature void createRightArrayBinaryTree(rightArrayBinaryTree, text);
 *
 * @param[in] rightArrayBinaryTree A wavelet tree structure.
 * @param[in] text                 A @link TextConcept text @endlink.
 */

// This function computes the wavelet tree structure.
template <typename TChar, typename TSpec, typename TIterSpec, typename TBorderString, typename TPrefixSums>
inline void _createRightArrayBinaryTreeImpl(Iter<RightArrayBinaryTree<TChar, TSpec>, TIterSpec> & it,
                                            TBorderString & borderString,
                                            TPrefixSums const & sums)
{
    do
    {
        // TODO (singer): Comment this
        if (back(borderString).i2 - back(borderString).i1 + 1 < 3 ||
            sums[back(borderString).i1] == sums[back(borderString).i2 + 1])
        {
            setCharacter(it, back(borderString).i1 + 1);
            SEQAN_ASSERT_MSG(isLeaf(it), "You just deleted a subtree.");
        }
        else
            _setChildVertices(it, borderString, sums);

        if (!_goDownConstruction(it) && !_setAndGoRight(it, borderString))
            while (_goUpStructureConstruction(it, borderString) && !_setAndGoRight(it, borderString)) ;
    }
    while (!isRoot(it));
}

// This function computes the wavelet tree structure.
template <typename TChar, typename TSpec, typename TIterSpec, typename TPrefixSums>
inline void _createRightArrayBinaryTreeImpl(Iter<RightArrayBinaryTree<TChar, TSpec>, TIterSpec> & it,
                                            TPrefixSums const & sums)
{
    typedef RightArrayBinaryTree<TChar, TSpec>                       TRightArrayBinaryTree;
    typedef typename Size<RightArrayBinaryTree<TChar, TSpec> >::Type TSize;

    TRightArrayBinaryTree & waveletTreeStructure = container(it);

    String<Pair<TSize> > borderString;
    appendValue(borderString, Pair<TSize>(0, ValueSize<TChar>::VALUE - 1));
    _resize(waveletTreeStructure, 1, Exact());
    _createRightArrayBinaryTreeImpl(it, borderString, sums);
}

// ----------------------------------------------------------------------------
// Function createRightArrayBinaryTree()
// ----------------------------------------------------------------------------

template <typename TChar, typename TSpec, typename TPrefixSums>
inline void
createRightArrayBinaryTree(RightArrayBinaryTree<TChar, TSpec> & waveletTreeStructure, TPrefixSums const & sums)
{
    typename Iterator<RightArrayBinaryTree<TChar, TSpec>, TopDown<ParentLinks<> > >::Type it(waveletTreeStructure, 0u);
    _createRightArrayBinaryTreeImpl(it, sums);
}

// ----------------------------------------------------------------------------
// Function empty()
// ----------------------------------------------------------------------------

/*!
 * @fn RightArrayBinaryTree#empty
 * @headerfile <seqan/index.h>
 * @brief Checks whether or not a right-array-binary tree contains any elements.
 *
 * @signature bool empty(rightArrayBinaryTree);
 *
 * @param[in] rightArrayBinaryTree The right-array-binary tree to be queried.
 *
 * @return bool Returns <tt>true</tt> if the rank-support-bit string is empty and <tt>false</tt> otherwise.
 */

template <typename TChar, typename TSpec>
inline bool empty(RightArrayBinaryTree<TChar, TSpec> const & treeStructure)
{
    return empty(getFibre(treeStructure, FibreTreeStructureEncoding()));
}

// ----------------------------------------------------------------------------
// Function getFibre()
// ----------------------------------------------------------------------------
/*!
 * @fn RightArrayBinaryTree#getFibre
 * @headerfile <seqan/index.h>
 * @brief Returns a specific fibre of a right-array-binary tree.
 *
 * @signature TFibre getFibre(rightArrayBinaryTree, fibreTag);
 *
 * @param[in] rightArrayBinaryTree
 *                      The container holding the fibre.
 * @param[in] fibreTag  A tag that identifies the @link Fibre @endlink.  Types: @link RightArrayBinaryTreeFibres
 *                      RightArrayBinaryTree Fibres @endlink.
 *
 * @return TFibre A reference to the @link Fibre @endlink object.
 */
template <typename TChar, typename TSpec>
inline typename Fibre<RightArrayBinaryTree<TChar, TSpec>, FibreTreeStructureEncoding>::Type &
getFibre(RightArrayBinaryTree<TChar, TSpec> & treeStructure, FibreTreeStructureEncoding)
{
    return treeStructure.treeVertices;
}

template <typename TChar, typename TSpec>
inline typename Fibre<RightArrayBinaryTree<TChar, TSpec>, FibreTreeStructureEncoding>::Type const &
getFibre(RightArrayBinaryTree<TChar, TSpec> const & treeStructure, FibreTreeStructureEncoding)
{
    return treeStructure.treeVertices;
}

// ----------------------------------------------------------------------------
// Function length()
// ----------------------------------------------------------------------------

// This function returns the number of different entries in the wavelet tree structure.
/*!
 * @fn RightArrayBinaryTree#length
 * @brief Returns the number of nodes in the right-array-binary-tree.
 *
 * @signature TSize length(rightArrayBinaryTree);
 *
 * @param[in] tree The right-array-binary-tree to query for its size.
 *
 * @return TSize The number of nodes in the right-array-binary-tree.
 */
template <typename TChar, typename TSpec>
inline typename Size<RightArrayBinaryTree<TChar, TSpec> >::Type
length(RightArrayBinaryTree<TChar, TSpec> const & tree)
{
    return length(tree.treeVertices);
}

// ----------------------------------------------------------------------------
// Function _resize()
// ----------------------------------------------------------------------------

// This function resizes the string holding the nodes of the wavelet tree structure.
template <typename TChar, typename TSpec, typename TSize, typename TExpand>
inline typename Size<RightArrayBinaryTree<TChar, TSpec> >::Type
_resize(RightArrayBinaryTree<TChar, TSpec> & treeStructure, TSize size, Tag<TExpand> tag)
{
    return resize(treeStructure.treeVertices, size, tag);
}

// This function resizes the string holding the nodes of the wavelet tree structure.
template <typename TChar, typename TSpec, typename TSize, typename TExpand>
inline typename Size<RightArrayBinaryTree<TChar, TSpec> >::Type
_resize(RightArrayBinaryTree<TChar, TSpec> & treeStructure, TSize size,
        typename Value<typename Fibre<RightArrayBinaryTree<TChar, TSpec>, FibreTreeStructureEncoding>::Type>::Type value,
        Tag<TExpand> tag)
{
    return resize(treeStructure.treeVertices, size, value, tag);
}

// ----------------------------------------------------------------------------
// Function open()
// ----------------------------------------------------------------------------
/*!
 * @fn RightArrayBinaryTree#open
 *
 * @headerfile <seqan/index.h>
 *
 * @brief This functions loads a @link RightArrayBinaryTree @endlink from disk.
 *
 * @signature bool open(rightArrayBinaryTree, fileName [, openMode])
 *
 * @param[in,out] rightArrayBinaryTree
 *                               The RightArrayBinaryTree.
 * @param[in]     fileName       C-style character string containing the file name.
 * @param[in]     openMode       The combination of flags defining how the file should be
 *                               opened.  To open a file read-only, write-only or to read and
 *                               write use <tt>OPEN_RDONLY</tt>, <tt>OPEN_WRONLY</tt>, or
 *                               <tt>OPEN_RDWR</tt>.  To create or overwrite a file add
 *                               <tt>OPEN_CREATE</tt>.  To append a file if existing add
 *                               <tt>OPEN_APPEND</tt>.  To circumvent problems, files are always
 *                               opened in binary mode.  Default: <tt>OPEN_RDWR | OPEN_CREATE |
 *                               OPEN_APPEND</tt>.
 *
 * @return bool A <tt>bool</tt> which is <tt>true</tt> on success.
 */

template <typename TChar, typename TSpec>
inline bool open(RightArrayBinaryTree<TChar, TSpec> & treeStructure, const char * fileName, int openMode)
{
    String<char> name;

    name = fileName;
    append(name, ".rtv");
    if (!open(getFibre(treeStructure, FibreTreeStructureEncoding()), toCString(name), openMode)) return false;

    name = fileName;
    append(name, ".rtm");
    if (!open(treeStructure.minCharValue, toCString(name), openMode)) return false;

    return true;
}

template <typename TChar, typename TSpec>
inline bool open(RightArrayBinaryTree<TChar, TSpec> & treeStructure, const char * fileName)
{
    return open(treeStructure, fileName, DefaultOpenMode<RightArrayBinaryTree<TChar, TSpec> >::VALUE);
}

// ----------------------------------------------------------------------------
// Function save()
// ----------------------------------------------------------------------------
/*!
 * @fn RightArrayBinaryTree#save
 *
 * @headerfile <seqan/index.h>
 *
 * @brief This functions saves a @link RightArrayBinaryTree @endlink to disk.
 *
 * @signature bool save(rightArrayBinaryTree, fileName [, openMode])
 *
 * @param[in,out] rightArrayBinaryTree
 *                               The RightArrayBinaryTree.
 * @param[in]     fileName       C-style character string containing the file name.
 * @param[in]     openMode       The combination of flags defining how the file should be
 *                               opened.  To open a file read-only, write-only or to read and
 *                               write use <tt>OPEN_RDONLY</tt>, <tt>OPEN_WRONLY</tt>, or
 *                               <tt>OPEN_RDWR</tt>.  To create or overwrite a file add
 *                               <tt>OPEN_CREATE</tt>.  To append a file if existing add
 *                               <tt>OPEN_APPEND</tt>.  To circumvent problems, files are always
 *                               opened in binary mode.  Default: <tt>OPEN_RDWR | OPEN_CREATE |
 *                               OPEN_APPEND</tt>.
 *
 * @return bool A <tt>bool</tt> which is <tt>true</tt> on success.
 */

template <typename TChar, typename TSpec>
inline bool save(RightArrayBinaryTree<TChar, TSpec> const & treeStructure, const char * fileName, int openMode)
{
    String<char> name;

    name = fileName;
    append(name, ".rtv");
    if (!save(getFibre(treeStructure, FibreTreeStructureEncoding()), toCString(name), openMode)) return false;

    name = fileName;
    append(name, ".rtm");
    if (!save(treeStructure.minCharValue, toCString(name), openMode)) return false;

    return true;
}

template <typename TChar, typename TSpec>
inline bool save(RightArrayBinaryTree<TChar, TSpec> const & treeStructure, const char * fileName)
{
    return save(treeStructure, fileName, DefaultOpenMode<RightArrayBinaryTree<TChar, TSpec> >::VALUE);
}

}

#endif // INDEX_FM_RIGHT_ARRAY_BINARY_TREE_H
