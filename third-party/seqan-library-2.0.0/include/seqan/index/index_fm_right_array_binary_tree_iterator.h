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
//       from this software withoFIut specific prior written permission.
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

#ifndef INDEX_FM_RIGHT_ARRAY_BINARY_TREE_ITERATOR_H
#define INDEX_FM_RIGHT_ARRAY_BINARY_TREE_ITERATOR_H

//SEQAN_NO_DDDOC:do not generate documentation for this file

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

template <typename TSpec>
struct RightArrayBinaryTreeIterator;


// ============================================================================
// Metafunctions
// ============================================================================

// ----------------------------------------------------------------------------
// Metafunction Iterator
// ----------------------------------------------------------------------------

template <typename TChar, typename TSpec, typename TIterSpec>
struct Iterator<RightArrayBinaryTree<TChar, TSpec>, TopDown<TIterSpec> >
{
    typedef Iter<RightArrayBinaryTree<TChar, TSpec>, RightArrayBinaryTreeIterator<TopDown<TIterSpec> > > Type;
};

template <typename TChar, typename TSpec, typename TIterSpec>
struct Iterator<RightArrayBinaryTree<TChar, TSpec> const, TopDown<TIterSpec> >
{
    typedef Iter<RightArrayBinaryTree<TChar, TSpec> const, RightArrayBinaryTreeIterator<TopDown<TIterSpec> > > Type;
};

// ----------------------------------------------------------------------------
// Metafunction Spec
// ----------------------------------------------------------------------------

template <typename TChar, typename TSpec, typename TIterSpec>
struct Spec<Iter<RightArrayBinaryTree<TChar, TSpec>, RightArrayBinaryTreeIterator<TIterSpec> > >
{
    typedef TIterSpec Type;
};

template <typename TChar, typename TSpec, typename TIterSpec>
struct Spec<Iter<RightArrayBinaryTree<TChar, TSpec> const, RightArrayBinaryTreeIterator<TIterSpec> > >
{
    typedef TIterSpec Type;
};

// ============================================================================
// Classes
// ============================================================================
/*!
 * @class RightArrayBinaryTreeIterator RightArrayBinaryTree Iterator
 * @extends Iter
 * @headerfile <seqan/index.h>
 * @brief An iterator for @link RightArrayBinaryTree @endlink.
 *
 * @signature template <typename TSpec>
 *            class Iter<TRightArrayBinaryTree, TSpec >;
 *
 * @tparam TSpec                 Specialisation Tag. Types: TopDownIterator
 * @tparam TRightArrayBinaryTree The @link RightArrayBinaryTree @endlink.
 */
template <typename TTree, typename TIterSpec>
class Iter<TTree, RightArrayBinaryTreeIterator<TopDown<TIterSpec> > >
{
    typedef typename Fibre<TTree, FibreTreeStructureEncoding>::Type TWaveletTreeVertices;
    typedef typename Value<TWaveletTreeVertices>::Type TWaveletTreeVertex;
    typedef typename Value<TWaveletTreeVertex, 2>::Type TPos;

public:
    TPos position;
    TTree * waveletTreeStructure;

    Iter() :
        position(),
        waveletTreeStructure()
    {}

    Iter(TTree & treeStructure, TPos pos = 0) :
        position(pos),
        waveletTreeStructure(&treeStructure)
    {}
};

template <typename TTree, typename TSpec>
class Iter<TTree, RightArrayBinaryTreeIterator<TopDown<ParentLinks<TSpec> > > > :
    public Iter<TTree, RightArrayBinaryTreeIterator<TopDown<> > >
{
    typedef Iter<TTree, RightArrayBinaryTreeIterator<TopDown<> > > TBase;
    typedef typename Fibre<TTree, FibreTreeStructureEncoding>::Type TWaveletTreeVertices;
    typedef typename Value<TWaveletTreeVertices>::Type TWaveletTreeVertex;
    typedef typename Value<TWaveletTreeVertex, 2>::Type TPos;

public:
    String<TPos, Block<> > history;

    Iter() :
        TBase(),
        history()
    {}

    Iter(TTree & treeStructure, TPos pos = 0) :
        TBase(treeStructure,  pos),
        history()
    {
        appendValue(history, pos);
    }

};

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function begin()
// ----------------------------------------------------------------------------
/*!
 * @fn RightArrayBinaryTree#begin
 * @headerfile <seqan/index.h>
 * @brief The begin (root) of a @link RightArrayBinaryTree @endlink.
 *
 * @signature TIterator begin(rightArrayBinaryTree, iterSpec);
 *
 * @param[in] rightArrayBinaryTree The right-array-binary tree.
 * @param[in] iterSpec             A specialisation tag. Types: TopDown&lt;&gt;, TopDown&lt;ParentLinks&lt;&gt; &gt;.
 *
 * @return TIterator An iterator to the first item in <tt>object</tt>. Metafunctions: Metafunction.Iterator
 */
template <typename TChar, typename TSpec, typename TIterSpec>
inline typename Iterator<RightArrayBinaryTree<TChar, TSpec> const, TIterSpec>::Type
begin(RightArrayBinaryTree<TChar, TSpec> const & waveletTreeStructure, TIterSpec const &)
{
    return typename Iterator<RightArrayBinaryTree<TChar, TSpec> const, TIterSpec>::Type(waveletTreeStructure);
}

template <typename TChar, typename TSpec, typename TIterSpec>
inline typename Iterator<RightArrayBinaryTree<TChar, TSpec>, TIterSpec>::Type
begin(RightArrayBinaryTree<TChar, TSpec> & waveletTreeStructure, TIterSpec const &)
{
    return typename Iterator<RightArrayBinaryTree<TChar, TSpec>, TIterSpec>::Type(waveletTreeStructure);
}

// ----------------------------------------------------------------------------
// Function container()
// ----------------------------------------------------------------------------
/*!
 * @fn RightArrayBinaryTreeIterator#container
 * @headerfile <seqan/index.h>
 * @brief Container of an iterator.
 * @signature TContainer container(iterator)
 *
 * @param[in] iterator An iterator.
 *
 * @return TContainer The container that <tt>iterator</tt> traverses.
 */
template <typename TTree, typename TIterSpec>
inline TTree &
container(Iter<TTree, RightArrayBinaryTreeIterator<TopDown<TIterSpec> > > & it)
{
    return *it.waveletTreeStructure;
}

template <typename TTree, typename TIterSpec>
inline TTree &
container(Iter<TTree, RightArrayBinaryTreeIterator<TopDown<TIterSpec> > > const & it)
{
    return *it.waveletTreeStructure;
}

// ----------------------------------------------------------------------------
// Function end()
// ----------------------------------------------------------------------------
/*!
 * @fn RightArrayBinaryTree#end
 * @headerfile <seqan/index.h>
 * @brief The end (rigthmost leaf) of a @link RightArrayBinaryTree @endlink.
 *
 * @signature TIterator end(rightArrayBinaryTree, iterSpec);
 *
 * @param[in] rightArrayBinaryTree The right-array-binary tree.
 * @param[in] iterSpec             A specialisation tag. Types: TopDown&lt;&gt;, TopDown&lt;ParentLinks&lt;&gt; &gt;.
 *
 * @return TIterator An iterator to the first item in <tt>object</tt>.  Metafunctions: Metafunction.Iterator
 */
template <typename TChar, typename TSpec, typename TIterSpec>
inline typename Iterator<RightArrayBinaryTree<TChar, TSpec> const, TIterSpec>::Type
end(RightArrayBinaryTree<TChar, TSpec> const & waveletTreeStructure, TIterSpec const &)
{
    return typename Iterator<RightArrayBinaryTree<TChar, TSpec> const, TIterSpec>::Type(waveletTreeStructure, length(waveletTreeStructure));
}

template <typename TChar, typename TSpec, typename TIterSpec>
inline typename Iterator<RightArrayBinaryTree<TChar, TSpec>, TIterSpec>::Type
end(RightArrayBinaryTree<TChar, TSpec> & waveletTreeStructure, TIterSpec const &)
{
    return typename Iterator<RightArrayBinaryTree<TChar, TSpec>, TIterSpec>::Type(waveletTreeStructure, length(waveletTreeStructure));
}

// ----------------------------------------------------------------------------
// Function getCharacter()
// ----------------------------------------------------------------------------
/*!
 * @fn RightArrayBinaryTreeIterator#getCharacter
 * @headerfile <seqan/index.h>
 * @brief This function returns the pivot character of the node the iterator currently points to.
 *
 * @signature TChar getCharacter(iterator);
 *
 * @param[in] iterator The iterator.
 *
 * @return TChar The resulting character.
 */
template <typename TTree, typename TIterSpec>
inline typename Value<typename Value<TTree>::Type, 1>::Type
getCharacter(Iter<TTree, RightArrayBinaryTreeIterator<TopDown<TIterSpec> > > const & iter)
{
    return iter.waveletTreeStructure->treeVertices[getPosition(iter)].i1;
}

// ----------------------------------------------------------------------------
// Function getLeftChildPos()
// ----------------------------------------------------------------------------
/*!
 * @fn RightArrayBinaryTreeIterator#getLeftChildPos
 * @headerfile <seqan/index.h>
 * @brief Returns the position in @link RightArrayBinaryTree @endlink of the left child node.
 *
 * @signature unsigned getLeftChildPos(iterator);
 *
 * @param[in] iterator The iterator.
 *
 * @return unsigned The left child position.
 */
template <typename TTree, typename TIterSpec>
inline unsigned int getLeftChildPos(Iter<TTree, RightArrayBinaryTreeIterator<TopDown<TIterSpec> > > const & iter)
{
    if (iter.waveletTreeStructure->treeVertices[getPosition(iter)].i2 > 1)
    {
        return getPosition(iter) + 1;
    }
    return 0;
}

// ----------------------------------------------------------------------------
// Function getSubTreeSize()
// ----------------------------------------------------------------------------
/*!
 * @fn RightArrayBinaryTreeIterator#getSubTreeSize
 * @headerfile <seqan/index.h>
 * @brief Returns the number of vertices in the subtree starting at the position
 *        an iterator points to.
 *
 * @signature unsigned getSubTreeSize(iterator);
 *
 * @param[in] iterator The iterator.
 *
 * @return unsigned The subtree size.
 */
template <typename TTree, typename TIterSpec>
inline typename Size<TTree>::Type
getSubTreeSize(Iter<TTree, RightArrayBinaryTreeIterator<TopDown<TIterSpec> > > const & it)
{
    typedef typename Size<TTree>::Type TSize;

    Iter<TTree, RightArrayBinaryTreeIterator<TopDown<> > > _it(container(it));
    TSize originalPos = getPosition(it);
    goToPosition(_it, originalPos);
    while (goRightChild(_it) || goLeftChild(_it))
        continue;

    return getPosition(_it) - originalPos;
}

// ----------------------------------------------------------------------------
// Function getPosition()
// ----------------------------------------------------------------------------
/*!
 * @fn RightArrayBinaryTreeIterator#getPosition
 * @headerfile <seqan/index.h>
 * @brief Returns the position of the iterator in the host.
 *
 * @signature unsigned getPosition(iterator);
 *
 * @param[in] iterator The iterator.
 *
 * @return unsigned The position.
 */
template <typename TTree, typename TIterSpec>
inline typename Size<TTree>::Type
getPosition(Iter<TTree, RightArrayBinaryTreeIterator<TopDown<TIterSpec> > > const & it)
{
    return it.position;
}

// ----------------------------------------------------------------------------
// Function getRightChildPos()
// ----------------------------------------------------------------------------
/*!
 * @fn RightArrayBinaryTreeIterator#getRightChildPos
 * @headerfile <seqan/index.h>
 * @brief Returns the position in @link RightArrayBinaryTree @endlink of the right child node.
 *
 * @signature unsigned getLeftChildPos(iterator);
 *
 * @param[in] iterator The iterator.
 *
 * @return unsigned The left child position.
 */
template <typename TTree, typename TIterSpec>
inline typename Size<TTree>::Type
getRightChildPos(Iter<TTree, RightArrayBinaryTreeIterator<TopDown<TIterSpec> > > const & it)
{
    if (it.waveletTreeStructure->treeVertices[getPosition(it)].i2 > 2)
    {
        return it.waveletTreeStructure->treeVertices[getPosition(it)].i2 - 2;
    }
    if (it.waveletTreeStructure->treeVertices[getPosition(it)].i2 == 1)
    {
        return getPosition(it) + 1;
    }
    return 0;
}

// ----------------------------------------------------------------------------
// Function _historyPush()
// ----------------------------------------------------------------------------

template <typename TTree, typename TIterSpec, typename TPos>
inline void _historyPush(Iter<TTree, RightArrayBinaryTreeIterator<TopDown<TIterSpec> > > & /*tag*/ , TPos /*tag*/)
{}

template <typename TTree, typename TIterSpec, typename TPos>
inline void _historyPush(Iter<TTree, RightArrayBinaryTreeIterator<TopDown<ParentLinks<TIterSpec> > > > & it, TPos pos)
{
    appendValue(it.history, pos);
}

// ----------------------------------------------------------------------------
// Function goDown()
// ----------------------------------------------------------------------------
/*!
 * @fn RightArrayBinaryTreeIterator#goDown
 * @headerfile <seqan/index.h>
 * @brief Iterates down the leftmost edge in a @link RightArrayBinaryTree @endlink.
 *
 * @signature bool goDown(iterator);
 *
 * @param[in] iterator The iterator
 *
 * @return bool <tt>true</tt> if an edge to go down exists, otherwise <tt>false</tt>.
 */
template <typename TTree, typename TIterSpec>
inline bool goDown(Iter<TTree, RightArrayBinaryTreeIterator<TopDown<TIterSpec> > > & iter)
{
    if (goLeftChild(iter)) return true;
    if (goRightChild(iter)) return true;
    return false;
}

// ----------------------------------------------------------------------------
// Function _goDownConstruction()
// ----------------------------------------------------------------------------

template <typename TTree, typename TIterSpec>
inline bool _goDownConstruction(Iter<TTree, RightArrayBinaryTreeIterator<TopDown<TIterSpec> > > & it)
{
    if (goDown(it))
    {
        resize(container(it).treeVertices, length(container(it).treeVertices) + 1);
        return true;
    }
    return false;
}

// ----------------------------------------------------------------------------
// Function goLeftChild()
// ----------------------------------------------------------------------------

/*!
 * @fn RightArrayBinaryTreeIterator#goLeftChild
 * @headerfile <seqan/index.h>
 * @brief Sets the iterator to the left child of the current node if it exists and returns true, otherwise the
 *        iterator does not change position and the function returns false.
 *
 * @signature bool goLeftChild(iterator);
 *
 * @param[in] iterator The iterator
 *
 * @return bool <tt>true</tt> if the edge to go down exists, otherwise <tt>false</tt>.
 */
template <typename TTree, typename TIterSpec>
inline bool goLeftChild(Iter<TTree, RightArrayBinaryTreeIterator<TopDown<TIterSpec> > > & it)
{
    typedef typename Size<TTree>::Type TSize;

    TSize leftChildPos = getLeftChildPos(it);
    if (leftChildPos == 0)
        return false;

    if (!goToPosition(it, leftChildPos))
        return false;

    _historyPush(it, leftChildPos);
    return true;
}

// ----------------------------------------------------------------------------
// Function goRight()
// ----------------------------------------------------------------------------
/*!
 * @fn RightArrayBinaryTreeIterator#goRight
 * @headerfile <seqan/index.h>
 * @brief Iterates to the next sibling in a @link RightArrayBinaryTree @endlink.
 *
 * @signature bool goRight(iterator);
 *
 * @param[in] iterator The iterator
 *
 * @return bool <tt>true</tt> if the iterator could be moved, otherwise <tt>false</tt>.
 */
template <typename TTree, typename TIterSpec>
inline bool goRight(Iter<TTree, RightArrayBinaryTreeIterator<TopDown<TIterSpec> > > & it)
{
    typedef typename Size<TTree>::Type TSize;

    TSize pos = getPosition(it);
    if (goUp(it))
    {
        if (goRightChild(it))
        {
            if (pos != getPosition(it))
                return true;
        }
        else
            goLeftChild(it);
    }

    return false;
}

// ----------------------------------------------------------------------------
// Function goRightChild()
// ----------------------------------------------------------------------------
/*!
 * @fn RightArrayBinaryTreeIterator#goRightChild
 * @headerfile <seqan/index.h>
 * @brief Sets the iterator to the right child of the current node if it exists and returns true, otherwise the
 *        iterator does not change position and the function returns false.
 *
 * @signature bool goRightChild(iterator);
 *
 * @param[in] iterator The iterator.
 *
 * @return bool <tt>true</tt> if the edge to go down exists, otherwise <tt>false</tt>.
 */
template <typename TTree, typename TIterSpec>
inline bool goRightChild(Iter<TTree, RightArrayBinaryTreeIterator<TopDown<TIterSpec> > > & it)
{
    typedef typename Size<TTree>::Type TSize;

    TSize rightChildPos = getRightChildPos(it);
    if (rightChildPos == 0)
        return false;

    if (!goToPosition(it, rightChildPos))
        return false;

    _historyPush(it, rightChildPos);
    return true;
}

// template <typename TTree, typename TIterSpec>
// inline bool goRightChild(Iter<TTree, RightArrayBinaryTreeIterator<TopDown<ParentLinks<TIterSpec> > > > & iter)
// {
//     unsigned rightChildPos = getRightChildPos(iter);
//     if (rightChildPos == 0)
//         return false;
//
//     appendValue(iter.position, rightChildPos);
//     return true;
// }

// ----------------------------------------------------------------------------
// Function goToPosition()
// ----------------------------------------------------------------------------

// TODO(singer): Make this work!

template <typename TTree, typename TIterSpec, typename TPos>
inline bool goToPosition(Iter<TTree, RightArrayBinaryTreeIterator<TopDown<TIterSpec> > > & it, TPos pos)
{
    it.position = pos;
    return true;
}

// ----------------------------------------------------------------------------
// Function goUp()
// ----------------------------------------------------------------------------
/*!
 * @fn RightArrayBinaryTreeIterator#goUp
 * @headerfile <seqan/index.h>
 * @brief Iterates up one edge to the parent in a @link RightArrayBinaryTree @endlink.
 *
 * @signature bool goUp(iterator);
 *
 * @param[in] iterator The iterator.
 *
 * @return bool <tt>true</tt> if the iterator could be moved, otherwise <tt>false</tt>.
 */

template <typename TTree, typename TIterSpec>
inline bool goUp(Iter<TTree, RightArrayBinaryTreeIterator<TopDown<ParentLinks<TIterSpec> > > > & it)
{
    typedef typename Size<TTree>::Type TSize;

    TSize treeLevel = length(it.history);

    if (isRoot(it))
        return false;

    resize(it.history, treeLevel - 1);
    goToPosition(it, back(it.history));

    return true;
}

// ----------------------------------------------------------------------------
// Function _goUpStructureConstruction()
// ----------------------------------------------------------------------------

// This function implements the functionality of go up and
// resizes the borderString of the structure construction.
template <typename TTree, typename TIterSpec, typename TBorderString>
inline bool _goUpStructureConstruction(Iter<TTree, RightArrayBinaryTreeIterator<TopDown<ParentLinks<TIterSpec> > > > & it, TBorderString & borderString)
{
    if (goUp(it))
    {
        resize(borderString, length(it.history));
        return true;
    }
    return false;
}

// ----------------------------------------------------------------------------
// Function isLeaf()
// ----------------------------------------------------------------------------
/*!
 * @fn RightArrayBinaryTreeIterator#isLeaf
 * @headerfile <seqan/index.h>
 * @brief Tests whether a given node is a leaf or not.
 *
 * @signature bool isLeaf(iterator);
 *
 * @param[in] iterator The iterator.
 *
 * @return bool <tt>true</tt> if the node is a leaf.
 */
template <typename TTree, typename TIterSpec>
inline bool isLeaf(Iter<TTree, RightArrayBinaryTreeIterator<TopDown<TIterSpec> > > & iter)
{
    return container(iter).treeVertices[getPosition(iter)].i2 == 0;
}

// ----------------------------------------------------------------------------
// Function _setAndGoRight()
// ----------------------------------------------------------------------------

// This function creates the right sibling of the current node
// and goes to that one.
// Note: It acn only be called, if the right sibling really exists!
template <typename TTree, typename TIterSpec, typename TBorderString>
inline bool _setAndGoRight(Iter<TTree, RightArrayBinaryTreeIterator<TopDown<TIterSpec> > > & it, TBorderString & borderString)
{
    typedef typename Value<typename Value<TTree>::Type, 1>::Type TChar;

    if (isRoot(it) || (back(borderString).i2 == borderString[length(borderString) - 2].i2))
        return false;

    goUp(it);

    if (borderString[length(borderString) - 2].i2 == ordValue(getCharacter(it)))
    {
        goLeftChild(it);
        return false;
    }

    resize(container(it).treeVertices, length(container(it).treeVertices) + 1);
    TChar pivot = getCharacter(it);
    _setRightChildPos(it, length(container(it).treeVertices) - 1);
    goRightChild(it);
    back(borderString).i1 = ordValue(pivot);
    back(borderString).i2 = borderString[length(borderString) - 2].i2;

    return true;
}

// ----------------------------------------------------------------------------
// Function setCharacter()
// ----------------------------------------------------------------------------
/*!
 * @fn RightArrayBinaryTreeIterator#setCharacter
 * @headerfile <seqan/index.h>
 * @brief The function sets the character of the node the iterator points to to character.
 *
 * @signature void setCharacter(iterator, character);
 *
 * @param[in,out] iterator  The iterator.
 * @param[in]     character The character to be assigned to a node.
 */
template <typename TTree, typename TIterSpec, typename TChar2>
inline void setCharacter(Iter<TTree, RightArrayBinaryTreeIterator<TopDown<TIterSpec> > > & iter,
                         TChar2 character)
{
    container(iter).treeVertices[getPosition(iter)].i1 = character;
}

// ----------------------------------------------------------------------------
// Function _getPivotPosition()
// ----------------------------------------------------------------------------

// This function returns the position of the character which ensures that the sum of occurrences of the characters from
// beginPos to the computed pos and the sum of occurrences from the computed pos to endPos are about the same.
template <typename TPrefixSums, typename TBeginPos, typename TEndPos>
TBeginPos
_getPivotPosition(TPrefixSums const & sums, TBeginPos beginPos, TEndPos endPos)
{
    typedef typename Value< TPrefixSums >::Type TValue;

    TBeginPos realBeginPos = beginPos + 1;
    TEndPos realEndPos = endPos + 1;
    TBeginPos lengthRange = realEndPos - realBeginPos + 1;
    TBeginPos pivotPos = realBeginPos + lengthRange / 2 - 1;

    TValue tooSmallValues = sums[beginPos];
    long currentMin = sums[realEndPos] + 1;

    if (sums[pivotPos] - tooSmallValues >= sums[realEndPos] - sums[pivotPos])
    {
        while ((pivotPos >= realBeginPos) && std::abs((long)(sums[pivotPos] - tooSmallValues) - (long)((sums[realEndPos] - sums[pivotPos]))) <= currentMin)
        {
            currentMin = std::abs((long)((sums[pivotPos] - tooSmallValues)) - (long)((sums[realEndPos] - sums[pivotPos])));
            --pivotPos;
        }
        ++pivotPos;
    }
    else
    {
        while (std::abs((long)((sums[pivotPos] - tooSmallValues)) - (long)((sums[realEndPos] - sums[pivotPos]))) < currentMin && (pivotPos < realEndPos))
        {
            currentMin = std::abs((long)((sums[pivotPos] - tooSmallValues)) - (long)((sums[realEndPos] - sums[pivotPos])));
            ++pivotPos;
        }
        --pivotPos;
    }

    return pivotPos;
}

// ----------------------------------------------------------------------------
// Function _setChildVertices()
// ----------------------------------------------------------------------------

// This function sets the left child of the current node, or the right if there is no left child.
template <typename TTree, typename TIterSpec, typename TBorderString, typename TPrefixSums>
void _setChildVertices(Iter<TTree, RightArrayBinaryTreeIterator<TopDown<TIterSpec> > > & it,
                       TBorderString & borderString,
                       TPrefixSums & sums)
{
    typedef typename Size<TTree>::Type TSize;
    typedef typename Value<TBorderString>::Type TBorderStringValue;

    TSize leftBorder = back(borderString).i1;
    TSize rightBorder = back(borderString).i2;
    TSize pivotPosition = _getPivotPosition(sums, leftBorder, rightBorder);

    setCharacter(it, pivotPosition);

    if (leftBorder == pivotPosition - 1)
    {
        // set the right child to be the only one
        container(it).treeVertices[getPosition(it)].i2 = 1;
        appendValue(borderString, TBorderStringValue(pivotPosition, back(borderString).i2));
        return;
    }

    _setLeftChildPos(it);

    appendValue(borderString, TBorderStringValue(back(borderString).i1, pivotPosition - 1));
}

// ----------------------------------------------------------------------------
// Function _setLeftChildPos()
// ----------------------------------------------------------------------------

// This functions sets the pointer to the left child.
template <typename TTree, typename TIterSpec>
inline bool _setLeftChildPos(Iter<TTree, RightArrayBinaryTreeIterator<TopDown<TIterSpec> > > & iter)
{
    switch (iter.waveletTreeStructure->treeVertices[getPosition(iter)].i2)
    {
    case 0:
        iter.waveletTreeStructure->treeVertices[getPosition(iter)].i2 = 2;
        return true;

    case 2:
        return true;

    default:
        return false;
    }
}

// ----------------------------------------------------------------------------
// Function _setRightChildPos()
// ----------------------------------------------------------------------------

// This functions sets the pointer to the left child.
template <typename TTree, typename TPos, typename TIterSpec>
inline bool _setRightChildPos(Iter<TTree, RightArrayBinaryTreeIterator<TopDown<TIterSpec> > > & iter, TPos rightChildPosition)
{
    switch (iter.waveletTreeStructure->treeVertices[getPosition(iter)].i2)
    {
    case 0:
        SEQAN_ASSERT_EQ_MSG(rightChildPosition, 0u, "Wrong right child position!");
        iter.waveletTreeStructure->treeVertices[getPosition(iter)].i2 = 1;
        return true;

    case 2:
        iter.waveletTreeStructure->treeVertices[getPosition(iter)].i2 = rightChildPosition + 2;
        return true;

    case 1:
        SEQAN_ASSERT_MSG(rightChildPosition == 0u, "Wrong right child position!");
        return true;

    default:
        return false;
    }
}

// ----------------------------------------------------------------------------
// Function isRoot()
// ----------------------------------------------------------------------------
/*!
 * @fn RightArrayBinaryTreeIterator#isRoot
 * @headerfile <seqan/index.h>
 * @brief Test whether a iterator points to the root node.
 *
 * @signature bool isRoot(iterator);
 *
 * @param[in] iterator The iterator.
 *
 * @return bool <tt>true</tt> if <tt>iterator</tt> points to the root of the tree, otherwise <tt>false</tt>.
 */
template <typename TTree, typename TIterSpec>
inline bool isRoot(Iter<TTree, RightArrayBinaryTreeIterator<TopDown<TIterSpec> > > const & it)
{
    return getPosition(it) == 0;
}

/*
template <typename TTree, typename TIterSpec, typename TString>
inline void _writeGraphImpl(Iter<TTree, RightArrayBinaryTreeIterator<TopDown<TIterSpec> > > & iter, TString name)
{
    //typedef typename BitVector_<BitsPerValue<typename Value<TText>::Type>::VALUE>::Type TValue;
    Iter<TTree, RightArrayBinaryTreeIterator<TopDown<TIterSpec> > > iter2 = iter;
    std::ofstream stream(toCString(name), std::ios::app);
    unsigned pos = getLeftChildPos(iter);
    if (pos)
    {
        stream << (unsigned)ordValue(iter.waveletTreeStructure->treeVertices[getPosition(iter)].i1) << " -> " << (unsigned)ordValue(iter.waveletTreeStructure->treeVertices[pos].i1) << ";" << std::endl;
        goLeftChild(iter);
        writeGraphImpl(iter, name);
    }
    else
    {
        stream << (unsigned)ordValue(iter.waveletTreeStructure->treeVertices[getPosition(iter)].i1) << " -> " << "leave1" << (unsigned)ordValue(getPosition(iter)) << ";" << std::endl;
    }

    pos = getRightChildPos(iter2);
    if (pos)
    {
        stream << (unsigned) ordValue(iter2.waveletTreeStructure->treeVertices[getPosition(iter2)].i1) << " -> " << (unsigned)ordValue(iter2.waveletTreeStructure->treeVertices[pos].i1) << ";" << std::endl;

        goRightChild(iter2);
        writeGraphImpl(iter2, name);
    }
    else
    {
        stream << (unsigned)ordValue(iter2.waveletTreeStructure->treeVertices[getPosition(iter2)].i1) << " -> " << "leave2" << (unsigned)ordValue(getPosition(iter2)) << ";" << std::endl;
    }
    stream.close();
}

template <typename TTree>
inline void _writeGraph(RightArrayBinaryTree<TChar, TSpec> & treeStructure)
{

    typename Iterator<RightArrayBinaryTree<TChar, TSpec>, TopDown<ParentLinks<> > >::Type iter(treeStructure, 0);

    String<char> name = "testfile.dot";
    std::ofstream stream(toCString(name), std::ios::out);
    stream << "digraph G {" << std::endl;
    stream.close();
    writeGraphImpl(iter, name);

    stream.open(toCString(name), std::ios::app);
    stream << "}" << std::endl;
    stream.close();
}
*/

}
#endif // INDEX_FM_RIGHT_ARRAY_BINARY_TREE_ITERATOR_H
