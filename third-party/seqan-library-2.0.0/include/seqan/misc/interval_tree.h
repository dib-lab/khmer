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
// Author: Anne-Katrin Emde <emde@fu-berlin.de>
// ==========================================================================


// TODO(holtgrew): Bring back random center computation using new random module?

// (weese:) This header needs improvement
//
// documentation lacks:
//   - is createIntervalTree clearing the tree first (what happens after a second call of createIntervalTree?)
//   - where is clear(itree)?
//   - some flags are missing in the dddoc entries, e.g. ..cat:, class:, or ..type: of arguments
//
// interface flaws:
//   - there are 2 interfaces for createIntervalTree
//      - based on IntervalTree
//      - based on Graph and PropertyMap (should be removed)
//   - addInterval


#ifndef SEQAN_HEADER_MISC_INTERVAL_TREE_H
#define SEQAN_HEADER_MISC_INTERVAL_TREE_H

#include <seqan/graph_types.h>


namespace SEQAN_NAMESPACE_MAIN {
//////////////////////////////////////////////////////////////////////////////
// Graph - Interval Tree Types
//////////////////////////////////////////////////////////////////////////////

///---------------------------------------------------------------///

//////////////////// Interval and ID type ///////////////////

// TODO(holtgrew): Are these actual the first/last or begin/end entries?

/*!
 * @class IntervalAndCargo
 * @headerfile <seqan/misc/interval_tree.h>
 * @brief A simple record type that stores an interval and a cargo value.
 *
 * @signature template <[typename TValue[, typename TCargo]]>
 *            class IntervalAndCargo;
 *
 * @tparam TValue The value type.  Default: <tt>int</tt>.
 * @tparam TCargo The cargo type.  Default: <tt>int</tt>.
 *
 * @fn IntervalAndCargo::IntervalAndCargo
 * @brief Constructor
 *
 * @signature IntervalAndCargo::IntervalAndCargo();
 * @signature IntervalAndCargo::IntervalAndCargo(i1, i2, cargo);
 *
 * @param[in] i1    The first element in the interval.
 * @param[in] i2    The last element in the interval.
 * @param[in] cargo The cargo to store together with the interval.
 *
 * @var TValue IntervalAndCargo::i1;
 * @brief The first element in the interval.
 *
 * @var TValue IntervalAndCargo::i2;
 * @brief The last element in the interval.
 *
 * @var TCargo IntervalAndCargo::cargo;
 * @brief The stored cargo.
 */

template <typename TValue = int, typename TCargo = int>
class IntervalAndCargo
{
public:
    TValue i1;

    TValue i2;

    TCargo cargo;

    IntervalAndCargo() :
        i1(), i2(), cargo()
    {}

    IntervalAndCargo(TValue i1, TValue i2, TCargo cargo) :
        i1(i1), i2(i2), cargo(cargo)
    {
        SEQAN_CHECKPOINT
    }

};



/////////////////////// Point and ID type ////////////////

/*!
 * @class PointAndCargo
 * @headerfile <seqan/misc/interval_tree.h>
 * @brief Simple record class storing a point (one-value) interval and cargo.
 *
 * @signature template <[typename TValue[, typename TCargo]]>
 *            class PointAndCargo;
 *
 * @tparam TValue The value type.
 * @tparam TCargo The cargo type.
 *
 * @fn PointAndCargo::PointAndCargo
 * @brief Constructor
 *
 * @signature PointAndCargo::PointAndCargo();
 * @signature PointAndCargo::PointAndCargo(point, cargo);
 *
 * @param[in] point The point to store.
 * @param[in] cargo The cargo to store.
 *
 * @var TValue PointAndCargo::point
 * @brief The point to store.
 *
 * @var TCargo PointAndCargo::cargo
 * @brief The cargo to store.
 */

template <typename TValue = int, typename TCargo = int>
class PointAndCargo
{
public:
    TValue point;

    TCargo cargo;

    PointAndCargo() :
        point(), cargo()
    {}

    PointAndCargo(TValue point, TCargo cargo) :
        point(point), cargo(cargo)
    {}
};

///////////////////////////////////////////////////////////////////////////
/////////////////////////// IntervalTreeNode    ///////////////////////////

/*!
 * @defgroup IntervalTreeNodeTypeTags IntervalTree Node Types Tags
 * @brief Tags to select the node type for @link IntervalTree @endlink.
 *
 * @tag IntervalTreeNodeTypeTags#StorePointsOnly
 * @headerfile <seqan/misc/interval_tree.h>
 * @signature struct StorePointsOnly {};
 * @brief The tree nodes store points.
 *
 * @tag IntervalTreeNodeTypeTags#StoreIntervals
 * @headerfile <seqan/misc/interval_tree.h>
 * @signature struct StoreIntervals {};
 * @brief The tree nodes store intervals.
 */

struct StorePointsOnly {};


struct StoreIntervals {};

/*!
 * @class IntervalTreeNode
 * @headerfile <seqan/misc/interval_tree.h>
 * @brief Element of @link IntervalTree @endlink.
 *
 * @signature template <typename TInterval[, typename TSpec]>
 *            class IntervalTreeNode;
 *
 * @tparam TInterval The type to use for intervals.
 * @tparam TSpec     The specializing tag.  Default: @link IntervalTreeNodeTypeTags#StorePointsOnly @endlink.
 */


template <typename TInterval, typename TSpec = StorePointsOnly>
class IntervalTreeNode;

/*!
 * @class StoreIntervalsIntervalTreeNode
 * @extends IntervalTreeNode
 * @headerfile <seqan/misc/interval_tree.h>
 * @brief An IntervalTreeNode that stores intervals explicitely in each node.
 *
 * @signature template <typename TInterval>
 *            class IntervalTreeNode<TInterval, StoreIntervals>;
 *
 * @tparam TInterval The Interval type to use.
 */

/*!
 * @var TValue StoreIntervalsIntervalTreeNode::center;
 * @brief Center of the interval tree node.
 *
 * @var TString StoreIntervalsIntervalTreeNode::list1;
 * @brief @link AllocString @endlink of intervals sorted by begin point.
 *
 * @var TString StoreIntervalsIntervalTreeNode::list;
 * @brief @link AllocString @endlink of intervals sorted by end point.
 */

template <typename TInterval>
class IntervalTreeNode<TInterval, StoreIntervals>
{
public:
    typedef typename Value<TInterval>::Type TValue;

    TValue center;
    String<TInterval> list1;
    String<TInterval> list2;

    IntervalTreeNode() :
        center()
    {}
};

/*!
 * @class StorePointsOnlyIntervalTreeNode
 * @extends IntervalTreeNode
 * @headerfile <seqan/misc/interval_tree.h>
 * @brief An IntervalTreeNode that stores only the relevant points in each node.
 *
 * Only the end points of the intervals in the list sorted by endpoints (list2) and only the begin point of the interval
 * list sorted by begin points (list1) are stored.
 *
 * @signature template <typename TInterval>
 *            class IntervalTreeNode<TInterval, StoreIntervals>;
 *
 * @tparam TInterval The Interval type to use.
 */

/*!
 * @var TValue StorePointsOnlyIntervalTreeNode::center;
 * @brief Center of the interval.
 *
 * @var TString StorePointsOnlyIntervalTreeNode::list1;
 * @brief Points with cargo sorted by the begin points.
 *
 * @var TString StorePointsOnlyIntervalTreeNode::list2;
 * @brief Points with cargo sorted by the end points.
 */

template <typename TInterval>
class IntervalTreeNode<TInterval, StorePointsOnly>
{
public:
    typedef typename Cargo<TInterval>::Type TCargo;
    typedef typename Value<TInterval>::Type TValue;

    TValue center;
    String<PointAndCargo<TValue, TCargo> > list1;
    String<PointAndCargo<TValue, TCargo> > list2;

    IntervalTreeNode() :
        center()
    {
        SEQAN_CHECKPOINT
    }

    IntervalTreeNode(IntervalTreeNode const & other) :
        center(other.center),
        list1(other.list1),
        list2(other.list2)
    {
        SEQAN_CHECKPOINT
    }

};





//////////////////////////////////////////////////////////////////////////////
// Graph - Interval Tree
//////////////////////////////////////////////////////////////////////////////

/*!
 * @class IntervalTree
 * @headerfile <seqan/misc/interval_tree.h>
 * @brief Data structure for efficient interval storage.
 *
 * @signature template <[typename TValue[, typename TCargo]]>
 *            class IntervalTree;
 *
 * @tparam TValue The type to use for coordinates.  Default: <tt>int</tt>.
 * @tparam TCargo The type to use for cargo.  Default: <tt>unsigned</tt>.
 *
 * @section Remarks
 *
 * If the intervals are not associated with cargos/IDs, they will be numbered consecutively.
 *
 * @section Example
 *
 * The following example creates an integer interval tree with string keys.  This tree is quired for keys of intervals
 * that overlap the interval <tt>[550, 990)</tt>.
 *
 * @include demos/misc/interval_tree_example.cpp
 *
 * The resulting keys are:
 *
 * @code{.console}
 * gene
 * exon2
 * coding2
 * @endcode
 *
 *
 * @fn IntervalTree::IntervalTree
 * @brief Constructor
 *
 * @signature IntervalTree::IntervalTree();
 * @signature IntervalTree::IntervalTree(intervals);
 * @signature IntervalTree::IntervalTree(intervals[, center]);
 * @signature IntervalTree::IntervalTree(intervals[, tag]);
 * @signature IntervalTree::IntervalTree(intervalBegins, intervalEnds, [intervalCargos,] len);
 *
 * @param[in] intervals Container of intervals.  A strin gof <tt>IntervalAndCargo&lt;Value, TCargo&gt;</tt>
 *                      objects, see @link IntervalAndCargo @endlink.
 * @param[in] intervalBegins
 *                      Iterator pointing to begin position of first interval.
 * @param[in] intervalEnds
 *                      Iterator pointing to end position of first interval.
 * @param[in] intervalCargos
 *                      Iterator pointing to cargo/ids for intervals.
 * @param[in] len       Number of intervals to store in tree.
 * @param[in] tag       Tag for tree construction method.
 */

template <typename TValue = int, typename TCargo = unsigned int>
class IntervalTree
{
public:
    typedef Graph<Directed<void, WithoutEdgeId> > TGraph;
    typedef IntervalAndCargo<TValue, TCargo> TInterval;
    typedef IntervalTreeNode<TInterval> TNode;
    typedef String<TNode> TPropertyMap;

    TGraph g;
    TPropertyMap pm;
    size_t interval_counter;

    IntervalTree()
    {
        SEQAN_CHECKPOINT
            interval_counter = 0;
    }

    template <typename TIterator, typename TCargoIterator>
    IntervalTree(TIterator interval_begins,
                 TIterator interval_ends,
                 TCargoIterator interval_cargos,
                 size_t len)
    {
        SEQAN_CHECKPOINT
        String<TInterval> intervals;
        resize(intervals, len);
        size_t i = 0;
        while (i < len)
        {
            intervals[i].i1 = value(interval_begins);
            ++interval_begins;
            intervals[i].i2 = value(interval_ends);
            ++interval_ends;
            intervals[i].cargo = value(interval_cargos);
            ++interval_cargos;
            ++i;
        }
        createIntervalTree(*this, intervals);
    }

    template <typename TIterator>
    IntervalTree(TIterator interval_begins,
                 TIterator interval_ends,
                 size_t len)
    {
        SEQAN_CHECKPOINT
        String<TInterval> intervals;
        resize(intervals, len);
        size_t i = 0;
        while (i < len)
        {
            intervals[i].i1 = value(interval_begins);
            ++interval_begins;
            intervals[i].i2 = value(interval_ends);
            ++interval_ends;
            intervals[i].cargo = i;
            ++i;
        }
        createIntervalTree(*this, intervals);
    }

    IntervalTree(String<TInterval> intervals)
    {
        SEQAN_CHECKPOINT
            createIntervalTree(*this, intervals);
    }

    template <typename TTagSpec>
    IntervalTree(String<TInterval> intervals, Tag<TTagSpec> const tag)
    {
        SEQAN_CHECKPOINT
            interval_counter = length(intervals);
        createIntervalTree(g, pm, intervals, tag);
    }

    IntervalTree(String<TInterval> intervals, TValue center)
    {
        SEQAN_CHECKPOINT
            interval_counter = length(intervals);
        createIntervalTree(g, pm, intervals, center);
    }

};



///////Specs for the way interval centers are determined

//template <typename TSpec = SpecPointAndCargo>
struct TagComputeCenter_;
typedef Tag<TagComputeCenter_> const ComputeCenter;


///////////////////////////////////////////////////////////////////////////
///////////////////// IntervalAndCargo functions //////////////////////////
///////////////////////////////////////////////////////////////////////////

/*!
 * @fn IntervalAndCargo#leftBoundary
 * @brief Access to left boundary.
 *
 * @signature TBoundary leftBoundary(interval);
 *
 * @param[in] interval The IntervalAndCargo to query for its left boundary.
 *
 * @return TBoundary Reference to the left boundary value.
 */

template <typename TValue, typename TCargo>
TValue &
leftBoundary(IntervalAndCargo<TValue, TCargo> & interval)
{
    return interval.i1;
}

template <typename TValue, typename TCargo>
TValue const &
leftBoundary(IntervalAndCargo<TValue, TCargo> const & interval)
{
    return interval.i1;
}

/*!
 * @fn IntervalAndCargo#rightBoundary
 * @brief Access to right boundary.
 *
 * @signature TBoundary rightBoundary(interval);
 *
 * @param[in] interval The IntervalAndCargo to query for its right boundary.
 *
 * @return TBoundary Reference to the right boundary value.
 */

template <typename TValue, typename TCargo>
TValue &
rightBoundary(IntervalAndCargo<TValue, TCargo> & interval)
{
    return interval.i2;
}

template <typename TValue, typename TCargo>
TValue const &
rightBoundary(IntervalAndCargo<TValue, TCargo> const & interval)
{
    return interval.i2;
}

/*!
 * @fn IntervalAndCargo#getLeftBoundary
 * @brief Access to getLeft boundary.
 *
 * @signature TBoundary getLeftBoundary(interval);
 *
 * @param[in] interval The IntervalAndCargo to query for its left boundary.
 *
 * @return TBoundary Copy of the left boundary value.
 */

template <typename TValue, typename TCargo>
TValue
getLeftBoundary(IntervalAndCargo<TValue, TCargo> const & interval)
{
    return interval.i1;
}

/*!
 * @fn IntervalAndCargo#getRightBoundary
 * @brief Access to getRight boundary.
 *
 * @signature TBoundary getRightBoundary(interval);
 *
 * @param[in] interval The IntervalAndCargo to query for its right boundary.
 *
 * @return TBoundary Copy of the right boundary value.
 */

template <typename TValue, typename TCargo>
TValue
getRightBoundary(IntervalAndCargo<TValue, TCargo> const & interval)
{
    return interval.i2;
}

/*!
 * @fn IntervalAndCargo#cargo
 * @brief Access to the cargo.
 *
 * @signature TCargo cargo(interval);
 *
 * @param[in] interval The IntervalAndCargo to query for its cargo.
 *
 * @return TCargo Reference to the cargo member.
 */

template <typename TValue, typename TCargo>
TCargo const &
cargo(IntervalAndCargo<TValue, TCargo> const & interval)
{
    return interval.cargo;
}

template <typename TValue, typename TCargo>
TCargo &
cargo(IntervalAndCargo<TValue, TCargo> & interval)
{
    return interval.cargo;
}

/*!
 * @fn IntervalAndCargo#getCargo
 * @brief Access to the cargo.
 *
 * @signature TCargo getCargo(interval);
 *
 * @param[in] interval The IntervalAndCargo to query for its cargo.
 *
 * @return TCargo Copy of the cargo member.
 */

template <typename TValue, typename TCargo>
TCargo
getCargo(IntervalAndCargo<TValue, TCargo> const & interval)
{
    return interval.cargo;
}

/////////////////// Metafunctions //////////////////////

/*!
 * @mfn IntervalAndCargo#Value
 * @brief Return the value type.
 *
 * @signature Value<TIntervalAndCargo>::Type;
 */

template <typename TValue, typename TCargo>
struct Value<IntervalAndCargo<TValue, TCargo> >
{
    typedef TValue Type;
};

/*!
 * @mfn IntervalAndCargo#Cargo
 * @brief Return the cargo type.
 *
 * @signature Cargo<TIntervalAndCargo>::Type;
 */

template <typename TValue, typename TCargo>
struct Cargo<IntervalAndCargo<TValue, TCargo> >
{
    typedef TCargo Type;
};


///////////////////////////////////////////////////////////////////////////
///////////////////// PointAndCargo functions /////////////////////////////
///////////////////////////////////////////////////////////////////////////

/*!
 * @fn PointAndCargo#leftBoundary
 * @brief Access to left boundary.
 *
 * @signature TBoundary leftBoundary(point);
 *
 * @param[in] point The PointAndCargo to query for its left boundary.
 *
 * @return TBoundary Reference to the left boundary value.
 */

template <typename TValue, typename TCargo>
TValue const &
leftBoundary(PointAndCargo<TValue, TCargo> const & point)
{
    return point.point;
}

template <typename TValue, typename TCargo>
TValue &
leftBoundary(PointAndCargo<TValue, TCargo> & point)
{
    return point.point;
}

/*!
 * @fn PointAndCargo#rightBoundary
 * @brief Access to right boundary.
 *
 * @signature TBoundary rightBoundary(point);
 *
 * @param[in] point The PointAndCargo to query for its right boundary.
 *
 * @return TBoundary Reference to the right boundary value.
 */

template <typename TValue, typename TCargo>
TValue const &
rightBoundary(PointAndCargo<TValue, TCargo> const & point)
{
    return point.point;
}

template <typename TValue, typename TCargo>
TValue &
rightBoundary(PointAndCargo<TValue, TCargo> & point)
{
    return point.point;
}

template <typename TValue, typename TCargo>
TValue
getLeftBoundary(PointAndCargo<TValue, TCargo> const & point)
{
    return point.point;
}

/*!
 * @fn PointAndCargo#getLeftBoundary
 * @brief Access to getLeft boundary.
 *
 * @signature TBoundary getLeftBoundary(point);
 *
 * @param[in] point The PointAndCargo to query for its left boundary.
 *
 * @return TBoundary Copy of the left boundary value.
 */

template <typename TValue, typename TCargo>
TValue
getRightBoundary(PointAndCargo<TValue, TCargo> const & point)
{
    return point.point;
}

/*!
 * @fn PointAndCargo#cargo
 * @brief Access to the cargo.
 *
 * @signature TCargo cargo(point);
 *
 * @param[in] point The PointAndCargo to query for its cargo.
 *
 * @return TCargo Reference to the cargo member.
 */


template <typename TValue, typename TCargo>
TCargo const &
cargo(PointAndCargo<TValue, TCargo> const & point)
{
    return point.cargo;
}

template <typename TValue, typename TCargo>
TCargo &
cargo(PointAndCargo<TValue, TCargo> & point)
{
    return point.cargo;
}

/*!
 * @fn PointAndCargo#getCargo
 * @brief Access to the cargo.
 *
 * @signature TCargo getCargo(point);
 *
 * @param[in] point The PointAndCargo to query for its cargo.
 *
 * @return TCargo Copy of the cargo member.
 */

template <typename TValue, typename TCargo>
TCargo
getCargo(PointAndCargo<TValue, TCargo> const & point)
{
    return point.cargo;
}

////////////////// Metafunctions //////////////////

/*!
 * @mfn PointAndCargo#Value
 * @brief Return the value type.
 *
 * @signature Value<TPointAndCargo>::Type;
 */

template <typename TValue, typename TCargo>
struct Value<PointAndCargo<TValue, TCargo> >
{
    typedef TValue Type;
};

/*!
 * @mfn PointAndCargo#Cargo
 * @brief Return the cargo type.
 *
 * @signature Cargo<TPointAndCargo>::Type;
 */

template <typename TValue, typename TCargo>
struct Cargo<PointAndCargo<TValue, TCargo> >
{
    typedef TCargo Type;
};


//// Comparators
template <typename TPair>
bool _less_compI1_ITree(TPair const & p1, TPair const & p2)
{
    return leftBoundary(p1) < leftBoundary(p2);
}

template <typename TPair>
bool _greater_compI2_ITree(TPair const & p1, TPair const & p2)
{
    return rightBoundary(p1) > rightBoundary(p2);
}

///////////////////////////////////////////////////////////////////////////
///////////////////// IntervalTreeNode functions //////////////////////////
///////////////////////////////////////////////////////////////////////////





// internal set node functions
template <typename TValue, typename TInterval>
void
_setIntervalTreeNode(IntervalTreeNode<TInterval, StoreIntervals> & knot, TValue center, TInterval const & interval)
{
    knot.center = center;
    appendValue(knot.list1, interval);
    appendValue(knot.list2, interval);
}

// append intervals to lists in node knot
template <typename TInterval>
void
_appendIntervalTreeNodeLists(IntervalTreeNode<TInterval, StoreIntervals> & knot, TInterval const & interval)
{
    appendValue(knot.list1, interval);
    appendValue(knot.list2, interval);
}

//internal set node functions
template <typename TValue, typename TInterval>
void
_setIntervalTreeNode(IntervalTreeNode<TInterval, StorePointsOnly> & knot, TValue center, TInterval const & interval)
{
    knot.center = center;
    appendValue(knot.list1, PointAndCargo<TValue, typename Cargo<TInterval>::Type>(leftBoundary(interval), cargo(interval)));
    appendValue(knot.list2, PointAndCargo<TValue, typename Cargo<TInterval>::Type>(rightBoundary(interval), cargo(interval)));
}

template <typename TInterval>
void
_appendIntervalTreeNodeLists(IntervalTreeNode<TInterval, StorePointsOnly> & knot, TInterval const & interval)
{
    appendValue(knot.list1, PointAndCargo<typename Value<TInterval>::Type, typename Cargo<TInterval>::Type>(leftBoundary(interval), cargo(interval)));
    appendValue(knot.list2, PointAndCargo<typename Value<TInterval>::Type, typename Cargo<TInterval>::Type>(rightBoundary(interval), cargo(interval)));


}

/////////////////// Metafunctions ///////////////////////

/*!
 * @mfn IntervalTreeNode#Value
 * @brief Return value type.
 *
 * @signature Value<TNode>::Type;
 */

template <typename TInterval, typename TSpec>
struct Value<IntervalTreeNode<TInterval, TSpec> >
{
    typedef typename Value<TInterval>::Type Type;
};

/*!
 * @mfn IntervalTreeNode#Cargo
 * @brief Return cargo type.
 *
 * @signature Cargo<TNode>::Type;
 */

template <typename TInterval, typename TSpec>
struct Cargo<IntervalTreeNode<TInterval, TSpec> >
{
    typedef typename Cargo<TInterval>::Type Type;
};


/*!
 * @mfn IntervalTreeNode#ListType
 * @brief Type of the lists in tree nodes.
 *
 * @signature ListType<T>::Type;
 */

template <typename T>
struct ListType;


template <typename TInterval>
struct ListType<IntervalTreeNode<TInterval, StorePointsOnly> >
{
    typedef String<PointAndCargo<typename Value<TInterval>::Type, typename Cargo<TInterval>::Type> > Type;

};


template <typename TInterval>
struct ListType<IntervalTreeNode<TInterval, StoreIntervals> >
{
    typedef String<IntervalAndCargo<typename Value<TInterval>::Type, typename Cargo<TInterval>::Type> > Type;

};



///////////////////////////////////////////////////////////////////////////
/////////////////////// IntervalTree functions ////////////////////////////
///////////////////////////////////////////////////////////////////////////

/*!
 * @fn IntervalTree#createIntervalTree
 * @brief Create an interval tree.
 *
 * @signature void createIntervalTree(intervalTree, intervals[, tag]);
 * @signature void createIntervalTree(g, pm, intervals[, tag]);
 * @signature void createIntervalTree(g, pm, intervals, center[, tag]]);
 *
 * @param[in,out] intervalTree An interval tree Types: IntervalTree
 * @param[in,out] g            DirectedGraph to create interval tree in. Types: @link Graph @endlink.
 * @param[in,out] pm           Property map to use for the created interval tree.
 * @param[in]     tag          Tag for tree construction method;
 * @param[in]     intervals    Container of intervals.  A string of <tt>IntervalAndCargo&lt;TValue, TCargo&gt;</tt>
 *                             objects, see @link IntervalAndCargo @endlink. Types: @link AllocString @endlink.
 */

template <typename TGraph, typename TPropertyMap, typename TIntervals, typename TSpec>
inline void
createIntervalTree(TGraph & g,
                   TPropertyMap & pm,
                   TIntervals & intervals,
                   Tag<TSpec> const tag)
{
    typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
    typedef typename Value<TIntervals>::Type TInterval;
    typedef typename Value<TInterval>::Type TValue;

    reserve(g.data_vertex, length(intervals));
    reserve(pm, length(intervals));

    TVertexDescriptor root = addVertex(g);
    resizeVertexMap(pm, g);

    if (length(intervals) > 0u)
    {
        TValue center = _calcIntervalTreeRootCenter(intervals);

        std::sort(begin(intervals, Standard()), end(intervals, Standard()), _less_compI1_ITree<TInterval>);

        String<TInterval *> interval_pointers;
        // interval tree stores pointers to intervals, not original intervals
        _makePointerInterval(intervals, interval_pointers);

        _createIntervalTree(g, pm, interval_pointers, root, (TValue)0.0, center, length(intervals), tag);
        reserve(pm, length(pm), Exact());
        reserve(g.data_vertex, length(g.data_vertex), Exact());
    }
}

template <typename TGraph, typename TPropertyMap, typename TIntervals>
inline void
createIntervalTree(TGraph & g,
                   TPropertyMap & pm,
                   TIntervals & intervals)
{
    createIntervalTree(g, pm, intervals, ComputeCenter());
}

template <typename TGraph, typename TPropertyMap, typename TIntervals, typename TSpec>
inline void
createIntervalTree(TGraph & g,
                   TPropertyMap & pm,
                   TIntervals & intervals,
                   typename Value<typename Value<TIntervals>::Type>::Type center,
                   Tag<TSpec> const tag)
{
    typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
    typedef typename Value<TIntervals>::Type TInterval;
    typedef typename Value<typename Value<TIntervals>::Type>::Type TValue;

    reserve(g.data_vertex, length(intervals));
    reserve(pm, length(intervals));

    TVertexDescriptor root = addVertex(g);
    resizeVertexMap(pm, g);

    if (length(intervals) > 0u)
    {
        TInterval a;
        typename Iterator<TIntervals, Standard>::Type begin_ = begin(intervals, Standard());
        typename Iterator<TIntervals, Standard>::Type end_ = end(intervals, Standard());
        std::sort(begin_, end_, _less_compI1_ITree<TInterval>);

        String<TInterval const *> interval_pointers;
        _makePointerInterval(intervals, interval_pointers);

        if (length(intervals) == 1) // if there is just one interval ->  center = center of this interval
            center = (rightBoundary(intervals[0]) - leftBoundary(intervals[0])) / (TValue)2.0;

        _createIntervalTree(g, pm, interval_pointers, root, (TValue)0.0, center, length(intervals), tag);

        reserve(pm, length(pm), Exact());
        reserve(g.data_vertex, length(g.data_vertex), Exact());
    }
}

// ComputeCenter tag as default construction method
template <typename TGraph, typename TPropertyMap, typename TIntervals>
inline void
createIntervalTree(TGraph & g, TPropertyMap & pm,
                   TIntervals & intervals,
                   typename Value<typename Value<TIntervals>::Type>::Type center)
{
    createIntervalTree(g, pm, intervals, center, ComputeCenter());
}

template <typename TValue, typename TCargo, typename TIntervals, typename TSpec>
inline void
createIntervalTree(IntervalTree<TValue, TCargo> & it,
                   TIntervals & intervals,
                   Tag<TSpec> const tag)
{
    it.interval_counter = length(intervals);
    createIntervalTree(it.g, it.pm, intervals, tag);
}

template <typename TValue, typename TCargo, typename TIntervals>
inline void
createIntervalTree(IntervalTree<TValue, TCargo> & it,
                   TIntervals & intervals)
{
    createIntervalTree(it, intervals, ComputeCenter());
}

//////////////////////////////////////////////////////////////////////////////
//remembers minimum and maximum of point values in intervals and sets the center
//of each node to min+(max-min)/2
template <typename TGraph, typename TPropertyMap, typename TIntervalPointer, typename TValue>
inline void
_createIntervalTree(TGraph & g, TPropertyMap & pm,
                    String<TIntervalPointer *> & intervals,
                    typename VertexDescriptor<TGraph>::Type & knot,
                    TValue,
                    TValue center,
                    typename VertexDescriptor<TGraph>::Type len,
                    Tag<TagComputeCenter_> const tag)
{
    //  Rekursionsanker
    if (len == 1)
    {
        _setIntervalTreeNode(value(pm, knot), center, *intervals[0]);
        return;
    }

    typedef typename Value<TPropertyMap>::Type TNode;
    typedef typename ListType<TNode>::Type TList;
    typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
    typedef String<TIntervalPointer *> TIntervalPointers;

    // one list of interval pointers for the intervals to the left of center
    TIntervalPointers S_left;
    // one list of interval pointers for the intervals to the right of center
    TIntervalPointers S_right;

    TValue min1 = maxValue<TValue>();
    TValue min2 = maxValue<TValue>();
    TValue max1 = minValue<TValue>();
    TValue max2 = minValue<TValue>();

    value(pm, knot).center = center;


    typedef typename Iterator<TIntervalPointers, Standard>::Type TIntervalIterator;
    TIntervalIterator it = begin(intervals, Standard());
    TIntervalIterator it_end = end(intervals, Standard());

    // walk through intervals
    while (it != it_end)
    {
        // interval belongs to the left list
        if ((**it).i2 <= center)
        {
            appendValue(S_left, *it, Generous());
            //remember right most and left most point in left list
            if ((**it).i2 > max1)
                max1 = (**it).i2;
            if ((**it).i1 < min1)
                min1 = (**it).i1;
        }
        else
        {
            // interval belongs to the right list
            if ((**it).i1 > center)
            {
                appendValue(S_right, (*it), Generous());
                //remember right most and left most point in right list
                if ((**it).i2 > max2)
                    max2 = (**it).i2;
                if ((**it).i1 < min2)
                    min2 = (**it).i1;
            }
            else // interval belongs to this node
            {
                _appendIntervalTreeNodeLists(value(pm, knot), **it);
            }
        }
        ++it;
    }

//    std::sort(begin(value(pm,knot).list1),end(value(pm,knot).list1),_less_compI1_ITree<typename Value<TList>::Type>);
    std::sort(begin(value(pm, knot).list2), end(value(pm, knot).list2), _greater_compI2_ITree<typename Value<TList>::Type>);

    // build subtree to the left
    if (!empty(S_left))
    {
        TVertexDescriptor vd = addVertex(g);
        resize(pm, vd + 1);
        addEdge(g, knot, vd);
        _createIntervalTree(g, pm, S_left, vd, center, min1 + (max1 - min1) / 2, length(S_left), tag);
    }
    // build subtree to the right
    if (!empty(S_right))
    {
        TVertexDescriptor vd = addVertex(g);
        resize(pm, vd + 1);
        addEdge(g, knot, vd);
        _createIntervalTree(g, pm, S_right, vd, center, min2 + (max2 - min2) / 2, length(S_right), tag);
    }
}

//////////////////////////////////////////////////////////////////////////////
//createIntervalTree for all specs except CompCenter, the center value of each
//node is determined by functions _calcIntervalTreeNodeCenterLeft and
//_calcIntervalTreeNodeCenterRight
template <typename TGraph, typename TPropertyMap, typename TSpec, typename TInterval, typename TValue>
inline void
_createIntervalTree(TGraph & g, TPropertyMap & pm,
                    String<TInterval *> & intervals,
                    typename VertexDescriptor<TGraph>::Type & knot,
                    TValue last_center, TValue center,
                    typename VertexDescriptor<TGraph>::Type len,
                    Tag<TSpec> const tag)
{
    SEQAN_CHECKPOINT
    // Rekursionsanker
    if (len == 1)
    {
        _setIntervalTreeNode(value(pm, knot), center, *value(intervals, 0));
        return;
    }

    typedef typename Value<TPropertyMap>::Type TNode;
    typedef typename ListType<TNode>::Type TList;
    typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
    typedef String<TInterval *> TIntervalPointers;

    // one list of interval pointers for the intervals to the left of center
    TIntervalPointers S_left;
    // one list of interval pointers for the intervals to the right of center
    TIntervalPointers S_right;

    value(pm, knot).center = center;

    typedef typename Iterator<TIntervalPointers, Standard>::Type TIntervalIterator;
    TIntervalIterator it = begin(intervals, Standard());
    TIntervalIterator it_end = end(intervals, Standard());

    // walk through intervals
    while (it != it_end)
    {
        // interval belongs to the left list
        if ((**it).i2 <= center)
        {
            appendValue(S_left, *it, Generous());
        }
        else
        {   // interval belongs to the right list
            if ((**it).i1 > center)
                appendValue(S_right, (*it), Generous());
            else
                // interval belongs to the current node
                _appendIntervalTreeNodeLists(value(pm, knot), **it);
        }
        ++it;
    }

//    std::sort(begin(value(pm,knot).list1),end(value(pm,knot).list1),_less_compI1_ITree<typename Value<TList>::Type>);
    std::sort(begin(value(pm, knot).list2), end(value(pm, knot).list2), _greater_compI2_ITree<typename Value<TList>::Type>);

    // build subtree to the left
    if (!empty(S_left))
    {
        TVertexDescriptor vd = addVertex(g);
        resizeVertexMap(pm, g);
        addEdge(g, knot, vd);
        TValue next_center = _calcIntervalTreeNodeCenterLeft(S_left, last_center, center, tag);
        _createIntervalTree(g, pm, S_left, vd, center, next_center, length(S_left), tag);
    }
    // build subtree to the right
    if (!empty(S_right))
    {
        TVertexDescriptor vd = addVertex(g);
        resizeVertexMap(pm, g);
        addEdge(g, knot, vd);
        TValue next_center = _calcIntervalTreeNodeCenterRight(S_right, last_center, center, tag);
        _createIntervalTree(g, pm, S_right, vd, center, next_center, length(S_right), tag);
    }
}

// fill the container interval_pointers with pointers to the corresponding objects in intervals.
// this is done to avoid copying and passing the whole IntervalAndCargo objects during interval tree construction
template <typename TIntervals, typename TIntervalPointers>
void
_makePointerInterval(TIntervals & intervals, TIntervalPointers & interval_pointers)
{
    typedef typename Iterator<TIntervals, Standard>::Type TIntervalIterator;
    typedef typename Iterator<TIntervalPointers, Standard>::Type TIntervalPointerIterator;

    resize(interval_pointers, length(intervals));
    if (empty(intervals))
        return;

    TIntervalIterator it = begin(intervals, Standard());
    TIntervalIterator itEnd = end(intervals, Standard());
    TIntervalPointerIterator iit = begin(interval_pointers, Standard());

    for (; it != itEnd; ++it, ++iit)
        *iit = it;
}

// if the center of the root is not given, it is placed in the "ComputeCenter way": in the middle of minValue and maxValue
// where minValue is the minimum left boundary and maxValue is the maximum right boundary of all intervals
template <typename TIntervals>
typename Value<typename Value<TIntervals>::Type>::Type
_calcIntervalTreeRootCenter(TIntervals & intervals)
{
    SEQAN_CHECKPOINT;

    SEQAN_ASSERT_GT(length(intervals), 0u);

    typedef typename Value<typename Value<TIntervals>::Type>::Type TValue;
    typedef typename Iterator<TIntervals, Standard>::Type TIntervalIterator;

    TIntervalIterator it = begin(intervals);
    TIntervalIterator it_end = end(intervals);

    TValue min = maxValue<TValue>();
    TValue max = minValue<TValue>();

    // get min and max
    while (it != it_end)
    {
        if (leftBoundary(*it) < min)
            min = leftBoundary(*it);
        if (rightBoundary(*it) > max)
            max = rightBoundary(*it);
        SEQAN_ASSERT_LEQ(min, max);
        ++it;
    }

    SEQAN_ASSERT_LEQ(min, max);

    // return middle between max and min
    return min + (max - min) / (TValue)2.0;

}

/*!
 * @fn IntervalTree#addInterval
 *
 * @headerfile <seqan/misc/interval_tree.h>
 *
 * @brief Adds an interval to an interval tree.
 *
 * @signature void addInterval(intervalTree, interval);
 * @signature void addInterval(intervalTree, begin, end[, cargo]);
 * @signature void addInterval(graph, propertyMap, interval);
 *
 * @param[in,out] intervalTree The interval tree to add the interval to. Types: @link IntervalTree @endlink.
 * @param[in]     interval     The interval to be added to the interval tree.
 * @param[in]     begin        Begin position of interval of type TValue.
 * @param[in]     end          End position of interval of type TValue.
 * @param[in]     cargo        Cargo to attach to the interval. Types: @link IntervalAndCargo @endlink.
 * @param[in,out] graph        The directed graph that contains the topography of the interval tree.
 * @param[in,out] propertyMap  The property map containing the node properties of the interval tree.
 */


template <typename TGraph, typename TPropertyMap, typename TInterval>
void
addInterval(TGraph & g, TPropertyMap & pm, TInterval interval)
{
    SEQAN_CHECKPOINT

    typedef typename Iterator<TGraph, OutEdgeIterator>::Type TOutEdgeIterator;
    typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
    typedef typename Value<TPropertyMap>::Type TProperty;
    typedef typename Value<TInterval>::Type TValue;
    typedef typename ListType<TProperty>::Type TList;


    if (empty(pm))
    {
        TVertexDescriptor vd = addVertex(g);
        resizeVertexMap(pm, g);
        _setIntervalTreeNode(property(pm, vd), (rightBoundary(interval) + leftBoundary(interval)) / 2, interval);
        return;

    }
    // start at root
    TVertexDescriptor act_knot = 0;
    TProperty act_prop = property(pm, act_knot);
    TProperty next_prop;

    // look for the right node to add interval to
    while (true)
    {
        TOutEdgeIterator it(g, act_knot);
        act_prop = property(pm, act_knot);
        if (act_prop.center < leftBoundary(interval)) // interval to the left of current node?
        {
            if (atEnd(it))
            {
                TVertexDescriptor vd = addVertex(g);
                resizeVertexMap(pm, g);
                addEdge(g, act_knot, vd);
                _setIntervalTreeNode(property(pm, vd), (rightBoundary(interval) + leftBoundary(interval)) / (TValue)2.0, interval);
                break;
            }
            else
            {
                next_prop = property(pm, targetVertex(it));
                if (next_prop.center <= act_prop.center)
                {
                    goNext(it);
                    if (atEnd(it))
                    {
                        TVertexDescriptor vd = addVertex(g);
                        resizeVertexMap(pm, g);
                        addEdge(g, act_knot, vd);
                        _setIntervalTreeNode(property(pm, vd), (rightBoundary(interval) + leftBoundary(interval)) / (TValue)2.0, interval);
                        break;
                    }
                }
            }
            act_knot = targetVertex(it);
        }
        else
        {
            if (rightBoundary(interval) <= act_prop.center) // interval to the right of current node?
            {
                if (atEnd(it))
                {
                    TVertexDescriptor vd = addVertex(g);
                    resizeVertexMap(pm, g);
                    addEdge(g, act_knot, vd);
                    _setIntervalTreeNode(property(pm, vd), (rightBoundary(interval) + leftBoundary(interval)) / 2, interval);
                    break;
                }
                else
                {
                    next_prop = property(pm, targetVertex(it));
                    if (next_prop.center >= act_prop.center)
                    {
                        goNext(it);
                        if (atEnd(it))
                        {
                            TVertexDescriptor vd = addVertex(g);
                            resizeVertexMap(pm, g);
                            addEdge(g, act_knot, vd);
                            _setIntervalTreeNode(property(pm, vd), (rightBoundary(interval) + leftBoundary(interval)) / 2, interval);
                            break;
                        }
                    }
                }
                act_knot = targetVertex(it);
            }
            else  // need to create new node for interval
            {
                _appendIntervalTreeNodeLists(property(pm, act_knot), interval);
                std::sort(begin(property(pm, act_knot).list1), end(property(pm, act_knot).list1), _less_compI1_ITree<typename Value<TList>::Type>);
                std::sort(begin(property(pm, act_knot).list2), end(property(pm, act_knot).list2), _greater_compI2_ITree<typename Value<TList>::Type>);
                break;
            }
        }
    }

}

template <typename TValue, typename TCargo, typename TInterval>
void
addInterval(IntervalTree<TValue, TCargo> & itree, TInterval interval)
{
    SEQAN_CHECKPOINT

    ++ itree.interval_counter;
    addInterval(itree.g, itree.pm, interval);

}

template <typename TValue, typename TCargo>
void
addInterval(IntervalTree<TValue, TCargo> & itree, TValue begin, TValue end, TCargo cargo)
{
    SEQAN_CHECKPOINT

    IntervalAndCargo<TValue, TCargo> interval;
    interval.i1 = begin;
    interval.i2 = end;
    interval.cargo = cargo;
    ++itree.interval_counter;
    addInterval(itree.g, itree.pm, interval);

}

template <typename TValue, typename TCargo>
void
addInterval(IntervalTree<TValue, TCargo> & itree, TValue begin, TValue end)
{
    SEQAN_CHECKPOINT

    IntervalAndCargo<TValue, TCargo> interval;
    interval.i1 = begin;
    interval.i2 = end;
    interval.cargo = itree.interval_counter;
    ++itree.interval_counter;
    addInterval(itree.g, itree.pm, interval);

}

/*!
 * @fn IntervalTree#findIntervals
 * @brief Find all intervals that contain the query point or overlap with the query interval.
 *
 * @signature void findIntervals(result, intervalTree, query);
 * @signature void findIntervals(result, intervalTree, queryBegin, queryEnd);
 * @signature void findIntervals(result, graph, propertyMap, query);
 *
 * @param[out] result       A reference to the result string of <tt>TCargo</tt> objects. Types: @link String @endlink.
 * @param[in]  intervalTree An IntervalTree.
 * @param[in]  graph        The directed @link Graph graph @endlink that contains the topography of the interval tree.
 * @param[in]  propertyMap  The property map containing the node properties of the interval tree.
 * @param[in]  query        A query point.
 * @param[in]  queryBegin   The begin position of the query interval.
 * @param[in]  queryEnd     The end position of the query interval.
 */

template <typename TSpec, typename TPropertyMap, typename TValue, typename TCargo>
inline void
findIntervals(
        String<TCargo> & result,
        Graph<TSpec> const & g,
        TPropertyMap const & pm,
        TValue query)
{
    typedef Graph<TSpec> const TGraph;
    typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
    typedef typename Value<TPropertyMap>::Type TProperty;
    typedef typename Value<TProperty>::Type TPropertyValue;
    typedef typename Iterator<TGraph, OutEdgeIterator>::Type TOutEdgeIterator;

    resize(result, 0);
    if (empty(g))
        return;

    // start at root
    TVertexDescriptor act_knot = 0;
    TProperty act_prop = property(pm, act_knot);
    TProperty next_prop;

    while (true)
    {
        typename Iterator<Graph<TSpec>, OutEdgeIterator>::Type it7;
        Iter<Graph<TSpec>, GraphIterator<InternalOutEdgeIterator<OutEdgeIterator> > > it5(g, act_knot);
        TOutEdgeIterator it4;
        TOutEdgeIterator it(g, act_knot);
        act_prop = property(pm, act_knot);
        if (act_prop.center < (TPropertyValue)query) // look in current node and right subtree
        {
            unsigned int i = 0;
            while (i<length(act_prop.list2) && rightBoundary(value(act_prop.list2, i))>(TPropertyValue) query)
            {
                appendValue(result, cargo(value(act_prop.list2, i)), Generous());
                ++i;
            }
            if (atEnd(it))
                break;

            next_prop = property(pm, targetVertex(it));
            if (next_prop.center <= act_prop.center)
            {
                goNext(it);
                if (atEnd(it))
                    break;
            }
            act_knot = targetVertex(it);
        }
        else
        {
            if ((TPropertyValue)query < act_prop.center) // look in current node and left subtree
            {
                unsigned int i = 0;
                while (i < length(act_prop.list1) && leftBoundary(value(act_prop.list1, i)) <= (TPropertyValue)query)
                {
                    appendValue(result, cargo(value(act_prop.list1, i)), Generous());
                    ++i;
                }
                if (atEnd(it))
                    break;

                next_prop = property(pm, targetVertex(it));
                if (next_prop.center >= act_prop.center)
                {
                    goNext(it);
                    if (atEnd(it))
                        break;
                }
                act_knot = targetVertex(it);
            }
            else  // look in current node only, as query is center
            {
                for (unsigned int i = 0; i < length(act_prop.list1); ++i)
                    appendValue(result, cargo(value(act_prop.list1, i)), Generous());
                break;
            }
        }
    }

}

template <typename TValue, typename TCargo, typename TValue2>
inline void
findIntervals(
        String<TCargo> & result,
        IntervalTree<TValue, TCargo> const & it,
        TValue2 query)
{
    findIntervals(result, it.g, it.pm, query);
}

template <typename TValue, typename TCargo, typename TValue2>
inline void
findIntervals(
        String<TCargo> & result,
        IntervalTree<TValue, TCargo> const & tree,
        TValue2 query_begin,
        TValue2 query_end)
{
    findIntervals(result, tree.g, tree.pm, query_begin, query_end);
}

template <typename TSpec, typename TPropertyMap, typename TValue, typename TCargo>
inline void
findIntervals(
        String<TCargo> & result,
        Graph<TSpec> const & g,
        TPropertyMap const & pm,
        TValue query_begin,
        TValue query_end)
{
    typedef Graph<TSpec> const TGraph;
    typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;

    resize(result, 0);

    // start at root
    TVertexDescriptor act_knot = 0;
    findIntervals(result, g, pm, act_knot, query_begin, query_end);
}

template <
    typename TSpec,
    typename TPropertyMap,
    typename TVertexDescriptor,
    typename TValue,
    typename TCargo>
inline void
findIntervals(
        String<TCargo> & result,
        Graph<TSpec> const & g,
        TPropertyMap const & pm,
        TVertexDescriptor & act_knot,
        TValue query_begin,
        TValue query_end)
{
    typedef typename Value<TPropertyMap>::Type TProperty;
    typedef typename Iterator<Graph<TSpec>, OutEdgeIterator>::Type TOutEdgeIterator;

    if (empty(g))
        return;

    TProperty act_prop = property(pm, act_knot);
    TProperty next_prop;

    while (true)
    {
        TOutEdgeIterator it(g, act_knot);
        act_prop = property(pm, act_knot);
        //
        if (act_prop.center < query_begin) // query interval is to the right of node center
        {
            unsigned int i = 0;
            while (i<length(act_prop.list2) && rightBoundary(value(act_prop.list2, i))> query_begin)
            {
                appendValue(result, cargo(value(act_prop.list2, i)), Generous());
                ++i;
            }
            if (atEnd(it))
                break;

            next_prop = property(pm, targetVertex(it));
            if (next_prop.center <= act_prop.center)
            {
                goNext(it);
                if (atEnd(it))
                    break;
            }
            act_knot = targetVertex(it);
        }
        else
        {
            if (query_end <= act_prop.center) // query interval is to the left of node center
            {
                unsigned int i = 0;
                while (i < length(act_prop.list1) && leftBoundary(value(act_prop.list1, i)) < query_end)
                {
                    appendValue(result, cargo(value(act_prop.list1, i)), Generous());
                    ++i;
                }

                if (atEnd(it))
                    break;

                next_prop = property(pm, targetVertex(it));
                if (next_prop.center >= act_prop.center)
                {
                    goNext(it);
                    if (atEnd(it))
                        break;
                }
                act_knot = targetVertex(it);
            }
            else
            { //node center is contained in query interval
                for (unsigned int i = 0; i < length(act_prop.list1); ++i)
                    appendValue(result, cargo(value(act_prop.list1, i)), Generous());

                while (!atEnd(it))
                {
                    TVertexDescriptor next_knot = targetVertex(it);
                    findIntervals(result, g, pm, next_knot, query_begin, query_end);
                    goNext(it);
                }
                break;

                //break; //dont break! continue in both subtrees!!
            }
        }
    }
}

/*!
 * @fn IntervalTree#findIntervalsExcludeTouching
 * @brief Find all intervals that contain the query point, exclude intervals that touch the query, i.e. where the query
 *        point equals the start or end point.
 *
 * @signature void findIntervalsExcludeTouching(result, intervalTree, query);
 * @signature void findIntervalsExcludeTouching(result, graph, propertyMap, query,);
 *
 * @param[out] result      The resulting string of cargos/ids of the intervals that contain the query point.  Should
 *                         be a string of TCargo. Types: String
 * @param[in] intervalTree An interval tree Types: IntervalTree
 * @param[in] graph        The directed graph that contains the topography of the interval tree.
 * @param[in] query        The TValue to query here.
 * @param[in] propertyMap  The property map containing the node properties of the interval tree
 */

template <typename TSpec, typename TPropertyMap, typename TValue, typename TCargo>
inline void
findIntervalsExcludeTouching(
        String<TCargo> & result,
        Graph<TSpec> const & g,
        TPropertyMap const & pm,
        TValue query)
{
    typedef Graph<TSpec> const TGraph;
    typedef typename Iterator<TGraph, OutEdgeIterator>::Type TOutEdgeIterator;
    typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
    typedef typename Value<TPropertyMap>::Type TProperty;

    resize(result, 0);
    if (empty(g))
        return;

    // start at root
    TVertexDescriptor act_knot = 0;
    TProperty act_prop = property(pm, act_knot);
    TProperty next_prop;

    while (true)
    {
        TOutEdgeIterator it(g, act_knot);
        act_prop = property(pm, act_knot);
        if ((TValue) act_prop.center < query) // look in current node and right subtree
        {
            int i = 0;
            while (i < (int)length(act_prop.list2) && (TValue)rightBoundary(value(act_prop.list2, i)) > query)
            {
                appendValue(result, cargo(value(act_prop.list2, i)), Generous());
                ++i;
            }
            if (atEnd(it))
                break;

            next_prop = property(pm, targetVertex(it));
            if (next_prop.center <= act_prop.center)
            {
                goNext(it);
                if (atEnd(it))
                    break;
            }
            act_knot = targetVertex(it);
        }
        else
        {
            if (query < (TValue) act_prop.center) // look in current node and left subtree
            {
                int i = 0;
                while (i < (int)length(act_prop.list1) && (TValue)leftBoundary(value(act_prop.list1, i)) < query)
                {
                    appendValue(result, cargo(value(act_prop.list1, i)), Generous());
                    ++i;
                }
                if (atEnd(it))
                    break;

                next_prop = property(pm, targetVertex(it));
                if (next_prop.center >= act_prop.center)
                {
                    goNext(it);
                    if (atEnd(it))
                        break;
                }
                act_knot = targetVertex(it);
            }
            else  // look in current node only
            {
                int i = 0;
                while (i < (int)length(act_prop.list1) && (TValue)leftBoundary(value(act_prop.list1, i)) < query)
                {
                    appendValue(result, cargo(value(act_prop.list1, i)), Generous());
                    ++i;
                }
                break;
            }
        }
    }

}

template <typename TValue, typename TCargo>
inline void
findIntervalsExcludeTouching(
        String<TCargo> & result,
        IntervalTree<TValue, TCargo> const & tree,
        TValue query)
{
    findIntervalsExcludeTouching(result, tree.g, tree.pm, query);
}

/*!
 * @fn IntervalTree#removeInterval
 * @brief Removes an interval from the interval tree.
 *
 * @signature bool removeInterval(intervalTree, iBegin, iEnd, iId);
 *
 * @param[in,out] intervalTree An interval tree Types: IntervalTree
 * @param[in]     iBegin       The begin position of the interval to be removed.
 * @param[in]     iEnd         The end position of the interval to be removed.
 * @param[in]     iId          The ID of the interval to be removed.
 *
 * @return bool <tt>true</tt> on success, <tt>false</tt> on failure.
 */

template <
    typename TSpec,
    typename TPropertyMap,
    typename TVertexDescriptor,
    typename TValue,
    typename TCargo>
inline bool
removeInterval(
    Graph<TSpec> & g,
    TPropertyMap & pm,
    TVertexDescriptor & act_knot,
    TValue i_begin,
    TValue i_end,
    TCargo i_id)
{
    SEQAN_CHECKPOINT

    typedef typename Value<TPropertyMap>::Type TProperty;
    typedef typename ListType<TProperty>::Type TList;
    typedef typename Iterator<TList, Standard>::Type TListIterator;
    typedef typename Iterator<Graph<TSpec>, OutEdgeIterator>::Type TOutEdgeIterator;

    if (empty(g))
        return false;

    TProperty act_prop = property(pm, act_knot);
    TProperty next_prop;

    while (true)
    {
        TOutEdgeIterator it(g, act_knot);
        act_prop = property(pm, act_knot);
        //
        if (act_prop.center < i_begin) // interval is to the right of node center
        {
            if (atEnd(it))
                break;

            next_prop = property(pm, targetVertex(it));
            if (next_prop.center <= act_prop.center)
            {
                goNext(it);
                if (atEnd(it))
                    break;
            }
            act_knot = targetVertex(it);
        }
        else
        {
            if (i_end < act_prop.center) // interval is to the left of node center
            {
                if (atEnd(it))
                    break;

                next_prop = property(pm, targetVertex(it));
                if (next_prop.center >= act_prop.center)
                {
                    goNext(it);
                    if (atEnd(it))
                        break;
                }
                act_knot = targetVertex(it);
            }
            else //node center is contained in interval, this is where we should find the interval to be removed
            {   // remove from list1
                TProperty & change_prop = property(pm, act_knot);
                bool foundInLeft = false;
                TListIterator list_it_keep = begin(change_prop.list1, Standard());
                TListIterator list_it = begin(change_prop.list1, Standard());
                while (list_it != end(change_prop.list1, Standard()))
                {
                    //std::cout << "Element: " << getLeftBoundary(*list_it) << ".." << getRightBoundary(*list_it) << " " << cargo(*list_it) << std::endl;
                    if (getLeftBoundary(*list_it) == i_begin && cargo(*list_it) == i_id)
                    {
                        foundInLeft = true;
                        ++list_it;
                        continue;
                    }
                    *list_it_keep = *list_it;
                    //std::cout << "Element keep: " << getLeftBoundary(*list_it_keep) << ".." << getRightBoundary(*list_it_keep) << " " << cargo(*list_it_keep) << std::endl;
                    ++list_it; ++list_it_keep;
                }

                bool foundInRight = false;
                list_it_keep = begin(change_prop.list2, Standard());
                list_it = begin(change_prop.list2, Standard());
                while (list_it != end(change_prop.list2, Standard()))
                {
                    //std::cout << "Element: " << getLeftBoundary(*list_it) << ".." << getRightBoundary(*list_it) << " " << cargo(*list_it) << std::endl;
                    if (getRightBoundary(*list_it) == i_end && cargo(*list_it) == i_id)
                    {
                        foundInRight = true;
                        ++list_it;
                        continue;
                    }
                    *list_it_keep = *list_it;
                    //std::cout << "Element keep: " << getLeftBoundary(*list_it_keep) << ".." << getRightBoundary(*list_it_keep) << " " << cargo(*list_it_keep) << std::endl;
                    ++list_it; ++list_it_keep;

                }


                // TODO: if node is empty and does not have child nodes --> remove node.
                // keeping these empty leaf nodes just takes space unnecessarily
                if (foundInRight && foundInLeft)
                {
                    resize(change_prop.list2, length(change_prop.list2) - 1);
                    resize(change_prop.list1, length(change_prop.list1) - 1);
                    return true;
                }

            }
        }
    }
    return false;
}

template <
    typename TSpec,
    typename TPropertyMap,
    typename TValue,
    typename TCargo>
inline bool
removeInterval(
    Graph<TSpec> & g,
    TPropertyMap & pm,
    TValue i_begin,
    TValue i_end,
    TCargo i_id)
{
    typedef typename VertexDescriptor<Graph<TSpec> >::Type TVertexDescriptor;

    // start looking at root
    TVertexDescriptor act_knot = 0;
    return removeInterval(g, pm, act_knot, i_begin, i_end, i_id);
}

template <typename TValue, typename TCargo>
inline bool
removeInterval(
    IntervalTree<TValue, TCargo> & tree,
    TValue i_begin,
    TValue i_end,
    TCargo i_id)
{
    return removeInterval(tree.g, tree.pm, i_begin, i_end, i_id);
    // we do not decrease the interval_counter of tree, as it would mix up interval IDs
}

/////////////////// Metafunctions ///////////////////////

template <typename TValue, typename TCargo>
struct Value<IntervalTree<TValue, TCargo> >
{
    typedef TValue Type;
};


template <typename TValue, typename TCargo>
struct Cargo<IntervalTree<TValue, TCargo> >
{
    typedef TCargo Type;
};

}  // namespace SEQAN_NAMESPACE_MAIN

#endif  //#ifndef SEQAN_MISC_INTERVAL_TREE_H
