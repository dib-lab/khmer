// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2013, Knut Reinert, FU Berlin
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


namespace SEQAN_NAMESPACE_MAIN
{
//////////////////////////////////////////////////////////////////////////////
// Graph - Interval Tree Types
//////////////////////////////////////////////////////////////////////////////

///---------------------------------------------------------------///

//////////////////// Interval and ID type ///////////////////
/**
.Class.IntervalAndCargo:
..cat:Miscellaneous
..summary:A simple record type that stores an interval and a cargo value.
..signature:IntervalAndCargo<TValue, TCargo>
..param.TValue:The value type, that is the type of the interval borders.
...default:int.
...metafunction:Metafunction.Value
..param.TCargo:The cargo type.
...default:int.
...metafunction:Metafunction.Cargo
..include:seqan/misc/misc_interval_tree.h
*/
template<typename TValue = int, typename TCargo = int>
class IntervalAndCargo
{
public:
    /**
.Memvar.IntervalAndCargo#i1:
..class:Class.PointAndCargo
..summary:The first element in the interval of type i1.
     */
	TValue i1;

    /**
.Memvar.IntervalAndCargo#i2:
..class:Class.PointAndCargo
..summary:The last element in the interval of type i2.
     */
	TValue i2;

    /**
.Memvar.IntervalAndCargo#cargo:
..class:Class.PointAndCargo
..summary:The stored cargo of type TCargo.
     */
	TCargo cargo;

    /**
.Memfunc.IntervalAndCargo#IntervalAndCargo:
..signature:IntervalAndCargo()
     */
    IntervalAndCargo()
    {
SEQAN_CHECKPOINT
    }

    /**
.Memfunc.IntervalAndCargo#IntervalAndCargo:
..class:Class.IntervalAndCargo
..summary:Constructor.
..signature:IntervalAndCargo(i1, i2, cargo)
..param.i1:The first element in the interval, of type TValue.
..param.i2:The last element in the interval of type TValue.
..param.cargo:The cargo value of type TCargo.
     */
	IntervalAndCargo(TValue i1, TValue i2, TCargo cargo):
		i1(i1), i2(i2), cargo(cargo)
	{
SEQAN_CHECKPOINT
	}
};



/////////////////////// Point and ID type ////////////////
/**
.Class.PointAndCargo:
..cat:Miscellaneous
..summary:Simple record class storing a point (one-value interval) and a cargo.
..signature:PointAndCargo<TValue, TCargo>
..param.TValue:
...default:int.
...metafunction:Metafunction.Value
..param.TCargo:
...default:int.
...metafunction:Metafunction.Value
..include:seqan/misc/misc_interval_tree.h
*/
template<typename TValue=int, typename TCargo=int>
class PointAndCargo {
public:
    /**
.Memvar.PointAndCargo#point:
..class:Class.PointAndCargo
..summary:The stored point of type TValue.
     */
	TValue point;

    /**
.Memvar.PointAndCargo#cargo:
..class:Class.PointAndCargo
..summary:The stored cargo of type TCargo.
     */
	TCargo cargo;

    /**
.Memfunc.PointAndCargo#PointAndCargo
..class:Class.PointAndCargo
..signature:PointAndCargo(point, cargo)
    */
	PointAndCargo() {
SEQAN_CHECKPOINT
	}

    /**
.Memfunc.PointAndCargo#PointAndCargo
..class:Class.PointAndCargo
..summary:Constructor.
..signature:PointAndCargo(point, cargo)
..param.point:
...summary:The point to store of type TValue.
..param.cargo:
...summary:The cargo to store of type TCargo.
    */
	PointAndCargo(TValue point, TCargo cargo):
		point(point), cargo(cargo)
	{
SEQAN_CHECKPOINT
	}
};

///////////////////////////////////////////////////////////////////////////
/////////////////////////// IntervalTreeNode	///////////////////////////

/**
.Tag.IntervalTree Node Types
..summary:Tags to select the node type for @Class.IntervalTree@.
..cat:Miscellaneous
..see:Class.IntervalTree

..tag.StorePointsOnly:The tree nodes store points.
..include:seqan/misc/misc_interval_tree.h
*/
struct StorePointsOnly {};


///..tag.StoreIntervals:The tree nodes store intervals.
struct StoreIntervals {};


/**
.Class.IntervalTreeNode:
..cat:Miscellaneous
..summary:Element of @Class.IntervalTree@.
..signature:IntervalTreeNode<TInterval, TSpec>
..param.TInterval:The type of interval to store.
..param.TSpec:The type of interval to store.
...default:StorePointsOnly.
...metafunction:Metafunction.Spec
..include: seqan/misc/misc_interval_tree.h

.Memvar.IntervalTreeNode#center:
..class:Class.IntervalTreeNode
..summary:The center of the interval of type TValue.

.Memvar.IntervalTreeNode#list1
..class:Class.IntervalTreeNode
..summary:Sorted list of pointers to intervals, sorted in ascending according to their left boundary points.

.Memvar.IntervalTreeNode#list2
..class:Class.IntervalTreeNode
..summary:Sorted list of pointers to intervals, sorted in descending according to their right boundary points.
 */
template<typename TInterval, typename TSpec=StorePointsOnly>
class IntervalTreeNode;


/**
.Spec.Interval Tree Node
..cat:Miscellaneous
..general:Class.IntervalTreeNode
..summary:An Interval Tree Node that stores intervals explicitely in each node.
..signature:IntervalTreeNode<TInterval, StoreIntervals>
..param.TInterval:The interval type to store in the node.
..include:seqan/misc/misc_interval_tree.h
 */
template<typename TInterval>
class IntervalTreeNode<TInterval, StoreIntervals> {
public:
	typedef typename Value<TInterval>::Type TValue;

	TValue center;
	String<TInterval> list1;
	String<TInterval> list2;

    IntervalTreeNode() : center()
    {}
};


/**
.Spec.Points Only Tree Node
..cat:Miscellaneous
..general:Class.IntervalTreeNode
..summary:Spec for IntervalTreeNode that stores only the relevant point in each node meaning the endpoint of the interval in the list sorted by endpoints (list2) and only the beginpoint of the interval in the list sorted by beginpoints (list1).
..signature:IntervalTreeNode<TInterval, StorePointsOnly>
..param.TInterval:The interval type to store in the node.
..include:seqan/misc/misc_interval_tree.h
 */
template<typename TInterval>
class IntervalTreeNode<TInterval, StorePointsOnly> {
public:
	typedef typename Cargo<TInterval>::Type TCargo;
	typedef typename Value<TInterval>::Type TValue;

	TValue center;
	String<PointAndCargo<TValue,TCargo> > list1;
	String<PointAndCargo<TValue,TCargo> > list2;

    /**
.Memfunc.IntervalTreeNode#IntervalTreeNode:
..class:Class.IntervalTreeNode
..summary:Default constructor.
..signature:IntervalTreeNode()
     */
    IntervalTreeNode() : center()
    {
SEQAN_CHECKPOINT
    }

	IntervalTreeNode(IntervalTreeNode const & other):
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

/**
.Class.IntervalTree:
..cat:Miscellaneous
..summary:A datastructure that efficiently stores intervals.
..signature:IntervalTree<TValue, TCargo>
..param.TValue:The value type.
...default:int
..param.TCargo:The cargo/id type.
...default:int
...remarks:If the intervals are not associated with cargos/IDs, they will be numbered consecutively.
..example:The following example creates an integer interval tree with string keys.
This tree is queried for keys of intervals that overlap the interval [550,900).
...file:demos/misc/interval_tree_example.cpp
...text:The resulting keys are:
...output:
gene
exon2
coding2
..include:seqan/misc/misc_interval_tree.h
*/
template<typename TValue=int, typename TCargo=unsigned int>
class IntervalTree
{
public:
	typedef Graph<Directed<void,WithoutEdgeId> > TGraph;
	typedef IntervalAndCargo<TValue,TCargo> TInterval;
	typedef IntervalTreeNode<TInterval> TNode;
	typedef String<TNode> TPropertyMap;

	TGraph g;
	TPropertyMap pm;
	size_t interval_counter;
	
/**
.Memfunc.IntervalTree#IntervalTree:
..class:Class.IntervalTree
..summary:Constructor
..signature:IntervalTree()
..signature:IntervalTree(String<TInterval> intervals)
..signature:IntervalTree(String<TInterval> intervals, TValue center)
..signature:IntervalTree(String<TInterval> intervals, Tag<TTagSpec> const tag)
..signature:IntervalTree(intervalBegins, intervalEnds, len)
..signature:IntervalTree(intervalBegins, intervalEnds, intervalCargos, len)
..param.intervals:Container of intervals.
...type:Spec.Alloc String
...remarks:A string of $IntervalAndCargo<TValue, TCargo>$ objects, see @Class.IntervalAndCargo@.
..param.intervalBegins:Iterator pointing to begin position of first interval.
..param.intervalEnds:Iterator pointing to end position of first interval.
..param.intervalCargos:Iterator pointing to cargos/ids for intervals.
..param.len:Number of intervals to store in tree.
..param.tag:Tag for tree construction method; @Tag.IntervalTree Centers.tag.ComputeCenter@
...default:@Tag.IntervalTree Centers.tag.ComputeCenter@
..remarks:center of root node is computed by _calcIntervalTreeRootCenter
*/
	
	IntervalTree()
	{
SEQAN_CHECKPOINT
		interval_counter = 0;
	}
	
	template<typename TIterator,typename TCargoIterator>
	IntervalTree(TIterator interval_begins,
				 TIterator interval_ends, 
				 TCargoIterator interval_cargos, 
				 size_t len)	
	{
SEQAN_CHECKPOINT
		String<TInterval> intervals;
		resize(intervals,len);
		size_t i = 0;
		while(i<len)
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

	template<typename TIterator>
	IntervalTree(TIterator interval_begins,
				 TIterator interval_ends,
				 size_t len)
	{
SEQAN_CHECKPOINT
		String<TInterval> intervals;
		resize(intervals,len);
		size_t i = 0;
		while(i<len)
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
		createIntervalTree(g,pm,intervals,tag);
	}

	IntervalTree(String<TInterval> intervals, TValue center)	
	{
SEQAN_CHECKPOINT
		interval_counter = length(intervals);
		createIntervalTree(g,pm,intervals,center);
	}
};



///////Specs for the way interval centers are determined
/**
.Tag.IntervalTree Centers
..cat:Miscellaneous
..summary:Tag to select a specific way to compute the center of an interval tree node.
..see:Class.IntervalTree
..include:seqan/misc/misc_interval_tree.h
 */


/**
..tag.ComputeCenter
...summary:For intervals that are more or less uniformly distributed in the value range, using the ComputeCenter tag may result in a more balanced tree compared to using a random approach.
...signature:ComputeCenter
...remarks:center = minbegin + (maxend-minbegin)/2
 */
//template <typename TSpec = SpecPointAndCargo>
struct TagComputeCenter_;
typedef Tag<TagComputeCenter_> const ComputeCenter;


///////////////////////////////////////////////////////////////////////////
///////////////////// IntervalAndCargo functions //////////////////////////
///////////////////////////////////////////////////////////////////////////




/**
.Function.leftBoundary
..cat:Miscellaneous
..summary:Access to the left boundary.
..signature:leftBoundary(interval)
..class:Class.IntervalAndCargo
..param.interval:The interval to return the left boundary for.
...type:Class.IntervalAndCargo
..returns:The reference to the left boundary of the interval of type TValue&.
..see:Function.getLeftBoundary
..see:Function.rightBoundary
..see:Function.getRightBoundary
*/

template<typename TValue, typename TCargo>
TValue &
leftBoundary(IntervalAndCargo<TValue, TCargo> & interval)
{
	return interval.i1;
}

template<typename TValue, typename TCargo>
TValue const &
leftBoundary(IntervalAndCargo<TValue, TCargo> const & interval)
{
	return interval.i1;
}


/**
.Function.rightBoundary
..cat:Miscellaneous
..summary:Access to the right boundary.
..signature:leftBoundary(interval)
..class:Class.IntervalAndCargo
..param.interval:The interval to return the left boundary for.
...type:Class.IntervalAndCargo
..returns:The reference to the right boundary of the interval of type TValue&.
..see:Function.getRightBoundary
..see:Function.leftBoundary
..see:Function.getLeftBoundary
*/

template<typename TValue, typename TCargo>
TValue &
rightBoundary(IntervalAndCargo<TValue, TCargo> & interval)
{
	return interval.i2;
}

template<typename TValue, typename TCargo>
TValue const &
rightBoundary(IntervalAndCargo<TValue, TCargo> const & interval)
{
	return interval.i2;
}


/**
.Function.getLeftBoundary
..cat:Miscellaneous
..summary:Get method for the left boundary.
..signature:leftBoundary(interval)
..class:Class.IntervalAndCargo
..param.interval:The interval to return the left boundary for.
...type:Class.IntervalAndCargo
..returns:The copy of the left boundary of the interval of type TValue.
..see:Function.leftBoundary
..see:Function.getRightBoundary
..see:Function.rightBoundary
*/

template<typename TValue, typename TCargo>
TValue
getLeftBoundary(IntervalAndCargo<TValue, TCargo> const & interval)
{
	return interval.i1;
}


/**
.Function.getRightBoundary
..cat:Miscellaneous
..summary:Get method for the right boundary.
..signature:leftBoundary(interval)
..class:Class.IntervalAndCargo
..param.interval:The interval to return the right boundary for.
...type:Class.IntervalAndCargo
..returns:The copy of the right boundary of the interval of type TValue.
..see:Function.rightBoundary
..see:Function.getLeftBoundary
..see:Function.leftBoundary
*/

template<typename TValue, typename TCargo>
TValue
getRightBoundary(IntervalAndCargo<TValue, TCargo> const & interval)
{
	return interval.i2;
}


/**
.Function.cargo
..signature:cargo(me)
..class:Class.IntervalAndCargo
..param.me:
...type:Class.IntervalAndCargo
..see:Function.getCargo
*/

template<typename TValue, typename TCargo>
TCargo const &
cargo(IntervalAndCargo<TValue, TCargo> const & interval)
{
	return interval.cargo;
}

template<typename TValue, typename TCargo>
TCargo &
cargo(IntervalAndCargo<TValue, TCargo> & interval)
{
	return interval.cargo;
}

/**
.Function.getCargo
..signature:getCargo(me)
..class:Class.IntervalAndCargo
..param.me:
...type:Class.IntervalAndCargo
..see:Function.cargo
*/

template<typename TValue, typename TCargo>
TCargo
getCargo(IntervalAndCargo<TValue,TCargo> const & interval)
{
	return interval.cargo;
}


/////////////////// Metafunctions //////////////////////
    
///.Metafunction.Value.param.T.type:Class.IntervalAndCargo
///.Metafunction.Value.class:Class.IntervalAndCargo
template<typename TValue,typename TCargo>
struct Value<IntervalAndCargo<TValue,TCargo> >
{
	typedef TValue Type;
};


///.Metafunction.Cargo.param.T.type:Class.IntervalAndCargo
///.Metafunction.Cargo.class:Class.IntervalAndCargo
template<typename TValue,typename TCargo>
struct Cargo<IntervalAndCargo<TValue, TCargo> >
{
	typedef TCargo Type;
};


///////////////////////////////////////////////////////////////////////////
///////////////////// PointAndCargo functions /////////////////////////////
///////////////////////////////////////////////////////////////////////////

/**
.Function.leftBoundary
..signature:leftBoundary(point)
..class:Class.PointAndCargo
..param.point.type:Class.PointAndCargo
 */

template<typename TValue, typename TCargo>
TValue const &
leftBoundary(PointAndCargo<TValue, TCargo> const & point)
{
	return point.point;
}

template<typename TValue, typename TCargo>
TValue &
leftBoundary(PointAndCargo<TValue, TCargo> & point)
{
	return point.point;
}


/**
.Function.rightBoundary
..signature:rightBoundary(point)
..class:Class.PointAndCargo
..param.point.type:Class.PointAndCargo
 */

template<typename TValue, typename TCargo>
TValue const &
rightBoundary(PointAndCargo<TValue, TCargo> const & point)
{
	return point.point;
}

template<typename TValue, typename TCargo>
TValue &
rightBoundary(PointAndCargo<TValue, TCargo> & point)
{
	return point.point;
}


/**
.Function.getLeftBoundary
..signature:getLeftBoundary(point)
..class:Class.PointAndCargo
..param.point.type:Class.PointAndCargo
 */

template<typename TValue, typename TCargo>
TValue
getLeftBoundary(PointAndCargo<TValue, TCargo> const & point)
{
    return point.point;
}


/**
.Function.getRightBoundary
..signature:getRightBoundary(point)
..class:Class.PointAndCargo
..param.point.type:Class.PointAndCargo
 */

template<typename TValue, typename TCargo>
TValue
getRightBoundary(PointAndCargo<TValue,TCargo> const & point)
{
	return point.point;
}


/**
.Function.cargo
..signature:cargo(point)
..class:Class.PointAndCargo
..param.point.type:Class.PointAndCargo
 */

template<typename TValue, typename TCargo>
TCargo const &
cargo(PointAndCargo<TValue,TCargo> const & point)
{
	return point.cargo;
}

template<typename TValue, typename TCargo>
TCargo &
cargo(PointAndCargo<TValue,TCargo> & point)
{
	return point.cargo;
}


/**
.Function.cargo
..signature:getCargo(point)
..class:Class.IntervalAndCargo
..param.point:
...type:Class.IntervalAndCargo
..see:Function.getCargo
*/

template<typename TValue, typename TCargo>
TCargo
getCargo(PointAndCargo<TValue,TCargo> const & point)
{
	return point.cargo;
}

////////////////// Metafunctions //////////////////
///.Metafunction.Value.param.T.type:Class.PointAndCargo
///.Metafunction.Value.class:Class.PointAndCargo
template<typename TValue,typename TCargo>
struct Value<PointAndCargo<TValue,TCargo> >
{
	typedef TValue Type;
};


///.Metafunction.Cargo.param.T.type:Class.PointAndCargo
///.Metafunction.Cargo.class:Class.PointAndCargo
template<typename TValue,typename TCargo>
struct Cargo<PointAndCargo<TValue,TCargo> >
{
	typedef TCargo Type;
};


//// Comparators
template <typename TPair>
bool _less_compI1_ITree(TPair const & p1, TPair const & p2)
{
    return (leftBoundary(p1) < leftBoundary(p2));
}


template <typename TPair>
bool _greater_compI2_ITree(TPair const & p1, TPair const & p2)
{
    return (rightBoundary(p1) > rightBoundary(p2));
}



///////////////////////////////////////////////////////////////////////////
///////////////////// IntervalTreeNode functions //////////////////////////
///////////////////////////////////////////////////////////////////////////





// internal set node functions
template<typename TValue, typename TInterval>
void
_setIntervalTreeNode(IntervalTreeNode<TInterval, StoreIntervals> & knot, TValue center, TInterval const & interval)
{
	knot.center = center;
	appendValue(knot.list1, interval);
	appendValue(knot.list2, interval);
}

// append intervals to lists in node knot
template<typename TInterval>
void
_appendIntervalTreeNodeLists(IntervalTreeNode<TInterval, StoreIntervals> & knot,TInterval const & interval)
{
	appendValue(knot.list1, interval);
	appendValue(knot.list2, interval);
}


//internal set node functions
template<typename TValue, typename TInterval>
void
_setIntervalTreeNode(IntervalTreeNode<TInterval, StorePointsOnly> & knot, TValue center, TInterval const & interval)
{
	knot.center = center;
	appendValue(knot.list1,PointAndCargo<TValue,typename Cargo<TInterval>::Type>(leftBoundary(interval),cargo(interval)));
	appendValue(knot.list2,PointAndCargo<TValue,typename Cargo<TInterval>::Type>(rightBoundary(interval),cargo(interval)));
}


template<typename TInterval>
void
_appendIntervalTreeNodeLists(IntervalTreeNode<TInterval, StorePointsOnly> & knot, TInterval const & interval)
{
	appendValue(knot.list1,PointAndCargo<typename Value<TInterval>::Type,typename Cargo<TInterval>::Type>(leftBoundary(interval),cargo(interval)));
	appendValue(knot.list2,PointAndCargo<typename Value<TInterval>::Type,typename Cargo<TInterval>::Type>(rightBoundary(interval),cargo(interval)));
	

}

/////////////////// Metafunctions ///////////////////////
///.Metafunction.Value.param.T.type:Class.IntervalTreeNode
///.Metafunction.Value.class:Class.IntervalTreeNode
template<typename TInterval, typename TSpec>
struct Value<IntervalTreeNode<TInterval,TSpec> >
{
	typedef typename Value<TInterval>::Type Type;
};


///.Metafunction.Cargo.param.T.type:Class.IntervalTreeNode
///.Metafunction.Cargo.class:Class.IntervalTreeNode
template<typename TInterval, typename TSpec>
struct Cargo<IntervalTreeNode<TInterval,TSpec> >
{
	typedef typename Cargo<TInterval>::Type Type;
};


/**
.Metafunction.ListType:
..cat:Miscellaneous
..signature:ListType<T>::ListType
..class.Class.IntervalTreeNode
..summary:Type of lists in tree nodes.
..param.T:The type to retrieve the list type for.
..returns:Returns the type of the lists in @Class.IntervalTreeNode@ objects.
..include:seqan/misc/misc_interval_tree.h
 */
template<typename T>
struct ListType;


///.Metafunction.ListType.param.T.type:Class.IntervalTreeNode
///.Metafunction.ListType.class:Class.IntervalTreeNode
template<typename TInterval>
struct ListType<IntervalTreeNode<TInterval,StorePointsOnly> >
{
	typedef String<PointAndCargo<typename Value<TInterval>::Type,typename Cargo<TInterval>::Type> > Type;

};


///.Metafunction.ListType.param.T.type:Class.IntervalTreeNode
///.Metafunction.ListType.class:Class.IntervalTreeNode
template<typename TInterval>
struct ListType<IntervalTreeNode<TInterval,StoreIntervals> >
{
	typedef String<IntervalAndCargo<typename Value<TInterval>::Type,typename Cargo<TInterval>::Type> > Type;

};



///////////////////////////////////////////////////////////////////////////
/////////////////////// IntervalTree functions ////////////////////////////
///////////////////////////////////////////////////////////////////////////

/**
.Function.createIntervalTree
..summary:Create an interval tree.
..cat:Miscellaneous
..class:Class.IntervalTree
..signature:createIntervalTree(intervalTree, intervals [, tag])
..signature:createIntervalTree(g, pm, intervals [, tag])
..signature:createIntervalTree(g, pm, intervals, center [, tag]])
..param.intervalTree:An interval tree
...type:Class.IntervalTree
..param.g:DirectedGraph to create interval tree in.
...type:Class.Graph
..param.pm:Property map to use for the created interval tree.
..param.intervals:Container of intervals.
...type:Spec.Alloc String
...remarks:A string of $IntervalAndCargo<TValue, TCargo>$ objects, see @Class.IntervalAndCargo@.
..param.tag:Tag for tree construction method; @Tag.IntervalTree Centers.tag.ComputeCenter@
...default:@Tag.IntervalTree Centers.tag.ComputeCenter@
..remarks:center of root node is computed by _calcIntervalTreeRootCenter
..include:seqan/misc/misc_interval_tree.h
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

	reserve(g.data_vertex,length(intervals));
	reserve(pm,length(intervals));

	TVertexDescriptor root = addVertex(g);
	resizeVertexMap(g,pm);

    if (length(intervals) > 0u) {
        TValue center =	_calcIntervalTreeRootCenter(intervals);

        std::sort(begin(intervals, Standard()),end(intervals, Standard()),_less_compI1_ITree<TInterval>);

        String<TInterval *> interval_pointers;
        // interval tree stores pointers to intervals, not original intervals
        _makePointerInterval(intervals, interval_pointers);

        _createIntervalTree(g,pm,interval_pointers,root,(TValue)0.0,center,length(intervals),tag);
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
	
	reserve(g.data_vertex,length(intervals));
	reserve(pm,length(intervals));
	
	TVertexDescriptor root = addVertex(g);
	resizeVertexMap(g,pm);
	
    if (length(intervals) > 0u) {
	    TInterval a;
	    typename Iterator<TIntervals, Standard>::Type begin_ = begin(intervals, Standard());
	    typename Iterator<TIntervals, Standard>::Type end_ = end(intervals, Standard());
	    std::sort(begin_, end_ ,_less_compI1_ITree<TInterval>);

	    String<TInterval const*> interval_pointers;
	    _makePointerInterval(intervals,interval_pointers);

	    if(length(intervals)==1) // if there is just one interval ->  center = center of this interval
		    center = (rightBoundary(intervals[0])-leftBoundary(intervals[0]))/(TValue)2.0;

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


template <typename TValue,typename TCargo, typename TIntervals, typename TSpec>
inline void
createIntervalTree(IntervalTree<TValue,TCargo> & it,
				   TIntervals & intervals,
				   Tag<TSpec> const tag)
{
	it.interval_counter = length(intervals);
	createIntervalTree(it.g, it.pm, intervals, tag);
}

template <typename TValue,typename TCargo, typename TIntervals>
inline void 
createIntervalTree(IntervalTree<TValue, TCargo> & it,
				   TIntervals & intervals)
{
	createIntervalTree(it, intervals, ComputeCenter());
}


//////////////////////////////////////////////////////////////////////////////
//remembers minimum and maximum of point values in intervals and sets the center
//of each node to min+(max-min)/2
template<typename TGraph, typename TPropertyMap, typename TIntervalPointer, typename TValue>
inline void
_createIntervalTree(TGraph & g, TPropertyMap & pm, 
				   String<TIntervalPointer*> & intervals,
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
	typedef String<TIntervalPointer*> TIntervalPointers;

	// one list of interval pointers for the intervals to the left of center
	TIntervalPointers S_left;
	// one list of interval pointers for the intervals to the right of center
	TIntervalPointers S_right;

	TValue min1 = maxValue<TValue>();
	TValue min2 = maxValue<TValue>();
	TValue max1 = minValue<TValue>();
	TValue max2 = minValue<TValue>();

	value(pm,knot).center = center;
	
 
	typedef typename Iterator<TIntervalPointers,Standard>::Type TIntervalIterator;
	TIntervalIterator it = begin(intervals,Standard());
	TIntervalIterator it_end = end(intervals,Standard());
	
	// walk through intervals
	while(it != it_end)
	{
		// interval belongs to the left list
		if((**it).i2<=center)
		{
			appendValue(S_left,*it, Generous());
			 //remember right most and left most point in left list
			if((**it).i2 > max1)
				max1 = (**it).i2;
			if((**it).i1 < min1)
				min1 = (**it).i1;
		}
		else
		{
			// interval belongs to the right list
			if((**it).i1>center)
			{
				appendValue(S_right,(*it), Generous());
				 //remember right most and left most point in right list
				if((**it).i2 > max2)
					max2 = (**it).i2;
				if ((**it).i1 < min2)
					min2 = (**it).i1;
			}
			else // interval belongs to this node
			{
				_appendIntervalTreeNodeLists(value(pm,knot), **it);
			}
		}
        ++it;
	}

//	std::sort(begin(value(pm,knot).list1),end(value(pm,knot).list1),_less_compI1_ITree<typename Value<TList>::Type>);
	std::sort(begin(value(pm,knot).list2),end(value(pm,knot).list2),_greater_compI2_ITree<typename Value<TList>::Type>);

	// build subtree to the left
	if(!empty(S_left))
	{
		TVertexDescriptor vd = addVertex(g);
		resize(pm, vd+1); 
		addEdge(g,knot,vd);
		_createIntervalTree(g,pm,S_left,vd,center,min1+(max1-min1)/2,length(S_left),tag);
	}
	// build subtree to the right
	if(!empty(S_right))
	{
		TVertexDescriptor vd = addVertex(g);
		resize(pm, vd+1); 
		addEdge(g,knot,vd);
		_createIntervalTree(g,pm,S_right,vd,center,min2+(max2-min2)/2,length(S_right),tag);
	}
}




//////////////////////////////////////////////////////////////////////////////
//createIntervalTree for all specs except CompCenter, the center value of each 
//node is determined by functions _calcIntervalTreeNodeCenterLeft and 
//_calcIntervalTreeNodeCenterRight
template<typename TGraph, typename TPropertyMap, typename TSpec, typename TInterval, typename TValue>
inline void
_createIntervalTree(TGraph & g, TPropertyMap & pm, 
				   String<TInterval*> & intervals, 
				   typename VertexDescriptor<TGraph>::Type & knot, 
				   TValue last_center, TValue center, 
				   typename VertexDescriptor<TGraph>::Type len,
				   Tag<TSpec> const tag)
{
SEQAN_CHECKPOINT
	// Rekursionsanker
	if(len==1){
		_setIntervalTreeNode(value(pm,knot),center,*value(intervals,0));
		return;
	}

	typedef typename Value<TPropertyMap>::Type TNode;
	typedef typename ListType<TNode>::Type TList;
	typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
	typedef String<TInterval*> TIntervalPointers;
	
	// one list of interval pointers for the intervals to the left of center
	TIntervalPointers S_left;
	// one list of interval pointers for the intervals to the right of center
	TIntervalPointers S_right;
		
	value(pm,knot).center = center;
	
	typedef typename Iterator<TIntervalPointers,Standard>::Type TIntervalIterator;
	TIntervalIterator it = begin(intervals,Standard());
	TIntervalIterator it_end = end(intervals,Standard());
	
	// walk through intervals
	while(it != it_end)
	{
		// interval belongs to the left list
		if((**it).i2<=center)
		{
			appendValue(S_left,*it, Generous());
		}
		else
		{	// interval belongs to the right list
			if((**it).i1>center)
			{
				appendValue(S_right,(*it), Generous());
			}
			else
			{
				// interval belongs to the current node
				_appendIntervalTreeNodeLists(value(pm,knot),**it);
			}
		}
        ++it;
	}

//	std::sort(begin(value(pm,knot).list1),end(value(pm,knot).list1),_less_compI1_ITree<typename Value<TList>::Type>);
	std::sort(begin(value(pm,knot).list2),end(value(pm,knot).list2),_greater_compI2_ITree<typename Value<TList>::Type>);

	// build subtree to the left
	if(!empty(S_left))
	{
		TVertexDescriptor vd = addVertex(g);
		resizeVertexMap(g,pm);
		addEdge(g,knot,vd);
		TValue next_center = _calcIntervalTreeNodeCenterLeft(S_left,last_center,center,tag);
		_createIntervalTree(g,pm,S_left,vd,center,next_center,length(S_left),tag);
	}
	// build subtree to the right
	if(!empty(S_right))
	{
		TVertexDescriptor vd = addVertex(g);
		resizeVertexMap(g,pm);
		addEdge(g,knot,vd);
		TValue next_center = _calcIntervalTreeNodeCenterRight(S_right,last_center,center,tag);
		_createIntervalTree(g,pm,S_right,vd,center,next_center,length(S_right),tag);
	}
}


// fill the container interval_pointers with pointers to the corresponding objects in intervals.
// this is done to avoid copying and passing the whole IntervalAndCargo objects during interval tree construction
template<typename TIntervals, typename TIntervalPointers>
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
template<typename TIntervals>
typename Value<typename Value<TIntervals>::Type>::Type
_calcIntervalTreeRootCenter(TIntervals & intervals)
{
    SEQAN_CHECKPOINT;

    SEQAN_ASSERT_GT(length(intervals), 0u);
	
	typedef typename Value<typename Value<TIntervals>::Type>::Type TValue;
	typedef typename Iterator<TIntervals,Standard>::Type TIntervalIterator;

	TIntervalIterator it = begin(intervals);
	TIntervalIterator it_end = end(intervals);

	TValue min = maxValue<TValue>();
	TValue max = minValue<TValue>();

    // get min and max
	while(it != it_end)
	{
		if(leftBoundary(*it)<min) min = leftBoundary(*it);
		if(rightBoundary(*it)>max) max = rightBoundary(*it);
	  SEQAN_ASSERT_LEQ(min, max);
		++it;
	}

	SEQAN_ASSERT_LEQ(min, max);
	
    // return middle between max and min
	return (min+(max-min)/(TValue)2.0);

}



/**
.Function.addInterval
..summary:Adds an interval to an interval tree.
..cat:Miscellaneous
..class:Class.IntervalTree
..signature:addInterval(intervalTree, interval)
..signature:addInterval(intervalTree, begin, end)
..signature:addInterval(intervalTree, begin, end, cargo)
..signature:addInterval(graph, propertyMap, interval)
..param.intervalTree:The interval tree to add the interval to.
...type:Class.IntervalTree
..param.interval:The interval to be added to the interval tree.
..param.begin:Begin position of interval of type TValue.
..param.end:End position of interval of type TValue.
..param.cargo:Cargo to attach to the interval.
...type:Class.IntervalAndCargo
..param.graph:The directed graph that contains the topography of the interval tree.
..param.propertyMap:The property map containing the node properties of the interval tree.
..include:seqan/misc/misc_interval_tree.h
*/

template<typename TGraph, typename TPropertyMap, typename TInterval>
void
addInterval(TGraph & g, TPropertyMap & pm, TInterval interval)
{
SEQAN_CHECKPOINT

	typedef typename Iterator<TGraph, OutEdgeIterator>::Type TOutEdgeIterator;
	typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
	typedef typename Value<TPropertyMap>::Type TProperty;
	typedef typename Value<TInterval>::Type TValue;
	typedef typename ListType<TProperty>::Type TList;
	

	if(empty(pm))
	{
		TVertexDescriptor vd = addVertex(g);
		resizeVertexMap(g,pm);
		_setIntervalTreeNode(property(pm,vd),(rightBoundary(interval)+leftBoundary(interval))/2,interval);
		return;
		
	}
	// start at root
	TVertexDescriptor act_knot = 0;
	TProperty act_prop = property(pm,act_knot);
	TProperty next_prop;
	
    // look for the right node to add interval to
	while(true)
	{
		TOutEdgeIterator it(g, act_knot);
		act_prop = property(pm,act_knot);
		if(act_prop.center < leftBoundary(interval)) // interval to the left of current node?
		{
			if(atEnd(it)){
				TVertexDescriptor vd = addVertex(g);
				resizeVertexMap(g,pm);
				addEdge(g,act_knot,vd);
				_setIntervalTreeNode(property(pm,vd),(rightBoundary(interval)+leftBoundary(interval))/(TValue)2.0,interval);
				break;
			}
			else{
				next_prop = property(pm,targetVertex(it));
				if(next_prop.center <= act_prop.center)
				{
					goNext(it);
					if(atEnd(it)){
						TVertexDescriptor vd = addVertex(g);
						resizeVertexMap(g,pm);
						addEdge(g,act_knot,vd);
						_setIntervalTreeNode(property(pm,vd),(rightBoundary(interval)+leftBoundary(interval))/(TValue)2.0,interval);
						break;
					}
				}
			}
			act_knot = targetVertex(it);
		}
		else{
			if(rightBoundary(interval) <= act_prop.center) // interval to the right of current node?
			{
				if(atEnd(it)){
					TVertexDescriptor vd = addVertex(g);
					resizeVertexMap(g,pm);
					addEdge(g,act_knot,vd);
					_setIntervalTreeNode(property(pm,vd),(rightBoundary(interval)+leftBoundary(interval))/2,interval);
					break;
				}
				else
				{
					next_prop = property(pm,targetVertex(it));
					if(next_prop.center >= act_prop.center)
					{
						goNext(it);
						if(atEnd(it)){
							TVertexDescriptor vd = addVertex(g);
							resizeVertexMap(g,pm);
							addEdge(g,act_knot,vd);
							_setIntervalTreeNode(property(pm,vd),(rightBoundary(interval)+leftBoundary(interval))/2,interval);
							break;
						}
					}
				}
				act_knot = targetVertex(it);
			}
			else{ // need to create new node for interval
				_appendIntervalTreeNodeLists(property(pm, act_knot),interval);
				std::sort(begin(property(pm,act_knot).list1),end(property(pm,act_knot).list1),_less_compI1_ITree<typename Value<TList>::Type>);
				std::sort(begin(property(pm,act_knot).list2),end(property(pm,act_knot).list2),_greater_compI2_ITree<typename Value<TList>::Type>);
				break;
			}
		}
	}

}

template<typename TValue, typename TCargo, typename TInterval>
void
addInterval(IntervalTree<TValue,TCargo> & itree, TInterval interval)
{
SEQAN_CHECKPOINT

	++itree.interval_counter;
	addInterval(itree.g,itree.pm,interval);

}

template<typename TValue, typename TCargo>
void
addInterval(IntervalTree<TValue,TCargo> & itree, TValue begin, TValue end, TCargo cargo)
{
SEQAN_CHECKPOINT

	IntervalAndCargo<TValue,TCargo> interval;
	interval.i1 = begin;
	interval.i2 = end;
	interval.cargo = cargo;
	++itree.interval_counter;
	addInterval(itree.g,itree.pm,interval);

}

template<typename TValue, typename TCargo>
void
addInterval(IntervalTree<TValue,TCargo> & itree, TValue begin, TValue end)
{
SEQAN_CHECKPOINT

	IntervalAndCargo<TValue,TCargo> interval;
	interval.i1 = begin;
	interval.i2 = end;
	interval.cargo = itree.interval_counter;
	++itree.interval_counter;
	addInterval(itree.g,itree.pm,interval);

}

/**
.Function.findIntervals
..summary:Find all intervals that contain the query point or overlap with the query interval.
..cat:Miscellaneous
..class:Class.IntervalTree
..signature:findIntervals(intervalTree, query, result)
..signature:findIntervals(intervalTree, query_begin, query_end, result)
..signature:findIntervals(graph, propertyMap, query, result)
..param.intervalTree:An interval tree
...type:Class.IntervalTree
..param.query:A query point.
..param.query_begin:The begin position of the query interval.
..param.query_end:The end position of the query interval.
..param.result:A reference to the result string of $TCargo$ objects
...type:Class.String
..include:seqan/misc/misc_interval_tree.h
*/
template <typename TSpec, typename TPropertyMap, typename TValue, typename TCargo>
inline void
findIntervals(
    Graph<TSpec> const & g,
    TPropertyMap const & pm,
    TValue query,
    String<TCargo> & result)
{
SEQAN_CHECKPOINT

    typedef Graph<TSpec> const TGraph;
	typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
	typedef typename Value<TPropertyMap>::Type TProperty;
	typedef typename Value<TProperty>::Type TPropertyValue;
	typedef typename Iterator<TGraph, OutEdgeIterator>::Type TOutEdgeIterator;
    
    resize(result,0);
    if (empty(g)) return;

	// start at root
	TVertexDescriptor act_knot = 0;
	TProperty act_prop = property(pm,act_knot);
	TProperty next_prop;
		
	while(true)
	{
        typename Iterator<Graph<TSpec> , OutEdgeIterator>::Type it7;
        Iter<Graph<TSpec>, GraphIterator<InternalOutEdgeIterator<OutEdgeIterator> > > it5(g,act_knot);
        TOutEdgeIterator it4;
		TOutEdgeIterator it(g, act_knot);
		act_prop = property(pm,act_knot);
		if(act_prop.center < (TPropertyValue)query) // look in current node and right subtree
		{
			unsigned int i = 0;
			while(i < length(act_prop.list2) && rightBoundary(value(act_prop.list2,i)) > (TPropertyValue)query)
			{
				appendValue(result,cargo(value(act_prop.list2,i)), Generous());
				++i;	
			}
			if(atEnd(it)) break;
			else{
				next_prop = property(pm,targetVertex(it));
				if(next_prop.center <= act_prop.center)
				{
					goNext(it);
					if(atEnd(it)) break;
				}
			}
			act_knot = targetVertex(it);
		}
		else{
			if((TPropertyValue)query < act_prop.center) // look in current node and left subtree
			{
				unsigned int i = 0;
				while(i < length(act_prop.list1) && leftBoundary(value(act_prop.list1,i)) <= (TPropertyValue)query)
				{
					appendValue(result,cargo(value(act_prop.list1,i)), Generous());
					++i;
				}
				if(atEnd(it)) break;
				else
				{
					next_prop = property(pm,targetVertex(it));
					if(next_prop.center >= act_prop.center)
					{
						goNext(it);
						if(atEnd(it)) break;
					}
				}
				act_knot = targetVertex(it);
			}
			else{ // look in current node only, as query is center
				for(unsigned int i = 0; i < length(act_prop.list1); ++i)
                    appendValue(result, cargo(value(act_prop.list1,i)), Generous());
				break;
			}
		}
	}

}

template <typename TValue, typename TCargo>
inline void
findIntervals(
    IntervalTree<TValue,TCargo> const & it,
    TValue query,
    String<TCargo> & result)
{
SEQAN_CHECKPOINT
	findIntervals(it.g,it.pm,query,result);
}

template <typename TValue, typename TCargo>
inline void
findIntervals(
    IntervalTree<TValue,TCargo> const & tree,
    TValue query_begin,
    TValue query_end,
    String<TCargo> & result)
{
SEQAN_CHECKPOINT
	findIntervals(tree.g,tree.pm,query_begin,query_end,result);
}

template <typename TSpec, typename TPropertyMap, typename TValue, typename TCargo>
inline void
findIntervals(
    Graph<TSpec> const & g,
    TPropertyMap const & pm,
    TValue query_begin,
    TValue query_end,
    String<TCargo> & result)
{
    typedef Graph<TSpec> const TGraph;
	typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;

    resize(result,0);

	// start at root
	TVertexDescriptor act_knot = 0;
	findIntervals(g, pm, act_knot, query_begin, query_end, result);
}


template <
    typename TSpec,
    typename TPropertyMap,
    typename TVertexDescriptor,
    typename TValue,
    typename TCargo >
inline void
findIntervals(
    Graph<TSpec> const & g,
    TPropertyMap const & pm,
    TVertexDescriptor & act_knot,
    TValue query_begin, 
    TValue query_end, 
    String<TCargo> & result)
{
SEQAN_CHECKPOINT

	typedef typename Value<TPropertyMap>::Type TProperty;
	typedef typename Iterator<Graph<TSpec>, OutEdgeIterator>::Type TOutEdgeIterator;

    if (empty(g)) return;

    TProperty act_prop = property(pm,act_knot);
	TProperty next_prop;
		
	while(true)
	{
		TOutEdgeIterator it(g, act_knot);
		act_prop = property(pm,act_knot);
		//
		if(act_prop.center < query_begin) // query interval is to the right of node center
		{
			unsigned int i = 0;
			while(i < length(act_prop.list2) && rightBoundary(value(act_prop.list2,i)) > query_begin)
			{
				appendValue(result,cargo(value(act_prop.list2,i)), Generous());
				++i;	
			}
			if(atEnd(it)) break;
			else{
				next_prop = property(pm,targetVertex(it));
				if(next_prop.center <= act_prop.center)
				{
					goNext(it);
					if(atEnd(it)) break;
				}
			}
			act_knot = targetVertex(it);
		}
		else{
			if(query_end < act_prop.center) // query interval is to the left of node center
			{
				unsigned int i = 0;
				while(i < length(act_prop.list1) && leftBoundary(value(act_prop.list1,i)) < query_end)
				{
					appendValue(result,cargo(value(act_prop.list1,i)), Generous());
					++i;
				}
				if(atEnd(it)) break;
				else
				{
					next_prop = property(pm,targetVertex(it));
					if(next_prop.center >= act_prop.center)
					{
						goNext(it);
						if(atEnd(it)) break;
					}
				}
				act_knot = targetVertex(it);
			}
			else{//node center is contained in query interval
				for(unsigned int i = 0; i < length(act_prop.list1); ++i)
                    appendValue(result, cargo(value(act_prop.list1,i)), Generous());
				
				while(!atEnd(it))
				{
					TVertexDescriptor next_knot = targetVertex(it);
					findIntervals(g,pm, next_knot, query_begin, query_end, result);
					goNext(it);
				}
				break;

				//break; //dont break! continue in both subtrees!!
			}
		}
	}
}


/**
.Function.findIntervalsExcludeTouching
..summary:Find all intervals that contain the query point, exclude intervals that touch the query, i.e. where the query point equals the start or end point.
..signature:findIntervalsExcludeTouching(intervalTree, query, result)
..signature:findIntervalsExcludeTouching(graph, propertyMap, query, result)
..cat:Miscellaneous
..class:Class.IntervalTree
..param.intervalTree:An interval tree
...type:Class.IntervalTree
..param.graph:The directed graph that contains the topography of the interval tree.
..param.propertyMap:The property map containing the node properties of the interval tree
..param.query:The TValue to query here.
..param.result:The resulting string of cargos/ids of the intervals that contain the query point.
...type:Class.String
...remarks:Should be a string of TCargo.
..include:seqan/misc/misc_interval_tree.h
*/
template <typename TSpec, typename TPropertyMap, typename TValue, typename TCargo>
inline void
findIntervalsExcludeTouching(
    Graph<TSpec> const & g,
    TPropertyMap const & pm,
    TValue query,
    String<TCargo> & result)
{
SEQAN_CHECKPOINT

    typedef Graph<TSpec> const TGraph;
	typedef typename Iterator<TGraph, OutEdgeIterator >::Type TOutEdgeIterator;
	typedef typename VertexDescriptor<TGraph>::Type TVertexDescriptor;
	typedef typename Value<TPropertyMap>::Type TProperty;
	
	resize(result,0);
    if (empty(g)) return;

	// start at root
	TVertexDescriptor act_knot = 0;
	TProperty act_prop = property(pm,act_knot);
	TProperty next_prop;
		
	while(true)
	{
		TOutEdgeIterator it(g, act_knot);
		act_prop = property(pm,act_knot);
		if( (TValue) act_prop.center < query) // look in current node and right subtree
		{
			int i = 0;
			while(i < (int) length(act_prop.list2) && (TValue) rightBoundary(value(act_prop.list2,i)) > query)
			{
				appendValue(result,cargo(value(act_prop.list2,i)), Generous());
				++i;	
			}
			if(atEnd(it)) break;
			else{
				next_prop = property(pm,targetVertex(it));
				if(next_prop.center <= act_prop.center)
				{
					goNext(it);
					if(atEnd(it)) break;
				}
			}
			act_knot = targetVertex(it);
		}
		else{
			if(query < (TValue) act_prop.center) // look in current node and left subtree
			{
				int i = 0;
				while(i < (int) length(act_prop.list1) && (TValue) leftBoundary(value(act_prop.list1,i)) < query)
				{
					appendValue(result,cargo(value(act_prop.list1,i)), Generous());
					++i;
				}
				if(atEnd(it)) break;
				else
				{
					next_prop = property(pm,targetVertex(it));
					if(next_prop.center >= act_prop.center)
					{
						goNext(it);
						if(atEnd(it)) break;
					}
				}
				act_knot = targetVertex(it);
			}
			else{ // look in current node only
				int i = 0;
				while(i < (int) length(act_prop.list1) && (TValue) leftBoundary(value(act_prop.list1,i)) < query)
				{
					appendValue(result,cargo(value(act_prop.list1,i)), Generous());
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
    IntervalTree<TValue,TCargo> const & tree,
    TValue query,
    String<TCargo> & result)
{
SEQAN_CHECKPOINT
	findIntervalsExcludeTouching(tree.g,tree.pm,query,result);
}


/**
.Function.removeInterval
..summary:Removes an interval from the interval tree.
..signature:removeInterval(intervalTree, i_begin, i_end, i_id)
..cat:Miscellaneous
..class:Class.IntervalTree
..param.intervalTree:An interval tree
...type:Class.IntervalTree
..param.i_begin:The begin position of the interval to be removed.
..param.i_end:The end position of the interval to be removed.
..param.i_id:The ID of the interval to be removed.
*/

template <
    typename TSpec,
    typename TPropertyMap,
    typename TVertexDescriptor,
    typename TValue,
    typename TCargo >
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

    if (empty(g)) return false;

    TProperty act_prop = property(pm,act_knot);
	TProperty next_prop;
		
	while(true)
	{
		TOutEdgeIterator it(g, act_knot);
		act_prop = property(pm,act_knot);
		//
		if(act_prop.center < i_begin) // interval is to the right of node center
		{
			if(atEnd(it)) break;
			else{
				next_prop = property(pm,targetVertex(it));
				if(next_prop.center <= act_prop.center)
				{
					goNext(it);
					if(atEnd(it)) break;
				}
			}
			act_knot = targetVertex(it);
		}
		else{
			if(i_end < act_prop.center) // interval is to the left of node center
			{
				if(atEnd(it)) break;
				else
				{
					next_prop = property(pm,targetVertex(it));
					if(next_prop.center >= act_prop.center)
					{
						goNext(it);
						if(atEnd(it)) break;
					}
				}
				act_knot = targetVertex(it);
			}
			else{//node center is contained in interval, this is where we should find the interval to be removed
                // remove from list1
                TProperty& change_prop = property(pm,act_knot);
                bool foundInLeft = false;
                TListIterator list_it_keep = begin(change_prop.list1,Standard());
                TListIterator list_it = begin(change_prop.list1,Standard());
				while(list_it != end(change_prop.list1,Standard()))
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
                list_it_keep = begin(change_prop.list2,Standard());
                list_it = begin(change_prop.list2,Standard());
				while(list_it != end(change_prop.list2,Standard()))
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
				if(foundInRight && foundInLeft)
                {
                    resize(change_prop.list2, length(change_prop.list2)-1);
                    resize(change_prop.list1, length(change_prop.list1)-1);
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
    typename TCargo >
inline bool
removeInterval(
    Graph<TSpec> & g,
    TPropertyMap & pm,
    TValue i_begin,
    TValue i_end,
    TCargo i_id)
{
SEQAN_CHECKPOINT

	typedef typename VertexDescriptor<Graph<TSpec> >::Type TVertexDescriptor;

	// start looking at root
	TVertexDescriptor act_knot = 0;
	return removeInterval(g, pm, act_knot, i_begin, i_end, i_id);
}

template <typename TValue, typename TCargo>
inline bool
removeInterval(
    IntervalTree<TValue,TCargo> & tree,
    TValue i_begin,
    TValue i_end,
    TCargo i_id)
{
SEQAN_CHECKPOINT

	return removeInterval(tree.g, tree.pm, i_begin, i_end, i_id);
    // we do not decrease the interval_counter of tree, as it would mix up interval IDs
}




/////////////////// Metafunctions ///////////////////////

///.Metafunction.Value.param.T.type:Class.IntervalTree
template<typename TValue, typename TCargo>
struct Value<IntervalTree<TValue,TCargo> >
{
	typedef TValue Type;
};


///.Metafunction.Cargo.param.T.type:Class.IntervalTree
template<typename TValue, typename TCargo>
struct Cargo<IntervalTree<TValue,TCargo> >
{
	typedef TCargo Type;
};

}  // namespace SEQAN_NAMESPACE_MAIN

#endif  //#ifndef SEQAN_MISC_INTERVAL_TREE_H
