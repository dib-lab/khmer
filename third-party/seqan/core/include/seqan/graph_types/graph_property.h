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

#ifndef SEQAN_HEADER_GRAPH_PROPERTY_H
#define SEQAN_HEADER_GRAPH_PROPERTY_H

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////
//	Graph - External property map
//////////////////////////////////////////////////////////////////////////////

/**
.Class.External Property Map
..cat:Graph
..summary:An external property map.
..remarks:The external property map is assumed to be an instance of @Class.String@.
It is indexed via VertexDescriptors or EdgeDescriptors.
..signature:String<TValue, TSpec>
..param.TValue:The value type. That is the type of information stored in the property map.
...metafunction:Metafunction.Value
..param.TSpec:The specializing type.
...metafunction:Metafunction.Spec
...default:$Alloc<>$, see @Spec.Alloc String@.
..include:seqan/graph_types.h
*/


//////////////////////////////////////////////////////////////////////////////
// FUNCTIONS
//////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////

/**
.Function.resizeVertexMap
..class:Class.Graph
..cat:Graph
..summary:Initializes a vertex map. 
..signature:resizeVertexMap(g, pm [, prototype])
..param.g:A Graph.
...type:Class.Graph
..param.pm:An External Property Map.
...type:Class.External Property Map
..param.prototype:An optional prototype that is used for initializing the property map.
..returns:void
..see:Function.resizeEdgeMap
..include:seqan/graph_types.h
*/

template<typename TSpec, typename TPropertyMap>
inline void
resizeVertexMap(Graph<TSpec> const & g,
                TPropertyMap & pm)
{
	SEQAN_CHECKPOINT;
    typedef typename Value<TPropertyMap>::Type TValue;
	resize(pm, getIdUpperBound(_getVertexIdManager(g)), TValue(), Generous());
}

template<typename TSpec, typename TPropertyMap, typename TPrototype>
inline void
resizeVertexMap(Graph<TSpec> const & g,
                TPropertyMap & pm,
                TPrototype const & prototype)
{
	SEQAN_CHECKPOINT;
	resize(pm, getIdUpperBound(_getVertexIdManager(g)), prototype, Generous());
}

//////////////////////////////////////////////////////////////////////////////

/**
.Function.resizeEdgeMap:
..cat:Graph
..summary:Initializes an edge map
..signature:resizeEdgeMap(g, pm [, prototype])
..param.g:A Graph.
...type:Class.Graph
..param.pm:An External or Internal Property Map.
...type:Class.External Property Map
...type:Class.InternalMap
...type:Class.InternalPointerMap
...type:Class.InternalRawMap
..param.prototype:An optional prototype that is used for initializing the property map.
..returns:void
..see:Function.resizeVertexMap
..include:seqan/graph_types.h
*/

template<typename TSpec, typename TPropertyMap>
inline void
resizeEdgeMap(Graph<TSpec> const & g,
			  TPropertyMap & pm)
{
	SEQAN_CHECKPOINT;
    typedef typename Value<TPropertyMap>::Type TValue;
	resize(pm, getIdUpperBound(_getEdgeIdManager(g)), TValue(), Generous());
}

template<typename TSpec, typename TPropertyMap, typename TPrototype>
inline void
resizeEdgeMap(Graph<TSpec> const & g,
              TPropertyMap & pm,
              TPrototype const & prototype)
{
	SEQAN_CHECKPOINT;
	resize(pm, getIdUpperBound(_getEdgeIdManager(g)), prototype, Generous());
}

//////////////////////////////////////////////////////////////////////////////

/**
.Function.assignProperty:
..class:Class.External Property Map
..class:Class.InternalMap
..class:Class.InternalPointerMap
..class:Class.InternalRawMap
..cat:Graph
..summary:Assigns a property to an item in the property map.
..signature:assignProperty(pm, d, val)
..param.pm:An External or Internal Property Map.
...type:Class.External Property Map
...type:Class.InternalMap
...type:Class.InternalPointerMap
...type:Class.InternalRawMap
..param.d:A vertex or edge descriptor.
...remarks:Identifies the item in the property map.
...type:Metafunction.VertexDescriptor
...type:Metafunction.EdgeDescriptor
..param.val:The new value.
...remarks:Type of the new value must match the value type of the property map.
See @Metafunction.Value@.
..returns:void
..see:Function.getProperty
..see:Function.property
..include:seqan/graph_types.h
*/

template<typename TPropertyMap, typename TDescriptor, typename TValue>
inline void
assignProperty(TPropertyMap& pm,
			   TDescriptor const d,
			   TValue const val)
{
	SEQAN_CHECKPOINT
	assignValue(pm, _getId(d), val);
}

//////////////////////////////////////////////////////////////////////////////

/**
.Function.property
..class:Class.External Property Map
..class:Class.InternalMap
..class:Class.InternalPointerMap
..class:Class.InternalRawMap
..cat:Graph
..summary:Accesses the property of an item in the property map.
..signature:property(pm, d)
..param.pm:An External or Internal Property Map.
...type:Class.External Property Map
...type:Class.InternalMap
...type:Class.InternalPointerMap
...type:Class.InternalRawMap
..param.d:A vertex or edge descriptor.
...remarks:Identifies the item in the property map.
...type:Metafunction.VertexDescriptor
...type:Metafunction.EdgeDescriptor
..returns:Reference to the item in the property map.
..see:Function.getProperty
..see:Function.assignProperty
..include:seqan/graph_types.h
*/

template<typename TPropertyMap, typename TDescriptor>
inline typename Reference<TPropertyMap>::Type
property(TPropertyMap& pm,
		TDescriptor const d)
{
	SEQAN_CHECKPOINT
	return value(pm, _getId(d));
}

//////////////////////////////////////////////////////////////////////////////

template<typename TPropertyMap, typename TDescriptor>
inline typename Reference<TPropertyMap const>::Type
property(TPropertyMap const& pm,
		TDescriptor const d)
{
	SEQAN_CHECKPOINT
	return value(pm, _getId(d));
}

//////////////////////////////////////////////////////////////////////////////

/**
.Function.getProperty
..class:Class.External Property Map
..class:Class.InternalMap
..class:Class.InternalPointerMap
..class:Class.InternalRawMap
..cat:Graph
..summary:Get method for an item's property.
..signature:getProperty(pm, d)
..param.pm:An External or Internal Property Map.
...type:Class.External Property Map
...type:Class.InternalMap
...type:Class.InternalPointerMap
...type:Class.InternalRawMap
..param.d:A vertex or edge descriptor.
...remarks:Identifies the item in the property map.
...type:Metafunction.VertexDescriptor
...type:Metafunction.EdgeDescriptor
..returns:Value of the item in the property map.
..see:Function.property
..see:Function.assignProperty
..include:seqan/graph_types.h
*/

template<typename TPropertyMap, typename TDescriptor>
inline typename GetValue<TPropertyMap const>::Type
getProperty(TPropertyMap const& pm,
			TDescriptor const d)
{
	SEQAN_CHECKPOINT
	return getValue(pm, _getId(d));
}






//////////////////////////////////////////////////////////////////////////////
// Graph - Internal Property Manager using member ids (only for edges!!!)
//////////////////////////////////////////////////////////////////////////////

/**
.Class.InternalMap:
..cat:Graph
..summary:An internal property map using member ids.
..remarks:Internal property maps are used to access internal edge cargos.
..signature:InternalMap<TContainer, MemberId>
..param.TContainer:The cargo type.
...metafunction:Metafunction.Cargo
..param.MemberId:An unsigned int.
...remarks:Specifies the position of the member in the cargo.
Note: If zero it is assumed that the cargo is a simple type (e.g., int).
...default:$0$.
..include:seqan/graph_types.h
..see:Class.External Property Map
..see:Class.InternalPointerMap
..see:Class.InternalRawMap
*/
template<typename TContainer, unsigned int const MemberId = 0>
struct InternalMap 
{
};


//////////////////////////////////////////////////////////////////////////////
//	Internal Property Manager using member Ids - Metafunctions
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

///.Metafunction.Value.param.T.type:Class.InternalMap
///.Metafunction.Value.class:Class.InternalMap

template<typename T1, typename T2>
struct Value<InternalMap<Pair<T1, T2>, 1> const> {
	typedef T1 const Type;
};

//////////////////////////////////////////////////////////////////////////////

template<typename T1, typename T2>
struct Value<InternalMap<Pair<T1, T2>, 1> > {
	typedef T1 Type;
};

//////////////////////////////////////////////////////////////////////////////

template<typename T1, typename T2>
struct Value<InternalMap<Pair<T1, T2>, 2> const> {
	typedef T2 const Type;
};

//////////////////////////////////////////////////////////////////////////////

template<typename T1, typename T2>
struct Value<InternalMap<Pair<T1, T2>, 2> > {
	typedef T2 Type;
};

//////////////////////////////////////////////////////////////////////////////

template<typename T>
struct Value<InternalMap<T, 0> const> {
	typedef T const Type;
};

//////////////////////////////////////////////////////////////////////////////

template<typename T>
struct Value<InternalMap<T, 0> > {
	typedef T Type;
};


//////////////////////////////////////////////////////////////////////////////
// Internal Property Manager using member Ids - FUNCTIONS
//////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////

template<typename TSpec, typename TContainer, unsigned int const MemberId>
inline void
resizeEdgeMap(Graph<TSpec> const&,
			  InternalMap<TContainer, MemberId>&)
{
}

//////////////////////////////////////////////////////////////////////////////

template<typename TSpec, typename TContainer, unsigned int const MemberId>
inline void
resizeEdgeMap(Graph<TSpec>&,
			  InternalMap<TContainer, MemberId>&)
{
}

//////////////////////////////////////////////////////////////////////////////

template<typename TSpec, typename TContainer, unsigned int const MemberId, typename TSource>
inline void
assignEdgeMap(Graph<TSpec> const &,
			  InternalMap<TContainer, MemberId> &,
              TSource const &)
{
}

//////////////////////////////////////////////////////////////////////////////

template<typename TSpec, typename TContainer, unsigned int const MemberId, typename TSource>
inline void
assignEdgeMap(Graph<TSpec> &,
			  InternalMap<TContainer, MemberId> &,
              TSource const &)
{
}

//////////////////////////////////////////////////////////////////////////////

template<typename T1, typename T2, typename TEdgeDescriptor, typename TValue>
inline void
assignProperty(InternalMap<Pair<T1, T2>, 1>&,
			   TEdgeDescriptor const e,
			   TValue const val)
{
	SEQAN_CHECKPOINT
	(cargo(e)).i1 = val;
}

//////////////////////////////////////////////////////////////////////////////

template<typename T1, typename T2, typename TEdgeDescriptor, typename TValue>
inline void
assignProperty(InternalMap<Pair<T1, T2>, 2>&,
			   TEdgeDescriptor const e,
			   TValue const val)
{
	SEQAN_CHECKPOINT
	(cargo(e)).i2 = val;
}

//////////////////////////////////////////////////////////////////////////////

template<typename T, typename TEdgeDescriptor, typename TValue>
inline void
assignProperty(InternalMap<T, 0>&,
			   TEdgeDescriptor const e,
			   TValue const val)
{
	SEQAN_CHECKPOINT
	assignCargo(e, val);
}

//////////////////////////////////////////////////////////////////////////////

template<typename T1, typename T2, typename TEdgeDescriptor>
inline typename Value<InternalMap<Pair<T1, T2>, 2> >::Type&
property(InternalMap<Pair<T1, T2>, 2>&,
		 TEdgeDescriptor e)
{
	SEQAN_CHECKPOINT
	return (cargo(e)).i2;
}

//////////////////////////////////////////////////////////////////////////////

template<typename T1, typename T2, typename TEdgeDescriptor>
inline typename Value<InternalMap<Pair<T1, T2>, 2> const>::Type&
property(InternalMap<Pair<T1, T2>, 2> const&,
		 TEdgeDescriptor e)
{
	SEQAN_CHECKPOINT
	return (cargo(e)).i2;
}

//////////////////////////////////////////////////////////////////////////////

template<typename T1, typename T2, typename TEdgeDescriptor>
inline typename Value<InternalMap<Pair<T1, T2>, 1> >::Type&
property(InternalMap<Pair<T1, T2>, 1>&,
		 TEdgeDescriptor e)
{
	SEQAN_CHECKPOINT
	return (cargo(e)).i1;
}

//////////////////////////////////////////////////////////////////////////////

template<typename T1, typename T2, typename TEdgeDescriptor>
inline typename Value<InternalMap<Pair<T1, T2>, 1> const>::Type&
property(InternalMap<Pair<T1, T2>, 1> const&,
		 TEdgeDescriptor e)
{
	SEQAN_CHECKPOINT
	return (cargo(e)).i1;
}

//////////////////////////////////////////////////////////////////////////////

template<typename T, typename TEdgeDescriptor>
inline typename Value<InternalMap<T, 0> >::Type&
property(InternalMap<T, 0>&,
		 TEdgeDescriptor e)
{
	SEQAN_CHECKPOINT
	return cargo(e);
}

//////////////////////////////////////////////////////////////////////////////

template<typename T, typename TEdgeDescriptor>
inline typename Value<InternalMap<T, 0> const>::Type&
property(InternalMap<T, 0> const&,
		 TEdgeDescriptor e)
{
	SEQAN_CHECKPOINT
	return cargo(e);
}

//////////////////////////////////////////////////////////////////////////////

template<typename T1, typename T2, typename TEdgeDescriptor>
inline typename Value<InternalMap<Pair<T1, T2>, 1> >::Type
getProperty(InternalMap<Pair<T1, T2>, 1> const&,
			TEdgeDescriptor e)
{
	SEQAN_CHECKPOINT
	return (getCargo(e)).i1;
}

//////////////////////////////////////////////////////////////////////////////

template<typename T1, typename T2, typename TEdgeDescriptor>
inline typename Value<InternalMap<Pair<T1, T2>, 1> >::Type
getProperty(InternalMap<Pair<T1, T2>, 1>&,
			TEdgeDescriptor e)
{
	SEQAN_CHECKPOINT
	return (getCargo(e)).i1;
}

//////////////////////////////////////////////////////////////////////////////

template<typename T1, typename T2, typename TEdgeDescriptor>
inline typename Value<InternalMap<Pair<T1, T2>, 2> >::Type
getProperty(InternalMap<Pair<T1, T2>, 2> const&,
			TEdgeDescriptor e)
{
	SEQAN_CHECKPOINT
	return (getCargo(e)).i2;
}

//////////////////////////////////////////////////////////////////////////////

template<typename T1, typename T2, typename TEdgeDescriptor>
inline typename Value<InternalMap<Pair<T1, T2>, 2> >::Type
getProperty(InternalMap<Pair<T1, T2>, 2>&,
			TEdgeDescriptor e)
{
	SEQAN_CHECKPOINT
	return (getCargo(e)).i2;
}

//////////////////////////////////////////////////////////////////////////////

template<typename T, typename TEdgeDescriptor>
inline typename Value<InternalMap<T, 0> >::Type
getProperty(InternalMap<T, 0> const&,
			TEdgeDescriptor e)
{
	SEQAN_CHECKPOINT
	return getCargo(e);
}

//////////////////////////////////////////////////////////////////////////////

template<typename T, typename TEdgeDescriptor>
inline typename Value<InternalMap<T, 0> >::Type
getProperty(InternalMap<T, 0>&,
			TEdgeDescriptor e)
{
	SEQAN_CHECKPOINT
	return getCargo(e);
}


//////////////////////////////////////////////////////////////////////////////
// Graph - Internal Property Manager using pointer to members (only for edges!!!)
//////////////////////////////////////////////////////////////////////////////

/**
.Class.InternalPointerMap:
..cat:Graph
..summary:An internal property map using pointer to members.
..remarks:Internal property maps are used to access internal edge cargos.
..signature:InternalPointerMap<TPropmap, Instance>
..param.TPropmap:A pointer to member type.
..param.Instance:A pointer to a member of type TPropmap.
..include:seqan/graph_types.h
..see:Class.External Property Map
..see:Class.InternalMap
..see:Class.InternalRawMap
*/
template <typename TPropmap, TPropmap const Instance> 
struct InternalPointerMap 
{
}; 

//////////////////////////////////////////////////////////////////////////////
//	Internal Property Manager using member Ids - Metafunctions
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

///.Metafunction.Value.param.T.type:Class.InternalPointerMap
///.Metafunction.Value.class:Class.InternalMap

template<typename TClass, typename TValue, TValue TClass:: * TPMember>
struct Value<InternalPointerMap<TValue TClass::*, TPMember> const> {
	typedef TValue const Type;
};

//////////////////////////////////////////////////////////////////////////////

template<typename TClass, typename TValue, TValue TClass:: * TPMember>
struct Value<InternalPointerMap<TValue TClass::*, TPMember> > {
	typedef TValue Type;
};


//////////////////////////////////////////////////////////////////////////////
// Internal Property Manager using pointer to members - FUNCTIONS
//////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////

template<typename TSpec, typename TPropmap, TPropmap const Instance>
inline void
resizeEdgeMap(Graph<TSpec>&,
			  InternalPointerMap<TPropmap, Instance>&)
{
	SEQAN_CHECKPOINT
}

//////////////////////////////////////////////////////////////////////////////

template<typename TSpec, typename TPropmap, TPropmap const Instance>
inline void
resizeEdgeMap(Graph<TSpec> const&,
			  InternalPointerMap<TPropmap, Instance>&)
{
}

//////////////////////////////////////////////////////////////////////////////

template<typename TClass, typename TValue, TValue TClass:: * TPMember, typename TEdgeDescriptor>
inline void
assignProperty(InternalPointerMap<TValue TClass::*, TPMember>&,
			TEdgeDescriptor const e,
			TValue const val)
{
	SEQAN_CHECKPOINT
	(cargo(e)).*TPMember = val;
}

//////////////////////////////////////////////////////////////////////////////

template<typename TClass, typename TValue, TValue TClass:: * TPMember, typename TEdgeDescriptor>
inline typename Value<InternalPointerMap<TValue TClass::*, TPMember> >::Type&
property(InternalPointerMap<TValue TClass::*, TPMember>&,
		 TEdgeDescriptor const e)
{
	SEQAN_CHECKPOINT
	return (cargo(e)).*TPMember;
}

//////////////////////////////////////////////////////////////////////////////

template<typename TClass, typename TValue, TValue TClass:: * TPMember, typename TEdgeDescriptor>
inline typename Value<InternalPointerMap<TValue TClass::*, TPMember> const>::Type&
property(InternalPointerMap<TValue TClass::*, TPMember> const&,
		 TEdgeDescriptor const e)
{
	SEQAN_CHECKPOINT
	return (cargo(e)).*TPMember;
}

//////////////////////////////////////////////////////////////////////////////

template<typename TClass, typename TValue, TValue TClass:: * TPMember, typename TEdgeDescriptor>
inline typename Value<InternalPointerMap<TValue TClass::*, TPMember> >::Type
getProperty(InternalPointerMap<TValue TClass::*, TPMember> const&,
			TEdgeDescriptor const e)
{
	SEQAN_CHECKPOINT
	return (getCargo(e)).*TPMember;
}

//////////////////////////////////////////////////////////////////////////////


template<typename TClass, typename TValue, TValue TClass:: * TPMember, typename TEdgeDescriptor>
inline typename Value<InternalPointerMap<TValue TClass::*, TPMember> >::Type
getProperty(InternalPointerMap<TValue TClass::*, TPMember>&,
			TEdgeDescriptor const e)
{
	SEQAN_CHECKPOINT
	return (getCargo(e)).*TPMember;
}



//////////////////////////////////////////////////////////////////////////////
// Graph - Internal Property Manager using raw pointer to member (only for edges!!!)
//////////////////////////////////////////////////////////////////////////////

/**
.Class.InternalRawMap:
..cat:Graph
..summary:An internal property map using raw pointer to members.
..remarks:Internal property maps are used to access internal edge cargos.
..include:seqan/graph_types.h
..see:Class.External Property Map
..see:Class.InternalMap
..see:Class.InternalPointerMap
*/

//////////////////////////////////////////////////////////////////////////////
// Internal Property Manager using raw pointer to member - Metafunctions
//////////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////

///.Metafunction.Value.param.T.type:Class.InternalRawMap
///.Metafunction.Value.class:Class.InternalRawMap

template <typename TClass, typename TValue> 
struct Value<TValue TClass:: *> {
	typedef TValue Type;
};

//////////////////////////////////////////////////////////////////////////////

template <typename TClass, typename TValue> 
struct Value<TValue TClass:: * const> {
	typedef TValue const Type;
};


//////////////////////////////////////////////////////////////////////////////
// Raw pointer to member - FUNCTIONS
//////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////

template <typename TSpec, typename TClass, typename TValue> 
inline void
resizeEdgeMap(Graph<TSpec> const&,
			  TValue TClass:: *)
{
}

//////////////////////////////////////////////////////////////////////////////

template <typename TSpec, typename TClass, typename TValue> 
inline void
resizeEdgeMap(Graph<TSpec>&,
			  TValue TClass:: *)
{
}

//////////////////////////////////////////////////////////////////////////////


template <typename TClass, typename TValue, typename TEdgeDescriptor> 
inline void 
assignProperty(TValue TClass:: * ptr_to_member, 
			TEdgeDescriptor const e, 
			TValue const val) 
{
	SEQAN_CHECKPOINT
	(cargo(e)).*ptr_to_member=val; 
}


//////////////////////////////////////////////////////////////////////////////

template <typename TClass, typename TValue, typename TEdgeDescriptor> 
inline TValue& 
property(TValue TClass:: * const ptr_to_member, 
		 TEdgeDescriptor const e) 
{
	SEQAN_CHECKPOINT
	return (cargo(e)).*ptr_to_member; 
}

//////////////////////////////////////////////////////////////////////////////

template <typename TClass, typename TValue, typename TEdgeDescriptor> 
inline TValue
getProperty(TValue TClass:: * const ptr_to_member, 
			TEdgeDescriptor const e) 
{
	SEQAN_CHECKPOINT
	return (getCargo(e)).*ptr_to_member; 
} 



//////////////////////////////////////////////////////////////////////////////
// Init functions for all Maps - FUNCTIONS
//////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////

/**
.Function.assignVertexMap
..class:Class.External Property Map
..class:Class.InternalMap
..class:Class.InternalPointerMap
..class:Class.InternalRawMap
..cat:Graph
..summary:Initializes a vertex map with values of an array.
..signature:assignVertexMap(g, pm, prop)
..param.g:A Graph.
...type:Class.Graph
..param.pm:An External Property Map.
...type:Class.External Property Map
..param.prop:An array with properties that are to be assigned to the items in the property map.
...remarks:For every vertex descriptor there must be an entry in the array.
..returns:void
..see:Function.assignEdgeMap
..include:seqan/graph_types.h
*/

template<typename TSpec, typename TPropertyMap, typename TProperties>
inline void
assignVertexMap(Graph<TSpec> const & g,
                TPropertyMap & pm,
                TProperties const & prop)
{
	SEQAN_CHECKPOINT
	resize(pm, getIdUpperBound(_getVertexIdManager(g)), Generous());
	typedef Graph<TSpec> TGraph;
	typedef typename Iterator<TGraph, VertexIterator>::Type TVertexIterator;
	TVertexIterator it(g);
	for (; !atEnd(it); goNext(it))
		assignProperty(pm,getValue(it), getValue(prop, _getId(value(it))));
}

//////////////////////////////////////////////////////////////////////////////

/**
.Function.assignEdgeMap
..class:Class.External Property Map
..class:Class.InternalMap
..class:Class.InternalPointerMap
..class:Class.InternalRawMap
..cat:Graph
..summary:Initializes a vertex map with values of an array.
..signature:assignEdgeMap(g, pm, prop)
..param.g:A Graph.
...type:Class.Graph
..param.pm:An External or Internal Property Map.
...type:Class.External Property Map
...type:Class.InternalMap
...type:Class.InternalPointerMap
...type:Class.InternalRawMap
..param.prop:An array with properties that are to be assigned to the items in the property map.
...remarks:For every edge id there must be an entry in the array.
..returns:void
..see:Function.assignVertexMap
..include:seqan/graph_types.h
*/

template<typename TSpec, typename TPropertyMap, typename TProperties>
inline void
assignEdgeMap(Graph<TSpec> const & g,
              TPropertyMap & pm,
              TProperties const & prop)
{
	SEQAN_CHECKPOINT
	resize(pm, getIdUpperBound(_getEdgeIdManager(g)), Generous());
	typedef Graph<TSpec> TGraph;
	typedef typename Iterator<TGraph, EdgeIterator>::Type TEdgeIterator;
	TEdgeIterator it(g);
	for (; !atEnd(it); goNext(it))
		assignProperty(pm,*it,prop[_getId(*it)]);
}

}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
