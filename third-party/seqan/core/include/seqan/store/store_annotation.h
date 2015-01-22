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
// Author: David Weese <david.weese@fu-berlin.de>
// ==========================================================================

#ifndef SEQAN_HEADER_STORE_ANNOTATION_H
#define SEQAN_HEADER_STORE_ANNOTATION_H

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////
// Annotation Store
//////////////////////////////////////////////////////////////////////////////

/**
.Class.AnnotationStoreElement
..summary:Represents an annotation of a contig feature.
..cat:Fragment Store
..signature:AnnotationStoreElement<>
..signature:AnnotationStoreElement<TPos[, TSpec]>
..param.TPos:Type to store (gap-space) positions.
..param.TSpec:The specialization type.
...default:$void$
..remarks:Value type of the @Memvar.FragmentStore#annotationStore@ string.
The id of an annotation (aka annotationId) is not stored explicitly, as it is implicitly given by the position in the @Memvar.FragmentStore#annotationStore@.
..include:seqan/store.h

.Typedef.AnnotationStoreElement#TId
..summary:Type of annotationId and @Memvar.AnnotationStoreElement#contigId@.
..remarks:$TId$ equals the result of $Id<AnnotationStoreElement<> >::Type$, see @Metafunction.Id@.
..class:Class.AnnotationStoreElement
..include:seqan/store.h

.Typedef.AnnotationStoreElement#TPos
..summary:Type of the @Memvar.AnnotationStoreElement#beginPos@ and @Memvar.AnnotationStoreElement#endPos@ members.
..class:Class.AnnotationStoreElement
..include:seqan/store.h

.Typedef.AnnotationStoreElement#TValues
..summary:@Class.StringSet@ type of the @Memvar.AnnotationStoreElement#values@ member.
..class:Class.AnnotationStoreElement
..include:seqan/store.h


.Memfunc.AnnotationStoreElement#AnnotationStoreElement
..summary:Constructor
..signature:AnnotationStoreElement()
..remarks:The default constructor sets all members to @Memvar.AnnotationStoreElement#INVALID_ID@ and 
@Memvar.AnnotationStoreElement#beginPos@ and @Memvar.AnnotationStoreElement#endPos@ to @Memvar.AnnotationStoreElement#INVALID_POS@.
..class:Class.AnnotationStoreElement
..include:seqan/store.h

.Memvar.AnnotationStoreElement#contigId
..summary:Refers to the contig in the @Memvar.FragmentStore#contigStore@ the annotation is part of.
..type:Metafunction.Id
..class:Class.AnnotationStoreElement
..include:seqan/store.h

.Memvar.AnnotationStoreElement#typeId
..summary:Refers to an entry in the @Memvar.FragmentStore#annotationTypeStore@. 
There are some type ids predefined for commonly used types, e.g. $ANNO_GENE$. See @Enum.Predefined Annotation Types@.
..type:Metafunction.Id
..class:Class.AnnotationStoreElement
..include:seqan/store.h

.Memvar.AnnotationStoreElement#beginPos
..summary:Begin position of the annotation in gap-space.
..type:Typedef.AnnotationStoreElement#TPos
..class:Class.AnnotationStoreElement
..include:seqan/store.h

.Memvar.AnnotationStoreElement#endPos
..summary:End position of the annotation in gap-space. If @Memvar.AnnotationStoreElement#endPos@ < @Memvar.AnnotationStoreElement#beginPos@, 
the annotated feature is located on the reverse strand, where @Memvar.AnnotationStoreElement#beginPos@ and @Memvar.AnnotationStoreElement#endPos@
are the corresponding positions on the forward strand.
..type:Typedef.AnnotationStoreElement#TPos
..class:Class.AnnotationStoreElement
..include:seqan/store.h

.Memvar.AnnotationStoreElement#values
..summary:@Class.StringSet@ that stores additional annotation values addressed by $keyId$. The GFF/GTF file format allows to define user-specific key-value pairs. The set of all keys addressed by $keyId$ are stored in the @Memvar.FragmentStore#annotationKeyStore@.
..type:Typedef.AnnotationStoreElement#TValues
..class:Class.AnnotationStoreElement
..include:seqan/store.h

.Memvar.AnnotationStoreElement#parentId
..summary:The id of the parent annotation.
..type:Metafunction.Id
..class:Class.AnnotationStoreElement
..include:seqan/store.h

.Memvar.AnnotationStoreElement#nextSiblingId
..summary:The id of the right sibling annotation.
..type:Metafunction.Id
..class:Class.AnnotationStoreElement
..include:seqan/store.h

.Memvar.AnnotationStoreElement#lastChildId
..summary:The id of the rightmost child annotation.
..type:Metafunction.Id
..class:Class.AnnotationStoreElement
..include:seqan/store.h

.Memvar.AnnotationStoreElement#INVALID_ID
..summary:Constant to represent an invalid id.
..type:Metafunction.Id
..class:Class.AnnotationStoreElement
..include:seqan/store.h

.Memvar.AnnotationStoreElement#INVALID_POS
..summary:Constant to represent an invalid position.
..type:Typedef.AnnotationStoreElement#TPos
..class:Class.AnnotationStoreElement
..include:seqan/store.h
*/

template <typename TPos_, typename TSpec = void>
struct AnnotationStoreElement
{
	typedef typename Id<AnnotationStoreElement>::Type		TId;
	typedef TPos_											TPos;
	typedef StringSet<CharString, Owner< ConcatDirect<> > >	TValues;

	static const TId  INVALID_ID;
	static const TPos INVALID_POS;
	
	TId					parentId;
	TId					contigId;
	TId					countId;
	TId					typeId;			// gene, intron, ...

	TPos				beginPos;		// begin position of the annotation in the gapped contig sequence (i.e. in gap-space)
	TPos				endPos;			// end position of ..., for annotations on the reverse strand holds end < begin
	
	TId					lastChildId;	// generated back links to child
	TId					nextSiblingId;	// and sibling
	
	TValues				values;			// stores values for each keyId of (key,value) pairs

	AnnotationStoreElement() : 
		parentId(INVALID_ID), contigId(INVALID_ID), countId(INVALID_ID), typeId(INVALID_ID), 
		beginPos(INVALID_POS), endPos(INVALID_POS),
		lastChildId(INVALID_ID), nextSiblingId(INVALID_ID) {}
};

//////////////////////////////////////////////////////////////////////////////

template <typename TPos, typename TSpec>
const typename Id<AnnotationStoreElement<TPos, TSpec> >::Type
AnnotationStoreElement<TPos, TSpec>::INVALID_ID = MaxValue<typename Id<AnnotationStoreElement<TPos, TSpec> >::Type>::VALUE;

template <typename TPos, typename TSpec>
const TPos
AnnotationStoreElement<TPos, TSpec>::INVALID_POS = MaxValue<TPos>::VALUE;

//////////////////////////////////////////////////////////////////////////////

template <typename TSpec = void>
struct AnnotationTree {};

/**
.Spec.AnnotationTree Iterator:
..cat:FragmentStore
..summary:Iterator of the annotation tree represented by a @Class.FragmentStore@.
..remarks:This iterator can move @Function.AnnotationTree Iterator#goDown|down@, @Function.AnnotationTree Iterator#goRight|right@, and @Function.AnnotationTree Iterator#goUp|up@ in the tree and supports a preorder dfs traversal via the functions @Function.goBegin@, @Function.goNext@, and @Function.atEnd@.
Preorder means that the iterator visits a node before its children.
..remarks:To access the annotation, the iterator points to, use @Function.getAnnotation@. The annotation id is returned by @Function.value@.
..signature:Iter<TFragmentStore, AnnotationTree<> >
..signature:Iterator<TFragmentStore, AnnotationTree<> >::Type
..general:Class.Iter
..implements:Concept.RootedIteratorConcept
..param.TFragmentStore:A FragmentStore class.
...type:Class.FragmentStore
..include:seqan/store.h
..example:
...image:AnnotationTree|Typical annotation tree hierarchy.
..example:
...text:A new annotation tree iterator can be instantiated as follows:
...code:
Iterator<FragmentStore<>, AnnotationTree<> >::Type it;
it = begin(store, AnnotationTree<>());
...text:Or shorter (see @Memfunc.AnnotationTree Iterator#AnnotationTree Iterator|AnnotationTree Iterator constructor@):
...code:
Iterator<FragmentStore<>, AnnotationTree<> >::Type it(store);

.Memfunc.AnnotationTree Iterator#AnnotationTree Iterator
..summary:Constructor
..signature:Iter()
..signature:Iter(store [, startInNode])
..class:Spec.AnnotationTree Iterator
..param.store:A @Class.FragmentStore@ object.
...type:Class.FragmentStore
..param.startInNode:Annotation id of the node the iterator should start at.
...default:$0$, the id of the root node.
..include:seqan/store.h
..remarks:The @Function.begin@ function can also be used to create a tree iterator that starts in the root node:
...code:
Iterator<FragmentStore<>, AnnotationTree<> >::Type it;
it = begin(store, AnnotationTree<>());
*/

template <typename TFragmentStore, typename TSpec>
class Iter<TFragmentStore, AnnotationTree<TSpec> >
{
public:
	typedef typename TFragmentStore::TAnnotationStore	TAnnotationStore;
	typedef typename Value<TAnnotationStore>::Type		TAnnotation;
	typedef typename TAnnotation::TId					TId;

	TFragmentStore *store;
	TId _id;
	
	Iter():
		store(NULL),
		_id(TAnnotation::INVALID_ID) {}

	Iter(TFragmentStore &_store):
		store(&_store),
		_id(0) {}

	Iter(TFragmentStore &_store, TId startInNode):
		store(&_store),
		_id(startInNode) {}

	Iter(TFragmentStore &_store, MinimalCtor):
		store(&_store),
		_id(TAnnotation::INVALID_ID) {}

	inline Iter const &
	operator = (Iter const &_origin)
	{
		store = &container(_origin);
		_id = _origin._id;
		return *this;
	}
};

//////////////////////////////////////////////////////////////////////////////

template <typename TFragmentStore, typename TSpec>
struct Iterator< TFragmentStore, AnnotationTree<TSpec> > {
	typedef Iter< TFragmentStore, AnnotationTree<TSpec> > Type;
};

//////////////////////////////////////////////////////////////////////////////

template <typename TFragmentStore, typename TSpec>
struct Value< Iter< TFragmentStore, AnnotationTree<TSpec> > >:
	VertexDescriptor<TFragmentStore> {};

template <typename TFragmentStore, typename TSpec>
struct Size< Iter< TFragmentStore, AnnotationTree<TSpec> > > :
	Size<TFragmentStore> {};

template <typename TFragmentStore, typename TSpec>
struct Position< Iter< TFragmentStore, AnnotationTree<TSpec> > > :
	Position<TFragmentStore> {};

//////////////////////////////////////////////////////////////////////////////

template <typename TFragmentStore, typename TSpec>
inline typename VertexDescriptor<TFragmentStore>::Type & 
value(Iter< TFragmentStore, AnnotationTree<TSpec> > &it) { 
	return it._id;
}

template <typename TFragmentStore, typename TSpec>
inline typename VertexDescriptor<TFragmentStore>::Type const & 
value(Iter< TFragmentStore, AnnotationTree<TSpec> > const &it) { 
	return it._id;
}

//////////////////////////////////////////////////////////////////////////////

template <typename TFragmentStore, typename TSpec>
inline TFragmentStore & 
container(Iter< TFragmentStore, AnnotationTree<TSpec> > &it) { 
	return *it.store;
}

template <typename TFragmentStore, typename TSpec>
inline TFragmentStore & 
container(Iter< TFragmentStore, AnnotationTree<TSpec> > const &it) { 
	return *it.store;
}

//////////////////////////////////////////////////////////////////////////////

/**
.Function.getAnnotation
..class:Spec.AnnotationTree Iterator
..summary:Returns the current annotation.
..cat:Fragment Store
..signature:getAnnotation(iter)
..param.iter:An annotation tree iterator.
...type:Spec.AnnotationTree Iterator
..returns:A reference to the @Class.AnnotationStoreElement@ the iterator points at.
..include:seqan/store.h
*/

template <typename TFragmentStore, typename TSpec>
inline typename GetValue<typename TFragmentStore::TAnnotationStore>::Type
getAnnotation(Iter<TFragmentStore const, AnnotationTree<TSpec> > const & it)
{
	return getValue(it.store->annotationStore, it._id);
}

template <typename TFragmentStore, typename TSpec>
inline typename Reference<typename TFragmentStore::TAnnotationStore>::Type
getAnnotation(Iter<TFragmentStore, AnnotationTree<TSpec> > const & it)
{
	return value(it.store->annotationStore, it._id);
}

/**
.Function.AnnotationTree Iterator#getName
..class:Spec.AnnotationTree Iterator
..summary:Returns the identifier of the current annotation.
..cat:Fragment Store
..signature:getName(iter)
..param.iter:An annotation tree iterator.
...type:Spec.AnnotationTree Iterator
..returns:The name of the current annotation.
...remarks:This a reference to the corresponding position in the @Memvar.FragmentStore#annotationNameStore@.
...type:Concept.SequenceConcept
..see:Function.AnnotationTree Iterator#setName
..include:seqan/store.h
*/

template <typename TFragmentStore, typename TSpec>
inline typename GetValue<typename TFragmentStore::TAnnotationNameStore>::Type
getName(Iter<TFragmentStore, AnnotationTree<TSpec> > const & it)
{
	return getAnnoName(*it.store, it._id);
}

/**
.Function.AnnotationTree Iterator#setName
..class:Spec.AnnotationTree Iterator
..summary:Sets the identifier of the current annotation.
..cat:Fragment Store
..signature:setName(iter, name)
..param.iter:An annotation tree iterator.
...type:Spec.AnnotationTree Iterator
..param.name:The new identifier of the current annotation element.
...type:Concept.SequenceConcept
..see:Function.AnnotationTree Iterator#getName
..include:seqan/store.h
*/

template <typename TFragmentStore, typename TSpec, typename TName>
inline void
setName(Iter<TFragmentStore, AnnotationTree<TSpec> > const & it, TName const & name)
{
	if (length(it.store->annotationNameStore) <= it._id)
		resize(it.store->annotationNameStore, it._id + 1);
	it.store->annotationNameStore[it._id] = name;
}

/**
.Function.getParentName
..class:Spec.AnnotationTree Iterator
..summary:Returns the identifier of the parent node in the annotation tree of the current annotation.
..cat:Fragment Store
..signature:getParentName(iter)
..param.iter:An annotation tree iterator.
...type:Spec.AnnotationTree Iterator
..returns:The name of the parent of the current annotation.
...remarks:This a reference to the corresponding position in the @Memvar.FragmentStore#annotationNameStore@.
...type:Concept.SequenceConcept
..see:Function.AnnotationTree Iterator#getName
..include:seqan/store.h
*/

template <typename TFragmentStore, typename TSpec>
inline typename GetValue<typename TFragmentStore::TAnnotationNameStore>::Type
getParentName(Iter<TFragmentStore, AnnotationTree<TSpec> > & it)
{
	typedef typename TFragmentStore::TAnnotationStore	TAnnotationStore;
	typedef typename Value<TAnnotationStore>::Type		TAnnotation;
	typedef typename TAnnotation::TId					TId;

	TId parentId = it.store->annotationStore[it._id].parentId;
	if (parentId == TAnnotation::INVALID_ID) parentId = it._id;
	return getAnnoName(*it.store, parentId);
}

//////////////////////////////////////////////////////////////////////////////
/**
.Function.getType
..class:Spec.AnnotationTree Iterator
..summary:Returns the type name of the current annotation.
..cat:Fragment Store
..signature:getType(iter)
..param.iter:An annotation tree iterator.
...type:Spec.AnnotationTree Iterator
..returns:The type name of the current annotation, e.g. "exon" or "mRNA".
...remarks:This a reference to an entry in the @Memvar.FragmentStore#annotationTypeStore@.
...type:Concept.SequenceConcept
..see:Function.setType
..include:seqan/store.h
*/

template <typename TFragmentStore, typename TSpec>
inline typename GetValue<typename TFragmentStore::TAnnotationTypeStore>::Type
getType(Iter<TFragmentStore, AnnotationTree<TSpec> > const & it)
{
	return getAnnoType(*it.store, getAnnotation(it).typeId);
}

/**
.Function.setType
..class:Spec.AnnotationTree Iterator
..summary:Sets the type name of the current annotation.
..cat:Fragment Store
..signature:setType(iter, typeName)
..param.iter:An annotation tree iterator.
...type:Spec.AnnotationTree Iterator
..param.typeName:The new type name, e.g. "exon" or "mRNA".
...type:Concept.SequenceConcept
..see:Function.getType
..include:seqan/store.h
*/

template <typename TFragmentStore, typename TSpec, typename TTypeName>
inline void
setType(Iter<TFragmentStore, AnnotationTree<TSpec> > const & it, TTypeName & typeName)
{
	_storeAppendType(*it.store, getAnnotation(it).typeId, typeName);
}

//////////////////////////////////////////////////////////////////////////////
/**
.Function.getUniqueName
..class:Spec.AnnotationTree Iterator
..summary:Returns a unique name of the current annotation.
..cat:Fragment Store
..signature:getUniqueName(iter)
..param.iter:An annotation tree iterator.
...type:Spec.AnnotationTree Iterator
..returns:A unique name of the current annotation.
...type:Shortcut.CharString
..remarks:Some annotation file formats doesn't require that every annotation has a non-empty name.
This function returns the name if non-empty and otherwise generates one using the type and id.
..see:Function.AnnotationTree Iterator#getName
..include:seqan/store.h
*/

template <typename TFragmentStore, typename TSpec>
inline CharString
getUniqueName(Iter<TFragmentStore, AnnotationTree<TSpec> > const & it)
{
	return getAnnoUniqueName(*it.store, it._id);
}

//////////////////////////////////////////////////////////////////////////////
/**
.Function.clearValues
..class:Spec.AnnotationTree Iterator
..summary:Clear all key-value pairs of the current annotation.
..cat:Fragment Store
..signature:clearValues(iter)
..param.iter:An annotation tree iterator.
...type:Spec.AnnotationTree Iterator
..see:Function.assignValueByKey
..see:Function.getValueByKey
..include:seqan/store.h
*/

template <typename TFragmentStore, typename TSpec>
inline void 
clearValues(Iter<TFragmentStore, AnnotationTree<TSpec> > & it)
{
	clear(getAnnotation(it).values);
}

/**
.Function.assignValueByKey
..class:Spec.AnnotationTree Iterator
..summary:Add or update a key-value pair of the current annotation.
..cat:Fragment Store
..signature:assignValueByKey(iter, key, value)
..param.iter:An annotation tree iterator.
...type:Spec.AnnotationTree Iterator
..param.key:The key whose value should be changed.
...type:Concept.SequenceConcept
..param.value:The new value of the key.
...type:Concept.SequenceConcept
..see:Function.clearValues
..see:Function.getValueByKey
..include:seqan/store.h
*/

template <typename TFragmentStore, typename TSpec, typename TKey, typename TValue>
inline void 
assignValueByKey(
	Iter<TFragmentStore, AnnotationTree<TSpec> > & it,
	TKey const & key,
	TValue const & value)
{
	annotationAssignValueByKey(*it.store, getAnnotation(it), key, value);
}

/**
.Function.getValueByKey
..class:Spec.AnnotationTree Iterator
..summary:Given a key, retrieve its value of the current annotation.
..cat:Fragment Store
..signature:getValueByKey(iter, key)
..signature:getValueByKey(iter, key, value)
..param.iter:An annotation tree iterator.
...type:Spec.AnnotationTree Iterator
..param.key:The key to retrieve the value from.
...type:Concept.SequenceConcept
..param.value:The sequence where the return value should be wrote to.
...type:Concept.SequenceConcept
..returns:The three-parameter function returns a boolean indicating whether a non-empty value was returned.
..returns:The two-parameter function returns the value for the given key.
...type:nolink:bool
...type:Shortcut.CharString
..see:Function.clearValues
..see:Function.assignValueByKey
..include:seqan/store.h
*/

template <typename TFragmentStore, typename TSpec, typename TKey, typename TValue>
inline bool 
getValueByKey(
	Iter<TFragmentStore, AnnotationTree<TSpec> > const & it,
	TKey const & key,
	TValue & value)
{
	return annotationGetValueByKey(*it.store, getAnnotation(it), key, value);
}

template <typename TFragmentStore, typename TSpec, typename TKey>
inline CharString
getValueByKey(
	Iter<TFragmentStore, AnnotationTree<TSpec> > const & it,
	TKey const & key)
{
	return annotationGetValueByKey(*it.store, getAnnotation(it), key);
}

//////////////////////////////////////////////////////////////////////////////

///.Function.goBegin.param.iterator.type:Spec.AnnotationTree Iterator
template <typename TFragmentStore, typename TSpec>
inline void
goBegin(Iter<TFragmentStore, AnnotationTree<TSpec> > & it)
{
	it._id = 0;
}

///.Function.goEnd.param.iterator.type:Spec.AnnotationTree Iterator
template <typename TFragmentStore, typename TSpec>
inline void
goEnd(Iter<TFragmentStore, AnnotationTree<TSpec> > & it)
{
	typedef typename TFragmentStore::TAnnotationStore	TAnnotationStore;
	typedef typename Value<TAnnotationStore>::Type		TAnnotation;

	it._id = TAnnotation::INVALID_ID;
}

/**
.Function.goTo
..class:Spec.AnnotationTree Iterator
..summary:Moves the iterator to an arbitrary node given its annotationId.
..cat:Fragment Store
..signature:goTo(iter, annotationId)
..param.iter:An annotation tree iterator.
...type:Spec.AnnotationTree Iterator
..param.annotationId:The id of the new annotation.
..returns: Iterator to the new node.
...type:Spec.AnnotationTree Iterator
..include:seqan/store.h
*/

template <typename TFragmentStore, typename TSpec, typename TId>
inline void
goTo(Iter<TFragmentStore, AnnotationTree<TSpec> > & it, TId _id)
{
	it._id = _id;
}

template <typename TFragmentStore, typename TSpec>
inline void
clear(Iter<TFragmentStore, AnnotationTree<TSpec> > & it)
{
	typedef typename TFragmentStore::TAnnotationStore	TAnnotationStore;
	typedef typename Value<TAnnotationStore>::Type		TAnnotation;

	it._id = TAnnotation::INVALID_ID;
}

//////////////////////////////////////////////////////////////////////////////

///.Function.atBegin.param.iterator.type:Spec.AnnotationTree Iterator
template <typename TFragmentStore, typename TSpec>
inline bool
atBegin(Iter<TFragmentStore, AnnotationTree<TSpec> > & it)
{
	return it._id == 0;
}

///.Function.atEnd.param.iterator.type:Spec.AnnotationTree Iterator
template <typename TFragmentStore, typename TSpec>
inline bool
atEnd(Iter<TFragmentStore, AnnotationTree<TSpec> > & it)
{
	typedef typename TFragmentStore::TAnnotationStore	TAnnotationStore;
	typedef typename Value<TAnnotationStore>::Type		TAnnotation;

	return it._id == TAnnotation::INVALID_ID;
}

//////////////////////////////////////////////////////////////////////////////

///.Function.goNext.param.iterator.type:Spec.AnnotationTree Iterator
template <typename TFragmentStore, typename TSpec>
inline void
goNext(Iter<TFragmentStore, AnnotationTree<TSpec> > & it)
{
	// preorder dfs
	if (!goDown(it) && !goRight(it))
		while (goUp(it) && !goRight(it)) ;
	if (isRoot(it)) {
		clear(it);
		return;
	}
}

/**
.Function.AnnotationTree Iterator#goNextRight
..class:Spec.AnnotationTree Iterator
..summary:Go to the next node in preorder DFS skipping the current node's subtree.
..cat:Fragment Store
..signature:goNextRight(iter)
..param.iter:An annotation tree iterator.
...type:Spec.AnnotationTree Iterator
..see:Function.goNext
..see:Function.AnnotationTree Iterator#goNextUp
..include:seqan/store.h
*/

template <typename TFragmentStore, typename TSpec>
inline void
goNextRight(Iter<TFragmentStore, AnnotationTree<TSpec> > & it)
{
	// preorder dfs
	if (!goRight(it))
		while (goUp(it) && !goRight(it)) ;
	if (isRoot(it)) {
		clear(it);
		return;
	}
}

/**
.Function.AnnotationTree Iterator#goNextUp
..class:Spec.AnnotationTree Iterator
..summary:Go to the next node in preorder DFS skipping the subtrees of the current node and of all its siblings.
..cat:Fragment Store
..signature:goNextUp(iter)
..param.iter:An annotation tree iterator.
...type:Spec.AnnotationTree Iterator
..see:Function.goNext
..see:Function.AnnotationTree Iterator#goNextRight
..include:seqan/store.h
*/

template <typename TFragmentStore, typename TSpec>
inline void
goNextUp(Iter<TFragmentStore, AnnotationTree<TSpec> > & it)
{
	// preorder dfs
	while (goUp(it) && !goRight(it)) ;
	if (isRoot(it)) {
		clear(it);
		return;
	}
}

//////////////////////////////////////////////////////////////////////////////

/**
.Function.AnnotationTree Iterator#goRoot
..class:Spec.AnnotationTree Iterator
..summary:Go to the root node in the annotation tree.
..cat:Fragment Store
..signature:goRoot(iter)
..param.iter:An annotation tree iterator.
...type:Spec.AnnotationTree Iterator
..see:Function.goNext
..see:Function.AnnotationTree Iterator#goNextRight
..include:seqan/store.h
*/

template <typename TFragmentStore, typename TSpec>
inline void
goRoot(Iter<TFragmentStore, AnnotationTree<TSpec> > & it)
{
	it._id = 0;
}

/**
.Function.AnnotationTree Iterator#goUp
..class:Spec.AnnotationTree Iterator
..summary:Move the iterator up in the annotation tree.
..cat:Fragment Store
..signature:goUp(iter)
..param.iter:An annotation tree iterator.
...type:Spec.AnnotationTree Iterator
..returns:$true$ if the iterator could be moved, otherwise $false$.
...type:nolink:bool
..see:Function.AnnotationTree Iterator#goDown
..see:Function.AnnotationTree Iterator#goRight
..include:seqan/store.h
*/

template <typename TFragmentStore, typename TSpec>
inline bool
goUp(Iter<TFragmentStore, AnnotationTree<TSpec> > & it)
{
	typedef typename TFragmentStore::TAnnotationStore	TAnnotationStore;
	typedef typename Value<TAnnotationStore>::Type		TAnnotation;
	typedef typename TAnnotation::TId					TId;
	
	TId parentId = getAnnotation(it).parentId;
	if (parentId != TAnnotation::INVALID_ID)
	{
		it._id = parentId;
		return true;
	}
	return false;
}

/**
.Function.AnnotationTree Iterator#goDown
..class:Spec.AnnotationTree Iterator
..summary:Move the iterator down to the left-most child in the annotation tree.
..cat:Fragment Store
..signature:goDown(iter)
..param.iter:An annotation tree iterator.
...type:Spec.AnnotationTree Iterator
..returns:$true$ if the iterator could be moved, otherwise $false$.
...type:nolink:bool
..see:Function.AnnotationTree Iterator#goUp
..see:Function.AnnotationTree Iterator#goRight
..include:seqan/store.h
*/

template <typename TFragmentStore, typename TSpec>
inline bool
goDown(Iter<TFragmentStore, AnnotationTree<TSpec> > & it)
{
	typedef typename TFragmentStore::TAnnotationStore	TAnnotationStore;
	typedef typename Value<TAnnotationStore>::Type		TAnnotation;
	typedef typename TAnnotation::TId					TId;
	
	TId lastChildId = getAnnotation(it).lastChildId;
	if (lastChildId != TAnnotation::INVALID_ID)
	{
		it._id = it.store->annotationStore[lastChildId].nextSiblingId;
		return true;
	}
	return false;
}

/**
.Function.AnnotationTree Iterator#goRight
..class:Spec.AnnotationTree Iterator
..summary:Move the iterator right to the next sibling in the annotation tree.
..cat:Fragment Store
..signature:goRight(iter)
..param.iter:An annotation tree iterator.
...type:Spec.AnnotationTree Iterator
..returns:$true$ if the iterator could be moved, otherwise $false$.
...type:nolink:bool
..see:Function.AnnotationTree Iterator#goUp
..see:Function.AnnotationTree Iterator#goDown
..include:seqan/store.h
*/

template <typename TFragmentStore, typename TSpec>
inline bool
goRight(Iter<TFragmentStore, AnnotationTree<TSpec> > & it)
{
	typedef typename TFragmentStore::TAnnotationStore	TAnnotationStore;
	typedef typename Value<TAnnotationStore>::Type		TAnnotation;
	typedef typename TAnnotation::TId					TId;
	
	TAnnotation const &anno = getAnnotation(it);
	TId nextSiblingId = anno.nextSiblingId;
	if (nextSiblingId != TAnnotation::INVALID_ID)
	{
		TId lastChildId = it.store->annotationStore[anno.parentId].lastChildId;
		if (it._id != lastChildId)
		{
			it._id = nextSiblingId;
			return true;
		}
	}
	return false;
}

//////////////////////////////////////////////////////////////////////////////
/**
.Function.AnnotationTree Iterator#nodeUp
..class:Spec.AnnotationTree Iterator
..summary:Returns a new iterator to the parent node of the current annotation in the annotation tree.  
..cat:Fragment Store
..signature:nodeUp(iter)
..param.iter:An annotation tree iterator.
...type:Spec.AnnotationTree Iterator
..returns: A new iterator to the parent node.
....type:Spec.AnnotationTree Iterator
..see:Function.AnnotationTree Iterator#nodeDown
..see:Function.AnnotationTree Iterator#nodeRight
..include:seqan/store.h
*/

template <typename TFragmentStore, typename TSpec>
inline Iter<TFragmentStore, AnnotationTree<TSpec> >
nodeUp(Iter<TFragmentStore, AnnotationTree<TSpec> > const & it)
{
	Iter<TFragmentStore, AnnotationTree<TSpec> > tmp(it);
	goUp(tmp);
	return tmp;
}

/**
.Function.AnnotationTree Iterator#nodeDown
..class:Spec.AnnotationTree Iterator
..summary:Returns a new iterator to the first child node of the current annotation in the annotation tree.  
..cat:Fragment Store
..signature:nodeDown(iter)
..param.iter:An annotation tree iterator.
...type:Spec.AnnotationTree Iterator
..returns: A new iterator to the parent node.
....type:Spec.AnnotationTree Iterator
..see:Function.AnnotationTree Iterator#nodeUp
..see:Function.AnnotationTree Iterator#nodeRight
..include:seqan/store.h
*/

template <typename TFragmentStore, typename TSpec>
inline Iter<TFragmentStore, AnnotationTree<TSpec> >
nodeDown(Iter<TFragmentStore, AnnotationTree<TSpec> > const & it)
{
	Iter<TFragmentStore, AnnotationTree<TSpec> > tmp(it);
	goDown(tmp);
	return tmp;
}

/**
.Function.AnnotationTree Iterator#nodeRight
..class:Spec.AnnotationTree Iterator
..summary:Returns a new iterator to the right sibling of the current annotation in the annotation tree.  
..cat:Fragment Store
..signature:nodeRight(iter)
..param.iter:An annotation tree iterator.
...type:Spec.AnnotationTree Iterator
..returns: A new iterator to the right sibling.
....type:Spec.AnnotationTree Iterator
..see:Function.AnnotationTree Iterator#nodeUp
..see:Function.AnnotationTree Iterator#nodeDown
..include:seqan/store.h
*/

template <typename TFragmentStore, typename TSpec>
inline Iter<TFragmentStore, AnnotationTree<TSpec> >
nodeRight(Iter<TFragmentStore, AnnotationTree<TSpec> > const & it)
{
	Iter<TFragmentStore, AnnotationTree<TSpec> > tmp(it);
	goRight(tmp);
	return tmp;
}

//////////////////////////////////////////////////////////////////////////////

// insert a new _id into a cyclic list and returns new last child _id
template <typename TAnnotationStore, typename TId>
inline TId
_cyclicListFrontInsert(TAnnotationStore & annotationStore, TId newId, TId lastChildId)
{
	typedef typename Value<TAnnotationStore>::Type TAnnotation;

	TId nextId, newLastId;
	if (lastChildId != TAnnotation::INVALID_ID)
	{
		// get last node in the cycle
		TAnnotation &lastChild = annotationStore[lastChildId];
		// last child points to first child
		nextId = lastChild.nextSiblingId;
		// insert new node between last and first
		lastChild.nextSiblingId = newId;
		// last child remains the same
		newLastId = lastChildId;
	} else
		// cyclic list was empty
		newLastId = nextId = newId;
	
	// link new node to former first node
	annotationStore[newId].nextSiblingId = nextId;
	
	return newLastId;
}

// delete an _id from a cyclic list and returns new last child _id
template <typename TAnnotationStore, typename TId>
inline TId
_cyclicListSearchPrev(TAnnotationStore & annotationStore, TId _id, TId lastChildId)
{
	typedef typename Value<TAnnotationStore>::Type TAnnotation;

	if (lastChildId == TAnnotation::INVALID_ID)
		return TAnnotation::INVALID_ID;
	
	TId prevId, i = lastChildId;
	do {
		prevId = i;
		i = annotationStore[i].nextSiblingId;
		if (i == _id) break;
	} while (i != lastChildId);

	if (i == _id)
		return prevId;
	else
		return TAnnotation::INVALID_ID;
}

// delete an _id from a cyclic list and returns new last child _id
template <typename TAnnotationStore, typename TId>
inline TId
_cyclicListRemove(TAnnotationStore & annotationStore, TId _id, TId lastChildId)
{
	typedef typename Value<TAnnotationStore>::Type TAnnotation;

	TId prevId = _cyclicListSearchPrev(annotationStore, _id, lastChildId);
	
	if (prevId != TAnnotation::INVALID_ID)
	{
		annotationStore[prevId].nextSiblingId = annotationStore[_id].nextSiblingId;
		
		if (_id == lastChildId)
		{
			if (prevId != _id)
				return prevId;
			else
				return TAnnotation::INVALID_ID;
		} else
			return lastChildId;
	}
	return lastChildId;
}

/**
.Function.createLeftChild
..class:Spec.AnnotationTree Iterator
..summary:Creates a new left-most child of the current node and returns an iterator to it.
..cat:Fragment Store
..signature:createLeftChild(iter)
..param.iter:An annotation tree iterator.
...type:Spec.AnnotationTree Iterator
..returns: Iterator to the new left-most child.
....type:Spec.AnnotationTree Iterator
..see:Function.createRightChild
..see:Function.createSibling
..include:seqan/store.h
*/

template <typename TFragmentStore, typename TSpec>
inline Iter<TFragmentStore, AnnotationTree<TSpec> >
createLeftChild(Iter<TFragmentStore, AnnotationTree<TSpec> > & it)
{
	typedef typename TFragmentStore::TAnnotationStore	TAnnotationStore;
	typedef typename Value<TAnnotationStore>::Type		TAnnotation;
	typedef typename TAnnotation::TId					TId;
	
	appendValue(it.store->annotationStore, getAnnotation(it));
	TAnnotation &parentAnno = getAnnotation(it);

	TId childId = length(it.store->annotationStore) - 1;
	TAnnotation &childAnno = it.store->annotationStore[childId];
	
	parentAnno.lastChildId = _cyclicListFrontInsert(it.store->annotationStore, childId, parentAnno.lastChildId);
	childAnno.parentId = it._id;
	childAnno.lastChildId = TAnnotation::INVALID_ID;
	
	Iter<TFragmentStore, AnnotationTree<TSpec> > childIter(it);
	childIter._id = childId;
	return childIter;
}

/**
.Function.createRightChild
..class:Spec.AnnotationTree Iterator
..summary:Creates a new right-most child of the current node and returns an iterator to it.
..cat:Fragment Store
..signature:createRightChild(iter)
..param.iter:An annotation tree iterator.
...type:Spec.AnnotationTree Iterator
..returns: Iterator to the new right-most child.
....type:Spec.AnnotationTree Iterator
..see:Function.createLeftChild
..see:Function.createSibling
..include:seqan/store.h
*/

template <typename TFragmentStore, typename TSpec>
inline Iter<TFragmentStore, AnnotationTree<TSpec> >
createRightChild(Iter<TFragmentStore, AnnotationTree<TSpec> > & it)
{
	typedef typename TFragmentStore::TAnnotationStore	TAnnotationStore;
	typedef typename Value<TAnnotationStore>::Type		TAnnotation;
	typedef typename TAnnotation::TId					TId;
	
	appendValue(it.store->annotationStore, getAnnotation(it));
	TAnnotation &parentAnno = getAnnotation(it);

	TId childId = length(it.store->annotationStore) - 1;
	TAnnotation &childAnno = it.store->annotationStore[childId];
	
	_cyclicListFrontInsert(it.store->annotationStore, childId, parentAnno.lastChildId);
	parentAnno.lastChildId = childId;
	childAnno.parentId = it._id;
	childAnno.lastChildId = TAnnotation::INVALID_ID;
	
	Iter<TFragmentStore, AnnotationTree<TSpec> > childIter(it);
	childIter._id = childId;
	return childIter;
}

/**
.Function.createSibling
..class:Spec.AnnotationTree Iterator
..summary:Creates a new right sibling of the current node and returns an iterator to it.
..cat:Fragment Store
..signature:createSibling(iter)
..param.iter:An annotation tree iterator.
...type:Spec.AnnotationTree Iterator
..returns: Iterator to the new right sibling.
....type:Spec.AnnotationTree Iterator
..see:Function.createLeftChild
..see:Function.createRightChild
..include:seqan/store.h
*/

template <typename TFragmentStore, typename TSpec>
inline Iter<TFragmentStore, AnnotationTree<TSpec> >
createSibling(Iter<TFragmentStore, AnnotationTree<TSpec> > & it)
{
	typedef typename TFragmentStore::TAnnotationStore	TAnnotationStore;
	typedef typename Value<TAnnotationStore>::Type		TAnnotation;
	typedef typename TAnnotation::TId					TId;

	appendValue(it.store->annotationStore, getAnnotation(it));
	TAnnotation &anno = getAnnotation(it);

	TId siblingId = length(it.store->annotationStore) - 1;

	TAnnotation &parentAnno = it.store->annotationStore[anno.parentId];
	if (parentAnno.lastChildId == it._id)
		parentAnno.lastChildId = siblingId;

	TAnnotation &siblingAnno = it.store->annotationStore[siblingId];
	siblingAnno.nextSiblingId = anno.nextSiblingId;
	siblingAnno.parentId = anno.parentId;
	siblingAnno.lastChildId = TAnnotation::INVALID_ID;
	anno.nextSiblingId = siblingId;
	
	Iter<TFragmentStore, AnnotationTree<TSpec> > siblingIter(it);
	siblingIter._id = siblingId;
	return siblingIter;
}

//////////////////////////////////////////////////////////////////////////////

/**
.Function.isRoot:
..class:Spec.AnnotationTree Iterator
..param.iterator.type:Spec.AnnotationTree Iterator
*/

template <typename TFragmentStore, typename TSpec>
inline bool
isRoot(Iter<TFragmentStore, AnnotationTree<TSpec> > const & it)
{
	typedef typename TFragmentStore::TAnnotationStore	TAnnotationStore;
	typedef typename Value<TAnnotationStore>::Type		TAnnotation;

//	if (it._id >= length(it.store->annotationStore)) return false;
	return it.store->annotationStore[it._id].parentId == TAnnotation::INVALID_ID;
}

/**
.Function.isLeaf:
..class:Spec.AnnotationTree Iterator
..param.iterator.type:Spec.AnnotationTree Iterator
*/

template <typename TFragmentStore, typename TSpec>
inline bool
isLeaf(Iter<TFragmentStore, AnnotationTree<TSpec> > const & it)
{
	typedef typename TFragmentStore::TAnnotationStore	TAnnotationStore;
	typedef typename Value<TAnnotationStore>::Type		TAnnotation;

//	if (it._id >= length(it.store->annotationStore)) return false;
	return it.store->annotationStore[it._id].lastChildId == TAnnotation::INVALID_ID;
}

/**
.Function.isLastChild
..class:Spec.AnnotationTree Iterator
..summary:Returns a boolean value that indicates whether the current node is the last child.
..cat:Fragment Store
..signature:isLastChild(iter)
..param.iter:An annotation tree iterator.
...type:Spec.AnnotationTree Iterator
..returns: $true$ if the iterator is the last child, otherwise $false$.
...type:nolink:bool
..include:seqan/store.h
*/

template <typename TFragmentStore, typename TSpec>
inline bool
isLastChild(Iter<TFragmentStore, AnnotationTree<TSpec> > const & it)
{
	typedef typename TFragmentStore::TAnnotationStore	TAnnotationStore;
	typedef typename Value<TAnnotationStore>::Type		TAnnotation;
	typedef typename TAnnotation::TId					TId;

	TAnnotation const &anno = getAnnotation(it);
	TId nextSiblingId = anno.nextSiblingId;
	if (nextSiblingId != TAnnotation::INVALID_ID)
	{
		TId lastChildId = it.store->annotationStore[anno.parentId].lastChildId;
		return it._id == lastChildId;
	}
	return true;
}

//////////////////////////////////////////////////////////////////////////////

template <typename TAnnotationStore>
inline void
_storeClearAnnoBackLinks(TAnnotationStore & me)
{
	typedef typename Value<TAnnotationStore>::Type				TAnnotation;
	typedef typename Iterator<TAnnotationStore, Standard>::Type TAnnoIter;

	TAnnoIter it = begin(me, Standard());
	TAnnoIter itEnd = end(me, Standard());
	
	for (; it != itEnd; ++it)
	{
		(*it).lastChildId = TAnnotation::INVALID_ID;
		(*it).nextSiblingId = TAnnotation::INVALID_ID;
	}
}

template <typename TAnnotationStore>
inline void
_storeCreateAnnoBackLinks(TAnnotationStore & me)
{
	typedef typename Value<TAnnotationStore>::Type				TAnnotation;
	typedef typename TAnnotation::TId							TId;
	typedef typename Iterator<TAnnotationStore, Standard>::Type TAnnoIter;
	
	TAnnoIter itBegin = begin(me, Standard());
	TAnnoIter itEnd = end(me, Standard());
	TId _id = (itEnd - itBegin) - 1;
	TAnnoIter it = itBegin + _id;
	
	for (; itBegin <= it; --it, --_id)
	{
		if ((*it).parentId != TAnnotation::INVALID_ID)
		{
			TAnnoIter parent = itBegin + (*it).parentId;
			if ((*parent).lastChildId == TAnnotation::INVALID_ID)
			{
				(*parent).lastChildId = _id;
				(*it).nextSiblingId = _id;
			}

			if ((*it).nextSiblingId == TAnnotation::INVALID_ID)
			{
				TAnnoIter lastChild = itBegin + (*parent).lastChildId;
				(*it).nextSiblingId = (*lastChild).nextSiblingId;
				(*lastChild).nextSiblingId = _id;
			}
		}
		else 
			(*it).nextSiblingId = TAnnotation::INVALID_ID;
	}
}

template <typename TPos, typename TSpec>
inline std::ostream &
operator << (std::ostream & out, AnnotationStoreElement<TPos, TSpec> const & anno)
{
    out << "parentId:     \t" << anno.parentId << std::endl;
    out << "contigId:     \t" << anno.contigId << std::endl;
    out << "countId:      \t" << anno.countId << std::endl;
    out << "typeId:       \t" << anno.typeId << std::endl;
    out << "beginPos:     \t" << anno.beginPos << std::endl;
    out << "endPos:       \t" << anno.endPos << std::endl;
    out << "lastChildId:  \t" << anno.lastChildId << std::endl;
    out << "nextSiblingId:\t" << anno.nextSiblingId << std::endl;
    
    return out;
}

}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
