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
// Author: Andreas Gogol-Döring <andreas.doering@mdc-berlin.de>
// ==========================================================================
// Definition of basic Metafunctions.
// ==========================================================================

// TODO(holtgrew): Rename to "shared metafunctions.h"?
// TODO(holtgrew): This could use some cleanup.

#ifndef SEQAN_BASIC_BASIC_TYPE_H_
#define SEQAN_BASIC_BASIC_TYPE_H_

namespace seqan {

//////////////////////////////////////////////////////////////////////////////


// TODO(holtgrew): Remove default implementation.
// SEQREV: elements-are-containers should not have a default implementation

template <typename T, const int i>
struct Value
{
    typedef T Type;
};

template <typename T>
struct Value<T const>
{
    typedef T const Type;
};

//____________________________________________________________________________


// TODO(holtgrew): Remove default implementation.
template <typename T>
struct GetValue
{
    typedef typename Value<T>::Type const & Type;
};

template <typename T>
struct GetValue<T const>:
    public GetValue<T>
{
};

//____________________________________________________________________________

// TODO(holtgrew): Remove default implementation.
template <typename T>
struct Reference
{
    typedef typename Value<T>::Type & Type;
};

template <typename T>
struct Reference<T const>
{
    typedef typename Value<T>::Type const & Type;
};

//____________________________________________________________________________


// TODO(holtgrew): Remove default implementation.
template <typename T>
struct Size
{
    typedef size_t Type;
};

template <typename T>
struct Size<T const>:
    Size<T>
{
};

//____________________________________________________________________________



// TODO(holtgrew): Remove default implementation.
template <typename T>
struct Difference
{
    typedef ptrdiff_t Type;
};

template <typename T>
struct Difference<T const>:
    Difference<T>
{
};

//____________________________________________________________________________


// TODO(holtgrew): Remove default implementation.
template <typename T>
struct Position
{
    typedef typename Size<T>::Type Type;
};

template <typename T>
struct Position<T const>:
    Position<T>
{
};

//____________________________________________________________________________

// SEQREV: elements-are-containers should not have a default implementation, if a type has no host, do not return self

template <typename T>
struct Host
{
    typedef T Type;
};

//____________________________________________________________________________

//____________________________________________________________________________

/**
.Metafunction.Cargo:
..cat:Basic
..summary:Type of additional data stored in an object. 
..signature:Cargo<T>::Type
..param.T:Type for which the cargo tyoe is determined.
..returns.param.Type:Cargo of $T$.
..remarks:The definition of Cargo allows the addition of user specific data to existing data structures.
..include:seqan/basic.h
*/

// TODO(holtgrew): Should this have a default implementation?

template <typename T>
struct Cargo {
    typedef Nothing Type;
};
template <typename T>
struct Cargo<T const> {
    typedef typename Cargo<T>::Type const Type;
};

//____________________________________________________________________________

/**
.Metafunction.VertexDescriptor:
..cat:Graph
..summary:Type of an object that represents a vertex descriptor.
..signature:VertexDescriptor<T>::Type
..param.T:Type T must be a graph. All graphs currently use ids as vertex descriptors.
..returns.param.Type:VertexDescriptor type.
..remarks.text:The vertex descriptor is a unique handle to a vertex in a graph.
It is used in various graph functions, e.g., to add edges, to create OutEdge Iterators or to remove a vertex.
It is also used to attach properties to vertices.
..example.code:VertexDescriptor<Graph<> >::Type vD; //vD is a vertex descriptor
..include:seqan/basic.h
*/

// TODO(holtgrew): Should this have a default implementation? For all graphs -- OK but for all types?

template <typename T>
struct VertexDescriptor {
    typedef void* Type;
};
template <typename T>
struct VertexDescriptor<T const>:
    public VertexDescriptor<T> {};


//____________________________________________________________________________

    
/**
.Metafunction.Id:
..cat:Graph
..summary:Type of an object that represents an id.
..signature:Id<T>::Type
..param.T:Type for which a suitable id type is determined.
..returns.param.Type:Id type.
..remarks.text:The id type of a container is the type that is used to uniquely identify its elements.
In most cases this type is unsigned int.
..example.code:Id<Graph<> >::Type id; //id has type unsigned int
..include:seqan/basic.h
*/

// TODO(holtgrew): Should this have a default implementation?

template<typename T>
struct Id {
    typedef unsigned int Type;
};

template<typename T>
struct Id<T const> : Id<T> {};

//____________________________________________________________________________

/**
.Metafunction.Key:
..cat:Graph
..summary:Key type of a key to cargo mapping.
..signature:Key<T>::Type
..param.T:Type for which a key type is determined.
..returns.param.Type:Key type.
...default:The type $T$ itself.
..include:seqan/basic.h
*/

// TODO(holtgrew): Should this have a default implementation?

template< typename T >
struct Key
{
    typedef T Type;
};

template <typename T>
struct Key<T const>:
    Key<T> {};

//____________________________________________________________________________

/*VERALTET
.Metafunction.Object:
..summary:Object type of a key to object mapping.
..signature:Object<T>::Type
..param.T:Type for which a object type is determined.
..returns.param.Type:Object type.
..include:seqan/basic.h
*/

template<typename T>
struct Object; 

template <typename T>
struct Object<T const>:
    Object<T> {};


//____________________________________________________________________________

// TODO(holtgrew): Move to alignments?
// TODO(holtgrew): Is this default implementation what we want?

/**
.Metafunction.Source
..cat:Alignments
*/

template < typename TSpec = void >
struct Source
{
    typedef TSpec Type;
};

template <typename T>
struct Source<T const>:
    Source<T>
{
};

//____________________________________________________________________________

/**
.Internal.Parameter_:
..cat:Metafunctions
..summary:Type for function parameters and return values.
..signature:Parameter_<T>::Type
..param.T:A type.
..returns.param.Type:The parameter type for arguments of type $T$.
...text:If $T$ is a pointer or array type, then $Parameter_<T>::Type$ is $T$, 
otherwise $Parameter_<T>::Type$ is $T &$.
*/

// TODO(holtgrew): Really required?

template <typename T>
struct Parameter_
{
    typedef T & Type;
};

template <typename T>
struct Parameter_<T *>
{
    typedef T * Type;
};

template <typename T, size_t I>
struct Parameter_<T [I]>
{
    typedef T * Type;
};

/**
.Internal._toParameter:
..cat:Functions
..summary:Transforms pointers to parameter types.
..signature:_toParameter<T>(pointer)
..param.pointer:A pointer.
..param.T:A Type.
...text:$object$ is transformed into the parameter type of $T$ that is given by @Internal.Parameter_@.
...note:This type must be explicitely specified.
..returns:To $TParameter$ transformed $object$.
..see:Internal.Parameter_
*/

// TODO(holtgrew): Really required?

template <typename T>
typename Parameter_<T>::Type
inline _toParameter(T * _object)
{
    return * _object;
}

template <typename T>
typename Parameter_<T>::Type
inline _toParameter(T & _object)
{
    return _object;
}

template <typename T>
typename Parameter_<T const>::Type
inline _toParameter(T const & _object)
{
    return _object;
}

//____________________________________________________________________________

/**
.Internal.ConstParameter_:
..cat:Metafunctions
..summary:Type for constant function parameters and return values.
..signature:ConstParameter_<T>::Type
..param.T:A type.
..returns.param.Type:The const parameter type for arguments of type $T$.
...text:If $T$ is a pointer or array type, then $Parameter_<T>::Type$ is a pointer to a const array, 
otherwise $Parameter_<T>::Type$ is $T const &$.
..see:Internal.Parameter_
*/

// TODO(holtgrew): Really required?

template <typename T>
struct ConstParameter_
{
    typedef T const & Type;
};

template <typename T>
struct ConstParameter_<T const>:
    public ConstParameter_<T> {};

template <typename T>
struct ConstParameter_<T *>
{
    typedef T const * Type;
};

template <typename T>
struct ConstParameter_<T const *>
{
    typedef T const * Type;
};

template <typename T, size_t I>
struct ConstParameter_<T [I]>
{
    typedef T const * Type;
};

template <typename T, size_t I>
struct ConstParameter_<T const [I]>
{
    typedef T const * Type;
};

//____________________________________________________________________________

/**
.Internal.Pointer_:
..cat:Metafunctions
..summary:The associated pointer type.
..signature:Pointer_<T>::Type
..param.T:A type.
..returns.param.Type:A pointer type for $T$.
...text:if $T$ is already a pointer type, then $Pointer_<T>::Type$ is $T$,
otherwise $Pointer_<T>::Type$ is $T *$.
..see:Internal.Parameter_
..see:Internal._toParameter
*/

// TODO(holtgrew): Really required?

template <typename T>
struct Pointer_
{
    typedef T * Type;
};

template <typename T>
struct Pointer_<T *>
{
    typedef T * Type;
};
template <typename T>
struct Pointer_<T * const>
{
    typedef T * Type;
};

template <typename T, size_t I>
struct Pointer_<T [I]>
{
    typedef T * Type;
};

//non const version of Pointer_ for return values

template <typename T>
struct NonConstPointer_:
    Pointer_<T>
{
};
template <typename T>
struct NonConstPointer_<T * const>
{
    typedef T * Type;
};

/**
.Internal._toPointer:
..cat:Functions
..summary:Transforms types into pointers.
..signature:_toPointer(object)
..param.object:An object.
..returns:$object$, transformed to a pointer. 
...text:The type of the returned pointer is given by @Internal.Pointer_@.
..see:Internal.Pointer_
*/

// TODO(holtgrew): Really required?

template <typename T>
typename NonConstPointer_<T>::Type
_toPointer(T & _object)
{
SEQAN_CHECKPOINT
    return & _object;
}
template <typename T>
typename NonConstPointer_<T const>::Type
_toPointer(T const & _object)
{
SEQAN_CHECKPOINT
    return & _object;
}

template <typename T>
typename NonConstPointer_<T *>::Type
_toPointer(T * _object)
{
SEQAN_CHECKPOINT
    return _object;
}

//____________________________________________________________________________


/**
.Metafunction.LENGTH:
..cat:Basic
..summary:Number of elements in a fixed-size container.
..signature:LENGTH<T>::VALUE
..param.T:Type for which the number of elements is determined.
..returns.param.VALUE:Number of elements.
..remarks.text:The default return value is 1 for dynamic-size containers.
..include:seqan/basic.h
*/

// SEQREV: elements-are-containers should probably not have a default implementation
// TODO(holtgrew): Rather switch to static const unsigned VALUE = ?

template <typename T>
struct LENGTH
{
    enum { VALUE = 1 };
};

template <typename T>
struct LENGTH<T const>:
    LENGTH<T>
{
};

/**
.Metafunction.WEIGHT:
..cat:Index
..summary:Number of relevant positions in a shape.
..signature:WEIGHT<T>::Type
..param.T:Shape type for which the number of relevant positions is determined.
...type:Class.Shape
..returns.param.VALUE:Number of relevant positions.
..remarks.text:The default return value is the result of the @Metafunction.LENGTH@ function.
For gapped shapes this is the number of '1's.
..include:seqan/basic.h
*/

// TODO(holtgrew): Should probably go to wherever shapes are defined.

template <typename T>
struct WEIGHT:
    LENGTH<T>
{
};
template <typename T>
struct WEIGHT<T const>:
    WEIGHT<T>
{
};

//////////////////////////////////////////////////////////////////////////////

//Iterator: see basic_iterator.h

//////////////////////////////////////////////////////////////////////////////


}  // namespace seqan

#endif  // #ifndef SEQAN_BASIC_BASIC_TYPE_H_
