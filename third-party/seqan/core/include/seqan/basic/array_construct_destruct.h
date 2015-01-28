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
// Functions for destructing and several ways of constructing (e.g. copy,
// move) values or arrays of values.
// ==========================================================================

// TODO(holtgrew): Order of parameters should be (target1, target2, ..., source1, source2, ...).
// TODO(holtgrew): Can we maybe replace at least part with http://www.cplusplus.com/reference/std/memory/?
// TODO(holtgrew): The metafunction should go into the alphabet submodule, the functions into the sequence/string module.

#include <new>

#ifndef SEQAN_CORE_INCLUDE_SEQAN_BASIC_ARRAY_CONSTRUCT_DESTRUCT_H_
#define SEQAN_CORE_INCLUDE_SEQAN_BASIC_ARRAY_CONSTRUCT_DESTRUCT_H_

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// ============================================================================
// Metafunctions
// ============================================================================

// ----------------------------------------------------------------------------
// Metafunction IsSimple
// ----------------------------------------------------------------------------

// TODO(holtgrew): This has to go to alphabet sub module, storage section, adaption.

/**
.Metafunction.IsSimple
..cat:Basic
..summary:Tests type to be simple.
..signature:IsSimple<T>::Type
..param.T:Type that is tested.
..returns.param.Type:@Tag.Logical Values.True@, if $T$ is a simple type, @Tag.Logical Values.False@ otherwise.
...default:@Tag.Logical Values.False@
..remarks:A simple type is a type that does not need constructors to be created,
a destructor to be destroyed, and copy assignment operators or copy constructors
to be copied. All POD ("plain old data") types are simple, but some
non-POD types could be simple too, e.g. some specializations of @Class.SimpleType@.
..see:Class.SimpleType
..include:seqan/basic.h
*/

template <typename T>
struct IsSimple_
{
    typedef False Type;
    enum { VALUE = 0 };
};

// ----------------------------------------------------------------------------
// Metafunction IsSimple
// ----------------------------------------------------------------------------

template <> struct IsSimple_<bool> { typedef True Type; };
template <> struct IsSimple_<char> { typedef True Type; };

template <> struct IsSimple_<unsigned char> { typedef True Type; };
template <> struct IsSimple_<unsigned short> { typedef True Type; };
template <> struct IsSimple_<unsigned int> { typedef True Type; };
template <> struct IsSimple_<unsigned long> { typedef True Type; };

template <> struct IsSimple_<signed char> { typedef True Type; };
template <> struct IsSimple_<signed short> { typedef True Type; };
template <> struct IsSimple_<signed int> { typedef True Type; };
template <> struct IsSimple_<signed long> { typedef True Type; };

template <> struct IsSimple_<float> { typedef True Type; };
template <> struct IsSimple_<double> { typedef True Type; };
template <> struct IsSimple_<long double> { typedef True Type; };

template <typename T>
struct IsSimple
{
    typedef typename IsSimple_<T>::Type Type;
    enum { VALUE = Type::VALUE };
};

// user defined types (re-specializations are allowed here)
template <> struct IsSimple<wchar_t> { typedef True Type; };
template <> struct IsSimple<__int64> { typedef True Type; };
template <> struct IsSimple<__uint64> { typedef True Type; };

template <typename T>
struct IsSimple<T const> : public IsSimple<T> {};

// ----------------------------------------------------------------------------
// Metafunction Value
// ----------------------------------------------------------------------------

// TODO(holtgrew): This should probably to into sequence module along with this header.

///.Metafunction.Value.param.T.type:Adaption.char array
///.Metafunction.Value.class:Adaption.char array

template <typename TValue>
struct Value<TValue *>
{
    typedef TValue Type;
};

template <typename TValue>
struct Value<TValue * const>
{
    typedef TValue Type;
};

// TODO(holtgrew): Is this still a problem with dropped 2003 support of VC++?

//The next two metafunctions dont work in VC++ due to a compiler bug.
//(the default implementation in common_type.h is called instead)
//work-around: convert arrays to pointers.
template <typename TValue, size_t SIZE>
struct Value<TValue [SIZE]>
{
    typedef TValue Type;
};

template <typename TValue, size_t SIZE>
struct Value<TValue const [SIZE]>
{
    typedef TValue Type;
};

// ----------------------------------------------------------------------------
// Metafunction Reference
// ----------------------------------------------------------------------------

// TODO(holtgrew): This should probably to into sequence module along with this header.

///.Metafunction.Reference.param.T.type:Adaption.char array
///.Metafunction.Reference.class:Adaption.char array

template <typename TValue>
struct Reference<TValue *>
{
    typedef TValue & Type;
};

template <typename TValue>
struct Reference<TValue * const>
{
    typedef TValue const & Type;
};

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function value() for pointers.
// ----------------------------------------------------------------------------

// TODO(holtgrew): This has to go to iterator module, adaption of pointers to iterators.

template <typename T>
inline T &
value(T * me)
{
    SEQAN_CHECKPOINT;
    return *me;
}

// ----------------------------------------------------------------------------
// Function getValue() for pointers.
// ----------------------------------------------------------------------------

// TODO(holtgrew): This has to go to iterator module, adaption of pointers to iterators.

template <typename T>
inline T &
getValue(T * me)
{
    SEQAN_CHECKPOINT;
    return value(me);
}

// TODO(holtgrew): All of the helper structs could be replaced by global functions.

// TODO(holtgrew): First, the generic versions for iterators are defined.  Below are the versions for pointers.

// ----------------------------------------------------------------------------
// Function valueConstruct() using iterators
// ----------------------------------------------------------------------------

/**
.Function.valueConstruct
..cat:Content Manipulation
..summary:Constructs an object at specified position.
..signature:valueConstruct(iterator [, param [, move_tag] ])
..param.iterator:Pointer or iterator to position where the object should be constructed.
..param.param:Parameter that is forwarded to constructor. (optional)
..param.move_tag:Instance of the @Tag.Move Switch.move switch tag@. (optional)
...remarks:If the @Tag.Move Switch.move switch tag@ is specified, it is forwarded to the constructor,
so the constructed object must support move construction.
..remarks:The type of the destructed object is the @Metafunction.Value.value type@ of $iterator$.
..include:seqan/basic.h
*/

// Helper code for constructing values behind iterators that do not return
// proxies from their value() functions but references.
struct ValueConstructor_ 
{
    template <typename TIterator>
    static inline void
    construct(TIterator it)
    {
        typedef typename Value<TIterator>::Type     TValue;
        typedef typename RemoveConst_<TValue>::Type TNonConstValue;
        new( (void*) & value(it) ) TNonConstValue;
    }

    template <typename TIterator, typename TParam>
    static inline void
    construct(TIterator it,
              TParam const & param_)
    {
        typedef typename Value<TIterator>::Type     TValue;
        typedef typename RemoveConst_<TValue>::Type TNonConstValue;
        new( (void*) & value(it) ) TNonConstValue(param_);
    }

    template <typename TIterator, typename TParam>
    static inline void
    construct(TIterator it,
              TParam & param_,
              Move const & tag)
    {
        typedef typename Value<TIterator>::Type     TValue;
        typedef typename RemoveConst_<TValue>::Type TNonConstValue;
        new( (void*) & value(it) ) TNonConstValue(param_, tag);
    }
};

// Helper code for constructing values behind iterators that return proxies
// from their value() function.
//
// TODO(holtgrew): These implementations are empty and to be overwritten. Should we have dynamic/static asserstions here?
struct ValueConstructorProxy_ 
{
    template <typename TIterator>
    static inline void construct(TIterator) {}

    template <typename TIterator, typename TParam>
    static inline void construct(TIterator, TParam const &) {}

    template <typename TIterator, typename TParam>
    static inline void construct(TIterator, TParam &, Move const & ) {}
};

template <typename TIterator>
inline void
valueConstruct(TIterator it)
{
    SEQAN_CHECKPOINT;
    typedef typename IfC<
        IsSameType<
            typename Value<TIterator>::Type &,
            typename Reference<TIterator>::Type
        >::VALUE,
    // THEN
        ValueConstructor_,          // true,  types are equal
    // ELSE
        ValueConstructorProxy_      // false, types differ -> value() returns a proxy
    >::Type TConstructor;

    TConstructor::construct(it);
}

template <typename TIterator, typename TParam>
inline void
valueConstruct(TIterator it,
               TParam const & param_)
{
    SEQAN_CHECKPOINT;
    typedef typename IfC<
        IsSameType<
            typename Value<TIterator>::Type &,
            typename Reference<TIterator>::Type
        >::VALUE,
    // THEN
        ValueConstructor_,          // true,  types are equal
    // ELSE
        ValueConstructorProxy_      // false, types differ -> value() returns a proxy
    >::Type TConstructor;

    TConstructor::construct(it, param_);
}

template <typename TIterator, typename TParam>
inline void
valueConstruct(TIterator it,
               TParam & param_,
               Move const & tag)
{
    SEQAN_CHECKPOINT;
    typedef typename IfC<
        IsSameType<
            typename Value<TIterator>::Type &,
            typename Reference<TIterator>::Type
        >::VALUE,
    // THEN
        ValueConstructor_,          // true,  types are equal
    // ELSE
        ValueConstructorProxy_      // false, types differ -> value() returns a proxy
    >::Type TConstructor;

    TConstructor::construct(it, param_, tag);
}

// ----------------------------------------------------------------------------
// Function valueDestruct() using iterators
// ----------------------------------------------------------------------------

// Helper code for destructing values behind iterators that do not return
// proxies from their value() function but references.
struct ValueDestructor_ 
{
    template <typename TValue>
    static inline void
    _destruct(TValue * p)
    {
        p->~TValue();
    }

    template <typename TIterator>
    static inline void
    destruct(TIterator it)
    {
        _destruct(&value(it));
    }
};

// Helper code for destructing values behind iterators that return proxies
// from their value() function.
//
// TODO(holtgrew): These implementations are empty and to be overwritten. Should we have dynamic/static asserstions here?
struct ValueDestructorProxy_ 
{
    template <typename TIterator>
    static inline void destruct(TIterator) {}
};

/**
.Function.valueDestruct
..cat:Content Manipulation
..summary:Destoys an object at specified position.
..signature:valueDestruct(iterator)
..param.iterator:Pointer or iterator to position where the object should be destructed.
..remarks:The type of the constructed object is the @Metafunction.Value.value type@ of $iterator$.
..see:Function.valueConstruct
..include:seqan/basic.h
*/

template <typename TIterator>
inline void
valueDestruct(TIterator it)
{
    SEQAN_CHECKPOINT;
    typedef typename IfC<
        IsSameType<
            typename Value<TIterator>::Type &,
            typename Reference<TIterator>::Type
        >::VALUE,
    // THEN
        ValueDestructor_,           // true,  types are equal
    // ELSE
        ValueDestructorProxy_       // false, types differ -> value() returns a proxy
    >::Type TDestructor;

    TDestructor::destruct(it);
}

// ----------------------------------------------------------------------------
// Function arrayConstruct() using iterators
// ----------------------------------------------------------------------------

/**
.Function.arrayConstruct
..cat:Array Handling
..summary:Construct objects in a given memory buffer.
..signature:arrayConstruct(begin, end [, value])
..param.begin:Iterator to the begin of the range that is to be constructed.
..param.end:Iterator behind the end of the range.
..param.value:Argument that is forwarded to the constructor. (optional)
...text:An appropriate constructor is required. 
If $value$ is not specified, the default constructor is used. 
..remarks:The type of the constructed Objects is the @Metafunction.Value.value type@
of $begin$ and $end$.
..see:Function.arrayDestruct
..see:Function.arrayConstructCopy
..see:Function.arrayFill
..see:Class.SimpleType
..see:Function.valueConstruct
..include:seqan/basic.h
*/

// NOTE(holtgrew): Of course, it does not make sense to declare this in a move version!

template<typename TIterator1, typename TIterator2>
inline void 
_arrayConstructDefault(TIterator1 begin_, 
                       TIterator2 end_)
{
    SEQAN_CHECKPOINT;
    while (begin_ != end_)
    {
        valueConstruct(begin_);
        ++begin_;
    }
}

template<typename TIterator1, typename TIterator2>
inline void 
arrayConstruct(TIterator1 begin_, 
               TIterator2 end_)
{
    SEQAN_CHECKPOINT;
    _arrayConstructDefault(begin_, end_);
}

template<typename TIterator1, typename TIterator2, typename TParam>
inline void 
_arrayConstructDefault(TIterator1 begin_, 
                       TIterator2 end_, 
                       TParam const & param_)
{
    SEQAN_CHECKPOINT;
    while (begin_ != end_)
    {
        valueConstruct(begin_, param_);
        ++begin_;
    }
}

template<typename TIterator1, typename TIterator2, typename TParam>
inline void 
arrayConstruct(TIterator1 begin_, 
               TIterator2 end_, 
               TParam const & param_)
{
    SEQAN_CHECKPOINT;
    _arrayConstructDefault(begin_, end_, param_);
}

// ----------------------------------------------------------------------------
// Function arrayConstructCopy() using iterators
// ----------------------------------------------------------------------------

/**
.Function.arrayConstructCopy
..cat:Array Handling
..summary:Copy constructs an array of objects into in a given memory buffer.
..signature:arrayConstructCopy(source_begin, source_end, target)
..param.source_begin:Iterator to the first element of the source range.
..param.source_end:Iterator behind the last element of the source range.
...text:$source_end$ should have the same type as $source_begin$.
..param.target:Pointer to the memory block the new objects will be constructed in.
...text:The type of $target$ specifies the type of the constructed objects:
If $T*$ is the type of $target$, then the function constructs objects of type $T$. 
...text:The memory buffer should be large enough to store $source_end$ - $source_begin$ objects.
An appropriate (copy-) constructor that constructs an target objects given a source object is required.
..see:Function.arrayDestruct
..see:Function.arrayCopyForward
..see:Function.arrayCopy
..see:Function.valueConstruct
..include:seqan/basic.h
*/

template<typename TTarget, typename TSource1, typename TSource2>
inline void 
_arrayConstructCopyDefault(TSource1 source_begin, 
                           TSource2 source_end, 
                           TTarget target_begin)
{
    SEQAN_CHECKPOINT;
    while (source_begin != source_end)
    {
        // NOTE(holtgrew): getValue() is used here since value() could return
        // a proxy!
        valueConstruct(target_begin, getValue(source_begin));
        ++source_begin;
        ++target_begin;
    }
}

template<typename TTarget, typename TSource1, typename TSource2>
inline void 
arrayConstructCopy(TSource1 source_begin, 
                   TSource2 source_end, 
                   TTarget target_begin)
{
    SEQAN_CHECKPOINT;
    _arrayConstructCopyDefault(source_begin, source_end, target_begin);
}

// ----------------------------------------------------------------------------
// Function arrayConstructMove() using iterators
// ----------------------------------------------------------------------------

/**
.Function.arrayConstructMove
..cat:Array Handling
..summary:Move constructs an array of objects into in a given memory buffer.
..signature:arrayConstructMove(source_begin, source_end, target)
..param.source_begin:Iterator to the first element of the source range.
..param.source_end:Iterator behind the last element of the source range.
...text:$source_end$ should have the same type as $source_begin$.
..param.target:Pointer to the memory block the new objects will be constructed in.
...text:The type of $target$ specifies the type of the constructed objects:
If $T*$ is the type of $target$, then the function constructs objects of type $T$. 
...text:The memory buffer should be large enough to store $source_end$ - $source_begin$ objects.
An appropriate move constructor that constructs an target objects given a source object is required.
..see:Function.arrayDestruct
..see:Function.arrayConstructCopy
..see:Function.arrayMoveForward
..see:Function.arrayMove
..see:Function.valueConstruct
..include:seqan/basic.h
*/

template<typename TTarget, typename TSource1, typename TSource2>
inline void 
_arrayConstructMoveDefault(TSource1 source_begin, 
                           TSource2 source_end, 
                           TTarget target_begin)
{
    SEQAN_CHECKPOINT;
    while (source_begin < source_end)
    {
        // NOTE(holtgrew): Using value() here, used to be getValue() but
        // cannot move from const reference or proxy.
        // valueConstruct(target_begin, value(source_begin), Move());
        // TODO(holtgrew): We need a "has move constructor" metafunction to switch between move/copy constructing before we can use the line here.
        valueConstruct(target_begin, value(source_begin));
        ++source_begin;
        ++target_begin;
    }
}

template<typename TTarget, typename TSource1, typename TSource2>
inline void 
arrayConstructMove(TSource1 source_begin, 
                   TSource2 source_end, 
                   TTarget target_begin)
{
    SEQAN_CHECKPOINT;
    _arrayConstructMoveDefault(source_begin, source_end, target_begin);
}

// ----------------------------------------------------------------------------
// Function arrayDestruct() using iterators
// ----------------------------------------------------------------------------

/**
.Function.arrayDestruct
..cat:Array Handling
..summary:Destroys an array of objects.
..signature:arrayDestruct(begin, end)
..param.begin:Iterator to the begin of the range that is to be destructed.
..param.end:Iterator behind the end of the range.
..remarks:This function does not deallocates the memory.
..see:Class.SimpleType
..see:Function.valueDestruct
..include:seqan/basic.h
*/

template<typename TIterator1, typename TIterator2>
inline void 
_arrayDestructDefault(TIterator1 begin_, 
                      TIterator2 end_)
{
    SEQAN_CHECKPOINT;
    while (begin_ != end_)
    {
        valueDestruct(begin_);
        ++begin_;
    }
}

template<typename TIterator1, typename TIterator2>
inline void 
arrayDestruct(TIterator1 begin_, 
              TIterator2 end_)
{
    SEQAN_CHECKPOINT;
    _arrayDestructDefault(begin_, end_);
}

// ----------------------------------------------------------------------------
// Function arrayFill() using iterators
// ----------------------------------------------------------------------------

// TODO(holtgrew): What is the advantage over arrayConstruct() with prototype?

/**
.Function.arrayFill
..cat:Array Handling
..summary:Assigns one object to each element of a range.
..signature:arrayFill(begin, end, value)
..param.begin:Iterator to the begin of the range that is to be filled.
..param.end:Iterator behind the end of the range.
..param.value:Argument that is assigned to all $count$ objects in $array$.
..remarks:All objects $target_begin[0]$ to $target_begin[count-1]$ are set to $value$.
..see:Function.arrayCopy
..see:Function.arrayCopyForward
..include:seqan/basic.h
*/

// TODO(holtgrew): Redirects to fill_n. What are the exact semantics here? Do the array elements have to be initialized already? fill_n uses assignment, not copy construction!

template <typename TIterator, typename TValue>
inline void 
arrayFill(TIterator begin_,
          TIterator end_,
          TValue const & value)
{
    SEQAN_CHECKPOINT;
    ::std::fill_n(begin_, end_ - begin_, value);
}

template <typename TIterator, typename TValue>
inline void 
arrayFill(TIterator begin_,
          TIterator end_, 
          TValue const & value,
          Serial)
{
    arrayFill(begin_, end_, value);
}

// ----------------------------------------------------------------------------
// Function arrayCopyForward() using iterators
// ----------------------------------------------------------------------------

/**
.Function.arrayCopyForward
..cat:Array Handling
..summary:Copies a range of objects into another range of objects starting from the first element.
..signature:arrayCopyForward(source_begin, source_end, target)
..param.source_begin:Iterator to the first element of the source array.
..param.source_end:Iterator behind the last element of the source array.
...text:$source_end$ must have the same type as $source_begin$.
..param.target:Iterator to the first element of the target array.
...text:The target capacity should be at least as long as the source range.
..remarks.note:Be careful if source and target range overlap, because in this case some source elements could be accidently overwritten before they are copied.
..remarks:If there is no need for the source elements to persist, consider to use 
@Function.arrayMoveForward@ instead to improve performance.
..see:Class.SimpleType
..include:seqan/basic.h
*/

template<typename TTarget, typename TSource1, typename TSource2>
inline void 
_arrayCopyForwardDefault(TSource1 source_begin, 
                         TSource2 source_end, 
                         TTarget target_begin)
{
    SEQAN_CHECKPOINT;
    ::std::copy(source_begin, source_end, target_begin);
}

template<typename TTarget, typename TSource1, typename TSource2>
inline void 
arrayCopyForward(TSource1 source_begin, 
                 TSource2 source_end, 
                 TTarget target_begin)
{
    SEQAN_CHECKPOINT;
    _arrayCopyForwardDefault(source_begin, source_end, target_begin);   
}

// ----------------------------------------------------------------------------
// Function arrayCopyBackward() using iterators
// ----------------------------------------------------------------------------

/**
.Function.arrayCopyBackward
..cat:Array Handling
..summary:Copies a range of objects into another range of objects starting from the last element.
..signature:arrayCopyBackward(source_begin, source_end, target)
..param.source_begin:Iterator to the first element of the source array.
..param.source_end:Iterator behind the last element of the source array.
...text:$source_end$ must have the same type as $source_begin$.
..param.target:Iterator to the first element of the target array.
...text:The target capacity should be at least as long as the source range.
..remarks.note:Be careful if source and target range overlap, because in this case
    some source elements could be accidently overwritten before they are moved.
..remarks.text:If source and target do not overlap, consider to use the function
@Function.arrayCopyForward@ instead that is faster in some cases.
..remarks:If there is no need for the source elements to persist, consider to use 
@Function.arrayMoveBackward@ instead to improve performance.
..remarks.note:The semantic of this function's argument $target$ differ from the arguments of $::std::copy_backward$.
..see:Function.arrayCopyForward
..see:Class.SimpleType
..include:seqan/basic.h
*/

template<typename TTarget, typename TSource1, typename TSource2>
inline void 
_arrayCopyBackwardDefault(TSource1 source_begin, 
                           TSource2 source_end, 
                           TTarget target_begin)
{
    SEQAN_CHECKPOINT;
    ::std::copy_backward(source_begin, source_end, target_begin + (source_end - source_begin));
}

template<typename TTarget, typename TSource1, typename TSource2>
inline void 
arrayCopyBackward(TSource1 source_begin, 
                  TSource2 source_end, 
                  TTarget target_begin)
{
    SEQAN_CHECKPOINT;
    _arrayCopyBackwardDefault(source_begin, source_end, target_begin);
}

// ----------------------------------------------------------------------------
// Function arrayCopy() using iterators
// ----------------------------------------------------------------------------

/**
.Function.arrayCopy
..cat:Array Handling
..summary:Copies a range of objects into another range of objects.
..signature:arrayCopy(source_begin, source_end, target)
..param.source_begin:Iterator to the first element of the source range.
..param.source_end:Iterator behind the last element of the source range.
...text:$source_end$ must have the same type as $source_begin$.
..param.target:Iterator to the first element of the target range.
...text:The target capacity should be at least as long as the source range.
..remarks.text:If source and target range do not overlap, consider to use
    @Function.arrayCopyForward@ instead to improve performance.
..remarks:If there is no need for the source elements to persist, consider to use 
    @Function.arrayMoveForward@ instead to improve performance.
..DISABLED.remarks.note:Be careful if source and target range overlap and the size of the
    source elements differ from the size of target elements, because in this case
    some source elements could be accidently overwritten before they are moved.
..see:Function.arrayCopyForward
..see:Function.arrayCopyBackward
..see:Class.SimpleType
..include:seqan/basic.h
*/

template<typename TTarget, typename TSource1, typename TSource2>
inline void arrayCopy(TSource1 source_begin, 
                      TSource2 source_end, 
                      TTarget target_begin)
{
    if ((void *) source_begin >= (void *) target_begin) {
        SEQAN_CHECKPOINT;
        arrayCopyForward(source_begin, source_end, target_begin);
    } else {
        SEQAN_CHECKPOINT;
        arrayCopyBackward(source_begin, source_end, target_begin);
    }
}

// ----------------------------------------------------------------------------
// Function arrayMoveForward() using iterators
// ----------------------------------------------------------------------------

/**
.Function.arrayMoveForward
..cat:Array Handling
..summary:Moves a range of objects into another range of objects starting from the first element.
..signature:arrayMoveForward(source_begin, source_end, target)
..param.source_begin:Iterator to the first element of the source array.
..param.source_end:Iterator behind the last element of the source array.
...text:$source_end$ must have the same type as $source_begin$.
..param.target:Iterator to the first element of the target array.
...text:The target capacity should be at least as long as the source range.
..remarks:The function possibly clears (but does not destroy) the source elements.
    If source elements must persist, consider to use @Function.arrayCopyForward@ instead.
..remarks.note:Be careful if source and target range overlap, because in this case
    some source elements could be accidently overwritten before they are moved.
..see:Function.arrayCopyForward
..see:Class.SimpleType
..include:seqan/basic.h
*/

template<typename TTarget, typename TSource1, typename TSource2>
inline void 
_arrayMoveForwardDefault(TSource1 source_begin, 
                          TSource2 source_end, 
                          TTarget target_begin)
{
    SEQAN_CHECKPOINT;
    while (source_begin != source_end)
    {
        move(*target_begin, *source_begin);
        ++source_begin;
        ++target_begin;
    }
}

template<typename TTarget, typename TSource1, typename TSource2>
inline void 
arrayMoveForward(TSource1 source_begin, 
                 TSource2 source_end, 
                 TTarget target_begin)
{
    SEQAN_CHECKPOINT;
    _arrayMoveForwardDefault(source_begin, source_end, target_begin);   
}

// ----------------------------------------------------------------------------
// Function arrayMoveBackward() using iterators
// ----------------------------------------------------------------------------

/**
.Function.arrayMoveBackward
..cat:Array Handling
..summary:Moves a range of objects into another range of objects starting from the last element.
..signature:arrayMoveBackward(source_begin, source_end, target)
..param.source_begin:Iterator to the first element of the source array.
..param.source_end:Iterator behind the last element of the source array.
...text:$source_end$ must have the same type as $source_begin$.
..param.target:Iterator to the first element of the target array.
...text:The target capacity should be at least as long as the source range.
..remarks:The function possibly clears (but does not destroy) the source elements.
    If source elements must persist, consider to use @Function.arrayCopyBackward@ instead.
..remarks.note:Be careful if source and target range overlap, because in this case
    some source elements could be accidently overwritten before they are moved.
..remarks.text:If source and target do not overlap, consider to use the function
@Function.arrayMoveForward@ instead that is faster in some cases.
..remarks.note:The semantic of this function's argument $target$ differ from the arguments of $::std::copy_backward$.
..see:Function.arrayMoveForward
..see:Function.arrayCopyBackward
..see:Class.SimpleType
..include:seqan/basic.h
*/

template<typename TTarget, typename TSource1, typename TSource2>
inline void 
_arrayMoveBackwardDefault(TSource1 source_begin, 
                          TSource2 source_end, 
                          TTarget target_begin)
{
    SEQAN_CHECKPOINT;
    target_begin += (source_end - source_begin);
    while (source_end != source_begin) {
        --source_end;
        --target_begin;
        move(*target_begin, *source_end);
    }
}

template<typename TTarget, typename TSource1, typename TSource2>
inline void 
arrayMoveBackward(TSource1 source_begin, 
                  TSource2 source_end, 
                  TTarget target_begin)
{
    SEQAN_CHECKPOINT;
    _arrayMoveBackwardDefault(source_begin, source_end, target_begin);
}

// ----------------------------------------------------------------------------
// Function arrayMove using iterators
// ----------------------------------------------------------------------------

/**
.Function.arrayMove
..cat:Array Handling
..summary:Moves a range of objects into another range of objects.
..signature:arrayMove(source_begin, source_end, target)
..param.source_begin:Iterator to the first element of the source range.
..param.source_end:Iterator behind the last element of the source range.
...text:$source_end$ must have the same type as $source_begin$.
..param.target:Iterator to the first element of the target range.
...text:The target capacity should be at least as long as the source range.
..remarks:The function possibly clears (but does not destroy) the source elements.
    If source elements must persist, consider to use @Function.arrayCopy@ instead.
..remarks.text:If source and target range do not overlap, consider to use
    @Function.arrayMoveForward@ instead to improve performance.
..DISABLED.remarks.note:Be careful if source and target range overlap and the size of the
    source elements differ from the size of target elements, because in this case
    some source elements could be accidently overwritten before they are moved.
..remarks.note:Don't confuse this function with the standard $move$ function that
resembles @Function.arrayCopy@.
..see:Function.arrayMoveForward
..see:Function.arrayMoveBackward
..see:Function.arrayCopy
..see:Class.SimpleType
..include:seqan/basic.h
*/

template<typename TTarget, typename TSource1, typename TSource2>
inline void 
arrayMove(TSource1 source_begin, 
          TSource2 source_end,
          TTarget target_begin)
{
    if ((void *) source_begin >= (void *) target_begin) {
        SEQAN_CHECKPOINT;
        arrayMoveForward(source_begin, source_end, target_begin);
    } else {
        SEQAN_CHECKPOINT;
        arrayMoveBackward(source_begin, source_end, target_begin);
    }
}

// ----------------------------------------------------------------------------
// Function arrayClearSpace() using iterators
// ----------------------------------------------------------------------------

/**
.Function.arrayClearSpace
..cat:Array Handling
..summary:Destroys the begin of an array and keeps the rest.
..signature:arrayClearSpace(arr_begin, arr_length, keep_from, move_to)
..param.arr_begin:Pointer to the first element of the array.
..param.arr_length:Length of the array.
..param.keep_from:Offset of the first object that will be kept.
..param.move_to:Offset the first kept object will get at the end of the function. 
..remarks.text:The objects $arr[keep_from]$ to $arr[arr_length-1]$
are moved to the area beginning at positions $move_to$. 
All objects in $arr[0]$ to $arr[keep_from-1]$ are destroyed.
After this function, the first $move_to$ positions of the array
are free and dont contain objects. 
..remarks.text:The array must have at least enough space to store $arr_length + move_to - keep_from$ objects.
..remarks.text:The objects from $arr[0]$ to $arr[array_length-1]$ have to be initialized/constructed, arrays beyond $arr[array_length-1]$ are assumed not to be constructed. If this assumption is violated, memory might leak.
..see:Function.arrayCopy
..see:Function.arrayDestruct
..see:Function.arrayCopyForward
..see:Class.SimpleType
..include:seqan/basic.h
*/

// TODO(holtgrew): The feature that the range [0, array_begin) is deleted is used nowhere. Can this be removed to simplify behaviour?

template <typename TIterator>
void _arrayClearSpaceDefault(TIterator array_begin, 
                             size_t array_length, 
                             size_t keep_from, 
                             size_t move_to)
{
    if (keep_from == array_length) {
        // In the simplest case, we only destruct elements.
        arrayDestruct(array_begin, array_begin + array_length);
        return;
    }

    // Otherwise, we will perform the destruction & movement.
    SEQAN_ASSERT_LT(keep_from, array_length);
    if (keep_from == move_to) {
        // Case 1: No movement, just destroy elements.
        arrayDestruct(array_begin, array_begin + move_to);
        return;
    } else if (keep_from < move_to) {
        // Case 2: Move to the right.
        if (array_length > move_to) {
            // Case 2a: Moving right of array_length, i.e. we can move a part
            // of the objects and have to move-construct the rest.
            SEQAN_CHECKPOINT;
            size_t middle = array_length - (move_to - keep_from);
            arrayConstructMove(array_begin + middle, array_begin + array_length, array_begin + array_length);
            arrayMove(array_begin + keep_from, array_begin + middle, array_begin + move_to);
            arrayDestruct(array_begin, array_begin + move_to);
        } else {
            // Case 2b: We have to move-construct all target objects.
            SEQAN_CHECKPOINT;
            arrayConstructMove(array_begin + keep_from, array_begin + array_length, array_begin + move_to);
            arrayDestruct(array_begin, array_begin + array_length);
        }
    } else {
        // Case 3: Move to the left.
        SEQAN_CHECKPOINT;
        arrayMove(array_begin + keep_from, array_begin + array_length, array_begin + move_to);
        arrayDestruct(array_begin, array_begin + move_to);
        arrayDestruct(array_begin + array_length - (keep_from - move_to), array_begin + array_length);
    }
}

template <typename TIterator>
void arrayClearSpace(TIterator array_begin, 
                     size_t array_length, 
                     size_t keep_from, 
                     size_t move_to)
{
    _arrayClearSpaceDefault(array_begin, array_length, keep_from, move_to);
}

// ----------------------------------------------------------------------------
// Function arrayConstruct() using pointers
// ----------------------------------------------------------------------------

template<typename TIterator>
inline void 
_arrayConstructPointer(TIterator, 
                       TIterator,
                       True)
{
    SEQAN_CHECKPOINT;
    //nothing to do
}

template<typename TIterator>
inline void 
_arrayConstructPointer(TIterator begin_, 
                       TIterator end_,
                       False)
{
    SEQAN_CHECKPOINT;
    _arrayConstructDefault(begin_, end_);
}

template<typename TValue>
inline void 
arrayConstruct(TValue * begin_, 
               TValue * end_)
{
    SEQAN_CHECKPOINT;
    _arrayConstructPointer(begin_, end_, typename IsSimple<TValue>::Type() );
}

template<typename TIterator, typename TParam>
inline void 
_arrayConstructPointer(TIterator begin_, 
                        TIterator end_, 
                        TParam const & param_,
                        True)
{
    SEQAN_CHECKPOINT;
    arrayFill(begin_, end_, param_);
}

template<typename TIterator, typename TParam>
inline void 
_arrayConstructPointer(TIterator begin_, 
                        TIterator end_, 
                        TParam const & param_,
                        False)
{
    SEQAN_CHECKPOINT;
    _arrayConstructDefault(begin_, end_, param_);
}

template<typename TValue, typename TParam>
inline void 
arrayConstruct(TValue * begin_, 
               TValue * end_, 
               TParam const & param_)
{
    SEQAN_CHECKPOINT;
    _arrayConstructPointer(begin_, end_, param_, typename IsSimple<TValue>::Type());
}

// ----------------------------------------------------------------------------
// Function arrayConstructCopy() using pointers
// ----------------------------------------------------------------------------

template<typename TValueSource, typename TValueTarget>
inline void 
_arrayConstructCopyPointer(TValueSource * source_begin, 
                            TValueSource * source_end, 
                            TValueTarget * target_begin,
                            True)
{
    SEQAN_CHECKPOINT;
    arrayCopyForward(source_begin, source_end, target_begin);
}

template<typename TValueSource, typename TValueTarget>
inline void 
_arrayConstructCopyPointer(TValueSource * source_begin, 
                            TValueSource * source_end, 
                            TValueTarget const* target_begin,
                            True)
{
    SEQAN_CHECKPOINT;
    arrayCopyForward(source_begin, source_end, const_cast<TValueTarget *>(target_begin));
}

template<typename TValueSource, typename TValueTarget>
inline void 
_arrayConstructCopyPointer(TValueSource * source_begin, 
                            TValueSource * source_end, 
                            TValueTarget * target_begin,
                            False)
{
    SEQAN_CHECKPOINT;
    _arrayConstructCopyDefault(source_begin, source_end, target_begin);
}
template<typename TValueSource, typename TValueTarget>
inline void 
arrayConstructCopy(TValueSource * source_begin, 
                   TValueSource * source_end, 
                   TValueTarget * target_begin)
{
    SEQAN_CHECKPOINT;
    _arrayConstructCopyPointer(source_begin, source_end, target_begin, typename IsSimple<TValueTarget>::Type() );
}

// ----------------------------------------------------------------------------
// Function arrayConstructMove() using pointers
// ----------------------------------------------------------------------------

template<typename TValue>
inline void 
_arrayConstructMovePointer(TValue * source_begin, 
                            TValue * source_end, 
                            TValue * target_begin,
                            True)
{
    SEQAN_CHECKPOINT;
    arrayMoveForward(source_begin, source_end, target_begin);
}

template<typename TValue>
inline void 
_arrayConstructMovePointer(TValue * source_begin, 
                            TValue * source_end, 
                            TValue * target_begin,
                            False)
{
    SEQAN_CHECKPOINT;
    _arrayConstructMoveDefault(source_begin, source_end, target_begin);
}

template<typename TValue>
inline void 
arrayConstructMove(TValue * source_begin, 
                   TValue * source_end, 
                   TValue * target_begin)
{
    SEQAN_CHECKPOINT;
    _arrayConstructMovePointer(source_begin, source_end, target_begin, typename IsSimple<TValue>::Type() );
}

// ----------------------------------------------------------------------------
// Function arrayDestruct() using pointers
// ----------------------------------------------------------------------------

template<typename TValue>
inline void 
_arrayDestructPointer(TValue * /*begin_*/, 
                       TValue * /*end_*/,
                       True)
{
    SEQAN_CHECKPOINT;
    //do nothing
}

template<typename TValue>
inline void 
_arrayDestructPointer(TValue * begin_, 
                       TValue * end_,
                       False)
{
    SEQAN_CHECKPOINT;
    _arrayDestructDefault(begin_, end_);
}

template<typename TValue>
inline void 
arrayDestruct(TValue * begin_, 
              TValue * end_)
{
    SEQAN_CHECKPOINT;
    _arrayDestructPointer(begin_, end_, typename IsSimple<TValue>::Type() );
}

// ----------------------------------------------------------------------------
// Function arrayFill() using pointers
// ----------------------------------------------------------------------------

// TODO(holtgrew): Missing?

//no specializiation for pointer to simple

// ----------------------------------------------------------------------------
// Function arrayCopyBackward() using pointers
// ----------------------------------------------------------------------------

template<typename TValue>
inline void 
_arrayCopyForwardPointer(TValue * source_begin, 
                          TValue * source_end, 
                          TValue * target_begin,
                          True)
{
    SEQAN_CHECKPOINT;
    ::std::memmove(target_begin, source_begin, (source_end - source_begin) * sizeof(TValue));
}

template<typename TValue>
inline void 
_arrayCopyForwardPointer(TValue * source_begin, 
                          TValue * source_end, 
                          TValue * target_begin,
                          False)
{
    SEQAN_CHECKPOINT;
    _arrayCopyForwardDefault(source_begin, source_end, target_begin);
}

template<typename TValue>
inline void 
arrayCopyForward(TValue * source_begin, 
                 TValue * source_end, 
                 TValue * target_begin)
{
    SEQAN_CHECKPOINT;
    _arrayCopyForwardPointer(source_begin, source_end, target_begin, typename IsSimple<TValue>::Type() );
}

template <typename TValue>
inline void 
_arrayCopyBackwardPointer(TValue * source_begin, 
                           TValue * source_end, 
                           TValue * target_begin,
                           True)
{
    SEQAN_CHECKPOINT;
    ::std::memmove(target_begin, source_begin, (source_end - source_begin) * sizeof(TValue));
}

template <typename TValue>
inline void 
_arrayCopyBackwardPointer(TValue * source_begin, 
                           TValue * source_end, 
                           TValue * target_begin,
                           False)
{
    SEQAN_CHECKPOINT;
    _arrayCopyBackwardDefault(source_begin, source_end, target_begin); 
}

template<typename TValue>
inline void 
arrayCopyBackward(TValue * source_begin, 
                  TValue * source_end, 
                  TValue * target_begin)
{
    SEQAN_CHECKPOINT;
    _arrayCopyBackwardPointer(source_begin, source_end, target_begin, typename IsSimple<TValue>::Type() );
}

// ----------------------------------------------------------------------------
// Function arrayMoveBackward() using pointers
// ----------------------------------------------------------------------------

template<typename TValue>
inline void 
_arrayMoveForwardPointer(TValue * source_begin, 
                          TValue * source_end, 
                          TValue * target_begin,
                          True)
{
    SEQAN_CHECKPOINT;
    ::std::memmove(target_begin, source_begin, (source_end - source_begin) * sizeof(TValue));
}

template<typename TValue>
inline void 
_arrayMoveForwardPointer(TValue * source_begin, 
                          TValue * source_end, 
                          TValue * target_begin,
                          False)
{
    SEQAN_CHECKPOINT;
    _arrayMoveForwardDefault(source_begin, source_end, target_begin);
}

template<typename TValue>
inline void 
arrayMoveForward(TValue * source_begin, 
                 TValue * source_end, 
                 TValue * target_begin)
{
    SEQAN_CHECKPOINT;
    _arrayMoveForwardPointer(source_begin, source_end, target_begin, typename IsSimple<TValue>::Type() );
}

template <typename TValue>
inline void 
_arrayMoveBackwardPointer(TValue * source_begin, 
                           TValue * source_end, 
                           TValue * target_begin,
                           True)
{
    SEQAN_CHECKPOINT;
    ::std::memmove(target_begin, source_begin, (source_end - source_begin) * sizeof(TValue));
}
template <typename TValue>
inline void 
_arrayMoveBackwardPointer(TValue * source_begin, 
                           TValue * source_end, 
                           TValue * target_begin,
                           False)
{
    SEQAN_CHECKPOINT;
    _arrayMoveBackwardDefault(source_begin, source_end, target_begin); 
}

template<typename TValue>
inline void 
arrayMoveBackward(TValue * source_begin, 
                  TValue * source_end, 
                  TValue * target_begin)
{
    SEQAN_CHECKPOINT;
    _arrayMoveBackwardPointer(source_begin, source_end, target_begin, typename IsSimple<TValue>::Type() );
}

// ----------------------------------------------------------------------------
// Function arrayClearSpace() using pointers
// ----------------------------------------------------------------------------

// clearSpace() on simple type using pointers.
template <typename TValue>
inline void 
_arrayClearSpacePointer(TValue * array_begin, 
                        size_t array_length, 
                        size_t keep_from, 
                        size_t move_to,
                        True const & /*isSimple*/)
{
    if (keep_from == move_to) return;
    SEQAN_CHECKPOINT;
    // TODO(holtgrew): arrayCopy is more appropriate here since we are dealing with the IsSimple case.
    arrayMove(array_begin + keep_from, array_begin + array_length, array_begin + move_to);
}

// clearSpace() on non-simple type using pointers.
template <typename TValue>
inline void 
_arrayClearSpacePointer(TValue * array_begin, 
                        size_t array_length, 
                        size_t keep_from, 
                        size_t move_to,
                        False const & /*isSimple*/)
{
    _arrayClearSpaceDefault(array_begin, array_length, keep_from, move_to);
}

template <typename TValue>
void arrayClearSpace(TValue * array_begin, 
                     size_t array_length, 
                     size_t keep_from, 
                     size_t move_to)
{
    _arrayClearSpacePointer(array_begin, array_length, keep_from, move_to, typename IsSimple<TValue>::Type());
}


}  // namespace seqan

#endif  // #ifndef SEQAN_CORE_INCLUDE_SEQAN_BASIC_ARRAY_CONSTRUCT_DESTRUCT_H_
