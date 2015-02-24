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
// Author: Andreas Gogol-Döring <andreas.doering@mdc-berlin.de>
// ==========================================================================
// Functions for destructing and several ways of constructing (e.g. copy,
// move) values or arrays of values.
// ==========================================================================

// TODO(holtgrew): Order of parameters should be (target1, target2, ..., source1, source2, ...).
// TODO(holtgrew): Can we maybe replace at least part with http://www.cplusplus.com/reference/std/memory/?
// TODO(holtgrew): The metafunction should go into the alphabet submodule, the functions into the sequence/string module.

#include <new>

#ifndef SEQAN_INCLUDE_SEQAN_BASIC_ARRAY_CONSTRUCT_DESTRUCT_H_
#define SEQAN_INCLUDE_SEQAN_BASIC_ARRAY_CONSTRUCT_DESTRUCT_H_

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

/*!
 * @mfn IsSimple
 * @headerfile <seqan/basic.h>
 * @brief Tests type to be simple.
 *
 * @signature IsSimple<T>::Type;
 *
 * @tparam T Type that is tested.
 *
 * @return Type Either True or False, depending on T being a "POD" type.
 *
 * A simple type is a type that does not need constructors to be created, a destructor to be destroyed, and copy
 * assignment operators or copy constructors to be copied.  All POD ("plain old data") types are simple, but some
 * non-POD types could be simple too, e.g. some specializations of SimpleType.
 *
 * @see SimpleType
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
SEQAN_HOST_DEVICE inline T &
value(T * me)
{
    return *me;
}

// ----------------------------------------------------------------------------
// Function getValue() for pointers.
// ----------------------------------------------------------------------------

// TODO(holtgrew): This has to go to iterator module, adaption of pointers to iterators.

template <typename T>
SEQAN_HOST_DEVICE inline T &
getValue(T * me)
{
    return value(me);
}

// TODO(holtgrew): All of the helper structs could be replaced by global functions.

// TODO(holtgrew): First, the generic versions for iterators are defined.  Below are the versions for pointers.

// ----------------------------------------------------------------------------
// Function valueConstruct() using iterators
// ----------------------------------------------------------------------------

/*!
 * @fn valueConstruct
 * @headerfile <seqan/basic.h>
 * @brief Constructs an object at specified position.
 *
 * @signature void valueConstruct(iterator [, param [, move_tag]]);
 *
 * @param[in,out] iterator Pointer or iterator to position where the object should be constructed.
 * @param[in]     param    Parameter that is forwarded to constructor.
 * @param[in]     moveTag  Instance of the move tag.  If the tag is specified, it is forwarded to the constructor,
 *                         so the constructed object must support move construction.
 *
 * The type of the destructed object is the value type of <tt>iterator</tt>.
 */

// Helper code for constructing values behind iterators that do not return
// proxies from their value() functions but references.
struct ValueConstructor_
{
    template <typename TIterator>
    static inline void
    construct(TIterator it)
    {
        typedef typename Value<TIterator>::Type    TValue;
        typedef typename RemoveConst<TValue>::Type TNonConstValue;
        new( (void*) & value(it) ) TNonConstValue;
    }

    template <typename TIterator, typename TParam>
    static inline void
    construct(TIterator it,
              TParam SEQAN_FORWARD_CARG param_)
    {
        typedef typename Value<TIterator>::Type    TValue;
        typedef typename RemoveConst<TValue>::Type TNonConstValue;
        new( (void*) & value(it) ) TNonConstValue(SEQAN_FORWARD(TParam, param_));
    }

#ifndef SEQAN_CXX11_STANDARD
    template <typename TIterator, typename TParam>
    static inline void
    construct(TIterator it,
              TParam & param_,
              Move const & tag,
              True)
    {
        typedef typename Value<TIterator>::Type    TValue;
        typedef typename RemoveConst<TValue>::Type TNonConstValue;
        new( (void*) & value(it) ) TNonConstValue(param_, tag);
    }

    template <typename TIterator, typename TParam>
    static inline void
    construct(TIterator it,
              TParam & param_,
              Move const &,
              False)
    {
        typedef typename Value<TIterator>::Type    TValue;
        typedef typename RemoveConst<TValue>::Type TNonConstValue;
        new( (void*) & value(it) ) TNonConstValue(param_);
    }

    template <typename TIterator, typename TParam>
    static inline void
    construct(TIterator it,
              TParam & param_,
              Move const & tag)
    {
        typedef typename Value<TIterator>::Type    TValue;
        typedef typename RemoveConst<TValue>::Type TNonConstValue;
        construct(it, param_, tag, typename HasMoveConstructor<TNonConstValue>::Type());
    }
#endif
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
    static inline void construct(TIterator, TParam SEQAN_FORWARD_CARG) {}

#ifndef SEQAN_CXX11_STANDARD
    template <typename TIterator, typename TParam>
    static inline void construct(TIterator, TParam &, Move const & ) {}
#endif
};

template <typename TIterator>
inline void
valueConstruct(TIterator it)
{
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
               TParam SEQAN_FORWARD_CARG param_)
{
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

    TConstructor::construct(it, SEQAN_FORWARD(TParam, param_));
}

#ifndef SEQAN_CXX11_STANDARD
template <typename TIterator, typename TParam>
inline void
valueConstruct(TIterator it,
               TParam & param_,
               Move const & tag)
{
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
#endif

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

/*!
 * @fn valueDestruct
 * @headerfile <seqan/basic.h>
 * @brief Destroys an object at specified position.
 *
 * @signature void valueDestruct(iterator);
 *
 * @param[in,out] iterator Pointer or iterator to position where the object should be destructed.
 *
 * The type of the constructed object is the value type of <tt>iterator</tt>.
 */

template <typename TIterator>
inline void
valueDestruct(TIterator it)
{
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

/*!
 * @fn arrayConstruct
 * @headerfile <seqan/basic.h>
 * @brief Construct objects in a given memory buffer.
 *
 * @signature void arrayConstruct(begin, end[, value]);
 *
 * @param[in] begin Iterator to the begin of the range that is to be constructed.
 * @param[in] end   Iterator behind the end of the range.
 * @param[in] value Argument that is forwarded to the constructor.  An appropriate constructor is required.  If
 *                 <tt>value</tt> is not specified, the default constructor is used.
 *
 * The type of the constructed Objects is the value type of <tt>begin</tt> and <tt>end</tt>.
 */

// NOTE(holtgrew): Of course, it does not make sense to declare this in a move version!

template<typename TIterator1, typename TIterator2>
inline void
_arrayConstructDefault(TIterator1 begin_,
                       TIterator2 end_)
{
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
    _arrayConstructDefault(begin_, end_);
}

template<typename TIterator1, typename TIterator2, typename TParam>
inline void
_arrayConstructDefault(TIterator1 begin_,
                       TIterator2 end_,
                       TParam const & param_)
{
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
    _arrayConstructDefault(begin_, end_, param_);
}

// ----------------------------------------------------------------------------
// Function arrayConstructCopy() using iterators
// ----------------------------------------------------------------------------

/*!
 * @fn arrayConstructCopy
 * @headerfile <seqan/basic.h>
 * @brief Copy constructs an array of objects into in a given memory buffer.
 *
 * @signature void arrayConstructCopy(sourceBegin, sourceEnd, target);
 *
 * @param[in] sourceBegin Iterator to the first element of the source range.
 * @param[in] sourceEnd   Iterator behind the last element of the source range.  <tt>sourceEnd</tt> should have the same
 *                        type as <tt>sourceBegin</tt>.
 * @param[in] target      Pointer to the memory block the new objects will be constructed in.  The type of <tt>target</tt>
 *                        specifies the type of the constructed objects: If <tt>T*</tt> is the type of <tt>target</tt>, then
 *                        the function constructs objects of type <tt>T</tt>.  The memory buffer should be large enough to
 *                        store <tt>sourceEnd</tt> - <tt>sourceBegin</tt> objects.  An appropriate (copy-) constructor that
 *                        constructs an target objects given a source object is required.
 */

template<typename TTarget, typename TSource1, typename TSource2>
inline void
_arrayConstructCopyDefault(TSource1 source_begin,
                           TSource2 source_end,
                           TTarget target_begin)
{
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
    _arrayConstructCopyDefault(source_begin, source_end, target_begin);
}

// ----------------------------------------------------------------------------
// Function arrayConstructMove() using iterators
// ----------------------------------------------------------------------------

/*!
 * @fn arrayConstructMove
 * @headerfile <seqan/basic.h>
 * @brief Move constructs an array of objects into in a given memory buffer.
 *
 * @signature void arrayConstructMove(sourceBegin, sourceEnd, target);
 *
 * @param[in] sourceEnd   Iterator behind the last element of the source range.  <tt>sourceEnd</tt> should have the same
 *                        type as <tt>sourceBegin</tt>.
 * @param[in] sourceBegin Iterator to the first element of the source range.
 * @param[in] target      Pointer to the memory block the new objects will be constructed in.  The type of <tt>target</tt>
 *                        specifies the type of the constructed objects: If <tt>T*</tt> is the type of <tt>target</tt>, then
 *                        the function constructs objects of type <tt>T</tt>.  The memory buffer should be large enough to
 *                        store <tt>sourceEnd</tt> - <tt>sourceBegin</tt> objects.  An appropriate move constructor that
 *                        constructs an target objects given a source object is required.
 */

template<typename TTarget, typename TSource1, typename TSource2>
inline void
_arrayConstructMoveDefault(TSource1 source_begin,
                           TSource2 source_end,
                           TTarget target_begin)
{
    while (source_begin < source_end)
    {
        // NOTE(holtgrew): Using value() here, used to be getValue() but
        // cannot move from const reference or proxy.
        // valueConstruct(target_begin, value(source_begin), Move());
        // TODO(holtgrew): We need a "has move constructor" metafunction to switch between move/copy constructing before we can use the line here.
#ifdef SEQAN_CXX11_STANDARD
        valueConstruct(target_begin, std::move<decltype(*source_begin)>(*source_begin));
#else
        valueConstruct(target_begin, value(source_begin), Move());
#endif
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
    _arrayConstructMoveDefault(source_begin, source_end, target_begin);
}

// ----------------------------------------------------------------------------
// Function arrayDestruct() using iterators
// ----------------------------------------------------------------------------

/*!
 * @fn arrayDestruct
 * @headerfile <seqan/basic.h>
 * @brief Destroys an array of objects.
 *
 * @signature void arrayDestruct(begin, end);
 *
 * @param[in] begin Iterator to the begin of the range that is to be destructed.
 * @param[in] end   Iterator behind the end of the range.
 *
 * This function does not deallocates the memory.
 */

template<typename TIterator1, typename TIterator2>
inline void
_arrayDestructDefault(TIterator1 begin_,
                      TIterator2 end_)
{
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
    _arrayDestructDefault(begin_, end_);
}

// ----------------------------------------------------------------------------
// Function arrayFill() using iterators
// ----------------------------------------------------------------------------

// TODO(holtgrew): What is the advantage over arrayConstruct() with prototype?

/*!
 * @fn arrayFill
 * @headerfile <seqan/basic.h>
 * @brief Assigns one object to each element of a range.
 *
 * @signature void arrayFill(begin, end, value[, parallelTag]);
 *
 * @param[in] begin       Iterator to the begin of the range that is to be filled.
 * @param[in] end         Iterator behind the end of the range.
 * @param[in] value       Argument that is assigned to all <tt>count</tt> objects in <tt>array</tt>.
 * @param[in] parallelTag Tag to enable/disable parallelism. Types: Serial, Parallel
 *
 * All objects <tt>target_begin[0]</tt> to <tt>target_begin[count-1]</tt> are set to <tt>value</tt>.
 */

// TODO(holtgrew): Redirects to fill_n. What are the exact semantics here? Do the array elements have to be initialized already? fill_n uses assignment, not copy construction!

template <typename TIterator, typename TValue>
inline void
arrayFill(TIterator begin_,
          TIterator end_,
          TValue const & value)
{
    std::fill_n(begin_, end_ - begin_, value);
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

/*!
 * @fn arrayCopyForward
 * @headerfile <seqan/basic.h>
 * @brief Copies a range of objects into another range of objects starting from the first element.
 *
 * @signature void arrayCopyForward(sourceBegin, sourceEnd, target);
 *
 * @param[in] sourceEnd   Iterator behind the last element of the source array.  <tt>sourceEnd</tt> must have the same type
 *                        as <tt>sourceBegin</tt>.
 * @param[in] sourceBegin Iterator to the first element of the source array.
 * @param[in] target      Iterator to the first element of the target array.  The target capacity should be at least as
 *                        long as the source range.
 *
 * Be careful if source and target range overlap, because in this case some source elements could be accidently
 * overwritten before they are copied.
 *
 * If there is no need for the source elements to persist, consider to use arrayMoveForward instead to improve
 * performance.
 */

template<typename TTarget, typename TSource1, typename TSource2>
inline void
_arrayCopyForwardDefault(TSource1 source_begin,
                         TSource2 source_end,
                         TTarget target_begin)
{
    std::copy(source_begin, source_end, target_begin);
}

template<typename TTarget, typename TSource1, typename TSource2>
inline void
arrayCopyForward(TSource1 source_begin,
                 TSource2 source_end,
                 TTarget target_begin)
{
    _arrayCopyForwardDefault(source_begin, source_end, target_begin);
}

// ----------------------------------------------------------------------------
// Function arrayCopyBackward() using iterators
// ----------------------------------------------------------------------------

/*!
 * @fn arrayCopyBackward
 * @headerfile <seqan/basic.h>
 * @brief Copies a range of objects into another range of objects starting from the last element.
 *
 * @signature void arrayCopyBackward(source_begin, source_end, target);
 *
 * @param[in] sourceBegin Iterator to the first element of the source array.
 * @param[in] sourceEnd   Iterator behind the last element of the source array.  <tt>sourceEnd</tt> must have the same type
 *                        as <tt>source_begin</tt>.
 * @param[in] target      Iterator to the first element of the target array.  The target capacity should be at least as
 *                        long as the source range.
 *
 * Be careful if source and target range overlap, because in this case some source elements could be accidently
 * overwritten before they are moved.
 *
 * If source and target do not overlap, consider to use the function arrayCopyForward instead that is faster in some
 * cases.
 *
 * If there is no need for the source elements to persist, consider to use arrayMoveBackward instead to improve
 * performance.
 *
 * The semantic of this function's argument <tt>target</tt> differ from the arguments of <tt>std::copy_backward</tt>.
 */

template<typename TTarget, typename TSource1, typename TSource2>
inline void
_arrayCopyBackwardDefault(TSource1 source_begin,
                          TSource2 source_end,
                          TTarget target_begin)
{
    std::copy_backward(source_begin, source_end, target_begin + (source_end - source_begin));
}

template<typename TTarget, typename TSource1, typename TSource2>
inline void
arrayCopyBackward(TSource1 source_begin,
                  TSource2 source_end,
                  TTarget target_begin)
{
    _arrayCopyBackwardDefault(source_begin, source_end, target_begin);
}

// ----------------------------------------------------------------------------
// Function arrayCopy() using iterators
// ----------------------------------------------------------------------------

/*!
 * @fn arrayCopy
 *
 * @headerfile <seqan/basic.h>
 *
 * @brief Copies a range of objects into another range of objects.
 *
 * @signature void arrayCopy(sourceBegin, sourceEnd, target);
 *
 * @param[in] sourceEnd   Iterator behind the last element of the source range.  <tt>sourceEnd</tt> must have the same type
 *                        as <tt>sourceBegin</tt>.
 * @param[in] sourceBegin Iterator to the first element of the source range.
 * @param[in] target      Iterator to the first element of the target range.The target capacity should be at least as long
 *                        as the source range.
 *
 * If source and target range do not overlap, consider to use arrayCopyForward instead to improve performance.
 *
 * If there is no need for the source elements to persist, consider to use arrayMoveForward instead to improve
 * performance.
 */

template<typename TTarget, typename TSource1, typename TSource2>
inline void arrayCopy(TSource1 source_begin,
                      TSource2 source_end,
                      TTarget target_begin)
{
    if (target_begin <= source_begin)
        arrayCopyForward(source_begin, source_end, target_begin);
    else
        arrayCopyBackward(source_begin, source_end, target_begin);
}

// ----------------------------------------------------------------------------
// Function arrayMoveForward() using iterators
// ----------------------------------------------------------------------------

/*!
 * @fn arrayMoveForward
 * @headerfile <seqan/basic.h>
 * @brief Moves a range of objects into another range of objects starting from the first element.
 *
 * @signature void arrayMoveForward(sourceBegin, sourceEnd, target);
 *
 * @param[in] sourceEnd   Iterator behind the last element of the source array.  <tt>sourceEnd</tt> must have the same type
 *                        as <tt>sourceBegin</tt>.
 * @param[in] sourceBegin Iterator to the first element of the source array.
 * @param[in] target      Iterator to the first element of the target array.  The target capacity should be at least as
 *                        long as the source range.
 *
 * The function possibly clears (but does not destroy) the source elements.  If source elements must persist, consider
 * to use arrayCopyForward instead.
 *
 * Be careful if source and target range overlap, because in this case some source elements could be accidently
 * overwritten before they are moved.
 */

template<typename TTarget, typename TSource1, typename TSource2>
inline void
_arrayMoveForwardDefault(TSource1 source_begin,
                          TSource2 source_end,
                          TTarget target_begin)
{
#ifdef SEQAN_CXX11_STANDARD
    std::move(source_begin, source_end, target_begin);
#else
    while (source_begin != source_end)
    {
        move(*target_begin, *source_begin);
        ++source_begin;
        ++target_begin;
    }
#endif
}

template<typename TTarget, typename TSource1, typename TSource2>
inline void
arrayMoveForward(TSource1 source_begin,
                 TSource2 source_end,
                 TTarget target_begin)
{
    _arrayMoveForwardDefault(source_begin, source_end, target_begin);
}

// ----------------------------------------------------------------------------
// Function arrayMoveBackward() using iterators
// ----------------------------------------------------------------------------

/*!
 * @fn arrayMoveBackward
 * @headerfile <seqan/basic.h>
 * @brief Moves a range of objects into another range of objects starting from the last element.
 *
 * @signature void arrayMoveBackward(sourceBegin, sourceEnd, target);
 *
 * @param[in] sourceEnd   Iterator behind the last element of the source array.  <tt>sourceEnd</tt> must have the same type
 *                        as <tt>sourceBegin</tt>.
 * @param[in] sourceBegin Iterator to the first element of the source array.
 * @param[in] target      Iterator to the first element of the target array.The target capacity should be at least as long
 *                        as the source range.
 *
 * The function possibly clears (but does not destroy) the source elements.  If source elements must persist, consider
 * to use arrayCopyBackward instead.
 *
 * Be careful if source and target range overlap, because in this case some source elements could be accidently
 * overwritten before they are moved.
 *
 * If source and target do not overlap, consider to use the function arrayMoveForward instead that is faster in some
 * cases.
 *
 * The semantic of this function's argument <tt>target</tt> differ from the arguments of <tt>std::copy_backward</tt>.
 */

template<typename TTarget, typename TSource1, typename TSource2>
inline void
_arrayMoveBackwardDefault(TSource1 source_begin,
                          TSource2 source_end,
                          TTarget target_begin)
{
#ifdef SEQAN_CXX11_STANDARD
    std::move_backward(source_begin, source_end, target_begin + (source_end - source_begin));
#else
    target_begin += (source_end - source_begin);
    while (source_end != source_begin) {
        --source_end;
        --target_begin;
        move(*target_begin, *source_end);
    }
#endif
}

template<typename TTarget, typename TSource1, typename TSource2>
inline void
arrayMoveBackward(TSource1 source_begin,
                  TSource2 source_end,
                  TTarget target_begin)
{
    _arrayMoveBackwardDefault(source_begin, source_end, target_begin);
}

// ----------------------------------------------------------------------------
// Function arrayMove using iterators
// ----------------------------------------------------------------------------

/*!
 * @fn arrayMove
 * @headerfile <seqan/basic.h>
 * @brief Moves a range of objects into another range of objects.
 *
 * @signature void arrayMove(sourceBegin, sourceEnd, target);
 *
 * @param[in] sourceBegin Iterator to the first element of the source range.
 * @param[in] sourceEnd   Iterator behind the last element of the source range.  <tt>sourceEnd</tt> must have the same type
 *                        as <tt>sourceBegin</tt>.
 * @param[in] target      Iterator to the first element of the target range.  The target capacity should be at least as
 *                        long as the source range.
 *
 * The function possibly clears (but does not destroy) the source elements.  If source elements must persist, consider
 * to use arrayCopy instead.
 *
 * If source and target range do not overlap, consider to use arrayMoveForward instead to improve performance.
 *
 * Don't confuse this function with the standard <tt>move</tt> function that resembles arrayCopy.
 */

template<typename TTarget, typename TSource1, typename TSource2>
inline void
arrayMove(TSource1 source_begin,
          TSource2 source_end,
          TTarget target_begin)
{
    if (target_begin <= source_begin)
        arrayMoveForward(source_begin, source_end, target_begin);
    else
        arrayMoveBackward(source_begin, source_end, target_begin);
}

// ----------------------------------------------------------------------------
// Function arrayClearSpace() using iterators
// ----------------------------------------------------------------------------

/*!
 * @fn arrayClearSpace
 * @headerfile <seqan/basic.h>
 * @brief Destroys the begin of an array and keeps the rest.
 *
 * @signature void arrayClearSpace(arrBegin, arrLength, keepFrom, moveTo);
 *
 * @param[in] arrBegin  Pointer to the first element of the array.
 * @param[in] keepFrom  Offset of the first object that will be kept.
 * @param[in] arrLength Length of the array.
 * @param[in] moveTo    Offset the first kept object will get at the end of the function.
 *
 * The objects <tt>arr[keep_from]</tt> to <tt>arr[arr_length-1]</tt> are moved to the area beginning at positions
 * <tt>move_to</tt>. All objects in <tt>arr[0]</tt> to <tt>arr[keep_from-1]</tt> are destroyed. After this function, the
 * first <tt>move_to</tt> positions of the array are free and dont contain objects.
 *
 * The array must have at least enough space to store <tt>arr_length + move_to - keep_from</tt> objects.
 *
 * The objects from <tt>arr[0]</tt> to <tt>arr[array_length-1]</tt> have to be initialized/constructed, arrays beyond
 * <tt>arr[array_length-1]</tt> are assumed not to be constructed. If this assumption is violated, memory might leak.
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
            size_t middle = array_length - (move_to - keep_from);
            arrayConstructMove(array_begin + middle, array_begin + array_length, array_begin + array_length);
            arrayMove(array_begin + keep_from, array_begin + middle, array_begin + move_to);
            arrayDestruct(array_begin, array_begin + move_to);
        } else {
            // Case 2b: We have to move-construct all target objects.
            arrayConstructMove(array_begin + keep_from, array_begin + array_length, array_begin + move_to);
            arrayDestruct(array_begin, array_begin + array_length);
        }
    } else {
        // Case 3: Move to the left.
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
    //nothing to do
}

template<typename TIterator>
inline void
_arrayConstructPointer(TIterator begin_,
                       TIterator end_,
                       False)
{
    _arrayConstructDefault(begin_, end_);
}

template <typename TValue>
inline void
arrayConstruct(TValue * begin_,
               TValue * end_)
{
    _arrayConstructPointer(begin_, end_, typename IsSimple<TValue>::Type() );
}

template <typename TValue, typename TParam>
inline void
_arrayConstructPointer(TValue * begin_,
                       TValue * end_,
                       TParam const & param_,
                       True)
{
    arrayFill(begin_, end_, static_cast<TValue>(param_));
}

template <typename TValue, typename TParam>
inline void
_arrayConstructPointer(TValue * begin_,
                       TValue * end_,
                       TParam const & param_,
                       False)
{
    _arrayConstructDefault(begin_, end_, param_);
}

template<typename TValue, typename TParam>
inline void
arrayConstruct(TValue * begin_,
               TValue * end_,
               TParam const & param_)
{
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
    arrayCopyForward(source_begin, source_end, target_begin);
}

template<typename TValueSource, typename TValueTarget>
inline void
_arrayConstructCopyPointer(TValueSource * source_begin,
                            TValueSource * source_end,
                            TValueTarget const* target_begin,
                            True)
{
    arrayCopyForward(source_begin, source_end, const_cast<TValueTarget *>(target_begin));
}

template<typename TValueSource, typename TValueTarget>
inline void
_arrayConstructCopyPointer(TValueSource * source_begin,
                            TValueSource * source_end,
                            TValueTarget * target_begin,
                            False)
{
    _arrayConstructCopyDefault(source_begin, source_end, target_begin);
}
template<typename TValueSource, typename TValueTarget>
inline void
arrayConstructCopy(TValueSource * source_begin,
                   TValueSource * source_end,
                   TValueTarget * target_begin)
{
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
    arrayMoveForward(source_begin, source_end, target_begin);
}

template<typename TValue>
inline void
_arrayConstructMovePointer(TValue * source_begin,
                            TValue * source_end,
                            TValue * target_begin,
                            False)
{
    _arrayConstructMoveDefault(source_begin, source_end, target_begin);
}

template<typename TValue>
inline void
arrayConstructMove(TValue * source_begin,
                   TValue * source_end,
                   TValue * target_begin)
{
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
    //do nothing
}

template<typename TValue>
inline void
_arrayDestructPointer(TValue * begin_,
                       TValue * end_,
                       False)
{
    _arrayDestructDefault(begin_, end_);
}

template<typename TValue>
inline void
arrayDestruct(TValue * begin_,
              TValue * end_)
{
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
    std::memmove(target_begin, source_begin, (source_end - source_begin) * sizeof(TValue));
}

template<typename TValue>
inline void
_arrayCopyForwardPointer(TValue * source_begin,
                          TValue * source_end,
                          TValue * target_begin,
                          False)
{
    _arrayCopyForwardDefault(source_begin, source_end, target_begin);
}

template<typename TValue>
inline void
arrayCopyForward(TValue * source_begin,
                 TValue * source_end,
                 TValue * target_begin)
{
    _arrayCopyForwardPointer(source_begin, source_end, target_begin, typename IsSimple<TValue>::Type() );
}

template <typename TValue>
inline void
_arrayCopyBackwardPointer(TValue * source_begin,
                           TValue * source_end,
                           TValue * target_begin,
                           True)
{
    std::memmove(target_begin, source_begin, (source_end - source_begin) * sizeof(TValue));
}

template <typename TValue>
inline void
_arrayCopyBackwardPointer(TValue * source_begin,
                           TValue * source_end,
                           TValue * target_begin,
                           False)
{
    _arrayCopyBackwardDefault(source_begin, source_end, target_begin);
}

template<typename TValue>
inline void
arrayCopyBackward(TValue * source_begin,
                  TValue * source_end,
                  TValue * target_begin)
{
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
    std::memmove(target_begin, source_begin, (source_end - source_begin) * sizeof(TValue));
}

template<typename TValue>
inline void
_arrayMoveForwardPointer(TValue * source_begin,
                          TValue * source_end,
                          TValue * target_begin,
                          False)
{
    _arrayMoveForwardDefault(source_begin, source_end, target_begin);
}

template<typename TValue>
inline void
arrayMoveForward(TValue * source_begin,
                 TValue * source_end,
                 TValue * target_begin)
{
    _arrayMoveForwardPointer(source_begin, source_end, target_begin, typename IsSimple<TValue>::Type() );
}

template <typename TValue>
inline void
_arrayMoveBackwardPointer(TValue * source_begin,
                           TValue * source_end,
                           TValue * target_begin,
                           True)
{
    std::memmove(target_begin, source_begin, (source_end - source_begin) * sizeof(TValue));
}
template <typename TValue>
inline void
_arrayMoveBackwardPointer(TValue * source_begin,
                           TValue * source_end,
                           TValue * target_begin,
                           False)
{
    _arrayMoveBackwardDefault(source_begin, source_end, target_begin);
}

template<typename TValue>
inline void
arrayMoveBackward(TValue * source_begin,
                  TValue * source_end,
                  TValue * target_begin)
{
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

#endif  // #ifndef SEQAN_INCLUDE_SEQAN_BASIC_ARRAY_CONSTRUCT_DESTRUCT_H_
