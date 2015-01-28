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
// Author: Andreas Gogol-Doering <andreas.doering@mdc-berlin.de>
// ==========================================================================
// Allocator class definition and generic interface.
// ==========================================================================

// TODO(holtgrew): Perform some benchmarks and use a better malloc, e.g. tcmalloc and see whether our allocator infrastructure is worth keeping around.
// TODO(holtgrew): Rename to allocator_base.h?

#ifndef SEQAN_CORE_INCLUDE_SEQAN_BASIC_ALLOCATOR_INTERFACE_H_
#define SEQAN_CORE_INCLUDE_SEQAN_BASIC_ALLOCATOR_INTERFACE_H_

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

struct Tristate_;
typedef Tag<Tristate_> Tristate;
template <typename TValue, typename TSpec> struct Holder;

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

/**
.Tag.Allocator Usage:
..cat:Memory
..summary:The purpose of an allocated memory block.
..tag.TagAllocateUnspecified:Not specified.
..tag.TagAllocateTemp:Temporary memory. 
..tag.TagAllocateStorage:Memory for storing container content. 
..see:Function.allocate
..see:Function.deallocate
..include:seqan/basic.h
*/

// TODO(holtgrew): ANY use/difference?

struct AllocateUnspecified_;
typedef Tag<AllocateUnspecified_> TagAllocateUnspecified;

struct AllocateTemp_;
typedef Tag<AllocateTemp_> TagAllocateTemp;

struct AllocateStorage_;
typedef Tag<AllocateStorage_> TagAllocateStorage;

/**
.Class.Allocator:
..cat:Basic
..summary:Manager for allocated memory.
..signature:Allocator<TSpec>
..param.TSpec:The specializing type.
...metafunction:Metafunction.Spec
..include:basic.h
..remarks:There are two reasons for using non-trivial allocators:
...text:1. Allocators support the function @Function.Allocator#clear@ for a fast deallocation of all 
allocated memory blocks. 
...text:2. Some allocators are faster in allocating an deallocating memory.
Pool allocators like e.g. @Spec.Single Pool Allocator@ or @Spec.Multi Pool Allocator@
speed up @Function.allocate@, @Function.deallocate@, and @Function.Allocator#clear@ for
pooled memory blocks.
..include:seqan/basic.h
*/

template <typename TSpec>
struct Allocator;

///.Function.allocate.param.object.type:Class.Allocator
///.Function.allocate.class:Class.Allocator
///.Function.deallocate.param.object.type:Class.Allocator
///.Function.deallocate.class:Class.Allocator

// ============================================================================
// Metafunctions
// ============================================================================

//.Metafunction.Spec.param.T.type:Class.Allocator

template <typename TSpec>
struct Spec<Allocator<TSpec> >
{
    typedef TSpec Type;
};

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function allocate()
// ----------------------------------------------------------------------------

/**
.Function.allocate
..class:Class.Allocator
..cat:Memory
..summary:Allocates memory from heap.
..signature:allocate(object, data, count [, usage_tag])
..param.object:Allocator object.
...remarks:$object$ is conceptually the "owner" of the allocated memory.
 Objects of all types can be used as allocators. If no special behavior is implemented,
 default functions allocation/deallocation are applied that uses standard
 $new$ and $delete$ operators.
..param.count:Number of items that could be stored in the allocated memory.
...text:The type of the allocated items is given by the type of $data$.
..param.usage_tag:A tag the specifies the purpose for the allocated memory.
...value:@Tag.Allocator Usage@
..returns.param.data:Pointer to allocated memory.
...remarks:The value of this pointer is overwritten by the function.
..remarks:
...text:The function allocates at least $count*sizeof(data)$ bytes. 
 The allocated memory is large enough 
 to hold $count$ objects of type $T$, where $T *$ is type of $data$.
...note:These objects are not constructed by $allocate$.
...text:Use e.g. one of the functions @Function.valueConstruct@, @Function.arrayConstruct@, @Function.arrayConstructCopy@ or @Function.arrayFill@
to construct the objects.
A $new$ operator which is part of the C++ standard (defined in $<new>$)
 can also be used to construct objects at a given memory address.
..note:All allocated memory blocks should be deallocated by the corresponding function @Function.deallocate@.
..see:Function.deallocate
..see:Function.valueConstruct
..see:Function.arrayFill
..see:Function.arrayConstruct
..see:Function.arrayConstructCopy
..include:seqan/basic.h
*/

template <typename T, typename TValue, typename TSize>
inline void
allocate(T const & me,
         TValue * & data,
         TSize count)
{
    allocate(me, data, count, TagAllocateUnspecified());
}

template <typename T, typename TValue, typename TSize>
inline void
allocate(T & me,
         TValue * & data,
         TSize count)
{
    allocate(me, data, count, TagAllocateUnspecified());
}

template <typename T, typename TValue, typename TSize, typename TUsage>
inline void
allocate(T const &, 
         TValue * & data,
         TSize count,
         Tag<TUsage> const &)
{
//  data = (TValue *) operator new(count * sizeof(TValue));
#ifdef PLATFORM_WINDOWS_VS
    data = (TValue *) _aligned_malloc(count * sizeof(TValue), __alignof(TValue));
#else
/*#if _POSIX_C_SOURCE >= 200112L || _XOPEN_SOURCE >= 600
    const size_t align = (__alignof__(TValue) < sizeof(void*))? sizeof(void*): __alignof__(TValue);
    if (posix_memalign(&(void* &)data, align, count * sizeof(TValue)))
        data = NULL;
#else
    data = (TValue *) malloc(count * sizeof(TValue));
#endif*/
    data = (TValue *) operator new(count * sizeof(TValue));
#endif

#ifdef SEQAN_PROFILE
    if (data)
        SEQAN_PROADD(SEQAN_PROMEMORY, count * sizeof(TValue));
#endif
}

template <typename T, typename TValue, typename TSize, typename TUsage>
inline void
allocate(T &, 
         TValue * & data,
         TSize count,
         Tag<TUsage> const &)
{
//  data = (TValue *) operator new(count * sizeof(TValue));
#ifdef PLATFORM_WINDOWS_VS
    data = (TValue *) _aligned_malloc(count * sizeof(TValue), __alignof(TValue));
#else
/*#if _POSIX_C_SOURCE >= 200112L || _XOPEN_SOURCE >= 600
    const size_t align = (__alignof__(TValue) < sizeof(void*))? sizeof(void*): __alignof__(TValue);
    if (posix_memalign(&(void* &)data, align, count * sizeof(TValue)))
        data = NULL;
#else
    data = (TValue *) malloc(count * sizeof(TValue));
#endif
*/  data = (TValue *) operator new(count * sizeof(TValue));
#endif

#ifdef SEQAN_PROFILE
    if (data)
        SEQAN_PROADD(SEQAN_PROMEMORY, count * sizeof(TValue));
#endif
}

// ----------------------------------------------------------------------------
// Function deallocate()
// ----------------------------------------------------------------------------

/**
.Function.deallocate
..class:Class.Allocator
..cat:Memory
..summary:Deallocates memory.
..signature:deallocate(object, data, count [, usage_tag])
..param.object:Allocator object.
...remarks:$object$ is conceptually the "owner" of the allocated memory.
 Objects of all types can be used as allocators. If no special behavior is implemented,
 default functions allocation/deallocation are applied that uses standard
 $new$ and $delete$ operators.
..param.data:Pointer to allocated memory that was allocated by $allocate$.
..param.count:Number of items that could be stored in the allocated memory.
..param.usage_tag:A tag the specifies the purpose for the allocated memory.
...value:@Tag.Allocator Usage@
..remarks:
...text:The values for $object$, $count$ and $usage_tag$ should be the same that was 
used when $allocate$ was called. The value of $data$ should be the same that was
returned by $allocate$.
...note:$deallocate$ does not destruct objects.
...text:Use e.g. one of the functions @Function.valueDestruct@ or @Function.arrayDestruct@ to destruct the objects.
$delete$ and $delete []$ operators which are part of the C++ standard (defined in $<new>$)
 can also be used to destruct objects at a given memory address.
..see:Function.valueDestruct
..see:Function.arrayDestruct
..include:seqan/basic.h
*/

template <typename T, typename TValue, typename TSize>
inline void 
deallocate(T const & me, 
           TValue * data, 
           TSize const count)
{
    deallocate(me, data, count, TagAllocateUnspecified());
}

template <typename T, typename TValue, typename TSize>
inline void 
deallocate(T & me, 
           TValue * data, 
           TSize const count)
{
    deallocate(me, data, count, TagAllocateUnspecified());
}

template <typename T, typename TValue, typename TSize, typename TUsage>
inline void
deallocate(
    T const & /*me*/,
    TValue * data,
#ifdef SEQAN_PROFILE 
    TSize count,
#else
    TSize,
#endif
    Tag<TUsage> const)
{
#ifdef SEQAN_PROFILE
    if (data && count)  // .. to use count if SEQAN_PROFILE is not defined
        SEQAN_PROSUB(SEQAN_PROMEMORY, count * sizeof(TValue));
#endif
//  operator delete ((void *) data);
#ifdef PLATFORM_WINDOWS_VS
    _aligned_free((void *) data);
#else
//  free((void *) data);
    operator delete ((void *) data);
#endif
}

template <typename T, typename TValue, typename TSize, typename TUsage>
inline void 
deallocate(
    T & /*me*/,
    TValue * data,
#ifdef SEQAN_PROFILE 
    TSize count,
#else
    TSize,
#endif
    Tag<TUsage> const)
{
#ifdef SEQAN_PROFILE
    if (data && count)  // .. to use count if SEQAN_PROFILE is not defined
        SEQAN_PROSUB(SEQAN_PROMEMORY, count * sizeof(TValue));
#endif
//  operator delete ((void *) data);
#ifdef PLATFORM_WINDOWS_VS
    _aligned_free((void *) data);
#else
//  free((void *) data);
    operator delete ((void *) data);
#endif
}

}  // namespace seqan

#endif  // #ifndef SEQAN_CORE_INCLUDE_SEQAN_BASIC_ALLOCATOR_INTERFACE_H_
