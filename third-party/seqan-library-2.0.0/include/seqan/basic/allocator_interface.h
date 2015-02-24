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
// Author: Andreas Gogol-Doering <andreas.doering@mdc-berlin.de>
// ==========================================================================
// Allocator class definition and generic interface.
// ==========================================================================

// TODO(holtgrew): Perform some benchmarks and use a better malloc, e.g. tcmalloc and see whether our allocator infrastructure is worth keeping around.
// TODO(holtgrew): Rename to allocator_base.h?

#ifndef SEQAN_INCLUDE_SEQAN_BASIC_ALLOCATOR_INTERFACE_H_
#define SEQAN_INCLUDE_SEQAN_BASIC_ALLOCATOR_INTERFACE_H_

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

/*!
 * @defgroup AllocatorUsageTags Allocator Usage Tags
 * @brief The purpose of an allocated memory block.
 *
 * @tag AllocatorUsageTags#TagAllocateUnspecified
 * @headerfile <seqan/basic.h>
 * @brief Not specified.
 *
 * @tag AllocatorUsageTags#TagAllocateTemp
 * @headerfile <seqan/basic.h>
 * @brief Temporary memory.
 *
 * @tag AllocatorUsageTags#TagAllocateStorage
 * @headerfile <seqan/basic.h>
 * @brief Memory for storing container content.
 */

// TODO(holtgrew): ANY use/difference?

struct AllocateUnspecified_;
typedef Tag<AllocateUnspecified_> TagAllocateUnspecified;

struct AllocateTemp_;
typedef Tag<AllocateTemp_> TagAllocateTemp;

struct AllocateStorage_;
typedef Tag<AllocateStorage_> TagAllocateStorage;

/*!
 * @class Allocator
 * @headerfile <seqan/basic.h>
 * @brief Manager for allocated memory.
 *
 * @signature template <typename TSpec>
 *            class Allocator;
 *
 * @tparam TSpec The specializing type.
 *
 * @section Remarks
 *
 * There are two reasons for using non-trivial allocators:
 *
 * <ol>
 *   <li>Allocators support the function @link Allocator#clear @endlink for a fast deallocation of all allocated
 *       memory blocks.</li>
 *   <li>Some allocators are faster in allocating an deallocating memory.  Pool allocators like e.g.
 *       @link SinglePoolAllocator @endlink or @link MultiPoolAllocator @endlink speed up
 *       @link Allocator#allocate @endlink, @link Allocator#deallocate * @endlink, and
 *       @link Allocator#clear @endlink for pooled memory blocks.</li>
 * </ol>
 */

template <typename TSpec>
struct Allocator;

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

/*!
 * @fn Allocator#allocate
 * @headerfile <seqan/basic.h>
 * @brief Allocates memory from heap.
 *
 * @signature void allocate(allocator, data, count[, usageTag]);
 *
 * @param[in]     count     Number of items that could be stored in the allocated memory.  The type of the allocated
 *                          items is given by the type of <tt>data</tt>.
 * @param[in]     usageTag  A tag the specifies the purpose for the allocated memory.  Values:
 *                          @link AllocatorUsageTags @endlink.
 * @param[in,out] allocator Allocator object.  <tt>allocator</tt> is conceptually the "owner" of the allocated
 *                          memory.  Objects of all types can be used as allocators.  If no special behavior is
 *                          implemented,  default functions allocation/deallocation are applied that uses standard
 *                          <tt>new</tt> and <tt>delete</tt> operators. Types: Allocator
 *
 * @section Remarks
 *
 * The function allocates at least <tt>count*sizeof(data)</tt> bytes.  The allocated memory is large enough  to hold
 * <tt>count</tt> objects of type <tt>T</tt>, where <tt>T *</tt> is type of <tt>data</tt>.
 *
 * These objects are not constructed by <tt>allocate</tt>.
 *
 * Use e.g. one of the functions @link valueConstruct @endlink, @link arrayConstruct @endlink, @link arrayConstructCopy
 * @endlink or @link arrayFill @endlink to construct the objects. A <tt>new</tt> operator which is part of the C++
 * standard (defined in <tt>&lt;new&gt;</tt>)  can also be used to construct objects at a given memory address.
 *
 * @section Remarks
 *
 * All allocated memory blocks should be deallocated by the corresponding function @link Allocator#deallocate @endlink.
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

/*!
 * @fn Allocator#deallocate
 * @headerfile <seqan/basic.h>
 * @brief Deallocates memory.
 *
 * @signature void deallocate(object, data, count[, usageTag])
 *
 * @param[in,out] object Allocator object.<tt>object</tt> is conceptually the "owner" of the allocated memory.
 *                       Objects of all types can be used as allocators.  If no special behavior is implemented,
 *                       default functions allocation/deallocation are applied that uses standard  <tt>new</tt>
 *                       and <tt>delete</tt> operators.  Types: Allocator
 * @param[out]     data  Pointer to allocated memory that was allocated by <tt>allocate</tt>.
 * @param[in]     count    Number of items that could be stored in the allocated memory.
 * @param[in]     usageTag A tag the specifies the purpose for the allocated memory.
 *                         Values: @link AllocatorUsageTags @endlink.
 *
 * The values for <tt>object</tt>, <tt>count</tt> and <tt>usageTag</tt> should be the same that was used when
 * <tt>allocate</tt> was called. The value of <tt>data</tt> should be the same that was returned by <tt>allocate</tt>.
 *
 * <tt>deallocate</tt> does not destruct objects.
 *
 * Use e.g. one of the functions @link valueDestruct @endlink or @link arrayDestruct @endlink to destruct the objects.
 * <tt>delete</tt> and <tt>delete []</tt> operators which are part of the C++ standard (defined in <tt>&lt;new&gt;</tt>)
 * can also be used to destruct objects at a given memory address.
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

#endif  // #ifndef SEQAN_INCLUDE_SEQAN_BASIC_ALLOCATOR_INTERFACE_H_
