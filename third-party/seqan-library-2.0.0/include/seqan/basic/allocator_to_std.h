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
// Wrapper that adapts SeqAn allocators to STL allocators.
// ==========================================================================

// TODO(holtgrew): Rename STD to STL?
// TODO(holtgrew): Rename to allocator_to_stl.h, remove basic_ prefix of all other allocator headers.

#ifndef SEQAN_INCLUDE_SEQAN_BASIC_ALLOCATOR_TO_STD_H_
#define SEQAN_INCLUDE_SEQAN_BASIC_ALLOCATOR_TO_STD_H_

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

/*!
 * @class ToStdAllocator
 * @headerfile <seqan/basic.h>
 * @brief Emulates standard conform allocator.
 *
 * @signature template <typename THost, typename TValue>
 *            class ToStdAllocator;
 *
 * @tparam TValue Type of allocated items.
 * @tparam THost  Type of the host allocator object.This object is used to call @link Allocator#allocate @endlink and
 *                @link Allocator#deallocate @endlink.
 *
 * The member functions <tt>allocate</tt> and <tt>deallocate</tt> of <tt>ToStdAllocator</tt> call the (globale)
 * functions @link Allocator#allocate @endlink and @link Allocator#deallocate @endlink, respectively. The globale
 * functions get an allocator object as their first arguments. This allocator object is not the <tt>ToStdAllocator</tt>
 * object itself, but the host object that was given to the constructor.
 *
 * @fn ToStdAllocator#ToStdAllocator
 * @brief Constructor
 *
 * @signature ToStdAllocator::ToStdAllocator(host);
 *
 * @param[in] host The host object that is used as allocator for @link Allocator#allocate @endlink and @link
 *                 Allocator#deallocate @endlink.
 */

template <typename THost, typename TValue>
struct ToStdAllocator
{
    typedef TValue value_type; // nolint
    typedef value_type * pointer; // nolint
    typedef value_type & reference; // nolint
    typedef value_type const * const_pointer; // nolint
    typedef value_type const & const_reference; // nolint

    typedef size_t size_type; // nolint
    typedef ptrdiff_t difference_type; // nolint

    ToStdAllocator(THost & host): m_host(& host)
    {}

    template <typename TValue2>
    ToStdAllocator(ToStdAllocator<THost, TValue2> const & alloc)
            : m_host(alloc.m_host)
    {}

    ToStdAllocator & operator= (ToStdAllocator const & alloc)
    {
        m_host = alloc.m_host;
        return *this;
    }

    pointer allocate(size_type count)
    {
        value_type * ptr;
        seqan::allocate(*m_host, ptr, count);
        return pointer(ptr);
    }

    pointer allocate(size_type count, const void *)
    {
        value_type * ptr;
        seqan::allocate(*m_host, ptr, count);
        return pointer(ptr);
    }

    void deallocate(pointer data, size_type count)
    {
        seqan::deallocate(*m_host, data, count);
    }

    void construct(pointer ptr, const_reference data)
    {
        new(ptr) TValue(data);
    }

    void destroy(pointer ptr)
    {
        ptr->~TValue();
    }

    pointer address(reference value) const
    {
        return (&value);
    }

    const_pointer address(const_reference value) const
    {
        return (&value);
    }

    size_type max_size() const
    {
        return ~0UL / sizeof(value_type);
    }

    template <class TValue2>
    struct rebind // nolint
    {
        typedef ToStdAllocator<THost, TValue2> other; // nolint
    };

    template <typename THost2, typename TValue2>
    friend
    struct ToStdAllocator;

private:
    THost * m_host;
};

// ============================================================================
// Metafunctions
// ============================================================================

// ----------------------------------------------------------------------------
// Metafunction ToStdAllocator()
// ----------------------------------------------------------------------------

template <typename T, typename TData>
struct StdAllocator
{
    typedef ToStdAllocator<T, TData> Type;
};

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function host()
// ----------------------------------------------------------------------------

/*!
 * @fn ToStdAllocator#host
 * @headerfile <seqan/basic.h>
 * @brief The object a given object depends on.
 *
 * @signature THost host(allocator);
 *
 * @param[in] allocator The allocator to query.
 *
 * @return THost The host object.
 */

template <typename THost, typename TValue>
THost &
host(ToStdAllocator<THost, TValue> & me)
{
   return *me.m_host;
}

}  // namespace seqan

#endif  // #ifndef SEQAN_INCLUDE_SEQAN_BASIC_ALLOCATOR_TO_STD_H_
