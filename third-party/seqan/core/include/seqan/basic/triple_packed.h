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
// Packed triple specialization.
// ==========================================================================

#ifndef SEQAN_CORE_INCLUDE_SEQAN_TRIPLE_PACKED_H_
#define SEQAN_CORE_INCLUDE_SEQAN_TRIPLE_PACKED_H_

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

/*!
 * @class PackedTriple
 * @headerfile <seqan/basic.h>
 * @extends Triple
 * @brief Disable memory alignment to store memory.
 *
 * @signature template <typename T1, typename T2, typename T3>
 *            class Triple<T1, T2, T3, Pack>;
 *
 * @tparam T1 The first type of the Triple.
 * @tparam T2 The second type of the Triple.
 * @tparam T3 The third type of the Triple.
 */

/**
.Spec.Packed Triple:
..cat:Aggregates
..general:Class.Triple
..summary:Stores three arbitrary objects. Saves memory by disabling memory alignment.
..signature:Triple<T1, T2, T3, Pack>
..param.T1:The type of the first object.
..param.T2:The type of the second object.
..param.T3:The type of the third object.
..notes:Useful for external storage.
..remarks:Memory access could be slower. Direct access to members by pointers is not allowed on all platforms.
..include:seqan/basic.h

.Memfunc.Triple#Triple.class:Spec.Packed Triple
.Memvar.Triple#i1.class:Spec.Packed Triple
.Memvar.Triple#i2.class:Spec.Packed Triple
.Memvar.Triple#i3.class:Spec.Packed Triple
*/

#ifdef PLATFORM_WINDOWS
    #pragma pack(push,1)
#endif
template <typename T1, typename T2, typename T3>
struct Triple<T1, T2, T3, Pack>
{
    // -----------------------------------------------------------------------
    // Members
    // -----------------------------------------------------------------------

    T1 i1;
    T2 i2;
    T3 i3;

    // -----------------------------------------------------------------------
    // Constructors
    // -----------------------------------------------------------------------

    inline Triple() : i1(T1()), i2(T2()), i3(T3()) {}
    
    inline Triple(Triple const &_p)
            : i1(_p.i1), i2(_p.i2), i3(_p.i3) {}
    
    inline Triple(T1 const &_i1, T2 const &_i2, T3 const &_i3)
            : i1(_i1), i2(_i2), i3(_i3) {}
    
    template <typename T1_, typename T2_, typename T3_, typename TSpec__>
    inline Triple(Triple<T1_, T2_, T3_, TSpec__> const & _p)
            : i1(getValueI1(_p)), i2(getValueI2(_p)), i3(getValueI3(_p)) {}
}
#ifndef PLATFORM_WINDOWS
    __attribute__((packed))
#endif
    ;
#ifdef PLATFORM_WINDOWS
    #pragma pack(pop)
#endif

// ============================================================================
// Metafunctions
// ============================================================================

template <typename T1, typename T2, typename T3, typename TSpec>
struct MakePacked< Triple<T1, T2, T3, TSpec> >
{
    typedef Triple<T1, T2, T3, Pack> Type;
};

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function set().
// ----------------------------------------------------------------------------

// References to members to packed structs do not work.  Always copy.

template <typename T1, typename T2, typename T3, typename T>
inline void
set(Triple<T1, T2, T3, Pack> & t1, Triple<T1, T2, T3, Pack> & t2)
{
    t1 = t2;
}

// ----------------------------------------------------------------------------
// Function move().
// ----------------------------------------------------------------------------

// References to members to packed structs do not work.  Always copy.

template <typename T1, typename T2, typename T3, typename T>
inline void
move(Triple<T1, T2, T3, Pack> & t1, Triple<T1, T2, T3, Pack> & t2)
{
    t1 = t2;
}

// -----------------------------------------------------------------------
// Function setValueIX()
// -----------------------------------------------------------------------

// References to members to packed structs do not work.  Always copy.

template <typename T1, typename T2, typename T3, typename T>
inline void setValueI1(Triple<T1, T2, T3, Pack> & triple, T const & _i)
{
    triple.i1 = _i;
}

template <typename T1, typename T2, typename T3, typename T>
inline void setValueI2(Triple<T1, T2, T3, Pack> & triple, T const & _i)
{
    triple.i2 = _i;
}

template <typename T1, typename T2, typename T3, typename T>
inline void setValueI3(Triple<T1, T2, T3, Pack> & triple, T const & _i)
{
    triple.i3 = _i;
}

// -----------------------------------------------------------------------
// Function moveValueIX()
// -----------------------------------------------------------------------

// References to members to packed structs do not work.  Always copy.

template <typename T1, typename T2, typename T3, typename T>
inline void moveValueI1(Triple<T1, T2, T3, Pack> & triple, T const & _i)
{
    triple.i1 = _i;
}

template <typename T1, typename T2, typename T3, typename T>
inline void moveValueI2(Triple<T1, T2, T3, Pack> & triple, T const & _i)
{
    triple.i2 = _i;
}

template <typename T1, typename T2, typename T3, typename T>
inline void moveValueI3(Triple<T1, T2, T3, Pack> & triple, T const & _i)
{
    triple.i3 = _i;
}

}  // namespace seqan

#endif  // #ifndef SEQAN_CORE_INCLUDE_SEQAN_TRIPLE_PACKED_H_
