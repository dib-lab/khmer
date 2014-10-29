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
// Bit-packed pair specialization.
// ==========================================================================

#ifndef SEQAN_CORE_INCLUDE_SEQAN_BASIC_PAIR_BIT_PACKED_H_
#define SEQAN_CORE_INCLUDE_SEQAN_BASIC_PAIR_BIT_PACKED_H_

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

/**
.Spec.Bit Packed Pair:
..cat:Aggregates
..general:Class.Pair
..summary:Stores two arbitrary objects. Saves memory by packing bits with bit fields.
..signature:Pair<T1, T2, BitPacked<BITSIZE1, BITSIZE2> >
..param.T1:The type of the first object.
..param.T2:The type of the second object.
..param.BITSIZE1:Number of bits to store $T1$.
..param.BITSIZE2:Number of bits to store $T2$.
..notes:Useful for external storage.
..remarks:Memory access could be slower. Direct access to members by pointers is not allowed.
..remarks:Functions $value()$ is not implemented yet since there it would require using a proxy. Use $getValue()$, $assignValue()$, $moveValue()$, $setValue()$ instead.
..include:seqan/basic.h
.Memfunc.Pair#Pair.class:Spec.Bit Packed Pair
.Memvar.Pair#i1.class:Spec.Bit Packed Pair
.Memvar.Pair#i2.class:Spec.Bit Packed Pair
*/

#ifdef PLATFORM_WINDOWS
    #pragma pack(push,1)
#endif
template <typename T1, typename T2, unsigned BITSIZE1, unsigned BITSIZE2>
struct Pair<T1, T2, BitPacked<BITSIZE1, BITSIZE2> >
{
    // ------------------------------------------------------------------------
    // Members
    // ------------------------------------------------------------------------

    T1 i1:BITSIZE1;
    T2 i2:BITSIZE2;

    // ------------------------------------------------------------------------
    // Constructors
    // ------------------------------------------------------------------------

    inline Pair() : i1(T1()), i2(T2()) {}

    inline Pair(Pair const & _p) : i1(_p.i1), i2(_p.i2) {}

    inline Pair(T1 const & _i1, T2 const & _i2) : i1(_i1), i2(_i2) {}

    template <typename T1_, typename T2_, typename TSpec__>
    // TODO(holtgrew): explicit?
    inline Pair(Pair<T1_, T2_, TSpec__> const &_p)
            : i1(getValueI1(_p)), i2(getValueI2(_p)) {}
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

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function set().
// ----------------------------------------------------------------------------

// References to members to packed structs do not work.  Always copy.

template <typename T1, typename T2, unsigned BITSIZE1, unsigned BITSIZE2>
inline void
set(Pair<T1, T2, BitPacked<BITSIZE1, BITSIZE2> > & p1, Pair<T1, T2, BitPacked<BITSIZE1, BITSIZE2> > & p2)
{
    p1 = p2;
}

// ----------------------------------------------------------------------------
// Function move().
// ----------------------------------------------------------------------------

// References to members to packed structs do not work.  Always copy.

template <typename T1, typename T2, unsigned BITSIZE1, unsigned BITSIZE2>
inline void
move(Pair<T1, T2, BitPacked<BITSIZE1, BITSIZE2> > & p1, Pair<T1, T2, BitPacked<BITSIZE1, BITSIZE2> > & p2)
{
    p1 = p2;
}

// -----------------------------------------------------------------------
// Function setValueIX()
// -----------------------------------------------------------------------

// Cannot be setValue with index since T1 can be != T2.

template <typename T1, typename T2, typename T, unsigned BITSIZE1, unsigned BITSIZE2>
inline void setValueI1(Pair<T1, T2, BitPacked<BITSIZE1, BITSIZE2> > & pair, T const & _i)
{
    pair.i1 = _i;
}

template <typename T1, typename T2, typename T, unsigned BITSIZE1, unsigned BITSIZE2>
inline void setValueI2(Pair<T1, T2, BitPacked<BITSIZE1, BITSIZE2> > & pair, T const & _i)
{
    pair.i2 = _i;
}

// -----------------------------------------------------------------------
// Function moveValueIX()
// -----------------------------------------------------------------------

// Cannot be moveValue with index since T1 can be != T2.

template <typename T1, typename T2, typename T, unsigned BITSIZE1, unsigned BITSIZE2>
inline void moveValueI1(Pair<T1, T2, BitPacked<BITSIZE1, BITSIZE2> > & pair, T & _i)
{
    pair.i1 = _i;
}

template <typename T1, typename T2, typename T, unsigned BITSIZE1, unsigned BITSIZE2>
inline void moveValueI2(Pair<T1, T2, BitPacked<BITSIZE1, BITSIZE2> > & pair, T & _i)
{
    pair.i2 = _i;
}

}  // namespace seqan

#endif  // #ifndef SEQAN_CORE_INCLUDE_SEQAN_BASIC_PAIR_BIT_PACKED_H_
