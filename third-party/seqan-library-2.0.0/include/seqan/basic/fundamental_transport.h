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
// Default implementations of the transport functions assign(), set() and
// move().
// ==========================================================================

// TODO(holtgrew): Do we want to get rid of move() and HasMoveConstructor<>? Will get rrvalues in C++11 and for everything else, swap() would be better.

#ifndef SEQAN_INCLUDE_SEQAN_BASIC_FUNDAMENTAL_TRANSPORT_H_
#define SEQAN_INCLUDE_SEQAN_BASIC_FUNDAMENTAL_TRANSPORT_H_

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

template <typename TDest, typename TSource>
void assignValue(TDest &, TSource const &);

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// ============================================================================
// Metafunctions
// ============================================================================

/*!
 * @mfn HasMoveConstructor
 * @headerfile <seqan/basic.h>
 * @brief Query whether a class has a move constructor.
 *
 * @signature HasMoveConstructor<T>::Type;
 * @signature HasMoveConstructor<T>::VALUE;
 *
 * @tparam T Type to query for availability of move constructor.
 */

template <typename T>
struct HasMoveConstructor
{
    typedef False Type;
    enum { VALUE = 0 };
};

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function assign()
// ----------------------------------------------------------------------------

/*!
 * @fn AssignableConcept#assign
 * @headerfile <seqan/basic.h>
 * @brief Assigns one object to another object.
 *
 * @signature void assign(target, source);
 *
 * @param[out] target Reference to assign to.
 * @param[in]  source Value to assign.
 *
 * Assign value of source to target.
 */

template <typename TTarget, typename TSource>
inline void
assign(TTarget & target,
       TSource & source)
{
    SEQAN_CHECKPOINT;
    target = source;
}

template <typename TTarget, typename TSource>
inline void
assign(TTarget & target,
       TSource const & source)
{
    SEQAN_CHECKPOINT;
    target = source;
}

// TODO(holtgrew): assign() for proxies should be defined in basic_proxy.h.

template <typename TSpec> class Proxy;

template<typename TTargetSpec, typename TSource>
inline void
assign(Proxy<TTargetSpec> & target,
       TSource & source)
{
    SEQAN_CHECKPOINT;
    assignValue(iter(target), source);
}

template<typename TTargetSpec, typename TSource>
inline void
assign(Proxy<TTargetSpec> & target,
       TSource const & source)
{
    SEQAN_CHECKPOINT;
    assignValue(iter(target), source);
}

// ----------------------------------------------------------------------------
// Function set()
// ----------------------------------------------------------------------------

/*!
 * @fn AssignableConcept#set
 * @headerfile <seqan/basic.h>
 * @brief Assigns one object to another object avoiding to copy contents.
 *
 * @signature set(target, source);
 *
 * @param[out] target Reference to the set to source.
 * @param[in]  source Value to set to target.
 *
 * The default implementation copies.  Types implementing AssignableConcept can implement more efficient variants.
 */

template<typename TTarget, typename TSource>
inline void
set(TTarget & target,
    TSource & source)
{
    SEQAN_CHECKPOINT;
    assign(target, source);
}

template<typename TTarget, typename TSource>
inline void
set(TTarget const & target,
    TSource & source)
{
    SEQAN_CHECKPOINT;
    assign(target, source);
}

template<typename TTarget, typename TSource>
inline void
set(TTarget & target,
    TSource const & source)
{
    SEQAN_CHECKPOINT;
    assign(target, source);
}

template<typename TTarget, typename TSource>
inline void
set(TTarget const & target,
    TSource const & source)
{
    SEQAN_CHECKPOINT;
    assign(target, source);
}

// ----------------------------------------------------------------------------
// Function move()
// ----------------------------------------------------------------------------

/*!
 * @fn AssignableConcept#move
 * @headerfile <seqan/basic.h>
 * @brief Hands over content from one object to another object.
 *
 * @signature void move(target, source);
 *
 * @param[out]     target Where to move source to.
 * @param[in,out]  source What to move to target.
 *
 * The default implementation will call @link AssignableConcept#assign @endlink and classes implementing
 * AssignableConcept can override move to provide a more efficient implementation.
 */

// TODO(holtgrew): Are all specializations necessary?

template<typename TTarget, typename TSource>
inline void
move(TTarget & target,
     TSource & source)
{
    SEQAN_CHECKPOINT;
    assign(target, source);
}

template<typename TTarget, typename TSource>
inline void
move(TTarget const & target,
     TSource & source)
{
    SEQAN_CHECKPOINT;
    assign(target, source);
}

template<typename TTarget, typename TSource>
inline void
move(TTarget & target,
     TSource const & source)
{
    SEQAN_CHECKPOINT;
    assign(target, source);
}

template<typename TTarget, typename TSource>
inline void
move(TTarget const & target,
     TSource const & source)
{
    SEQAN_CHECKPOINT;
    assign(target, source);
}

}  // namespace seqan

#endif  // #ifndef SEQAN_INCLUDE_SEQAN_BASIC_FUNDAMENTAL_TRANSPORT_H_
