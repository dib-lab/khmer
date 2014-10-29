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
// Implementation of the Holder base class.
// ==========================================================================

#ifndef SEQAN_BASIC_HOLDER_BASE_H_
#define SEQAN_BASIC_HOLDER_BASE_H_

// By default, disable holders to pointers, this is used/tested nowhere and does probably not work.
#ifndef SEQAN_ENABLE_POINTER_HOLDER
#define SEQAN_ENABLE_POINTER_HOLDER 0
#endif  //#ifndef SEQAN_ENABLE_POINTER_HOLDER

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// ----------------------------------------------------------------------------
// Class Holder
// ----------------------------------------------------------------------------

/**
.Class.Holder:
..cat:Basic
..summary:Manages relationship to another object.
..signature:Holder<TValue, TSpec>
..param.TValue:Type of the managed object.
...metafunction:Metafunction.Value
..param.TSpec:The specializing type.
...metafunction:Metafunction.Spec
...default:$Tristate$
..remarks.text:The main purpose of this class is to facilitate the handling of
member objects. If we want class $A$ to be dependent on or the owner of another object of class $B$, 
then we add a data member of type $Holder<B>$ to $A$. 
$Holder$ offers some useful access functions and stores the kind of relationship between $A$ and $B$.
..include:seqan/basic.h

.Memfunc.Holder:
..class:Class.Holder
..summary:Constructor
..signature:Holder<TValue, TSpec>()
..signature:Holder<TValue, TSpec>(holder)
..signature:Holder<TValue, TSpec>(value)
..param.holder:Another holder object.
..param.value:An object of type $TValue$.
..remarks.text:
The default constructor creates a holder that is in state 'empty'.
If a $value$ is passed to the constructor, the holder will be in state 'dependent'.
 */

// Tag for default Holder specialization.
struct Tristate_;
typedef Tag<Tristate_> Tristate;

// Tag for default Simple specialization.
struct Simple_;
typedef Tag<Simple_> Simple;

template <typename TValue, typename TSpec = Tristate>
struct Holder;

// ============================================================================
// Metafunctions
// ============================================================================

// ----------------------------------------------------------------------------
// Metafunction Value
// ----------------------------------------------------------------------------

///.Metafunction.Value.param.T.type:Class.Holder

template <typename TValue, typename TSpec>
struct Value<Holder<TValue, TSpec> >
{
    typedef TValue Type;
};

template <typename TValue, typename TSpec>
struct Value<Holder<TValue, TSpec> const>
{
    typedef TValue Type;
};

#if SEQAN_ENABLE_POINTER_HOLDER
// TODO(holtgrew): What about holders on pointers?
template <typename TValue, typename TSpec>
struct Value<Holder<TValue * const, TSpec> >
{
    typedef TValue * Type;
};
#endif  // #if SEQAN_ENABLE_POINTER_HOLDER

// ----------------------------------------------------------------------------
// Metafunction Spec
// ----------------------------------------------------------------------------

///.Metafunction.Spec.param.T.type:Class.Holder
///.Metafunction.Spec.class:Class.Holder

template <typename TValue, typename TSpec>
struct Spec<Holder<TValue, TSpec> >
{
    typedef TSpec Type;
};

template <typename TValue, typename TSpec>
struct Spec<Holder<TValue, TSpec> const>
{
    typedef TSpec Type;
};

// ----------------------------------------------------------------------------
// Metafunction Reference
// ----------------------------------------------------------------------------

///.Metafunction.Reference.param.T.type:Class.Holder
///.Metafunction.Reference.class:Class.Holder

template <typename TValue, typename TSpec>
struct Reference<Holder<TValue, TSpec> >
{
    typedef typename Value<Holder<TValue, TSpec> >::Type & Type;
};

template <typename TValue, typename TSpec>
struct Reference< Holder<TValue, TSpec> const>
{
    typedef typename Value<Holder<TValue, TSpec> const>::Type & Type;
};

#if SEQAN_ENABLE_POINTER_HOLDER
template <typename TValue, typename TSpec>
struct Reference<Holder<TValue *, TSpec> const>
{
    typedef typename Value<Holder<TValue *, TSpec> const>::Type const & Type;
};
#endif  // #if SEQAN_ENABLE_POINTER_HOLDER

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function getValue()
// ----------------------------------------------------------------------------

// TODO(holtgrew): We would rather have only one here.

template <typename TValue, typename TSpec>
inline typename GetValue<Holder<TValue, TSpec> >::Type
getValue(Holder<TValue, TSpec> const & holder)
{
    return value(holder);
}

template <typename TValue, typename TSpec>
inline typename GetValue<Holder<TValue, TSpec> >::Type
getValue(Holder<TValue, TSpec> & holder)
{
    return value(holder);
}

}  // namespace seqan

#endif  // #ifndef SEQAN_BASIC_HOLDER_BASE_H_
