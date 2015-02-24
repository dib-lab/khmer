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

/*!
 * @class Holder
 * @headerfile <seqan/basic.h>
 * @brief Manages relationship to another object.
 *
 * @signature template <typename TValue[, typename TSpec>
 *            class Holder;
 *
 * @tparam TSpec  The specializing type. Default: <tt>Tristate</tt>
 * @tparam TValue Type of the managed object.
 *
 * @section Remarks
 *
 * The main purpose of this class is to facilitate the handling of member objects. If we want class <tt>A</tt> to be
 * dependent on or the owner of another object of class <tt>B</tt>, then we add a data member of type <tt>Holder&lt;B&gt;</tt>
 * to <tt>A</tt>.  <tt>Holder</tt> offers some useful access functions and stores the kind of relationship between
 * <tt>A</tt> and <tt>B</tt>.
 *
 * @see Holder#create
 * @see Holder#setValue
 */

/*
 * @fn Holder::Holder
 * @brief Constructor.
 *
 * You can construct a Holder using the default constructor, copy from another Holder, or directly with a value to hold.
 *
 * @signature Holder::Holder();
 * @signature Holder::Holder(other);
 * @signature Holder::Holder(value);
 *
 * @param[in] other Other Holder to copy from.
 * @param[in] value The value to hold.
 *
 * The main purpose of this class is to facilitate the handling of member objects.  If we want a class <tt>A</tt> to be
 * dependent on or the owner or another object of class <tt>B</tt> then we add a data member of type
 * <tt>Holder&lt;B&gt;</tt> to <tt>A</tt>.  Holder offers some useful access functions and stores the kind of
 * relationship between <tt>A</tt> and <tt>B</tt>.
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

/*!
 * @mfn Holder#GetValue
 * @brief Return get-value type of Holder.
 *
 * @signature GetValue<THolder>::Type;
 *
 * @tparam THolder The Holder to query for its value type.
 *
 * @return Type The get-value type for its holder.
 */

// ----------------------------------------------------------------------------
// Metafunction Value
// ----------------------------------------------------------------------------

/*!
 * @mfn Holder#Value
 * @brief Return value type of Holder.
 *
 * @signature Value<THolder>::Type;
 *
 * @tparam THolder The Holder to query for its value type.
 *
 * @return Type The value type for its holder.
 */

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

/*!
 * @mfn Holder#Spec
 * @brief Return the specialization tag for a Holder.
 *
 * @signature Spec<THolder>::Type;
 *
 * @tparam THolder The Holder to query for its value type.
 *
 * @return Type The resulting specialization tag.
 */

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

/*!
 * @mfn Holder#Reference
 * @brief Return the reference type of a Holder.
 *
 * @signature Reference<THolder>::Type;
 *
 * @tparam THolder The holder to query for its reference type.
 *
 * @return Type The resulting reference type.
 */

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
// Function create()
// ----------------------------------------------------------------------------

/*!
 * @fn Holder#create
 * @brief Makes an object to owner of its content.
 *
 * @signature void create(holder[, object]);
 *
 * @param[in,out] holder The Holder to create the object of.
 * @param[in]     object Object from which a copy is made and stored in <tt>holder</tt>.
 *
 * After this operation, <tt>holder</tt> will be in state 'owner'.  If <tt>object</tt> is specified, <tt>holder</tt>
 * will hold a copy of <tt>object</tt> at the end of this function.  If <tt>object</tt> is not specified, the action
 * depends on the former state of <tt>holder</tt>:
 *
 * <ul>
 *   <li>If the state of <tt>holder</tt> was 'empty', a new object is default constructed and stored into
 *       <tt>holder</tt>.</li>
 *   <li>If the state of <tt>holder</tt> was 'dependent', a copy of the former object is made and stored into
 *       <tt>holder</tt>.</li>
 *   <li>If the state of <tt>holder</tt> was already 'owner', nothing happens.</li>
 * </ul>
 *
 * It is guaranteed, that after calling this function <tt>source</tt> and <tt>target</tt> can be used independently.
 */

// ----------------------------------------------------------------------------
// Function detach()
// ----------------------------------------------------------------------------

/*!
 * @fn Holder#detach
 * @brief Makes an object independent from other objects.

 * @signature void detach(holder);
 *
 * @param[in,out] holder The Holder to detach.
 *
 * @section Remarks
 *
 * After this function, <tt>holder</tt> does not depends from any other entity outside of <tt>holder</tt>, like a source
 * or a host, and dependent(holer) returns <tt>false</tt>
 */

// ----------------------------------------------------------------------------
// Function setValue()
// ----------------------------------------------------------------------------

/*!
 * @fn Holder#setValue
 * @brief Makes holder dependent.
 *
 * @signature void setValue(holder, object);
 *
 * @param[in,out] holder A holder object. Types: Holder
 * @param[in]     object Object from which <tt>holder</tt> will be dependent.
 *
 * After this operation, <tt>holder</tt> will be dependent in state 'dependent'.
 */

// ----------------------------------------------------------------------------
// Function empty()
// ----------------------------------------------------------------------------

/*!
 * @fn Holder#empty
 * @brief Test a Holder for being empty.
 *
 * @signature bool empty(holder);
 *
 * @param[in] holder A Holder.
 *
 * @return bool <tt>true</tt> if <tt>holder</tt> contains no elements, otherwise <tt>false</tt>.
 *
 * @section Remarks
 *
 * <tt>empty(x)</tt> is guaranteed to be at least as fast as <tt>length(me) == 0</tt>, but can be significantly faster
 * in some cases.
 *
 * @see HostedConcept#emptyHost
 */

// ----------------------------------------------------------------------------
// Function assignValue()
// ----------------------------------------------------------------------------

/*!
 * @fn Holder#assignValue
 *
 * @headerfile <seqan/basic.h>
 * @headerfile <seqan/sequence.h>
 *
 * @brief Assigns value to item.
 *
 * @signature void assignValue(object, value);
 *
 * @param[in,out] object An object that holds a value or points to a value. Types:
 *                       Holder, Iter Concepts: Concept.BasicOutputIteratorConcept
 * @param[in]     value  A value that is assigned to the item <tt>object</tt> holds or points to.
 *
 * This function is similar to @link AssignableConcept#assign @endlink. The difference is, that
 * <tt>assignValue</tt> just changes a value stored in <tt>object</tt> or the
 * value <tt>object</tt> points to, while @link AssignableConcept#assign @endlink changes the
 * whole object.
 *
 * If <tt>object</tt> is a container (that is <tt>pos</tt> is not specified),
 * the whole content of <tt>object</tt> is replaced by <tt>value</tt>.
 *
 * If <tt>value</tt> is not used again after calling this function,     then
 * consider to use @link Holder#moveValue @endlink that could be faster in some cases
 * instead.
 *
 * @see AssignableConcept#assign
 * @see Holder#moveValue
 */

// ----------------------------------------------------------------------------
// Function dependent()
// ----------------------------------------------------------------------------

/*
 * @fn Holder#dependent
 * @headerfile <seqan/basic.h>
 * @brief Test whether Holder depends on other objects.
 *
 * @signature bool dependent(holder);
 *
 * @param[in] holder The Holder to query;
 *
 * @return bool <tt>true</tt> if <tt>object</tt> depends one some other object, <tt>false</tt> otherwise.
 *
 * An object "<tt>a</tt>" depends on another object "<tt>b</tt>", if changing "<tt>b</tt>" can invalidate "<tt>a</tt>";
 * especially the destruction of "<tt>b</tt>" invalidates "<tt>a</tt>".
 */

// ----------------------------------------------------------------------------
// Function clear()
// ----------------------------------------------------------------------------

/*!
 * @fn Holder#clear
 * @brief Clear/destruct the Holder's value.
 *
 * @signature void clear(holder);
 *
 * @param[in,out] holder The Holder to clear.
 */

// ----------------------------------------------------------------------------
// Function value()
// ----------------------------------------------------------------------------

/*!
 * @fn Holder#value
 * @brief Return a reference to the value of the holder.
 *
 * @signature TReference value(holder);
 *
 * @param[in] holder The Holder to query for its reference.
 *
 * @return TReference The reference of the Holder's value.
 */

// ----------------------------------------------------------------------------
// Function moveValue()
// ----------------------------------------------------------------------------

/*!
 * @fn Holder#moveValue
 *
 * @headerfile <seqan/sequence.h>
 *
 * @brief Move a value of into a holder.
 *
 * @signature void moveValue(holder, value);
 *
 * @param[in,out] holder  The Holder to manipulate.
 * @param[in,out] value   The value to move into <tt>holder</tt>.
 */

// ----------------------------------------------------------------------------
// Function getValue()
// ----------------------------------------------------------------------------

/*!
 * @fn Holder#getValue
 * @brief Return the get-value of the holder.
 *
 * @signature TGetValue getValue(holder);
 *
 * @param[in] holder The Holder to query for its get-value type.
 *
 * @return TGetValue The get-value of the Holder.
 */

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
