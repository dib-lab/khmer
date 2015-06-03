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
// Tristate Holder Implementation.
// ==========================================================================

#ifndef SEQAN_BASIC_HOLDER_TRISTATE_H_
#define SEQAN_BASIC_HOLDER_TRISTATE_H_

// TODO(holtgrew): What about const holders?
// TODO(holtgrew): Are holders on pointers used anywhere?

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

#if SEQAN_ENABLE_POINTER_HOLDER
template <typename T> struct IsSimple;

// Used in functions section.
// TODO(holtgrew): This will go away, right?
template <typename TValue> inline size_t length(TValue const * me);
template <typename TValue> inline size_t length(TValue * me);
#endif  // #if SEQAN_ENABLE_POINTER_HOLDER

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

/*!
 * @class TristateHolder
 * @extends Holder
 * @headerfile <seqan/basic.h>
 * @brief Holder that can be empty, dependent, or owner.
 *
 * @signature template <typename TValue>
 *            class Holder<TValue, Tristate>;
 *
 * @tparam TValue Type of the managed object.
 *
 * A tristate holder <tt>A</tt> that holds an object <tt>B</tt> has one of the following states:
 *
 * <ul>
 *   <li>owner: <tt>A</tt> is the owner of <tt>B</tt>. If <tt>A</tt> is destroyed, <tt>B</tt> will be destroyed
 *       automatically.</li>
 *   <li>dependent: <tt>A</tt> depends on <tt>B</tt>. <tt>B</tt> should not be destroyed as long as <tt>A</tt> is
 *       used.</li>
 *   <li>empty: there is currently no object reference stored in the holder <tt>A</tt>.</li>
 * </ul>
 *
 * The state of the holder can be determined by empty and dependent.
 *
 * If a holder object is in owner state when destructed, the owned object is destructed as well.
 */

// TODO(holtgrew): This is broken for const TValue since we use Value<Holder>::Type below.

template <typename TValue>
struct Holder<TValue, Tristate>
{
    enum EHolderState
    {
        EMPTY = 0,
        OWNER = 1,
        DEPENDENT = 2
    };

    typedef typename Value<Holder>::Type THostValue;

    // ------------------------------------------------------------------------
    // Members
    // ------------------------------------------------------------------------

    THostValue * data_value;
    EHolderState data_state;

    // ------------------------------------------------------------------------
    // Constructors; Destructor
    // ------------------------------------------------------------------------

    Holder() : data_value(NULL), data_state(EMPTY)
    {
        SEQAN_CHECKPOINT;
    }

    Holder(Holder const & source_) : data_value(NULL), data_state(EMPTY)
    {
        SEQAN_CHECKPOINT;
        assign(*this, source_);
    }

    explicit
    Holder(THostValue & value_) : data_value(NULL), data_state(EMPTY)
    {
        SEQAN_CHECKPOINT;
        setValue(*this, value_);
    }

    explicit
    Holder(THostValue const & value_) : data_value(NULL), data_state(EMPTY)
    {
        SEQAN_CHECKPOINT;
        assignValue(*this, value_);
    }

    ~Holder()
    {
        SEQAN_CHECKPOINT;
        clear(*this);
    }

    // ------------------------------------------------------------------------
    // Assignment Operators;  Must be defined in class.
    // ------------------------------------------------------------------------

    inline Holder const &
    operator=(Holder const & source_)
    {
        SEQAN_CHECKPOINT;
        assign(*this, source_);
        return *this;
    }

    inline Holder const &
    operator=(THostValue const & value_)
    {
        SEQAN_CHECKPOINT;
        assignValue(*this, value_);
        return *this;
    }

    // ------------------------------------------------------------------------
    // Conversion Operators;  Must be defined in class.
    // ------------------------------------------------------------------------

    inline operator THostValue()
    {
        SEQAN_CHECKPOINT;
        return _dataValue(*this);
    }
};

// TODO(holtgrew): Why does this not reliably work without copying where it can?
template <typename TValue>
struct Holder<TValue const, Tristate>
{
    enum EHolderState
    {
        EMPTY = 0,
        OWNER = 1,
        DEPENDENT = 2
    };

    typedef typename Value<Holder>::Type THostValue;

    // ------------------------------------------------------------------------
    // Members
    // ------------------------------------------------------------------------

    THostValue * data_value;
    EHolderState data_state;

    // ------------------------------------------------------------------------
    // Constructors; Destructor
    // ------------------------------------------------------------------------

    Holder() : data_value(), data_state(EMPTY)
    {}

    Holder(Holder const & source_) : data_value(), data_state(EMPTY)
    {
        assign(*this, source_);
    }

    Holder(TValue & value_) : data_value(), data_state(EMPTY)
    {
        setValue(*this, value_);
    }

    Holder(TValue const & value_) : data_value(), data_state(EMPTY)
    {
        setValue(*this, value_);
    }

    ~Holder()
    {
        clear(*this);
    }

    // ------------------------------------------------------------------------
    // Assignment Operators;  Has to be defined in class.
    // ------------------------------------------------------------------------

    inline Holder &
    operator=(Holder const & source_)
    {
        SEQAN_CHECKPOINT;
        assign(*this, source_);
        return *this;
    }

    inline Holder &
    operator=(THostValue const & value_)
    {
        SEQAN_CHECKPOINT;
        assignValue(*this, value_);
        return *this;
    }

    // ------------------------------------------------------------------------
    // Conversion Operators;  Has to be defined in class.
    // ------------------------------------------------------------------------

    inline operator THostValue()
    {
        SEQAN_CHECKPOINT;
        return _dataValue(*this);
    }
};

#if SEQAN_ENABLE_POINTER_HOLDER
template <typename TValue>
struct Holder<TValue *, Tristate>
{

    enum EHolderState
    {
        EMPTY = 0,
        OWNER = 1,
        DEPENDENT = 2
    };

    typedef typename Value<Holder>::Type THostValue;

    // ------------------------------------------------------------------------
    // Members
    // ------------------------------------------------------------------------

    THostValue data_value;
    EHolderState data_state;

    // ------------------------------------------------------------------------
    // Constructor; Destructors
    // ------------------------------------------------------------------------

    Holder() : data_state(EMPTY)
    {
        SEQAN_CHECKPOINT;
    }

    Holder(Holder const & source_) : data_state(EMPTY)
    {
        SEQAN_CHECKPOINT;
        assign(*this, source_);
    }

    Holder(TValue * value_) : data_state(EMPTY)
    {
        SEQAN_CHECKPOINT;
        setValue(*this, value_);
    }

    ~Holder()
    {
        SEQAN_CHECKPOINT;
        clear(*this);
    }

    // ------------------------------------------------------------------------
    // Assignment Operators;  Have to be defined in class.
    // ------------------------------------------------------------------------

    inline Holder const &
    operator = (Holder const & source_)
    {
        SEQAN_CHECKPOINT;
        assign(*this, source_);
        return *this;
    }

    inline Holder const &
    operator = (THostValue value_)
    {
        SEQAN_CHECKPOINT;
        setValue(*this, value_);
        return *this;
    }

    // ------------------------------------------------------------------------
    // Conversion Operators;  Have to be defined in class.
    // ------------------------------------------------------------------------

    inline operator THostValue()
    {
        SEQAN_CHECKPOINT;
        return _dataValue(*this);
    }
};

template <typename TValue>
struct Holder<TValue * const, Tristate>
{
    enum EHolderState
    {
        EMPTY = 0,
        OWNER = 1,
        DEPENDENT = 2
    };

    typedef typename Value<Holder>::Type THostValue;

    // ------------------------------------------------------------------------
    // Members
    // ------------------------------------------------------------------------

    THostValue data_value;
    EHolderState data_state;

    // ------------------------------------------------------------------------
    // Constructors; Destructor
    // ------------------------------------------------------------------------

    Holder() : data_state(EMPTY)
    {
        SEQAN_CHECKPOINT;
    }

    Holder(Holder const & source_) : data_state(EMPTY)
    {
        SEQAN_CHECKPOINT;
        assign(*this, source_);
    }

    Holder(TValue * value_) : data_state(EMPTY)
    {
        SEQAN_CHECKPOINT;
        setValue(*this, value_);
    }

    ~Holder()
    {
        SEQAN_CHECKPOINT;
        clear(*this);
    }

    // ------------------------------------------------------------------------
    // Assignment Operators;  Have to be defined in class.
    // ------------------------------------------------------------------------

    inline Holder &
    operator=(Holder const & source_)
    {
        SEQAN_CHECKPOINT;
        assign(*this, source_);
        return *this;
    }

    inline Holder &
    operator=(THostValue value_)
    {
        SEQAN_CHECKPOINT;
        setValue(*this, value_);
        return *this;
    }

    // ------------------------------------------------------------------------
    // Conversion Operators;  Have to be defined in class.
    // ------------------------------------------------------------------------

    inline operator THostValue()
    {
        SEQAN_CHECKPOINT;
        return _dataValue(*this);
    }
};
#endif  // #if SEQAN_ENABLE_POINTER_HOLDER

// ============================================================================
// Metafunctions
// ============================================================================

// ----------------------------------------------------------------------------
// Metafunction Spec
// ----------------------------------------------------------------------------

// Spec has to be given for all default specializations.

template <typename TValue>
struct Spec<Holder<TValue, Tristate> >
{
    typedef Tristate Type;
};

template <typename TValue>
struct Spec<Holder<TValue, Tristate> const>
{
    typedef Tristate Type;
};

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function _dataValue()
// ----------------------------------------------------------------------------

template <typename TValue>
inline typename Reference<Holder<TValue, Tristate> >::Type
_dataValue(Holder<TValue, Tristate> & me)
{
    return * me.data_value;
}

template <typename TValue>
inline typename Reference<Holder<TValue, Tristate> const>::Type
_dataValue(Holder<TValue, Tristate> const & me)
{
    return * me.data_value;
}

#if SEQAN_ENABLE_POINTER_HOLDER
template <typename TValue>
inline typename Reference<Holder<TValue *, Tristate> >::Type
_dataValue(Holder<TValue *, Tristate> & me)
{
    return me.data_value;
}

template <typename TValue>
inline typename Reference<Holder<TValue *, Tristate> const>::Type
_dataValue(Holder<TValue *, Tristate> const & me)
{
    return me.data_value;
}

template <typename TValue>
inline typename Reference<Holder<TValue * const, Tristate> >::Type
_dataValue(Holder<TValue * const, Tristate> & me)
{
    return me.data_value;
}

template <typename TValue>
inline typename Reference<Holder<TValue * const, Tristate> const>::Type
_dataValue(Holder<TValue * const, Tristate> const & me)
{
    return me.data_value;
}
#endif  // #if SEQAN_ENABLE_POINTER_HOLDER

// ----------------------------------------------------------------------------
// Function empty()
// ----------------------------------------------------------------------------

template <typename TValue>
inline bool
empty(Holder<TValue, Tristate> const & me)
{
SEQAN_CHECKPOINT
    return (me.data_state == Holder<TValue, Tristate>::EMPTY);
}

// ----------------------------------------------------------------------------
// Function dependent()
// ----------------------------------------------------------------------------

template <typename TValue>
inline bool
dependent(Holder<TValue, Tristate> const & me)
{
    SEQAN_CHECKPOINT;
    return me.data_state == Holder<TValue, Tristate>::DEPENDENT;
}

/// ----------------------------------------------------------------------------
// Function clear()
// ----------------------------------------------------------------------------

template <typename THolder, typename TValue>
inline void
_holderDeallocate(THolder & me, TValue const & data)
{
    valueDestruct(& data);
    deallocate(me, & data, 1);
}

template <typename THolder, typename TValue>
inline void
_holderDeallocate(THolder & me, TValue * data, True)                // is a pointer to an *array* of objects
{
    size_t len = length(data)+1;
    arrayDestruct(data, data+len);
    deallocate(me, data, len);
}

template <typename THolder, typename TValue>
inline void
_holderDeallocate(THolder & /*me*/, TValue * /*data*/, False)               // is a pointer to *one* object
{
}

template <typename THolder, typename TValue>
inline void
_holderDeallocate(THolder & me, TValue * data)
{
    return _holderDeallocate(me, data, IsSimple<TValue>());         // try to distinguish between a pointer to one/array of object(s)
}

///the state of this object is set to 'empty'.

template <typename TValue>
inline void
clear(Holder<TValue, Tristate> & me)
{
    switch (me.data_state) {
        case Holder<TValue, Tristate>::EMPTY:
            break;
        case Holder<TValue, Tristate>::DEPENDENT:
            SEQAN_CHECKPOINT;
            me.data_state = Holder<TValue, Tristate>::EMPTY;
            break;
        default:  // case Holder<TValue, TSpec>::OWNER
            SEQAN_CHECKPOINT;
            _holderDeallocate(me, _dataValue(me));
            me.data_state = Holder<TValue, Tristate>::EMPTY;
            break;
    }
}

// ----------------------------------------------------------------------------
// Function create()
// ----------------------------------------------------------------------------

template <typename THolder, typename TValue>
inline typename Value<THolder, 0>::Type *
_holderAllocateObject(THolder & me, TValue const & data)
{
    typename Value<THolder>::Type * ret;
    allocate(me, ret, 1);
    valueConstruct(ret, data);
    return ret;
}

#if SEQAN_ENABLE_POINTER_HOLDER
template <typename THolder, typename TValue>
inline typename Value<THolder, 0>::Type
_holderAllocatePointer(THolder & me, TValue * data, True)           // is a pointer to an *array* of objects
{
    typename Value<THolder>::Type ret;
    size_t len = length(data)+1;
    allocate(me, ret, len);
    arrayConstructCopy(data, data + len, ret);
    return ret;
}

template <typename THolder, typename TValue>
inline typename Value<THolder, 0>::Type
_holderAllocatePointer(THolder & /*me*/, TValue * data, False)          // is a pointer to *one* object
{
    return data;
}

template <typename THolder, typename TValue>
inline typename Value<THolder, 0>::Type
_holderAllocatePointer(THolder & me, TValue * data)
{
    return _holderAllocatePointer(me, data, IsSimple<TValue>());    // try to distinguish between a pointer to one/array of object(s)
}
#endif  // #if SEQAN_ENABLE_POINTER_HOLDER

template <typename TValue>
inline void
create(Holder<TValue, Tristate> & me)
{
    typedef Holder<TValue, Tristate> THolder;

    switch (me.data_state)
    {
        case Holder<TValue, Tristate>::EMPTY:
            SEQAN_CHECKPOINT;
            allocate(me, me.data_value, 1);
            valueConstruct(me.data_value);
            me.data_state = THolder::OWNER;
            break;

        case THolder::DEPENDENT:
            SEQAN_CHECKPOINT;
            create(me, _dataValue(me));
            break;
        default:;
    }
}

#if SEQAN_ENABLE_POINTER_HOLDER
template <typename TValue>
inline void
create(Holder<TValue *, Tristate> & me)
{
    typedef Holder<TValue *, Tristate> THolder;

    switch (me.data_state) {
        case Holder<TValue *, Tristate>::EMPTY:
            SEQAN_CHECKPOINT;
            valueConstruct(& me.data_value);
            me.data_state = THolder::OWNER;
            break;
        case THolder::DEPENDENT:
            SEQAN_CHECKPOINT;
            create(me, _dataValue(me));
            break;
        default:;
    }
}

template <typename TValue>
inline void
create(Holder<TValue * const, Tristate> & me)
{
    typedef Holder<TValue *, Tristate> THolder;

    switch (me.data_state) {
        case Holder<TValue *, Tristate>::EMPTY:
            SEQAN_CHECKPOINT;
            valueConstruct(& me.data_value);
            me.data_state = THolder::OWNER;
            break;
        case THolder::DEPENDENT:
            SEQAN_CHECKPOINT;
            create(me, _dataValue(me));
            break;
        default:;
    }
}
#endif  // #if SEQAN_ENABLE_POINTER_HOLDER

template <typename TValue, typename TValue2>
inline void
create(Holder<TValue, Tristate> & me,
       TValue2 & value_)
{
    SEQAN_CHECKPOINT;

    if (me.data_state == Holder<TValue, Tristate>::OWNER) {
        assign(_dataValue(me), value_);
        return;
    }

    clear(me);
    me.data_value = _holderAllocateObject(me, value_);
    me.data_state = Holder<TValue, Tristate>::OWNER;
}

template <typename TValue, typename TValue2>
inline void
create(Holder<TValue const, Tristate> & me,
       TValue2 & value_)
{
    SEQAN_CHECKPOINT;

    clear(me);
    me.data_value = _holderAllocateObject(me, value_);
    me.data_state = Holder<TValue const, Tristate>::OWNER;
}

#if SEQAN_ENABLE_POINTER_HOLDER
template <typename TValue, typename TValue2>
inline void
create(Holder<TValue *, Tristate> & me,
       TValue2 & value_)
{
SEQAN_CHECKPOINT

    clear(me);
    me.data_value = _holderAllocatePointer(me, value_);
    me.data_state = Holder<TValue *, Tristate>::OWNER;
}

template <typename TValue, typename TValue2>
inline void
create(Holder<TValue * const, Tristate> & me,
       TValue2 & value_)
{
SEQAN_CHECKPOINT

    clear(me);
    me.data_value = _holderAllocatePointer(me, value_);
    me.data_state = Holder<TValue *, Tristate>::OWNER;
}
#endif  // #if SEQAN_ENABLE_POINTER_HOLDER

template <typename TValue, typename TValue2>
inline void
create(Holder<TValue const, Tristate> & me,
       TValue2 & value_,
       Move const &)
{
    SEQAN_CHECKPOINT;
    // TODO(holtgrew): Real implementation once HasMoveConstructor metafunction is in place.
    me.data_value = value_;
}

template <typename TValue, typename TValue2>
inline void
create(Holder<TValue, Tristate> & me,
       TValue2 & value_,
       Move const &)
{
    SEQAN_CHECKPOINT;
    // TODO(holtgrew): Real implementation once HasMoveConstructor metafunction is in place.
    me.data_value = value_;
}

//////////////////////////////////////////////////////////////////////////////
template <typename TValue>
inline void
detach(Holder<TValue, Tristate> & me)
{
    SEQAN_CHECKPOINT;
    create(me);
}

// ----------------------------------------------------------------------------
// Function setValue()
// ----------------------------------------------------------------------------

template <typename TValue>
inline void
setValue(Holder<TValue, Tristate> & me,
         TValue & value_)
{
    SEQAN_CHECKPOINT;
    clear(me);
    me.data_value = & value_;
    me.data_state = Holder<TValue, Tristate>::DEPENDENT;
}

template <typename TValue>
inline void
setValue(Holder<TValue const, Tristate> & me,
         TValue & value_)
{
    SEQAN_CHECKPOINT;
    clear(me);
    me.data_value = & value_;
    me.data_state = Holder<TValue const, Tristate>::DEPENDENT;
}

// ----------------------------------------------------------------------------
// Function setValue()
// ----------------------------------------------------------------------------

#if SEQAN_ENABLE_POINTER_HOLDER
template <typename TValue>
inline void
setValue(Holder<TValue *, Tristate> & me,
         TValue * & value_)
{
    SEQAN_CHECKPOINT;
    clear(me);
    me.data_value = value_;
    me.data_state = Holder<TValue *, Tristate>::DEPENDENT;
}

template <typename TValue>
inline void
setValue(Holder<TValue *, Tristate> & me,
         TValue * const & value_)
{
    SEQAN_CHECKPOINT;
    clear(me);
    me.data_value = value_;
    me.data_state = Holder<TValue *, Tristate>::DEPENDENT;
}

template <typename TValue>
inline void
setValue(Holder<TValue * const, Tristate> & me,
         TValue * & value_)
{
    SEQAN_CHECKPOINT;
    clear(me);
    me.data_value = value_;
    me.data_state = Holder<TValue *, Tristate>::DEPENDENT;
}

template <typename TValue>
inline void
setValue(Holder<TValue * const, Tristate> & me,
         TValue * const & value_)
{
    SEQAN_CHECKPOINT;
    clear(me);
    me.data_value = value_;
    me.data_state = Holder<TValue *, Tristate>::DEPENDENT;
}

template <typename TValue, size_t I>
inline void
setValue(Holder<TValue *, Tristate> & me,
         TValue (& value_)[I])
{
    SEQAN_CHECKPOINT;
    clear(me);
    me.data_value = value_;
    me.data_state = Holder<TValue *, Tristate>::DEPENDENT;
}

template <typename TValue, size_t I>
inline void
setValue(Holder<TValue *, Tristate> & me,
         TValue const (& value_)[I])
{
    SEQAN_CHECKPOINT;
    clear(me);
    me.data_value = value_;
    me.data_state = Holder<TValue *, Tristate>::DEPENDENT;
}

template <typename TValue, size_t I>
inline void
setValue(Holder<TValue * const, Tristate> & me,
         TValue (& value_)[I])
{
    SEQAN_CHECKPOINT;
    clear(me);
    me.data_value = value_;
    me.data_state = Holder<TValue *, Tristate>::DEPENDENT;
}

template <typename TValue, size_t I>
inline void
setValue(Holder<TValue * const, Tristate> & me,
         TValue const (& value_)[I])
{
    SEQAN_CHECKPOINT;
    clear(me);
    me.data_value = value_;
    me.data_state = Holder<TValue *, Tristate>::DEPENDENT;
}
#endif  // #if SEQAN_ENABLE_POINTER_HOLDER

template <typename TValue, typename TValue2>
inline void
setValue(Holder<TValue, Tristate> & me,
         TValue2 & value_)
{
    SEQAN_CHECKPOINT;
    set(value(me), value_);
}

template <typename TValue, typename TValue2>
inline void
setValue(Holder<TValue, Tristate> & me,
         TValue2 const & value_)
{
    SEQAN_CHECKPOINT;
    set(value(me), value_);
}

// ----------------------------------------------------------------------------
// Function swap()
// ----------------------------------------------------------------------------

template <typename TValue>
void swap(Holder<TValue, Tristate> & lhs, Holder<TValue, Tristate> & rhs)
{
    std::swap(lhs.data_value, rhs.data_value);
    std::swap(lhs.data_state, rhs.data_state);
}

// ----------------------------------------------------------------------------
// Function value()
// ----------------------------------------------------------------------------

template <typename TValue>
inline typename Reference<Holder<TValue, Tristate> >::Type
value(Holder<TValue, Tristate> & me)
{
    if (empty(me))
        create(me);

    return _dataValue(me);
}

template <typename TValue>
inline typename Reference<Holder<TValue, Tristate> const>::Type
value(Holder<TValue, Tristate> const & me)
{
    SEQAN_CHECKPOINT;
    SEQAN_ASSERT_NOT(empty(me));

    return _dataValue(me);
}

// ----------------------------------------------------------------------------
// Function assignValue()
// ----------------------------------------------------------------------------

template <typename TValue, typename TSource>
inline void
assignValue(Holder<TValue, Tristate> & me,
            TSource const & value_)
{
    if (empty(me)) {
        create(me, value_);
    } else {
        assign(_dataValue(me), value_);
    }
}

// ----------------------------------------------------------------------------
// Function moveValue()
// ----------------------------------------------------------------------------

template <typename TValue, typename TSource>
inline void
moveValue(Holder<TValue, Tristate> & me,
          TSource const & value_)
{
    SEQAN_CHECKPOINT;
    if (empty(me)) {
        create(me, value_);
    } else {
        move(_dataValue(me), value_);
    }
}

// ----------------------------------------------------------------------------
// Function assign()
// ----------------------------------------------------------------------------

template <typename TValue>
inline void
assign(Holder<TValue, Tristate> & target_,
       Holder<TValue, Tristate> const & source_)
{
    SEQAN_CHECKPOINT;
    switch(source_.data_state) {
        case Holder<TValue, Tristate>::EMPTY:
            clear(target_);
            break;

        case Holder<TValue, Tristate>::OWNER:
            assignValue(target_, value(source_));
            break;

        default:  // case Holder<TValue, Tristate>::DEPENDENT
            setValue(target_, value(source_));
            break;
    }
}

template <typename TValue>
inline void
assign(Holder<TValue const, Tristate> & target_,
       Holder<TValue const, Tristate> const & source_)
{
    SEQAN_CHECKPOINT;
    switch(source_.data_state) {
        case Holder<TValue, Tristate>::EMPTY:
            clear(target_);
            break;

        case Holder<TValue, Tristate>::OWNER:
            create(target_, value(source_));
            break;

        default:  // case Holder<TValue, Tristate>::DEPENDENT
            setValue(target_, value(source_));
            break;
    }
}

}  // namespace seqan

#endif  // #ifndef SEQAN_BASIC_HOLDER_TRISTATE_H_
