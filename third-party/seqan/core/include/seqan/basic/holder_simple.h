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
// Simple Holder specialization.
// ==========================================================================

#ifndef SEQAN_BASIC_HOLDER_SIMPLE_H_
#define SEQAN_BASIC_HOLDER_SIMPLE_H_

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

#ifdef PLATFORM_WINDOWS_VS
// Disable warning C4521 locally (multiple copy constructors).
#pragma warning( push )
#pragma warning( disable: 4521 )
#endif  // PLATFORM_WINDOWS_VS

/**
.Spec.Simple Holder
..cat:Holders
..summary:Simple copying wrapper without any additional state.
..signature:Holder<TValue, Simple>
..param.TValue:Type of the managed object.
..general:Class.Holder
..remarks.text:This holder stores a copy of the value.
..include:seqan/basic.h
 */

template <typename TValue>
struct Holder<TValue, Simple>
{
    typedef typename Value<Holder>::Type THolderValue;
    typedef typename Parameter_<THolderValue>::Type THolderParameter;

    // ------------------------------------------------------------------------
    // Members
    // ------------------------------------------------------------------------

    mutable typename RemoveConst_<THolderValue>::Type data_value;

    // ------------------------------------------------------------------------
    // Constructors
    // ------------------------------------------------------------------------

    Holder() {}

    Holder(Holder & source_) : data_value(source_.data_value)
    {
        SEQAN_CHECKPOINT;
    }

    Holder(Holder const & source_) : data_value(source_.data_value)
    {
        SEQAN_CHECKPOINT;
    }

    template <typename TSource>
    explicit
    Holder(TSource & value_) : data_value(value_)
    {
        SEQAN_CHECKPOINT;
    }

    template <typename TSource>
    explicit
    Holder(TSource const & value_) : data_value(value_)
    {
        SEQAN_CHECKPOINT;
    }

    // ------------------------------------------------------------------------
    // Assignment Operators;  Have to be defined in class.
    // ------------------------------------------------------------------------

    Holder &
    operator=(Holder const & source_)
    {
        SEQAN_CHECKPOINT;
        data_value = source_.data_value;
        return *this;
    }

    Holder &
    operator=(THolderValue const & value_)
    {
        SEQAN_CHECKPOINT;
        data_value = value_;
        return *this;
    }

    // ------------------------------------------------------------------------
    // Conversion Operators;  Have to be defined in class.
    // ------------------------------------------------------------------------

    operator THolderParameter()
    {
        SEQAN_CHECKPOINT;
        return *data_value;
    }
};

#ifdef PLATFORM_WINDOWS_VS
// Set old warning C4521 state again (multiple copy constructors).
#pragma warning( pop )
#endif  // PLATFORM_WINDOWS_VS

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function empty()
// ----------------------------------------------------------------------------

template <typename TValue>
inline bool
empty(Holder<TValue, Simple> const & /*me*/)
{
    SEQAN_CHECKPOINT;
    return false;
}

// ----------------------------------------------------------------------------
// Function dependent()
// ----------------------------------------------------------------------------

template <typename TValue>
inline bool
dependent(Holder<TValue, Simple> const & /*me*/)
{
    SEQAN_CHECKPOINT;
    return false;
}

// ----------------------------------------------------------------------------
// Function clear()
// ----------------------------------------------------------------------------

template <typename TValue>
inline void
clear(Holder<TValue, Simple> & /*me*/)
{
    SEQAN_CHECKPOINT;
}

// ----------------------------------------------------------------------------
// Function create()
// ----------------------------------------------------------------------------

template <typename TValue>
inline void
create(Holder<TValue, Simple> & /*me*/)
{
    SEQAN_CHECKPOINT;
    // TODO(holtgrew): Should be create(me.data_value), right?
}

template <typename TValue>
inline void
create(Holder<TValue, Simple> & me,
       TValue const & value_)
{
    SEQAN_CHECKPOINT;
    me.data_value = value_;
}

template <typename TValue>
inline void
create(Holder<TValue, Simple> & me,
       TValue const & value_,
       Move const &)
{
    SEQAN_CHECKPOINT;
    // TODO(holtgrew): Real implementation once HasMoveConstructor metafunction is in place.
    me.data_value = value_;
}

// ----------------------------------------------------------------------------
// Function detach()
// ----------------------------------------------------------------------------

template <typename TValue>
inline void
detach(Holder<TValue, Simple> & /*me*/)
{
    SEQAN_CHECKPOINT;
}

// ----------------------------------------------------------------------------
// Function setValue()
// ----------------------------------------------------------------------------

template <typename TValue>
inline void
setValue(Holder<TValue, Simple> & me,
         TValue const & value_)
{
    SEQAN_CHECKPOINT;
    set(me.data_value, value_);
}

// ----------------------------------------------------------------------------
// Function value()
// ----------------------------------------------------------------------------

template <typename TValue>
inline typename Reference<Holder<TValue, Simple> >::Type
value(Holder<TValue, Simple> & me)
{
    SEQAN_CHECKPOINT;
    return me.data_value;
}

template <typename TValue>
inline typename Reference<Holder<TValue, Simple> const>::Type
value(Holder<TValue, Simple> const & me)
{
    SEQAN_CHECKPOINT;
    return me.data_value;
}

// ----------------------------------------------------------------------------
// Function getValue()
// ----------------------------------------------------------------------------

template <typename TValue>
inline typename GetValue<Holder<TValue, Simple> >::Type
getValue(Holder<TValue, Simple> & me)
{
    SEQAN_CHECKPOINT;
    return me.data_value;
}

template <typename TValue>
inline typename GetValue<Holder<TValue, Simple> const>::Type
getValue(Holder<TValue, Simple> const & me)
{
    SEQAN_CHECKPOINT;
    return me.data_value;
}

// ----------------------------------------------------------------------------
// Function assignValue()
// ----------------------------------------------------------------------------

template <typename TValue, typename TSource>
inline void
assignValue(Holder<TValue, Simple> & me,
            TSource const & value_)
{
    SEQAN_CHECKPOINT;
    assign(me.data_value, value_);
}

// ----------------------------------------------------------------------------
// Function moveValue()
// ----------------------------------------------------------------------------

template <typename TValue, typename TSource>
inline void
moveValue(Holder<TValue, Simple> & me,
          TSource const & value_)
{
    SEQAN_CHECKPOINT;
    move(me.data_value, value_);
}

// ----------------------------------------------------------------------------
// Function assign()
// ----------------------------------------------------------------------------

template <typename TValue>
inline void
assign(Holder<TValue, Simple> & target_,
       Holder<TValue, Simple> const & source_)
{
    SEQAN_CHECKPOINT;
    assignValue(target_, source_);
}

}  // namespace seqan

#endif  // #ifndef SEQAN_BASIC_HOLDER_SIMPLE_H_
