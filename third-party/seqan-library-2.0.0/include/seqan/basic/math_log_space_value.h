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
// Author: Tobias Rausch <rausch@embl.de>
// Author: David Weese <david.weese@fu-berlin.de>
// ==========================================================================
// Type to store double values in log space.  Multiplication with normal
// number corresponds to addtion there, such they are more stable for large
// values.
// ==========================================================================

// TODO(holtgrew): Rename to "LogSpaceValue"?

#ifndef SEQAN_INCLUDE_SEQAN_BASIC_MATH_LOG_SPACE_VALUE_H_
#define SEQAN_INCLUDE_SEQAN_BASIC_MATH_LOG_SPACE_VALUE_H_

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

/*!
 * @class LogProb
 * @headerfile <seqan/basic.h>
 * @brief Value type for computation in log-space.
 *
 * @signature template <[typename T]>
 *            class LogProb;
 *
 * @tparam T Floating number type to use as the basic, defaults to <tt>double</tt>.
 *
 * @section Remarks
 *
 * Internally, the logarithms of the original values are stored.  This is numerically more stable for multiplications
 * and large numbers.
 *
 * This type can be used like an ordinary <tt>double</tt> or <tt>float</tt> value.
 *
 * @section Examples
 *
 * @code{.cpp}
 * LogProb<double> x = 10;
 * x *= 3;
 * x += 4;
 * int y = x;
 * @endcode
 */

template<typename TValue = double, typename TSpec = Default>
class LogProb;


// ----------------------------------------------------------------------------
// Function log()
// ----------------------------------------------------------------------------

template <typename TValue, typename TSpec>
inline TValue
log(LogProb<TValue, TSpec> const & value)
{
    return value.data_value;
}

template <typename TValue>
inline double
log(TValue const & value)
{
    return std::log((double)value);
}


template<typename TValue, typename TSpec>
class LogProb
{
  public:
    TValue data_value;

    // ------------------------------------------------------------------------
    // Constructors
    // ------------------------------------------------------------------------

    LogProb() : data_value(std::log(0.0)) {}

    template <typename TValue2>
    LogProb(TValue2 const & _other)
    {
        data_value = log(_other);
    }

    // ------------------------------------------------------------------------
    // Type conversion operators;  Have to be defined in class.
    // ------------------------------------------------------------------------

    template <typename TResult>
    inline TResult _convert(False) const
    {
        return (TResult)std::exp(data_value);
    }

    template <typename TResult>
    inline TResult _convert(True) const
    {
        return (TResult)round(std::exp(data_value));
    }

    template <typename TValue2>
    inline
    operator TValue2() const
    {
        return _convert<TValue2>(typename IsInteger<TValue2>::Type());
    }
/*
    inline
    operator float() const
    {
        return (float) std::exp(data_value);
    }

    inline
    operator double() const
    {
        return (double) std::exp(data_value);
    }
*/
    // ------------------------------------------------------------------------
    // Function operator=();  Has to be defined in class.
    // ------------------------------------------------------------------------

    template<typename TRhs>
    inline
    LogProb &
    operator=(TRhs const& rhs)
    {
        data_value = std::log(rhs);
        return *this;
    }

    template<typename TValue2, typename TSpec2>
    inline
    LogProb &
    operator=(LogProb<TValue2, TSpec2> const & rhs)
    {
        data_value = rhs.data_value;
        return *this;
    }
};

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================


// ----------------------------------------------------------------------------
// Function operator*=()
// ----------------------------------------------------------------------------

template <typename TValue, typename TSpec, typename TRhs>
inline
LogProb<TValue, TSpec> &
operator*=(LogProb<TValue, TSpec> & lhs,
           TRhs const & rhs)
{
    lhs.data_value += log(rhs);
    return lhs;
}

// ----------------------------------------------------------------------------
// Function operator*()
// ----------------------------------------------------------------------------

template <typename TValue, typename TSpec, typename TRhs>
inline
LogProb<TValue, TSpec>
operator*(LogProb<TValue, TSpec> const & lhs,
          TRhs const & rhs)
{
    LogProb<TValue, TSpec> result = lhs;
    result *= LogProb<TValue, TSpec>(rhs);
    return result;
}

template <typename TValue, typename TSpec, typename TValue2, typename TSpec2>
inline
LogProb<TValue, TSpec>
operator*(LogProb<TValue, TSpec> const & lhs,
          LogProb<TValue2, TSpec2> const & rhs)
{
    LogProb<TValue, TSpec> result = lhs;
    result *= rhs;
    return result;
}

// ----------------------------------------------------------------------------
// Function operator/=()
// ----------------------------------------------------------------------------

template <typename TValue, typename TSpec, typename TRhs>
inline
LogProb<TValue, TSpec> &
operator/=(LogProb<TValue, TSpec> & lhs,
           TRhs const & rhs)
{
    lhs.data_value -= log(rhs);
    return lhs;
}

// ----------------------------------------------------------------------------
// Function operator/()
// ----------------------------------------------------------------------------

template <typename TValue, typename TSpec, typename TRhs>
inline
LogProb<TValue, TSpec>
operator/(LogProb<TValue, TSpec> const & lhs,
          TRhs const & rhs)
{
    LogProb<TValue, TSpec> result = lhs;
    result /= LogProb<TValue, TSpec>(rhs);
    return result;
}

template <typename TValue, typename TSpec, typename TValue2, typename TSpec2>
inline
LogProb<TValue, TSpec>
operator/(LogProb<TValue, TSpec> const & lhs,
          LogProb<TValue2, TSpec2> const & rhs)
{
    LogProb<TValue, TSpec> result = lhs;
    result /= rhs;
    return result;
}

// ----------------------------------------------------------------------------
// Function operator+=()
// ----------------------------------------------------------------------------

template <typename TValue, typename TSpec, typename TRhs>
inline
LogProb<TValue, TSpec> &
operator+=(LogProb<TValue, TSpec> & lhs,
           TRhs const & rhs) {
    lhs.data_value = std::log(std::exp(lhs.data_value) + rhs);
    return lhs;
}

template <typename TValue, typename TSpec, typename TValue2, typename TSpec2>
inline
LogProb<TValue, TSpec> &
operator+=(LogProb<TValue, TSpec> & lhs,
           LogProb<TValue2, TSpec2> const & rhs)
{
    if (lhs.data_value > rhs.data_value) {
        if ((rhs.data_value == std::log(0.0)) || (lhs.data_value - rhs.data_value > 100))
            return lhs;
        lhs.data_value = lhs.data_value + std::log(1 + std::exp(rhs.data_value - lhs.data_value));
    } else {
        if ((lhs.data_value == std::log(0.0)) || (rhs.data_value - lhs.data_value > 100)) {
            lhs.data_value = rhs.data_value;
            return lhs;
        }
        lhs.data_value = rhs.data_value + std::log(1 + std::exp(lhs.data_value - rhs.data_value));
    }
    return lhs;
}

// ----------------------------------------------------------------------------
// Function operator+()
// ----------------------------------------------------------------------------

template <typename TValue, typename TSpec, typename TRhs>
inline
LogProb<TValue, TSpec>
operator+(LogProb<TValue, TSpec> const & lhs,
          TRhs const & rhs) {
    LogProb<TValue, TSpec> result = lhs;
    result += LogProb<TValue, TSpec>(rhs);
    return result;
}

template <typename TValue, typename TSpec, typename TValue2, typename TSpec2>
inline
LogProb<TValue, TSpec>
operator+(LogProb<TValue, TSpec> const & lhs,
          LogProb<TValue2, TSpec2> const & rhs)
{
    LogProb<TValue, TSpec> result = lhs;
    result += rhs;
    return result;
}

// ----------------------------------------------------------------------------
// Function operator-=()
// ----------------------------------------------------------------------------

template <typename TValue, typename TSpec, typename TRhs>
inline
LogProb<TValue, TSpec> &
operator-=(LogProb<TValue, TSpec> & lhs,
           TRhs const & rhs)
{
    lhs.data_value = std::log(std::exp(lhs.data_value) - rhs);
    return lhs;
}

template <typename TValue, typename TSpec, typename TValue2, typename TSpec2>
inline
LogProb<TValue, TSpec> &
operator-=(LogProb<TValue, TSpec> & lhs,
           LogProb<TValue2, TSpec2> const & rhs)
{
    if (lhs.data_value > rhs.data_value) {
        if ((rhs.data_value == std::log(0.0)) || (lhs.data_value - rhs.data_value > 100))
            return lhs;
        lhs.data_value = lhs.data_value + std::log(1 - std::exp(rhs.data_value - lhs.data_value));
    } else {
        if ((lhs.data_value == std::log(0.0)) || (rhs.data_value - lhs.data_value > 100)) {
            lhs.data_value = rhs.data_value;
            return lhs;
        }
        lhs.data_value = rhs.data_value + std::log(1 - std::exp(lhs.data_value - rhs.data_value));
    }
    return lhs;
}

// ----------------------------------------------------------------------------
// Function operator-()
// ----------------------------------------------------------------------------

template <typename TValue, typename TSpec, typename TRhs>
inline
LogProb<TValue, TSpec>
operator-(LogProb<TValue, TSpec> const & lhs,
          TRhs const & rhs) {
    LogProb<TValue, TSpec> result = lhs;
    result -= LogProb<TValue, TSpec>(rhs);
    return result;
}

template <typename TValue, typename TSpec, typename TValue2, typename TSpec2>
inline
LogProb<TValue, TSpec>
operator-(LogProb<TValue, TSpec> const & lhs,
          LogProb<TValue2, TSpec2> const & rhs)
{
    LogProb<TValue, TSpec> result = lhs;
    result -= rhs;
    return result;
}

// ----------------------------------------------------------------------------
// Function operator==()
// ----------------------------------------------------------------------------

template <typename TValue, typename TSpec, typename TRhs>
inline
bool
operator==(LogProb<TValue, TSpec> const & lhs,
           TRhs const & rhs)
{
    return lhs.data_value == log(rhs);
}

// ----------------------------------------------------------------------------
// Function operator!=()
// ----------------------------------------------------------------------------

template <typename TValue, typename TSpec, typename TRhs>
inline
bool
operator!=(LogProb<TValue, TSpec> const & lhs,
           TRhs const & rhs)
{
    return lhs.data_value != log(rhs);
}

// ----------------------------------------------------------------------------
// Function operator<()
// ----------------------------------------------------------------------------

template <typename TValue, typename TSpec, typename TRhs>
inline
bool
operator<(LogProb<TValue, TSpec> const & lhs,
          TRhs const & rhs)
{
    return lhs.data_value < log(rhs);
}

// ----------------------------------------------------------------------------
// Function operator>()
// ----------------------------------------------------------------------------

template <typename TValue, typename TSpec, typename TRhs>
inline
bool
operator>(LogProb<TValue, TSpec> const & lhs,
          TRhs const & rhs)
{
    return lhs.data_value > log(rhs);
}

// ----------------------------------------------------------------------------
// Function operator<=()
// ----------------------------------------------------------------------------

template <typename TValue, typename TSpec, typename TRhs>
inline
bool
operator<=(LogProb<TValue, TSpec> const & lhs,
           TRhs const & rhs)
{
    return lhs.data_value <= log(rhs);
}

// ----------------------------------------------------------------------------
// Function operator>=()
// ----------------------------------------------------------------------------

template <typename TValue, typename TSpec, typename TRhs>
inline
bool
operator>=(LogProb<TValue, TSpec> const & lhs,
           TRhs const & rhs)
{
    return lhs.data_value >= log(rhs);
}

// ----------------------------------------------------------------------------
// Function operator<<
// ----------------------------------------------------------------------------

// Stream output for LogProb values.

template<typename TStream, typename TValue, typename TSpec>
inline
TStream &
operator<<(TStream & stream, LogProb<TValue, TSpec> const & rhs)
{
    return stream << std::exp(rhs.data_value);
}

}  // namespace seqan

#endif  // #ifndef SEQAN_INCLUDE_SEQAN_BASIC_MATH_LOG_SPACE_VALUE_H_
