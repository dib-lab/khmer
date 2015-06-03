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
// Author: David Weese <david.weese@fu-berlin.de>
// ==========================================================================
// Triple base class.
// ==========================================================================

// TODO(holtgrew): What about move construction? Useful for pairs of strings and such. Tricky to implement since ints have no move constructor, for example.

#ifndef SEQAN_INCLUDE_SEQAN_BASIC_TRIPLE_BASE_H_
#define SEQAN_INCLUDE_SEQAN_BASIC_TRIPLE_BASE_H_

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

/*!
 * @class Triple
 * @implements ComparableConcept
 * @implements LessThanComparableConcept
 * @headerfile <seqan/basic.h>
 * @brief Store three arbitrary object.
 *
 * @signature template <typename T1, typename T3, typename T3[, typename TSpec]>
 *            class Triple;
 *
 * @tparam T1 Type of first object.
 * @tparam T2 Type of second object.
 * @tparam T3 Type of third object.
 * @tparam TSpec Tag for specialization (Default: <tt>void</tt>).
 */

/*!
 * @fn Triple::Triple
 * @brief Default and copy construction and construction with three objects.
 *
 * @signature Triple::Triple()
 * @signature Triple::Triple(other)
 * @signature Triple::Triple(x1, x2, x3)
 *
 * @param[in] other Other Triple object to copy from.
 * @param[in] x1 First object.
 * @param[in] x2 Second object.
 * @param[in] x3 Third object.
 *
 * <tt>x1</tt> must be convertible to T1, <tt>x2</tt> to T2, <tt>x3</tt> to T3.  For example, a Triple of three
 * <tt>int</tt> values can be constructed with three <tt>double</tt> values.
 */

/*!
 * @var T1 Triple::i1
 * @brief First value of triple.
 *
 * signature T1 Triple::i1;
 */

/*!
 * @var T2 Triple::i2
 * @brief Second value of triple.
 *
 * signature T2 Triple::i2;
 */

/*!
 * @var T3 Triple::i3
 * @brief Third value of triple.
 *
 * signature T3 Triple::i3;
 */

template <typename T1, typename T2 = T1, typename T3 = T1, typename TSpec = void>
struct Triple
{
    // ------------------------------------------------------------------------
    // Members
    // ------------------------------------------------------------------------

    T1 i1;
    T2 i2;
    T3 i3;

    // ------------------------------------------------------------------------
    // Constructors
    // ------------------------------------------------------------------------

    inline Triple() : i1(T1()), i2(T2()), i3(T3()) {}

    inline Triple(Triple const & _p)
            : i1(_p.i1), i2(_p.i2), i3(_p.i3) {}

    inline Triple(T1 const & _i1, T2 const & _i2, T3 const & _i3)
            : i1(_i1), i2(_i2), i3(_i3) {}

    template <typename T1_, typename T2_, typename T3_, typename TSpec__>
    inline Triple(Triple<T1_, T2_, T3_, TSpec__> const & _p)
            : i1(getValueI1(_p)), i2(getValueI2(_p)), i3(getValueI3(_p)) {}

    // TODO(holtgrew): Move comparison operators to global functions?
    inline bool
    operator==(Triple const & other) const
    {
        return i1 == other.i1 && i2 == other.i2 && i3 == other.i3;
    }

    inline bool
    operator<(Triple const & other) const
    {
        if (i1 < other.i1)
            return true;
        if (i1 == other.i1 && i2 < other.i2)
            return true;
        if (i1 == other.i1 && i2 == other.i2 && i3 < other.i3)
                return true;
        return false;
    }
};

// ============================================================================
// Metafunctions
// ============================================================================

// -----------------------------------------------------------------------
// Metafunction LENGTH
// -----------------------------------------------------------------------

/*!
 * @mfn Triple#LENGTH
 * @brief Return (only type-depending) length of a triple: 3.
 *
 * @signature LENGTH<TTriple>::VALUE
 *
 * @tparam TTriple The Triple specialization to get the length of.
 *
 * @return VALUE Length of the triple (always 3).
 */

template <typename T1, typename T2, typename T3, typename TSpec>
struct LENGTH<Triple<T1, T2, T3, TSpec> >
{
    enum { VALUE = 3 };
};

// Const variant is mapped to non-const.

// -----------------------------------------------------------------------
// Metafunction Value
// -----------------------------------------------------------------------

/*!
 * @mfn Triple#Value
 * @brief Return i<sup>th</sup> type of the triple.
 *
 * @signature Value<TTriple, I>::Type;
 *
 * @tparam TTriple The Triple to return the <tt>I</tt>-th value of.
 * @tparam I       The index of the value to return, one of 1, 2, or 3.
 */

template <typename T1, typename T2, typename T3, typename TSpec>
struct Value<Triple<T1, T2, T3, TSpec>, 1>
{
    typedef T1 Type;
};

template <typename T1, typename T2, typename T3, typename TSpec>
struct Value<Triple<T1, T2, T3, TSpec>, 2>
{
    typedef T2 Type;
};

template <typename T1, typename T2, typename T3, typename TSpec>
struct Value<Triple<T1, T2, T3, TSpec>, 3 >
{
    typedef T3 Type;
};

// -----------------------------------------------------------------------
// Metafunction Spec
// -----------------------------------------------------------------------

/*!
 * @mfn Triple#Spec
 * @brief Return specialization tag.
 *
 * @signature Spec<TTriple>::Type;
 *
 * @tparam TTriple The Triple specialization to query for the specialization tag.
 *
 * @return Type The specialization type.
 */

template <typename T1, typename T2, typename T3, typename TSpec>
struct Spec<Triple<T1, T2, T3, TSpec> >
{
    typedef TSpec Type;
};

// ============================================================================
// Functions
// ============================================================================

// -----------------------------------------------------------------------
// Function operator<<();  Stream Output.
// -----------------------------------------------------------------------

template <typename TTarget, typename T1, typename T2, typename T3, typename TSpec>
inline void
write(TTarget &target, Triple<T1, T2, T3, TSpec> const & p)
{
    write(target, "< ");
    write(target, getValueI1(p));
    write(target, " , ");
    write(target, getValueI2(p));
    write(target, " , ");
    write(target, getValueI3(p));
    write(target, " >");
}

template <typename TStream, typename T1, typename T2, typename T3, typename TSpec>
inline TStream &
operator<<(TStream & target,
           Triple<T1, T2, T3, TSpec> const & source)
{
    typename DirectionIterator<TStream, Output>::Type it = directionIterator(target, Output());
    write(it, source);
    return target;
}

// -----------------------------------------------------------------------
// Function getValueIX()
// -----------------------------------------------------------------------

/*!
 * @fn Triple#getValueI1
 * @brief The get-value of the Triple's first entry.
 *
 * @signature T1 getValue(triple);
 *
 * @param[in] triple The triple to get entry from.
 *
 * @return T1 The first entry of the Triple.
 */

template <typename T1, typename T2, typename T3, typename TSpec>
inline T1
getValueI1(Triple<T1, T2, T3, TSpec> const & triple)
{
    return triple.i1;
}

/*!
 * @fn Triple#getValueI2
 * @brief The get-value of the Triple's second entry.
 *
 * @signature T2 getValue(triple);
 *
 * @param[in] triple The triple to get entry from.
 *
 * @return T2 The second entry of the Triple.
 */

template <typename T1, typename T2, typename T3, typename TSpec>
inline T2
getValueI2(Triple<T1, T2, T3, TSpec> const & triple)
{
    return triple.i2;
}

/*!
 * @fn Triple#getValueI3
 * @brief The get-value of the Triple's third entry.
 *
 * @signature T3 getValue(triple);
 *
 * @param[in] triple The triple to get entry from.
 *
 * @return T3 The third entry of the Triple.
 */

template <typename T1, typename T2, typename T3, typename TSpec>
inline T3
getValueI3(Triple<T1, T2, T3, TSpec> const & triple)
{
    return triple.i3;
}

// -----------------------------------------------------------------------
// Function assignValueIX()
// -----------------------------------------------------------------------

/*!
 * @fn Triple#assignValueI1
 * @brief Set first entry of a triple.
 *
 * @signature void assignValueI1(triple, val);
 *
 * @param[in] triple The triple to get entry from.
 * @param[in] val    Set the value of the Triple's first entry.
 */

template <typename T1, typename T2, typename T3, typename TSpec, typename T>
inline void assignValueI1(Triple<T1, T2, T3, TSpec> & triple, T const & _i)
{
    triple.i1 = _i;
}

/*!
 * @fn Triple#assignValueI2
 * @brief Set second entry of a triple.
 *
 * @signature void assignValueI2(triple, val);
 *
 * @param[in] triple The triple to get entry from.
 * @param[in] val    Set the value of the Triple's second entry.
 */

template <typename T1, typename T2, typename T3, typename TSpec, typename T>
inline void assignValueI2(Triple<T1, T2, T3, TSpec> & triple, T const & _i)
{
    triple.i2 = _i;
}

/*!
 * @fn Triple#assignValueI3
 * @brief Set third entry of a triple.
 *
 * @signature void assignValueI3(triple, val);
 *
 * @param[in] triple The triple to get entry from.
 * @param[in] val    Set the value of the Triple's third entry.
 */

template <typename T1, typename T2, typename T3, typename TSpec, typename T>
inline void assignValueI3(Triple<T1, T2, T3, TSpec> & triple, T const & _i)
{
    triple.i3 = _i;
}

// -----------------------------------------------------------------------
// Function setValueIX()
// -----------------------------------------------------------------------

/*!
 * @fn Triple#setValueI1
 * @brief Set first entry of a triple.
 *
 * @signature void setValueI1(triple, val);
 *
 * @param[in] triple The triple to get entry from.
 * @param[in] val    Set the value of the Triple's first entry.
 */

template <typename T1, typename T2, typename T3, typename TSpec, typename T>
inline void setValueI1(Triple<T1, T2, T3, TSpec> & triple, T const & _i)
{
    set(triple.i1, _i);
}

/*!
 * @fn Triple#setValueI2
 * @brief Set second entry of a triple.
 *
 * @signature void setValueI2(triple, val);
 *
 * @param[in] triple The triple to get entry from.
 * @param[in] val    Set the value of the Triple's second entry.
 */

template <typename T1, typename T2, typename T3, typename TSpec, typename T>
inline void setValueI2(Triple<T1, T2, T3, TSpec> & triple, T const & _i)
{
    set(triple.i2, _i);
}

/*!
 * @fn Triple#setValueI3
 * @brief Set third entry of a triple.
 *
 * @signature void setValueI3(triple, val);
 *
 * @param[in] triple The triple to get entry from.
 * @param[in] val    Set the value of the Triple's third entry.
 */

template <typename T1, typename T2, typename T3, typename TSpec, typename T>
inline void setValueI3(Triple<T1, T2, T3, TSpec> & triple, T const & _i)
{
    set(triple.i3, _i);
}

// -----------------------------------------------------------------------
// Function moveValueIX()
// -----------------------------------------------------------------------

template <typename T1, typename T2, typename T3, typename TSpec, typename T>
inline void moveValueI1(Triple<T1, T2, T3, TSpec> & triple, T const & _i)
{
    move(triple.i1, _i);
}

template <typename T1, typename T2, typename T3, typename TSpec, typename T>
inline void moveValueI2(Triple<T1, T2, T3, TSpec> & triple, T const & _i)
{
    move(triple.i2, _i);
}

template <typename T1, typename T2, typename T3, typename TSpec, typename T>
inline void moveValueI3(Triple<T1, T2, T3, TSpec> & triple, T const & _i)
{
    move(triple.i3, _i);
}

// -----------------------------------------------------------------------
// Function operator<()
// -----------------------------------------------------------------------

template <
    typename L1, typename L2, typename L3, typename LPack,
    typename R1, typename R2, typename R3, typename RPack>
inline bool
operator<(Triple<L1, L2, L3, LPack> const & _left,
          Triple<R1, R2, R3, RPack> const & _right)
{
    return _left.i1 < _right.i1 || (_left.i1 == _right.i1 && _left.i2 < _right.i2) || (_left.i1 == _right.i1 && _left.i2 == _right.i2 && _left.i3 < _right.i3);
}

// -----------------------------------------------------------------------
// Function operator>()
// -----------------------------------------------------------------------

template <
    typename L1, typename L2, typename L3, typename LPack,
    typename R1, typename R2, typename R3, typename RPack>
inline bool
operator>(Triple<L1, L2, L3, LPack> const & _left,
          Triple<R1, R2, R3, RPack> const & _right)
{
    return _left.i1 > _right.i1 || (_left.i1 == _right.i1 && _left.i2 > _right.i2) || (_left.i1 == _right.i1 && _left.i2 == _right.i2 && _left.i3 > _right.i3);
}

// -----------------------------------------------------------------------
// Function operator<=()
// -----------------------------------------------------------------------

template <
    typename L1, typename L2, typename L3, typename LPack,
    typename R1, typename R2, typename R3, typename RPack>
inline bool
operator<=(Triple<L1, L2, L3, LPack> const & _left,
           Triple<R1, R2, R3, RPack> const & _right)
{
    return !operator>(_left, _right);
}

// -----------------------------------------------------------------------
// Function operator==()
// -----------------------------------------------------------------------

template <
    typename L1, typename L2, typename L3, typename LPack,
    typename R1, typename R2, typename R3, typename RPack>
inline bool
operator==(Triple<L1, L2, L3, LPack> const & _left,
           Triple<R1, R2, R3, RPack> const & _right)
{
    return _left.i1 == _right.i1 && _left.i2 == _right.i2 && _left.i3 == _right.i3;
}

// -----------------------------------------------------------------------
// Function operator>()
// -----------------------------------------------------------------------

template <
    typename L1, typename L2, typename L3, typename LPack,
    typename R1, typename R2, typename R3, typename RPack>
inline bool
operator>=(Triple<L1, L2, L3, LPack> const & _left,
           Triple<R1, R2, R3, RPack> const & _right)
{
    return !operator<(_left, _right);
}

// -----------------------------------------------------------------------
// Function operator!=()
// -----------------------------------------------------------------------

template <
    typename L1, typename L2, typename L3, typename LPack,
    typename R1, typename R2, typename R3, typename RPack>
inline bool
operator!=(Triple<L1, L2, L3, LPack> const & _left,
           Triple<R1, R2, R3, RPack> const & _right)
{
    return !operator==(_left, _right);
}
}  // namespace seqan

#endif  // #ifndef SEQAN_INCLUDE_SEQAN_BASIC_TRIPLE_BASE_H_
