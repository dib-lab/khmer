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
// Pair base class.
// ==========================================================================

// TODO(holtgrew): What about move construction? Useful for pairs of strings and such. Tricky to implement since ints have no move constructor, for example.

#ifndef SEQAN_CORE_INCLUDE_SEQAN_BASIC_PAIR_BASE_H_
#define SEQAN_CORE_INCLUDE_SEQAN_BASIC_PAIR_BASE_H_

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

/*!
 * @class Pair
 * @extends ComparableConcept
 * @headerfile <seqan/basic.h>
 * @brief Store two arbitrary objects.
 *
 * @signature template <typename T1, typename T2, typename TSpec>
 *            class Pair;
 *
 * @tparam T1    The type of the first member.
 * @tparam T2    The type of the second member.
 * @tparam TSpec Tag used for the specialization.
 */

/*!
 * @fn Pair#Pair
 * @brief Default and copy construction and construction for two values.
 *
 * @signature Pair::Pair();
 * @signature Pair::Pair(other);
 * @signature Pair::Pair(x1, x2);
 *
 * @param other The other Pair object to copy from.
 * @param x1    Copied to first member.
 * @param x2    Copied to second member
 */

/*!
 * @var T1 Pair::i1
 * @brief First member
 */

/*!
 * @var T2 Pair::i2
 * @brief Second member
 */

/**
.Class.Pair:
..cat:Aggregates
..concept:Concept.AggregateConcept
..summary:Stores two arbitrary objects.
..signature:Pair<T1[, T2[, TSpec]]>
..param.T1:The type of the first object.
..param.T2:The type of the second object.
...default:$T1$
..param.TSpec:The specializing type.
...default:$void$, no packing (faster access).
.Memfunc.Pair#Pair:
..class:Class.Pair
..summary:Constructor
..signature:Pair<T1, T2[, TSpec]> ()    
..signature:Pair<T1, T2[, TSpec]> (pair)
..signature:Pair<T1, T2[, TSpec]> (i1, i2)
..param.pair:Other Pair object. (copy constructor)
..param.i1:T1 object.
..param.i2:T2 object.
.Memvar.Pair#i1:
..class:Class.Pair
..summary:T1 object
.Memvar.Pair#i2:
..class:Class.Pair
..summary:T2 object
..include:seqan/basic.h
*/

// TODO(holtgrew): Should default specs be specialized with void or Default?
// TODO(holtgrew): Move construction, will be a bit tricky, either with enable_if or with 4 base classes and all constructors are forwarded there.

template <typename T1, typename T2 = T1, typename TSpec = void>
struct Pair
{
    // ------------------------------------------------------------------------
    // Members
    // ------------------------------------------------------------------------

    T1 i1;
    T2 i2;

    // ------------------------------------------------------------------------
    // Constructors
    // ------------------------------------------------------------------------

    Pair() : i1(T1()), i2(T2()) {}

    template <typename T1_, typename T2_>
    Pair(Pair<T1_, T2_> const & _p) : i1(_p.i1), i2(_p.i2) {}

    inline
    Pair(T1 const & _i1, T2 const & _i2) : i1(_i1), i2(_i2) {}

    template <typename T1_, typename T2_, typename TSpec__>
    // TODO(holtgrew): explicit?
    inline Pair(Pair<T1_, T2_, TSpec__> const &_p)
            : i1(getValueI1(_p)), i2(getValueI2(_p))
    {}
};

// ============================================================================
// Metafunctions
// ============================================================================

// -----------------------------------------------------------------------
// Metafunction LENGTH
// -----------------------------------------------------------------------

/*!
 * @mfn Pair#LENGTH
 * @brief Return number of members in a Pair (2).
 *
 * @signature LENGTH<TPair>::VALUE;
 *
 * @tparam TPair The Pair specialization.
 *
 * @return VALUE The number of element in a Pair (2).
 */

///.Metafunction.LENGTH.param.T.type:Class.Pair
///.Metafunction.LENGTH.class:Class.Pair

template <typename T1, typename T2, typename TSpec>
struct LENGTH<Pair<T1, T2, TSpec> >
{
    enum { VALUE = 2 };
};

// Const variant is mapped to non-const by default implementation.

// ----------------------------------------------------------------------------
// Metafunction Value
// ----------------------------------------------------------------------------

/*!
 * @mfn Pair#Value
 * @brief Return type of the i-th value.
 *
 * @signature Value<TTuple, I>::Type;
 *
 * @tparam TTuple Tuple specialization to get the type of.
 * @tparam I      The index of the member to get (1 or 2).
 *
 * @return Type Result type.
 */

/**
.Metafunction.Value
..class:Class.Pair
..class:Class.Triple
..class:Class.Tuple
..signature:Value<TTuple, POSITION>::Type
..param.TTuple:@Class.Pair@, @Class.Triple@, or @Class.Tuple@ to return value from.
...type:Class.Pair
...type:Class.Triple
...type:Class.Tuple
..param.POSITION:Position of the type to query.
...type:nolink:$int$
 */

template <typename T1, typename T2, typename TSpec>
struct Value<Pair<T1, T2, TSpec>, 1>
{
    typedef T1 Type;
};

template <typename T1, typename T2, typename TSpec>
struct Value<Pair<T1, T2, TSpec>, 2>
{
        typedef T2 Type;
};

// ----------------------------------------------------------------------------
// Metafunction Spec
// ----------------------------------------------------------------------------

/*!
 * @mfn Pair#Spec
 * @brief Return specialization tag.
 *
 * @signature Spec<TPair>::Type;
 *
 * @tparam TPair The Pair specialization.
 *
 * @return Type The resulting type.
 */

///.Metafunction.Spec.param.T.type:Class.Pair
///.Metafunction.Spec.class:Class.Pair

template <typename T1, typename T2, typename TSpec>
struct Spec<Pair<T1, T2, TSpec> >
{
    typedef TSpec Type;
};

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function set().
// ----------------------------------------------------------------------------

template <typename T1, typename T2, typename TSpec>
inline void
set(Pair<T1, T2, TSpec> & p1, Pair<T1, T2, TSpec> & p2)
{
    set(p1.i1, p2.i1);
    set(p1.i2, p2.i2);
}

// ----------------------------------------------------------------------------
// Function move().
// ----------------------------------------------------------------------------

template <typename T1, typename T2, typename TSpec>
inline void
move(Pair<T1, T2, TSpec> & p1, Pair<T1, T2, TSpec> & p2)
{
    move(p1.i1, p2.i1);
    move(p1.i2, p2.i2);
}

// ----------------------------------------------------------------------------
// Function operator<<();  Stream Output.
// ----------------------------------------------------------------------------

template <typename T1, typename T2, typename TSpec>
inline
std::ostream & operator<<(std::ostream & out, Pair<T1, T2, TSpec> const & p)
{
    // TODO(holtgrew): Incorporate this into new stream concept? Adapt from stream module?
    out << "< " << getValueI1(p) << " , " << getValueI2(p) << " >";
    return out;
}

// -----------------------------------------------------------------------
// Function getValueIX()
// -----------------------------------------------------------------------

/*!
 * @fn Pair#getValueI1
 * @brief The get-value of the Pair's first entry.
 *
 * @signature T1 getValue(pair);
 *
 * @param pair The pair to get entry from.
 *
 * @return T1 The first entry of the Pair.
 */

// There can be no getValue with index since T1 can be != T2.

template <typename T1, typename T2, typename TSpec>
inline T1 getValueI1(Pair<T1, T2, TSpec> const & pair)
{
    return pair.i1;
}

/*!
 * @fn Pair#getValueI2
 * @brief The get-value of the Pair's second entry.
 *
 * @signature T2 getValue(pair);
 *
 * @param pair The pair to get entry from.
 *
 * @return T2 The second entry of the Pair.
 */

template <typename T1, typename T2, typename TSpec>
inline T2 getValueI2(Pair<T1, T2, TSpec> const & pair)
{
    return pair.i2;
}

// -----------------------------------------------------------------------
// Function assignValueIX()
// -----------------------------------------------------------------------

/*!
 * @fn Pair#assignValueI1
 * @brief Set first entry of a pair.
 *
 * @signature void assignValueI1(pair, val);
 *
 * @param pair The pair to get entry from.
 * @param val  Set the value of the Pair's first entry.
 */

// Cannot be assignValue with index since T1 can be != T2.

template <typename T1, typename T2, typename TSpec, typename T>
inline void assignValueI1(Pair<T1, T2, TSpec> & pair, T const & _i)
{
    pair.i1 = _i;
}

/*!
 * @fn Pair#assignValueI2
 * @brief Set second entry of a pair.
 *
 * @signature void assignValueI1(pair, val);
 *
 * @param pair The pair to get entry from.
 * @param val  Set the value of the Pair's second entry.
 */

template <typename T1, typename T2, typename TSpec, typename T>
inline void assignValueI2(Pair<T1, T2, TSpec> & pair, T const & _i)
{
    pair.i2 = _i;
}

// -----------------------------------------------------------------------
// Function setValueIX()
// -----------------------------------------------------------------------

/*!
 * @fn Pair#setValueI1
 * @brief Set first entry of a pair.
 *
 * @signature void setValueI1(pair, val);
 *
 * @param pair The pair to get entry from.
 * @param val  Set the value of the Pair's first entry.
 */

// Cannot be setValue with index since T1 can be != T2.

template <typename T1, typename T2, typename TSpec, typename T>
inline void setValueI1(Pair<T1, T2, TSpec> & pair, T const & _i)
{
    set(pair.i1, _i);
}

/*!
 * @fn Pair#setValueI2
 * @brief Set second entry of a pair.
 *
 * @signature void setValueI1(pair, val);
 *
 * @param pair The pair to get entry from.
 * @param val  Set the value of the Pair's second entry.
 */

template <typename T1, typename T2, typename TSpec, typename T>
inline void setValueI2(Pair<T1, T2, TSpec> & pair, T const & _i)
{
    set(pair.i2, _i);
}

// -----------------------------------------------------------------------
// Function moveValueIX()
// -----------------------------------------------------------------------

// Cannot be moveValue with index since T1 can be != T2.

template <typename T1, typename T2, typename TSpec, typename T>
inline void moveValueI1(Pair<T1, T2, TSpec> & pair, T & _i)
{
    move(pair.i1, _i);
}

template <typename T1, typename T2, typename TSpec, typename T>
inline void moveValueI2(Pair<T1, T2, TSpec> & pair, T & _i)
{
    move(pair.i2, _i);
}

// -----------------------------------------------------------------------
// Function operator<()
// -----------------------------------------------------------------------

template <typename L1, typename L2, typename LPack, typename R1, typename R2, typename RPack>
inline bool
operator<(Pair<L1, L2, LPack> const & _left,
          Pair<R1, R2, RPack> const & _right)
{
    return (_left.i1 < _right.i1) || (_left.i1 == _right.i1 && _left.i2 < _right.i2);
}

// -----------------------------------------------------------------------
// Function operator>()
// -----------------------------------------------------------------------

template <typename L1, typename L2, typename LPack, typename R1, typename R2, typename RPack>
inline bool
operator>(Pair<L1, L2, LPack> const & _left,
          Pair<R1, R2, RPack> const & _right)
{
    return (_left.i1 > _right.i1) || (_left.i1 == _right.i1 && _left.i2 > _right.i2);
}

// -----------------------------------------------------------------------
// Function operator==()
// -----------------------------------------------------------------------

template <typename L1, typename L2, typename LPack, typename R1, typename R2, typename RPack>
inline bool
operator==(Pair<L1, L2, LPack> const & _left,
           Pair<R1, R2, RPack> const & _right)
{
    return _left.i1 == _right.i1 && _left.i2 == _right.i2;
}

// -----------------------------------------------------------------------
// Function operator<=()
// -----------------------------------------------------------------------

template <typename L1, typename L2, typename LPack, typename R1, typename R2, typename RPack>
inline bool
operator<=(Pair<L1, L2, LPack> const & _left,
           Pair<R1, R2, RPack> const & _right)
{
    return !operator>(_left, _right);
}

// -----------------------------------------------------------------------
// Function operator>=()
// -----------------------------------------------------------------------

template <typename L1, typename L2, typename LPack, typename R1, typename R2, typename RPack>
inline bool
operator>=(Pair<L1, L2, LPack> const & _left,
           Pair<R1, R2, RPack> const & _right)
{
    return !operator<(_left, _right);
}

// -----------------------------------------------------------------------
// Function operator!=()
// -----------------------------------------------------------------------

template <typename L1, typename L2, typename LPack, typename R1, typename R2, typename RPack>
inline bool
operator!=(Pair<L1, L2, LPack> const & _left,
           Pair<R1, R2, RPack> const & _right)
{
    return !operator==(_left, _right);
}

}  // namespace seqan

#endif  // #ifndef SEQAN_CORE_INCLUDE_SEQAN_BASIC_PAIR_BASE_H_
