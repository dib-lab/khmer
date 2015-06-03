// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2015, Knut Reinert, FU Berlin
// Copyright (c) 2013 NVIDIA Corporation
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
// Tuple base class.
// ==========================================================================

// TODO(holtgrew): What about move construction? Useful for pairs of strings and such. Tricky to implement since ints have no move constructor, for example.

#ifndef SEQAN_INCLUDE_SEQAN_BASIC_TUPLE_BASE_H_
#define SEQAN_INCLUDE_SEQAN_BASIC_TUPLE_BASE_H_

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

template <typename TValue>
struct StoredTupleValue_
{
    typedef TValue Type;
};

template <typename TValue, typename TSpec>
struct StoredTupleValue_< SimpleType<TValue, TSpec> >
{
    typedef TValue Type;
};

/*!
 * @class Tuple
 * @headerfile <seqan/basic.h>
 * @brief A constant-size tuple of the same type.
 *
 * @signature template <typename TValue, unsigned SIZE[, typename TSpec]>
 *            class Tuple;
 *
 * @tparam TValue The value that the tuple should be over.
 * @tparam SIZE   The number of entries in the tuple.
 * @tparam TSpec  Specialization tag, optional;  defaults to void.
 */

/*!
 * @fn Tuple::i
 * @brief Array of tuple value.
 *
 * @signature TValue Tuple::i[SIZE];
 *
 * The array of the tuple's values.
 */

template <typename TValue, unsigned SIZE, typename TSpec = void>
struct Tuple
{
    // -----------------------------------------------------------------------
    // Members
    // -----------------------------------------------------------------------

    typename StoredTupleValue_<TValue>::Type i[SIZE];

    // -----------------------------------------------------------------------
    // Subscription Operators;  Have to be declared in class.
    // -----------------------------------------------------------------------

    // TODO(holtgrew): Return Value<>::Type?

    template <typename TPos>
    SEQAN_HOST_DEVICE inline
    typename StoredTupleValue_<TValue>::Type &
    operator[](TPos k)
    {
        SEQAN_ASSERT_GEQ(static_cast<__int64>(k), 0);
        SEQAN_ASSERT_LT(static_cast<__int64>(k), static_cast<__int64>(SIZE));
        return i[k];
    }

    template <typename TPos>
    SEQAN_HOST_DEVICE inline
    typename StoredTupleValue_<TValue>::Type const &
    operator[](TPos k) const
    {
        SEQAN_ASSERT_GEQ(static_cast<__int64>(k), 0);
        SEQAN_ASSERT_LT(static_cast<__int64>(k), static_cast<__int64>(SIZE));
        return i[k];
    }

    // This has to be inline because elements (like this tuple) of packed
    // structs can't be arguments.
    template <typename TPos, typename TValue2>
    inline TValue2
    assignValue(TPos k, TValue2 const source)
    {
        return i[k] = source;
    }
};


#ifdef PLATFORM_WINDOWS
    #pragma pack(push,1)
#endif
template <typename TValue, unsigned SIZE>
struct Tuple<TValue, SIZE, Pack>
{
    // -----------------------------------------------------------------------
    // Members
    // -----------------------------------------------------------------------

    typename StoredTupleValue_<TValue>::Type i[SIZE];

    // -----------------------------------------------------------------------
    // Subscription Operators;  Have to be declared in class.
    // -----------------------------------------------------------------------

    // TODO(holtgrew): Return Value<>::Type?

    template <typename TPos>
    inline typename StoredTupleValue_<TValue>::Type &
    operator[](TPos k)
    {
        SEQAN_ASSERT_GEQ(static_cast<__int64>(k), 0);
        SEQAN_ASSERT_LT(static_cast<__int64>(k), static_cast<__int64>(SIZE));
        return i[k];
    }

    template <typename TPos>
    inline typename StoredTupleValue_<TValue>::Type const &
    operator[](TPos k) const
    {
        SEQAN_ASSERT_GEQ(static_cast<__int64>(k), 0);
        SEQAN_ASSERT_LT(static_cast<__int64>(k), static_cast<__int64>(SIZE));
        return i[k];
    }

    // This has to be inline because elements (like this tuple) of packed
    // structs can't be arguments.
    template <typename TPos, typename TValue2>
    inline TValue2
    assignValue(TPos k, TValue2 const source)
    {
        return i[k] = source;
    }
}
#ifndef PLATFORM_WINDOWS
    __attribute__((packed))
#endif
    ;
#ifdef PLATFORM_WINDOWS
      #pragma pack(pop)
#endif

//template <typename TValue, unsigned SIZE>
//const unsigned Tuple<TValue, SIZE, Pack>::SIZE = SIZE;

// ============================================================================
// Metafunctions
// ============================================================================

// -----------------------------------------------------------------------
// Metafunction LENGTH
// -----------------------------------------------------------------------

/*!
 * @mfn Tuple#LENGTH
 * @brief Return the length of a tuple.
 *
 * @signature LENGTH<TTuple>::VALUE;
 *
 * @tparam TTuple The tuple to query for its length.
 *
 * @return VALUE The length of the tuple.
 */

template <typename TValue, unsigned SIZE, typename TSpec>
struct LENGTH<Tuple<TValue, SIZE, TSpec> >
{
    enum { VALUE = SIZE };
};

// -----------------------------------------------------------------------
// Metafunction Value
// -----------------------------------------------------------------------

/*!
 * @mfn Tuple#Value
 * @brief Return the value type of a tuple.
 *
 * @signature Value<TTuple>::Type;
 *
 * @tparam TTuple The tuple type to query for its value type.
 *
 * @return Type The resulting value type.
 */

template <typename TValue, unsigned SIZE, typename TSpec>
struct Value<Tuple<TValue, SIZE, TSpec> >
{
    typedef TValue Type;
};

// -----------------------------------------------------------------------
// Metafunction Spec
// -----------------------------------------------------------------------

template <typename TValue, unsigned SIZE, typename TSpec>
struct Spec<Tuple<TValue, SIZE, TSpec> >
{
    typedef TSpec Type;
};

// ============================================================================
// Functions
// ============================================================================

// -----------------------------------------------------------------------
// Function operator<<();  Stream Output.
// -----------------------------------------------------------------------

template <typename TTarget, typename TValue, unsigned SIZE, typename TSpec>
inline void
write(TTarget &target, Tuple<TValue, SIZE, TSpec> const &a)
{
    writeValue(target, '[');
    if (SIZE > 0)
        write(target, (TValue)a[0]);
    for (unsigned j = 1; j < SIZE; ++j)
    {
        writeValue(target, ' ');
        write(target, (TValue)a[j]);
    }
    writeValue(target, ']');
}

template <typename TStream, typename TValue, unsigned SIZE, typename TSpec>
inline TStream &
operator<<(TStream & target,
           Tuple<TValue, SIZE, TSpec> const & source)
{
    typename DirectionIterator<TStream, Output>::Type it = directionIterator(target, Output());
    write(it, source);
    return target;
}

// ----------------------------------------------------------------------------
// Function set().
// ----------------------------------------------------------------------------

template <typename TTuple1, typename TTuple2>
struct TupleMoveSetWorkerContext_
{
    TTuple1 & t1;
    TTuple2 & t2;

    TupleMoveSetWorkerContext_(TTuple1 & _t1, TTuple2 & _t2)
            : t1(_t1), t2(_t2)
    {}
};

struct TupleSetWorker_
{
    template <typename TArg>
    static inline void body(TArg & arg, unsigned I)
    {
        set(arg.t1.i[I - 1], arg.t2.i[I - 1]);
    }
};

template <typename TValue, unsigned SIZE, typename TSpec>
inline void
set(Tuple<TValue, SIZE, TSpec> & t1, Tuple<TValue, SIZE, TSpec> const & t2)
{
    typedef Tuple<TValue, SIZE, TSpec> TTuple1;
    typedef Tuple<TValue, SIZE, TSpec> const TTuple2;
    TupleMoveSetWorkerContext_<TTuple1, TTuple2> context(t1, t2);
    Loop<TupleSetWorker_, SIZE>::run(context);
}

template <typename TValue, unsigned SIZE, typename TSpec>
inline void
set(Tuple<TValue, SIZE, TSpec> & t1, Tuple<TValue, SIZE, TSpec> & t2)
{
    set(t1, const_cast<Tuple<TValue, SIZE, TSpec> const &>(t2));
}

// ----------------------------------------------------------------------------
// Function move().
// ----------------------------------------------------------------------------

struct TupleMoveWorker_
{
    template <typename TArg>
    static inline void body(TArg & arg, unsigned I)
    {
        move(arg.t1.i[I - 1], arg.t2.i[I - 1]);
    }
};

template <typename TValue, unsigned SIZE, typename TSpec>
inline void
move(Tuple<TValue, SIZE, TSpec> & t1, Tuple<TValue, SIZE, TSpec> & t2)
{
    typedef Tuple<TValue, SIZE, TSpec> TTuple1;
    typedef Tuple<TValue, SIZE, TSpec> TTuple2;
    TupleMoveSetWorkerContext_<TTuple1, TTuple2> context(t1, t2);
    Loop<TupleMoveWorker_, SIZE>::run(context);
}

// -----------------------------------------------------------------------
// Function assignValue()
// -----------------------------------------------------------------------

template <typename TValue, unsigned SIZE, typename TSpec, typename TPos, typename TValue2>
inline TValue2
assignValue(Tuple<TValue, SIZE, TSpec> & me, TPos k, TValue2 const source)
{
    SEQAN_CHECK((unsigned(k) < SIZE), "Invalid position, k = %u, SIZE = %u.", unsigned(k), unsigned(SIZE));
    return me.i[k] = source;
}

// -----------------------------------------------------------------------
// Function getValue()
// -----------------------------------------------------------------------

template <typename TValue, unsigned SIZE, typename TSpec, typename TPos>
inline TValue
getValue(Tuple<TValue, SIZE, TSpec> & me, TPos k)
{
    SEQAN_CHECK((unsigned(k) < SIZE), "Invalid position, k = %u, SIZE = %u.", unsigned(k), unsigned(SIZE));
    return me.i[k];
}

template <typename TValue, unsigned SIZE, typename TSpec, typename TPos>
inline TValue
getValue(Tuple<TValue, SIZE, TSpec> const & me, TPos k)
{
    SEQAN_CHECK((unsigned(k) < SIZE), "Invalid position, k = %u, SIZE = %u.", unsigned(k), unsigned(SIZE));
    return me.i[k];
}

// -----------------------------------------------------------------------
// Function setValue()
// -----------------------------------------------------------------------

template <typename TValue, unsigned SIZE, typename TSpec, typename TPos, typename TValue2>
inline void
setValue(Tuple<TValue, SIZE, TSpec> & me, TPos k, TValue2 const & source)
{
    SEQAN_CHECK((unsigned(k) < SIZE), "Invalid position, k = %u, SIZE = %u.", unsigned(k), unsigned(SIZE));
    set(me.i[k], source);
}

// -----------------------------------------------------------------------
// Function moveValue()
// -----------------------------------------------------------------------

template <typename TValue, unsigned SIZE, typename TSpec, typename TPos, typename TValue2>
inline void
moveValue(Tuple<TValue, SIZE, TSpec> & me, TPos k, TValue2 & source)
{
    SEQAN_CHECK((unsigned(k) < SIZE), "Invalid position, k = %u, SIZE = %u.", unsigned(k), unsigned(SIZE));
    move(me.i[k], source);
}

// -----------------------------------------------------------------------
// Function shiftLeft()
// -----------------------------------------------------------------------

// TODO(holtgrew): Document!

struct TupleShiftLeftWorker_
{
    template <typename TArg>
    static inline void body(TArg & arg, unsigned I)
    {
        arg[I-1] = arg[I];  // TODO(holtgrew): Do we really want assignment or movement here?
    }
};

template <typename TValue, unsigned SIZE, typename TSpec>
inline void shiftLeft(Tuple<TValue, SIZE, TSpec> &me)
{
    Loop<TupleShiftLeftWorker_, SIZE - 1>::run(me.i);
}

// -----------------------------------------------------------------------
// Function shiftRight()
// -----------------------------------------------------------------------

// TODO(holtgrew): Document!

struct TupleShiftRightWorker_
{
    template <typename TArg>
    static inline void body(TArg & arg, unsigned I)
    {
        arg[I] = arg[I - 1];  // TODO(holtgrew): Do we really want assignment or movement here?
    }
};

template <typename TValue, unsigned SIZE, typename TSpec>
inline void shiftRight(Tuple<TValue, SIZE, TSpec> & me)
{
    LoopReverse<TupleShiftRightWorker_, SIZE - 1>::run(me.i);
}

// -----------------------------------------------------------------------
// Function length()
// -----------------------------------------------------------------------

template <typename TValue, unsigned SIZE, typename TSpec>
inline unsigned length(Tuple<TValue, SIZE, TSpec> const &)
{
    return SIZE;
}

// -----------------------------------------------------------------------
// Function clear()
// -----------------------------------------------------------------------

template <typename TValue, unsigned SIZE, typename TSpec>
inline void clear(Tuple<TValue, SIZE, TSpec> & me)
{
   memset<sizeof(me.i), 0>(&(me.i));
}

// -----------------------------------------------------------------------
// Function operator==()
// -----------------------------------------------------------------------

template <typename TTupleL, typename TTupleR>
struct ComparisonWorkerContext_
{
    int result;
    TTupleL const & left;
    TTupleR const & right;

    ComparisonWorkerContext_(int b, TTupleL const & l, TTupleR const & r)
            : result(b), left(l), right(r)
    {}
};

struct TupleComparisonWorkerEq_
{
    template <typename TArg>
    static inline void body(TArg & arg, unsigned I)
    {
        if (arg.result != 1)
            return;
        if (getValue(arg.left, I - 1) != getValue(arg.right, I - 1))
            arg.result = 0;
    }
};

template <typename TValue, unsigned SIZE, typename TSpecL, typename TSpecR>
inline bool
operator==(Tuple<TValue, SIZE, TSpecL> const & left,
           Tuple<TValue, SIZE, TSpecR> const & right)
{
    typedef Tuple<TValue, SIZE, TSpecL> TTupleL;
    typedef Tuple<TValue, SIZE, TSpecR> TTupleR;
    ComparisonWorkerContext_<TTupleL, TTupleR> context(1, left, right);
    Loop<TupleComparisonWorkerEq_, SIZE>::run(context);
    return context.result == 1;
}

// -----------------------------------------------------------------------
// Function operator!=()
// -----------------------------------------------------------------------


template <typename TValue, unsigned SIZE, typename TSpecL, typename TSpecR>
inline bool
operator!=(Tuple<TValue, SIZE, TSpecL> const & left,
           Tuple<TValue, SIZE, TSpecR> const & right)
{
    return !operator==(left, right);
}

// -----------------------------------------------------------------------
// Function operator<()
// -----------------------------------------------------------------------

struct TupleComparisonWorkerLt_
{
    template <typename TArg>
    static inline void body(TArg & arg, unsigned I)
    {
        if (arg.result != -1)
            return;
        if (arg.left.i[I - 1] == arg.right.i[I - 1])
            return;
        if (arg.left.i[I - 1] < arg.right.i[I - 1])
            arg.result = 1;
        if (arg.left.i[I - 1] > arg.right.i[I - 1])
            arg.result = 0;
    }
};

template <typename TValue, unsigned SIZE, typename TSpecL, typename TSpecR>
inline bool
operator<(Tuple<TValue, SIZE, TSpecL> const & left,
          Tuple<TValue, SIZE, TSpecR> const & right)
{
    typedef Tuple<TValue, SIZE, TSpecL> TTupleL;
    typedef Tuple<TValue, SIZE, TSpecR> TTupleR;
    ComparisonWorkerContext_<TTupleL, TTupleR> context(-1, left, right);
    Loop<TupleComparisonWorkerLt_, SIZE>::run(context);
    return context.result == 1;
}

// -----------------------------------------------------------------------
// Function operator>()
// -----------------------------------------------------------------------

struct TupleComparisonWorkerGt_
{
    template <typename TArg>
    static inline void body(TArg & arg, unsigned I)
    {
        if (arg.result != -1)
            return;
        if (arg.left.i[I - 1] == arg.right.i[I - 1])
            return;
        if (arg.left.i[I - 1] > arg.right.i[I - 1])
            arg.result = 1;
        if (arg.left.i[I - 1] < arg.right.i[I - 1])
            arg.result = 0;
    }
};

template <typename TValue, unsigned SIZE, typename TSpecL, typename TSpecR>
inline bool
operator>(Tuple<TValue, SIZE, TSpecL> const & left,
          Tuple<TValue, SIZE, TSpecR> const & right)
{
    typedef Tuple<TValue, SIZE, TSpecL> TTupleL;
    typedef Tuple<TValue, SIZE, TSpecR> TTupleR;
    ComparisonWorkerContext_<TTupleL, TTupleR> context(-1, left, right);
    Loop<TupleComparisonWorkerGt_, SIZE>::run(context);
    return context.result == 1;
}

// -----------------------------------------------------------------------
// Function operator<=()
// -----------------------------------------------------------------------

template <typename TValue, unsigned SIZE, typename TSpecL, typename TSpecR>
inline bool
operator<=(Tuple<TValue, SIZE, TSpecL> const & left,
           Tuple<TValue, SIZE, TSpecR> const & right)
{
    return !operator>(left, right);
}

// -----------------------------------------------------------------------
// Function operator>=()
// -----------------------------------------------------------------------

template <typename TValue, unsigned SIZE, typename TSpecL, typename TSpecR>
inline bool
operator>=(Tuple<TValue, SIZE, TSpecL> const & left,
           Tuple<TValue, SIZE, TSpecR> const & right)
{
    return !operator<(left, right);
}

// -----------------------------------------------------------------------
// Function operator+()
// -----------------------------------------------------------------------

template <typename TValue, unsigned SIZE, typename TSpecL, typename TSpecR>
inline Tuple<TValue, SIZE, TSpecL>
operator+(Tuple<TValue, SIZE, TSpecL> const & left,
          Tuple<TValue, SIZE, TSpecR> const & right)
{
    Tuple<TValue, SIZE, TSpecL>  tuple;

    for (unsigned j = 0; j < SIZE; ++j)
        tuple[j] = left[j] + right[j];

    return tuple;
}

}  // namespace seqan

#endif  // #ifndef SEQAN_INCLUDE_SEQAN_BASIC_TUPLE_BASE_H_
