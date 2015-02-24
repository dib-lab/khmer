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
// Author: Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>
// ==========================================================================
// Implementation of the SeqAn accumulator code.
//
// While Boost also has a versatile accumulator library, we create our own
// since it is much simpler and focused on the requirements in sequence
// analysis algorithms than then Boost one.
// ==========================================================================

#ifndef INCLUDE_SEQAN_MISC_ACCUMULATORS_H_
#define INCLUDE_SEQAN_MISC_ACCUMULATORS_H_

#include <seqan/sequence.h>

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

template <typename T, typename TTag> struct Result;

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

struct AccuAverage_;
typedef Tag<AccuAverage_> AccuAverage;

struct AccuSum_;
typedef Tag<AccuSum_> AccuSum;

struct AccuCount_;
typedef Tag<AccuCount_> AccuCount;

/*!
 * @class Accumulator
 * @headerfile <seqan/misc/accumulators.h>
 * @brief Accumulator base class.
 *
 * @signature template <typename TValue, typename TSpec>
 *            struct Accumulator;
 *
 * @tparam TSpec  The specialization tag.
 * @tparam TValue The type of the values to accumulate.
 *
 * Accumulators are for computing statistics on streams of values.
 *
 * Currently, this is only meant for accumulating integers.
 */

template <typename TValue, typename TSpec>
struct Accumulator;

/*!
 * @class AverageAccumulator
 * @extends Accumulator
 * @headerfile <seqan/misc/accumulators.h>
 * @brief Accumulator for computing averages.
 *
 * @signature template <typename TValue>
 *            struct Accumulator<TValue, AccuAverage>;
 *
 * @tparam TValue The type of the values to compute the average of.
 *
 * The average of an empty sequence is defined to be 0.
 *
 * @section Examples
 *
 * This program shows how to use the Average Accumulator.
 *
 * @code{.cpp}
 * Accumulator<int, AccuAverage> acc;
 * push(acc, 1);
 * push(acc, 2);
 * push(acc, 3);
 * std::cout << "average: " << average(acc) << "\n"
 *           << "sum:     " << sum(acc) << "\n"
 *           << "count:   " << count(acc) << "\n";
 * @endcode
 *
 * The output is then:
 *
 * @code{.console}
 * average: 2
 * sum:     6
 * count:   3
 * @endcode
 */

template <typename TValue>
struct Accumulator<TValue, AccuAverage>
{
    typedef typename Result<Accumulator<TValue, AccuAverage>, AccuSum>::Type TAccuSum_;

    TAccuSum_ sum_;
    unsigned count_;

    Accumulator() : sum_(0), count_(0) {}
};

// ============================================================================
// Metafunctions
// ============================================================================

// ----------------------------------------------------------------------------
// Metafunction Value                                       Average Accumulator
// ----------------------------------------------------------------------------

/*!
 * @mfn Accumulator#Value
 * @brief Return the type of the values to accumulate.
 *
 * @signature Value<TAccumulator>::Type;
 *
 * @tparam TAccumulator The Accumulator type to query.
 *
 * @return Type The value type.
 */

template <typename TValue>
struct Value<Accumulator<TValue, AccuAverage> >
{
    typedef double Type;
};

template <typename TValue>
struct Value<Accumulator<TValue, AccuAverage> const > : Value<Accumulator<TValue, AccuAverage> >
{};

// ----------------------------------------------------------------------------
// Metafunction Result                                      Average Accumulator
// ----------------------------------------------------------------------------

// TODO(holtgrew): This could probably go to basic, as a part of an "AlgorithmState" concept?

/*!
 * @mfn Accumulator#Result
 * @brief Return the type for accumulation results.
 *
 * @signature Result<TAccumulator>::Type;
 *
 * @tparam TAccumulator The Accumulator type to query.
 *
 * @return Type The result type.
 */

template <typename T, typename TTag = void>
struct Result;

// ----------------------------------------------------------------------------
// Metafunction Result                                      Average Accumulator
// ----------------------------------------------------------------------------

template <typename TValue>
struct Result<Accumulator<TValue, AccuAverage>, AccuAverage>
{
    typedef double Type;
};

template <typename TValue>
struct Result<Accumulator<TValue, AccuAverage> const, AccuAverage> : Result<Accumulator<TValue, AccuAverage>, AccuAverage>
{};

template <typename TValue>
struct Result<Accumulator<TValue, AccuAverage>, AccuCount>
{
    typedef unsigned Type;
};

template <typename TValue>
struct Result<Accumulator<TValue, AccuAverage> const, AccuCount> : Result<Accumulator<TValue, AccuAverage>, AccuCount>
{};

template <typename TValue>
struct Result<Accumulator<TValue, AccuAverage>, AccuSum>
{
    typedef typename IfC<Is<IntegerConcept<TValue> >::VALUE,
                         __int64,
                         double>::Type Type;
};

template <typename TValue>
struct Result<Accumulator<TValue, AccuAverage> const, AccuSum> : Result<Accumulator<TValue, AccuAverage>, AccuSum>
{};

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function clear()                                                 Accumulator
// ----------------------------------------------------------------------------

// ----------------------------------------------------------------------------
// Function push()                                                  Accumulator
// ----------------------------------------------------------------------------

/*!
 * @fn Accumulator#push
 * @brief Include value into sequence of values to accumulate.
 *
 * @signature void push(acc, x);
 *
 * @param[in,out] acc The Accumulator to push the value to.
 * @param[in]     x   The value to include in the accumulation (@link IntegerConcept @endlink).
 */

// ----------------------------------------------------------------------------
// Function clear()                                                 Accumulator
// ----------------------------------------------------------------------------

/*!
 * @fn Accumulator#clear
 * @brief Clear the current accumulator state.
 *
 * @signature void clear(acc);
 *
 * @param[in,out] acc The Accumulator to clear.
 */

// ----------------------------------------------------------------------------
// Function clear()                                         Average Accumulator
// ----------------------------------------------------------------------------

template <typename TValue>
inline void
clear(Accumulator<TValue, AccuAverage> & accumulator)
{
    accumulator.sum_ = 0;
    accumulator.count_ = 0;
}

// ----------------------------------------------------------------------------
// Function push()                                          Average Accumulator
// ----------------------------------------------------------------------------

template <typename TValue, typename TValue2>
inline void
push(Accumulator<TValue, AccuAverage> & acc, TValue2 value)
{
    typedef typename Result<Accumulator<TValue, AccuAverage>, AccuSum>::Type TAccuSum;
    acc.sum_ += static_cast<TAccuSum>(value);
    acc.count_ += 1;
}

// ----------------------------------------------------------------------------
// Function average()                                       Average Accumulator
// ----------------------------------------------------------------------------

/*!
 * @fn AverageAccumulator#average
 * @brief Return the average of the included values.
 *
 * @signature TResult average(acc);
 *
 * @param[in] acc The Accumulator to compute the average for.
 *
 * @return TResult The average of the values (Metafunction: @link Accumulator#Result @endlink).
 */

template <typename TValue>
inline typename Result<Accumulator<TValue, AccuAverage>, AccuAverage>::Type
average(Accumulator<TValue, AccuAverage> const & acc)
{
    typedef typename Result<Accumulator<TValue, AccuAverage>, AccuAverage>::Type TResult;
    if (acc.count_ == 0u)
        return 0;
    return static_cast<TResult>(acc.sum_ / static_cast<TResult>(acc.count_));
}

// ----------------------------------------------------------------------------
// Function sum()                                           Average Accumulator
// ----------------------------------------------------------------------------

/*!
 * @fn AverageAccumulator#sum
 * @brief Return the sum of the included values.
 *
 * @signature TResult sum(acc);
 *
 * @param[in] acc The Accumulator to compute the sum for.
 *
 * @return TResult The sum of the values (Metafunction: @link Accumulator#Result @endlink).
 */

template <typename TValue>
inline typename Result<Accumulator<TValue, AccuAverage>, AccuSum>::Type
sum(Accumulator<TValue, AccuAverage> const & acc)
{
    return acc.sum_;
}

// ----------------------------------------------------------------------------
// Function count()                                         Average Accumulator
// ----------------------------------------------------------------------------

/*!
 * @fn AverageAccumulator#count
 * @brief Return the number of included values.
 *
 * @signature TResult count(acc);
 *
 * @param[in] count The number of values pushed to the accumulator.
 *
 * @return TResult The number of pushed values (Metafunction: @link Accumulator#Result @endlink).
 */

template <typename TValue>
inline typename Result<Accumulator<TValue, AccuAverage>, AccuCount>::Type
count(Accumulator<TValue, AccuAverage> const & acc)
{
    return acc.count_;
}

}  // namespace seqan

#endif  // #ifndef INCLUDE_SEQAN_MISC_ACCUMULATORS_H_
