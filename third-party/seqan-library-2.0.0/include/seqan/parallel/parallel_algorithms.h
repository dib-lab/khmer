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
// Basic parallel algorithms.
// ==========================================================================

#ifndef SEQAN_PARALLEL_PARALLEL_ALGORITHMS_H_
#define SEQAN_PARALLEL_PARALLEL_ALGORITHMS_H_

namespace seqan {

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function arrayFill
// ----------------------------------------------------------------------------

template <typename TIterator, typename TValue, typename TParallelTag>
inline void
arrayFill(TIterator begin_,
          TIterator end_,
          TValue const & value,
          Tag<TParallelTag> parallelTag)
{
    Splitter<TIterator> splitter(begin_, end_, parallelTag);

    SEQAN_OMP_PRAGMA(parallel for)
    for (int job = 0; job < (int)length(splitter); ++job)
        arrayFill(splitter[job], splitter[job + 1], value, Serial());
}

// ----------------------------------------------------------------------------
// Function sum
// ----------------------------------------------------------------------------

/*!
 * @fn sum
 * @headerfile <seqan/parallel.h>
 * @brief Returns the sum of all elements in a sequence.
 *
 * @signature TValue sum(seq[, parallelTag]);
 *
 * @param[in] seq         The sequence to sum up, o ftype <tt>TSequence</tt>.
 * @param[in] parallelTag Tag to enable/disable parallelism, one of <tt>Serial</tt> and <tt>Parallel</tt>,
 *                        default is <tt>Serial</tt>.
 *
 * @return TValue The sum of the elements in <tt>seq</tt>, of type <tt>Value&lt;TSequence&gt;::Type</tt>.
 *
 * The sequence alphabet must support the <tt>operator+</tt> and conversion from zero.
 *
 * @see partialSum
 */

template <typename TSequence>
inline typename Value<TSequence>::Type
sum(TSequence const &seq, Serial)
{
    typename Iterator<TSequence const>::Type it = begin(seq, Standard());
    typename Iterator<TSequence const>::Type itEnd = end(seq, Standard());
    typename Value<TSequence>::Type sum = 0;
    for (; it != itEnd; ++it)
        sum += *it;
    return sum;
}

template <typename TSequence, typename TParallelTag>
inline typename Value<TSequence>::Type
sum(TSequence const &seq, Tag<TParallelTag> parallelTag)
{
    Splitter<typename Size<TSequence>::Type> splitter(0, length(seq), parallelTag);

    typename Value<TSequence>::Type threadSum = 0;
    SEQAN_OMP_PRAGMA(parallel for reduction(+:threadSum))
    for (int job = 0; job < (int)length(splitter); ++job)
        threadSum += sum(infix(seq, splitter[job], splitter[job + 1]), Serial());
    return threadSum;
}

template <typename TSequence>
inline typename Value<TSequence>::Type
sum(TSequence const &seq)
{
    return sum(seq, Serial());
}

// ----------------------------------------------------------------------------
// Function partialSum
// ----------------------------------------------------------------------------

/*!
 * @fn partialSum
 * @headerfile <seqan/parallel.h>
 * @brief Computes the partial sum of a sequence.
 *
 * @signature TValue partialSum(target, source[, parallelTag]);
 *
 * @param[in]  source      A sequence of elements that should be partially summed.  The sequence alphabet must support
 *                         the <tt>operator+</tt> and conversion from zero, the type is <tt>TSource</tt>.
 * @param[in]  parallelTag Tag to enable/disable parallelism, one of <tt>Serial</tt>, <tt>Parallel</tt>, default is
 *                         <tt>Serial</tt>.
 * @param[out] target      The resulting partial sum.  This sequence will have the same length as <tt>source</tt> and
 *                         contains at position <tt>i</tt> the sum of elements <tt>source[0]</tt>, <tt>source[1]</tt>,
 *                         ..., <tt>source[i]</tt>.
 *
 * @return TValue The sum of all elements in <tt>source</tt>.  The returned value equals the last value in target.
 *                <tt>TValue</tt> is <tt>Value&lt;TSource&gt;::Type</tt>.
 *
 * @see sum
 */

template <typename TTarget, typename TSource, typename TParallelTag>
inline typename Value<TSource>::Type
partialSum(TTarget &target, TSource const &source, Tag<TParallelTag> parallelTag)
{
    typedef typename Value<TSource>::Type TValue;
    typedef typename Size<TSource>::Type TSize;
    typedef typename Iterator<TSource const, Standard>::Type TConstIterator;
    typedef typename Iterator<TTarget, Standard>::Type TIterator;

    resize(target, length(source), Exact());
    if (empty(target))
        return 0;

    Splitter<TSize> splitter(0, length(source), parallelTag);
    String<TValue> localSums;
    resize(localSums, length(splitter), Exact());
    localSums[0] = 0;

    // STEP 1: compute sums of all subintervals (in parallel)
    //
    SEQAN_OMP_PRAGMA(parallel for)
    for (int job = 0; job < (int)length(splitter) - 1; ++job)
        localSums[job + 1] = sum(infix(source, splitter[job], splitter[job + 1]), Serial());

    // STEP 2: compute partial sums (of subinterval sums) to get offsets for each subinterval (sequentially)
    //
    for (int job = 2; job < (int)length(splitter); ++job)
        localSums[job] += localSums[job - 1];

    // STEP 3: compute partial sums of each subinterval starting from offset (in parallel)
    //
    SEQAN_OMP_PRAGMA(parallel for)
    for (int job = 0; job < (int)length(splitter); ++job)
    {
        TConstIterator it = begin(source, Standard()) + splitter[job];
        TConstIterator itEnd = begin(source, Standard()) + splitter[job + 1];
        TIterator dstIt = begin(target, Standard()) + splitter[job];
        TValue sum = localSums[job];
        for (; it != itEnd; ++it, ++dstIt)
        {
            sum += *it;
            *dstIt = sum;
        }
        localSums[job] = sum;
    }

    return back(localSums);
}

template <typename TTarget, typename TParallelTag>
inline typename Value<TTarget>::Type
partialSum(TTarget & target, Tag<TParallelTag> parallelTag)
{
    return partialSum(target, target, parallelTag);
}

template <typename TTarget, typename TSource>
inline typename Value<TSource>::Type
partialSum(TTarget & target, TSource const & source)
{
    return partialSum(target, source, Serial());
}

template <typename TTarget>
inline typename Value<TTarget>::Type
partialSum(TTarget & target)
{
    return partialSum(target, target);
}

// ----------------------------------------------------------------------------
// Function iterate()
// ----------------------------------------------------------------------------

template <typename TContainer, typename TFunctor, typename TIterTag, typename TParallelTag>
inline void iterate(TContainer & c, TFunctor f, Tag<TIterTag> const & iterTag, Tag<TParallelTag> const & /* tag */)
{
    typedef Tag<TIterTag> const                                     TIterSpec;
    typedef typename Iterator<TContainer, TIterSpec>::Type          TIter;

    for (TIter it = begin(c, iterTag); !atEnd(it, c); ++it)
        f(it);
}

// ----------------------------------------------------------------------------
// Function iterate(Parallel)
// ----------------------------------------------------------------------------

template <typename TContainer, typename TFunctor, typename TIterTag>
inline void iterate(TContainer & c, TFunctor f, Tag<TIterTag> const & iterTag, Parallel)
{
    typedef Tag<TIterTag> const                                     TIterSpec;
    typedef typename Position<TContainer>::Type                     TPos;
    typedef typename Iterator<TContainer, TIterSpec>::Type          TIter;

    Splitter<TPos> splitter(0, length(c), Parallel());

    SEQAN_OMP_PRAGMA(parallel for firstprivate(f))
    for (TPos i = 0; i < length(splitter); ++i)
    {
       TIter it = begin(c, iterTag) + splitter[i];
       TIter itEnd = begin(c, iterTag) + splitter[i + 1];

       for (; it != itEnd; ++it)
            f(it);
    }
}

// ----------------------------------------------------------------------------
// Function iterate()
// ----------------------------------------------------------------------------

template <typename TContainer, typename TFunctor>
inline void iterate(TContainer & c, TFunctor f)
{
    iterate(c, f, typename DefaultIteratorSpec<TContainer>::Type(), Serial());
}

// ============================================================================
// STL Wrappers
// ============================================================================

// ----------------------------------------------------------------------------
// Function forEach()
// ----------------------------------------------------------------------------

template <typename TContainer, typename TFunctor, typename TParallelTag>
inline TFunctor
forEach(TContainer & c, TFunctor f, Tag<TParallelTag> const & /* tag */)
{
    return std::for_each(begin(c, Standard()), end(c, Standard()), f);
}

template <typename TContainer, typename TFunctor, typename TParallelTag>
inline TFunctor
forEach(TContainer const & c, TFunctor f, Tag<TParallelTag> const & /* tag */)
{
    return std::for_each(begin(c, Standard()), end(c, Standard()), f);
}

// ----------------------------------------------------------------------------
// Function transform()
// ----------------------------------------------------------------------------

template <typename TTarget, typename TSource, typename TUnaryOperator, typename TParallelTag>
inline void transform(TTarget & target, TSource & source, TUnaryOperator o, Tag<TParallelTag> const & /* tag */)
{
    SEQAN_ASSERT_GEQ(length(target), length(source));
    std::transform(begin(source, Standard()), end(source, Standard()), begin(target, Standard()), o);
}

template <typename TContainer, typename TUnaryOperator, typename TParallelTag>
inline void transform(TContainer & c, TUnaryOperator o, Tag<TParallelTag> const & tag)
{
    transform(c, c, o, tag);
}

// ----------------------------------------------------------------------------
// Function generate()
// ----------------------------------------------------------------------------

template <typename TTarget, typename TGenerator, typename TParallelTag>
inline void generate(TTarget & target, TGenerator g, Tag<TParallelTag> const & /* tag */)
{
    std::generate(begin(target, Standard()), end(target, Standard()), g);
}

// ----------------------------------------------------------------------------
// Function iota()
// ----------------------------------------------------------------------------

#ifdef SEQAN_CXX11_STANDARD
template <typename TTarget, typename TValue, typename TParallelTag>
inline void iota(TTarget & target, TValue val, Tag<TParallelTag> const & /* tag */)
{
    std::iota(begin(target, Standard()), end(target, Standard()), val);
}
#endif

// ----------------------------------------------------------------------------
// Function count()
// ----------------------------------------------------------------------------

template <typename TContainer, typename TValue, typename TParallelTag>
inline typename Difference<TContainer>::Type
count(TContainer const & c, TValue const & v, Tag<TParallelTag> const & /* tag */)
{
    return std::count(begin(c, Standard()), end(c, Standard()), v);
}

// ----------------------------------------------------------------------------
// Function countIf()
// ----------------------------------------------------------------------------

template <typename TContainer, typename TUnaryPredicate, typename TParallelTag>
inline typename Difference<TContainer>::Type
countIf(TContainer const & c, TUnaryPredicate p, Tag<TParallelTag> const & /* tag */)
{
    return std::count_if(begin(c, Standard()), end(c, Standard()), p);
}

// ----------------------------------------------------------------------------
// Function maxElement()
// ----------------------------------------------------------------------------

template <typename TContainer, typename TUnaryPredicate, typename TParallelTag>
inline typename Reference<TContainer const>::Type
maxElement(TContainer const & c, TUnaryPredicate p, Tag<TParallelTag> const & /* tag */)
{
    SEQAN_ASSERT_NOT(empty(c));
    return value(std::max_element(begin(c, Standard()), end(c, Standard()), p));
}

template <typename TContainer, typename TParallelTag>
inline typename Reference<TContainer const>::Type
maxElement(TContainer const & c, Tag<TParallelTag> const & /* tag */)
{
    SEQAN_ASSERT_NOT(empty(c));
    return value(std::max_element(begin(c, Standard()), end(c, Standard())));
}

// ----------------------------------------------------------------------------
// Function minElement()
// ----------------------------------------------------------------------------

template <typename TContainer, typename TUnaryPredicate, typename TParallelTag>
inline typename Reference<TContainer const>::Type
minElement(TContainer const & c, TUnaryPredicate p, Tag<TParallelTag> const & /* tag */)
{
    SEQAN_ASSERT_NOT(empty(c));
    return value(std::min_element(begin(c, Standard()), end(c, Standard()), p));
}

template <typename TContainer, typename TParallelTag>
inline typename Reference<TContainer const>::Type
minElement(TContainer const & c, Tag<TParallelTag> const & /* tag */)
{
    SEQAN_ASSERT_NOT(empty(c));
    return value(std::min_element(begin(c, Standard()), end(c, Standard())));
}

// ----------------------------------------------------------------------------
// Function sort()
// ----------------------------------------------------------------------------

template <typename TContainer, typename TBinaryPredicate, typename TParallelTag>
inline void sort(TContainer & c, TBinaryPredicate p, Tag<TParallelTag> const & /* tag */)
{
    return std::sort(begin(c, Standard()), end(c, Standard()), p);
}

template <typename TContainer, typename TParallelTag>
inline void sort(TContainer & c, Tag<TParallelTag> const & /* tag */)
{
    return std::sort(begin(c, Standard()), end(c, Standard()));
}

// ----------------------------------------------------------------------------
// Function stableSort()
// ----------------------------------------------------------------------------

template <typename TContainer, typename TBinaryPredicate, typename TParallelTag>
inline void stableSort(TContainer & c, TBinaryPredicate p, Tag<TParallelTag> const & /* tag */)
{
    return std::stable_sort(begin(c, Standard()), end(c, Standard()), p);
}

template <typename TContainer, typename TParallelTag>
inline void stableSort(TContainer & c, Tag<TParallelTag> const & /* tag */)
{
    return std::stable_sort(begin(c, Standard()), end(c, Standard()));
}

// ============================================================================
// MCSTL Wrappers
// ============================================================================

// use MCSTL which is part of the GCC since version 4.3
#if defined(_OPENMP) && defined(PLATFORM_GCC) && __GNUC__ >= 4 && __GNUC_MINOR__ >= 3

// ----------------------------------------------------------------------------
// Function forEach(Parallel)
// ----------------------------------------------------------------------------

template <typename TContainer, typename TFunctor>
inline TFunctor forEach(TContainer & c, TFunctor f, Parallel)
{
    return __gnu_parallel::for_each(begin(c, Standard()), end(c, Standard()), f);
}

template <typename TContainer, typename TFunctor>
inline TFunctor forEach(TContainer const & c, TFunctor f, Parallel)
{
    return __gnu_parallel::for_each(begin(c, Standard()), end(c, Standard()), f);
}

// ----------------------------------------------------------------------------
// Function transform(Parallel)
// ----------------------------------------------------------------------------

template <typename TTarget, typename TSource, typename TUnaryOperator>
inline void transform(TTarget & target, TSource & source, TUnaryOperator o, Parallel)
{
    SEQAN_ASSERT_GEQ(length(target), length(source));
    __gnu_parallel::transform(begin(source, Standard()), end(source, Standard()), begin(target, Standard()), o);
}

// ----------------------------------------------------------------------------
// Function generate(Parallel)
// ----------------------------------------------------------------------------

template <typename TTarget, typename TGenerator>
inline void generate(TTarget & target, TGenerator g, Parallel)
{
    __gnu_parallel::generate(begin(target, Standard()), end(target, Standard()), g);
}

// ----------------------------------------------------------------------------
// Function count(Parallel)
// ----------------------------------------------------------------------------

template <typename TContainer, typename TValue>
inline typename Difference<TContainer>::Type
count(TContainer const & c, TValue const & v, Parallel)
{
    return __gnu_parallel::count(begin(c, Standard()), end(c, Standard()), v);
}

// ----------------------------------------------------------------------------
// Function countIf(Parallel)
// ----------------------------------------------------------------------------

template <typename TContainer, typename TUnaryPredicate>
inline typename Difference<TContainer>::Type
countIf(TContainer const & c, TUnaryPredicate p, Parallel)
{
    return __gnu_parallel::count_if(begin(c, Standard()), end(c, Standard()), p);
}

// ----------------------------------------------------------------------------
// Function maxElement(Parallel)
// ----------------------------------------------------------------------------

template <typename TContainer, typename TUnaryPredicate>
inline typename Reference<TContainer const>::Type
maxElement(TContainer const & c, TUnaryPredicate p, Parallel)
{
    SEQAN_ASSERT_NOT(empty(c));
    return value(__gnu_parallel::max_element(begin(c, Standard()), end(c, Standard()), p));
}

template <typename TContainer>
inline typename Reference<TContainer const>::Type
maxElement(TContainer const & c, Parallel)
{
    SEQAN_ASSERT_NOT(empty(c));
    return value(__gnu_parallel::max_element(begin(c, Standard()), end(c, Standard())));
}

// ----------------------------------------------------------------------------
// Function minElement(Parallel)
// ----------------------------------------------------------------------------

template <typename TContainer, typename TUnaryPredicate>
inline typename Reference<TContainer const>::Type
minElement(TContainer const & c, TUnaryPredicate p, Parallel)
{
    SEQAN_ASSERT_NOT(empty(c));
    return value(__gnu_parallel::min_element(begin(c, Standard()), end(c, Standard()), p));
}

template <typename TContainer>
inline typename Reference<TContainer const>::Type
minElement(TContainer const & c, Parallel)
{
    SEQAN_ASSERT_NOT(empty(c));
    return value(__gnu_parallel::min_element(begin(c, Standard()), end(c, Standard())));
}

// ----------------------------------------------------------------------------
// Function sort(Parallel)
// ----------------------------------------------------------------------------

template <typename TContainer, typename TBinaryPredicate>
inline void sort(TContainer & c, TBinaryPredicate p, Parallel)
{
    return __gnu_parallel::sort(begin(c, Standard()), end(c, Standard()), p);
}

template <typename TContainer>
inline void sort(TContainer & c, Parallel)
{
    return __gnu_parallel::sort(begin(c, Standard()), end(c, Standard()));
}

// ----------------------------------------------------------------------------
// Function stableSort(Parallel)
// ----------------------------------------------------------------------------

template <typename TContainer, typename TBinaryPredicate>
inline void stableSort(TContainer & c, TBinaryPredicate p, Parallel)
{
    return __gnu_parallel::stable_sort(begin(c, Standard()), end(c, Standard()), p);
}

template <typename TContainer>
inline void stableSort(TContainer & c, Parallel)
{
    return __gnu_parallel::stable_sort(begin(c, Standard()), end(c, Standard()));
}

#endif  // #ifdef PLATFORM_GCC

// ============================================================================
// Shortcuts for STL Wrappers
// ============================================================================

// ----------------------------------------------------------------------------
// Function forEach()
// ----------------------------------------------------------------------------

template <typename TContainer, typename TFunctor>
inline TFunctor forEach(TContainer & c, TFunctor f)
{
    return forEach(c, f, Serial());
}

template <typename TContainer, typename TFunctor>
inline TFunctor forEach(TContainer const & c, TFunctor f)
{
    return forEach(c, f, Serial());
}

// ----------------------------------------------------------------------------
// Function transform()
// ----------------------------------------------------------------------------

template <typename TTarget, typename TSource, typename TUnaryOperator>
inline void transform(TTarget & target, TSource & source, TUnaryOperator o)
{
    transform(target, source, o, Serial());
}

template <typename TContainer, typename TUnaryOperator>
inline void transform(TContainer & c, TUnaryOperator o)
{
    transform(c, o, Serial());
}

// ----------------------------------------------------------------------------
// Function generate()
// ----------------------------------------------------------------------------

template <typename TTarget, typename TGenerator>
inline void generate(TTarget & target, TGenerator g)
{
    generate(target, g, Serial());
}

// ----------------------------------------------------------------------------
// Function iota()
// ----------------------------------------------------------------------------

#ifdef SEQAN_CXX11_STANDARD
template <typename TTarget, typename TValue>
inline void iota(TTarget & target, TValue val)
{
    iota(target, val, Serial());
}
#endif

// ----------------------------------------------------------------------------
// Function count()
// ----------------------------------------------------------------------------

template <typename TContainer, typename TValue>
inline typename Difference<TContainer>::Type
count(TContainer const & c, TValue const & v)
{
    return count(c, v, Serial());
}

// ----------------------------------------------------------------------------
// Function countIf()
// ----------------------------------------------------------------------------

template <typename TContainer, typename TUnaryPredicate>
inline typename Difference<TContainer>::Type
countIf(TContainer const & c, TUnaryPredicate p)
{
    return countIf(c, p, Serial());
}

// ----------------------------------------------------------------------------
// Function maxElement()
// ----------------------------------------------------------------------------

template <typename TContainer, typename TUnaryPredicate>
inline typename Reference<TContainer const>::Type
maxElement(TContainer const & c, TUnaryPredicate p)
{
    return maxElement(c, p, Serial());
}

template <typename TContainer>
inline typename Reference<TContainer const>::Type
maxElement(TContainer const & c)
{
    return maxElement(c, Serial());
}

// ----------------------------------------------------------------------------
// Function minElement()
// ----------------------------------------------------------------------------

template <typename TContainer, typename TUnaryPredicate>
inline typename Reference<TContainer const>::Type
minElement(TContainer const & c, TUnaryPredicate p)
{
    return minElement(c, p, Serial());
}

template <typename TContainer>
inline typename Reference<TContainer const>::Type
minElement(TContainer const & c)
{
    return minElement(c, Serial());
}

// ----------------------------------------------------------------------------
// Function sort()
// ----------------------------------------------------------------------------

template <typename TContainer, typename TBinaryPredicate>
inline void sort(TContainer & c, TBinaryPredicate p)
{
    sort(c, p, Serial());
}

template <typename TContainer>
inline void sort(TContainer & c)
{
    sort(c, Serial());
}

// ----------------------------------------------------------------------------
// Function stableSort()
// ----------------------------------------------------------------------------

template <typename TContainer, typename TBinaryPredicate>
inline void stableSort(TContainer & c, TBinaryPredicate p)
{
    stableSort(c, p, Serial());
}

template <typename TContainer>
inline void stableSort(TContainer & c)
{
    stableSort(c, Serial());
}

}  // namespace seqan

#endif  // #ifndef SEQAN_PARALLEL_PARALLEL_ALGORITHMS_H_
