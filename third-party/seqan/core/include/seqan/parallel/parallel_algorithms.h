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
// Basic parallel algorithms.
// ==========================================================================

#ifndef SEQAN_PARALLEL_PARALLEL_ALGORITHMS_H_
#define SEQAN_PARALLEL_PARALLEL_ALGORITHMS_H_

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

/**
.Function.arrayFill
..signature:arrayFill(begin, end, value, parallelTag)
..param.parallelTag:Tag to enable/disable parallelism.
...type:Tag.Serial
...type:Tag.Parallel
..include:seqan/basic.h
*/

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

/**
.Function.sum
..cat:Miscellaneous
..summary:Returns the sum of all elements in a sequence.
..signature:sum(seq[, parallelTag])
..param.seq:A sequence of elements that should be summed.
...remarks:The sequence alphabet must support the $operator+$ and conversion from zero.
..param.parallelTag:Tag to enable/disable parallelism.
...default:Tag.Serial
...type:Tag.Serial
...type:Tag.Parallel
..returns:The sum of all contained elements.
..see:Function.partialSum
..include:seqan/parallel.h
*/

template <typename TSequence>
inline typename Value<TSequence>::Type
sum(TSequence const &seq, Serial)
{
    register typename Iterator<TSequence const>::Type it = begin(seq, Standard());
    register typename Iterator<TSequence const>::Type itEnd = end(seq, Standard());
    register typename Value<TSequence>::Type sum = 0;
    for (; it != itEnd; ++it)
        sum += *it;
    return sum;
}

template <typename TSequence, typename TParallelTag>
inline typename Value<TSequence>::Type
sum(TSequence const &seq, Tag<TParallelTag> parallelTag)
{
    Splitter<typename Size<TSequence>::Type> splitter(0, length(seq), parallelTag);

    register typename Value<TSequence>::Type threadSum = 0;
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

/**
.Function.partialSum
..cat:Miscellaneous
..summary:Computes the partial sum of a sequence.
..signature:partialSum(target, source[, parallelTag])
..param.source:A sequence of elements that should be partially summed.
...remarks:The sequence alphabet must support the $operator+$ and conversion from zero.
..param.target:The resulting partial sum. This sequence will have the same length as $source$ and contains
at position $i$ the sum of elements $source[0]$, $source[1]$, ..., $source[i]$.
..param.parallelTag:Tag to enable/disable parallelism.
...default:@Tag.Serial@
...type:Tag.Serial
...type:Tag.Parallel
..returns:The sum of all elements in $source$.
...remarks:The returned value equals the last value in target.
..see:Function.sum
..include:seqan/parallel.h
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
        register TConstIterator it = begin(source, Standard()) + splitter[job];
        register TConstIterator itEnd = begin(source, Standard()) + splitter[job + 1];
        register TIterator dstIt = begin(target, Standard()) + splitter[job];
        register TValue sum = localSums[job];
        for (; it != itEnd; ++it, ++dstIt)
        {
            sum += *it;
            *dstIt = sum;
        }
        localSums[job] = sum;
    }
    
    return back(localSums);
}

template <typename TTarget, typename TSource>
inline typename Value<TSource>::Type
partialSum(TTarget &target, TSource const &source)
{
    return partialSum(target, source, Serial());
}

}  // namespace seqan

#endif  // #ifndef SEQAN_PARALLEL_PARALLEL_ALGORITHMS_H_
