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
// Author: Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>
// Author: David Weese <david.weese@fu-berlin.de>
// ==========================================================================
// Utility macros for parallelism.
// ==========================================================================

#ifndef SEQAN_PARALLEL_PARALLEL_SPLITTING_H_
#define SEQAN_PARALLEL_PARALLEL_SPLITTING_H_

namespace seqan {

struct Equidistant_;
typedef Tag<Equidistant_> Equidistant;

/**
.Class.Splitter
..cat:Parallelism
..summary:Splits an interval into subintervals.
..description:This class divides an interval into the disjoint union of subintervals and enumerates its boundaries.
It can be used to parallelize large for-loops that iterate over a contiguous range of elements.
The interval and the number of subintervals can be set in the constructor @Memfunc.Splitter#Splitter@.
@Function.length@ and @Function.resize@ can be used to retrieve or change the number of subintervals later.
In contrast to other containers the Splitter allows to access one more element than its length would imply to allow to retrieve the right boundary of each subinterval (see example code below).

..signature:Splitter<TValue, TSpec>
..param.TValue:Type of the interval boundaries.
..param.TSpec:Tag to select the way the values are sampled.
...default:Equidistant Splitter
..example:
...text:Simple example for equidistant (default) splitting.
...file:demos/parallel/splitter_example.cpp
...output:
[10,14)
[14,17)
[17,20)
..include:seqan/parallel.h

.Memfunc.Splitter#Splitter
..summary:Constructor
..signature:Splitter(beginPos, endPos[, subintervalCount])
..signature:Splitter(beginPos, endPos, parallelTag)
..param.beginPos:Left interval boundary.
..param.endPos:Right interval boundary.
..param.subintervalCount:Number of subintervals.
...default:The minimum of interval size and the number of available threads (returned by $omp_get_max_threads()$).
...remarks:@Function.length@ and @Function.resize@ can be used to retrieve or change the number of subintervals later.
..param.parallelTag:Tag to generically enable/disable parallelism. If its type is @Tag.Parallel@, the default number of subintervals is used.
If it is @Tag.Serial@, only one subinterval is used.
...type:Tag.Serial
...type:Tag.Parallel
...remarks:This tag should be used to write generic parallel algorithms and to switch between parallel and serial variants.
..class:Class.Splitter
 */

/**
.Spec.Equidistant Splitter
..cat:Parallelism
..general:Class.Splitter
..summary:Splits an interval into equal-sized subintervals.
..signature:Splitter<TValue[, Equidistant]>
..description:This @Class.Splitter@ specialization divides an interval into subintervals of (almost) equal length, i.e.
two subintervals differ by at most 1 in size.
..param.TValue:Type of the interval boundaries.
..example:
...text:Simple example for equidistant splitting.
...file:demos/parallel/splitter_example.cpp
...output:
[10,14)
[14,17)
[17,20)
..include:seqan/parallel.h */
 
template <typename TValue, typename TSpec = Equidistant>
class Splitter
{
public:
    typedef typename Size<Splitter>::Type TSize;

    TValue beginPos;
    TSize subintervalCount;
    TSize blockLength;
    TSize rest;

    Splitter(TValue beginPos_, TValue endPos):
        beginPos(beginPos_)
    {
        // we choose the counts automatically and don't want to have empty jobs
        _resize(*this, endPos - beginPos, _min((TSize)(endPos - beginPos), (TSize)omp_get_max_threads()));
    }

    Splitter(TValue beginPos_, TValue endPos, Parallel):
        beginPos(beginPos_)
    {
        // we choose the counts automatically and don't want to have empty jobs
        _resize(*this, endPos - beginPos, _min((TSize)(endPos - beginPos), (TSize)omp_get_max_threads()));
    }

    Splitter(TValue beginPos_, TValue endPos, Serial):
        beginPos(beginPos_)
    {
        // we produce at most 1 job (or none if interval is empty)
        _resize(*this, endPos - beginPos, _min((TSize)(endPos - beginPos), (TSize)1));
    }
    
    Splitter(TValue beginPos_, TValue endPos, TSize subintervalCount):
        beginPos(beginPos_)
    {
        _resize(*this, endPos - beginPos, subintervalCount);
    }

    TValue operator[] (TSize i) const
    {
        SEQAN_ASSERT_LEQ_MSG(i, subintervalCount, "Trying to access an element behind the last one!");
        return beginPos + blockLength * i + std::min(i, rest);
    }
};


///.Metafunction.Size.param.T.type:Class.Splitter
template <typename TValue, typename TSpec>
struct Size<Splitter<TValue, TSpec> >
{
    typedef typename MakeUnsigned<typename Difference<TValue>::Type>::Type Type;
};

///.Metafunction.Value.param.T.type:Class.Splitter
template <typename TValue, typename TSpec>
struct Value<Splitter<TValue, TSpec> >
{
    typedef TValue Type;
};

///.Function.length.param.object.type:Class.Splitter
///.Function.length.class:Class.Splitter
template <typename TValue, typename TSpec>
inline typename Size<Splitter<TValue, TSpec> >::Type
length(Splitter<TValue, TSpec> const &splitter)
{
    return splitter.subintervalCount;
}

template <typename TValue, typename TSpec, typename TSize1, typename TSize2>
inline void
_resize(Splitter<TValue, TSpec> &splitter, TSize1 intervalLen, TSize2 newCount)
{
    if (newCount != 0)
    {
        splitter.blockLength = intervalLen / newCount;
        splitter.rest = intervalLen % newCount;
    }
    else
    {
        splitter.blockLength = 0;
        splitter.rest = intervalLen;
    }
    splitter.subintervalCount = newCount;
}

///.Function.resize.param.object.type:Class.Splitter
///.Function.resize.class:Class.Splitter
template <typename TValue, typename TSpec, typename TSize>
inline typename Size<Splitter<TValue, TSpec> >::Type
resize(Splitter<TValue, TSpec> &splitter, TSize newCount)
{
    _resize(splitter, splitter.blockLength * splitter.subintervalCount + splitter.rest, newCount);
    return newCount;
}

/**
.Function.computeSplitters
..cat:Parallelism
..summary:Compute splitters for a sequence of objects.
..signature:computeSplitters(splitters, size, count)
..param.splitters:Resulting splitters, will be resized to contain $count + 1$ elements.
...type:Spec.Alloc String
..param.size:The number of objects to split.
..param.count:The number of chunks.
..remarks:The first $count - 1$ chunks will have the size $ceil(size / count)$, the last chunk will contain the rest.
..example.text:Most simple case for splitting.
..example.code:String<unsigned> splitters;
computeSplitters(splitters, 10, 5);
// splitters == {0, 5, 10}
..example.text:In this case, the last chunks will stay empty.
..example.code:computeSplitters(splitters, 3, 5);
// splitters == {0, 1, 2, 3, 3, 3}
..include:seqan/parallel.h
 */

template <typename TPosString, typename TSize, typename TCount>
void computeSplitters(TPosString & splitters, TSize size, TCount count)
{
    typedef typename Value<TPosString>::Type TPos;

    SEQAN_ASSERT_GEQ(count, (TCount)0);

    resize(splitters, count + 1);

    TSize blockLength = (size / count) + 1;
    TCount rest = size % count;
    TPos pos = 0;

    TCount i = 0;

    // the first (size % count) many blocks have length (size / count) + 1
    for (; i < rest; ++i, pos += blockLength)
        splitters[i] = pos;

    // the remaining blocks have length (size / count)
    --blockLength;
    for (; i <= count; ++i, pos += blockLength)
        splitters[i] = pos;

    SEQAN_ASSERT_EQ(back(splitters), static_cast<TPos>(size));
}

}  // namespace seqan

#endif  // SEQAN_PARALLEL_PARALLEL_SPLITTING_H_
