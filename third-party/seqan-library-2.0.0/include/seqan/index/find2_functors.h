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
//     * Neither the name of NVIDIA Corporation nor the names of
//       its contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL NVIDIA CORPORATION BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
// LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
// OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
// DAMAGE.
//
// ==========================================================================
// Author: Enrico Siragusa <enrico.siragusa@fu-berlin.de>
// ==========================================================================
// Some useful functors to be used along with the Finder class.
// ==========================================================================

#ifndef SEQAN_INDEX_FIND_FUNCTORS_H_
#define SEQAN_INDEX_FIND_FUNCTORS_H_

namespace seqan {

// ============================================================================
// Tags
// ============================================================================

// ----------------------------------------------------------------------------
// Tag Counts_
// ----------------------------------------------------------------------------

struct Counts_;

// ============================================================================
// Classes
// ============================================================================

// ----------------------------------------------------------------------------
// Class OccurrencesCounter_
// ----------------------------------------------------------------------------

template <typename TIndex, typename TSpec = typename ExecSpec<TIndex>::Type>
struct OccurrencesCounter_
{
    typename Member<OccurrencesCounter_, Counts_>::Type    counts;

    OccurrencesCounter_() {}

    template <typename TPattern>
    OccurrencesCounter_(TPattern const & pattern)
    {
        _init(*this, pattern);
    }

    template <typename TFinder>
    inline SEQAN_HOST_DEVICE void
    operator() (TFinder const & finder)
    {
        counts[getThreadId()] += countOccurrences(_textIterator(finder));
    }
};

// ----------------------------------------------------------------------------
// Member Counts_; Count
// ----------------------------------------------------------------------------

template <typename TIndex, typename TSpec>
struct Member<OccurrencesCounter_<TIndex, TSpec>, Counts_>
{
    typedef String<typename Size<TIndex>::Type> Type;
};

#ifdef PLATFORM_CUDA
template <typename TIndex, typename TSpec>
struct Member<OccurrencesCounter_<TIndex, Device<TSpec> >, Counts_>
{
    typedef thrust::device_vector<typename Size<TIndex>::Type>  Type;
};
#endif

template <typename TIndex, typename TSpec>
struct Member<OccurrencesCounter_<TIndex, View<Device<TSpec> > >, Counts_>
{
    typedef typename View<typename Member<OccurrencesCounter_<TIndex, Device<TSpec> >, Counts_>::Type>::Type Type;
};

// ----------------------------------------------------------------------------
// Metafunction View
// ----------------------------------------------------------------------------

template <typename TIndex, typename TSpec>
struct View<OccurrencesCounter_<TIndex, TSpec> >
{
    typedef OccurrencesCounter_<TIndex, View<TSpec> >  Type;
};

// ----------------------------------------------------------------------------
// Metafunction Device
// ----------------------------------------------------------------------------

template <typename TIndex, typename TSpec>
struct Device<OccurrencesCounter_<TIndex, TSpec> >
{
    typedef OccurrencesCounter_<TIndex, Device<TSpec> >  Type;
};

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function _init()
// ----------------------------------------------------------------------------

template <typename TIndex, typename TSpec, typename TPattern>
inline void
_init(OccurrencesCounter_<TIndex, TSpec> & counter, TPattern const & /* pattern */)
{
    resize(counter.counts, omp_get_max_threads(), 0, Exact());
}

template <typename TIndex, typename TSpec, typename TPattern>
inline void
_init(OccurrencesCounter_<TIndex, Device<TSpec> > & counter, TPattern const & pattern)
{
    resize(counter.counts, length(needle(pattern)), Exact());
}

// ----------------------------------------------------------------------------
// Function view()
// ----------------------------------------------------------------------------

template <typename TIndex, typename TSpec>
typename View<OccurrencesCounter_<TIndex, TSpec> >::Type
view(OccurrencesCounter_<TIndex, TSpec> & counter)
{
    typename View<OccurrencesCounter_<TIndex, TSpec> >::Type counterView;

    counterView.counts = view(counter.counts);

    return counterView;
}

// ----------------------------------------------------------------------------
// Function _getCount()
// ----------------------------------------------------------------------------

template <typename TIndex, typename TSpec>
inline typename Size<TIndex>::Type
_getCount(OccurrencesCounter_<TIndex, TSpec> & counter)
{
    return sum(counter.counts);
}

#ifdef PLATFORM_CUDA
template <typename TIndex, typename TSpec>
inline typename Size<TIndex>::Type
_getCount(OccurrencesCounter_<TIndex, Device<TSpec> > & counter)
{
    return thrust::reduce(begin(counter.counts, Standard()), end(counter.counts, Standard()));
}
#endif

// --------------------------------------------------------------------------
// Function countOccurrences()
// --------------------------------------------------------------------------
// Count the occurrences of a set of needles in a indexed haystack.

template <typename TText, typename TSpec, typename TNeedle, typename TSSetSpec>
typename Size<Index<TText, TSpec> >::Type
countOccurrences(Index<TText, TSpec> & index, StringSet<TNeedle, TSSetSpec> & needles)
{
    typedef Index<TText, TSpec>                         TIndex;
    typedef StringSet<TNeedle, TSSetSpec>               TNeedles;
    typedef Multiple<FinderSTree>                       TAlgorithmSpec;
    typedef Pattern<TNeedles, TAlgorithmSpec>           TPattern;
    typedef Finder_<TIndex, TPattern, TAlgorithmSpec>   TFinder;
    typedef OccurrencesCounter_<TIndex>                  TCounter;

    // Instantiate a finder object holding the context of the search algorithm.
    TFinder finder(index);

    // Instantiate a pattern object holding the needles.
    TPattern pattern(needles);

    // Instantiate a functor object counting the number of occurrences.
    TCounter counter(pattern);

    // Find all needles in haystack and call counter() on match.
    _find(finder, pattern, counter);

    // Return the number of occurrences.
    return _getCount(counter);
}

}

#endif  // #ifndef SEQAN_INDEX_FIND_FUNCTORS_H_
