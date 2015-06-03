// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
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

#ifndef SEQAN_INDEX_FIND_INDEX_MULTIPLE_H
#define SEQAN_INDEX_FIND_INDEX_MULTIPLE_H

#ifdef PLATFORM_CUDA
#include <thrust/sort.h>
#include <thrust/reduce.h>
#endif

namespace seqan {

// ============================================================================
// Tags
// ============================================================================

// ----------------------------------------------------------------------------
// Tag Multiple
// ----------------------------------------------------------------------------

template <typename TSpec>
struct Multiple;

// ----------------------------------------------------------------------------
// Finder Member Tags
// ----------------------------------------------------------------------------

struct Factory_;

// ----------------------------------------------------------------------------
// Pattern Member Tags
// ----------------------------------------------------------------------------

struct Needles_;
struct Hashes_;
struct Permutation_;

// ============================================================================
// Classes
// ============================================================================

// ----------------------------------------------------------------------------
// Class Finder; Multiple
// ----------------------------------------------------------------------------

template <typename TText, typename TPattern, typename TSpec>
struct Finder_<TText, TPattern, Multiple<TSpec> >
{
    typename Member<Finder_, Factory_>::Type    _factory;

    Finder_() {}

    Finder_(TText & text) :
        _factory(text)
    {}
};

// ----------------------------------------------------------------------------
// Class Finder; Index, Multiple
// ----------------------------------------------------------------------------

template <typename TText, typename TIndexSpec, typename TPattern, typename TSpec>
struct Finder_<Index<TText, TIndexSpec>, TPattern, Multiple<TSpec> >
{
    typedef Index<TText, TIndexSpec>            TIndex;

    typename Member<Finder_, Factory_>::Type    _factory;
    typename Score_<TSpec>::Type                _scoreThreshold;

    Finder_() :
        _scoreThreshold()
    {}

    Finder_(TIndex & index) :
        _factory(index),
        _scoreThreshold()
    {}
};

// ----------------------------------------------------------------------------
// Class Pattern; Multiple
// ----------------------------------------------------------------------------

template <typename TNeedles, typename TSpec>
class Pattern<TNeedles, Multiple<TSpec> >
{
public:
    typename Member<Pattern, Needles_>::Type        data_host;
    typename Member<Pattern, Hashes_>::Type         _hashes;
    typename Member<Pattern, Permutation_>::Type    _permutation;

    Pattern() {}

    Pattern(TNeedles const & needles) :
        data_host(needles)
    {
        _preprocess(*this);
    }

//    Pattern(TNeedles const & needles)
//    {
//        setHost(*this, needles);
//    }
};

// ----------------------------------------------------------------------------
// Class FinderContext_
// ----------------------------------------------------------------------------

template <typename TText, typename TPattern, typename TSpec, typename TDelegate>
struct FinderContext_;

template <typename TText, typename TPattern, typename TSpec, typename TDelegate>
struct FinderContext_<TText, TPattern, Multiple<TSpec>, TDelegate>
{
    typedef typename Needle<TPattern>::Type                 TNeedles;
    typedef typename Value<TPattern>::Type                  TNeedle;
    typedef Finder_<TText, TNeedle, TSpec>                  TBaseFinder;
    typedef Finder_<TText, TPattern, Multiple<TSpec> >      TFinder;
    typedef typename Iterator<TNeedles, Standard>::Type     TPatternIterator;

    TFinder &      finder;
    TDelegate &    delegate;
    TBaseFinder    baseFinder;
//    TPatternIterator _patternIt;
    unsigned        _patternIt;

    explicit SEQAN_HOST_DEVICE
    FinderContext_(TFinder & finder, TDelegate & delegate) :
        finder(finder),
        delegate(delegate),
        baseFinder(getObject(finder._factory, getThreadId()))
    {}

    // NOTE(esiragusa): This is called on firstprivate.
    FinderContext_(FinderContext_ & ctx) :
        finder(ctx.finder),
        delegate(ctx.delegate),
        baseFinder(getObject(finder._factory, getThreadId()))
    {}

    template <typename TOther>
    SEQAN_HOST_DEVICE inline void
    operator()(TOther & /* other */)
    {
        delegate(*this);
    }
};

// ============================================================================
// Metafunctions
// ============================================================================

// ----------------------------------------------------------------------------
// Metafunction Score_; Multiple
// ----------------------------------------------------------------------------

template <typename TSpec>
struct Score_<Multiple<TSpec> >
{
    typedef typename Score_<TSpec>::Type    Type;
};

// ----------------------------------------------------------------------------
// Metafunction Host                                                  [Pattern]
// ----------------------------------------------------------------------------

template <typename TNeedles, typename TSpec>
struct Host<Pattern<TNeedles, Multiple<TSpec> > >
{
    typedef TNeedles Type;
};

// ----------------------------------------------------------------------------
// Metafunction PatternShape_                                         [Pattern]
// ----------------------------------------------------------------------------
// TODO(esiragusa): Automatically select a good shape.

template <typename TPattern>
struct PatternShape_
{
    typedef typename Host<TPattern>::Type                               TNeedles_;
    typedef typename Value<TNeedles_>::Type                             TNeedle_;
    typedef typename Value<TNeedle_>::Type                              TAlphabet_;
    typedef UngappedShape<Max<MinLength<TNeedles_>::VALUE, 1>::VALUE>   TShapeSpec_;

    typedef Shape<TAlphabet_, TShapeSpec_>                              Type;
};

// ----------------------------------------------------------------------------
// Member Needles_                                                    [Pattern]
// ----------------------------------------------------------------------------

template <typename TNeedles, typename TSpec>
struct Member<Pattern<TNeedles, Multiple<TSpec> >, Needles_>
{
    typedef typename IfView<TNeedles, TNeedles, Holder<TNeedles> >::Type    Type;
};

// ----------------------------------------------------------------------------
// Member Hashes_                                                     [Pattern]
// ----------------------------------------------------------------------------

template <typename TNeedles, typename TSpec>
struct Member<Pattern<TNeedles, Multiple<TSpec> >, Hashes_>
{
    typedef Nothing Type;
};

#ifdef PLATFORM_CUDA
template <typename TValue, typename TAlloc, typename TSSetSpec, typename TSpec>
struct Member<Pattern<StringSet<thrust::device_vector<TValue, TAlloc>, TSSetSpec>, Multiple<TSpec> >, Hashes_>
{
    typedef thrust::device_vector<TValue, TAlloc>   TNeedle_;
    typedef StringSet<TNeedle_, TSSetSpec>          TNeedles_;
    typedef Pattern<TNeedles_, Multiple<TSpec> >    TPattern_;
    typedef typename PatternShape_<TPattern_>::Type TShape_;
    typedef typename Value<TShape_>::Type           THash_;

    typedef thrust::device_vector<THash_>           Type;
};
#endif

template <typename TNeedle, typename TViewSpec, typename TSSetSpec, typename TSpec>
struct Member<Pattern<StringSet<ContainerView<TNeedle, TViewSpec>, TSSetSpec>, Multiple<TSpec> >, Hashes_>
{
    typedef Pattern<StringSet<TNeedle, TSSetSpec>, Multiple<TSpec> >        TPattern_;
    typedef typename View<typename Member<TPattern_, Hashes_>::Type>::Type  Type;
};

// ----------------------------------------------------------------------------
// Member Permutation_                                                [Pattern]
// ----------------------------------------------------------------------------

template <typename TNeedles, typename TSpec>
struct Member<Pattern<TNeedles, Multiple<TSpec> >, Permutation_>
{
    typedef Nothing Type;
};

#ifdef PLATFORM_CUDA
template <typename TValue, typename TAlloc, typename TSSetSpec, typename TSpec>
struct Member<Pattern<StringSet<thrust::device_vector<TValue, TAlloc>, TSSetSpec>, Multiple<TSpec> >, Permutation_>
{
    typedef thrust::device_vector<TValue, TAlloc>   TNeedle_;
    typedef StringSet<TNeedle_, TSSetSpec>          TNeedles_;
    typedef typename Size<TNeedles_>::Type          TSize_;

    typedef thrust::device_vector<TSize_>           Type;
};
#endif

template <typename TNeedle, typename TViewSpec, typename TSSetSpec, typename TSpec>
struct Member<Pattern<StringSet<ContainerView<TNeedle, TViewSpec>, TSSetSpec>, Multiple<TSpec> >, Permutation_>
{
    typedef Pattern<StringSet<TNeedle, TSSetSpec>, Multiple<TSpec> >            TPattern_;
    typedef typename View<typename Member<TPattern_, Permutation_>::Type>::Type Type;
};

// ----------------------------------------------------------------------------
// Member Factory_                                                     [Finder]
// ----------------------------------------------------------------------------

template <typename TText, typename TPattern, typename TSpec>
struct Member<Finder_<TText, TPattern, Multiple<TSpec> >, Factory_>
{
    typedef Factory<typename TextIterator_<TText, TPattern, TSpec>::Type>   Type;
};

// ----------------------------------------------------------------------------
// Metafunction Delegated                                              [Finder]
// ----------------------------------------------------------------------------

//template <typename TText, typename TPattern, typename TSpec>
//struct Delegated<Finder_<TText, TPattern, Multiple<TSpec> > >
//{
//    typedef Proxy<Finder_<TText, TPattern, Multiple<TSpec> > > Type;
//};

// ============================================================================
// Kernels
// ============================================================================

// ----------------------------------------------------------------------------
// Kernel _preprocessKernel()                                         [Pattern]
// ----------------------------------------------------------------------------

#ifdef PLATFORM_CUDA
template <typename TNeedles, typename TSpec>
SEQAN_GLOBAL void
_preprocessKernel(Pattern<TNeedles, Multiple<TSpec> > pattern)
{
    typedef Pattern<TNeedles, Multiple<TSpec> >     TPattern;
    typedef typename PatternShape_<TPattern>::Type  TShape;

    unsigned threadId = getThreadId();

    // Return silently if there is no job left.
    if (threadId >= length(pattern.data_host)) return;

    // Compute the hash of a needle.
    TShape shape;
    SEQAN_ASSERT_LEQ(weight(shape), length(pattern.data_host[threadId]));
    pattern._hashes[threadId] = hash(shape, begin(pattern.data_host[threadId], Standard()));

    // Fill with the identity permutation.
    pattern._permutation[threadId] = threadId;
}
#endif

// ----------------------------------------------------------------------------
// Kernel _findKernel()                                                [Finder]
// ----------------------------------------------------------------------------

#ifdef PLATFORM_CUDA
template <typename TText, typename TPattern, typename TSpec, typename TDelegate>
SEQAN_GLOBAL void
_findKernel(Finder_<TText, TPattern, Multiple<TSpec> > finder, TPattern pattern, TDelegate delegate)
{
    typedef FinderContext_<TText, TPattern, Multiple<TSpec>, TDelegate> TFinderContext;

    unsigned threadId = getThreadId();

    // Return if there is no job left.
    if (threadId >= length(pattern.data_host)) return;

    // Instantiate a thread context.
    TFinderContext ctx(finder, delegate);

    // Get the sorted needle id.
    ctx._patternIt = pattern._permutation[threadId];

    // Find a single needle.
    clear(ctx.baseFinder);
    _setScoreThreshold(ctx.baseFinder, _getScoreThreshold(finder));
    _find(ctx.baseFinder, pattern.data_host[ctx._patternIt], ctx);
}
#endif

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function _preprocess()                                   [Pattern; ExecHost]
// ----------------------------------------------------------------------------

template <typename TNeedles, typename TSpec>
inline void
_preprocess(Pattern<TNeedles, Multiple<TSpec> > & /* pattern */, ExecHost const & /* tag */) {}

// ----------------------------------------------------------------------------
// Function _preprocess()                                 [Pattern; ExecDevice]
// ----------------------------------------------------------------------------

#ifdef PLATFORM_CUDA
template <typename TNeedles, typename TSpec>
inline void
_preprocess(Pattern<TNeedles, Multiple<TSpec> > & pattern, ExecDevice const & /* tag */)
{
    typedef Pattern<TNeedles, Multiple<TSpec> > TPattern;
    typedef typename Size<TNeedles>::Type       TSize;

    TSize needlesCount = length(needle(pattern));

    resize(pattern._hashes, needlesCount, Exact());
    resize(pattern._permutation, needlesCount, Exact());

    // Compute grid size.
    unsigned ctaSize = CtaSize<TPattern>::VALUE;
    unsigned activeBlocks = (needlesCount + ctaSize - 1) / ctaSize;

    // Launch the preprocessing kernel.
    _preprocessKernel<<<activeBlocks, ctaSize>>>(view(pattern));
    SEQAN_ASSERT_EQ(cudaGetLastError(), cudaSuccess);

    // Sort the pattern.
    thrust::sort_by_key(pattern._hashes.begin(), pattern._hashes.end(), pattern._permutation.begin());
}
#endif

// ----------------------------------------------------------------------------
// Function _preprocess()                                             [Pattern]
// ----------------------------------------------------------------------------

template <typename TNeedles, typename TSpec>
inline void
_preprocess(Pattern<TNeedles, Multiple<TSpec> > & pattern)
{
    typedef Pattern<TNeedles, Multiple<TSpec> > TPattern;

    _preprocess(pattern, typename ExecSpace<TPattern>::Type());
}

// ----------------------------------------------------------------------------
// Function setHost()                                                 [Pattern]
// ----------------------------------------------------------------------------

//template <typename TNeedles, typename TSpec, typename TOtherNeedles>
//SEQAN_HOST_DEVICE inline void
//setHost(Pattern<TNeedles, Multiple<TSpec> > & pattern, TOtherNeedles const & needles)
//{
//}

// ----------------------------------------------------------------------------
// Function _initFactory()                                             [Finder]
// ----------------------------------------------------------------------------

template <typename TText, typename TPattern, typename TSpec, typename THistorySize, typename TObjectsSize>
inline void
_initFactory(Finder_<TText, TPattern, Multiple<TSpec> > & finder,
             THistorySize maxHistoryLength,
             TObjectsSize maxObjects)
{
    // Initialize the iterator factory.
    setMaxHistoryLength(finder._factory, maxHistoryLength);
    setMaxObjects(finder._factory, maxObjects);
    build(finder._factory);
}

// ----------------------------------------------------------------------------
// Function _find()                                          [Finder; ExecHost]
// ----------------------------------------------------------------------------

template <typename TText, typename TPattern, typename TSpec, typename TDelegate>
inline void
_find(Finder_<TText, TPattern, Multiple<TSpec> > & finder,
      TPattern & pattern,
      TDelegate & delegate,
      ExecHost const & /* tag */)
{
    typedef Finder_<TText, TPattern,  Multiple<TSpec> >                         TFinder;
    typedef typename View<TText>::Type                                          TTextView;
    typedef typename View<TPattern>::Type                                       TPatternView;
    typedef typename View<TFinder>::Type                                        TFinderView;
    typedef typename Needle<TPattern>::Type                                     TNeedles;
    typedef typename View<TNeedles>::Type                                       TNeedlesView;
    typedef typename Size<TNeedles>::Type                                       TSize;
    typedef FinderContext_<TTextView, TPatternView, Multiple<TSpec>, TDelegate> TFinderContext;

    // Initialize the iterator factory.
    _initFactory(finder, maxLength(needle(pattern), Parallel()) + 1, omp_get_max_threads());

    TNeedles & needles = needle(pattern);
    TSize needlesCount = length(needles);

    // Use a finder view.
    TFinderView finderView = view(finder);
    TNeedlesView needlesView = view(needles);

    // Instantiate a thread context.
    // NOTE(esiragusa): Each thread initializes its private context on firstprivate.
    TFinderContext ctx(finderView, delegate);

    // Find all needles in parallel.
    SEQAN_OMP_PRAGMA(parallel for schedule(dynamic) firstprivate(ctx))
    for (TSize needleId = 0; needleId < needlesCount; ++needleId)
    {
        clear(ctx.baseFinder);
        _setScoreThreshold(ctx.baseFinder, _getScoreThreshold(finder));
        ctx._patternIt = needleId;
        _find(ctx.baseFinder, needlesView[ctx._patternIt], ctx);
    }
}

// ----------------------------------------------------------------------------
// Function _find()                                        [Finder; ExecDevice]
// ----------------------------------------------------------------------------

#ifdef PLATFORM_CUDA
template <typename TText, typename TPattern, typename TSpec, typename TDelegate>
inline void
_find(Finder_<TText, TPattern, Multiple<TSpec> > & finder,
      TPattern /* const */ & pattern,
      TDelegate & delegate,
      ExecDevice const & /* tag */)
{
    typedef Finder_<TText, TPattern, Multiple<TSpec> >  TFinder;
    typedef typename View<TFinder>::Type                TFinderView;

    // NOTE(esiragusa): Use this to switch to persistent threads.
//    unsigned activeBlocks = cudaMaxActiveBlocks(_findKernel<TTextView, TPatternsView, TSpec, TDelegateView>, ctaSize, 0);
//    setMaxObjects(finder._factory, activeBlocks * ctaSize);

    // Compute grid size.
    unsigned ctaSize = CtaSize<TFinderView>::VALUE;
    unsigned activeBlocks = (length(needle(pattern)) + ctaSize - 1) / ctaSize;

    // Initialize the iterator factory.
    _initFactory(finder, 33 + 1, length(needle(pattern)));
    // TODO(esiragusa): Compute the longest pattern.
    // NOTE(esiragusa): back/value do not work on thrust::device_vector.
//    _initFactory(finder, length(back(needle(pattern))) + 1, length(needle(pattern)));

    // Launch the find kernel.
    _findKernel<<<activeBlocks, ctaSize>>>(view(finder), view(pattern), view(delegate));
    SEQAN_ASSERT_EQ(cudaGetLastError(), cudaSuccess);
}
#endif

// ----------------------------------------------------------------------------
// Function _find()                                                    [Finder]
// ----------------------------------------------------------------------------

template <typename TText, typename TPattern, typename TSpec, typename TDelegate>
inline void
_find(Finder_<TText, TPattern, Multiple<TSpec> > & finder, TPattern & pattern, TDelegate & delegate)
{
    typedef Finder_<TText, TPattern,  Multiple<TSpec> > TFinder;

    if (empty(needle(pattern))) return;

    _find(finder, pattern, delegate, typename ExecSpace<TFinder>::Type());
}

// ----------------------------------------------------------------------------
// Function view()                                                    [Pattern]
// ----------------------------------------------------------------------------

template <typename TNeedles, typename TSpec>
inline typename View<Pattern<TNeedles, Multiple<TSpec> > >::Type
view(Pattern<TNeedles, Multiple<TSpec> > & pattern)
{
    typename View<Pattern<TNeedles, Multiple<TSpec> > >::Type   patternView;

    patternView.data_host = view(needle(pattern));
    patternView._hashes = view(pattern._hashes);
    patternView._permutation = view(pattern._permutation);

    return patternView;
}

// ----------------------------------------------------------------------------
// Function view()                                                     [Finder]
// ----------------------------------------------------------------------------

template <typename TText, typename TPattern, typename TSpec>
inline typename View<Finder_<TText, TPattern, Multiple<TSpec> > >::Type
view(Finder_<TText, TPattern, Multiple<TSpec> > & finder)
{
    typename View<Finder_<TText, TPattern, Multiple<TSpec> > >::Type finderView;

    finderView._factory = view(finder._factory);

    return finderView;
}

// ----------------------------------------------------------------------------
// Function _textIterator()                                     [Finder Context]
// ----------------------------------------------------------------------------

template <typename TText, typename TPattern, typename TSpec, typename TDelegate>
inline SEQAN_HOST_DEVICE typename TextIterator_<TText, TPattern, Multiple<TSpec> >::Type &
_textIterator(FinderContext_<TText, TPattern, Multiple<TSpec>, TDelegate> & ctx)
{
    return _textIterator(ctx.baseFinder);
}

template <typename TText, typename TPattern, typename TSpec, typename TDelegate>
inline SEQAN_HOST_DEVICE typename TextIterator_<TText, TPattern, Multiple<TSpec> >::Type const &
_textIterator(FinderContext_<TText, TPattern, Multiple<TSpec>, TDelegate> const & ctx)
{
    return _textIterator(ctx.baseFinder);
}

// ----------------------------------------------------------------------------
// Function _getScore()
// ----------------------------------------------------------------------------

template <typename TText, typename TPattern, typename TSpec, typename TDelegate>
SEQAN_HOST_DEVICE inline typename Score_<TSpec>::Type
_getScore(FinderContext_<TText, TPattern, Multiple<TSpec>, TDelegate> const & ctx)
{
    return _getScore(ctx.baseFinder);
}

}

#endif  // #ifndef SEQAN_INDEX_FIND_INDEX_MULTIPLE_H
