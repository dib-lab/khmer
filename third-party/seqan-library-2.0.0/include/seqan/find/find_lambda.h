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
// Author: Enrico Siragusa <enrico.siragusa@fu-berlin.de>
// ==========================================================================

#ifndef SEQAN_FIND_LAMBDA_H
#define SEQAN_FIND_LAMBDA_H

namespace seqan {

// ============================================================================
// Metafunctions
// ============================================================================

// ----------------------------------------------------------------------------
// Metafunction DefaultFind<>
// ----------------------------------------------------------------------------

template <typename TText, typename TPattern>
struct DefaultFind
{
    typedef void Type;
};

// ----------------------------------------------------------------------------
// Metafunction FindState_<>
// ----------------------------------------------------------------------------

template <typename TText, typename TPattern, typename TAlgorithm>
struct FindState_
{
    typedef Nothing Type;
};

// ----------------------------------------------------------------------------
// Metafunction HasStatesPool_<>
// ----------------------------------------------------------------------------

template <typename TState, typename TThreading>
struct HasStatesPool_
{
    typedef typename IsSameType<TThreading, Parallel>::Type         IsParallel;
    typedef typename Not<IsSameType<TState, Nothing> >::Type        IsStateful;
    typedef typename And<IsParallel, IsStateful>::Type              Type;
};

// ----------------------------------------------------------------------------
// Metafunction StatesPool_<>
// ----------------------------------------------------------------------------

template <typename TState, typename TThreading>
struct StatesPool_
{
    typedef typename HasStatesPool_<TState, TThreading>::Type           HasStatesPool;
    typedef typename If<HasStatesPool, String<TState>, TState>::Type    Type;
};

// ============================================================================
// Classes
// ============================================================================

// ----------------------------------------------------------------------------
// Class FindDelegator_
// ----------------------------------------------------------------------------

template <typename TDelegate, typename TNeedlesIt>
struct FindDelegator_
{
    TDelegate & delegate;
    TNeedlesIt & needlesIt;

    FindDelegator_(TDelegate & delegate, TNeedlesIt & needlesIt) :
        delegate(delegate),
        needlesIt(needlesIt)
    {}

    template <typename TTextIt, typename TScore>
    void operator()(TTextIt & textIt, TScore score)
    {
        delegate(textIt, needlesIt, score);
    }
};

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function find(text, pattern, errors, [](...){}, Algorithm());
// ----------------------------------------------------------------------------

template <typename TText, typename TPattern, typename TThreshold, typename TDelegate, typename TAlgorithm>
inline void find(TText & text, TPattern const & pattern, TThreshold threshold, TDelegate && delegate, TAlgorithm)
{
    typedef typename FindState_<TText, TPattern const, TAlgorithm>::Type    TState;

    TState state;
    _findStateInit(state, text, pattern, threshold, TAlgorithm());
    _findImpl(state, text, pattern, threshold, delegate, TAlgorithm());
}

// ----------------------------------------------------------------------------
// Function find(text, pattern, [](...){});
// ----------------------------------------------------------------------------

template <typename TText, typename TPattern, typename TDelegate>
inline void find(TText & text, TPattern const & pattern, TDelegate && delegate)
{
    typedef typename DefaultFind<TText, TPattern const>::Type   TAlgorithm;

    find(text, pattern, unsigned(), delegate, TAlgorithm());
}

// ----------------------------------------------------------------------------
// Function _findStateInit()
// ----------------------------------------------------------------------------
// Stateless algorithm.

template <typename TState, typename TText, typename TPattern, typename TThreshold, typename TAlgorithm>
inline void
_findStateInit(TState & /* s */, TText const & /* t */, TPattern const & /* p */, TThreshold /* t */, TAlgorithm) {}

// ----------------------------------------------------------------------------
// Function _findStatesPoolInit()
// ----------------------------------------------------------------------------

template <typename TState>
inline void _findStatesPoolInit(TState & /* pool */, False) {}

template <typename TStates>
inline void _findStatesPoolInit(TStates & pool, True)
{
    resize(pool, omp_get_num_threads(), Exact());
}

// ----------------------------------------------------------------------------
// Function _findPickState()
// ----------------------------------------------------------------------------

template <typename TState>
inline TState &
_findPickState(TState & state, False)
{
    return state;
}

template <typename TStates>
inline typename Value<TStates>::Type &
_findPickState(TStates & pool, True)
{
    return pool[omp_get_thread_num()];
}

// ----------------------------------------------------------------------------
// Function find(text, needles, errors, [](...){}, Algorithm(), Parallel());
// ----------------------------------------------------------------------------

template <typename TText, typename TNeedle_, typename TSSetSpec,
          typename TThreshold, typename TDelegate, typename TAlgorithm, typename TThreading>
inline void find(TText & text,
                 StringSet<TNeedle_, TSSetSpec> const & needles,
                 TThreshold threshold,
                 TDelegate && delegate,
                 TAlgorithm,
                 TThreading)
{
    typedef StringSet<TNeedle_, TSSetSpec> const                    TNeedles;
    typedef typename Value<TNeedles>::Type                          TNeedle;
    typedef typename Iterator<TNeedles, Rooted>::Type               TNeedlesIt;
    typedef typename Reference<TNeedlesIt>::Type                    TNeedleRef;

    typedef typename FindState_<TText, TNeedle, TAlgorithm>::Type   TState;
    typedef typename StatesPool_<TState, TThreading>::Type          TStatesPool;
    typedef typename HasStatesPool_<TState, TThreading>::Type       HasStatesPool;

    typedef FindDelegator_<TDelegate, TNeedlesIt const>             TDelegator;

    TStatesPool pool;
    _findStatesPoolInit(pool, HasStatesPool());

    iterate(needles, [&](TNeedlesIt const & needlesIt)
    {
        TState & state = _findPickState(pool, HasStatesPool());
        TNeedleRef needle = value(needlesIt);
        TDelegator delegator(delegate, needlesIt);

        _findStateInit(state, text, needle, threshold, TAlgorithm());
        _findImpl(state, text, needle, threshold, delegator, TAlgorithm());
    },
    Rooted(), TThreading());
}

// ----------------------------------------------------------------------------
// Function find(text, needles, errors, [](...){}, Algorithm());
// ----------------------------------------------------------------------------

template <typename TText, typename TNeedle, typename TSSetSpec,
          typename TThreshold, typename TDelegate, typename TAlgorithm>
inline void find(TText & text,
                 StringSet<TNeedle, TSSetSpec> const & needles,
                 TThreshold threshold,
                 TDelegate && delegate,
                 TAlgorithm)
{
    find(text, needles, threshold, delegate, TAlgorithm(), Serial());
}

}

#endif  // #ifndef SEQAN_FIND_LAMBDA_H
