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
// Author: Enrico Siragusa <enrico.siragusa@fu-berlin.de>
// ==========================================================================

#ifndef SEQAN_BASIC_FUNCTORS_H_
#define SEQAN_BASIC_FUNCTORS_H_

namespace seqan {

// ============================================================================
// Functors
// ============================================================================

// ----------------------------------------------------------------------------
// Functor OrFunctor
// ----------------------------------------------------------------------------

template <typename TFunctor1, typename TFunctor2>
struct OrFunctor
{
    TFunctor1 func1;
    TFunctor2 func2;

    OrFunctor()
    {}

    OrFunctor(TFunctor1 const &func1, TFunctor2 const &func2):
        func1(func1), func2(func2)
    {}

    template <typename TValue>
    bool operator() (TValue const & val)
    {
        return func1(val) || func2(val);
    }
};

// ----------------------------------------------------------------------------
// Functor AndFunctor
// ----------------------------------------------------------------------------

template <typename TFunctor1, typename TFunctor2>
struct AndFunctor
{
    TFunctor1 func1;
    TFunctor2 func2;

    AndFunctor()
    {}

    AndFunctor(TFunctor1 const &func1, TFunctor2 const &func2):
        func1(func1), func2(func2)
    {}

    template <typename TValue>
    bool operator() (TValue const & val)
    {
        return func1(val) && func2(val);
    }
};

// ----------------------------------------------------------------------------
// Functor NotFunctor
// ----------------------------------------------------------------------------

template <typename TFunctor>
struct NotFunctor
{
    TFunctor func;

    NotFunctor()
    {}

    NotFunctor(TFunctor const &func):
        func(func)
    {}

    template <typename TValue>
    bool operator() (TValue const & val)
    {
        return !func(val);
    }
};

// ----------------------------------------------------------------------------
// Functor CountDownFunctor
// ----------------------------------------------------------------------------

template <typename TFunctor = True, __uint64 REMAINING = 0>
struct CountDownFunctor
{
    __uint64 remaining;
    TFunctor func;

    CountDownFunctor(__uint64 remaining = REMAINING):
        remaining(remaining)
    {}

    CountDownFunctor(__uint64 remaining, TFunctor const &func):
        remaining(remaining),
        func(func)
    {}

    template <typename TValue>
    bool operator() (TValue const & val)
    {
        if (remaining == 0)
            return true;
        if (func(val))
            --remaining;
        return false;
    }

    operator bool()
    {
        return remaining == 0;
    }
};

// ----------------------------------------------------------------------------
// Functor CountFunctor
// ----------------------------------------------------------------------------

template <typename TFunctor = True>
struct CountFunctor
{
    __uint64 count;
    TFunctor func;

    CountFunctor() : count(0)
    {}

    CountFunctor(TFunctor const & func) : count(0), func(func)
    {}

    template <typename TValue>
    bool operator() (TValue const & val)
    {
        if (func(val))
            ++count;
        return false;
    }

    operator __uint64() const
    {
        return count;
    }
};

template <typename TFunctor>
inline void clear(CountFunctor<TFunctor> &func)
{
    func.count = 0;
}

template <typename TFunctor>
inline __uint64 & value(CountFunctor<TFunctor> &func)
{
    return func.count;
}

}   // namespace seqan

#endif // SEQAN_BASIC_FUNCTORS_H_
