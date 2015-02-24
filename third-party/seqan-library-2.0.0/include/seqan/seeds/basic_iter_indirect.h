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
// A wrapper for an iterator that dereferences the wrapped iterator twice
// on dereferentiation.
// ==========================================================================

#ifndef SEQAN_SEEDS_BASIC_ITER_INDIRECT_H_
#define SEQAN_SEEDS_BASIC_ITER_INDIRECT_H_

namespace seqan {

// ===========================================================================
// Enums, Tags, Classes, Specializations
// ===========================================================================

template <typename TWrappedIter>
struct Indirect {};

// TODO(holtgrew): Make this an Iter specialization?
//
// Iterator that dereferences twice on derefentiation.
template <typename TContainer, typename TWrappedIter>
class Iter<TContainer, Indirect<TWrappedIter> >
{
public:
    TWrappedIter _wrappedIter;

    Iter() : _wrappedIter(0)
    {}

    Iter(Iter const & other) : _wrappedIter(other._wrappedIter)
    {}

    Iter(TWrappedIter const & wrappedIter) : _wrappedIter(wrappedIter)
    {}
};

// ===========================================================================
// Metafunctions
// ===========================================================================

template <typename TContainer, typename TWrappedIter>
struct Reference<Iter<TContainer, Indirect<TWrappedIter> > >
{
    typedef typename Reference<TContainer>::Type Type;
};

// ===========================================================================
// Functions
// ===========================================================================

template <typename TContainer, typename TWrappedIter>
inline bool
operator==(Iter<TContainer, Indirect<TWrappedIter> > const & a, Iter<TContainer, Indirect<TWrappedIter> > const & b)
{
    SEQAN_CHECKPOINT;
    return a._wrappedIter == b._wrappedIter;
}


template <typename TContainer, typename TWrappedIter>
inline bool
operator!=(Iter<TContainer, Indirect<TWrappedIter> > const & a, Iter<TContainer, Indirect<TWrappedIter> > const & b)
{
    SEQAN_CHECKPOINT;
    return a._wrappedIter != b._wrappedIter;
}

template <typename TContainer, typename TWrappedIter, typename TDiff>
Iter<TContainer, Indirect<TWrappedIter> >
operator+(Iter<TContainer, Indirect<TWrappedIter> > const & it, TDiff const diff)
{
    SEQAN_CHECKPOINT;
    typedef Iter<TContainer, Indirect<TWrappedIter> > TIter;
    return TIter(it._wrappedIter + diff);
}


template <typename TContainer, typename TWrappedIter>
inline Iter<TContainer, Indirect<TWrappedIter> > &
operator++(Iter<TContainer, Indirect<TWrappedIter> > & iter)
{
    SEQAN_CHECKPOINT;
    ++iter._wrappedIter;
    return iter;
}


template <typename TContainer, typename TWrappedIter>
inline Iter<TContainer, Indirect<TWrappedIter> >
operator++(Iter<TContainer, Indirect<TWrappedIter> > & iter, int /*postfix*/)
{
    SEQAN_CHECKPOINT;
    typedef Iter<TContainer, Indirect<TWrappedIter> > TIter;
    TIter tmp(iter);
    ++iter;
    return tmp;
}


template <typename TContainer, typename TWrappedIter>
inline Iter<TContainer, Indirect<TWrappedIter> > &
operator--(Iter<TContainer, Indirect<TWrappedIter> > & iter)
{
    SEQAN_CHECKPOINT;
    --iter._wrappedIter;
    return iter;
}


template <typename TContainer, typename TWrappedIter>
inline Iter<TContainer, Indirect<TWrappedIter> >
operator--(Iter<TContainer, Indirect<TWrappedIter> > & iter, int /*postfix*/)
{
    SEQAN_CHECKPOINT;
    typedef Iter<TContainer, Indirect<TWrappedIter> > TIter;
    TIter tmp(iter);
    --iter;
    return tmp;
}


template <typename TContainer, typename TWrappedIter>
inline typename Value<TWrappedIter>::Type &
operator*(Iter<TContainer, Indirect<TWrappedIter> > & iter)
{
    SEQAN_CHECKPOINT;
    return **iter._wrappedIter;
}


template <typename TContainer, typename TWrappedIter>
inline typename Value<TWrappedIter>::Type &
operator*(Iter<TContainer, Indirect<TWrappedIter> > const & iter)
{
    SEQAN_CHECKPOINT;
    return **iter._wrappedIter;
}


template <typename TContainer, typename TWrappedIter>
inline typename Reference<Iter<TContainer, Indirect<TWrappedIter> > >::Type
value(Iter<TContainer, Indirect<TWrappedIter> > & iter)
{
    SEQAN_CHECKPOINT;
    return **iter._wrappedIter;
}


template <typename TContainer, typename TWrappedIter>
inline typename Reference<Iter<TContainer, Indirect<TWrappedIter> > >::Type
value(Iter<TContainer, Indirect<TWrappedIter> > const & iter)
{
    SEQAN_CHECKPOINT;
    return **iter._wrappedIter;
}

}  // namespace seqan

#endif  // #ifndef SEQAN_SEEDS_BASIC_ITER_INDIRECT_H_
