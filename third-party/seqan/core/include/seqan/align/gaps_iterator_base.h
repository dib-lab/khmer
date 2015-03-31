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
// Author: Andreas Gogol-Doering <andreas.doering@mdc-berlin.de>
// Author: Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>
// ==========================================================================

// TODO(holtgrew): Switch to Host interface.

#ifndef SEQAN_CORE_INCLUDE_SEQAN_ALIGN_GAPS_ITERATOR_BASE_H_
#define SEQAN_CORE_INCLUDE_SEQAN_ALIGN_GAPS_ITERATOR_BASE_H_

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// Internally used tag for creating iterators at the begin of containers.
struct Begin__;
typedef Tag<Begin__> Begin_;

// Internally used tag for creating iterators at the end of containers.
struct End__;
typedef Tag<End__> End_;

// Internally used tag for creating iterators inside of containers.
struct Position__;
typedef Tag<Position__> Position_;

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

template <typename TSpec>
struct GapsIterator;

// ============================================================================
// Metafunctions
// ============================================================================

// ----------------------------------------------------------------------------
// Metafunction Position
// ----------------------------------------------------------------------------

template <typename TGaps, typename TSpec>
struct Position<Iter<TGaps, GapsIterator<TSpec> > > :
            Position<TGaps>
{};

template <typename TGaps, typename TSpec>
struct Position<Iter<TGaps, GapsIterator<TSpec> > const> :
            Position<Iter<TGaps, GapsIterator<TSpec> > >
{};

// ----------------------------------------------------------------------------
// Metafunction Difference
// ----------------------------------------------------------------------------

template <typename TGaps, typename TSpec>
struct Difference<Iter<TGaps, GapsIterator<TSpec> > > :
            Difference<TGaps>
{};

template <typename TGaps, typename TSpec>
struct Difference<Iter<TGaps, GapsIterator<TSpec> > const> :
            Difference<Iter<TGaps, GapsIterator<TSpec> > >
{};

// ----------------------------------------------------------------------------
// Metafunction Source
// ----------------------------------------------------------------------------

// TODO(holtgrew): Should this be Host? or SourceIterator?

template <typename TGaps, typename TSpec>
struct Source<Iter<TGaps, GapsIterator<TSpec> > >
{
	typedef typename Source<TGaps>::Type TSource_;
	typedef typename Iterator<TSource_, Rooted>::Type Type;
};

template <typename TGaps, typename TSpec>
struct Source<Iter<TGaps, GapsIterator<TSpec> > const>
{
	typedef typename Source<TGaps>::Type TSource_;
	typedef typename Iterator<TSource_, Rooted>::Type Type;
};

// ----------------------------------------------------------------------------
// Metafunction Value
// ----------------------------------------------------------------------------

template <typename TGaps, typename TSpec>
struct Value<Iter<TGaps, GapsIterator<TSpec> > >
{
	typedef typename Source<Iter<TGaps, GapsIterator<TSpec> > >::Type TSource_;
	typedef typename Value<TSource_>::Type TSourceValue_;
    //typedef TSourceValue_ Type;
    // TODO(holtgrew): We really want gapped values here but there are issues...
	typedef typename GappedValueType<TSourceValue_>::Type Type;
};

template <typename TGaps, typename TSpec>
struct Value<Iter<TGaps, GapsIterator<TSpec> > const> :
            Value<Iter<TGaps, GapsIterator<TSpec> > > {};

// ----------------------------------------------------------------------------
// Metafunction GetValue
// ----------------------------------------------------------------------------

template <typename TGaps, typename TSpec>
struct GetValue<Iter<TGaps, GapsIterator<TSpec> > > :
	Value<Iter<TGaps, GapsIterator<TSpec> > >
{
};

template <typename TGaps, typename TSpec>
struct GetValue<Iter<TGaps, GapsIterator<TSpec> > const> :
	Value<Iter<TGaps, GapsIterator<TSpec> > const>
{
};

// ----------------------------------------------------------------------------
// Metafunction Reference
// ----------------------------------------------------------------------------

template <typename TGaps, typename TSpec>
struct Reference<Iter<TGaps, GapsIterator<TSpec> > >
{
	typedef Iter<TGaps, GapsIterator<TSpec> > TIterator_;
	typedef Proxy<IteratorProxy<TIterator_> > Type;
};

template <typename TGaps, typename TSpec>
struct Reference<Iter<TGaps, GapsIterator<TSpec> > const>
{
	typedef Iter<TGaps, GapsIterator<TSpec> const > TIterator_;
	typedef Proxy<IteratorProxy<TIterator_> > Type;
};

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function operator++
// ----------------------------------------------------------------------------

// TODO(holtgrew): Could be general forward

template <typename TGaps, typename TSpec>
inline Iter<TGaps, GapsIterator<TSpec> > & 
operator++(Iter<TGaps, GapsIterator<TSpec> > & it)
{
	goNext(it);
	return it;
}

template <typename TGaps, typename TSpec>
inline Iter<TGaps, GapsIterator<TSpec> >
operator++(Iter<TGaps, GapsIterator<TSpec> > & it, int)
{
	Iter<TGaps, GapsIterator<TSpec> > ret = it;
	goNext(it);
	return ret;
}

// ----------------------------------------------------------------------------
// Function operator--
// ----------------------------------------------------------------------------

template <typename TGaps, typename TSpec>
inline Iter<TGaps, GapsIterator<TSpec> > & 
operator--(Iter<TGaps, GapsIterator<TSpec> > & it)
{
	goPrevious(it);
	return it;
}

template <typename TGaps, typename TSpec>
inline Iter<TGaps, GapsIterator<TSpec> >
operator--(Iter<TGaps, GapsIterator<TSpec> > & it, int)
{
	Iter<TGaps, GapsIterator<TSpec> > ret = it;
	goPrevious(it);
	return ret;
}

// ----------------------------------------------------------------------------
// Function insertGap()
// ----------------------------------------------------------------------------

// Forward to insertGaps() which has to be implemented by the specific gap
// iterator.

template <typename TGaps, typename TSpec>
inline void
insertGap(Iter<TGaps, GapsIterator<TSpec> > & it)
{
	insertGaps(it, 1);
}

// ----------------------------------------------------------------------------
// Function removeGap()
// ----------------------------------------------------------------------------

// Forward to removeGaps() which has to be implemented by the specific gap
// iterator.

template <typename TGaps, typename TSpec>
inline typename Size<TGaps>::Type
removeGap(Iter<TGaps, GapsIterator<TSpec> > & it)
{
	return removeGaps(it, 1);
}

// ----------------------------------------------------------------------------
// Function assignValue()
// ----------------------------------------------------------------------------

// TODO(holtgrew): Const consistency problems.

template <typename TGaps, typename TSpec, typename TValue>
inline void
assignValue(Iter<TGaps, GapsIterator<TSpec> > & me,
			TValue const & val)
{
	if (!isGap(me)) 
	{
		assignValue(source(me), val);
	}
    // TODO(holtgrew): Else, inserting gaps is problematic...
}

template <typename TGaps, typename TSpec, typename TValue>
inline void
assignValue(Iter<TGaps, GapsIterator<TSpec> > const & me,
			TValue const & val)
{
	if (!isGap(me)) 
	{
		assignValue(source(me), val);
	}
}

// ----------------------------------------------------------------------------
// Function container()
// ----------------------------------------------------------------------------

template <typename TGaps, typename TSpec>
inline TGaps &
container(Iter<TGaps, GapsIterator<TSpec> > & me)
{
	return *me._container;
}

template <typename TGaps, typename TSpec>
inline TGaps &
container(Iter<TGaps, GapsIterator<TSpec> > const & me)
{
	return *me._container;
}

// ----------------------------------------------------------------------------
// Function source
// ----------------------------------------------------------------------------

// Returns host iterator.

// TODO(holtgrew): Non-const version is superflous.
template <typename TGaps, typename TSpec>
inline typename Source<Iter<TGaps, GapsIterator<TSpec> > >::Type /*returns copy*/
source(Iter<TGaps, GapsIterator<TSpec> > & it)
{
    return iter(container(it), toSourcePosition(container(it), position(it)));
}

template <typename TGaps, typename TSpec>
inline typename Source<Iter<TGaps, GapsIterator<TSpec> > const>::Type /*returns copy*/
source(Iter<TGaps, GapsIterator<TSpec> > const & it)
{
    return iter(container(source(it)), toSourcePosition(container(it), position(it)));
}

// TODO(holtgrew): setSource? setContainer?

// ----------------------------------------------------------------------------
// Function operator+=
// ----------------------------------------------------------------------------

template <typename TGaps, typename TSpec, typename TDiff>
inline Iter<TGaps, GapsIterator<TSpec> > &
operator+=(Iter<TGaps, GapsIterator<TSpec> > & it, TDiff diff)
{
    goFurther(it, diff);
    return it;
}

// ----------------------------------------------------------------------------
// Function operator-=
// ----------------------------------------------------------------------------

template <typename TGaps, typename TSpec, typename TDiff>
inline Iter<TGaps, GapsIterator<TSpec> > &
operator-=(Iter<TGaps, GapsIterator<TSpec> > & it, TDiff diff)
{
    goFurther(it, -(__int64)(diff));
    return it;
}

// ----------------------------------------------------------------------------
// Function goFurther
// ----------------------------------------------------------------------------

// TODO(holtgrew): Implementation could be faster.
template <typename TGaps, typename TSpec, typename TDifference>
inline void
goFurther(Iter<TGaps, GapsIterator<TSpec> > & it,
		  TDifference steps)
{
	typedef typename MakeSigned<TDifference>::Type TSignedDifference;
    if (steps > TDifference(0))
        for (; steps; --steps)
            goNext(it);
    else
        for (; -static_cast<TSignedDifference>(steps); ++steps)
            goPrevious(it);
}

// ----------------------------------------------------------------------------
// Function isClipped()
// ----------------------------------------------------------------------------

template <typename TGaps, typename TSpec>
inline bool
isClipped(Iter<TGaps, GapsIterator<TSpec> > const &)
{
    return false;
}

}  // namespace seqan

#endif  // SEQAN_CORE_INCLUDE_SEQAN_ALIGN_GAPS_ITERATOR_BASE_H_
