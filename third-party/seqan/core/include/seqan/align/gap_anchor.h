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

#ifndef SEQAN_CORE_INCLUDE_SEQAN_ALIGN_GAP_ANCHOR_H_
#define SEQAN_CORE_INCLUDE_SEQAN_ALIGN_GAP_ANCHOR_H_

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// TODO(holtgrew): Document?
// Sorting tags (just for lower_bound and upper_bound, positions are always sorted)

struct SortSeqPos_;
typedef Tag<SortSeqPos_> const SortSeqPos;

struct SortGapPos_;
typedef Tag<SortGapPos_> const SortGapPos;

// ----------------------------------------------------------------------------
// Specialization GapAnchor
// ----------------------------------------------------------------------------

// TODO(holtgrew): Make a class instead of a struct?

/*!
 * @class GapAnchor
 * @headerfile <seqan/align.h>
 * @brief Stores the position of an alignment character in sequence-space and in gap-space.
 * 
 * @signature template <typename TPosition>
 *            struct GapAnchor;
 * 
 * @tparam TPos Type to store gapped/ungapped positions.
 * 
 * @section Remarks
 * 
 * Value types of the <tt>gaps</tt> strings in @link ReadStoreElement @endlink and @link ContigStoreElement @endlink.
 */

/**
.Class.GapAnchor
..summary:Stores the position of an alignment character in sequence-space and in gap-space.
..cat:Alignments
..signature:GapAnchor<TPos>
..param.TPos:Type to store gapped/ungapped positions.
..remarks:Value types of the $gaps$ strings in @Class.ReadStoreElement@ and @Class.ContigStoreElement@.

.Memfunc.GapAnchor#GapAnchor
..summary:Constructor
..signature:GapAnchor<TPos> ()
..signature:GapAnchor<TPos> (TPos seqPos, TPos gapPos)
..param.seqPos:Sequence character position in the ungapped sequence.
..param.gapPos:Sequence character position in the gapped sequence.
..remarks:Default constructor sets both positions to $0$.
..class:Class.GapAnchor
.Memvar.GapAnchor#seqPos
..summary:Sequence character position in the ungapped sequence.
..class:Class.GapAnchor
.Memvar.GapAnchor#gapPos
..summary:Sequence character position in the gapped sequence.
..class:Class.GapAnchor
..include:seqan/store.h
*/

// We store gap anchors only for the first text character behind a gap or a clipped sequence character

template <typename TPos>
struct GapAnchor
{
/*!
 * @var VariableType GapAnchor::seqPos
 * @brief Sequence character position in the ungapped sequence.
 */
	TPos	seqPos;			// sequence character position in the ungapped sequence

/*!
 * @var VariableType GapAnchor::gapPos
 * @brief Sequence character position in the gapped sequence.
 */
	TPos	gapPos;			// sequence character position in the gapped sequence

/*!
 * @fn GapAnchor::GapAnchor
 * 
 * @brief Constructor
 * 
 * @signature GapAnchor::GapAnchor([other])
 * @signature GapAnchor::GapAnchor(seqPos, gapPos)
 *
 * @param other  GapAnchor object to copy from.
 * @param seqPos Sequence character position in the ungapped sequence.
 * @param gapPos Sequence character position in the gapped sequence.
 * 
 * @section Remarks
 * 
 * Default constructor sets both positions to <tt>0</tt>.
 */

	GapAnchor() : seqPos(0), gapPos(0) {}
	GapAnchor(TPos sP, TPos gP) : seqPos(sP), gapPos(gP) {}

	template <typename TPos_>
	GapAnchor(GapAnchor<TPos_> const &other)
	{
		seqPos = other.seqPos;
		gapPos = other.gapPos;
	}

	template <typename TPos_>
	inline GapAnchor const &
	operator = (GapAnchor<TPos_> const &other)
	{
		seqPos = other.seqPos;
		gapPos = other.gapPos;
		return *this;
	} 

	template <typename TOther>
	inline bool
	operator == (TOther const &other) const
	{
		return seqPos == other.seqPos && gapPos == other.gapPos;
	} 

	template <typename TOther>
	inline bool
	operator != (TOther const &other) const
	{
		return !(*this == other);
	} 

	template <typename TOther>
	inline bool
	operator < (TOther const &other) const
	{
		return seqPos < other.seqPos || gapPos < other.gapPos;
	} 

	template <typename TOther>
	inline bool
	operator > (TOther const &other) const
	{
		return seqPos > other.seqPos || gapPos > other.gapPos;
	} 

	template <typename TOther>
	inline bool
	operator <= (TOther const &other) const
	{
		return seqPos < other.seqPos || gapPos <= other.gapPos;
	} 

	template <typename TOther>
	inline bool
	operator >= (TOther const &other) const
	{
		return seqPos > other.seqPos || gapPos >= other.gapPos;
	}
};

// ============================================================================
// Metafunctions
// ============================================================================

// ----------------------------------------------------------------------------
// Metafunction Size                                                [GapAnchor]
// ----------------------------------------------------------------------------

template <typename TPos>
struct Size<GapAnchor<TPos> >
{
    typedef TPos Type;
};

template <typename TPos>
struct Size<GapAnchor<TPos> const> : public Size<GapAnchor<TPos> >
{};

// ----------------------------------------------------------------------------
// Metafunction Position                                            [GapAnchor]
// ----------------------------------------------------------------------------

template <typename TPos>
struct Position<GapAnchor<TPos> >
{
    typedef TPos Type;
};

template <typename TPos>
struct Position<GapAnchor<TPos> const> : public Position<GapAnchor<TPos> >
{};

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Functor _LessGapAnchor
// ----------------------------------------------------------------------------

// TODO(holtgrew): Leading underscore should become a trailing one.

template <typename TGapAnchor, typename TTag>
struct _LessGapAnchor;

template <typename TGapAnchor>
struct _LessGapAnchor<TGapAnchor, SortSeqPos> :
	public ::std::binary_function<TGapAnchor, TGapAnchor, bool>
{
	inline bool 
	operator() (TGapAnchor const& a1, TGapAnchor const& a2) const {
		return (a1.seqPos) < (a2.seqPos);
	}
};

template <typename TGapAnchor>
struct _LessGapAnchor<TGapAnchor, SortGapPos> :
	public ::std::binary_function<TGapAnchor, TGapAnchor, bool>
{
	inline bool 
	operator() (TGapAnchor const& a1, TGapAnchor const& a2) const {
		return (a1.gapPos) < (a2.gapPos);
	}
};

// ----------------------------------------------------------------------------
// Function lowerBoundGapAnchor()
// ----------------------------------------------------------------------------

template <typename TGapAnchor, typename TSearchValue>
inline typename Iterator<TGapAnchor const, Standard>::Type
lowerBoundGapAnchor(TGapAnchor const & gaps, 
					TSearchValue const val,
					SortSeqPos) 
{
	typedef typename Value<TGapAnchor>::Type TGapAnchorElement;
	TGapAnchorElement el;
	el.seqPos = val;
	return ::std::lower_bound(
		begin(gaps, Standard()), 
		end(gaps, Standard()), 
		el,
		_LessGapAnchor<typename Value<TGapAnchor>::Type, SortSeqPos const>() );
}

template <typename TGapAnchor, typename TSearchValue>
inline typename Iterator<TGapAnchor, Standard>::Type
lowerBoundGapAnchor(TGapAnchor & gaps, 
					TSearchValue const val,
					SortSeqPos) 
{
	typedef typename Value<TGapAnchor>::Type TGapAnchorElement;
	TGapAnchorElement el;
	el.seqPos = val;
	return ::std::lower_bound(
		begin(gaps, Standard()), 
		end(gaps, Standard()), 
		el,
		_LessGapAnchor<typename Value<TGapAnchor>::Type, SortSeqPos const>() );
}

template <typename TGapAnchor, typename TSearchValue>
inline typename Iterator<TGapAnchor const, Standard>::Type
lowerBoundGapAnchor(TGapAnchor const & gaps, 
					TSearchValue const val,
					SortGapPos) 
{
	typedef typename Value<TGapAnchor>::Type TGapAnchorElement;
	TGapAnchorElement el;
	el.gapPos = val;
	return ::std::lower_bound(
		begin(gaps, Standard()), 
		end(gaps, Standard()), 
		el,
		_LessGapAnchor<typename Value<TGapAnchor>::Type, SortGapPos const>() );
}

template <typename TGapAnchor, typename TSearchValue>
inline typename Iterator<TGapAnchor, Standard>::Type
lowerBoundGapAnchor(TGapAnchor & gaps, 
					TSearchValue const val,
					SortGapPos) 
{
	typedef typename Value<TGapAnchor>::Type TGapAnchorElement;
	TGapAnchorElement el;
	el.gapPos = val;
	return ::std::lower_bound(
		begin(gaps, Standard()), 
		end(gaps, Standard()), 
		el,
		_LessGapAnchor<typename Value<TGapAnchor>::Type, SortGapPos const>() );
}

// ----------------------------------------------------------------------------
// Function upperBoundGapAnchor()
// ----------------------------------------------------------------------------

template <typename TGapAnchors, typename TSearchValue>
inline typename Iterator<TGapAnchors const, Standard>::Type
upperBoundGapAnchor(TGapAnchors const & gaps,
					TSearchValue const val,
					SortSeqPos) 
{
	typedef typename Value<TGapAnchors>::Type TGapAnchorElement;
	TGapAnchorElement el;
	el.seqPos = val;
	return ::std::upper_bound(
		begin(gaps, Standard()), 
		end(gaps, Standard()), 
		el,
		_LessGapAnchor<typename Value<TGapAnchors>::Type, SortSeqPos const>() );
}

template <typename TGapAnchors, typename TSearchValue>
inline typename Iterator<TGapAnchors, Standard>::Type
upperBoundGapAnchor(TGapAnchors & gaps,
					TSearchValue const val,
					SortSeqPos) 
{
	typedef typename Value<TGapAnchors>::Type TGapAnchorElement;
	TGapAnchorElement el;
	el.seqPos = val;
	return ::std::upper_bound(
		begin(gaps, Standard()), 
		end(gaps, Standard()), 
		el,
		_LessGapAnchor<typename Value<TGapAnchors>::Type, SortSeqPos const>() );
}

template <typename TGapAnchors, typename TSearchValue>
inline typename Iterator<TGapAnchors const, Standard>::Type
upperBoundGapAnchor(TGapAnchors const & gaps, 
					TSearchValue const val,
					SortGapPos) 
{
	typedef typename Value<TGapAnchors>::Type TGapAnchorElement;
	TGapAnchorElement el;
	el.gapPos = val;
	return ::std::upper_bound(
		begin(gaps, Standard()), 
		end(gaps, Standard()), 
		el,
		_LessGapAnchor<typename Value<TGapAnchors>::Type, SortGapPos const>() );
}

template <typename TGapAnchors, typename TSearchValue>
inline typename Iterator<TGapAnchors, Standard>::Type
upperBoundGapAnchor(TGapAnchors & gaps, 
					TSearchValue const val,
					SortGapPos) 
{
	typedef typename Value<TGapAnchors>::Type TGapAnchorElement;
	TGapAnchorElement el;
	el.gapPos = val;
	return ::std::upper_bound(
		begin(gaps, Standard()), 
		end(gaps, Standard()), 
		el,
		_LessGapAnchor<typename Value<TGapAnchors>::Type, SortGapPos const>() );
}

}  // namespace seqan

#endif  // #ifndef SEQAN_CORE_INCLUDE_SEQAN_ALIGN_GAP_ANCHOR_H_
