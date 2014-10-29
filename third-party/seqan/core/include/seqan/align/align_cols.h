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
// ==========================================================================
// AlignCols is a wrapper around Align that allows the virtual access to the
// rows of an alignment.
// ==========================================================================

#ifndef SEQAN_CORE_INCLUDE_SEQAN_ALIGN_ALIGN_COLS_H_
#define SEQAN_CORE_INCLUDE_SEQAN_ALIGN_ALIGN_COLS_H_

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// ----------------------------------------------------------------------------
// Class AlignCols
// ----------------------------------------------------------------------------

/*!
 * @class AlignCols
 * @implements EqualityComparableConcept
 * @implements RandomAccessContainerConcept
 * @headerfile <seqan/align.h>
 * @brief Pseudo columns container for row-based Align class.
 *
 * @signature template <typename TAlign>
 *            class AlignCols<TAlign>;
 *
 * @tparam TAlign The Align type.
 */

/**
.Class.AlignCols:
..cat:Alignments
..summary:Pseudo columns container for row-based alignment classes.
..signature:AlignCols<TAlign>
..param.TAlign:Alignment type.
...metafunction:Metafunction.Host
..remarks:
This class emulates a container of columns on alignment classes that store the alignment in a container of rows.
Note that accessing a row-based alignment column-wise can be significantly slower than accessing the alignment row-wise.
..see:Class.Align
..include:seqan/align.h
 */

template <typename TAlign>
struct AlignCols
{
    // TODO(holtgrew): Do we need this mutable?
    mutable TAlign * data_align;

    AlignCols() :
        data_align(0)
    {}


    AlignCols(TAlign & align) : data_align(&align)
    {}

    template <typename TPosition>
    inline typename Value<AlignCols>::Type
    operator[](TPosition _pos)
    {
        return value(*this, _pos);
    }

    template <typename TPosition>
    inline typename Value<AlignCols const>::Type
    operator[](TPosition _pos) const
    {
        return value(*this, _pos);
    }
};

// ----------------------------------------------------------------------------
// Specialization AlignCols
// ----------------------------------------------------------------------------

/**
.Spec.AlignColIterator:
..cat:Iterators
..summary:Iterator for @Class.AlignCols@ pseudo container.
..signature:Iter< TAlign, AlignColIterator<TSpec> >
..param.TSpec:Specialization tag.
..general:Class.Iter
..see:Class.AlignCols
..include:seqan/align.h
*/

template <typename TSpec>
struct AlignColIterator;

// ============================================================================
// Metafunctions
// ============================================================================

// ----------------------------------------------------------------------------
// Metafunction Host
// ----------------------------------------------------------------------------

// TODO(holtgrew): Add HostedTypeConcept and make AlignCols object implement the concept.

///.Metafunction.Host.param.T.type:Class.AlignCols

template <typename TAlign>
struct Host<AlignCols<TAlign> >
{
    typedef TAlign Type;
};
template <typename TAlign>
struct Host<AlignCols<TAlign> const>
{
    typedef TAlign Type;
};

// ----------------------------------------------------------------------------
// Metafunction AlignColIterator
// ----------------------------------------------------------------------------

///.Metafunction.Iterator.param.T.type:Class.AlignCols

template <typename TAlign, typename TIteratorSpec>
struct Iterator<AlignCols<TAlign>, TIteratorSpec>
{
    typedef Iter<TAlign, AlignColIterator<void> > Type;
};
template <typename TAlign, typename TIteratorSpec>
struct Iterator<AlignCols<TAlign> const, TIteratorSpec>
{
    typedef Iter<TAlign, AlignColIterator<void> > Type;
};

// ----------------------------------------------------------------------------
// Metafunction Value
// ----------------------------------------------------------------------------

// Iterator is also used as Value

///.Metafunction.Value.param.T.type:Class.AlignCols

template <typename TAlign>
struct Value<AlignCols<TAlign> >:
    Iterator<AlignCols<TAlign>, Standard>
{};

template <typename TAlign>
struct Value<AlignCols<TAlign> const>:
    Iterator<AlignCols<TAlign> const, Standard>
{};

// ----------------------------------------------------------------------------
// Metafunction Size
// ----------------------------------------------------------------------------

///.Metafunction.Size.param.T.type:Class.AlignCols

template <typename TAlign>
struct Size<AlignCols<TAlign> >:
    Size<typename Row<TAlign>::Type>
{};

template <typename TAlign>
struct Size<AlignCols<TAlign> const>:
    Size<typename Row<TAlign const>::Type>
{};

// ----------------------------------------------------------------------------
// Metafunction Position
// ----------------------------------------------------------------------------

///.Metafunction.Position.param.T.type:Class.AlignCols

template <typename TAlign>
struct Position<AlignCols<TAlign> >:
    Position<typename Row<TAlign>::Type>
{};

template <typename TAlign>
struct Position<AlignCols<TAlign> const>:
    Position<typename Row<TAlign const>::Type>
{};

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function host()
// ----------------------------------------------------------------------------

///.Function.host.param.object.type:Class.AlignCols

template <typename TAlign>
inline typename Host<AlignCols<TAlign> >::Type &
host(AlignCols<TAlign> & me)
{
    SEQAN_ASSERT(me.data_align != 0);
    return *me.data_align;
}

template <typename TAlign>
inline typename Host<AlignCols<TAlign> const>::Type &
host(AlignCols<TAlign> const & me)
{
    SEQAN_ASSERT(me.data_align != 0);
    return *me.data_align;
}

// ----------------------------------------------------------------------------
// Function iter()
// ----------------------------------------------------------------------------

///.Function.iter.param.object.type:Class.AlignCols

template <typename TAlign, typename TPosition, typename TTag>
inline typename Iterator<AlignCols<TAlign>, Tag<TTag> const>::Type
iter(AlignCols<TAlign> & me,
     TPosition pos_,
     Tag<TTag> const)
{
    return typename Iterator<AlignCols<TAlign>, Tag<TTag> const>::Type(host(me), pos_);
}

template <typename TAlign, typename TPosition, typename TTag>
inline typename Iterator<AlignCols<TAlign> const, Tag<TTag> const>::Type
iter(AlignCols<TAlign> const & me,
     TPosition pos_,
     Tag<TTag> const)
{
    return typename Iterator<AlignCols<TAlign> const, Tag<TTag> const>::Type(host(me), pos_);
}

// ----------------------------------------------------------------------------
// Function value()
// ----------------------------------------------------------------------------

///.Function.value.param.container.type:Class.AlignCols

template <typename TAlign, typename TPosition>
inline typename Value<AlignCols<TAlign> >::Type
value(AlignCols<TAlign> & me,
      TPosition _pos)
{
    return iter(me, _pos);
}

template <typename TAlign, typename TPosition>
inline typename Value<AlignCols<TAlign> const>::Type
value(AlignCols<TAlign> const & me,
      TPosition _pos)
{
    return iter(me, _pos);
}

// ----------------------------------------------------------------------------
// Function beginPosition()
// ----------------------------------------------------------------------------

///.Function.beginPosition.param.object.type:Class.AlignCols

template <typename TAlignCols>
inline typename Position<TAlignCols>::Type
_beginPositionAlignCols(TAlignCols const & me)
{
    typedef typename Host<TAlignCols>::Type TAlign;
    typename Position<typename Rows<TAlign>::Type>::Type _i = length(rows(host(me)));

    if (!_i)
    {
        return 0;
    }

    --_i;
    typename Position<TAlignCols>::Type _pos = beginPosition(row(host(me), _i));

    while (_i > 0)
    {
        --_i;
        typename Position<TAlignCols>::Type _pos2 = beginPosition(row(host(me), _i));
        if (_pos2 < _pos)
        {
            _pos = _pos2;
        }
    }
    return _pos;
}

template <typename TAlign>
inline typename Position<AlignCols<TAlign> >::Type
beginPosition(AlignCols<TAlign> const & me)
{
    return _beginPositionAlignCols(me);
}

template <typename TAlign>
inline typename Position<AlignCols<TAlign> >::Type
beginPosition(AlignCols<TAlign> & me)
{
    return _beginPositionAlignCols(me);
}

// ----------------------------------------------------------------------------
// Function begin()
// ----------------------------------------------------------------------------

///.Function.begin.param.object.type:Class.AlignCols

template <typename TAlign, typename TTag>
inline typename Iterator<AlignCols<TAlign>, Tag<TTag> const>::Type
begin(AlignCols<TAlign> & me,
      Tag<TTag> const tag_)
{
    return iter(me, beginPosition(me), tag_);
}

template <typename TAlign, typename TTag>
inline typename Iterator<AlignCols<TAlign> const, Tag<TTag> const>::Type
begin(AlignCols<TAlign> const & me,
      Tag<TTag> const tag_)
{
    return iter(me, beginPosition(me), tag_);
}

// ----------------------------------------------------------------------------
// Function endPosition()
// ----------------------------------------------------------------------------

///.Function.endPosition.param.object.type:Class.AlignCols

template <typename TAlignCols>
inline typename Position<TAlignCols>::Type
_endPositionAlignCols(TAlignCols const & me)
{
    typedef typename Host<TAlignCols>::Type TAlign;

    typename Position<typename Rows<TAlign>::Type>::Type _i = length(rows(host(me)));
    typename Position<TAlignCols>::Type _pos = 0;

    while (_i > 0)
    {
        --_i;
        typename Position<TAlignCols>::Type _pos2 = endPosition(row(host(me), _i));
        if (_pos2 > _pos)
        {
            _pos = _pos2;
        }
    }
    return _pos;
}

template <typename TAlign>
inline typename Position<AlignCols<TAlign> >::Type
endPosition(AlignCols<TAlign> & me)
{
    return _endPositionAlignCols(me);
}

template <typename TAlign>
inline typename Position<AlignCols<TAlign> const>::Type
endPosition(AlignCols<TAlign> const & me)
{
    return _endPositionAlignCols(me);
}

// ----------------------------------------------------------------------------
// Function end()
// ----------------------------------------------------------------------------

///.Function.end.param.object.type:Class.AlignCols

template <typename TAlign, typename TTag>
inline typename Iterator<AlignCols<TAlign>, Tag<TTag> const>::Type
end(AlignCols<TAlign> & me,
    Tag<TTag> const tag_)
{
    return iter(me, endPosition(me), tag_);
}

template <typename TAlign, typename TTag>
inline typename Iterator<AlignCols<TAlign> const, Tag<TTag> const>::Type
end(AlignCols<TAlign> const & me,
    Tag<TTag> const tag_)
{
    return iter(me, endPosition(me), tag_);
}

// ----------------------------------------------------------------------------
// Function length()
// ----------------------------------------------------------------------------

template <typename TAlign>
inline typename Size<AlignCols<TAlign> >::Type
length(AlignCols<TAlign> const & me)
{
    return endPosition(me) - beginPosition(me);
}

// ----------------------------------------------------------------------------
// Function operator==()
// ----------------------------------------------------------------------------

template <typename TAlign>
inline bool
operator==(AlignCols<TAlign> const & left,
           AlignCols<TAlign> const & right)
{
    return left.data_align == right.data_align;
}

}  // namespace seqan

#endif  // #ifndef SEQAN_CORE_INCLUDE_SEQAN_ALIGN_ALIGN_COLS_H_
