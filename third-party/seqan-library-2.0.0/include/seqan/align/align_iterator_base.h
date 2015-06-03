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
//  Author: Andreas Gogol-Doering <andreas.doering@mdc-berlin.de>
// ==========================================================================

#ifndef SEQAN_HEADER_ALIGN_ITERATOR_BASE_H
#define SEQAN_HEADER_ALIGN_ITERATOR_BASE_H

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////
// Align Iterator for Gaps alignment
//////////////////////////////////////////////////////////////////////////////

// TODO(holtgrew): Extend class Iter?
/*!
 * @class AlignColIterator
 * @extends Iter
 * @headerfile <seqan/align.h>
 * @brief Iterator for alignment columns.
 *
 * @signature template <typename TAlign, typename TSpec>
 *            class Iter<TAlign, AlignColIterator<TSpec> >;
 *
 * @tparam TAlign Align object to iterate columns of.
 * @tparam TSpec  Tag for specializing the class further.
 */

template <typename TAlign, typename TSpec>
class Iter<TAlign, AlignColIterator<TSpec> >
{
public:
    typedef typename Rows<TAlign>::Type TRows;
    typedef typename Row<TAlign>::Type TRow;
    typedef typename Iterator<TRow, Standard>::Type TRowIterator;
    typedef typename Position<TRow>::Type TRowPosition;
    typedef String<TRowIterator> TIterators;

    TAlign * data_host;
    TIterators data_iterators;

public:
    Iter()
    {
SEQAN_CHECKPOINT
    }
    Iter(TAlign & _align):
        data_host(& _align)
    {
SEQAN_CHECKPOINT
        typename Position<TRows>::Type _i = length(rows(_align));
        resize(data_iterators, _i, Exact());
    }
    Iter(TAlign & _align, TRowPosition _pos):
        data_host(& _align)
    {
SEQAN_CHECKPOINT
        typename Position<TRows>::Type _i = length(rows(_align));
        resize(data_iterators, _i, Exact());

        while (_i > 0)
        {
            --_i;
            data_iterators[_i] = iter(row(_align, _i), _pos);
        }
    }
    Iter(Iter const & _other):
        data_host(_other.data_host),
        data_iterators(_other.data_iterators)
    {
    }
    ~Iter()
    {
SEQAN_CHECKPOINT
    }

    Iter const &
    operator = (Iter const & _other)
    {
SEQAN_CHECKPOINT
        data_host = _other.data_host;
        data_iterators = _other.data_iterators;
        return *this;
    }
//____________________________________________________________________________
};

//////////////////////////////////////////////////////////////////////////////

// TODO(holtgrew); Document as dox/hosted?

template <typename TAlign, typename TSpec>
inline TAlign &
host(Iter<TAlign, AlignColIterator<TSpec> > & me)
{
SEQAN_CHECKPOINT
    return *me.data_host;
}
template <typename TAlign, typename TSpec>
inline TAlign &
host(Iter<TAlign, AlignColIterator<TSpec> > const & me)
{
SEQAN_CHECKPOINT
    return *me.data_host;
}

//////////////////////////////////////////////////////////////////////////////

template <typename TAlign, typename TSpec>
inline void
setHost(Iter<TAlign, AlignColIterator<TSpec> > & me, TAlign & _host)
{
SEQAN_CHECKPOINT
    me.data_host = & _host;
}

//////////////////////////////////////////////////////////////////////////////

template <typename TAlign, typename TSpec>
inline typename Cols<TAlign>::Type
container(Iter<TAlign, AlignColIterator<TSpec> > & me)
{
SEQAN_CHECKPOINT
    return cols(*me.data_host);
}
template <typename TAlign, typename TSpec>
inline typename Cols<TAlign>::Type
container(Iter<TAlign, AlignColIterator<TSpec> > const & me)
{
SEQAN_CHECKPOINT
    return cols(*me.data_host);
}

//////////////////////////////////////////////////////////////////////////////

template <typename TAlign, typename TSpec>
inline void
goNext(Iter<TAlign, AlignColIterator<TSpec> > & me)
{
SEQAN_CHECKPOINT
    typedef typename Row<TAlign>::Type TRow;
    typedef typename Iterator<TRow, Standard>::Type TRowIterator;
    typedef String<TRowIterator> TIterators;
    typedef typename Iterator<TIterators, Standard>::Type TIteratorsIterator;

    TIteratorsIterator _it = begin(me.data_iterators);
    TIteratorsIterator _it_end = end(me.data_iterators);

    while (_it != _it_end)
    {
        goNext(*_it);
        ++_it;
    }
}
//____________________________________________________________________________

template <typename TAlign, typename TSpec>
inline Iter<TAlign, AlignColIterator<TSpec> > &
operator ++(Iter<TAlign, AlignColIterator<TSpec> > & me)
{
SEQAN_CHECKPOINT
    goNext(me);
    return me;
}
//____________________________________________________________________________

template <typename TAlign, typename TSpec>
inline Iter<TAlign, AlignColIterator<TSpec> >
operator ++(Iter<TAlign, AlignColIterator<TSpec> > & me, int)
{
SEQAN_CHECKPOINT
    Iter<TAlign, AlignColIterator<TSpec> > ret = me;
    goNext(me);
    return ret;
}

//////////////////////////////////////////////////////////////////////////////

template <typename TAlign, typename TSpec>
inline void
goPrevious(Iter<TAlign, AlignColIterator<TSpec> > & me)
{
SEQAN_CHECKPOINT
    typedef typename Row<TAlign>::Type TRow;
    typedef typename Iterator<TRow, Standard>::Type TRowIterator;
    typedef String<TRowIterator> TIterators;
    typedef typename Iterator<TIterators, Standard>::Type TIteratorsIterator;

    TIteratorsIterator _it = begin(me.data_iterators);
    TIteratorsIterator _it_end = end(me.data_iterators);

    while (_it != _it_end)
    {
        goPrevious(*_it);
        ++_it;
    }
}
//____________________________________________________________________________

template <typename TAlign, typename TSpec>
inline Iter<TAlign, AlignColIterator<TSpec> > &
operator --(Iter<TAlign, AlignColIterator<TSpec> > & me)
{
SEQAN_CHECKPOINT
    goPrevious(me);
    return me;
}
//____________________________________________________________________________

template <typename TAlign, typename TSpec>
inline Iter<TAlign, AlignColIterator<TSpec> >
operator --(Iter<TAlign, AlignColIterator<TSpec> > & me, int)
{
SEQAN_CHECKPOINT
    Iter<TAlign, AlignColIterator<TSpec> > ret = me;
    goPrevious(me);
    return ret;
}

//////////////////////////////////////////////////////////////////////////////

template <typename TAlign1, typename TAlign2, typename TSpec>
inline bool
operator ==(Iter<TAlign1, AlignColIterator<TSpec> > & _left,
            Iter<TAlign2, AlignColIterator<TSpec> > & _right)
{
SEQAN_CHECKPOINT
    return getValue(_left.data_iterators, 0) == getValue(_right.data_iterators, 0);
}
template <typename TAlign1, typename TAlign2, typename TSpec>
inline bool
operator ==(Iter<TAlign1, AlignColIterator<TSpec> > const & _left,
            Iter<TAlign2, AlignColIterator<TSpec> > & _right)
{
SEQAN_CHECKPOINT
    return value(_left.data_iterators, 0) == value(_right.data_iterators, 0);
}
template <typename TAlign1, typename TAlign2, typename TSpec>
inline bool
operator ==(Iter<TAlign1, AlignColIterator<TSpec> > & _left,
            Iter<TAlign2, AlignColIterator<TSpec> > const & _right)
{
SEQAN_CHECKPOINT
    return value(_left.data_iterators, 0) == value(_right.data_iterators, 0);
}
template <typename TAlign1, typename TAlign2, typename TSpec>
inline bool
operator ==(Iter<TAlign1, AlignColIterator<TSpec> > const & _left,
            Iter<TAlign2, AlignColIterator<TSpec> > const & _right)
{
SEQAN_CHECKPOINT
    return value(_left.data_iterators, 0) == value(_right.data_iterators, 0);
}

//////////////////////////////////////////////////////////////////////////////

template <typename TAlign1, typename TAlign2, typename TSpec>
inline bool
operator !=(Iter<TAlign1, AlignColIterator<TSpec> > & _left,
            Iter<TAlign2, AlignColIterator<TSpec> > & _right)
{
SEQAN_CHECKPOINT
    return value(_left.data_iterators, 0) != value(_right.data_iterators, 0);
}
template <typename TAlign1, typename TAlign2, typename TSpec>
inline bool
operator !=(Iter<TAlign1, AlignColIterator<TSpec> > const & _left,
            Iter<TAlign2, AlignColIterator<TSpec> > & _right)
{
SEQAN_CHECKPOINT
    return value(_left.data_iterators, 0) != value(_right.data_iterators, 0);
}
template <typename TAlign1, typename TAlign2, typename TSpec>
inline bool
operator !=(Iter<TAlign1, AlignColIterator<TSpec> > & _left,
            Iter<TAlign2, AlignColIterator<TSpec> > const & _right)
{
SEQAN_CHECKPOINT
    return value(_left.data_iterators, 0) != value(_right.data_iterators, 0);
}
template <typename TAlign1, typename TAlign2, typename TSpec>
inline bool
operator !=(Iter<TAlign1, AlignColIterator<TSpec> > const & _left,
            Iter<TAlign2, AlignColIterator<TSpec> > const & _right)
{
SEQAN_CHECKPOINT
    return value(_left.data_iterators, 0) != value(_right.data_iterators, 0);
}

//////////////////////////////////////////////////////////////////////////////

template <typename TAlign, typename TSpec, typename TPosition>
inline typename Reference<TAlign>::Type
value(Iter<TAlign, AlignColIterator<TSpec> > & me,
      TPosition pos_)
{
SEQAN_CHECKPOINT
    return value(me.data_iterators[pos_]);
}
template <typename TAlign, typename TSpec, typename TPosition>
inline typename Reference<TAlign>::Type
value(Iter<TAlign, AlignColIterator<TSpec> > const & me,
      TPosition pos_)
{
SEQAN_CHECKPOINT
    return value(me.data_iterators[pos_]);
}
//////////////////////////////////////////////////////////////////////////////

template <typename TAlign, typename TSpec, typename TPosition>
inline typename GetValue<TAlign>::Type
getValue(Iter<TAlign, AlignColIterator<TSpec> > & me,
         TPosition pos_)
{
SEQAN_CHECKPOINT
    return getValue(me.data_iterators[pos_]);
}
template <typename TAlign, typename TSpec, typename TPosition>
inline typename GetValue<TAlign>::Type
getValue(Iter<TAlign, AlignColIterator<TSpec> > const & me,
         TPosition pos_)
{
SEQAN_CHECKPOINT
    return getValue(me.data_iterators[pos_]);
}

//////////////////////////////////////////////////////////////////////////////

template <typename TAlign, typename TSpec, typename TPosition, typename TValue>
inline void
assignValue(Iter<TAlign, AlignColIterator<TSpec> > & me,
            TPosition pos_,
            TValue & val)
{
SEQAN_CHECKPOINT
    return assignValue(me.data_iterators[pos_], val);
}
template <typename TAlign, typename TSpec, typename TPosition, typename TValue>
inline void
assignValue(Iter<TAlign, AlignColIterator<TSpec> > & me,
            TPosition pos_,
            TValue const & val)
{
SEQAN_CHECKPOINT
    return assignValue(me.data_iterators[pos_], val);
}
template <typename TAlign, typename TSpec, typename TPosition, typename TValue>
inline void
assignValue(Iter<TAlign, AlignColIterator<TSpec> > const & me,
            TPosition pos_,
            TValue & val)
{
SEQAN_CHECKPOINT
    return assignValue(me.data_iterators[pos_], val);
}
template <typename TAlign, typename TSpec, typename TPosition, typename TValue>
inline void
assignValue(Iter<TAlign, AlignColIterator<TSpec> > const & me,
            TPosition pos_,
            TValue const & val)
{
SEQAN_CHECKPOINT
    return assignValue(me.data_iterators[pos_], val);
}

//////////////////////////////////////////////////////////////////////////////

template <typename TAlign, typename TSpec, typename TPosition, typename TValue>
inline void
moveValue(Iter<TAlign, AlignColIterator<TSpec> > & me,
          TPosition pos_,
          TValue & val)
{
SEQAN_CHECKPOINT
    return moveValue(me.data_iterators[pos_], val);
}
template <typename TAlign, typename TSpec, typename TPosition, typename TValue>
inline void
moveValue(Iter<TAlign, AlignColIterator<TSpec> > & me,
          TPosition pos_,
          TValue const & val)
{
SEQAN_CHECKPOINT
    return moveValue(me.data_iterators[pos_], val);
}
template <typename TAlign, typename TSpec, typename TPosition, typename TValue>
inline void
moveValue(Iter<TAlign, AlignColIterator<TSpec> > const & me,
          TPosition pos_,
          TValue & val)
{
SEQAN_CHECKPOINT
    return moveValue(me.data_iterators[pos_], val);
}
template <typename TAlign, typename TSpec, typename TPosition, typename TValue>
inline void
moveValue(Iter<TAlign, AlignColIterator<TSpec> > const & me,
          TPosition pos_,
          TValue const & val)
{
SEQAN_CHECKPOINT
    return moveValue(me.data_iterators[pos_], val);
}

//////////////////////////////////////////////////////////////////////////////

//??? TODO
//disabled since GapsIterator has no operator - and +
/*
template <typename TAlign, typename TSpec, typename TSize>
inline Iter<TAlign, AlignColIterator<TSpec> > &
operator +=(Iter<TAlign, AlignColIterator<TSpec> > & me,
            TSize size)
{
SEQAN_CHECKPOINT
    typedef typename Row<TAlign>::Type TRow;
    typedef typename Iterator<TRow>::Type TRowIterator;
    typedef String<TRowIterator> TIterators;
    typedef typename Iterator<TIterators>::Type TIteratorsIterator;

    TIteratorsIterator _it = begin(me.data_iterators);
    TIteratorsIterator _it_end = end(me.data_iterators);

    while (_it != _it_end)
    {
        *_it += size;
        ++_it;
    }
    return me;
}


//////////////////////////////////////////////////////////////////////////////

template <typename TAlign, typename TSpec, typename TSize>
inline Iter<TAlign, AlignColIterator<TSpec> >
operator +(Iter<TAlign, AlignColIterator<TSpec> > & me,
           TSize size)
{
SEQAN_CHECKPOINT
    Iter<TAlign, AlignColIterator<TSpec> > ret = me;
    me += size;
    return me;
}
template <typename TAlign, typename TSpec, typename TSize>
inline Iter<TAlign, AlignColIterator<TSpec> >
operator +(Iter<TAlign, AlignColIterator<TSpec> > const & me,
           TSize size)
{
SEQAN_CHECKPOINT
    Iter<TAlign, AlignColIterator<TSpec> > ret = me;
    me += size;
    return me;
}

//////////////////////////////////////////////////////////////////////////////

template <typename TAlign, typename TSpec, typename TSize>
inline Iter<TAlign, AlignColIterator<TSpec> > &
operator -=(Iter<TAlign, AlignColIterator<TSpec> > & me,
            TSize size)
{
SEQAN_CHECKPOINT
    typedef typename Row<TAlign>::Type TRow;
    typedef typename Iterator<TRow>::Type TRowIterator;
    typedef String<TRowIterator> TIterators;
    typedef typename Iterator<TIterators>::Type TIteratorsIterator;

    TIteratorsIterator _it = begin(me.data_iterators);
    TIteratorsIterator _it_end = end(me.data_iterators);

    while (_it != _it_end)
    {
        *_it -= size;
        ++_it;
    }
    return me;
}

//////////////////////////////////////////////////////////////////////////////

template <typename TAlign, typename TSpec, typename TSize>
inline Iter<TAlign, AlignColIterator<TSpec> >
operator -(Iter<TAlign, AlignColIterator<TSpec> > & me,
           TSize size)
{
SEQAN_CHECKPOINT
    Iter<TAlign, AlignColIterator<TSpec> > ret = me;
    me -= size;
    return me;
}
template <typename TAlign, typename TSpec, typename TSize>
inline Iter<TAlign, AlignColIterator<TSpec> >
operator -(Iter<TAlign, AlignColIterator<TSpec> > const & me,
           TSize size)
{
SEQAN_CHECKPOINT
    Iter<TAlign, AlignColIterator<TSpec> > ret = me;
    me -= size;
    return me;
}

//____________________________________________________________________________

template <typename TAlign, typename TSpec>
inline typename Difference<TAlign>::Type
operator -(Iter<TAlign, AlignColIterator<TSpec> > const & left,
           Iter<TAlign, AlignColIterator<TSpec> > const & right)
{
SEQAN_CHECKPOINT
    SEQAN_ASSERT_GT(length(left.data_iterators), 0u);
    SEQAN_ASSERT_GT(length(right.data_iterators), 0u);

    return (left.data_iterators[0] - right.data_iterators[0]);
}

//////////////////////////////////////////////////////////////////////////////

template <typename TAlign, typename TSpec>
inline typename Position<TAlign>::Type
position(Iter<TAlign, AlignColIterator<TSpec> > & me)
{
SEQAN_CHECKPOINT
    return position(me.data_iterators[0], row(host(me), 0));
}
template <typename TAlign, typename TSpec>
inline typename Position<TAlign>::Type
position(Iter<TAlign, AlignColIterator<TSpec> > const & me)
{
SEQAN_CHECKPOINT
    return position(me.data_iterators[0], row(host(me), 0));
}
*/
//////////////////////////////////////////////////////////////////////////////



//////////////////////////////////////////////////////////////////////////////

}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
