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

#ifndef SEQAN_HEADER_MISC_DEQUEUE_H
#define SEQAN_HEADER_MISC_DEQUEUE_H

#include <algorithm>
#include <seqan/sequence.h>

//////////////////////////////////////////////////////////////////////////////

namespace SEQAN_NAMESPACE_MAIN
{

/*!
 * @class Deque
 * @implements ContainerConcept
 * @headerfile <seqan/basic.h>
 * @brief A double-ended queue implementation on top of a @link String @endlink.
 *
 * @signature template <typename TValue[, typename TSpec]>
 *            class Deque;
 *
 * @tparam TValue The value type.
 * @tparam TSpec  Specialization tag for the String.  Default: <tt>Alloc&lt;&gt;</tt>.
 */

template <typename TValue, typename TSpec = Alloc<> >
class Dequeue
{
public:
    typedef String<TValue, TSpec>                        TString;
    typedef typename Iterator<TString, Standard>::Type    TIter;

    String<TValue, TSpec> data_string;

    TIter data_begin;    // string beginning
    TIter data_end;        // string end

    TIter data_front;    // front fifo character
    TIter data_back;    // back fifo character
    bool data_empty;    // fifo is empty

//____________________________________________________________________________

public:
    inline Dequeue()
    {
        clear(*this);
    }

    // ----------------------------------------------------------------------
    // Subscription operators; have to be defined in class def.
    // ----------------------------------------------------------------------

    template <typename TPos>
    inline typename Reference<Dequeue>::Type
    operator[] (TPos pos)
    {
        SEQAN_CHECKPOINT;
        return value(*this, pos);
    }

    template <typename TPos>
    inline typename Reference<Dequeue const>::Type
    operator[] (TPos pos) const
    {
        SEQAN_CHECKPOINT;
        return value(*this, pos);
    }
};

//////////////////////////////////////////////////////////////////////////////

template <typename TValue, typename TSpec>
struct Reference<Dequeue<TValue, TSpec> >
{
    typedef typename Reference<TValue>::Type Type;
};

template <typename TValue, typename TSpec>
struct Reference<Dequeue<TValue, TSpec> const >
{
    typedef typename Reference<TValue>::Type Type;
};


//////////////////////////////////////////////////////////////////////////////
// Iterators
//////////////////////////////////////////////////////////////////////////////

template<typename TValue, typename TSpec>
struct Iterator<Dequeue<TValue, TSpec>, Standard>
{
    typedef Iter<Dequeue<TValue, TSpec>, PositionIterator> Type;
};

template<typename TValue, typename TSpec>
struct Iterator<Dequeue<TValue, TSpec> const, Standard>
{
    typedef Iter<Dequeue<TValue, TSpec> const, PositionIterator> Type;
};

template<typename TValue, typename TSpec>
struct Iterator<Dequeue<TValue, TSpec>, Rooted>
{
    typedef Iter<Dequeue<TValue, TSpec>, PositionIterator> Type;
};

template<typename TValue, typename TSpec>
struct Iterator<Dequeue<TValue, TSpec> const, Rooted>
{
    typedef Iter<Dequeue<TValue, TSpec> const, PositionIterator> Type;
};


//////////////////////////////////////////////////////////////////////////////

template <typename TValue, typename TSpec>
inline bool
empty(Dequeue<TValue, TSpec> const &me)
{
    return me.data_empty;
}

//////////////////////////////////////////////////////////////////////////////

template <typename TValue, typename TSpec>
inline void
clear(Dequeue<TValue, TSpec> &me)
{
    clear(me.data_string);
    me.data_begin = begin(me.data_string, Standard());
    me.data_end = end(me.data_string, Standard());

    me.data_front = me.data_back = me.data_begin;
    me.data_empty = true;
}

//////////////////////////////////////////////////////////////////////////////

template <typename TValue, typename TSpec, typename TPos>
inline TValue &
value(Dequeue<TValue, TSpec> &me, TPos pos)
{
    typedef typename Size<Dequeue<TValue, TSpec> >::Type TSize;
    TSize wrap = length(me.data_string) - (me.data_front - me.data_begin);

    if ((TSize)pos < wrap)
        return value(me.data_front + pos);
    else
        return value(me.data_begin + (pos - wrap));
}

template <typename TValue, typename TSpec, typename TPos>
inline TValue const &
value(Dequeue<TValue, TSpec> const &me, TPos pos)
{
    typedef typename Size<Dequeue<TValue, TSpec> >::Type TSize;
    TSize wrap = length(me.data_string) - (me.data_front - me.data_begin);

    if ((TSize)pos < wrap)
        return value(me.data_front + pos);
    else
        return value(me.data_begin + (pos - wrap));
}

//////////////////////////////////////////////////////////////////////////////

template <typename TValue, typename TSpec>
inline TValue &
front(Dequeue<TValue, TSpec> &me)
{
    return *me.data_front;
}

template <typename TValue, typename TSpec>
inline TValue const &
front(Dequeue<TValue, TSpec> const &me)
{
    return *me.data_front;
}

//////////////////////////////////////////////////////////////////////////////

template <typename TValue, typename TSpec>
inline TValue &
back(Dequeue<TValue, TSpec> &me)
{
    return *me.data_back;
}

template <typename TValue, typename TSpec>
inline TValue const &
back(Dequeue<TValue, TSpec> const &me)
{
    return *me.data_back;
}

//////////////////////////////////////////////////////////////////////////////

template <typename TValue, typename TSpec>
inline bool
popFront(Dequeue<TValue, TSpec> &me)
{
    if (me.data_empty) return false;

    if (me.data_front == me.data_back)
        me.data_empty = true;
    else
    {
        if (++me.data_front == me.data_end)
            me.data_front = me.data_begin;
    }

    return true;
}

template <typename TValue, typename TSpec>
inline bool
popBack(Dequeue<TValue, TSpec> &me)
{
    if (me.data_empty) return false;

    if (me.data_front == me.data_back)
        me.data_empty = true;
    else
    {
        if (me.data_back == me.data_begin)
            me.data_back = me.data_end;
        --me.data_back;
    }

    return true;
}

//////////////////////////////////////////////////////////////////////////////

template <typename TValue, typename TSpec>
inline void
pushFront(Dequeue<TValue, TSpec> &me, TValue const & _value)
{
    typedef typename Dequeue<TValue, TSpec>::TIter TIter;

    if (me.data_empty)
    {
        if (me.data_begin == me.data_end)
            reserve(me, computeGenerousCapacity(me, length(me.data_string) + 1), Generous());
        me.data_empty = false;
    }
    else
    {
        TIter new_front = me.data_front;
        if (new_front == me.data_begin)
            new_front = me.data_end;
        --new_front;

        if (new_front == me.data_back)
        {
            reserve(me, computeGenerousCapacity(me, length(me.data_string) + 1), Generous());

            if (me.data_front == me.data_begin)
                me.data_front = me.data_end;
            --me.data_front;
        } else
            me.data_front = new_front;
    }
    assign(*me.data_front, _value);
}

template <typename TValue, typename TSpec>
inline void
pushBack(Dequeue<TValue, TSpec> &me, TValue const & _value)
{
    typedef typename Dequeue<TValue, TSpec>::TIter TIter;

    if (me.data_empty)
    {
        if (me.data_begin == me.data_end)
            reserve(me, computeGenerousCapacity(me, length(me.data_string) + 1), Generous());
        me.data_empty = false;
    }
    else
    {
        TIter new_back = me.data_back;
        if (++new_back == me.data_end)
            new_back = me.data_begin;

        if (new_back == me.data_front)
        {
            reserve(me, computeGenerousCapacity(me, length(me.data_string) + 1), Generous());
            // in this case reserve adds new space behind data_back
            ++me.data_back;
        } else
            me.data_back = new_back;
    }
    assign(*me.data_back, _value);
}

//////////////////////////////////////////////////////////////////////////////

template <typename TValue, typename TSpec>
inline typename Size<Dequeue<TValue, TSpec> >::Type
length(Dequeue<TValue, TSpec> const &me)
{
    if (empty(me)) return 0;

    if (me.data_front <= me.data_back)
        return (me.data_back - me.data_front) + 1;
    else
        return (me.data_end - me.data_begin) - (me.data_front - me.data_back) + 1;
}

//////////////////////////////////////////////////////////////////////////////

template <typename TValue, typename TSpec, typename TSize_, typename TExpand>
inline typename Size<Dequeue<TValue, TSpec> >::Type
reserve(Dequeue<TValue, TSpec> &me, TSize_ new_capacity, Tag<TExpand> tag)
{
    typedef typename Size<Dequeue<TValue, TSpec> >::Type TSize;
//    std::cout << "resize to "<<new_capacity<<std::endl;
    TSize len = length(me);
    if (len < new_capacity && length(me.data_string) != new_capacity)
    {
        TSize pos_front = me.data_front - me.data_begin;
        TSize pos_back  = me.data_back  - me.data_begin;
        TSize new_freeSpace = new_capacity - len;

        if (pos_front <= pos_back)
        {
            // |empty|data|empty|
            // 0
            TSize freeSpace = length(me.data_string) - len;
            if (new_freeSpace > freeSpace)
                resize(me.data_string, new_capacity, tag);
            else
            {
                freeSpace -= new_freeSpace;    // reduce the free space by <freeSpace>
                if (pos_front >= freeSpace)
                {
                    resizeSpace(me.data_string, pos_front - freeSpace, (TSize)0, pos_front, tag);
                    pos_back -= freeSpace;
                    pos_front -= freeSpace;
                }
                else
                {
                    freeSpace -= pos_front;
                    resizeSpace(me.data_string, length(me.data_string) - freeSpace, pos_back + 1, length(me.data_string), tag);
                    resizeSpace(me.data_string, (TSize)0, (TSize)0, pos_front, tag);
                    pos_back -= pos_front;
                    pos_front = 0;
                }
            }
        }
        else
        {
            // |data|empty|data|
            // 0
            resizeSpace(me.data_string, new_freeSpace, pos_back + 1, pos_front, tag);
            pos_front += new_freeSpace;
        }

        me.data_begin = begin(me.data_string, Standard());
        me.data_end = end(me.data_string, Standard());
        me.data_front = me.data_begin + pos_front;
        me.data_back = me.data_begin + pos_back;
    }
    return length(me.data_string);
}

}

#endif

