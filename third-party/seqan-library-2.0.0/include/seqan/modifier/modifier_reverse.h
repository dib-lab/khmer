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
// Author: Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>
// ==========================================================================

#ifndef SEQAN_HEADER_MODIFIER_REVERSE_H
#define SEQAN_HEADER_MODIFIER_REVERSE_H

namespace seqan
{

// ==========================================================================
// Forwards
// ==========================================================================

// ==========================================================================
// Classes, Enums, Typedefs
// ==========================================================================

// --------------------------------------------------------------------------
// Class ModReverse Iterator
// --------------------------------------------------------------------------

/*!
 * @class ModReverseIterator
 * @extends ModifiedIterator
 * @headerfile <seqan/modifier.h>
 * @brief Mirror the characters from begin to end.
 *
 * @signature template <typename THost>
 *            class ModifiedIterator<THost, ModReverse>;
 *
 * @tparam THost original iterator.
 */

/*!
 * @class ModReverseString
 * @extends ModifiedString
 * @headerfile <seqan/modifier.h>
 * @brief Mirror the characters from begin to end.
 *
 * @signature template <typename THost>
 *            class ModifiedString<THost, ModReverse>;
 *
 * @tparam THost original string.
 */

struct ModReverse_;
typedef Tag<ModReverse_> ModReverse;

// ==========================================================================
// Metafunctions
// ==========================================================================

// --------------------------------------------------------------------------
// Metafunction Cargo                           [ModReverse ModifiedIterator]
// --------------------------------------------------------------------------

template <typename THost>
struct Cargo<ModifiedIterator<THost, ModReverse> >
{
    typedef Cargo Type;        // to reduce namespace pollution
    bool _atEnd;

    Cargo() : _atEnd(false)
    {}
};

// --------------------------------------------------------------------------
// Metafunction Iterator                          [ModReverse ModifiedString]
// --------------------------------------------------------------------------

template <typename THost>
struct Iterator<ModifiedString<THost, ModReverse>, Standard>
{
    typedef ModifiedIterator<typename Iterator<THost const, Rooted>::Type, ModReverse> Type;
};

template <typename THost>
struct Iterator<ModifiedString<THost, ModReverse> const, Standard>
{
    typedef ModifiedIterator<typename Iterator<THost const, Rooted>::Type, ModReverse> Type;
};

// --------------------------------------------------------------------------
// Metafunction DefaultIteratorSpec               [ModReverse ModifiedString]
// --------------------------------------------------------------------------

template <typename THost>
struct DefaultIteratorSpec< ModifiedString<THost, ModReverse> >
{
    typedef Rooted Type;
};

// ==========================================================================
// Functions
// ==========================================================================

// --------------------------------------------------------------------------
// Function goNext()                            [ModReverse ModifiedIterator]
// --------------------------------------------------------------------------

template <typename THost>
inline void
goNext(ModifiedIterator<THost, ModReverse> & me)
{
    if (atBegin(host(me)))
        cargo(me)._atEnd = true;
    else
        goPrevious(host(me));
}

// --------------------------------------------------------------------------
// Function goPrevious()                        [ModReverse ModifiedIterator]
// --------------------------------------------------------------------------

template <typename THost>
inline void
goPrevious(ModifiedIterator<THost, ModReverse> & me)
{
    if (cargo(me)._atEnd)
        cargo(me)._atEnd = false;
    else
        goNext(host(me));
}

// --------------------------------------------------------------------------
// Function goEnd()                             [ModReverse ModifiedIterator]
// --------------------------------------------------------------------------

template <typename THost, typename TContainer>
inline void
goEnd(ModifiedIterator<THost, ModReverse> & me,
      TContainer & container)
{
    goBegin(host(me), host(container));
    cargo(me)._atEnd = true;
}

// --------------------------------------------------------------------------
// Function goBegin()                           [ModReverse ModifiedIterator]
// --------------------------------------------------------------------------

template <typename THost, typename TContainer>
inline void
goBegin(ModifiedIterator<THost, ModReverse> & me,
        TContainer & container)
{
    typedef typename Host<TContainer>::Type THostContainer;
    typename Parameter_<THostContainer>::Type hostContainer = host(container);
    goEnd(host(me), hostContainer);
    if (atBegin(host(me), hostContainer))
    {
        cargo(me)._atEnd = true;
    }
    else
    {
        cargo(me)._atEnd = false;
        goPrevious(host(me));
    }
}

//template <typename THost>
//inline void
//goBegin(ModifiedIterator<THost, ModReverse> & me)
//{
//    goEnd(host(me));
//    if (atBegin(host(me)))
//    {
//        cargo(me)._atEnd = true;
//    }
//    else
//    {
//        cargo(me)._atEnd = false;
//        goPrevious(host(me));
//    }
//}

// --------------------------------------------------------------------------
// Function operator+=()                        [ModReverse ModifiedIterator]
// --------------------------------------------------------------------------

template <typename THost, typename TDelta>
inline ModifiedIterator<THost, ModReverse> &
operator+=(ModifiedIterator<THost, ModReverse> & me, TDelta delta_)
{
    typedef ModifiedIterator<THost, ModReverse> TIterator;
    typedef typename Position<TIterator>::Type TPosition;
    TPosition delta = delta_;

    if (delta == 0)
    {
        return me;
    }
    if (delta > 0)
    {
        if (position(host(me)) < delta)
        {
            cargo(me)._atEnd = true;
            --delta;
        }
        host(me) -= delta;
    }
    else
    {
        if (cargo(me)._atEnd)
        {
            cargo(me)._atEnd = false;
            ++delta;
        }
        host(me) -= delta;
    }
    return me;
}

// --------------------------------------------------------------------------
// Function operator-=()                        [ModReverse ModifiedIterator]
// --------------------------------------------------------------------------

template <typename THost, typename TDelta>
inline ModifiedIterator<THost, ModReverse> &
operator-=(ModifiedIterator<THost, ModReverse> & me, TDelta delta)
{
    typedef typename Position<THost>::Type TPos;
    typedef typename MakeSigned<TPos>::Type TSignedPos;

    if (delta > 0)
    {
        if (cargo(me)._atEnd)
        {
            cargo(me)._atEnd = false;
            --delta;
        }
        host(me) += delta;
    }
    else
    {
        if ((TSignedPos)position(host(me)) < -(TSignedPos)delta)
        {
            cargo(me)._atEnd = true;
            ++delta;
        }
        host(me) -= -delta;
    }
    return me;
}

// --------------------------------------------------------------------------
// Function operator-()                         [ModReverse ModifiedIterator]
// --------------------------------------------------------------------------

template <typename THost>
inline typename Difference< ModifiedIterator<THost, ModReverse> >::Type
operator-(ModifiedIterator<THost, ModReverse> const & a,
          ModifiedIterator<THost, ModReverse> const & b)
{
    typename Difference< ModifiedIterator<THost, ModReverse> >::Type diff = host(b) - host(a);
    if (cargo(a)._atEnd)
        ++diff;
    if (cargo(b)._atEnd)
        --diff;
    return diff;
}

// --------------------------------------------------------------------------
// Function position()                          [ModReverse ModifiedIterator]
// --------------------------------------------------------------------------

template <typename THost>
inline typename Position<ModifiedIterator<THost, ModReverse> const>::Type
position(ModifiedIterator<THost, ModReverse> const & me)
{
    if (cargo(me)._atEnd)
        return length(container(host(me)));
    else
        return length(container(host(me))) - 1 - position(host(me));
}

template <typename THost, typename TContainer>
inline typename Position<ModifiedIterator<THost, ModReverse> const>::Type
position(ModifiedIterator<THost, ModReverse> const & me, TContainer const &cont)
{
    if (cargo(me)._atEnd)
        return length(cont);
    else
        return length(cont) - 1 - position(host(me), cont);
}

// --------------------------------------------------------------------------
// Function setPosition()                       [ModReverse ModifiedIterator]
// --------------------------------------------------------------------------

template <typename THost, typename TPosition>
inline void
setPosition(ModifiedIterator<THost, ModReverse> const & me, TPosition pos)
{
    setPosition(host(me), length(container(host(me))) - 1 - pos);
}

// --------------------------------------------------------------------------
// Function operator==()                        [ModReverse ModifiedIterator]
// --------------------------------------------------------------------------

template <typename THost>
inline bool
operator==(ModifiedIterator<THost, ModReverse> const & a,
           ModifiedIterator<THost, ModReverse> const & b)
{
    return cargo(a)._atEnd == cargo(b)._atEnd && host(a) == host(b);
}

// --------------------------------------------------------------------------
// Function operator<()                         [ModReverse ModifiedIterator]
// --------------------------------------------------------------------------

template <typename THost>
inline bool
operator<(ModifiedIterator<THost, ModReverse> const & a,
          ModifiedIterator<THost, ModReverse> const & b)
 {
    return (!cargo(a)._atEnd && cargo(b)._atEnd) ||
            (!cargo(a)._atEnd && !cargo(b)._atEnd && host(a) > host(b));
}

// --------------------------------------------------------------------------
// Function atEnd()                             [ModReverse ModifiedIterator]
// --------------------------------------------------------------------------

template <typename THost, typename TContainer>
inline bool
atBegin(ModifiedIterator<THost, ModReverse> const & me,
        TContainer const & container)
{
    return position(me, container) == 0;
}

template <typename THost>
inline bool
atBegin(ModifiedIterator<THost, ModReverse> const & me)
{
    return position(me) == 0;
}

// --------------------------------------------------------------------------
// Function atEnd()                             [ModReverse ModifiedIterator]
// --------------------------------------------------------------------------

template <typename THost, typename TContainer>
inline bool
atEnd(ModifiedIterator<THost, ModReverse> const & me,
      TContainer const & /*container*/)
{
            return cargo(me)._atEnd;
}

template <typename THost>
inline bool
atEnd(ModifiedIterator<THost, ModReverse> const & me)
{
            return cargo(me)._atEnd;
}

// --------------------------------------------------------------------------
// Function value()                               [ModReverse ModifiedString]
// --------------------------------------------------------------------------

template <typename THost, typename TPos>
inline typename Reference<ModifiedString<THost, ModReverse> >::Type
value(ModifiedString<THost, ModReverse> & me, TPos pos)
{
    return value(host(me), (length(host(me)) - 1) - pos);
}

template <typename THost, typename TPos>
inline typename Reference<ModifiedString<THost, ModReverse> const>::Type
value(ModifiedString<THost, ModReverse> const & me, TPos pos)
{
    return value(host(me), (length(host(me)) - 1) - pos);
}

// --------------------------------------------------------------------------
// Function begin()                               [ModReverse ModifiedString]
// --------------------------------------------------------------------------

template < typename THost>
inline typename Iterator< ModifiedString<THost, ModReverse> const >::Type
begin(ModifiedString<THost, ModReverse> const & me)
{
    typename Iterator< ModifiedString<THost, ModReverse> const >::Type temp_(end(host(me), Rooted()));
    _copyCargo(temp_, me);
    goNext(temp_);
    return temp_;
}

template < typename THost >
inline typename Iterator< ModifiedString<THost, ModReverse> >::Type
begin(ModifiedString<THost, ModReverse> & me)
{
    typename Iterator< ModifiedString<THost, ModReverse> >::Type temp_(end(host(me), Rooted()));
    _copyCargo(temp_, me);
    goNext(temp_);
    return temp_;
}

template < typename THost, typename TTagSpec >
inline typename Iterator< ModifiedString<THost, ModReverse> const, Tag<TTagSpec> const >::Type
begin(ModifiedString<THost, ModReverse> const & me, Tag<TTagSpec> const)
{
    typename Iterator< ModifiedString<THost, ModReverse> const, Tag<TTagSpec> const >::Type temp_(end(host(me), Rooted()));
    _copyCargo(temp_, me);
    goNext(temp_);
    return temp_;
}

template < typename THost, typename TTagSpec >
inline typename Iterator< ModifiedString<THost, ModReverse>, Tag<TTagSpec> const >::Type
begin(ModifiedString<THost, ModReverse> & me, Tag<TTagSpec> const)
{
    typedef typename Iterator< ModifiedString<THost, ModReverse>, Tag<TTagSpec> const >::Type TIterator;
    TIterator temp_(end(host(me), Rooted()));
    _copyCargo(temp_, me);
    goNext(temp_);
    return temp_;
}

// --------------------------------------------------------------------------
// Function end()                                 [ModReverse ModifiedString]
// --------------------------------------------------------------------------

template <typename THost >
inline typename Iterator<ModifiedString<THost, ModReverse> const >::Type
end(ModifiedString<THost, ModReverse> const & me)
{
    typename Iterator<ModifiedString<THost, ModReverse> const >::Type temp_(begin(host(me), Rooted()));
    _copyCargo(temp_, me);
    goNext(temp_);
    return temp_;
}

template <typename THost >
inline typename Iterator<ModifiedString<THost, ModReverse> >::Type
end(ModifiedString<THost, ModReverse> & me)
{
    typename Iterator<ModifiedString<THost, ModReverse> >::Type temp_(begin(host(me), Rooted()));
    _copyCargo(temp_, me);
    goNext(temp_);
    return temp_;
}

template <typename THost, typename TTagSpec >
inline typename Iterator<ModifiedString<THost, ModReverse> const, Tag<TTagSpec> const>::Type
end(ModifiedString<THost, ModReverse> const & me, Tag<TTagSpec> const)
{
    typename Iterator<ModifiedString<THost, ModReverse> const, Tag<TTagSpec> const >::Type temp_(begin(host(me), Rooted()));
    _copyCargo(temp_, me);
    goNext(temp_);
    return temp_;
}

template <typename THost, typename TTagSpec >
inline typename Iterator<ModifiedString<THost, ModReverse>, Tag<TTagSpec> const>::Type
end(ModifiedString<THost, ModReverse> & me, Tag<TTagSpec> const)
{
    typename Iterator<ModifiedString<THost, ModReverse>, Tag<TTagSpec> const >::Type temp_(begin(host(me), Rooted()));
    _copyCargo(temp_, me);
    goNext(temp_);
    return temp_;
}

// --------------------------------------------------------------------------
// Function reverse()
// --------------------------------------------------------------------------

/*!
 * @fn reverse
 * @headerfile <seqan/modifier.h>
 * @brief Reverse a container in-place.
 *
 * @signature void reverse(sequence);
 * @signature void reverse(stringSet);
 *
 * @param[in,out] sequence  The sequence to reverse in-place.
 * @param[in,out] stringSet The StringSet to reverse in-place.
 *
 * @section Remarks
 *
 * StringSet objects are reverse element-wise, i.e. the entries are reverse-complemented but their order itself
 * remains the same.
 */

template <typename TValue>
inline bool
_reverseDoSequential(TValue size)
{
    return size < (TValue)10000;
}

inline bool
_reverseDoSequential(unsigned char)
{
    return true;
}

template < typename TSequence, typename TParallelTag >
inline void
reverse(TSequence & sequence, Tag<TParallelTag> parallelTag)
{
    typedef typename Position<TSequence>::Type              TPos;
    typedef typename Iterator<TSequence, Standard>::Type    TIter;

    TIter itBeg = begin(sequence, Standard());
    TIter itEnd = end(sequence, Standard());
    Splitter<TPos> splitter(0, length(sequence) / 2, parallelTag);

    // disable multi-threading if sequence is too small
    // __uint64 cast is for 8bit size types for which comparison would be always true
    if (IsSameType<Tag<TParallelTag>, Parallel>::VALUE && _reverseDoSequential(length(sequence)))
        resize(splitter, 1);

    // (weese:) We have to cast the result of length to int to circumvent an internal gcc compiler error
    SEQAN_OMP_PRAGMA(parallel for num_threads((int)length(splitter)) schedule(static))
    for (int job = 0; job < (int)length(splitter); ++job)
    {
        TIter it1 = itBeg + splitter[job];
        TIter it2 = itEnd - (splitter[job] + 1);
        TIter it1End = itBeg + splitter[job + 1];
        for (; it1 != it1End; ++it1, --it2)
            std::swap(*it1, *it2);
    }
}

template < typename TSequence, typename TSpec, typename TParallelTag >
inline void
reverse(StringSet<TSequence, TSpec> & stringSet, Tag<TParallelTag>)
{
    typedef typename Position<StringSet<TSequence, TSpec> >::Type   TPos;
    typedef typename MakeSigned<TPos>::Type                         TSPos;

    TSPos seqCount = length(stringSet);
    SEQAN_OMP_PRAGMA(parallel for if (IsSameType<Tag<TParallelTag>, Parallel>::VALUE))
    for (TSPos seqNo = 0; seqNo < seqCount; ++seqNo)
        reverse(stringSet[seqNo], Serial());
}

template <typename TValue>
inline void
reverse(std::list<TValue> & list)
{
    list.reverse();
}

template < typename TText >
inline void
reverse(TText & text)
{
    reverse(text, Parallel());
}

// const variants for segments/modifiers

template < typename TText >
inline void
reverse(TText const & text)
{
    reverse(const_cast<TText &>(text), Parallel());
}

template < typename TText, typename TParallelTag >
inline void
reverse(TText const & text, Tag<TParallelTag> parallelTag)
{
    reverse(const_cast<TText &>(text), parallelTag);
}


// --------------------------------------------------------------------------
// Function reverseString()
// --------------------------------------------------------------------------

template <typename THost>
inline ModifiedString<THost, ModReverse>
reverseString(THost & host)
{
    return ModifiedString<THost, ModReverse>(host);
}

}

#endif
