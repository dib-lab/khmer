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
// Adaptions for pointer and arrays to SeqAn strings.
//
// TODO(holtgrew): Break out into adapt_pointer.h and adapt_array.h? The important main distinction is the fixed size at compile time.
// ==========================================================================

#ifndef SEQAN_SEQUENCE_ADAPT_ARRAY_POINTER_H_
#define SEQAN_SEQUENCE_ADAPT_ARRAY_POINTER_H_

namespace seqan {

// ===========================================================================
// Enums, Tags, Classes, Specializations
// ===========================================================================

/**
.Adaption.char array:
..summary:Zero terminated $char[]$ or $wchar_t[]$.
..remarks:Char arrays only support the Insist @Tag.Overflow Strategy.overflow strategy@.
*/

/**
.Adaption.char array.remarks:The default overflow strategy
(both @Metafunction.DefaultOverflowImplicit@ and @Metafunction.DefaultOverflowExplicit@)
for all operations on char arrays is @Tag.Overflow Strategy.insist@.
*/

// ===========================================================================
// Metafunctions
// ===========================================================================

template <typename TValue>
struct DefaultOverflowImplicit;

template <typename TValue>
struct DefaultOverflowImplicit<TValue *>
{
    typedef Insist Type;
};

template <typename TValue>
struct DefaultOverflowImplicit<TValue * const>
{
    typedef Insist Type;
};

template <typename TValue, size_t SIZE>
struct DefaultOverflowImplicit<TValue [SIZE]>
{
    typedef Insist Type;
};

template <typename TValue, size_t SIZE>
struct DefaultOverflowImplicit<TValue const [SIZE]>
{
    typedef Insist Type;
};

template <typename TValue>
struct DefaultOverflowExplicit;

template <typename TValue>
struct DefaultOverflowExplicit< TValue * >
{
    typedef Insist Type;
};
template <typename TValue>
struct DefaultOverflowExplicit< TValue * const>
{
    typedef Insist Type;
};
template <typename TValue, size_t SIZE>
struct DefaultOverflowExplicit< TValue [SIZE] >
{
    typedef Insist Type;
};
template <typename TValue, size_t SIZE>
struct DefaultOverflowExplicit< TValue const [SIZE] >
{
    typedef Insist Type;
};

///.Metafunction.IsContiguous.param.T.type:Adaption.char array
///.Metafunction.IsContiguous.class:Adaption.char array

template <typename TValue>
struct IsContiguous;

template <typename TValue>
struct IsContiguous< TValue * >
{
    typedef True Type;
    enum { VALUE = true };
};

template <typename TValue, size_t SIZE>
struct IsContiguous< TValue [SIZE] >
{
    typedef True Type;
    enum { VALUE = true };
};

template <typename TValue, size_t SIZE>
struct IsContiguous< TValue const [SIZE] >
{
    typedef True Type;
    enum { VALUE = true };
};

///.Metafunction.IsSequence.param.T.type:Adaption.char array
///.Metafunction.IsSequence.class:Adaption.char array

template <typename TValue>
struct IsSequence< TValue * >
{
    typedef True Type;
    enum { VALUE = true };
};
template <typename TValue, size_t SIZE>
struct IsSequence< TValue [SIZE] >
{
    typedef True Type;
    enum { VALUE = true };
};
template <typename TValue, size_t SIZE>
struct IsSequence< TValue const [SIZE] >
{
    typedef True Type;
    enum { VALUE = true };
};

// ----------------------------------------------------------------------------
// Metafunction Iterator
// ----------------------------------------------------------------------------

///.Metafunction.Iterator.param.T.type:Adaption.char array
///.Metafunction.Iterator.class:Adaption.char array

template <typename TValue>
struct Iterator<TValue *, Standard>
{
    typedef TValue * Type;
};

template <typename TValue>
struct Iterator<TValue * const, Standard>
{
    typedef TValue * Type;
};

template <typename TValue, size_t SIZE>
struct Iterator<TValue [SIZE], Standard>
        : Iterator<TValue *, Standard>
{
};

template <typename TValue, size_t SIZE>
struct Iterator<TValue const [SIZE], Standard>
        : Iterator<TValue const *, Standard>
{
};

template <typename TValue, size_t SIZE>
struct Iterator<TValue [SIZE], Rooted>
        : Iterator<TValue *, Rooted>
{
};

template <typename TValue, size_t SIZE>
struct Iterator<TValue const [SIZE], Rooted>
        : Iterator<TValue const *, Rooted>
{
};

// ===========================================================================
// Functions
// ===========================================================================

template <typename T>
inline typename Iterator<T *, typename DefaultGetIteratorSpec<T>::Type>::Type
begin(T * me)
{
SEQAN_CHECKPOINT
    return begin(me, typename DefaultGetIteratorSpec<T>::Type()) ;
}

///.Function.begin.param.object.type:Adaption.char array
///.Function.begin.class:Adaption.char array

template <typename TValue>
inline typename Iterator<TValue *, Standard>::Type
begin(TValue * me,
      Standard)
{
SEQAN_CHECKPOINT
    return me;
}

// TODO(holtgrew): Is the following still required since we dropped support for VC++ 2003?
//folgende Versionen wurde wegen seltsamer Phaenomene bei VC++ 2003 hinzugenommen
template <typename TValue>
inline typename Iterator<TValue const *, Standard>::Type
begin(TValue const * me,
      Standard)
{
SEQAN_CHECKPOINT
    return me;
}

template <typename TValue, typename TSpec>
inline typename Iterator<TValue *, Tag<TSpec> const>::Type
begin(TValue * me,
      Tag<TSpec> const)
{
SEQAN_CHECKPOINT
    typedef typename Iterator<TValue *, Tag<TSpec> const>::Type TIterator;
    return TIterator(me, begin(me, Standard()));
}

template <typename TValue, typename TSpec>
inline typename Iterator<TValue const *, Tag<TSpec> const>::Type
begin(TValue const * me,
      Tag<TSpec> const)
{
SEQAN_CHECKPOINT
    typedef typename Iterator<TValue const *, Tag<TSpec> const>::Type TIterator;
    return TIterator(me, begin(me, Standard()));
}

///.Function.end.param.object.type:Adaption.char array
///.Function.end.class:Adaption.char array

template <typename TValue>
inline typename Iterator<TValue *, Standard>::Type
end(TValue * me,
    Standard)
{
SEQAN_CHECKPOINT
    return begin(me, Standard()) + length(me);
}

//folgende Version wurde wegen eines seltsamen Phaenomens bei VC++ hinzugenommen
template <typename TValue>
inline typename Iterator<TValue const *, Standard>::Type
end(TValue const * me,
    Standard)
{
SEQAN_CHECKPOINT
    return begin(me, Standard()) + length(me);
}

template <typename TValue, typename TSpec>
inline typename Iterator<TValue *, Tag<TSpec> const>::Type
end(TValue * me,
      Tag<TSpec> const tag_)
{
SEQAN_CHECKPOINT
    return begin(me, tag_) + length(me);
}

template <typename TValue, typename TSpec>
inline typename Iterator<TValue const *, Tag<TSpec> const>::Type
end(TValue const * me,
      Tag<TSpec> const tag_)
{
SEQAN_CHECKPOINT
    return begin(me, tag_) + length(me);
}

///.Function.value.param.container.type:Adaption.char array
///.Function.value.class:Adaption.char array

template <typename TValue, typename TPos>
inline TValue &
value(TValue * me,
      TPos pos)
{
SEQAN_CHECKPOINT
    return me[pos];
}

template <typename TValue, typename TPos>
inline TValue const &
value(TValue const * me,
      TPos pos)
{
SEQAN_CHECKPOINT
    return me[pos];
}

///.Function.assignValue.param.container.type:Adaption.char array
///.Function.assignValue.class:Adaption.char array

template <typename TValue, typename TPos>
inline void
assignValue(TValue * me,
            TPos pos,
            TValue const & _value)
{
SEQAN_CHECKPOINT
    assign(value(me, pos), _value);
}

///.Function.moveValue.param.container.type:Adaption.char array
///.Function.moveValue.class:Adaption.char array

template <typename TValue, typename TPos>
inline void
moveValue(TValue * me,
          TPos pos,
          TValue const & _value)
{
SEQAN_CHECKPOINT
    move(value(me, pos), _value);
}

// Function atEnd for pointers / array iterators.

template <typename TValue>
inline bool
atEnd(TValue * pos)
{
SEQAN_CHECKPOINT
    return *pos == 0;
}

template <typename TValue>
inline bool
atEnd(TValue const * pos,
      TValue const * /*container*/)
{
SEQAN_CHECKPOINT
    return *pos == 0;
}

///.Function.length.param.object.type:Adaption.char array
///.Function.length.class:Adaption.char array

template <typename TValue>
inline size_t
length(TValue * me)
{
SEQAN_CHECKPOINT
    if (!me) return 0;
    TValue * it = me;
    TValue zero = TValue(0);
    while ( *it != zero) ++it;
    return it - me;
}

template <typename TValue>
inline size_t
length(TValue const * me)
{
SEQAN_CHECKPOINT
    if (!me) return 0;
    TValue const * it = me;
    TValue const zero = TValue();
    while ( *it != zero) ++it;
    return it - me;
}

inline size_t
length(char * me)
{
SEQAN_CHECKPOINT
    return ::std::strlen(me);
}

inline size_t
length(char const * me)
{
SEQAN_CHECKPOINT
    return ::std::strlen(me);
}

template <typename TValue>
inline void
_setLength(TValue * me,
           size_t new_length)
{
SEQAN_CHECKPOINT
    me[new_length] = 0;
}

///.Function.clear.param.object.type:Adaption.char array

template <typename TValue>
inline void
clear(TValue * me)
{
SEQAN_CHECKPOINT
    // TODO(holtgrew): Review this.
    //arrayDestruct(begin(me), length(me)); //??? Die Laengenbestimmung ist meistens nutzlos, braucht man sowieso nur fuer non-pod
    _setLength(me, 0);
}

///.Function.empty.param.object.type:Adaption.char array

template <typename TValue>
inline bool
empty(TValue * me)
{
SEQAN_CHECKPOINT
    return !me || (*me == TValue());
}

template<typename TValue, typename TExpand>
inline size_t
_clearSpace(TValue * me,
           size_t size,
           Tag<TExpand>)
{
SEQAN_CHECKPOINT
    return ClearSpaceStringBase_<Tag<TExpand> >::_clearSpace_(me, size);
}

template<typename TValue, typename TExpand>
inline size_t
_clearSpace(TValue * me,
           size_t size,
           size_t limit,
           Tag<TExpand>)
{
SEQAN_CHECKPOINT
    return ClearSpaceStringBase_<Tag<TExpand> >::_clearSpace_(me, size, limit);
}

template<typename TValue, typename TPosition, typename TExpand>
inline size_t
_clearSpace(TValue * me,
           size_t size,
           TPosition pos_begin,
           TPosition pos_end,
           Tag<TExpand>)
{
SEQAN_CHECKPOINT
    return ClearSpaceStringBase_<Tag<TExpand> >::_clearSpace_(me, size, pos_begin, pos_end);
}

template<typename TValue, typename TPosition, typename TExpand>
inline size_t
_clearSpace(TValue * me,
           size_t size,
           TPosition pos_begin,
           TPosition pos_end,
           size_t limit,
           Tag<TExpand>)
{
SEQAN_CHECKPOINT
    return ClearSpaceStringBase_<Tag<TExpand> >::_clearSpace_(me, size, pos_begin, pos_end, limit);
}

///.Function.assign.param.target.type:.Adaption.char array
///.Function.assign.param.source.type:.Adaption.char array
///.Function.assign.class:Adaption.char array

//overload of binary version for strings:

template<typename TTargetValue, typename TSource>
inline typename EnableIf<IsCharType<TTargetValue> >::Type
assign(TTargetValue * target,
       TSource & source)
{
SEQAN_CHECKPOINT
    typedef TTargetValue * TTarget;
    assign(target, source, typename DefaultOverflowImplicit<TTarget>::Type());
}
template<typename TTargetValue, typename TSource>
inline typename EnableIf<IsCharType<TTargetValue> >::Type
assign(TTargetValue * target,
       TSource const & source)
{
SEQAN_CHECKPOINT
    typedef TTargetValue * TTarget;
    assign(target, source, typename DefaultOverflowImplicit<TTarget>::Type());
}

template<typename TTargetValue, typename TSource, typename TExpand>
inline void
assign(TTargetValue * target,
       TSource const & source,
       Tag<TExpand>)
{
SEQAN_CHECKPOINT
    AssignString_<Tag<TExpand> >::assign_(target, source);
}

template<typename TTargetValue, typename TSource, typename TExpand>
inline void
assign(TTargetValue * target,
       TSource const & source,
       size_t limit,
       Tag<TExpand>)
{
SEQAN_CHECKPOINT
    AssignString_<Tag<TExpand> >::assign_(target, source, limit);
}

//____________________________________________________________________________
//this variant is a workaround for the "const array"-bug of VC++

template<typename TTargetValue, typename TSourceValue, typename TExpand>
inline void
assign(TTargetValue * target,
       TSourceValue const * source,
       Tag<TExpand>)
{
SEQAN_CHECKPOINT
    AssignString_<Tag<TExpand> >::assign_(target, source);
}

template<typename TTargetValue, typename TSourceValue, typename TExpand>
inline void
assign(TTargetValue * target,
       TSourceValue const * source,
       size_t limit,
       Tag<TExpand>)
{
SEQAN_CHECKPOINT
    AssignString_<Tag<TExpand> >::assign_(target, source, limit);
}

///.Function.move.param.target.type:Adaption.char array
///.Function.move.class:Adaption.char array

//overload of binary version for strings:

template<typename TTargetValue, typename TSource>
inline void
move(TTargetValue * & target,
     TSource & source)
{
SEQAN_CHECKPOINT
    target = source;
}
template<typename TTargetValue, typename TSource>
inline void
move(TTargetValue * & target,
     TSource const & source)
{
SEQAN_CHECKPOINT
    target = source;
}


//////////////////////////////////////////////////////////////////////////////
// append
//////////////////////////////////////////////////////////////////////////////

///.Function.append.param.target.type:.Adaption.char array
///.Function.append.param.source.type:.Adaption.char array
///.Function.append.class:Adaption.char array

template<typename TTargetValue, typename TSource, typename TExpand>
inline void
append(TTargetValue * target,
       TSource const & source,
       Tag<TExpand>)
{
SEQAN_CHECKPOINT
    AppendString_<Tag<TExpand> >::append_(target, source);
}

template<typename TTargetValue, typename TSource, typename TExpand>
inline void
append(TTargetValue * target,
       TSource const & source,
       size_t limit,
       Tag<TExpand>)
{
SEQAN_CHECKPOINT
    AppendString_<Tag<TExpand> >::append_(target, source, limit);
}

//____________________________________________________________________________
//this variant is a workaround for the "const array"-bug of VC++

template<typename TTargetValue, typename TSourceValue, typename TExpand>
inline void
append(TTargetValue * target,
       TSourceValue const * source,
       Tag<TExpand>)
{
SEQAN_CHECKPOINT
    AppendString_<Tag<TExpand> >::append_(target, source);
}

template<typename TTargetValue, typename TSourceValue, typename TExpand>
inline void
append(TTargetValue * target,
       TSourceValue const * source,
       size_t limit,
       Tag<TExpand>)
{
SEQAN_CHECKPOINT
    AppendString_<Tag<TExpand> >::append_(target, source, limit);
}

//////////////////////////////////////////////////////////////////////////////
// replace
//////////////////////////////////////////////////////////////////////////////

///.Function.replace.param.target.type:.Adaption.char array
///.Function.replace.param.source.type:.Adaption.char array

template<typename TTargetValue, typename TSource, typename TExpand>
inline void
replace(TTargetValue * target,
        size_t pos_begin,
        size_t pos_end,
        TSource const & source,
        Tag<TExpand>)
{
SEQAN_CHECKPOINT
    ReplaceString_<Tag<TExpand> >::replace_(target, pos_begin, pos_end, source);
}

template<typename TTargetValue, typename TSource, typename TExpand>
inline void
replace(TTargetValue * target,
        size_t pos_begin,
        size_t pos_end,
        TSource const & source,
        size_t limit,
        Tag<TExpand>)
{
SEQAN_CHECKPOINT
    ReplaceString_<Tag<TExpand> >::replace_(target, pos_begin, pos_end, source, limit);
}
//____________________________________________________________________________
//this variant is a workaround for the "const array"-bug of VC++

template<typename TTargetValue, typename TSourceValue, typename TExpand>
inline void
replace(TTargetValue * target,
        size_t pos_begin,
        size_t pos_end,
        TSourceValue const * source,
        Tag<TExpand>)
{
SEQAN_CHECKPOINT
    ReplaceString_<Tag<TExpand> >::replace_(target, pos_begin, pos_end, source);
}

template<typename TTargetValue, typename TSourceValue, typename TExpand>
inline void
replace(TTargetValue * target,
        size_t pos_begin,
        size_t pos_end,
        TSourceValue const * source,
        size_t limit,
        Tag<TExpand>)
{
SEQAN_CHECKPOINT
    ReplaceString_<Tag<TExpand> >::replace_(target, pos_begin, pos_end, source, limit);
}

//////////////////////////////////////////////////////////////////////////////
// handling of iterators as begin and and
/*
template<typename TTargetValue, typename TSource, typename TExpand>
inline void
replace(TTargetValue * target,
        typename Iterator<TTargetValue *, Rooted>::Type pos_begin,
        typename Iterator<TTargetValue *, Rooted>::Type pos_end,
        TSource const & source,
        Tag<TExpand> tag)
{
    replace(target, position(pos_begin), position(pos_end), source, tag);
}
template<typename TTargetValue, typename TSource, typename TExpand>
inline void
replace(TTargetValue * target,
        typename Iterator<TTargetValue *, Rooted>::Type pos_begin,
        typename Iterator<TTargetValue *, Rooted>::Type pos_end,
        TSource const & source,
        size_t limit,
        Tag<TExpand> tag)
{
    replace(target, position(pos_begin), position(pos_end), source, limit, tag);
}
*/
//////////////////////////////////////////////////////////////////////////////
///.Function.resize.param.object.type:Adaption.char array

template <typename TValue, typename TSize, typename TExpand>
inline size_t
resize(
    TValue * me,
    TSize new_length,
    Tag<TExpand>)
{
SEQAN_CHECKPOINT
    return _Resize_String<Tag<TExpand> >::resize_(me, new_length);
}

template <typename TValue, typename TSize, typename TExpand>
inline size_t
resize(
    TValue * me,
    TSize new_length,
    TValue const & val,
    Tag<TExpand>)
{
SEQAN_CHECKPOINT
    return _Resize_String<Tag<TExpand> >::resize_(me, new_length, val);
}

//////////////////////////////////////////////////////////////////////////////
// TODO(holtgrew): Review this problem, and document in ticket system.
//PROBLEM: ambiguitiy "pointer/iterator" and "c-style string"
//workaround: disable all operators
/*
template <typename TLeftValue, typename TRight >
TLeftValue const *
operator += (TLeftValue * left,
             TRight const & right)
{
SEQAN_CHECKPOINT
    append(left, right);
    return left;
}
*/
//////////////////////////////////////////////////////////////////////////////

template <typename TLeftValue, typename TRight >
inline bool
isEqual(TLeftValue * left,
        TRight const & right)
{
SEQAN_CHECKPOINT
    typename Comparator<TLeftValue *>::Type _lex(left, right);
    return isEqual(_lex);
}
/*
template <typename TLeftValue, typename TRight >
inline bool
operator == (TLeftValue * left,
            TRight const & right)
{
SEQAN_CHECKPOINT
    typename Comparator<TLeftValue *>::Type _lex(left, right);
    return isEqual(_lex);
}
*/
//////////////////////////////////////////////////////////////////////////////

template <typename TLeftValue, typename TRight >
inline bool
isNotEqual(TLeftValue * left,
           TRight const & right)
{
SEQAN_CHECKPOINT
    typename Comparator<TLeftValue *>::Type _lex(left, right);
    return isNotEqual(_lex);
}
/*
template <typename TLeftValue, typename TRight >
inline bool
operator != (TLeftValue * left,
             TRight const & right)
{
SEQAN_CHECKPOINT
    typename Comparator<TLeftValue *>::Type _lex(left, right);
    return isNotEqual(_lex);
}
*/

//////////////////////////////////////////////////////////////////////////////

template <typename TLeftValue, typename TRight>
inline bool
isLess(TLeftValue * left,
       TRight const & right)
{
SEQAN_CHECKPOINT
    return isLess(left, right, typename DefaultPrefixOrder<TLeftValue *>::Type());
}
/*
template <typename TLeftValue, typename TRight>
inline bool
operator < (TLeftValue * left,
            TRight const & right)
{
SEQAN_CHECKPOINT
    return isLess(left, right, typename DefaultPrefixOrder<TLeftValue *>::Type());
}
*/
//////////////////////////////////////////////////////////////////////////////

template <typename TLeftValue, typename TRight>
inline bool
isLessOrEqual(TLeftValue * left,
             TRight const & right)
{
SEQAN_CHECKPOINT
    return isLessOrEqual(left, right, typename DefaultPrefixOrder<TLeftValue *>::Type());
}
/*
template <typename TLeftValue, typename TRight>
inline bool
operator <= (TLeftValue * left,
             TRight const & right)
{
SEQAN_CHECKPOINT
    return isLessOrEqual(left, right, typename DefaultPrefixOrder<TLeftValue *>::Type());
}
*/
//////////////////////////////////////////////////////////////////////////////

template <typename TLeftValue, typename TRight>
inline bool
isGreater(TLeftValue * left,
        TRight const & right)
{
SEQAN_CHECKPOINT
    return isGreater(left, right, typename DefaultPrefixOrder<TLeftValue *>::Type());
}
/*
template <typename TLeftValue, typename TRight>
inline bool
operator > (TLeftValue * left,
        TRight const & right)
{
SEQAN_CHECKPOINT
    return isGreater(left, right, typename DefaultPrefixOrder<TLeftValue *>::Type());
}
*/
//////////////////////////////////////////////////////////////////////////////

template <typename TLeftValue, typename TRight>
inline bool
isGreaterOrEqual(TLeftValue * left,
        TRight const & right)
{
SEQAN_CHECKPOINT
    return isGreaterOrEqual(left, right, typename DefaultPrefixOrder<TLeftValue *>::Type());
}
/*
template <typename TLeftValue, typename TRight>
inline bool
operator >= (TLeftValue * left,
        TRight const & right)
{
SEQAN_CHECKPOINT
    return isGreaterOrEqual(left, right, typename DefaultPrefixOrder<TLeftValue *>::Type());
}
*/
//////////////////////////////////////////////////////////////////////////////

}  // namespace seqan

//____________________________________________________________________________

#endif  // #ifndef SEQAN_SEQUENCE_ADAPT_ARRAY_POINTER_H_
