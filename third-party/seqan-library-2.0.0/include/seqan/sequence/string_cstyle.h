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
// Author: Andreas Gogol-Doering <andreas.doering@mdc-berlin.de>
// ==========================================================================
// Implementation of the CStyle String specialization,
// ==========================================================================

#ifndef SEQAN_SEQUENCE_STRING_CSTYLE_H_
#define SEQAN_SEQUENCE_STRING_CSTYLE_H_

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// Forward for tag.
struct CStyle_;
typedef Tag<CStyle_> CStyle;

// Used in constructors.
template <typename TValue>
inline void
clear(String<TValue, CStyle> & me);

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// TODO(holtgrew): Put example into demo.

/*!
 * @class CStyleString CStyle String
 * @extends String
 * @headerfile <seqan/sequence.h>
 * @brief Allows adaption of strings to C-style strings.
 *
 * @signature template <typename TValue, typename TCStyle>
 *            class String<TValue, CStyle>;
 *
 * @tparam TValue The value type, that is the type of the items/characters stored in the string.Use @link Value @endlink
 *                to get the value type for a given class.
 *
 * Assigning a string <tt>TValue *</tt> to a CStyle String will not create a
 * copy of the string but just copy pointers.
 *
 * @section Remarks
 *
 * The purpose of this class is to access to the content of a sequence in a "zero terminated string" style.  This can be
 * useful if SeqAn classes has to be integrated in programs that use <tt>char</tt> arrays to store strings.  Instances
 * of <tt>String<TValue, CStyle></tt> can implicitely converted to a <tt>TValue *</tt> that points to a zero terminated
 * CStyle of <tt>TValue</tt>.
 *
 * The stored c-style string object can be set by constructors or assignment.  The content of a c-style string can
 * eighter be stored in a separate buffer, that is the source string is copied.  Or the buffer of the source string
 * itself is used instead, in this case the c-style string depends on the source string and gets invalid as soon as the
 * buffer of the source string is destroyed.
 *
 * Hence, this class is a kind of adaptor from an arbitrary SeqAn string to char arrays. Of course, the opposite way is
 * possible too.
 *
 * @section Examples
 *
 * @code{.cpp}
 * // Create a string str:
 * String<char> str = "this is a test string";
 *
 * // Create a c-style string object for str:
 * String<char, CStyle> cStyle = str;
 *
 * // Now use cStyle as char array:
 * strcmp(cStyle, "compare it to this string");
 * @endcode
 * If the c-style string is needed only temporarily, the function
 * <tt>toCString</tt> can be used:
 *
 * @code{.cpp}
 * String<char> str = "this is a test string";
 * strcmp(toCString(str), "compare it to this string");
 * @endcode
 */

struct CStyle_;
typedef Tag<CStyle_> CStyle;

#ifdef PLATFORM_WINDOWS_VS
#pragma warning( push )
// Disable warning C4521 locally (multiple copy constructors).
#pragma warning( disable: 4521 )
// Disable warning C4522 locally (multiple assignment operators).
#pragma warning( disable: 4522 )
#endif  // PLATFORM_WINDOWS_VS

template <typename TValue>
class String <TValue, CStyle >
{
public:
    TValue * data_begin;
    TValue * data_end;
    // If data_size > 0, then the buffer is owned by me and must be deallocated.
    size_t data_size;

    // TODO(holtgrew): Maybe better point to 0?
    static TValue EMPTY_STRING;

    String()
        : data_begin(&EMPTY_STRING),
          data_end(&EMPTY_STRING),
          data_size(0)
    {
        SEQAN_CHECKPOINT;
    }

    template <typename TString>
    String(TString & str) : data_begin(0), data_end(0), data_size(0)
    {
        SEQAN_CHECKPOINT;
        assign(*this, str);
    }

    template <typename TString>
    String(TString const & str) : data_begin(0), data_end(0), data_size(0)
    {
        SEQAN_CHECKPOINT;
        assign(*this, str);
    }

    String(String & str) : data_begin(0), data_end(0), data_size(0)
    {
        SEQAN_CHECKPOINT;
        assign(*this, str);
    }

    String(String const & str) : data_begin(0), data_end(0), data_size(0)
    {
        SEQAN_CHECKPOINT;
        assign(*this, str);
    }

    String(TValue * str)
        : data_begin(str),
          data_end(end(str)),
          data_size(0)
    {
        SEQAN_CHECKPOINT;
    }

    ~String()
    {
        SEQAN_CHECKPOINT;
        clear(*this);
    }

    template <typename TString>
    inline
    String & operator=(TString & str)
    {
        SEQAN_CHECKPOINT;
        assign(*this, str);
        return *this;
    }

    template <typename TString>
    inline
    String & operator=(TString const & str)
    {
        SEQAN_CHECKPOINT;
        assign(*this, str);
        return *this;
    }

    inline
    String & operator=(String & str)
    {
        SEQAN_CHECKPOINT;
        assign(*this, str);
        return *this;
    }

    inline
    String & operator=(String const & str)
    {
        SEQAN_CHECKPOINT;
        assign(*this, str);
        return *this;
    }

    inline
    operator TValue * ()
    {
        SEQAN_CHECKPOINT;
        return data_begin;
    }

    inline
    operator TValue const * () const
    {
        SEQAN_CHECKPOINT;
        return data_begin;
    }
};

// Define the static member
template <typename TValue>
TValue String<TValue, CStyle >::EMPTY_STRING = TValue();

#ifdef PLATFORM_WINDOWS_VS
// Reset warning state to previous one for C4521, C4522.
#pragma warning( pop )
#endif  // PLATFORM_WINDOWS_VS

// ============================================================================
// Metafunctions
// ============================================================================

// --------------------------------------------------------------------------
// Metafunction DefaultOverflowImplicit
// --------------------------------------------------------------------------

template <typename TValue>
struct DefaultOverflowImplicit<String<TValue, CStyle> >
{
    typedef Exact Type;
};

// --------------------------------------------------------------------------
// Metafunction IsContiguous
// --------------------------------------------------------------------------

template <typename TValue>
struct IsContiguous< String<TValue, CStyle > >
{
    typedef True Type;
    enum { VALUE = true };
};

// ============================================================================
// Functions
// ============================================================================

// --------------------------------------------------------------------------
// Function move()
// --------------------------------------------------------------------------

template <typename TValue>
inline void
move(String<TValue, CStyle> & target,
     String<TValue, CStyle> & source)
{
    SEQAN_CHECKPOINT;
    clear(target);

    target.data_begin = source.data_begin;
    target.data_end = source.data_end;
    target.data_size = source.data_size;

    source.data_begin = 0;
    source.data_end = 0;
    source.data_size = 0;
}

template <typename TValue>
inline void
move(String<TValue, CStyle> & target,
     String<TValue, CStyle> const & source)
{
    SEQAN_CHECKPOINT;
    move(target, const_cast<String<TValue, CStyle> &>(source));
}

// --------------------------------------------------------------------------
// Function begin()
// --------------------------------------------------------------------------

template <typename TValue>
inline typename Iterator<String<TValue, CStyle >, Standard>::Type
begin(String <TValue, CStyle > & me,
      Standard const &)
{
    SEQAN_CHECKPOINT;
    return me.data_begin;
}

template <typename TValue>
inline typename Iterator<String<TValue, CStyle > const, Standard>::Type
begin(String <TValue, CStyle > const & me,
      Standard const &)
{
    SEQAN_CHECKPOINT;
    return me.data_begin;
}

// --------------------------------------------------------------------------
// Internal Function _setBegin()
// --------------------------------------------------------------------------

template <typename TValue, typename TValue2>
inline void
_setBegin(String <TValue, CStyle > & me, TValue2 new_begin)
{
    SEQAN_CHECKPOINT;
    me.data_begin = new_begin;
}

// --------------------------------------------------------------------------
// Function end()
// --------------------------------------------------------------------------

template <typename TValue>
inline typename Iterator<String <TValue, CStyle >, Standard>::Type
end(String <TValue, CStyle > & me,
    Standard const &)
{
    SEQAN_CHECKPOINT;
    return me.data_end;
}

template <typename TValue>
inline typename Iterator<String <TValue, CStyle > const, Standard>::Type
end(String <TValue, CStyle > const & me,
    Standard const &)
{
    SEQAN_CHECKPOINT;
    return me.data_end;
}

// --------------------------------------------------------------------------
// Internal Function _setEnd()
// --------------------------------------------------------------------------

template <typename TValue, typename TValue2>
inline void
_setEnd(String <TValue, CStyle > & me, TValue2 new_end)
{
    SEQAN_CHECKPOINT;
    me.data_end = new_end;
    if (new_end != NULL)
        *new_end = TValue(); //??? ist das wirklich sinnvoll fuer typen, die weder char noch wchar_t sind?
}

// --------------------------------------------------------------------------
// Function capacity()
// --------------------------------------------------------------------------

template <typename TValue>
inline size_t
capacity(String <TValue, CStyle > const & me)
{
    SEQAN_CHECKPOINT;
    if (me.data_size) return me.data_size -1;
    else return me.data_end - me.data_begin;
}

// --------------------------------------------------------------------------
// Internal Function _reallocateStorage()
// --------------------------------------------------------------------------

//this function works also for dependent buffers
template <typename TValue>
inline TValue *
_reallocateStorage(
    String <TValue, CStyle > & me,
    size_t new_capacity,
    Exact)
{
    SEQAN_CHECKPOINT;
    TValue * _returnValue;
    if (me.data_size)
    {//dependent
        _returnValue = me.data_begin;
    }
    else
    {//not dependent
        _returnValue = 0;
    }

    me.data_size = new_capacity + 1; //+1 for zero termination
    allocate(me, me.data_begin, me.data_size, TagAllocateStorage());
    return _returnValue;
}

// --------------------------------------------------------------------------
// Internal Function _deallocateStorage()
// --------------------------------------------------------------------------

template <typename TValue>
inline void
_deallocateStorage(
    String <TValue, CStyle > & me,
    TValue * ptr,
    size_t capacity)
{
    SEQAN_CHECKPOINT;
    size_t size = capacity + 1;
    deallocate(me, ptr, size, TagAllocateStorage());
}

// --------------------------------------------------------------------------
// Function dependent()
// --------------------------------------------------------------------------

template <typename TValue>
inline bool
dependent(String <TValue, CStyle > & me)
{
    SEQAN_CHECKPOINT;
    return (me.data_size == 0);
}

// --------------------------------------------------------------------------
// Function assign()
// --------------------------------------------------------------------------

//special implementation for char array sources
template <typename TValue>
inline void
assign(String <TValue, CStyle > & target,
       TValue * source)
{
    SEQAN_CHECKPOINT;
    clear(target);
    target.data_begin = source;
    target.data_end = end(source);
}

// --------------------------------------------------------------------------
// Function assign()
// --------------------------------------------------------------------------

template <typename TTargetValue, typename TSource, typename TExpand>
inline void
assign(String<TTargetValue, CStyle> & target,
       TSource & source,
       Tag<TExpand> tag)
{
    SEQAN_CHECKPOINT;
    create(target, source, tag);
}

template <typename TTargetValue, typename TSource, typename TExpand>
inline void
assign(String<TTargetValue, CStyle> & target,
       TSource const & source,
       Tag<TExpand> tag)
{
    SEQAN_CHECKPOINT;
    create(target, source, tag);
}

template <typename TTargetValue, typename TSource, typename TSize, typename TExpand>
inline void
assign(String<TTargetValue, CStyle> & target,
       TSource & source,
       TSize /*limit*/,
       Tag<TExpand> tag)
{
    SEQAN_CHECKPOINT;
    create(target, source, tag);
}

template <typename TTargetValue, typename TSource, typename TSize, typename TExpand>
inline void
assign(String<TTargetValue, CStyle> & target,
       TSource const & source,
       TSize limit,
       Tag<TExpand> tag)
{
    SEQAN_CHECKPOINT;
    create(target, source, limit, tag);
}

// TODO(holtgrew): Still needed with dropped VC++ 2003 support?
//this variant is a workaround for the "const array"-bug of VC++

template <typename TTargetValue, typename TSourceValue, typename TExpand>
inline void
assign(String<TTargetValue, CStyle> & target,
       TSourceValue const * source,
       Tag<TExpand> tag)
{
    SEQAN_CHECKPOINT;
    create(target, source, tag);
}

template <typename TTargetValue, typename TSourceValue, typename TSize, typename TExpand>
inline void
assign(String<TTargetValue, CStyle> & target,
       TSourceValue const * source,
       TSize limit,
       Tag<TExpand> tag)
{
    SEQAN_CHECKPOINT;
    create(target, source, limit, tag);
}

//If source is non-const String, then there could be the possibility
//to use the source buffer

template <typename TExpand, bool IS_CONTIGUOUS>
struct AssignStringToStringArray_;

template <typename TExpand>
struct AssignStringToStringArray_<TExpand, true>
{
    template <typename TValue, typename TSourceSpec>
    static inline void
    assign_(String<TValue, CStyle> & target,
        String<TValue, TSourceSpec> & source)
    {
        if (capacity(source) > length(source))
        {//use source's buffer
    SEQAN_CHECKPOINT;
            clear(target);
            _setBegin(target, begin(source));
            _setEnd(target, end(source));
        }
        else
        {
            create(target, source, TExpand());
        }
    }

//special treatment of char:
//_computeSizeForCapacity is specialized for char such that there
//is enough place for the zero termination

    template <typename TSourceSpec>
    static inline void
    assign_(String<char, CStyle> & target,
        String<char, TSourceSpec> & source)
    {
    SEQAN_CHECKPOINT;
        clear(target);
        typedef String<char, CStyle> TTarget;
        typedef typename Iterator<TTarget>::Type TIterator;
        _setBegin(target, TIterator(begin(source)));
        _setEnd(target, TIterator(end(source)));
    }
};

template <typename TExpand>
struct AssignStringToStringArray_<TExpand, false>
{
    template <typename TValue, typename TSourceSpec>
    static inline void
    assign_(String<TValue, CStyle> & target,
        String<TValue, TSourceSpec> & source)
    {
    SEQAN_CHECKPOINT;
        create(target, source, TExpand());
    }
};

template <typename TValue, typename TSourceSpec, typename TExpand>
inline void
assign(String<TValue, CStyle> & target,
    String<TValue, TSourceSpec> & source,
    Tag<TExpand>)
{
    typedef String<TValue, TSourceSpec> TSource;
    AssignStringToStringArray_<Tag<TExpand>, IsContiguous<TSource>::VALUE>::assign_(target, source);
}

// --------------------------------------------------------------------------
// Function clear()
// --------------------------------------------------------------------------

template <typename TValue>
inline void
clear(String<TValue, CStyle> & me)
{
    if (me.data_size)
    {
        SEQAN_CHECKPOINT
        //          arrayDestruct(me, length(me));
        deallocate(me, me.data_begin, me.data_size);
        me.data_size = 0;
    }
    me.data_begin = me.data_end = &me.EMPTY_STRING;
}

// --------------------------------------------------------------------------
// Function create()
// --------------------------------------------------------------------------

//see basic_holder
template <typename TExpand>
struct CreateArrayStringExpand_
{
    template <typename TTarget, typename TSource>
    static inline void
    create_(TTarget & target,
        TSource & source)
    {
        typename Size<TTarget>::Type source_length = length(source);
        if (dependent(target) || (capacity(target) < source_length))
        {
    SEQAN_CHECKPOINT;
            typename Size<TTarget>::Type old_target_capacity = capacity(target);
            typename Value<TTarget>::Type * buf = _reallocateStorage(target, source_length, TExpand());
            if (buf)
            {
                _deallocateStorage(target, buf, old_target_capacity);
            }
        }
        if (length(source) > 0)
        {
            assignValue(begin(target, Standard()), 0); //set target length to 0
            assign(begin(target, Standard()), source, Insist());
            typedef typename Iterator<TTarget>::Type TTargetIterator;
            _setEnd(target, TTargetIterator( begin(target) + source_length));
        }
    }

    template <typename TTarget, typename TSource, typename TLimit>
    static inline void
    create_(TTarget & target,
        TSource & source,
        TLimit limit)
    {
        typename Size<TTarget>::Type copy_length = length(source);
        if (limit < copy_length)
        {
            copy_length = limit;
        }
        if (dependent(target) || (capacity(target) < copy_length))
        {
    SEQAN_CHECKPOINT;
            typename Size<TTarget>::Type old_target_capacity = capacity(target);
            TTarget * buf = _reallocateStorage(target, copy_length, TExpand());
            if (buf)
            {
                _deallocateStorage(target, buf, old_target_capacity);
            }
        }
        assign(begin(target, Standard()), source, copy_length, Insist());
        _setEnd(target, begin(target, Standard()) + copy_length);
    }
};

template <typename TExpand>
struct CreateArrayString_;

template <>
struct CreateArrayString_<Insist>
{
    template <typename TTarget, typename TSource>
    static inline void
    create_(TTarget & target,
        TSource & source)
    {
    SEQAN_CHECKPOINT;
        typename Size<TTarget>::Type source_length = length(source);
        if (dependent(target))
            _reallocateStorage(target, source_length, Exact());
        assign(begin(target, Standard()), source, source_length, Insist());
        _setEnd(target, begin(target, Standard()) + source_length);
    }

    template <typename TTarget, typename TSource, typename TSize>
    static inline void
    create_(TTarget & target,
        TSource & source,
        TSize limit)
    {
    SEQAN_CHECKPOINT;
        typename Size<TTarget>::Type copy_size = length(source);
        if (limit < copy_size)
        {
            copy_size = limit;
        }
        if (dependent(target))
            _reallocateStorage(target, copy_size, Exact());
        assign(begin(target, Standard()), source, copy_size, Insist());
        _setEnd(target, begin(target, Standard()) + copy_size);
    }
};

template <>
struct CreateArrayString_<Limit>
{
    template <typename TTarget, typename TSource>
    static inline void
    create_(TTarget & target,
        TSource & source)
    {
    SEQAN_CHECKPOINT;
        CreateArrayString_<Insist>::create_(target, source, capacity(target));
    }

    template <typename TTarget, typename TSource, typename TSize>
    static inline void
    create_(TTarget & target,
        TSource & source,
        TSize & limit)
    {
    SEQAN_CHECKPOINT;
        typename Size<TTarget>::Type copy_size = capacity(target);
        if (copy_size > limit)
        {
            copy_size = limit;
        }
        CreateArrayString_<Insist>::create_(target, source, copy_size);
    }
};

template <>
struct CreateArrayString_<Exact>:
    CreateArrayStringExpand_<Exact>
{
};

template <>
struct CreateArrayString_<Generous>:
    CreateArrayStringExpand_<Generous>
{
};

template <typename TTargetValue, typename TSource>
inline void
create(String<TTargetValue, CStyle> & target,
       TSource & source)
{
    SEQAN_CHECKPOINT;
    typedef String<TTargetValue, CStyle> TTarget;
    create(target, source, typename DefaultOverflowImplicit<TTarget>::Type());
}

template <typename TTargetValue, typename TSource, typename TSize>
inline void
create(String<TTargetValue, CStyle> & target,
       TSource & source,
       TSize limit)
{
    SEQAN_CHECKPOINT;
    typedef String<TTargetValue, CStyle> TTarget;
    create(target, source, limit, typename DefaultOverflowImplicit<TTarget>::Type());
}

template <typename TTargetValue, typename TSource, typename TExpand>
inline void
create(String<TTargetValue, CStyle> & target,
       TSource & source,
       Tag<TExpand>)
{
    SEQAN_CHECKPOINT;
    CreateArrayString_<Tag<TExpand> >::create_(target, source);
}

template <typename TTargetValue, typename TSource, typename TSize, typename TExpand>
inline void
create(String<TTargetValue, CStyle> & target,
       TSource & source,
       TSize limit,
       Tag<TExpand>)
{
    SEQAN_CHECKPOINT;
    CreateArrayString_<Tag<TExpand> >::create_(target, source, limit);
}

template <typename TTargetValue, typename TSource, typename TExpand>
inline void
create(String<TTargetValue, CStyle> & target,
       TSource const & source,
       Tag<TExpand>)
{
    SEQAN_CHECKPOINT;
    CreateArrayString_<Tag<TExpand> >::create_(target, source);
}

template <typename TTargetValue, typename TSource, typename TSize, typename TExpand>
inline void
create(String<TTargetValue, CStyle> & target,
       TSource const & source,
       TSize limit,
       Tag<TExpand>)
{
    SEQAN_CHECKPOINT;
    CreateArrayString_<Tag<TExpand> >::create_(target, source, limit);
}

//this variant is a workaround for the "const array"-bug of VC++

template <typename TTargetValue, typename TSourceValue, typename TExpand>
inline void
create(String<TTargetValue, CStyle> & target,
       TSourceValue const * source,
       Tag<TExpand>)
{
    SEQAN_CHECKPOINT;
    CreateArrayString_<Tag<TExpand> >::create_(target, source);
}

template <typename TTargetValue, typename TSourceValue, typename TSize, typename TExpand>
inline void
create(String<TTargetValue, CStyle> & target,
       TSourceValue const * source,
       TSize limit,
       Tag<TExpand>)
{
    SEQAN_CHECKPOINT;
    CreateArrayString_<Tag<TExpand> >::create_(target, source, limit);
}

// --------------------------------------------------------------------------
// Function toCString()
// --------------------------------------------------------------------------

/*!
 * @fn String#toCString
 * @brief Access sequence as c-style string.
 * @signature TValue* toCString(seq)
 * @param seq The sequence to be accessed. Type: @link String @endlink
 * @return TValue* For strings that store their elements in a contiguous block (see @link IsContiguous @endlink)
                   a pointer to first element of $object$ is returned.
 * @section Remarks
 *
 * If the alphabet of $object$ is $char$ or $wchar_t$ the return value is a c-style string representing the
 * contents of <tt>object<tt/>.
 *
 * Calling this function for non-contiguous containers will raise a compilation error.  To create
 * c-style strings for non-contiguous strings or strings with different alphabets, use a @link CStyleString @endlink as an
 * intermediate.
 */

template <typename TValue>
inline TValue *
toCString(TValue * me)
{
    SEQAN_CHECKPOINT;
    return me;
}

template <typename TValue>
inline TValue const *
toCString(TValue const * me)
{
    SEQAN_CHECKPOINT;
    return me;
}

template <typename TValue>
inline TValue *
toCString(String<TValue, CStyle> & me)
{
    SEQAN_CHECKPOINT;
    return me;
}

template <typename TValue>
inline TValue const *
toCString(String<TValue, CStyle> const & me)
{
    SEQAN_CHECKPOINT;
    return me;
}

template <typename TValue, typename TSpec>
inline TValue *
_toCStringImpl(String<TValue, TSpec> & me, True)
{
    SEQAN_CHECKPOINT;
    typename Size< String<TValue, TSpec> >::Type len = length(me);
    if (len >= capacity(me))
        reserve(me, len + 1);
    if (end(me) != NULL)
        *end(me) = TValue();
    return begin(me);
}

// You called toCString with non-contiguous strings.
// Convert to String<..., CStyle> and try again!
template <typename TValue, typename TSpec>
inline TValue *
_toCStringImpl(String<TValue, TSpec> & me, False);


template <typename TValue, typename TSpec>
inline TValue *
toCString(String<TValue, TSpec> & me)
{
    SEQAN_CHECKPOINT;
    return _toCStringImpl(me, typename IsContiguous<String<TValue, TSpec> >::Type());
}

template <typename TValue, typename TSpec>
inline TValue *
toCString(String<TValue, TSpec> const & me)
{
    SEQAN_CHECKPOINT;
    return toCString(const_cast<String<TValue, TSpec> &>(me));
}

}  // namespace seqan

#endif  // #ifndef SEQAN_SEQUENCE_STRING_CSTYLE_H_
