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

#ifndef SEQAN_HEADER_MODIFIER_ALPHABET_EXPANSION_H
#define SEQAN_HEADER_MODIFIER_ALPHABET_EXPANSION_H

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////

// TODO(holtgrew): Rename ot ExpandedAlphabet?

/*!
 * @class AlphabetExpansion
 * @headerfile <seqan/modifier.h>
 * @brief Modifier that adds a character to an alphabet.
 *
 * @signature template <typename TAlphabet, char CHAR[, typename TSpec]>
 *            class ModifiedAlphabet<TAlphabet, ModExpand<CHAR,TSpec> >;
 *
 * @tparam TAlphabet Original value type.
 * @tparam CHAR      <tt>char</tt> character that specifies, what value should added to the alphabet.  <tt>CHAR</tt>
 *                   should not be a <tt>char</tt> that already stands for a value in <tt>TAlphabet</tt>.  For
 *                   example, do not use <tt>'A'</tt> or <tt>'a'</tt> as <tt>CHAR</tt> when expanding Dna.
 * @tparam TSpec     Optional specialization tag.  This modifier is intended to expand SimpleType classes, default is <tt>Default</tt>.
 *
 * @section Special Characters
 *
 * Some values of <tt>CHAR</tt> have special meaning:
 *
 * <dl>
 *   <dt><tt>'-'</tt></dt>
 *   <dd>A gap character.  The value in the expanded alphabet that corresponds to <tt>'-'</tt> will be returned
 *       by the function gapValue.</dd>
 *   <dt><tt>'$'</tt></dt>
 *   <dd>An end of string character.</dd>
 * </dl>
 */

template <char CHAR, typename TSpec = Default>
struct ModExpand;


template <typename THost, char CHAR, typename TSpec>
class ModifiedAlphabet<THost, ModExpand<CHAR, TSpec> >
{
public:
    typedef typename IntegralForValue<ModifiedAlphabet>::Type TData;
    TData data;

    ModifiedAlphabet() : data()
    {}

    ModifiedAlphabet(ModifiedAlphabet const & other)
        : data(_internalOrdValue(other))
    {}

    template <typename TOther>
    ModifiedAlphabet(TOther const & other_data)
        : data(_internalOrdValue(convert<ModifiedAlphabet>(other_data)))
    {}

    ModifiedAlphabet const &
    operator = (ModifiedAlphabet const & other)
    {
        data = other.data;
        return *this;
    }

    template <typename TOther>
    ModifiedAlphabet const &
    operator = (TOther const & other_data)
    {
        data = _internalOrdValue(convert<ModifiedAlphabet>(other_data));
        return *this;
    }

/*    operator TData ()
    {
        return data;
    }

    operator THost()
    {
        return convert<THost>(data);
    }
*/
//____________________________________________________________________________

    //TODO(weese): investigate the issue below
    //this cannot be a template since a template would be in conflict to
    //the template c'tor

/*
    // weese: I tried the template below without problems in the tests
    template <typename TValue>
    operator TValue() const
    {
SEQAN_CHECKPOINT
        return convert<TValue>(*this);
    }
*/
    operator long() const
    {
SEQAN_CHECKPOINT
        return convert<long>(*this);
    }
    operator unsigned long() const
    {
SEQAN_CHECKPOINT
        return convert<unsigned long>(*this);
    }
    operator int() const
    {
SEQAN_CHECKPOINT
        return convert<int>(*this);
    }
    operator unsigned int() const
    {
SEQAN_CHECKPOINT
        return convert<unsigned int>(*this);
    }
    operator short() const
    {
SEQAN_CHECKPOINT
        return convert<short>(*this);
    }
    operator unsigned short() const
    {
SEQAN_CHECKPOINT
        return convert<unsigned short>(*this);
    }
    operator char() const
    {
SEQAN_CHECKPOINT
        return convert<char>(*this);
    }
    operator signed char() const
    {
SEQAN_CHECKPOINT
        return convert<signed char>(*this);
    }
    operator unsigned char() const
    {
SEQAN_CHECKPOINT
        return convert<unsigned char>(*this);
    }
};


//////////////////////////////////////////////////////////////////////////////
// unknownValueImpl()
//////////////////////////////////////////////////////////////////////////////

template <typename THost, char CHAR, typename TSpec>
inline ModifiedAlphabet<THost, ModExpand<CHAR, TSpec> >
unknownValueImpl(ModifiedAlphabet<THost, ModExpand<CHAR, TSpec> > *)
{
    static const ModifiedAlphabet<THost, ModExpand<CHAR, TSpec> > n = 'N';
    return n;
}

//////////////////////////////////////////////////////////////////////////////
// sizes
//////////////////////////////////////////////////////////////////////////////

template <typename THost, char CHAR, typename TSpec>
struct BitsPerValue<ModifiedAlphabet<THost, ModExpand<CHAR, TSpec> > >
{
    typedef ModifiedAlphabet<THost, ModExpand<CHAR, TSpec> > TValue;
    enum { VALUE = Log2< ValueSize<TValue>::VALUE >::VALUE };
};

//////////////////////////////////////////////////////////////////////////////

template <typename THost, char CHAR, typename TSpec>
struct ValueSize<ModifiedAlphabet<THost, ModExpand<CHAR, TSpec> > >
{
    enum { VALUE = ValueSize<THost>::VALUE + 1 };
};

template <typename THost, char CHAR, typename TSpec>
struct InternalValueSize_<ModifiedAlphabet<THost, ModExpand<CHAR, TSpec> > >
{
    enum { VALUE = InternalValueSize_<THost>::VALUE + 1 };
};

//////////////////////////////////////////////////////////////////////////////
// _internalCreateChar
//////////////////////////////////////////////////////////////////////////////

template <typename TValue, typename TNum>
inline TValue
_internalCreateChar(TValue const &, TNum i)
{
    return convert<TValue>(i);
}

template <typename TValue, typename TSpec, typename TNum>
inline SimpleType<TValue, TSpec>
_internalCreateChar(SimpleType<TValue, TSpec> const &, TNum i)
{
    SimpleType<TValue, TSpec> s;
    s.value = i;
    return s;
}

//////////////////////////////////////////////////////////////////////////////
// gapValueImpl()
//////////////////////////////////////////////////////////////////////////////

template <typename THost, char CHAR, typename TSpec>
inline ModifiedAlphabet<THost, ModExpand<CHAR, TSpec> >
gapValueImpl(ModifiedAlphabet<THost, ModExpand<CHAR, TSpec> > /*const*/ *)
{
    typedef ModifiedAlphabet<THost, ModExpand<CHAR, TSpec> > T;
    if (CHAR == '-')
    {
        return _internalCreateChar(T(), (unsigned)ValueSize<T>::VALUE - 1);
    }
    else
    {
        THost * ptr = 0;
        return gapValueImpl(ptr);
    }
}


//////////////////////////////////////////////////////////////////////////////
// conversions
//////////////////////////////////////////////////////////////////////////////


// some type => ModExpand
template <typename THost, char CHAR, typename TSpec, typename TSource>
inline void
_initializeAlphabetConversionTable(ModifiedAlphabet<THost, ModExpand<CHAR, TSpec> > * buf,
                                   TSource const &)
{
    //assure that the conversion from TSource to THost is possible
//    AlphabetConversionTable_<THost, TSource>::initialize();

    //copy the conversion table for converting TSouce => THost
    //maybe, if there is no CHAR in TSource, the entry for CHAR is overwritten now
    for (int i = InternalValueSize_<TSource>::VALUE; i > 0; )
    {
        --i;
        buf[i].data = _internalOrdValue(convert<THost>(_internalCreateChar(TSource(), i)));
    }

    //add the new character CHAR to the table
    buf[_internalOrdValue(convert<TSource>(CHAR))].data = InternalValueSize_<THost>::VALUE;
}

template <int SIZE_OF_SOURCE>
struct ConvertImplModExpand_
{
    //default implementation for large source types
    template <typename THost, char CHAR, typename TSpec, typename T, typename TSource>
    inline
    static typename Convert<ModifiedAlphabet<THost, ModExpand<CHAR, TSpec> >, TSource>::Type
    _convertImpl(Convert<ModifiedAlphabet<THost, ModExpand<CHAR, TSpec> >, T> const,
        TSource const & source_)
    {
        typedef ModifiedAlphabet<THost, ModExpand<CHAR, TSpec> > TTarget;
        if (source_ == ValueSize<THost>::VALUE)
        {// the extra character
            TTarget tmp;
            tmp.data = InternalValueSize_<THost>::VALUE;
            return tmp;
        }
        return convert<TTarget>(convert<THost>(source_));
    }
};

//for 1 byte source: use translation table
template <>
struct ConvertImplModExpand_<1>
{
    template <typename THost, char CHAR, typename TSpec, typename T, typename TSource>
    inline
    static typename Convert<ModifiedAlphabet<THost, ModExpand<CHAR, TSpec> >, TSource>::Type
    _convertImpl(Convert<ModifiedAlphabet<THost, ModExpand<CHAR, TSpec> >, T> const,
        TSource const & source_)
    {
        typedef ModifiedAlphabet<THost, ModExpand<CHAR, TSpec> > TTarget;
        TTarget * table = AlphabetConversionTable_<TTarget, TSource>::table;
        return table[_internalOrdValue(source_)];
    }
};

//generic source: dispatch for size of BytesPerValue
template <typename THost, char CHAR, typename TSpec, typename T, typename TSource>
inline typename Convert<ModifiedAlphabet<THost, ModExpand<CHAR, TSpec> >, TSource>::Type
convertImpl(Convert<ModifiedAlphabet<THost, ModExpand<CHAR, TSpec> >, T> const convert_,
            TSource const & source_)
{
    return ConvertImplModExpand_<BytesPerValue<TSource>::VALUE>::_convertImpl(convert_, source_);
}

//for SimpleType sources
template <typename THost, char CHAR, typename TSpec, typename T, typename TSourceValue, typename TSourceSpec>
inline typename Convert<ModifiedAlphabet<THost, ModExpand<CHAR, TSpec> >, SimpleType<TSourceValue, TSourceSpec> >::Type
convertImpl(Convert<ModifiedAlphabet<THost, ModExpand<CHAR, TSpec> >, T> const,
            SimpleType<TSourceValue, TSourceSpec> const & source_)
{
SEQAN_CHECKPOINT
    typedef ModifiedAlphabet<THost, ModExpand<CHAR, TSpec> > TTarget;
    typedef SimpleType<TSourceValue, TSourceSpec> TSource;
    return AlphabetConversionTable_<TTarget, TSource>::table[_internalOrdValue(source_)];
}


// For SimpleType sources with the same underlying alphabet.
template <char CHAR, typename TSpec, typename T, typename TSourceValue, typename TSourceSpec>
inline typename Convert<ModifiedAlphabet<SimpleType<TSourceValue, TSourceSpec>, ModExpand<CHAR, TSpec> >, SimpleType<TSourceValue, TSourceSpec> >::Type
convertImpl(Convert<ModifiedAlphabet<SimpleType<TSourceValue, TSourceSpec>, ModExpand<CHAR, TSpec> >, T> const,
            SimpleType<TSourceValue, TSourceSpec> const & source_)
{
    SEQAN_CHECKPOINT;
    typedef SimpleType<TSourceValue, TSourceSpec> TSource;
    typedef ModifiedAlphabet<TSource, ModExpand<CHAR, TSpec> > TTarget;

    TTarget tmp;
    tmp.data = _internalOrdValue(source_);
    return tmp;
}


//for Proxy sources
template <typename THost, char CHAR, typename TSpec, typename T, typename TSpec2>
inline typename Convert<ModifiedAlphabet<THost, ModExpand<CHAR, TSpec> >, Proxy<TSpec2> >::Type
convertImpl(Convert<ModifiedAlphabet<THost, ModExpand<CHAR, TSpec> >, T> const,
            Proxy<TSpec2> const & source_)
{
SEQAN_CHECKPOINT
    typedef ModifiedAlphabet<THost, ModExpand<CHAR, TSpec> > TTarget;
    return convert<TTarget>(getValue(source_));
}


// ModExpand => some type

template <typename TTarget, typename THost, char CHAR, typename TSpec>
inline void
_initializeAlphabetConversionTable(TTarget * buf,
                                   ModifiedAlphabet<THost, ModExpand<CHAR, TSpec> > const &)
{
    //assure that the conversion from THost to TTarget is possible
    AlphabetConversionTable_<TTarget, THost>::initialize();

    //copy the conversion table for converting THost => TTarget
    for (int i = InternalValueSize_<THost>::VALUE; i > 0; )
    {
        --i;
        buf[i] = convert<TTarget>(_internalCreateChar(THost(), i));
    }

    //add the new character CHAR to the table
    buf[InternalValueSize_<THost>::VALUE] = convert<TTarget, char>(CHAR);
}

template <typename TTarget, typename THost, char CHAR, typename TSpec>
inline void
_initializeAlphabetOrdTable(TTarget * buf,
                            ModifiedAlphabet<THost, ModExpand<CHAR, TSpec> > const &)
{
    //assure that the conversion from THost to TTarget is possible
    AlphabetOrdTable_<THost>::initialize();

    //copy the conversion table for converting THost => TTarget
    for (int i = InternalValueSize_<THost>::VALUE; i > 0; )
    {
        --i;
        buf[i] = ordValue(_internalCreateChar(THost(), i));
    }

    //add the new character CHAR to the table
    buf[InternalValueSize_<THost>::VALUE] = ValueSize<THost>::VALUE;
}


// Conversion from modified alphabet to non-modified alphabet with an
// arbitrary type.  The conversion is done through an alphabet
// conversion table.  See below for a specialization where the
// underlying type of the modified alphabet and the target type are
// the same SimpleType specialization.
template <typename TTarget, typename T, typename THost, char CHAR, typename TSpec>
inline typename Convert<TTarget, ModifiedAlphabet<THost, ModExpand<CHAR, TSpec> > >::Type
convertImpl(Convert<TTarget, T> const,
            ModifiedAlphabet<THost, ModExpand<CHAR, TSpec> > const & source_)
{
    SEQAN_CHECKPOINT;
    typedef ModifiedAlphabet<THost, ModExpand<CHAR, TSpec> > TSource;
    return AlphabetConversionTable_<TTarget, TSource>::table[_internalOrdValue(source_)];
}


// Conversion from modified alphabet to non-modified with the same underlying type.
template <typename TTargetValue, typename TTargetSpec, char CHAR, typename TSpec>
inline typename Convert<SimpleType<TTargetValue, TTargetSpec>, ModifiedAlphabet<SimpleType<TTargetValue, TTargetSpec>, ModExpand<CHAR, TSpec> > >::Type
convertImpl(Convert<SimpleType<TTargetValue, TTargetSpec>, ModifiedAlphabet<SimpleType<TTargetValue, TTargetSpec>, ModExpand<CHAR, TSpec> > > const,
            ModifiedAlphabet<SimpleType<TTargetValue, TTargetSpec>, ModExpand<CHAR, TSpec> > const & source_)
{
    SEQAN_CHECKPOINT;
    typedef SimpleType<TTargetValue, TTargetSpec> TTarget;
    TTarget target;
    if (source_.data == InternalValueSize_<TTarget>::VALUE)
        assign(target, CHAR);
    else
        target.value = source_.data;
    return target;
}


//////////////////////////////////////////////////////////////////////////////


/*
template
<
    typename TTargetHost, char TARGET_CHAR, typename TTargetSpec, typename T,
    typename TSourceHost, char SOURCE_CHAR, typename TSourceSpec
>
inline typename Convert<ModifiedAlphabet<TTargetHost, ModExpand<TARGET_CHAR, TTargetSpec> > , ModifiedAlphabet<TSourceHost, ModExpand<SOURCE_CHAR, TSourceSpec> > >::Type
convertImpl(Convert<ModifiedAlphabet<TTargetHost, ModExpand<TARGET_CHAR, TTargetSpec> >, T> const,
            ModifiedAlphabet<TSourceHost, ModExpand<SOURCE_CHAR, TSourceSpec> > const & source_)
{
    ModifiedAlphabet<TTargetHost, ModExpand<TARGET_CHAR, TTargetSpec> > TTarget;
    ModifiedAlphabet<TSourceHost, ModExpand<SOURCE_CHAR, TSourceSpec> > TSource;
    return convert<TTarget>(convert<TTargetHost>(convert<TSourceHost>(source_)));
}
*/

//no conversion
template <typename THost, char CHAR, typename TSpec, typename T>
inline ModifiedAlphabet<THost, ModExpand<CHAR, TSpec> >
convertImpl(Convert<ModifiedAlphabet<THost, ModExpand<CHAR, TSpec> >, T> const,
            ModifiedAlphabet<THost, ModExpand<CHAR, TSpec> > const & source_)
{
    return source_;
}


//////////////////////////////////////////////////////////////////////////////

template <typename THost, char CHAR, typename TSpec>
inline unsigned
_internalOrdValue(ModifiedAlphabet<THost, ModExpand<CHAR, TSpec> > const & c)
{
    return c.data;
}

template <typename THost, char CHAR, typename TSpec>
inline unsigned
ordValue(ModifiedAlphabet<THost, ModExpand<CHAR, TSpec> > const & c)
{
    SEQAN_CHECKPOINT;
    typedef ModifiedAlphabet<THost, ModExpand<CHAR, TSpec> > TSource;
    return AlphabetOrdTable_<TSource>::table[_internalOrdValue(c)];
}

//////////////////////////////////////////////////////////////////////////////
// comparisons
//////////////////////////////////////////////////////////////////////////////

template <typename TModExpand, typename THost, typename TRight, typename TCompareHostRight>
struct CompareTypeModExpandImpl_
{
    typedef TCompareHostRight Type; //fallback
};
template <typename TModExpand, typename THost, typename TRight>
struct CompareTypeModExpandImpl_<TModExpand, THost, TRight, THost>
{
    typedef TModExpand Type;
};
template <typename TModExpand, typename THost, typename TRight>
struct CompareTypeModExpandImpl_<TModExpand, THost, TRight, TRight>
{
    typedef TRight Type;
};


template <typename THost, char CHAR, typename TSpec, typename TRight>
struct CompareType<ModifiedAlphabet<THost, ModExpand<CHAR, TSpec> >, TRight>
{
    typedef ModifiedAlphabet<THost, ModExpand<CHAR, TSpec> > TModExpand;
    typedef typename CompareType<THost, TRight>::Type TCompareHostRight;
    typedef typename CompareTypeModExpandImpl_<TModExpand, THost, TRight, TCompareHostRight>::Type Type;
};

//////////////////////////////////////////////////////////////////////////////

}// namespace

#endif
