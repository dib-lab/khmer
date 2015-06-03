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
// Author: Andreas Gogol-Döring <andreas.doering@mdc-berlin.de>
// ==========================================================================
// Iterator Proxy specialization through an iterator.
//
// The proxy stores an iterator to the referenced value.  Assignment is done
// using assignValue() on the iterator, reading via getValue().
// ==========================================================================

#ifndef SEQAN_INCLUDE_SEQAN_BASIC_PROXY_ITERATOR_H_
#define SEQAN_INCLUDE_SEQAN_BASIC_PROXY_ITERATOR_H_

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

/*!
 * @class IteratorProxy
 * @extends Proxy
 * @headerfile <seqan/basic.h>
 * @brief Proxy that is implemented by an iterator.
 *
 * @signature template <typename TIterator>
 *            class Proxy<IteratorProxy<TIterator> >;
 *
 * @tparam TIterator The iterator to store internally.
 */

// TODO(holtgrew): Assignment works through operator=() but we also need set() and move()!
// TODO(holtgrew): Proxy cannot work correctly for const containers, but should! RemoveConst_ removes the const of the values inside pointers and references.

template <typename TIterator>
struct IteratorProxy;

template <typename TIterator>
class Proxy<IteratorProxy<TIterator> >
{
public:
    typedef typename Value<Proxy>::Type TValue_;
    typedef typename GetValue<Proxy>::Type TAccessor_;

    // TODO(holtgrew): Is removing this reference necessary or does this only hide errors in the definition of GetValue<>::Type?
    //typedef typename RemoveReference_<typename RemoveConst_<TAccessor_>::Type>::Type TAccessorNotConst_;
    typedef typename RemoveConst_<TAccessor_>::Type TAccessorNotConst_;

    TIterator data_iterator;

    // ------------------------------------------------------------------------
    // Constructors
    // ------------------------------------------------------------------------

    Proxy(TIterator const _it)
            : data_iterator(_it)
    {
        SEQAN_CHECKPOINT;
    }

    Proxy(Proxy const & _other)
            : data_iterator(_other.data_iterator)
    {
        SEQAN_CHECKPOINT;
    }

    // ------------------------------------------------------------------------
    // Assignment operators;  Have to be defined in class.
    // ------------------------------------------------------------------------

    Proxy const &
    operator=(Proxy const & _other) const
    {
        SEQAN_CHECKPOINT;
        assignValue(data_iterator, getValue(_other.data_iterator));
        return *this;
    }

    Proxy const &
    operator=(TValue_ const & _value) const
    {
        SEQAN_CHECKPOINT;
        assignValue(data_iterator, _value);
        return *this;
    }

    // ------------------------------------------------------------------------
    // Type conversion operators;  Have to be defined in class.
    // ------------------------------------------------------------------------

    operator TAccessorNotConst_() const
    {
        SEQAN_CHECKPOINT;
        return getValue(data_iterator);
    }

//    operator TAccessor_() const
//    {
//        return getValue(data_iterator);
//    }

//    operator unsigned int() const
//    {
//        return static_cast<unsigned int>(getValue(data_iterator));
//    }
};

// ============================================================================
// Metafunctions
// ============================================================================

// ----------------------------------------------------------------------------
// Metafunction Value
// ----------------------------------------------------------------------------

template <typename TIterator>
struct Value<Proxy<IteratorProxy<TIterator> > >
        : Value<TIterator>
{
};

template <typename TIterator>
struct Value<Proxy<IteratorProxy<TIterator> > const>
{
    typedef typename Value<TIterator>::Type const Type;
};

// ----------------------------------------------------------------------------
// Metafunction GetValue
// ----------------------------------------------------------------------------

template <typename TIterator>
struct GetValue<Proxy<IteratorProxy<TIterator> > >
        : GetValue<TIterator>
{
};

template <typename TIterator>
struct GetValue<Proxy<IteratorProxy<TIterator> > const>
{
    typedef typename GetValue<TIterator const>::Type Type;
};

// ----------------------------------------------------------------------------
// Metafunction Reference
// ----------------------------------------------------------------------------

template <typename TIterator>
struct Reference<Proxy<IteratorProxy<TIterator> > >
{
    typedef Proxy<IteratorProxy<TIterator> > Type;
};

template <typename TIterator>
struct Reference<Proxy<IteratorProxy<TIterator> > const>
{
    typedef Proxy<IteratorProxy<TIterator> > const Type;
};

// ----------------------------------------------------------------------------
// Metafunction Size
// ----------------------------------------------------------------------------

template <typename TIterator>
struct Size<Proxy<IteratorProxy<TIterator> > >
        : Size<TIterator>
{
};

template <typename TIterator>
struct Size<Proxy<IteratorProxy<TIterator> > const>
        : Size<TIterator>
{
};

// ----------------------------------------------------------------------------
// Metafunction Difference
// ----------------------------------------------------------------------------

template <typename TIterator>
struct Difference<Proxy<IteratorProxy<TIterator> > >
        : Difference<TIterator>
{
};

template <typename TIterator>
struct Difference<Proxy<IteratorProxy<TIterator> > const>
        : Difference<TIterator>
{
};

// ----------------------------------------------------------------------------
// Metafunction CompareType
// ----------------------------------------------------------------------------

template <typename TIterator, typename T>
struct CompareTypeRemoveProxy<Proxy<IteratorProxy<TIterator> >, T>:
    CompareTypeRemoveProxy<typename Value<TIterator>::Type, T> {};

template <typename TIterator, typename T>
struct CompareTypeRemoveProxy<T, Proxy<IteratorProxy<TIterator> > >:
    CompareTypeRemoveProxy<T, typename Value<TIterator>::Type> {};

template <typename TIterator1, typename TIterator2>
struct CompareTypeRemoveProxy<Proxy<IteratorProxy<TIterator1> >, Proxy<IteratorProxy<TIterator2> > >:
    CompareTypeRemoveProxy<typename Value<TIterator1>::Type, typename Value<TIterator2>::Type> {};

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function iter()
// ----------------------------------------------------------------------------

template <typename TIterator>
inline TIterator &
iter(Proxy<IteratorProxy<TIterator> > & me)
{
    return me.data_iterator;
}

template <typename TIterator>
inline TIterator const &
iter(Proxy<IteratorProxy<TIterator> > const & me)
{
    return me.data_iterator;
}

// ----------------------------------------------------------------------------
// Function getValue()
// ----------------------------------------------------------------------------

// TODO(holtgrew): Deletion candidate, why should a proxy provide the iterator interface?

template <typename TIterator>
inline typename GetValue<Proxy<IteratorProxy<TIterator> > >::Type
getValue(Proxy<IteratorProxy<TIterator> > & me)
{
    return getValue(iter(me));
}

template <typename TIterator>
inline typename GetValue<Proxy<IteratorProxy<TIterator> > const>::Type
getValue(Proxy<IteratorProxy<TIterator> > const & me)
{
    return getValue(iter(me));
}

// ------------------------------------------------------------------------
// Function assign()
// ------------------------------------------------------------------------

// additional const definitions (that are not covered by default assign() definition)
template <typename TIterator, typename TValue>
inline void
assign(Proxy<IteratorProxy<TIterator> > const & me, TValue & value)
{
    me = value;
}

template <typename TIterator, typename TValue>
inline void
assign(Proxy<IteratorProxy<TIterator> > const & me, TValue const & value)
{
    me = value;
}

// ----------------------------------------------------------------------------
// Function operator<<()
// ----------------------------------------------------------------------------

template <typename TStream, typename TIterator>
inline TStream &
operator<<(TStream & stream, Proxy<IteratorProxy<TIterator> > const & it)
{
    stream << getValue(it);
    return stream;
}

}  // namespace seqan



namespace std {

// ------------------------------------------------------------------------
// Function swap()
// ------------------------------------------------------------------------

// std::swap doesn't cover the case of const-ref proxies, so we provide it here
template <typename TIterator>
inline void
swap(seqan::Proxy<seqan::IteratorProxy<TIterator> > const & left, seqan::Proxy<seqan::IteratorProxy<TIterator> > const & right)
{
    typedef seqan::Proxy<seqan::IteratorProxy<TIterator> >  TProxy;
    typedef typename seqan::GetValue<TProxy>::Type          TGetValue;

    TGetValue tmp = getValue(iter(left));
    assignValue(iter(left), getValue(iter(right)));
    assignValue(iter(right), tmp);
}

}  // namespace std

#endif  // #ifndef SEQAN_INCLUDE_SEQAN_BASIC_PROXY_ITERATOR_H_

