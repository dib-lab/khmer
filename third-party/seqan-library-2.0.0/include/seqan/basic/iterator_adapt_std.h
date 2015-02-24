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
// Adaption of STL iterators to the SeqAn Iter class and vice versa through
// iterator traits.
// ==========================================================================

#ifndef SEQAN_INCLUDE_SEQAN_BASIC_ITERATOR_ADAPT_STD_H_
#define SEQAN_INCLUDE_SEQAN_BASIC_ITERATOR_ADAPT_STD_H_

// ============================================================================
// Adaption of SeqAn Iterators to STL Iterators.
// ============================================================================

namespace std
{
    template<typename TContainer, typename TSpec>
    struct iterator_traits<seqan::Iter<TContainer, TSpec> > // nolint
    {
        typedef seqan::Iter<TContainer, TSpec> TIter; // nolint

        typedef random_access_iterator_tag iterator_category; // nolint
        typedef typename seqan::Value<TIter>::Type value_type; // nolint
        typedef typename seqan::Difference<TIter>::Type difference_type; // nolint
        typedef typename seqan::Value<TIter>::Type * pointer; // nolint
        typedef typename seqan::Reference<TIter>::Type reference; // nolint
    };
}

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// Used in the STD Iterator Adapter specialization.
template <typename TStdContainer>
struct StdContainerIterator;

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// TODO(holtgrew): Rename ot STL Adaptor Iterator?

/*!
 * @class StdAdaptorIterator
 * @extends Iter
 * @headerfile <seqan/basic.h>
 *
 * @brief Adapt STL iterators to SeqAn iterators.
 *
 * @signature template <typename TContainer>
 *            class Iter<TContaier, StdIteratorAdaptor>;
 *
 * @tparam TContainer The container to iterate over.
 *
 * This type is a wrapper around STL iterators that has a conversion operator back to the underlying iterator.
 */

struct StdIteratorAdaptor_;
typedef Tag<StdIteratorAdaptor_> StdIteratorAdaptor;

template <typename TContainer>
class Iter<TContainer, StdIteratorAdaptor>
{
public:
    typedef typename StdContainerIterator<TContainer>::Type TIterator;
    TIterator data_iterator;

    // ------------------------------------------------------------------------
    // Constructors
    // ------------------------------------------------------------------------

    Iter() {}

    Iter(Iter const & other_)
            : data_iterator(other_.data_iterator)
    {}

    Iter(TIterator const & iter_)
            : data_iterator(iter_)
    {}

    template <typename TContainer_>
    Iter(Iter<TContainer_, StdIteratorAdaptor> const & other,
         SEQAN_CTOR_ENABLE_IF(IsSameType<TContainer, TContainer_ const>)) :
            data_iterator(other.data_iterator)
    {
        ignoreUnusedVariableWarning(dummy);
    }

    // TODO(holtgrew): Do we want this? Besides being non-essential, this creates a cross dependency on the begin function!
    Iter(TContainer const & cont_)
            : data_iterator(begin(cont_))
    {}

    // ------------------------------------------------------------------------
    // Assignment Operators;  Have to be defined in the class.
    // ------------------------------------------------------------------------

    Iter &
    operator = (Iter const & other_)
    {
        data_iterator = other_.data_iterator;
        return *this;
    }

    Iter &
    operator = (TIterator const & iter_)
    {
        data_iterator = iter_;
        return *this;
    }

    // ------------------------------------------------------------------------
    // Pointer Operators;  Have to be defined within class.
    // ------------------------------------------------------------------------

    typename Value<Iter>::Type *
    operator->()
    {
        return &*data_iterator;
    }

    typename Value<Iter>::Type const *
    operator->() const
    {
        return &*data_iterator;
    }

    // ------------------------------------------------------------------------
    // Conversion Operators;  Have to be defined within class.
    // ------------------------------------------------------------------------

    operator TIterator & ()
    {
        return data_iterator;
    }
};

// ============================================================================
// Metafunctions
// ============================================================================

// ----------------------------------------------------------------------------
// Metafunction Value
// ----------------------------------------------------------------------------

template <typename TContainer>
struct Value<Iter<TContainer, StdIteratorAdaptor> >
{
    typedef typename TContainer::iterator TIterator_;
    typedef typename std::iterator_traits<TIterator_>::value_type Type;
};

template <typename TContainer>
struct Value<Iter<TContainer const, StdIteratorAdaptor> >
{
    typedef TContainer const TContainer_;
    typedef typename TContainer::const_iterator TIterator_;
    typedef typename std::iterator_traits<TIterator_>::value_type Type;
};

// ----------------------------------------------------------------------------
// Metafunction GetValue
// ----------------------------------------------------------------------------

template <typename TContainer>
struct GetValue<Iter<TContainer, StdIteratorAdaptor> >
{
    typedef typename TContainer::const_iterator TIterator_;
    typedef typename std::iterator_traits<TIterator_>::reference Type;
};

template <typename TContainer>
struct GetValue<Iter<TContainer const, StdIteratorAdaptor> >
{
    typedef TContainer const TContainer_;
    typedef typename TContainer::const_iterator TIterator_;
    typedef typename std::iterator_traits<TIterator_>::reference Type;
};

// ----------------------------------------------------------------------------
// Metafunction Reference
// ----------------------------------------------------------------------------

template <typename TContainer>
struct Reference<Iter<TContainer, StdIteratorAdaptor> >
{
    typedef typename TContainer::iterator TIterator_;
    typedef typename std::iterator_traits<TIterator_>::reference Type;
};

template <typename TContainer>
struct Reference<Iter<TContainer, StdIteratorAdaptor> const> :
        Reference<Iter<TContainer, StdIteratorAdaptor> > {};

template <typename TContainer>
struct Reference<Iter<TContainer const, StdIteratorAdaptor> >
{
    typedef TContainer const TContainer_;
    typedef typename TContainer::const_iterator TIterator_;
    typedef typename std::iterator_traits<TIterator_>::reference Type;
};

template <typename TContainer>
struct Reference<Iter<TContainer const, StdIteratorAdaptor> const> :
        Reference<Iter<TContainer const, StdIteratorAdaptor> > {};

// ----------------------------------------------------------------------------
// Metafunction StdContainerIterator
// ----------------------------------------------------------------------------

// TODO(holtgrew): This is a candidate for not beging publically documented

template <typename TStdContainer>
struct StdContainerIterator;

template <typename TStdContainer>
struct StdContainerIterator
{
     typedef typename TStdContainer::iterator Type;
};

template <typename TStdContainer>
struct StdContainerIterator<TStdContainer const>
{
     typedef typename TStdContainer::const_iterator Type;
};

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function hostIterator()
// ----------------------------------------------------------------------------

template <typename TContainer>
inline typename StdContainerIterator<TContainer>::Type &
hostIterator(Iter<TContainer, StdIteratorAdaptor> & me)
{
    return me.data_iterator;
}

template <typename TContainer>
inline typename StdContainerIterator<TContainer>::Type const &
hostIterator(Iter<TContainer, StdIteratorAdaptor> const & me)
{
    return me.data_iterator;
}

// ----------------------------------------------------------------------------
// Function value()
// ----------------------------------------------------------------------------

template <typename TContainer>
inline typename Reference<Iter<TContainer, StdIteratorAdaptor> >::Type
value(Iter<TContainer, StdIteratorAdaptor> & me)
{
    return *(me.data_iterator);
}

template <typename TContainer>
inline typename Reference<Iter<TContainer, StdIteratorAdaptor> const>::Type
value(Iter<TContainer, StdIteratorAdaptor> const & me)
{
    return *(me.data_iterator);
}

// ----------------------------------------------------------------------------
// Function operator*()
// ----------------------------------------------------------------------------

template <typename TContainer>
inline typename Reference<Iter<TContainer, StdIteratorAdaptor> >::Type
operator*(Iter<TContainer, StdIteratorAdaptor> & me)
{
    return *(me.data_iterator);
}

template <typename TContainer>
inline typename Reference<Iter<TContainer, StdIteratorAdaptor> const>::Type
operator*(Iter<TContainer, StdIteratorAdaptor> const & me)
{
    return *(me.data_iterator);
}

// ----------------------------------------------------------------------------
// Function assignValue()
// ----------------------------------------------------------------------------

template <typename TContainer, typename TValue>
inline void
assignValue(Iter<TContainer, StdIteratorAdaptor> & me,
            TValue & val)
{
    *(me.data_iterator) = val;
}

template <typename TContainer, typename TValue>
inline void
assignValue(Iter<TContainer, StdIteratorAdaptor> & me,
            TValue const & val)
{
    *(me.data_iterator) = val;
}

// ----------------------------------------------------------------------------
// Function moveValue()
// ----------------------------------------------------------------------------

template <typename TContainer, typename TValue>
inline void
moveValue(Iter<TContainer, StdIteratorAdaptor> & me,
          TValue & val)
{
    move(*(me.data_iterator), val);
}
template <typename TContainer, typename TValue>
inline void
moveValue(Iter<TContainer, StdIteratorAdaptor> & me,
          TValue const & val)
{
    move(*(me.data_iterator), val);
}

// ----------------------------------------------------------------------------
// Function operator==()
// ----------------------------------------------------------------------------

template <typename TContainer>
inline bool
operator==(Iter<TContainer, StdIteratorAdaptor> const & left,
           Iter<TContainer, StdIteratorAdaptor> const & right)
{
    SEQAN_CHECKPOINT;
    return hostIterator(left) == hostIterator(right);
}

// ----------------------------------------------------------------------------
// Function operator!=()
// ----------------------------------------------------------------------------

template <typename TContainer>
inline bool
operator!=(Iter<TContainer, StdIteratorAdaptor> const & left,
           Iter<TContainer, StdIteratorAdaptor> const & right)
{
    SEQAN_CHECKPOINT;
    return hostIterator(left) != hostIterator(right);
}

// ----------------------------------------------------------------------------
// Function operator<()
// ----------------------------------------------------------------------------

template <typename TContainer>
inline bool
operator<(Iter<TContainer, StdIteratorAdaptor> const & left,
          Iter<TContainer, StdIteratorAdaptor> const & right)
{
    SEQAN_CHECKPOINT;
    return hostIterator(left) < hostIterator(right);
}

// ----------------------------------------------------------------------------
// Function operator>()
// ----------------------------------------------------------------------------

template <typename TContainer>
inline bool
operator>(Iter<TContainer, StdIteratorAdaptor> const & left,
          Iter<TContainer, StdIteratorAdaptor> const & right)
{
    SEQAN_CHECKPOINT;
    return hostIterator(left) > hostIterator(right);
}

// ----------------------------------------------------------------------------
// Function operator<=()
// ----------------------------------------------------------------------------

template <typename TContainer>
inline bool
operator<=(Iter<TContainer, StdIteratorAdaptor> const & left,
           Iter<TContainer, StdIteratorAdaptor> const & right)
{
    return hostIterator(left) <= hostIterator(right);
}

// ----------------------------------------------------------------------------
// Function operator>=()
// ----------------------------------------------------------------------------

template <typename TContainer>
inline bool
operator>=(Iter<TContainer, StdIteratorAdaptor> const & left,
           Iter<TContainer, StdIteratorAdaptor> const & right)
{
    SEQAN_CHECKPOINT;
    return hostIterator(left) >= hostIterator(right);
}

// ----------------------------------------------------------------------------
// Function goNext()
// ----------------------------------------------------------------------------

template <typename TContainer>
inline void
goNext(Iter<TContainer, StdIteratorAdaptor> & me)
{
    SEQAN_CHECKPOINT;
    goNext(hostIterator(me));
}

// ----------------------------------------------------------------------------
// Function goPrevious()
// ----------------------------------------------------------------------------

template <typename TContainer>
inline void
goPrevious(Iter<TContainer, StdIteratorAdaptor> & me)
{
    SEQAN_CHECKPOINT;
    goPrevious(hostIterator(me));
}

// ----------------------------------------------------------------------------
// Function operator+()
// ----------------------------------------------------------------------------

template <typename TContainer, typename TIntegral>
inline Iter<TContainer, StdIteratorAdaptor>
operator+(Iter<TContainer, StdIteratorAdaptor> const & left,
          TIntegral right)
{
    SEQAN_CHECKPOINT;
    return Iter<TContainer, StdIteratorAdaptor>(hostIterator(left) + right);
}

// for <anonymous enum> types
template <typename TContainer>
inline Iter<TContainer, StdIteratorAdaptor>
operator+(Iter<TContainer, StdIteratorAdaptor> const & left,
          int right)
{
    SEQAN_CHECKPOINT;
    return Iter<TContainer, StdIteratorAdaptor>(hostIterator(left) + right);
}

template <typename TContainer, typename TIntegral>
inline Iter<TContainer, StdIteratorAdaptor>
operator+(TIntegral left,
          Iter<TContainer, StdIteratorAdaptor> const & right)
{
    SEQAN_CHECKPOINT;
    return Iter<TContainer, StdIteratorAdaptor>(hostIterator(right) + left);
}

// for <anonymous enum> types
template <typename TContainer>
inline Iter<TContainer, StdIteratorAdaptor>
operator+(int left,
          Iter<TContainer, StdIteratorAdaptor> const & right)
{
    SEQAN_CHECKPOINT;
    return Iter<TContainer, StdIteratorAdaptor>(hostIterator(right) + left);
}

// ----------------------------------------------------------------------------
// Function operator+=()
// ----------------------------------------------------------------------------

template <typename TContainer, typename TIntegral>
inline Iter<TContainer, StdIteratorAdaptor> &
operator+=(Iter<TContainer, StdIteratorAdaptor> & left,
           TIntegral right)
{
    SEQAN_CHECKPOINT;
    hostIterator(left) += right;
    return left;
}

// for <anonymous enum> types
template <typename TContainer>
inline Iter<TContainer, StdIteratorAdaptor> &
operator+=(Iter<TContainer, StdIteratorAdaptor> & left,
           int right)
{
    SEQAN_CHECKPOINT;
    hostIterator(left) += right;
    return left;
}

// ----------------------------------------------------------------------------
// Function operator-()
// ----------------------------------------------------------------------------

template <typename TContainer, typename TIntegral>
inline Iter<TContainer, StdIteratorAdaptor>
operator-(Iter<TContainer, StdIteratorAdaptor> const & left,
          TIntegral right)
{
    SEQAN_CHECKPOINT;
    return Iter<TContainer, StdIteratorAdaptor>(hostIterator(left) - right);
}

// for <anonymous enum> types
template <typename TContainer>
inline Iter<TContainer, StdIteratorAdaptor>
operator-(Iter<TContainer, StdIteratorAdaptor> const & left,
          int right)
{
SEQAN_CHECKPOINT
    return Iter<TContainer, StdIteratorAdaptor>(hostIterator(left) - right);
}

template <typename TContainer>
inline typename Difference<Iter<TContainer, StdIteratorAdaptor> >::Type
operator-(Iter<TContainer, StdIteratorAdaptor> const & left,
          Iter<TContainer, StdIteratorAdaptor> const & right)
{
    SEQAN_CHECKPOINT;
    return hostIterator(left) - hostIterator(right);
}

// ----------------------------------------------------------------------------
// Function operator-=()
// ----------------------------------------------------------------------------

template <typename TContainer, typename TIntegral>
inline Iter<TContainer, StdIteratorAdaptor> &
operator-=(Iter<TContainer, StdIteratorAdaptor> & left,
           TIntegral right)
{
    SEQAN_CHECKPOINT;
    hostIterator(left) -= right;
    return left;
}

// for <anonymous enum> types
template <typename TContainer>
inline Iter<TContainer, StdIteratorAdaptor> &
operator -= (Iter<TContainer, StdIteratorAdaptor> & left,
             int right)
{
    SEQAN_CHECKPOINT;
    hostIterator(left) -= right;
    return left;
}

// ----------------------------------------------------------------------------
// Function assign()
// ----------------------------------------------------------------------------

// Conversion assignment.
template <typename TTargetContainer, typename TSource>
inline void
assign(Iter<TTargetContainer, StdIteratorAdaptor> & target,
       TSource const & source)
{
    SEQAN_CHECKPOINT;
    target.data_iterator = begin(container(source)) + position(source);
}

}  // namespace seqan

#endif  // #ifndef SEQAN_INCLUDE_SEQAN_BASIC_ITERATOR_ADAPT_STD_H_
