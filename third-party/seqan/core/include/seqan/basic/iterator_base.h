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
// Author: Andreas Gogol-Döring <andreas.doering@mdc-berlin.de>
// ==========================================================================
// Basic declarations for the Iter class and generic implementations.
// ==========================================================================

// TODO(holtgrew): I think the interface is not completely specified here. Also, we could maybe have more generic implementations for operators?

#ifndef SEQAN_CORE_INCLUDE_SEQAN_BASIC_ITERATOR_BASE_H_
#define SEQAN_CORE_INCLUDE_SEQAN_BASIC_ITERATOR_BASE_H_

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

/*!
 * @class Iter
 * @extends IteratorAssociatedTypesConcept
 * @headerfile <seqan/basic.h>
 * @brief Base class for iterators to traverse containers.
 *
 * @signature template <typename TContainer, typename TSpec>
 *            class Iter;
 *
 * @tparam TContainer The type of the container to iterate.
 * @tparam TSpec      Type to use for specializing the <tt>Iter</tt> class.
 */

/**
.Class.Iter:
..cat:Basic
..summary:Iterator that is used to traverse containers.
..signature:Iter<TContainer, TSpec>
..param.TContainer:Type of the container that can be iterated by $Iter$.
...metafunction:Metafunction.Container
..param.TSpec:The specializing type.
...metafunction:Metafunction.Spec
..implements:Concept.IteratorAssociatedTypesConcept
..include:seqan/basic.h
*/

template <typename TContainer, typename TSpec>
class Iter;

// ============================================================================
// Metafunctions
// ============================================================================

// ----------------------------------------------------------------------------
// Metafunction IterComplementConst
// ----------------------------------------------------------------------------

/*!
 * @mfn Iter#IterComplementConst
 * @brief Metafunction that complements the const-ness of the container of an iterator.
 *
 * @signature IterComplementConst<TIter>::Type
 *
 * @tparam TIter The @link Iter @endlink to complement the container constness of.
 *
 * @return Type The type of the iterator that is the same as <tt>TIter</tt> except that the const-ness of the
 *              container is complemented.
 */

/**
.Metafunction.IterComplementConst:
..cat:Basic
..summary:Complement the "constness" of the container of an iterator.
..signature:IterComplementConst<TIterator>
..param.TIterator:The iterator to toggle the constness of its container for.
..include:seqan/basic.h
*/

template <typename TIterator>
struct IterComplementConst;

template <typename TContainer, typename TSpec>
struct IterComplementConst<Iter<TContainer, TSpec> >
{
    typedef Iter<
        typename IfC<
            IsSameType<typename RemoveConst_<TContainer>::Type, TContainer>::VALUE,
            TContainer const,
            typename RemoveConst_<TContainer>::Type>::Type,
        TSpec> Type;
};

template <typename TContainer, typename TSpec>
struct IterComplementConst<Iter<TContainer, TSpec> const>
        : IterComplementConst<Iter<TContainer, TSpec> > {};

// ----------------------------------------------------------------------------
// Metafunction IterMakeConst
// ----------------------------------------------------------------------------

/*!
 * @mfn Iter#IterMakeConst
 * @brief Metafunction to make enforce const-ness of the container of an iterator.
 *
 * @signature IterMakeConst<TIter>::Type
 *
 * @tparam TIter The iterator type to make the container const of.
 *
 * @return Type The resulting Iter type with a const container.
 */

/**
.Metafunction.IterMakeConst:
..cat:Basic
..summary:Make the container of an Iter const.
..signature:IterMakeConst<TIterator>
..param.TIterator:The iterator make the container const for.
..include:seqan/basic.h
*/

template <typename TIterator>
struct IterMakeConst;

template <typename TContainer, typename TSpec>
struct IterMakeConst<Iter<TContainer, TSpec> >
{
    typedef Iter<typename RemoveConst_<TContainer>::Type const, TSpec> Type;
};

template <typename TContainer, typename TSpec>
struct IterMakeConst<Iter<TContainer, TSpec> const>
        : IterMakeConst<Iter<TContainer, TSpec> > {};

// ----------------------------------------------------------------------------
// Metafunction Spec
// ----------------------------------------------------------------------------

/*!
 * @mfn Iter#Spec
 * @brief Return specialization tag of the <tt>Iter</tt> specialization.
 *
 * @signature Spec<TIter>::Type
 *
 * @tparam TIter The <tt>Iter</tt> class to get specialization tag of.
 *
 * @return Type The specialization tag used for the <tt>Iter</tt>.
 */

///.Metafunction.Spec.param.T.type:Class.Iter
///.Metafunction.Spec.class:Class.Iter

template <typename TContainer, typename TSpec>
struct Spec<Iter<TContainer, TSpec> >
{
    typedef TSpec Type;
};

template <typename TContainer, typename TSpec>
struct Spec<Iter<TContainer, TSpec> const>
{
    typedef TSpec Type;
};

// ----------------------------------------------------------------------------
// Metafunction Value
// ----------------------------------------------------------------------------

///.Metafunction.Value.param.T.type:Class.Iter
///.Metafunction.Value.class:Class.Iter

template <typename TContainer, typename TSpec>
struct Value<Iter<TContainer, TSpec> >:
    Value<TContainer>
{
};

template <typename TContainer, typename TSpec>
struct Value<Iter<TContainer, TSpec> const>:
    Value<TContainer>
{
};

// ----------------------------------------------------------------------------
// Metafunction GetValue
// ----------------------------------------------------------------------------

///.Metafunction.GetValue.param.T.type:Class.Iter
///.Metafunction.GetValue.class:Class.Iter

template <typename TContainer, typename TSpec>
struct GetValue<Iter<TContainer, TSpec> >:
    GetValue<TContainer>
{
};

template <typename TContainer, typename TSpec>
struct GetValue<Iter<TContainer, TSpec> const>:
    GetValue<TContainer>
{
};

// ----------------------------------------------------------------------------
// Metafunction Reference
// ----------------------------------------------------------------------------

///.Metafunction.Reference.param.T.type:Class.Iter
///.Metafunction.Reference.class:Class.Iter

template <typename TContainer, typename TSpec>
struct Reference<Iter<TContainer, TSpec> >:
    Reference<TContainer>
{
};

template <typename TContainer, typename TSpec>
struct Reference<Iter<TContainer, TSpec> const>:
    Reference<TContainer>
{
};

// ----------------------------------------------------------------------------
// Metafunction Container
// ----------------------------------------------------------------------------

/*!
 * @mfn Iter#Container
 * @brief The container type of the iterator.
 *
 * @signature Container<TIter>::Type
 *
 * @tparam TIter The <tt>TIter</tt> class to query for its container type.
 *
 * @return Type The container type of <tt>TIter</tt>
 */

///.Metafunction.Container.param.T.type:Class.Iter
///.Metafunction.Container.class:Class.Iter

template <typename T> struct Container;

template <typename TContainer, typename TSpec>
struct Container<Iter<TContainer, TSpec> >
{
    typedef TContainer Type;
};

template <typename TContainer, typename TSpec>
struct Container<Iter<TContainer, TSpec> const>
{
    typedef TContainer Type;
};

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function operator*()
// ----------------------------------------------------------------------------

template <typename TContainer, typename TSpec>
inline typename Reference<Iter<TContainer, TSpec> >::Type
operator*(Iter<TContainer, TSpec> & me)
{
    SEQAN_CHECKPOINT;
    return value(me);
}

template <typename TContainer, typename TSpec>
inline typename Reference<Iter<TContainer, TSpec> const>::Type
operator*(Iter<TContainer, TSpec> const & me)
{
    SEQAN_CHECKPOINT;
    return value(me);
}

// ----------------------------------------------------------------------------
// Function operator++()
// ----------------------------------------------------------------------------

template <typename TContainer, typename TSpec>
inline Iter<TContainer, TSpec> const &
operator++(Iter<TContainer, TSpec> & me)
{
    SEQAN_CHECKPOINT;
    goNext(me);
    return me;
}

template <typename TContainer, typename TSpec>
inline Iter<TContainer, TSpec>
operator++(Iter<TContainer, TSpec> & me, int)
{
    SEQAN_CHECKPOINT;
    Iter<TContainer, TSpec> temp_(me);
    goNext(me);
    return temp_;
}

// ----------------------------------------------------------------------------
// Function operator--()
// ----------------------------------------------------------------------------

template <typename TContainer, typename TSpec>
inline Iter<TContainer, TSpec> const &
operator--(Iter<TContainer, TSpec> & me)
{
    SEQAN_CHECKPOINT;
    goPrevious(me);
    return me;
}

template <typename TContainer, typename TSpec>
inline Iter<TContainer, TSpec>
operator--(Iter<TContainer, TSpec> & me, int)
{
    SEQAN_CHECKPOINT;
    Iter<TContainer, TSpec> temp_(me);
    goPrevious(me);
    return temp_;
}

// ----------------------------------------------------------------------------
// Function operator+()
// ----------------------------------------------------------------------------

template <typename TContainer, typename TSpec, typename TSize>
inline Iter<TContainer, TSpec>
operator+(Iter<TContainer, TSpec> const & me, TSize size)
{
    SEQAN_CHECKPOINT;
    Iter<TContainer, TSpec> temp_(me);
    goFurther(temp_, size);
    return temp_;
}

// ----------------------------------------------------------------------------
// Function operator+=()
// ----------------------------------------------------------------------------

template <typename TContainer, typename TSpec, typename TSize>
inline Iter<TContainer, TSpec> const &
operator+=(Iter<TContainer, TSpec> & me, TSize size)
{
    SEQAN_CHECKPOINT;
    goFurther(me, size);
    return me;
}

// ----------------------------------------------------------------------------
// Function operator-=()
// ----------------------------------------------------------------------------

/*
// TODO(doering): collides with Iter-Iter
// TODO(holtgrew): Try to reproduce error.
template <typename TContainer, typename TSpec, typename TSize>
inline Iter<TContainer, TSpec>
operator - (Iter<TContainer, TSpec> & me, TSize size)
{
SEQAN_CHECKPOINT
    Iter<TContainer, TSpec> temp_(me);
    goFurther(temp_, -size);
    return temp_;
}

template <typename TContainer, typename TSpec, typename TSize>
inline Iter<TContainer, TSpec>
operator - (Iter<TContainer, TSpec> const & me, TSize size)
{
SEQAN_CHECKPOINT
    Iter<TContainer, TSpec> temp_(me);
    goFurther(temp_, -size);
    return temp_;
}

template <typename TContainer, typename TSpec, typename TSize>
inline Iter<TContainer, TSpec> const &
operator -= (Iter<TContainer, TSpec> & me, TSize size)
{
SEQAN_CHECKPOINT
    goFurther(me, -size);
    return me;
}
*/

// ----------------------------------------------------------------------------
// Function position()
// ----------------------------------------------------------------------------

template <typename TContainer, typename TSpec, typename TContainer2>
inline typename Position<Iter<TContainer, TSpec> const>::Type
position(Iter<TContainer, TSpec> const & me,
         TContainer2 const &)
{
    SEQAN_CHECKPOINT;
    return position(me);
}

}  // namespace seqan

#endif  // #ifndef SEQAN_CORE_INCLUDE_SEQAN_BASIC_ITERATOR_BASE_H_
