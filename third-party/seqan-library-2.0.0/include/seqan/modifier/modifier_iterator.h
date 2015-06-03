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

#ifndef SEQAN_HEADER_MODIFIER_ITERATOR_H
#define SEQAN_HEADER_MODIFIER_ITERATOR_H

namespace SEQAN_NAMESPACE_MAIN
{

// ==========================================================================
// Forwards
// ==========================================================================

template <typename THost, typename TSpec> class ModifiedString;

// ==========================================================================
// Classes
// ==========================================================================

// --------------------------------------------------------------------------
// Class ModifiedIterator
// --------------------------------------------------------------------------

/*!
 * @class ModifiedIterator
 * @implements RandomAccessIteratorConcept
 * @headerfile <seqan/modifier.h>
 * @brief Allows you to modify arbitrary iterators by specializing what differs from an origin.
 *
 * @signature template <typename THost[, typename TSpec]>
 *            class ModifiedIterator;
 *
 * @tparam THost The host iterator type.
 * @tparam TSpec Tag used for the specialization, defaults to <tt>void</tt>.
 *
 * <tt>THost</tt> can also be a modified iterator, so you can create custom iterators by combining predefined ones.
 */

template <typename THost, typename TSpec = void>
class ModifiedIterator
{
public:
    typedef typename Cargo<ModifiedIterator>::Type TCargo_;

    THost _host;
    TCargo_ _cargo;

    ModifiedIterator()
    {}

    template <typename TOtherHost>
    ModifiedIterator(ModifiedIterator<TOtherHost, TSpec> const & origin):
            _host(origin._host), _cargo(origin._cargo)
    {}

    explicit
    ModifiedIterator(THost const & host): _host(host)
    {}
};

// ==========================================================================
// Metafunctions
// ==========================================================================

// --------------------------------------------------------------------------
// Metafunction Spec
// --------------------------------------------------------------------------

template < typename THost, typename TSpec >
struct Spec< ModifiedIterator<THost, TSpec> > {
    typedef TSpec Type;
};

template < typename THost, typename TSpec >
struct Spec< ModifiedIterator<THost, TSpec> const > {
    typedef TSpec Type;
};

// --------------------------------------------------------------------------
// Metafunction Value
// --------------------------------------------------------------------------

// an iterator is not the owner of the values pointing at
// it can be constant while
// - pointing to an alterable object
// - returning an non-constant value
// - being an iterator of an alterable container

template <typename THost, typename TSpec>
struct Value<ModifiedIterator<THost, TSpec> > : Value<THost>
{};

template <typename THost, typename TSpec>
struct Value<ModifiedIterator<THost, TSpec> const> : Value<ModifiedIterator<THost, TSpec> >
{};

// --------------------------------------------------------------------------
// Metafunction GetValue
// --------------------------------------------------------------------------

template <typename THost, typename TSpec>
struct GetValue< ModifiedIterator<THost, TSpec> > : GetValue<THost>
{};

template <typename THost, typename TSpec>
struct GetValue<ModifiedIterator<THost, TSpec> const> : GetValue<ModifiedIterator<THost, TSpec> >
{};

// --------------------------------------------------------------------------
// Metafunction Reference
// --------------------------------------------------------------------------

template <typename THost, typename TSpec>
struct Reference<ModifiedIterator<THost, TSpec> > : Reference<THost>
{};

template <typename THost, typename TSpec>
struct Reference<ModifiedIterator<THost, TSpec> const> : Reference< ModifiedIterator<THost, TSpec> >
{};

// --------------------------------------------------------------------------
// Metafunction Size
// --------------------------------------------------------------------------

template <typename THost, typename TSpec>
struct Size<ModifiedIterator<THost, TSpec> > : Size<THost>
{};

// --------------------------------------------------------------------------
// Metafunction Positions
// --------------------------------------------------------------------------

template <typename THost, typename TSpec>
struct Position<ModifiedIterator<THost, TSpec> > : Position<THost>
{};

// --------------------------------------------------------------------------
// Metafunction Difference
// --------------------------------------------------------------------------

template <typename THost, typename TSpec>
struct Difference<ModifiedIterator<THost, TSpec> > : Difference<THost>
{};

// --------------------------------------------------------------------------
// Metafunction Host
// --------------------------------------------------------------------------

template <typename THost, typename TSpec>
struct Host<ModifiedIterator<THost, TSpec> >
{
    typedef THost Type;
};

template <typename THost, typename TSpec>
struct Host<ModifiedIterator<THost, TSpec> const>
{
    typedef THost const Type;
};

// --------------------------------------------------------------------------
// Metafunction Container
// --------------------------------------------------------------------------

template <typename THost, typename TSpec >
struct Container<ModifiedIterator<THost, TSpec> >
{
    typedef typename Container<THost>::Type THostContainer;
    typedef ModifiedString<THostContainer, TSpec> Type;
};

template <typename THost, typename TSpec >
struct Container<ModifiedIterator<THost, TSpec> const>
{
    typedef typename Container<THost>::Type THostContainer;
    typedef ModifiedString<THostContainer, TSpec> Type;
};

// ==========================================================================
// Functions
// ==========================================================================

// --------------------------------------------------------------------------
// Function host()                                         [ModifiedIterator]
// --------------------------------------------------------------------------

template <typename THost, typename TSpec>
inline typename Host<ModifiedIterator<THost, TSpec> >::Type &
host(ModifiedIterator<THost, TSpec> & me)
{
    return me._host;
}

template <typename THost, typename TSpec>
inline typename Host<ModifiedIterator<THost, TSpec> const>::Type &
host(ModifiedIterator<THost, TSpec> const & me)
{
    return me._host;
}

// // --------------------------------------------------------------------------
// // Function host()
// // --------------------------------------------------------------------------

// template <typename THost, typename TSpec>
// inline THost &
// host(ModifiedIterator<THost, TSpec> & me)
// {
//     return value(me.data_host);
// }

// template <typename THost, typename TSpec>
// inline THost const &
// host(ModifiedIterator<THost, TSpec> const & me)
// {
//     return value(me._host);
// }

// --------------------------------------------------------------------------
// Function cargo()
// --------------------------------------------------------------------------

template <typename THost, typename TSpec>
inline typename Reference<typename Cargo<ModifiedIterator<THost, TSpec> >::Type>::Type
cargo(ModifiedIterator<THost, TSpec> & me)
{
    return me._cargo;
}

template <typename THost, typename TSpec>
inline typename Reference<typename Cargo<ModifiedIterator<THost, TSpec> const>::Type>::Type
cargo(ModifiedIterator<THost, TSpec> const & me)
{
    return me._cargo;
}

// --------------------------------------------------------------------------
// Function container()
// --------------------------------------------------------------------------

template <typename THost, typename TSpec>
inline typename Container<ModifiedIterator<THost, TSpec> >::Type //no reference
container(ModifiedIterator<THost, TSpec> & me)
{
    typedef typename Container<ModifiedIterator<THost, TSpec> >::Type TContainer;
    TContainer cont(container(host(me)));
    _copyCargo(cont, me);
    return cont;
}

template <typename THost, typename TSpec>
inline typename Container<ModifiedIterator<THost, TSpec> const>::Type //no reference
container(ModifiedIterator<THost, TSpec> const & me)
{
    typedef typename Container<ModifiedIterator<THost, TSpec> const>::Type TContainer;
    TContainer cont(container(host(me)));
    _copyCargo(cont, me);
    return cont;
}

// --------------------------------------------------------------------------
// Function setContainer()
// --------------------------------------------------------------------------

template <typename TIteratorHost, typename TSpec, typename TStringHost>
inline void
setContainer(
        ModifiedIterator<TIteratorHost, TSpec> & me,
        ModifiedString<TStringHost, TSpec> & cont)
{
    setContainer(host(me), host(cont));
    _copyCargo(me, cont);
}

template <typename TIteratorHost, typename TSpec, typename TStringHost>
inline void
setContainer(
        ModifiedIterator<TIteratorHost, TSpec> & me,
        ModifiedString<TStringHost, TSpec> const & cont)
{
    setContainer(host(me), host(const_cast<ModifiedString<TStringHost, TSpec> &>(cont)));
    _copyCargo(me, cont);
}

// --------------------------------------------------------------------------
// Function assign()
// --------------------------------------------------------------------------

// TODO(holtgrew): Do!

// --------------------------------------------------------------------------
// Function value()
// --------------------------------------------------------------------------

template <typename THost, typename TSpec>
inline typename Reference<ModifiedIterator<THost, TSpec> >::Type
value(ModifiedIterator<THost, TSpec> & me)
{
    return value(host(me));
}

template <typename THost, typename TSpec>
inline typename Reference<ModifiedIterator<THost, TSpec> const>::Type
value(ModifiedIterator<THost, TSpec> const & me)
{
    return value(host(me));
}

// --------------------------------------------------------------------------
// Function operator*()
// --------------------------------------------------------------------------

template <typename THost, typename TSpec>
inline typename Reference<ModifiedIterator<THost, TSpec> >::Type
operator*(ModifiedIterator<THost, TSpec> & me)
{
    return value(me);
}

template <typename THost, typename TSpec>
inline typename Reference<ModifiedIterator<THost, TSpec> const>::Type
operator*(ModifiedIterator<THost, TSpec> const & me)
{
    return value(me);
}

// --------------------------------------------------------------------------
// Function goNext()
// --------------------------------------------------------------------------

// redefinition candidate
template <typename THost, typename TSpec>
inline void
goNext(ModifiedIterator<THost, TSpec> & me)
{
    // goNext(host(me));
    ++host(me);
}

// --------------------------------------------------------------------------
// Function operator++()
// --------------------------------------------------------------------------

template <typename THost, typename TSpec>
inline ModifiedIterator<THost, TSpec> const &
operator++(ModifiedIterator<THost, TSpec> & me)
{
    goNext(me);
    return me;
}

template <typename THost, typename TSpec>
inline ModifiedIterator<THost, TSpec>
operator++(ModifiedIterator<THost, TSpec> & me, int)
{
    ModifiedIterator<THost, TSpec> temp(me);
    goNext(me);
    return temp;
}

// --------------------------------------------------------------------------
// Function goPrevious()
// --------------------------------------------------------------------------

template <typename THost, typename TSpec>
inline void
goPrevious(ModifiedIterator<THost, TSpec> & me)
{
    goPrevious(host(me));
}

// --------------------------------------------------------------------------
// Function operator--()
// --------------------------------------------------------------------------

template <typename THost, typename TSpec>
inline ModifiedIterator<THost, TSpec> const &
operator--(ModifiedIterator<THost, TSpec> & me)
{
    goPrevious(me);
    return me;
}

template <typename THost, typename TSpec>
inline ModifiedIterator<THost, TSpec>
operator--(ModifiedIterator<THost, TSpec> & me, int)
{
    ModifiedIterator<THost, TSpec> temp(me);
    goPrevious(me);
    return temp;
}

// --------------------------------------------------------------------------
// Function operator+=()
// --------------------------------------------------------------------------

template <typename THost, typename TSpec, typename TDelta>
inline ModifiedIterator<THost, TSpec> &
operator += (ModifiedIterator<THost, TSpec> & me, TDelta delta)
{
    host(me) += delta;
    return me;
}

// --------------------------------------------------------------------------
// Function operator+()
// --------------------------------------------------------------------------

template <typename THost, typename TSpec, typename TDelta>
inline ModifiedIterator<THost, TSpec>
operator+(ModifiedIterator<THost, TSpec> const & me, TDelta delta)
{
    ModifiedIterator<THost, TSpec> temp_(me);
    temp_ += delta;
    return temp_;
}

// --------------------------------------------------------------------------
// Function operator-=()
// --------------------------------------------------------------------------

template <typename THost, typename TSpec, typename TDelta>
inline ModifiedIterator<THost, TSpec> &
operator-=(ModifiedIterator<THost, TSpec> & me, TDelta delta)
{
    host(me) -= delta;
    return me;
}

// --------------------------------------------------------------------------
// Function operator-()
// --------------------------------------------------------------------------

template <typename THost, typename TSpec, typename TDelta>
inline ModifiedIterator<THost, TSpec>
operator-(ModifiedIterator<THost, TSpec> const & me, TDelta delta)
{
    ModifiedIterator<THost, TSpec> temp_(me);
    temp_ -= delta;
    return temp_;
}

template <typename THost, typename TSpec>
inline typename Difference< ModifiedIterator<THost, TSpec> >::Type
operator-(ModifiedIterator<THost, TSpec> const & a, ModifiedIterator<THost, TSpec> const & b)
{
    return host(a) - host(b);
}

// --------------------------------------------------------------------------
// Function goBegin()
// --------------------------------------------------------------------------

// (weese:) the default implementations should do the same

//template <typename THost, typename TSpec, typename TContainer>
//inline void
//goBegin(ModifiedIterator<THost, TSpec> & me,
//        TContainer & container)
//{
//    host(me) = begin(host(container));
//}
//
//template <typename THost, typename TSpec>
//inline void
//goBegin(ModifiedIterator<THost, TSpec> & me)
//{
//    goBegin(me, container(me));
//}

// --------------------------------------------------------------------------
// Function goEnd()
// --------------------------------------------------------------------------

//template <typename THost, typename TSpec, typename TContainer>
//inline void
//goEnd(ModifiedIterator<THost, TSpec> & me,
//      TContainer & container)
//{
//    host(me) = end(host(container));
//}
//
//template <typename THost, typename TSpec>
//inline void
//goEnd(ModifiedIterator<THost, TSpec> & me)
//{
//    goEnd(me, container(me));
//}

// --------------------------------------------------------------------------
// Function position()
// --------------------------------------------------------------------------

template <typename THost, typename TSpec>
inline typename Position<ModifiedIterator<THost, TSpec> const>::Type
position(ModifiedIterator<THost, TSpec> const & me)
{
    return position(host(me));
}

template <typename THost, typename TSpec, typename TContainer>
inline typename Position<ModifiedIterator<THost, TSpec> const>::Type
position(ModifiedIterator<THost, TSpec> const & me, TContainer const &cont)
{
    return position(host(me), cont);
}

// --------------------------------------------------------------------------
// Function operator==()
// --------------------------------------------------------------------------

template <typename THost, typename TSpec>
inline bool
operator == (ModifiedIterator<THost, TSpec> const & a, ModifiedIterator<THost, TSpec> const & b)
{
    return host(a) == host(b);
}

// --------------------------------------------------------------------------
// Function operator!=()
// --------------------------------------------------------------------------

template <typename THost, typename TSpec>
inline bool
operator != (ModifiedIterator<THost, TSpec> const & a, ModifiedIterator<THost, TSpec> const & b)
{
    return !(a == b);
}

// --------------------------------------------------------------------------
// Function operator<()
// --------------------------------------------------------------------------

template <typename THost, typename TSpec>
inline bool
operator<(ModifiedIterator<THost, TSpec> const & a, ModifiedIterator<THost, TSpec> const & b)
{
    return host(a) < host(b);
}

// --------------------------------------------------------------------------
// Function operator>()
// --------------------------------------------------------------------------

template <typename THost, typename TSpec>
inline bool
operator>(ModifiedIterator<THost, TSpec> const & a, ModifiedIterator<THost, TSpec> const & b)
{
    return b < a;
}

// --------------------------------------------------------------------------
// Function atBegin()
// --------------------------------------------------------------------------

template <typename THost, typename TSpec, typename TContainer>
inline bool
atBegin(ModifiedIterator<THost, TSpec> & me,
        TContainer const & container)
{
    return atBegin(const_cast<ModifiedIterator<THost, TSpec> const &>(me), container);
}

template <typename THost, typename TSpec, typename TContainer>
inline bool
atBegin(ModifiedIterator<THost, TSpec> const & me,
        TContainer const & container)
{
    return atBegin(host(me), container);
}

template <typename THost, typename TSpec>
inline bool
atBegin(ModifiedIterator<THost, TSpec> & me)
{
    return atBegin(const_cast<ModifiedIterator<THost, TSpec> const &>(me));
}

template <typename THost, typename TSpec>
inline bool
atBegin(ModifiedIterator<THost, TSpec> const & me)
{
    return atBegin(host(me));
}

// --------------------------------------------------------------------------
// Function atEnd()
// --------------------------------------------------------------------------

template <typename THost, typename TSpec, typename TContainer>
inline bool
atEnd(ModifiedIterator<THost, TSpec> & me,
      TContainer const & container)
{
    return atEnd(const_cast<ModifiedIterator<THost, TSpec> const &>(me), container);
}

template <typename THost, typename TSpec, typename TContainer>
inline bool
atEnd(ModifiedIterator<THost, TSpec> const & me,
      TContainer const & container)
{
    return atEnd(host(me), container);
}

template <typename THost, typename TSpec>
inline bool
atEnd(ModifiedIterator<THost, TSpec> & me)
{
    return atEnd(const_cast<ModifiedIterator<THost, TSpec> const &>(me));
}

template <typename THost, typename TSpec>
inline bool
atEnd(ModifiedIterator<THost, TSpec> const & me)
{
    return atEnd(host(me));
}

}

// Adapt SeqAn modified to std.
namespace std
{
    template<typename THost, typename TSpec>
    struct iterator_traits<seqan::ModifiedIterator<THost, TSpec> >
    {
        typedef seqan::ModifiedIterator<THost, TSpec> TIter;

        typedef random_access_iterator_tag iterator_category;
        typedef typename seqan::Value<TIter>::Type value_type;
        typedef typename seqan::Difference<TIter>::Type difference_type;
        typedef typename seqan::Value<TIter>::Type * pointer;
        typedef typename seqan::Reference<TIter>::Type reference;
    };
}

#endif
