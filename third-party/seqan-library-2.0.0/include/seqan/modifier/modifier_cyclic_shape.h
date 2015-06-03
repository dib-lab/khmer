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
// Author: Sascha Meiers <meiers@inf.fu-berlin.de>
// ==========================================================================
// Modified Iterator and Modified String for CyclicShapes
// ==========================================================================

#ifndef SEQAN_HEADER_MODIFIER_SHAPE_H
#define SEQAN_HEADER_MODIFIER_SHAPE_H


namespace seqan {

// NOTE(meiers): Only a few functions are documented, the rest should be derived

// ==========================================================================
// Classes, Enums, Typedefs
// ==========================================================================

// --------------------------------------------------------------------------
// Class ModCyclicShape Iterator
// --------------------------------------------------------------------------

template<typename TSpec>
struct ModCyclicShape
{};

/*!
 * @class ModCyclicShapeModifiedIterator ModCyclicShape ModifiedIterator
 * @extends ModifiedIterator
 * @headerfile <seqan/modifier.h>
 *
 * @brief Iterates over a string leaving out don't-care-positions defined in
 *      CyclicShape.
 *
 * @signature template <typename THost, typename TSpec>
 *            class ModifiedIterator<THost, ModCyclicShape<CyclicShape<TSpec> > >;
 *
 * @tparam THost The host iterator.
 * @tparam TSpec The specialization of @link CyclicShape @endlink.
 * @tparam TCargo Cargo type of ModCyclicShape
 * @tparam TSize Size type of ModCyclicShape
 *
 * Using ModCyclicShape you can <i>mask</i> characters at certain positions of a text.
 * A CylcicShape defines these positions (so called don't-care-positions) and the
 * CyclicShape is repeated till the end of the string. The iterator behaves as if
 * the don't-care-positions had been deleted without actually copying the string.
 *
 * @snippet demos/cyclic_shape_snippets.cpp Define CyclicShape Modified Iterator
 *
 * @see ModCyclicShapeModifiedString
 * @see CyclicShape
 *
 * @note ModCyclicShape Modifier should only be used to <b>view/read a sequence</b>,
 *       never to change the underlying sequence. Some operations will work, e.g.
 *       you can use the <tt>operator*()</tt> of the ModifiedIterator to change
 *       a character, but other operations like <tt>replace</tt> or <tt>append</tt>
 *       are unsafe or even conceptionally wrong.
 */

template<typename THost, typename TSpec>
class ModifiedIterator<THost, ModCyclicShape<CyclicShape<TSpec> > >
{
public:
    typedef typename Cargo<ModifiedIterator>::Type TCargo_;
    typedef CyclicShape<TSpec>  TCyclicShape;

    THost _host;

    /*!
     * @var TCargo ModCyclicShapeModifiedIterator::_cargo;
     * @brief Copy of a Cargo object, which is a CyclicShape.
     *
     * The Cargo is a CyclicShape. This can be very heavy in case of GenericCyclicShape,
     * so prefer the lightweight FixedCyclicShape.
     */
    TCargo_ _cargo;

    /*!
     * @var TSize ModCyclicShapeModifiedIterator::_idx;
     * @brief Internal position inside <tt>diffs</tt> of the CylicShape.
     */
    typename Size<TCyclicShape>::Type _idx;

    /*!
     * @fn ModCyclicShapeModifiedIterator::ModifiedIterator
     * @brief The constructor.
     *
     * @signature ModifiedIterator::ModifiedIterator();
     * @signature ModifiedIterator::ModifiedIterator(host);
     * @signature ModifiedIterator::ModifiedIterator(otherModIter);
     * @param[in] host Host iterator to set.
     * @param[in] otherModIter ModifiedIterator which may differ in its host type.
     *            Copy construction.
     *
     * The safe way to create this ModifiedIterator is to use the <tt>begin</tt> and <tt>end</tt>
     * functions of the ModifiedString or to copy construct.
     * When directly initialized with a host iterator,
     * the CyclicShape is not yet taken into account and you have to set it up yourself,
     * for example a CyclicShape starting with a 01... will have skip the first position.
     * So you would have to reset the host iterator using <tt>++host(myCyclicShape)</tt>.
     * Use <tt>cargo()</tt> to assign a shape to the iterator.
     */
    ModifiedIterator() :
        _idx(0)
    {}

    template<typename TOtherHost>
    ModifiedIterator(ModifiedIterator<TOtherHost, ModCyclicShape<CyclicShape<TSpec> > > const & origin) :
        _host(origin._host), _cargo(origin._cargo), _idx(origin._idx)
    {}

    template<typename T>
    explicit
    ModifiedIterator(T const & host) :
        _host(host), _idx(0)
    {}
};


/*!
 * @class ModCyclicShapeModifiedString ModCyclicShape ModifiedString
 * @extends ModifiedString
 * @headerfile <seqan/modifier.h>
 *
 * @brief A string leaving out don't-care-positions defined in CyclicShape.
 *
 * @signature template <typename THost, typename TSpec>
 *            class ModifiedString<THost, ModCyclicShape<CyclicShape<TSpec> > >;
 *
 * @tparam THost The host string.
 * @tparam TSpec The specialization of CyclicShape.
 *
 * Applies a CyclicShape to a string. The resulting modified string skips don't-care positions
 * and can thus be shorter than the original string. Of course the underlying string is not changed
 * by the modifier.
 *
 * @section Examples
 *
 * Use case for GenericCyclicShapes:
 *
 * @snippet demos/cyclic_shape_snippets.cpp Define GenericCyclicShape Modified String
 *
 * Use case for FixedCyclicShapes:
 *
 * @snippet demos/cyclic_shape_snippets.cpp Define FixedCyclicShape Modified String
 *
 *
 * @see CyclicShape
 * @see ModCyclicShapeModifiedIterator
 *
 * @note ModCyclicShape Modifier should only be used to <b>view/read a sequence</b>,
 *       never to change the underlying sequence. Some operations will work, e.g.
 *       you can use the <tt>operator*()</tt> of the ModifiedIterator to change
 *       a character, but other operations like <tt>replace</tt> or <tt>append</tt>
 *       are unsafe or even conceptionally wrong.
 */

template<typename THost, typename TSpec>
class ModifiedString<THost, ModCyclicShape<CyclicShape<TSpec> > >
{
public:
    typedef typename Pointer_<THost>::Type                  THostPointer_;
    typedef typename Cargo<ModifiedString>::Type            TCargo_;
    typedef CyclicShape<TSpec>                              TCyclicShape;

    mutable THostPointer_ _host;
    TCargo_ _cargo;

    // TODO(meiers): Once the Modified String has been propperly documented, document
    //               the special constructors for ModCyclicShape Modified Strings

    // Default constructor.
    ModifiedString()
    {}

    // Construct with the actual host.
    explicit
    ModifiedString(typename Parameter_<THost>::Type host) :
        _host(_toPointer(host))
    {}

    // Construct with the actual host; variant with shape.
    explicit
    ModifiedString(typename Parameter_<THost>::Type host, TCyclicShape const & shape) :
        _host(_toPointer(host)), _cargo(shape)
    {}

    // Constructor for creating a ModifiedString with const host with a non-const host.
    template<typename THost_>
    explicit
    ModifiedString(THost_ & host,
                   SEQAN_CTOR_ENABLE_IF(IsConstructible<THost, THost_>)) :
        _host(_toPointer(host))
    {
        ignoreUnusedVariableWarning(dummy);
    }

    // Constructor for creating a ModifiedString with const host with a non-const host; variant with shape.
    template<typename THost_>
    explicit ModifiedString(THost_ & host,
                            TCyclicShape const & shape,
                            SEQAN_CTOR_ENABLE_IF(IsConstructible<THost, THost_>)) :
        _host(_toPointer(host)), _cargo(shape)
    {
        ignoreUnusedVariableWarning(dummy);
    }

#ifdef SEQAN_CXX11_STANDARD

    // Constructor for innermost type; hand down to _host which is a ModifiedString itself.  Non-const variant.
    template<typename THost_>
    explicit
    ModifiedString(THost_ && host,
                   SEQAN_CTOR_ENABLE_IF(IsAnInnerHost<
                                            typename RemoveReference<THost>::Type,
                                            typename RemoveReference<THost_>::Type>)) :
        _host(host)                     // TODO(meiers): need std::forward<THost_>(host) here?
    {
        ignoreUnusedVariableWarning(dummy);
    }

    // Constructor for innermost type; hand down to _host which is a ModifiedString itself.
    // Non-const variant with shape.
    template<typename THost_>
    explicit
    ModifiedString(THost_ && host,
                   TCyclicShape const & shape,
                   SEQAN_CTOR_ENABLE_IF(IsAnInnerHost<
                                            typename RemoveReference<THost>::Type,
                                            typename RemoveReference<THost_>::Type>)) :
        _host(host), _cargo(shape)
    {
        ignoreUnusedVariableWarning(dummy);
    }

#else // SEQAN_CXX11_STANDARD

    // Constructor for innermost type; hand down to _host which is a ModifiedString itself.  Non-const variant.
    template<typename THost_>
    explicit
    ModifiedString(THost_ & host,
                   SEQAN_CTOR_ENABLE_IF(IsAnInnerHost<THost, THost_>)) :
        _host(host)
    {
        ignoreUnusedVariableWarning(dummy);
    }

    // Constructor for innermost type; hand down to _host which is a ModifiedString itself.  Const variant.
    template<typename THost_>
    explicit
    ModifiedString(THost_ const & host,
                   SEQAN_CTOR_ENABLE_IF(IsAnInnerHost<THost, THost_ const>)) :
        _host(host)
    {
        ignoreUnusedVariableWarning(dummy);
    }

    // Constructor for innermost type; hand down to _host which is a ModifiedString itself.  Non-const variant with
    // shape.
    template<typename THost_>
    explicit
    ModifiedString(THost_ & host,
                   TCyclicShape const & shape,
                   SEQAN_CTOR_ENABLE_IF(IsAnInnerHost<THost, THost_>)) :
        _host(host), _cargo(shape)
    {
        ignoreUnusedVariableWarning(dummy);
    }

    // Constructor for innermost type; hand down to _host which is a ModifiedString itself.
    // Const variant with shape.
    template<typename THost_>
    explicit
    ModifiedString(THost_ const & host,
                   TCyclicShape const & shape,
                   SEQAN_CTOR_ENABLE_IF(IsAnInnerHost<THost, THost_ const>)) :
        _host(host), _cargo(shape)
    {
        ignoreUnusedVariableWarning(dummy);
    }

#endif // SEQAN_CXX11_STANDARD

    template<typename TPos>
    inline typename Reference<ModifiedString>::Type
        operator[] (TPos pos)
    {
        return value(*this, pos);
    }

    template<typename TPos>
    inline typename Reference<ModifiedString const>::Type
        operator[] (TPos pos) const
    {
        return value(*this, pos);
    }
};

// ==========================================================================
// Metafunctions
// ==========================================================================

// --------------------------------------------------------------------------
// Metafunction Cargo                       [ModCyclicShape ModifiedIterator]
// --------------------------------------------------------------------------

/*!
 * @mfn ModCyclicShapeModifiedIterator#Cargo
 * @headerfile <seqan/modifier.h>
 * @brief Cargo of ModCyclicShape ModCyclicShape Modified Iterator and ModCyclicShape Modified String.
 * @signature Cargo<ModifiedIterator<THost, ModCyclicShape<CyclicShape<TSpec> > > >::Type;
 * @tparam THost Host iterator of ModifiedIterator.
 * @tparam TSpec Specialization of CyclicShape.
 * @return Cargo type of ModCyclicShape Modified Iterator is the CyclicShape itself.
 */

template<typename THost, typename TSpec>
struct Cargo<ModifiedIterator<THost, ModCyclicShape<CyclicShape<TSpec> > > >
{
    typedef CyclicShape<TSpec> Type;
};

// --------------------------------------------------------------------------
// Metafunction Cargo                         [ModCyclicShape ModifiedString]
// --------------------------------------------------------------------------

/*!
 * @mfn ModCyclicShapeModifiedString#Cargo
 * @headerfile <seqan/modifier.h>
 * @brief Cargo of ModCyclicShape ModCyclicShape Modified Iterator and ModCyclicShape Modified String.
 * @signature Cargo<ModifiedString<THost, ModCyclicShape<CyclicShape<TSpec> > > >::Type;
 * @tparam THost Host container of ModifiedString.
 * @tparam TSpec Specialization of CyclicShape.
 * @return Cargo type of ModCyclicShape Modified String is the CyclicShape itself.
 */

template<typename THost, typename TSpec>
struct Cargo<ModifiedString<THost, ModCyclicShape<CyclicShape<TSpec> > > >
{
    typedef CyclicShape<TSpec> Type;
};

// --------------------------------------------------------------------------
// Metafunction Iterator                      [ModCyclicShape ModifiedString]
// --------------------------------------------------------------------------

template<typename THost, typename TSpec>
struct Iterator<ModifiedString<THost, ModCyclicShape<CyclicShape<TSpec> > >, Standard>
{
    typedef ModifiedIterator<typename Iterator<THost, Standard>::Type,
        ModCyclicShape<CyclicShape<TSpec> > > Type;
};

// TODO(meiers): Rooted ?

template<typename THost, typename TSpec>
struct Iterator<ModifiedString<THost, ModCyclicShape<CyclicShape<TSpec> > > const, Standard>
{
    typedef ModifiedIterator<typename Iterator<THost const, Standard>::Type,
        ModCyclicShape<CyclicShape<TSpec> > > Type;
};

// --------------------------------------------------------------------------
// Metafunction DefaultIteratorSpec           [ModCyclicShape ModifiedString]
// --------------------------------------------------------------------------

template<typename THost, typename TSpec>
struct DefaultIteratorSpec<ModifiedString<THost, ModCyclicShape<CyclicShape<TSpec> > > >
{
    typedef Standard Type;
};

// ==========================================================================
// Functions
// ==========================================================================

// --------------------------------------------------------------------------
// Function begin()                           [ModCyclicShape ModifiedString]
// --------------------------------------------------------------------------

template<typename THost, typename TSpec, typename TTagSpec>
inline typename Iterator<ModifiedString<THost, ModCyclicShape<CyclicShape<TSpec> > > const,
    Tag<TTagSpec> const>::Type
begin(ModifiedString<THost, ModCyclicShape<CyclicShape<TSpec> > > const & me, Tag<TTagSpec> const tag_)
{
    typedef typename Iterator<ModifiedString<THost, ModCyclicShape<CyclicShape<TSpec> > > const,
        Tag<TTagSpec> const>::Type TResult;
    typedef typename Size<CyclicShape<TSpec> >::Type TSize;

    TResult tmp(begin(host(me), tag_));
    _copyCargo(tmp, me);
    host(tmp) += (TSize)cargo(me).loffset;
    return tmp;
}

template<typename THost, typename TSpec, typename TTagSpec>
inline typename Iterator<ModifiedString<THost, ModCyclicShape<CyclicShape<TSpec> > >,
    Tag<TTagSpec> const>::Type
begin(ModifiedString<THost, ModCyclicShape<CyclicShape<TSpec> > > & me, Tag<TTagSpec> const tag_)
{
    typedef typename Iterator<ModifiedString<THost, ModCyclicShape<CyclicShape<TSpec> > >,
        Tag<TTagSpec> const>::Type TResult;
    typedef typename Size<CyclicShape<TSpec> >::Type TSize;

    TResult tmp(begin(host(me), tag_));
    _copyCargo(tmp, me);
    host(tmp) += (TSize)cargo(me).loffset;
    return tmp;
}

// --------------------------------------------------------------------------
// Function end()                             [ModCyclicShape ModifiedString]
// --------------------------------------------------------------------------

/*!
 * @fn ModCyclicShapeModifiedString#end
 * @headerfile <seqan/modifier.h>
 * @brief Returns an iterator to the end of the container.
 * @signature TIterator end(modStr[, tag]);
 * @tparam TIterator ModCyclicShape Modified Iterator to be returned.
 * @param[in] modStr ModCyclicShape Modified String.
 * @param[in] An optional tag for selecting the type of the iterator.
 *            One of Standard and Rooted. When left out,
 *            @link ContainerConcept#DefaultGetIteratorSpec @endlink of <tt>modStr</tt> is used.
 * @return Iterator to the end of the container, the type is selected by #Fn<>ContainerConcept#Iterator
 *         with the given (or default) tag.
 *
 * The end iterator is not guaranteed to be at position <tt>length(modStr)</tt>, where
 * a standard iterator ends. Due to the jumping-tequniche of ModCyclicShape
 * Modified Iterator, the end iterator is set to the first care position after
 * the end of the string. The invariant is: <tt>host(end(modStr)) >= end(host(modStr))<tt>
 * where the the larger or equal operator refers to the position in the sequence.
 *
 * If for some reason you use the host iterator to check an end condition, make sure
 * you use <tt> for(...; host(modIter) >= end(host(modStr)); ...) </tt> instead of
 * <tt> for(...; host(modIter) != end(host(modStr)); ...) </tt>.
 */

template<typename THost, typename TSpec, typename TTagSpec>
inline typename Iterator<ModifiedString<THost, ModCyclicShape<CyclicShape<TSpec> > > const,
    Tag<TTagSpec> const>::Type
end(ModifiedString<THost, ModCyclicShape<CyclicShape<TSpec> > > const & me, Tag<TTagSpec> const tag_)
{
    typedef typename Iterator<ModifiedString<THost, ModCyclicShape<CyclicShape<TSpec> > > const,
        Tag<TTagSpec> const>::Type TResult;
    TResult tmp(end(host(me), tag_));
    _copyCargo(tmp, me);
    goEnd(tmp, me);
    return tmp;
}

template<typename THost, typename TSpec, typename TTagSpec>
inline typename Iterator<ModifiedString<THost, ModCyclicShape<CyclicShape<TSpec> > >,
    Tag<TTagSpec> const>::Type
end(ModifiedString<THost, ModCyclicShape<CyclicShape<TSpec> > > & me, Tag<TTagSpec> const tag_)
{
    typedef typename Iterator<ModifiedString<THost, ModCyclicShape<CyclicShape<TSpec> > >,
        Tag<TTagSpec> const>::Type TResult;
    TResult tmp(end(host(me), tag_));
    _copyCargo(tmp, me);
    goEnd(tmp, me);
    return tmp;
}

// --------------------------------------------------------------------------
// Function goNext()                        [ModCyclicShape ModifiedIterator]
// --------------------------------------------------------------------------

template<typename THost, typename TSpec>
inline void
goNext(ModifiedIterator<THost, ModCyclicShape<CyclicShape<TSpec> > > & me)
{
    host(me) += cargo(me).diffs[me._idx];
    if (++me._idx == weight(cargo(me)))
        me._idx = 0;
}

// --------------------------------------------------------------------------
// Function goPrevious()                    [ModCyclicShape ModifiedIterator]
// --------------------------------------------------------------------------

template<typename THost, typename TSpec>
inline void
goPrevious(ModifiedIterator<THost, ModCyclicShape<CyclicShape<TSpec> > > & me)
{
    if (me._idx == 0)
        me._idx = weight(cargo(me));
    --me._idx;
    host(me) -= cargo(me).diffs[me._idx];
}

// --------------------------------------------------------------------------
// Function goBegin()                       [ModCyclicShape ModifiedIterator]
// --------------------------------------------------------------------------

template<typename THost, typename TSpec>
inline void
goBegin(ModifiedIterator<THost, ModCyclicShape<CyclicShape<TSpec> > > & me)
{
    goBegin(host(me));
    me._idx = 0;
    host(me) += cargo(me).loffset;
}

template<typename THost, typename TSpec, typename TContainer>
inline void
goBegin(ModifiedIterator<THost, ModCyclicShape<CyclicShape<TSpec> > > & me, TContainer & cont)
{
    goBegin(host(me), host(cont));
    me._idx = 0;
    host(me) += cargo(me).loffset;
}

template<typename THost, typename TSpec, typename TContainer>
inline void
goBegin(ModifiedIterator<THost, ModCyclicShape<CyclicShape<TSpec> > > & me, TContainer const & cont)
{
    goBegin(host(me), host(cont));
    me._idx = 0;
    host(me) += cargo(me).loffset;
}

// --------------------------------------------------------------------------
// Function goEnd()                         [ModCyclicShape ModifiedIterator]
// --------------------------------------------------------------------------

template<typename THost, typename TSpec>
inline void
goEnd(ModifiedIterator<THost, ModCyclicShape<CyclicShape<TSpec> > > & me)
{
    goBegin(host(me));

    THost _end(host(me));
    goEnd(_end);

    typedef typename Position<THost>::Type TPos;
    TPos len = _end - host(me);
    TPos pos = cargo(me).span * (len / cargo(me).span) + cargo(me).loffset;

    for (me._idx = 0; pos < len; me._idx = (me._idx + 1) % weight(cargo(me)))
        pos += cargo(me).diffs[me._idx];

    host(me) += pos;
}

template<typename THost, typename TSpec, typename TContainer>
inline void
goEnd(ModifiedIterator<THost, ModCyclicShape<CyclicShape<TSpec> > > & me, TContainer & cont)
{
    goBegin(host(me), host(cont));

    typedef typename Position<THost>::Type TPos;
    TPos len = length(host(cont));
    TPos pos = cargo(me).span * (len / cargo(me).span) + cargo(me).loffset;

    for (me._idx = 0; pos < len; me._idx = (me._idx + 1) % weight(cargo(me)))
        pos += cargo(me).diffs[me._idx];

    host(me) += pos;
}

template<typename THost, typename TSpec, typename TContainer>
inline void
goEnd(ModifiedIterator<THost, ModCyclicShape<CyclicShape<TSpec> > > & me, TContainer const & cont)
{
    goBegin(host(me), host(cont));

    typedef typename Position<THost>::Type TPos;
    TPos len = length(host(cont));
    TPos pos = cargo(me).span * (len / cargo(me).span) + cargo(me).loffset;

    for (me._idx = 0; pos < len; me._idx = (me._idx + 1) % weight(cargo(me)))
        pos += cargo(me).diffs[me._idx];

    host(me) += pos;
}

// --------------------------------------------------------------------------
// Function operator+=()                    [ModCyclicShape ModifiedIterator]
// --------------------------------------------------------------------------

template<typename THost, typename TSpec, typename TDelta>
inline ModifiedIterator<THost, ModCyclicShape<CyclicShape<TSpec> > > &
operator+=(ModifiedIterator<THost, ModCyclicShape<CyclicShape<TSpec> > >&me, TDelta delta)
{
    if (isNegative(delta))
    {
        me -= -(typename MakeSigned<TDelta>::Type)delta;
    }
    else
    {
        // jump full patterns
        host(me) += (delta / weight(cargo(me))) * cargo(me).span;

        // number of jumps in dist that remain
        for (delta = delta % weight(cargo(me)); delta != 0; --delta)
            goNext(me);
    }
    return me;
}

// --------------------------------------------------------------------------
// Function operator-=()                    [ModCyclicShape ModifiedIterator]
// --------------------------------------------------------------------------

template<typename THost, typename TSpec, typename TDelta>
inline ModifiedIterator<THost, ModCyclicShape<CyclicShape<TSpec> > > &
operator-=(ModifiedIterator<THost, ModCyclicShape<CyclicShape<TSpec> > >&me, TDelta delta)
{
    if (isNegative(delta))
    {
        me += -(typename MakeSigned<TDelta>::Type)delta;
    }
    else
    {
        // jump full patterns
        host(me) -= (delta / weight(cargo(me))) * cargo(me).span;

        // number of jumps in dist that remain
        for (delta = delta % weight(cargo(me)); delta != 0; --delta)
            goPrevious(me);
    }
    return me;
}

// --------------------------------------------------------------------------
// Function atBegin()                       [ModCyclicShape ModifiedIterator]
// --------------------------------------------------------------------------

template<typename THost, typename TSpec>
inline bool
atBegin(ModifiedIterator<THost, ModCyclicShape<CyclicShape<TSpec> > > & me)
{
    THost _beg(host(me));
    goBegin(_beg);
    _beg += cargo(me).DIFFS[0];
    return host(me) == _beg && me._idx == 0;
}

// --------------------------------------------------------------------------
// Function atEnd()                         [ModCyclicShape ModifiedIterator]
// --------------------------------------------------------------------------

template<typename THost, typename TSpec>
inline bool
atEnd(ModifiedIterator<THost, ModCyclicShape<CyclicShape<TSpec> > > & me)
{
    ModifiedIterator<THost, ModCyclicShape<CyclicShape<TSpec> > > _tmp(me);
    goEnd(_tmp);
    return _tmp == me;
}

// --------------------------------------------------------------------------
// Function operator==()                        [ModCyclicShape ModifiedIterator]
// --------------------------------------------------------------------------

template<typename THost, typename TSpec>
inline bool
    operator == (ModifiedIterator<THost, ModCyclicShape<CyclicShape<TSpec> > > const & a,
                 ModifiedIterator<THost, ModCyclicShape<CyclicShape<TSpec> > > const & b)
{
    // TODO(meiers): check whether cargo is equal (but there is no operator == for CyclicShape yet)

    return host(a) == host(b) && a._idx == b._idx;
}

// --------------------------------------------------------------------------
// Function operator<()                     [ModCyclicShape ModifiedIterator]
// --------------------------------------------------------------------------

// uses general version from modifier_iterator.h

// --------------------------------------------------------------------------
// Function operator-()                     [ModCyclicShape ModifiedIterator]
// --------------------------------------------------------------------------

template<typename THost, typename TSpec>
inline typename Difference<ModifiedIterator<THost, ModCyclicShape<CyclicShape<TSpec> > > >::Type
operator-(ModifiedIterator<THost, ModCyclicShape<CyclicShape<TSpec> > > const & a,
                ModifiedIterator<THost, ModCyclicShape<CyclicShape<TSpec> > > const & b)
{
    typedef typename Difference<ModifiedIterator<THost, ModCyclicShape<CyclicShape<TSpec> > > >::Type TDiff;
    if (host(a) == host(b))
        return 0;

    if (a > b)
    {
        THost _a(host(a));
        THost _b(host(b));
        TDiff diff = ((_a - _b) / cargo(a).span) * weight(cargo(a));
        diff += a._idx;
        diff -= b._idx;
        return diff;
    }
    else
    {
        return -(b - a);
    }
}

// --------------------------------------------------------------------------
// Function length()                          [ModCyclicShape ModifiedString]
// --------------------------------------------------------------------------

template<typename THost, typename TSpec>
inline typename Size<ModifiedString<THost, ModCyclicShape<TSpec> > >::Type
length(ModifiedString<THost, ModCyclicShape<TSpec> > const & me)
{
    return end(me, Standard()) - begin(me, Standard());
}

// --------------------------------------------------------------------------
// Function position()                      [ModCyclicShape ModifiedIterator]
// --------------------------------------------------------------------------

template<typename THost, typename TSpec>
inline typename Position<ModifiedIterator<THost, ModCyclicShape<CyclicShape<TSpec> > > const>::Type
position(ModifiedIterator<THost, ModCyclicShape<CyclicShape<TSpec> > > const & me)
{
    ModifiedIterator<THost, ModCyclicShape<CyclicShape<TSpec> > > _tmp(me);
    goEnd(_tmp);
    return _tmp - me;
}

// TODO(meiers): What about setPosition() ??

} // namespace

#endif
