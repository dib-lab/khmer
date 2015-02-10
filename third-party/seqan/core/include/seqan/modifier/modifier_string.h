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
// Author: David Weese <david.weese@fu-berlin.de>
// Author: Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>
// ==========================================================================

#ifndef SEQAN_MODIFIER_MODIFIER_STRING_H_
#define SEQAN_MODIFIER_MODIFIER_STRING_H_

namespace seqan
{

// ==========================================================================
// Forwards
// ==========================================================================

template <typename T> struct InnermostHost_;

// ==========================================================================
// Classes
// ==========================================================================

// --------------------------------------------------------------------------
// Class ModifierString
// --------------------------------------------------------------------------

/**
.Class.ModifiedString:
..summary:Allows to modify arbitrary strings by specializing what differs from an origin.
..cat:Modifier
..signature:ModifiedString<THost[, TSpec]>
..param.THost:Original sequence type.
...type:Concept.ContainerConcept
..param.TSpec:The modifier type.
...metafunction:Metafunction.Spec
..implements:Concept.ContainerConcept
..remarks:$THost$ can also be a modified string, so you can create custom strings by combining predefined ones.
..example.file:demos/modifier/modified_string.cpp
..example.text:The output is as follows:
..example.output:
TATACGCGAAAA
AAAAGCGCATAT

TATACGCGTTTT
TTTTGCGCATAT
..include:seqan/modifier.h
*/

template <typename THost, typename TSpec = void>
class ModifiedString
{
public:
    typedef typename Pointer_<THost>::Type       THostPointer_;
    typedef typename Cargo<ModifiedString>::Type TCargo_;

    typedef typename InnermostHost_<ModifiedString>::Type TInnermostHost_;

    mutable THostPointer_ _host;
    TCargo_ _cargo;

    // Default constructor.
    ModifiedString() : _host(), _cargo()
    {}

    // Construct with the actual host.
    explicit
    ModifiedString(THost & host) : _host(_toPointer(host)), _cargo()
    {}

    // Constructor for creating a ModifiedString with const host with a non-const host.
    template <typename THost_>
    explicit
    ModifiedString(THost_ const & host,
                   SEQAN_CTOR_ENABLE_IF(IsSameType<THost, THost_>)) :
            _host(_toPointer(host)), _cargo()
    {
        ignoreUnusedVariableWarning(dummy);
    }

    // Constructor for innermost type; hand down to _host which is a ModifiedString itself.  Non-const variant.
    template <typename THost_>
    explicit
    ModifiedString(THost_ & host,
                   SEQAN_CTOR_ENABLE_IF(And<Not<IsSameType<TInnermostHost_, THost> >,
                                            IsSameType<TInnermostHost_, THost_> >)) :
            _host(host), _cargo()
    {
        ignoreUnusedVariableWarning(dummy);
    }

    // Constructor for innermost type; hand down to _host which is a ModifiedString itself.  Const variant.
    template <typename THost_>
    explicit
    ModifiedString(THost_ const & host,
                   SEQAN_CTOR_ENABLE_IF(And<Not<IsSameType<TInnermostHost_, THost> >,
                                            IsSameType<TInnermostHost_, THost_> >)) :
            _host(host), _cargo()
    {
        ignoreUnusedVariableWarning(dummy);
    }
};

// ==========================================================================
// Metafunctions
// ==========================================================================

// --------------------------------------------------------------------------
// Metafunction InnermostHost_
// --------------------------------------------------------------------------

// TODO(holtgrew): Test me!

// This metafunction returns the innermost host type for a ModifiedString cascade.

// Recurse down const/non-const.

template <typename THost, typename TInnerSpec, typename TOuterSpec>
struct InnermostHost_<ModifiedString<ModifiedString<THost, TInnerSpec> const, TOuterSpec> >
{
    typedef THost Type;
};

template <typename THost, typename TInnerSpec, typename TOuterSpec>
struct InnermostHost_<ModifiedString<ModifiedString<THost, TInnerSpec>, TOuterSpec> >
{
    typedef THost Type;
};

// Recursion stop.

template <typename THost, typename TSpec>
struct InnermostHost_<ModifiedString<THost, TSpec> const>
{
    typedef THost Type;
};

template <typename THost, typename TSpec>
struct InnermostHost_<ModifiedString<THost, TSpec> >
{
    typedef THost Type;
};

// --------------------------------------------------------------------------
// Metafunction Spec
// --------------------------------------------------------------------------

template < typename THost, typename TSpec >
struct Spec< ModifiedString<THost, TSpec> >
{
    typedef TSpec Type;
};

template < typename THost, typename TSpec >
struct Spec< ModifiedString<THost, TSpec> const >
{
    typedef TSpec Type;
};

// --------------------------------------------------------------------------
// Metafunction Value
// --------------------------------------------------------------------------

// use Value, GetValue, Reference, Size, ... from corresponding iterator
template < typename THost, typename TSpec >
struct Value< ModifiedString<THost, TSpec> >:
    Value< typename Iterator< ModifiedString<THost, TSpec>, Rooted >::Type > {};

template < typename THost, typename TSpec >
struct Value< ModifiedString<THost, TSpec> const >:
    Value< typename Iterator< ModifiedString<THost, TSpec> const, Rooted >::Type > {};

// --------------------------------------------------------------------------
// Metafunction GetValue
// --------------------------------------------------------------------------

template < typename THost, typename TSpec >
struct GetValue< ModifiedString<THost, TSpec> >:
    GetValue< typename Iterator< ModifiedString<THost, TSpec>, Rooted >::Type > {};

template < typename THost, typename TSpec >
struct GetValue< ModifiedString<THost, TSpec> const >:
    GetValue< typename Iterator< ModifiedString<THost, TSpec> const, Rooted >::Type > {};

// --------------------------------------------------------------------------
// Metafunction Reference
// --------------------------------------------------------------------------

template < typename THost, typename TSpec >
struct Reference< ModifiedString<THost, TSpec> >:
    Reference< typename Iterator< ModifiedString<THost, TSpec>, Rooted >::Type > {};

template < typename THost, typename TSpec >
struct Reference< ModifiedString<THost, TSpec> const >:
    Reference< typename Iterator< ModifiedString<THost, TSpec> const, Rooted >::Type > {};

// --------------------------------------------------------------------------
// Metafunction Size
// --------------------------------------------------------------------------

///.Metafunction.Size.param.T.type:Class.ModifiedString
///.Metafunction.Size.class:Class.ModifiedString

template < typename THost, typename TSpec >
struct Size< ModifiedString<THost, TSpec> >:
    Size< typename Iterator< ModifiedString<THost, TSpec>, Rooted >::Type > {};

// --------------------------------------------------------------------------
// Metafunction Position
// --------------------------------------------------------------------------

template < typename THost, typename TSpec >
struct Position< ModifiedString<THost, TSpec> >:
    Position< typename Iterator< ModifiedString<THost, TSpec>, Rooted >::Type > {};

// --------------------------------------------------------------------------
// Metafunction Difference
// --------------------------------------------------------------------------

template < typename THost, typename TSpec >
struct Difference< ModifiedString<THost, TSpec> >:
    Difference< typename Iterator< ModifiedString<THost, TSpec>, Rooted >::Type > {};

// --------------------------------------------------------------------------
// Metafunction Iterator
// --------------------------------------------------------------------------

// TODO(holtgrew): Should the result of Iterator<> be a const iterator for const ModifiedString objects?

///.Metafunction.Iterator.param.T.type:Class.ModifiedString
///.Metafunction.Iterator.class:Class.ModifiedString

template <typename THost, typename TSpec>
struct Iterator<ModifiedString<THost, TSpec>, Standard>
{
    typedef ModifiedIterator<typename Iterator<THost, Standard>::Type, TSpec> Type;
};

template <typename THost, typename TSpec >
struct Iterator<ModifiedString<THost, TSpec> const, Standard>
{
    typedef ModifiedIterator<typename Iterator<THost, Standard>::Type, TSpec> Type;
};

template <typename THost, typename TSpec>
struct Iterator<ModifiedString<THost, TSpec>, Rooted>
{
    typedef ModifiedIterator<typename Iterator<THost, Rooted>::Type, TSpec> Type;
};

template <typename THost, typename TSpec >
struct Iterator<ModifiedString<THost, TSpec> const, Rooted>
{
    typedef ModifiedIterator<typename Iterator<THost, Rooted>::Type, TSpec> Type;
};

// --------------------------------------------------------------------------
// Metafunction Host
// --------------------------------------------------------------------------

///.Metafunction.Host.param.T.type:Class.ModifiedString
///.Metafunction.Host.class:Class.ModifiedString

template <typename THost, typename TSpec >
struct Host<ModifiedString<THost, TSpec> > {
    typedef THost Type;
};

template <typename THost, typename TSpec >
struct Host<ModifiedString<THost, TSpec> const > {
    typedef THost const Type;
};

// --------------------------------------------------------------------------
// Metafunction Parameter_
// --------------------------------------------------------------------------

// don't uncomment this, it would cause the Segment() c'tor to take the address
// of a temporary copy of the ModifiedString

template <typename THost, typename TSpec >
struct Parameter_<ModifiedString<THost, TSpec> > {
    typedef ModifiedString<THost, TSpec> Type;
};

template <typename THost, typename TSpec >
struct Parameter_<ModifiedString<THost, TSpec> const > {
    typedef ModifiedString<THost, TSpec> Type;
};

// --------------------------------------------------------------------------
// Metafunction Pointer_
// --------------------------------------------------------------------------

template <typename THost, typename TSpec >
struct Pointer_<ModifiedString<THost, TSpec> >
{
    typedef ModifiedString<THost, TSpec> Type;
};

template <typename THost, typename TSpec >
struct Pointer_<ModifiedString<THost, TSpec> const > : Pointer_<ModifiedString<THost, TSpec> >
{};

// --------------------------------------------------------------------------
// Metafunction IsSequence
// --------------------------------------------------------------------------

///.Metafunction.IsSequence.param.T.type:Class.ModifiedString

template <typename THost, typename TSpec >
struct IsSequence<ModifiedString<THost, TSpec> > : True
{};

// --------------------------------------------------------------------------
// Metafunction AllowsFastRandomAccess
// --------------------------------------------------------------------------

template <typename THost, typename TSpec >
struct AllowsFastRandomAccess<ModifiedString<THost, TSpec> > : AllowsFastRandomAccess<THost>
{};

// ==========================================================================
// Functions
// ==========================================================================

// --------------------------------------------------------------------------
// Function _copyCargo
// --------------------------------------------------------------------------

template <typename TDest, typename TSource>
inline void _copyCargoImpl(TDest &, TSource &, False const)
{}

template <typename TDest, typename TSource>
inline void _copyCargoImpl(TDest & me, TSource & _origin, True const)
{
    cargo(me) = cargo(_origin);
}

template <typename TDest, typename TSource>
inline void _copyCargo(TDest & me, TSource & _origin)
{
    _copyCargoImpl(me, _origin, typename IsSameType<
                   typename RemoveConst_<typename Cargo<TDest>::Type >::Type, 
                   typename RemoveConst_<typename Cargo<TSource>::Type>::Type >::Type());
}

// --------------------------------------------------------------------------
// Function _toPointer()
// --------------------------------------------------------------------------

template <typename THost, typename TSpec>
inline typename Pointer_<ModifiedString<THost, TSpec> >::Type
_toPointer(ModifiedString<THost, TSpec> & me)
{
    return me;
}

template <typename THost, typename TSpec>
inline typename Pointer_<ModifiedString<THost, TSpec> const >::Type
_toPointer(ModifiedString<THost, TSpec> const & me)
{
    return me;
}

// --------------------------------------------------------------------------
// Function _toParameter()
// --------------------------------------------------------------------------

template <typename THost, typename TSpec>
inline typename Parameter_<ModifiedString<THost, TSpec> >::Type
_toParameter(ModifiedString<THost, TSpec> & me)
{
    return me;
}

template <typename THost, typename TSpec>
inline typename Parameter_<ModifiedString<THost, TSpec> const >::Type
_toParameter(ModifiedString<THost, TSpec> const & me)
{
    return me;
}

// --------------------------------------------------------------------------
// Function _fromPointer()
// --------------------------------------------------------------------------

// TODO(holtgrew): Replace with _toParameter()?

template <typename T>
T &
_fromPointer(T * ptr)
{
    return *ptr;
}

template <typename T>
T const &
_fromPointer(T const * ptr)
{
    return *ptr;
}

template <typename THost, typename TSpec>
ModifiedString<THost, TSpec> _fromPointer(ModifiedString<THost, TSpec> & me)
{
    return me;
}

template <typename THost, typename TSpec>
ModifiedString<THost, TSpec> _fromPointer(ModifiedString<THost, TSpec> const & me)
{
    return me;
}

// --------------------------------------------------------------------------
// Function host()
// --------------------------------------------------------------------------

template <typename THost, typename TSpec>
inline typename Parameter_<THost>::Type
host(ModifiedString<THost, TSpec> & me)
{
    return _toParameter(_fromPointer(me._host));
}

template <typename THost, typename TSpec>
inline typename Parameter_<THost>::Type
host(ModifiedString<THost, TSpec> const & me)
{
    return _toParameter(_fromPointer(me._host));
}

// --------------------------------------------------------------------------
// Function setHost()
// --------------------------------------------------------------------------

template <typename THost, typename TSpec>
inline void setHost(ModifiedString<THost, TSpec> & me, THost & host)
{
    me._host = _toPointer(host);
}

template <typename THost, typename TSpec>
inline void setHost(ModifiedString<THost, TSpec> & me, THost const & host)
{
    me._host = _toPointer(host);
}

// --------------------------------------------------------------------------
// Function cargo()
// --------------------------------------------------------------------------

template <typename THost, typename TSpec>
inline typename Reference<typename Cargo<ModifiedString<THost, TSpec> >::Type >::Type
cargo(ModifiedString<THost, TSpec> & me) 
{
    return me._cargo;
}

template <typename THost, typename TSpec>
inline typename Reference<typename Cargo<ModifiedString<THost, TSpec> const>::Type >::Type
cargo(ModifiedString<THost, TSpec> const & me) 
{
    return me._cargo;
}

// --------------------------------------------------------------------------
// Function value()
// --------------------------------------------------------------------------

template <typename THost, typename TSpec, typename TPos>
inline typename Reference<ModifiedString<THost, TSpec> >::Type 
value(ModifiedString<THost, TSpec> & me, TPos pos)
{
    return value(begin(me, Standard()) + pos);
}

template <typename THost, typename TSpec, typename TPos>
inline typename Reference<ModifiedString<THost, TSpec> const >::Type 
value(ModifiedString<THost, TSpec> const & me, TPos pos)
{
    return value(begin(me, Standard()) + pos);
}

// --------------------------------------------------------------------------
// Function length()
// --------------------------------------------------------------------------

template <typename THost, typename TSpec>
inline typename Size<ModifiedString<THost, TSpec> >::Type 
length(ModifiedString<THost, TSpec> const & me)
{
    return length(host(me));
}

// --------------------------------------------------------------------------
// Function begin()
// --------------------------------------------------------------------------

template <typename THost, typename TSpec>
inline typename Iterator<ModifiedString<THost, TSpec> const>::Type 
begin(ModifiedString<THost, TSpec> const & me)
{
    typedef typename Iterator<ModifiedString<THost, TSpec> const>::Type TResult;
    TResult tmp(begin(host(me)));
    _copyCargo(tmp, me);
    return tmp;
}

template <typename THost, typename TSpec>
inline typename Iterator<ModifiedString<THost, TSpec> >::Type 
begin(ModifiedString<THost, TSpec> & me)
{
    typedef typename Iterator<ModifiedString<THost, TSpec> >::Type TResult;
    TResult tmp(begin(host(me)));
    _copyCargo(tmp, me);
    return tmp;
}

template <typename THost, typename TSpec, typename TTagSpec>
inline typename Iterator<ModifiedString<THost, TSpec> const, Tag<TTagSpec> const>::Type 
begin(ModifiedString<THost, TSpec> const & me, Tag<TTagSpec> const tag_)
{
    typedef typename Iterator<ModifiedString<THost, TSpec> const, Tag<TTagSpec> const>::Type TResult;
    TResult tmp(begin(host(me), tag_));
    _copyCargo(tmp, me);
    return tmp;
}

template <typename THost, typename TSpec, typename TTagSpec>
inline typename Iterator<ModifiedString<THost, TSpec>, Tag<TTagSpec> const>::Type 
begin(ModifiedString<THost, TSpec> & me, Tag<TTagSpec> const tag_)
{
    typedef typename Iterator<ModifiedString<THost, TSpec>, Tag<TTagSpec> const>::Type TResult;
    TResult tmp(begin(host(me), tag_));
    _copyCargo(tmp, me);
    return tmp;
}

// --------------------------------------------------------------------------
// Function end()
// --------------------------------------------------------------------------

template < typename THost, typename TSpec >
inline typename Iterator< ModifiedString<THost, TSpec> const >::Type 
end(ModifiedString<THost, TSpec> const & me)
{
    typedef typename Iterator<ModifiedString<THost, TSpec> >::Type TResult;
    TResult tmp(end(host(me)));
    _copyCargo(tmp, me);
    return tmp;
}

template < typename THost, typename TSpec >
inline typename Iterator< ModifiedString<THost, TSpec> >::Type 
end(ModifiedString<THost, TSpec> & me)
{
    typedef typename Iterator<ModifiedString<THost, TSpec> const>::Type TResult;
    TResult tmp(end(host(me)));
    _copyCargo(tmp, me);
    return tmp;
}

template < typename THost, typename TSpec, typename TTagSpec >
inline typename Iterator< ModifiedString<THost, TSpec> const, Tag<TTagSpec> const >::Type 
end(ModifiedString<THost, TSpec> const & me, Tag<TTagSpec> const tag_)
{
    typedef typename Iterator<ModifiedString<THost, TSpec> const, Tag<TTagSpec> const>::Type TResult;
    TResult tmp(end(host(me), tag_));
    _copyCargo(tmp, me);
    return tmp;
}

template < typename THost, typename TSpec, typename TTagSpec >
inline typename Iterator< ModifiedString<THost, TSpec>, Tag<TTagSpec> const >::Type 
end(ModifiedString<THost, TSpec> & me, Tag<TTagSpec> const tag_)
{
    typedef typename Iterator<ModifiedString<THost, TSpec>, Tag<TTagSpec> const>::Type TResult;
    TResult tmp(end(host(me), tag_));
    _copyCargo(tmp, me);
    return tmp;
}

// --------------------------------------------------------------------------
// Function operator==()
// --------------------------------------------------------------------------

template <typename THost, typename TSpec, typename TRight >
inline bool
operator==(ModifiedString<THost, TSpec> const & left, 
           TRight const & right)
{
	typename Comparator<ModifiedString<THost, TSpec> >::Type _lex(left, right);
    return isEqual(_lex);
}

template <typename TLeftValue, typename THost, typename TSpec >
inline bool
operator==(TLeftValue * left,
           ModifiedString<THost, TSpec> const & right)
{
	typename Comparator<ModifiedString<THost, TSpec> >::Type _lex(left, right);
    return isEqual(_lex);
}

// --------------------------------------------------------------------------
// Function operator!=()
// --------------------------------------------------------------------------

template <typename THost, typename TSpec, typename TRight >
inline bool
operator!=(ModifiedString<THost, TSpec> const & left, 
           TRight const & right)
{
	typename Comparator<ModifiedString<THost, TSpec> >::Type _lex(left, right);
    return isNotEqual(_lex);
}

template <typename TLeftValue, typename THost, typename TSpec >
inline bool
operator!= (TLeftValue * left,
            ModifiedString<THost, TSpec> const & right)
{
	typename Comparator<ModifiedString<THost, TSpec> >::Type _lex(left, right);
    return isNotEqual(_lex);
}

// --------------------------------------------------------------------------
// Function operator<()
// --------------------------------------------------------------------------

template <typename THost, typename TSpec, typename TRight>
inline bool
operator<(ModifiedString<THost, TSpec> const & left, 
          TRight const & right)
{
	return isLess(left, right, typename DefaultPrefixOrder<ModifiedString<THost, TSpec> >::Type());
}

template <typename TLeftValue, typename THost, typename TSpec >
inline bool
operator<(TLeftValue * left,
          ModifiedString<THost, TSpec> const & right)
{
	return isLess(left, right, typename DefaultPrefixOrder<TLeftValue *>::Type());
}

// --------------------------------------------------------------------------
// Function operator<=()
// --------------------------------------------------------------------------

template <typename THost, typename TSpec, typename TRight>
inline bool
operator<=(ModifiedString<THost, TSpec> const & left, 
           TRight const & right)
{
	return isLessOrEqual(left, right, typename DefaultPrefixOrder<ModifiedString<THost, TSpec> >::Type());
}

template <typename TLeftValue, typename THost, typename TSpec >
inline bool
operator<=(TLeftValue * left,
           ModifiedString<THost, TSpec> const & right)
{
	return isLessOrEqual(left, right, typename DefaultPrefixOrder<TLeftValue *>::Type());
}

// --------------------------------------------------------------------------
// Function operator>()
// --------------------------------------------------------------------------

template <typename THost, typename TSpec, typename TRight>
inline bool
operator>(ModifiedString<THost, TSpec> const & left, 
          TRight const & right)
{
	return isGreater(left, right, typename DefaultPrefixOrder<ModifiedString<THost, TSpec> >::Type());
}

template <typename TLeftValue, typename THost, typename TSpec >
inline bool
operator>(TLeftValue * left,
          ModifiedString<THost, TSpec> const & right)
{
	return isGreater(left, right, typename DefaultPrefixOrder<TLeftValue *>::Type());
}

// --------------------------------------------------------------------------
// Function operator>=()
// --------------------------------------------------------------------------

template <typename THost, typename TSpec, typename TRight>
inline bool
operator>=(ModifiedString<THost, TSpec> const & left, 
           TRight const & right)
{
	return isGreaterOrEqual(left, right, typename DefaultPrefixOrder<ModifiedString<THost, TSpec> >::Type());
}

template <typename TLeftValue, typename THost, typename TSpec >
inline bool
operator>=(TLeftValue * left,
           ModifiedString<THost, TSpec> const & right)
{
	return isGreaterOrEqual(left, right, typename DefaultPrefixOrder<TLeftValue *>::Type());
}

// --------------------------------------------------------------------------
// Function operator<<()
// --------------------------------------------------------------------------

template < typename TStream, typename THost, typename TSpec >
inline TStream &
operator<<(TStream & target, ModifiedString<THost, TSpec> const & source)
{
    write(target, source);
    return target;
}

// --------------------------------------------------------------------------
// Function operator>>()
// --------------------------------------------------------------------------

template < typename TStream, typename THost, typename TSpec >
inline TStream &
operator>>(TStream & source, ModifiedString<THost, TSpec> & target)
{
    read(source, target);
    return source;
}

// --------------------------------------------------------------------------
// Function getObjectId()
// --------------------------------------------------------------------------

template <typename THost, typename TSpec>
inline void const *
getObjectId(ModifiedString<THost, TSpec> & me) 
{
    return getObjectId(host(me));
}

template <typename THost, typename TSpec>
inline void const *
getObjectId(ModifiedString<THost, TSpec> const & me) 
{
    return getObjectId(host(me));
}

}  // namespace seqan

#endif  // SEQAN_MODIFIER_MODIFIER_STRING_H_
