/***************************************************************************
 *  include/stxxl/bits/compat/type_traits.h
 *
 *  Part of the STXXL. See http://stxxl.sourceforge.net
 *
 *  Copyright (C) 2009-2011 Andreas Beckmann <beckmann@cs.uni-frankfurt.de>
 *
 *  Distributed under the Boost Software License, Version 1.0.
 *  (See accompanying file LICENSE_1_0.txt or copy at
 *  http://www.boost.org/LICENSE_1_0.txt)
 **************************************************************************/

#ifndef STXXL_COMPAT_TYPE_TRAITS_HEADER
#define STXXL_COMPAT_TYPE_TRAITS_HEADER

#include <stxxl/bits/config.h>
#include <stxxl/bits/namespace.h>

#if __cplusplus >= 201103L
#include <type_traits>
#elif defined(__GNUG__) && (__GNUC__ >= 4)
#include <tr1/type_traits>
#elif STXXL_BOOST_CONFIG
#include <boost/type_traits/remove_const.hpp>
#endif

STXXL_BEGIN_NAMESPACE

namespace compat {

#if __cplusplus >= 201103L
using std::remove_const;
#elif defined(__GNUG__) && (__GNUC__ >= 4)
using std::tr1::remove_const;
#elif STXXL_BOOST_CONFIG
using boost::remove_const;
#else
template <typename Type>
struct remove_const
{
    typedef Type type;
};

template <typename Type>
struct remove_const<Type const>
{
    typedef Type type;
};
#endif

#if defined(__GNUG__) && ((__GNUC__ * 10000 + __GNUC_MINOR__ * 100) < 40300)
// That is a small subset of what GCC 4.3 does:

// Utility for finding the signed versions of unsigned integral types.
template <typename Type>
struct _make_signed
{
    typedef Type type;
};

template <>
struct _make_signed<char>
{
    typedef signed char type;
};

template <>
struct _make_signed<unsigned char>
{
    typedef signed char type;
};

template <>
struct _make_signed<unsigned short>
{
    typedef signed short type;
};

template <>
struct _make_signed<unsigned int>
{
    typedef signed int type;
};

template <>
struct _make_signed<unsigned long>
{
    typedef signed long type;
};

template <>
struct _make_signed<unsigned long long>
{
    typedef signed long long type;
};

template <typename Type>
struct make_signed
{
    typedef typename _make_signed<Type>::type type;
};
#endif

} // namespace compat

STXXL_END_NAMESPACE

#endif // !STXXL_COMPAT_TYPE_TRAITS_HEADER
