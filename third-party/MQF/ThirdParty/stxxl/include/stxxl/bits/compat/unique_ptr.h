/***************************************************************************
 *  include/stxxl/bits/compat/unique_ptr.h
 *
 *  compatibility interface to unique_ptr (C++0x), previously auto_ptr
 *
 *  Part of the STXXL. See http://stxxl.sourceforge.net
 *
 *  Copyright (C) 2008 Andreas Beckmann <beckmann@cs.uni-frankfurt.de>
 *
 *  Distributed under the Boost Software License, Version 1.0.
 *  (See accompanying file LICENSE_1_0.txt or copy at
 *  http://www.boost.org/LICENSE_1_0.txt)
 **************************************************************************/

#ifndef STXXL_COMPAT_UNIQUE_PTR_HEADER
#define STXXL_COMPAT_UNIQUE_PTR_HEADER

#include <memory>
#include <stxxl/bits/namespace.h>

STXXL_BEGIN_NAMESPACE

template <class Type>
struct compat_unique_ptr {
#if __cplusplus >= 201103L && ((__GNUC__ * 10000 + __GNUC_MINOR__ * 100) >= 40400)
    typedef std::unique_ptr<Type> result;
#else
    // auto_ptr is inherently broken and is deprecated by unique_ptr in c++0x
    typedef std::auto_ptr<Type> result;
#endif
};

STXXL_END_NAMESPACE

#if defined(__GNUG__) && ((__GNUC__ * 10000 + __GNUC_MINOR__ * 100) == 30400)

namespace workaround_gcc_3_4 {

// std::swap in gcc 3.4 is broken, tmp is declared const there
template <typename Type>
inline void
swap(Type& a, Type& b)
{
    // concept requirements
    __glibcxx_function_requires(_SGIAssignableConcept<Type>)

    Type tmp = a;
    a = b;
    b = tmp;
}

} // namespace workaround_gcc_3_4

namespace std {

// overload broken std::swap<auto_ptr> to call a working swap()
template <typename Type>
inline void swap(std::auto_ptr<Type>& a, std::auto_ptr<Type>& b)
{
    workaround_gcc_3_4::swap(a, b);
}

} // namespace std

#endif

#endif // !STXXL_COMPAT_UNIQUE_PTR_HEADER
