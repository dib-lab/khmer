/***************************************************************************
 *  include/stxxl/bits/deprecated.h
 *
 *  Part of the STXXL. See http://stxxl.sourceforge.net
 *
 *  Copyright (C) 2008-2009 Andreas Beckmann <beckmann@cs.uni-frankfurt.de>
 *  Copyright (C) 2013 Timo Bingmann <tb@panthema.net>
 *
 *  Distributed under the Boost Software License, Version 1.0.
 *  (See accompanying file LICENSE_1_0.txt or copy at
 *  http://www.boost.org/LICENSE_1_0.txt)
 **************************************************************************/

#ifndef STXXL_DEPRECATED_HEADER
#define STXXL_DEPRECATED_HEADER

#include <stxxl/bits/config.h>

// deprecated functions

#if STXXL_NO_DEPRECATED
// dont issue deprecated warnings for forced instantiation tests -tb
  #define STXXL_DEPRECATED(x) x
#elif STXXL_MSVC
  #define STXXL_DEPRECATED(x) __declspec(deprecated) x
#elif defined(__GNUG__) && ((__GNUC__ * 10000 + __GNUC_MINOR__ * 100) < 30400)
// no __attribute__ ((__deprecated__)) in GCC 3.3
  #define STXXL_DEPRECATED(x) x
#else
  #define STXXL_DEPRECATED(x) x __attribute__ ((__deprecated__))
#endif

#endif // !STXXL_DEPRECATED_HEADER
// vim: et:ts=4:sw=4
