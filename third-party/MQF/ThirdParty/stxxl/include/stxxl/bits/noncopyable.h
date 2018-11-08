/***************************************************************************
 *  include/stxxl/bits/noncopyable.h
 *
 *  Inspired by boost::noncopyable.
 *
 *  Part of the STXXL. See http://stxxl.sourceforge.net
 *
 *  Copyright (C) 2007 Andreas Beckmann <beckmann@mpi-inf.mpg.de>
 *
 *  Distributed under the Boost Software License, Version 1.0.
 *  (See accompanying file LICENSE_1_0.txt or copy at
 *  http://www.boost.org/LICENSE_1_0.txt)
 **************************************************************************/

#ifndef STXXL_NONCOPYABLE_HEADER
#define STXXL_NONCOPYABLE_HEADER

#include <stxxl/bits/config.h>
#include <stxxl/bits/namespace.h>

#if STXXL_BOOST_CONFIG
  #include <boost/noncopyable.hpp>
#endif

STXXL_BEGIN_NAMESPACE

#if STXXL_BOOST_CONFIG

typedef boost::noncopyable noncopyable;

#else

class noncopyable
{
protected:
    noncopyable() { }

private:
    // copying and assignment is not allowed
    noncopyable(const noncopyable&);
    const noncopyable& operator = (const noncopyable&);
};

#endif

STXXL_END_NAMESPACE

#endif // !STXXL_NONCOPYABLE_HEADER
