/***************************************************************************
 *  include/stxxl/bits/common/seed.h
 *
 *  Part of the STXXL. See http://stxxl.sourceforge.net
 *
 *  Copyright (C) 2007 Andreas Beckmann <beckmann@mpi-inf.mpg.de>
 *
 *  Distributed under the Boost Software License, Version 1.0.
 *  (See accompanying file LICENSE_1_0.txt or copy at
 *  http://www.boost.org/LICENSE_1_0.txt)
 **************************************************************************/

#ifndef STXXL_COMMON_SEED_HEADER
#define STXXL_COMMON_SEED_HEADER

#include <stxxl/bits/namespace.h>

STXXL_BEGIN_NAMESPACE

//! set the global stxxl seed value
void set_seed(unsigned seed);

//! get a seed value for prng initialization, subsequent calls return a
//! sequence of different values
unsigned get_next_seed();

STXXL_END_NAMESPACE

#endif // !STXXL_COMMON_SEED_HEADER
