/***************************************************************************
 *  lib/common/rand.cpp
 *
 *  Part of the STXXL. See http://stxxl.sourceforge.net
 *
 *  Copyright (C) 2002 Roman Dementiev <dementiev@mpi-sb.mpg.de>
 *  Copyright (C) 2007 Andreas Beckmann <beckmann@mpi-inf.mpg.de>
 *
 *  Distributed under the Boost Software License, Version 1.0.
 *  (See accompanying file LICENSE_1_0.txt or copy at
 *  http://www.boost.org/LICENSE_1_0.txt)
 **************************************************************************/

#include <stxxl/bits/common/rand.h>
#include <stxxl/bits/common/seed.h>
#include <stxxl/bits/namespace.h>

STXXL_BEGIN_NAMESPACE

unsigned ran32State = get_next_seed();

STXXL_END_NAMESPACE
