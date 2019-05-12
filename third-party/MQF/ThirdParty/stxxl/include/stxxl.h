/***************************************************************************
 *  include/stxxl.h
 *
 *  Part of the STXXL. See http://stxxl.sourceforge.net
 *
 *  Copyright (C) 2002-2007 Roman Dementiev <dementiev@mpi-sb.mpg.de>
 *
 *  Distributed under the Boost Software License, Version 1.0.
 *  (See accompanying file LICENSE_1_0.txt or copy at
 *  http://www.boost.org/LICENSE_1_0.txt)
 **************************************************************************/

#ifndef STXXL_MAIN_HEADER
#define STXXL_MAIN_HEADER

#include <stxxl/bits/common/utils.h>

#include <stxxl/io>

#include <stxxl/mng>

#include <stxxl/vector>
#include <stxxl/stack>
#include <stxxl/priority_queue>
#if ! defined(__GNUG__) || ((__GNUC__ * 10000 + __GNUC_MINOR__ * 100) >= 30400)
// map does not work with g++ 3.3
#include <stxxl/map>
#endif
#include <stxxl/unordered_map>
#include <stxxl/queue>
#include <stxxl/deque>

#include <stxxl/algorithm>

#include <stxxl/stream>

#include <stxxl/random>

#include <stxxl/timer>

#endif // STXXL_MAIN_HEADER
