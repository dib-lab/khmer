/***************************************************************************
 *  include/stxxl/bits/algo/sort_base.h
 *
 *  Part of the STXXL. See http://stxxl.sourceforge.net
 *
 *  Copyright (C) 2009 Andreas Beckmann <beckmann@cs.uni-frankfurt.de>
 *
 *  Distributed under the Boost Software License, Version 1.0.
 *  (See accompanying file LICENSE_1_0.txt or copy at
 *  http://www.boost.org/LICENSE_1_0.txt)
 **************************************************************************/

#ifndef STXXL_ALGO_SORT_BASE_HEADER
#define STXXL_ALGO_SORT_BASE_HEADER

#include <cmath>
#include <stxxl/bits/common/types.h>

#ifndef STXXL_NO_WARN_RECURSIVE_SORT
#define STXXL_WARNMSG_RECURSIVE_SORT STXXL_ERRMSG
#else
#define STXXL_WARNMSG_RECURSIVE_SORT STXXL_VERBOSE
#endif

#ifndef STXXL_SORT_OPTIMAL_PREFETCHING
#define STXXL_SORT_OPTIMAL_PREFETCHING 1
#endif

#ifndef STXXL_CHECK_ORDER_IN_SORTS
#define STXXL_CHECK_ORDER_IN_SORTS 0
#endif

#ifndef STXXL_L2_SIZE
#define STXXL_L2_SIZE  (512 * 1024)
#endif

STXXL_BEGIN_NAMESPACE

// Optimal merging: merge r = pow(nruns,1/ceil(log(nruns)/log(m))) runs at once
inline unsigned_type optimal_merge_factor(unsigned_type num_runs, unsigned_type max_concurrent_runs)
{
    return unsigned_type(ceil(pow(double(num_runs), 1. / ceil(log(double(num_runs)) / log(double(max_concurrent_runs))))));
}

STXXL_END_NAMESPACE

#endif // !STXXL_ALGO_SORT_BASE_HEADER
// vim: et:ts=4:sw=4
