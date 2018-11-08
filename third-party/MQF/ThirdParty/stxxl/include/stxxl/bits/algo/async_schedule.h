/***************************************************************************
 *  include/stxxl/bits/algo/async_schedule.h
 *
 *  Part of the STXXL. See http://stxxl.sourceforge.net
 *
 *  Copyright (C) 2002 Roman Dementiev <dementiev@mpi-sb.mpg.de>
 *  Copyright (C) 2009 Andreas Beckmann <beckmann@cs.uni-frankfurt.de>
 *
 *  Distributed under the Boost Software License, Version 1.0.
 *  (See accompanying file LICENSE_1_0.txt or copy at
 *  http://www.boost.org/LICENSE_1_0.txt)
 **************************************************************************/

#ifndef STXXL_ALGO_ASYNC_SCHEDULE_HEADER
#define STXXL_ALGO_ASYNC_SCHEDULE_HEADER

// Implements the "prudent prefetching" as described in
// D. Hutchinson, P. Sanders, J. S. Vitter: Duality between prefetching
// and queued writing on parallel disks, 2005
// DOI: 10.1137/S0097539703431573

#include <stxxl/bits/common/types.h>
#include <stxxl/bits/common/simple_vector.h>
#include <stxxl/bits/namespace.h>

STXXL_BEGIN_NAMESPACE

void compute_prefetch_schedule(
    const int_type* first,
    const int_type* last,
    int_type* out_first,
    int_type m,
    int_type D);

inline void compute_prefetch_schedule(
    int_type* first,
    int_type* last,
    int_type* out_first,
    int_type m,
    int_type D)
{
    compute_prefetch_schedule(static_cast<const int_type*>(first), last, out_first, m, D);
}

template <typename RunType>
void compute_prefetch_schedule(
    const RunType& input,
    int_type* out_first,
    int_type m,
    int_type D)
{
    const int_type L = input.size();
    simple_vector<int_type> disks(L);
    for (int_type i = 0; i < L; ++i)
        disks[i] = input[i].bid.storage->get_device_id();
    compute_prefetch_schedule(disks.begin(), disks.end(), out_first, m, D);
}

template <typename BidIteratorType>
void compute_prefetch_schedule(
    BidIteratorType input_begin,
    BidIteratorType input_end,
    int_type* out_first,
    int_type m,
    int_type D)
{
    const int_type L = input_end - input_begin;
    simple_vector<int_type> disks(L);
    int_type i = 0;
    for (BidIteratorType it = input_begin; it != input_end; ++it, ++i)
        disks[i] = it->storage->get_device_id();
    compute_prefetch_schedule(disks.begin(), disks.end(), out_first, m, D);
}

STXXL_END_NAMESPACE

#endif // !STXXL_ALGO_ASYNC_SCHEDULE_HEADER
// vim: et:ts=4:sw=4
