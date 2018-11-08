/***************************************************************************
 *  lib/algo/async_schedule.cpp
 *
 *  Part of the STXXL. See http://stxxl.sourceforge.net
 *
 *  Copyright (C) 2002, 2009 Roman Dementiev <dementiev@mpi-sb.mpg.de>
 *  Copyright (C) 2009, 2010 Andreas Beckmann <beckmann@cs.uni-frankfurt.de>
 *
 *  Distributed under the Boost Software License, Version 1.0.
 *  (See accompanying file LICENSE_1_0.txt or copy at
 *  http://www.boost.org/LICENSE_1_0.txt)
 **************************************************************************/

// Implements the "prudent prefetching" as described in
// D. Hutchinson, P. Sanders, J. S. Vitter: Duality between prefetching
// and queued writing on parallel disks, 2005
// DOI: 10.1137/S0097539703431573

#include <stxxl/bits/algo/async_schedule.h>
#include <stxxl/bits/common/simple_vector.h>
#include <stxxl/bits/common/types.h>
#include <stxxl/bits/io/file.h>
#include <stxxl/bits/namespace.h>
#include <stxxl/bits/parallel.h>
#include <stxxl/bits/unused.h>
#include <stxxl/bits/verbose.h>

#include <algorithm>
#include <cassert>
#include <functional>
#include <queue>
#include <utility>
#include <vector>

STXXL_BEGIN_NAMESPACE

namespace async_schedule_local {

// only one type of event: WRITE COMPLETED
struct sim_event
{
    int_type timestamp;
    int_type iblock;
    inline sim_event(int_type t, int_type b) : timestamp(t), iblock(b) { }
};

struct sim_event_cmp : public std::binary_function<sim_event, sim_event, bool>
{
    inline bool operator () (const sim_event& a, const sim_event& b) const
    {
        return a.timestamp > b.timestamp;
    }
};

typedef std::pair<int_type, int_type> write_time_pair;
struct write_time_cmp : public std::binary_function<write_time_pair, write_time_pair, bool>
{
    inline bool operator () (const write_time_pair& a, const write_time_pair& b) const
    {
        return a.second > b.second;
    }
};

static inline int_type get_disk(int_type i, const int_type* disks, int_type D)
{
    int_type disk = disks[i];
    if (disk == (int_type)file::DEFAULT_DEVICE_ID)
        disk = D;      // remap to sentinel
    assert(0 <= disk && disk <= D);
    return disk;
}

int_type simulate_async_write(
    const int_type* disks,
    const int_type L,
    const int_type m_init,
    const int_type D,
    std::pair<int_type, int_type>* o_time)
{
    typedef std::priority_queue<sim_event, std::vector<sim_event>, sim_event_cmp> event_queue_type;
    typedef std::queue<int_type> disk_queue_type;
    assert(L >= D);
    simple_vector<disk_queue_type> disk_queues(D + 1); // + sentinel for remapping NO_ALLOCATOR
    event_queue_type event_queue;

    int_type m = m_init;
    int_type i = L - 1;
    int_type oldtime = 0;
    simple_vector<bool> disk_busy(D + 1);

    while (m && (i >= 0))
    {
        int_type disk = get_disk(i, disks, D);
        disk_queues[disk].push(i);
        i--;
        m--;
    }

    for (int_type ii = 0; ii <= D; ii++)
        if (!disk_queues[ii].empty())
        {
            int_type j = disk_queues[ii].front();
            disk_queues[ii].pop();
            event_queue.push(sim_event(1, j));
            //STXXL_MSG("Block "<<j<<" scheduled");
        }

    while (!event_queue.empty())
    {
        sim_event cur = event_queue.top();
        event_queue.pop();
        if (oldtime != cur.timestamp)
        {
            // clear disk_busy
            for (int_type j = 0; j <= D; j++)
                disk_busy[j] = false;

            oldtime = cur.timestamp;
        }

        STXXL_VERBOSE1("Block " << cur.iblock << " put out, time " << cur.timestamp << " disk: " << disks[cur.iblock]);
        o_time[cur.iblock] = std::pair<int_type, int_type>(cur.iblock, cur.timestamp);

        if (i >= 0)
        {
            int_type disk = get_disk(i, disks, D);
            if (disk_busy[disk])
            {
                disk_queues[disk].push(i--);
            }
            else
            {
                if (!disk_queues[disk].empty())
                {
                    STXXL_VERBOSE1("c Block " << disk_queues[disk].front() << " scheduled for time " << cur.timestamp + 1);
                    event_queue.push(sim_event(cur.timestamp + 1, disk_queues[disk].front()));
                    disk_queues[disk].pop();
                }
                else
                {
                    STXXL_VERBOSE1("a Block " << i << " scheduled for time " << cur.timestamp + 1);
                    event_queue.push(sim_event(cur.timestamp + 1, i--));
                }
                disk_busy[disk] = true;
            }
        }

        // add next block to write
        int_type disk = get_disk(cur.iblock, disks, D);
        if (!disk_busy[disk] && !disk_queues[disk].empty())
        {
            STXXL_VERBOSE1("b Block " << disk_queues[disk].front() << " scheduled for time " << cur.timestamp + 1);
            event_queue.push(sim_event(cur.timestamp + 1, disk_queues[disk].front()));
            disk_queues[disk].pop();
            disk_busy[disk] = true;
        }
    }

    assert(i == -1);
    for (int_type j = 0; j <= D; j++)
        assert(disk_queues[j].empty());

    return (oldtime - 1);
}

} // namespace async_schedule_local

void compute_prefetch_schedule(
    const int_type* first,
    const int_type* last,
    int_type* out_first,
    int_type m,
    int_type D)
{
    typedef std::pair<int_type, int_type> pair_type;
    int_type L = last - first;
    if (L <= D)
    {
        for (int_type i = 0; i < L; ++i)
            out_first[i] = i;

        return;
    }
    pair_type* write_order = new pair_type[L];

    int_type w_steps = async_schedule_local::simulate_async_write(first, L, m, D, write_order);

    STXXL_VERBOSE1("Write steps: " << w_steps);

    for (int_type i = 0; i < L; i++)
        STXXL_VERBOSE1(first[i] << " " << write_order[i].first << " " << write_order[i].second);

    std::stable_sort(write_order, write_order + L, async_schedule_local::write_time_cmp() _STXXL_FORCE_SEQUENTIAL);

    for (int_type i = 0; i < L; i++)
    {
        out_first[i] = write_order[i].first;
        //if(out_first[i] != i)
        STXXL_VERBOSE1(i << " " << out_first[i]);
    }

    delete[] write_order;
    STXXL_UNUSED(w_steps);
}

STXXL_END_NAMESPACE

// vim: et:ts=4:sw=4
