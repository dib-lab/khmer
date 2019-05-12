/***************************************************************************
 *  lib/mng/disk_allocator.cpp
 *
 *  Part of the STXXL. See http://stxxl.sourceforge.net
 *
 *  Copyright (C) 2002 Roman Dementiev <dementiev@mpi-sb.mpg.de>
 *  Copyright (C) 2008, 2010 Andreas Beckmann <beckmann@cs.uni-frankfurt.de>
 *
 *  Distributed under the Boost Software License, Version 1.0.
 *  (See accompanying file LICENSE_1_0.txt or copy at
 *  http://www.boost.org/LICENSE_1_0.txt)
 **************************************************************************/

#include <stxxl/bits/common/error_handling.h>
#include <stxxl/bits/common/exceptions.h>
#include <stxxl/bits/common/types.h>
#include <stxxl/bits/mng/disk_allocator.h>
#include <stxxl/bits/namespace.h>
#include <stxxl/bits/verbose.h>

#include <cassert>
#include <map>
#include <ostream>
#include <utility>

STXXL_BEGIN_NAMESPACE

void disk_allocator::dump() const
{
    int64 total = 0;
    sortseq::const_iterator cur = free_space.begin();
    STXXL_ERRMSG("Free regions dump:");
    for ( ; cur != free_space.end(); ++cur)
    {
        STXXL_ERRMSG("Free chunk: begin: " << (cur->first) << " size: " << (cur->second));
        total += cur->second;
    }
    STXXL_ERRMSG("Total bytes: " << total);
}

void disk_allocator::deallocation_error(
    stxxl::int64 block_pos, stxxl::int64 block_size,
    const sortseq::iterator& pred, const sortseq::iterator& succ) const
{
    STXXL_ERRMSG("Error deallocating block at " << block_pos << " size " << block_size);
    STXXL_ERRMSG(((pred == succ) ? "pred==succ" : "pred!=succ"));
    if (pred == free_space.end()) {
        STXXL_ERRMSG("pred==free_space.end()");
    }
    else {
        if (pred == free_space.begin())
            STXXL_ERRMSG("pred==free_space.begin()");
        STXXL_ERRMSG("pred: begin=" << pred->first << " size=" << pred->second);
    }
    if (succ == free_space.end()) {
        STXXL_ERRMSG("succ==free_space.end()");
    }
    else {
        if (succ == free_space.begin())
            STXXL_ERRMSG("succ==free_space.begin()");
        STXXL_ERRMSG("succ: begin=" << succ->first << " size=" << succ->second);
    }
    dump();
}

void disk_allocator::add_free_region(stxxl::int64 block_pos, stxxl::int64 block_size)
{
    //assert(block_size);
    //dump();
    STXXL_VERBOSE2("Deallocating a block with size: " << block_size << " position: " << block_pos);
    stxxl::int64 region_pos = block_pos;
    stxxl::int64 region_size = block_size;
    if (!free_space.empty())
    {
        sortseq::iterator succ = free_space.upper_bound(region_pos);
        sortseq::iterator pred = succ;
        if (pred != free_space.begin())
            pred--;
        if (pred != free_space.end())
        {
            if (pred->first <= region_pos && pred->first + pred->second > region_pos)
            {
                STXXL_THROW2(bad_ext_alloc, "disk_allocator::check_corruption", "Error: double deallocation of external memory, trying to deallocate region " << region_pos << " + " << region_size << "  in empty space [" << pred->first << " + " << pred->second << "]");
            }
        }
        if (succ != free_space.end())
        {
            if (region_pos <= succ->first && region_pos + region_size > succ->first)
            {
                STXXL_THROW2(bad_ext_alloc, "disk_allocator::check_corruption", "Error: double deallocation of external memory, trying to deallocate region " << region_pos << " + " << region_size << "  which overlaps empty space [" << succ->first << " + " << succ->second << "]");
            }
        }
        if (succ == free_space.end())
        {
            if (pred == free_space.end())
            {
                deallocation_error(block_pos, block_size, pred, succ);
                assert(pred != free_space.end());
            }
            if ((*pred).first + (*pred).second == region_pos)
            {
                // coalesce with predecessor
                region_size += (*pred).second;
                region_pos = (*pred).first;
                free_space.erase(pred);
            }
        }
        else
        {
            if (free_space.size() > 1)
            {
#if 0
                if (pred == succ)
                {
                    deallocation_error(block_pos, block_size, pred, succ);
                    assert(pred != succ);
                }
#endif
                bool succ_is_not_the_first = (succ != free_space.begin());
                if ((*succ).first == region_pos + region_size)
                {
                    // coalesce with successor
                    region_size += (*succ).second;
                    free_space.erase(succ);
                    //-tb: set succ to pred afterwards due to iterator invalidation
                    succ = pred;
                }
                if (succ_is_not_the_first)
                {
                    if (pred == free_space.end())
                    {
                        deallocation_error(block_pos, block_size, pred, succ);
                        assert(pred != free_space.end());
                    }
                    if ((*pred).first + (*pred).second == region_pos)
                    {
                        // coalesce with predecessor
                        region_size += (*pred).second;
                        region_pos = (*pred).first;
                        free_space.erase(pred);
                    }
                }
            }
            else
            {
                if ((*succ).first == region_pos + region_size)
                {
                    // coalesce with successor
                    region_size += (*succ).second;
                    free_space.erase(succ);
                }
            }
        }
    }

    free_space[region_pos] = region_size;
    free_bytes += block_size;

    //dump();
}

STXXL_END_NAMESPACE
// vim: et:ts=4:sw=4
