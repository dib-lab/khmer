/***************************************************************************
 *  lib/io/wbtl_file.cpp
 *
 *  a write-buffered-translation-layer pseudo file
 *
 *  Part of the STXXL. See http://stxxl.sourceforge.net
 *
 *  Copyright (C) 2008-2009 Andreas Beckmann <beckmann@cs.uni-frankfurt.de>
 *  Copyright (C) 2009 Johannes Singler <singler@ira.uka.de>
 *
 *  Distributed under the Boost Software License, Version 1.0.
 *  (See accompanying file LICENSE_1_0.txt or copy at
 *  http://www.boost.org/LICENSE_1_0.txt)
 **************************************************************************/

#include <stxxl/bits/io/wbtl_file.h>
#include <stxxl/bits/common/error_handling.h>

#if STXXL_HAVE_WBTL_FILE

#include <algorithm>
#include <iomanip>
#include <stxxl/bits/io/io.h>
#include <stxxl/bits/parallel.h>
#include <stxxl/aligned_alloc>

#ifndef STXXL_VERBOSE_WBTL
#define STXXL_VERBOSE_WBTL STXXL_VERBOSE2
#endif

STXXL_BEGIN_NAMESPACE

wbtl_file::wbtl_file(
    file* backend_file,
    size_type write_buffer_size,
    int write_buffers,
    int queue_id, int allocator_id)
    : disk_queued_file(queue_id, allocator_id), storage(backend_file),
      sz(0), write_block_size(write_buffer_size),
      free_bytes(0), curbuf(1), curpos(write_block_size)
{
    STXXL_UNUSED(write_buffers);
    assert(write_buffers == 2); // currently hardcoded
    write_buffer[0] = static_cast<char*>(stxxl::aligned_alloc<STXXL_BLOCK_ALIGN>(write_block_size));
    write_buffer[1] = static_cast<char*>(stxxl::aligned_alloc<STXXL_BLOCK_ALIGN>(write_block_size));
    buffer_address[0] = offset_type(-1);
    buffer_address[1] = offset_type(-1);
}

wbtl_file::~wbtl_file()
{
    stxxl::aligned_dealloc<STXXL_BLOCK_ALIGN>(write_buffer[1]);
    stxxl::aligned_dealloc<STXXL_BLOCK_ALIGN>(write_buffer[0]);
    delete storage;
    storage = 0;
}

void wbtl_file::serve(void* buffer, offset_type offset, size_type bytes,
                      request::request_type type)
{
    if (type == request::READ)
    {
        //stats::scoped_read_timer read_timer(size());
        sread(buffer, offset, bytes);
    }
    else
    {
        //stats::scoped_write_timer write_timer(size());
        swrite(buffer, offset, bytes);
    }
}

void wbtl_file::lock()
{
    storage->lock();
}

wbtl_file::offset_type wbtl_file::size()
{
    return sz;
}

void wbtl_file::set_size(offset_type newsize)
{
    scoped_mutex_lock mapping_lock(mapping_mutex);
    assert(sz <= newsize); // may not shrink
    if (sz < newsize) {
        _add_free_region(sz, newsize - sz);
        storage->set_size(newsize);
        sz = newsize;
        assert(sz == storage->size());
    }
}

#define FMT_A_S(_addr_, _size_) "0x" << std::hex << std::setfill('0') << std::setw(8) << (_addr_) << "/0x" << std::setw(8) << (_size_)
        #define FMT_A_C(_addr_, _size_) "0x" << std::setw(8) << (_addr_) << "(" << std::dec << (_size_) << ")"
        #define FMT_A(_addr_) "0x" << std::setw(8) << (_addr_)

// logical address
void wbtl_file::discard(offset_type offset, offset_type size)
{
    scoped_mutex_lock mapping_lock(mapping_mutex);
    sortseq::iterator physical = address_mapping.find(offset);
    STXXL_VERBOSE_WBTL("wbtl:discard l" << FMT_A_S(offset, size) << " @    p" << FMT_A(physical != address_mapping.end() ? physical->second : 0xffffffff));
    if (physical == address_mapping.end()) {
        // could be OK if the block was never written ...
        //STXXL_ERRMSG("discard: mapping not found: " << FMT_A_S(offset, size) << " ==> " << "???");
    }
    else {
        offset_type physical_offset = physical->second;
        address_mapping.erase(physical);
        _add_free_region(physical_offset, size);
        place_map::iterator reverse = reverse_mapping.find(physical_offset);
        if (reverse == reverse_mapping.end()) {
            STXXL_ERRMSG("discard: reverse mapping not found: " << FMT_A_S(offset, size) << " ==> " << "???");
        }
        else {
            assert(offset == (reverse->second).first);
            reverse_mapping.erase(reverse);
        }
        storage->discard(physical_offset, size);
    }
}

// physical address
void wbtl_file::_add_free_region(offset_type offset, offset_type size)
{
    // mapping_lock has to be aquired by caller
    STXXL_VERBOSE_WBTL("wbtl:addfre  p" << FMT_A_S(offset, size) << " F <= f" << FMT_A_C(free_bytes, free_space.size()));
    offset_type region_pos = offset;
    offset_type region_size = size;
    if (!free_space.empty())
    {
        sortseq::iterator succ = free_space.upper_bound(region_pos);
        sortseq::iterator pred = succ;
        pred--;
        check_corruption(region_pos, region_size, pred, succ);
        if (succ == free_space.end())
        {
            if (pred == free_space.end())
            {
                //dump();
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
                bool succ_is_not_the_first = (succ != free_space.begin());
                if ((*succ).first == region_pos + region_size)
                {
                    // coalesce with successor
                    region_size += (*succ).second;
                    free_space.erase(succ);
                }
                if (succ_is_not_the_first)
                {
                    if (pred == free_space.end())
                    {
                        //dump();
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
    free_bytes += size;
    STXXL_VERBOSE_WBTL("wbtl:free    p" << FMT_A_S(region_pos, region_size) << " F => f" << FMT_A_C(free_bytes, free_space.size()));
}

void wbtl_file::sread(void* buffer, offset_type offset, size_type bytes)
{
    scoped_mutex_lock buffer_lock(buffer_mutex);
    int cached = -1;
    offset_type physical_offset;
    // map logical to physical address
    {
        scoped_mutex_lock mapping_lock(mapping_mutex);
        sortseq::iterator physical = address_mapping.find(offset);
        if (physical == address_mapping.end()) {
            STXXL_ERRMSG("wbtl_read: mapping not found: " << FMT_A_S(offset, bytes) << " ==> " << "???");
            //STXXL_THROW_ERRNO(io_error, "wbtl_read of unmapped memory");
            physical_offset = 0xffffffff;
        }
        else {
            physical_offset = physical->second;
        }
    }

    if (buffer_address[curbuf] <= physical_offset &&
        physical_offset < buffer_address[curbuf] + write_block_size)
    {
        // block is in current write buffer
        assert(physical_offset + bytes <= buffer_address[curbuf] + write_block_size);
        memcpy(buffer, write_buffer[curbuf] + (physical_offset - buffer_address[curbuf]), bytes);
        stats::get_instance()->read_cached(bytes);
        cached = curbuf;
    }
    else if (buffer_address[1 - curbuf] <= physical_offset &&
             physical_offset < buffer_address[1 - curbuf] + write_block_size)
    {
        // block is in previous write buffer
        assert(physical_offset + bytes <= buffer_address[1 - curbuf] + write_block_size);
        memcpy(buffer, write_buffer[1 - curbuf] + (physical_offset - buffer_address[1 - curbuf]), bytes);
        stats::get_instance()->read_cached(bytes);
        cached = curbuf;
    }
    else if (physical_offset == 0xffffffff) {
        // block was deleted or never written before
        char* uninitialized = (char*)malloc(sizeof(char));
        memset(buffer, *uninitialized, bytes);
        free(uninitialized);
    }
    else
    {
        // block is not cached
        request_ptr req = storage->aread(buffer, physical_offset, bytes);
        req->wait(false);
    }
    STXXL_VERBOSE_WBTL("wbtl:sread   l" << FMT_A_S(offset, bytes) << " @    p" << FMT_A(physical_offset) << " " << std::dec << cached);
    STXXL_UNUSED(cached);
}

void wbtl_file::swrite(void* buffer, offset_type offset, size_type bytes)
{
    scoped_mutex_lock buffer_lock(buffer_mutex);
    // is the block already mapped?
    {
        scoped_mutex_lock mapping_lock(mapping_mutex);
        sortseq::iterator physical = address_mapping.find(offset);
        STXXL_VERBOSE_WBTL("wbtl:swrite  l" << FMT_A_S(offset, bytes) << " @ <= p" <<
                           FMT_A_C(physical != address_mapping.end() ? physical->second : 0xffffffff, address_mapping.size()));
        if (physical != address_mapping.end()) {
            mapping_lock.unlock();
            // FIXME: special case if we can replace it in the current writing block
            discard(offset, bytes);
        }
    }

    if (bytes > write_block_size - curpos)
    {
        // not enough space in the current write buffer

        if (buffer_address[curbuf] != offset_type(-1)) {
            STXXL_VERBOSE_WBTL("wbtl:w2disk  p" << FMT_A_S(buffer_address[curbuf], write_block_size));

            // mark remaining part as free
            if (curpos < write_block_size)
                _add_free_region(buffer_address[curbuf] + curpos, write_block_size - curpos);

            if (backend_request.get()) {
                backend_request->wait(false);
            }

            backend_request = storage->awrite(write_buffer[curbuf], buffer_address[curbuf], write_block_size);
        }

        curbuf = 1 - curbuf;

        buffer_address[curbuf] = get_next_write_block();
        curpos = 0;
    }
    assert(bytes <= write_block_size - curpos);

    // write block into buffer
    memcpy(write_buffer[curbuf] + curpos, buffer, bytes);
    stats::get_instance()->write_cached(bytes);

    scoped_mutex_lock mapping_lock(mapping_mutex);
    address_mapping[offset] = buffer_address[curbuf] + curpos;
    reverse_mapping[buffer_address[curbuf] + curpos] = place(offset, bytes);
    STXXL_VERBOSE_WBTL("wbtl:swrite  l" << FMT_A_S(offset, bytes) << " @ => p" << FMT_A_C(buffer_address[curbuf] + curpos, address_mapping.size()));
    curpos += bytes;
}

wbtl_file::offset_type wbtl_file::get_next_write_block()
{
    // mapping_lock has to be aquired by caller
    sortseq::iterator space =
        std::find_if(free_space.begin(), free_space.end(),
                     bind2nd(FirstFit(), write_block_size) _STXXL_FORCE_SEQUENTIAL);

    if (space != free_space.end())
    {
        offset_type region_pos = (*space).first;
        offset_type region_size = (*space).second;
        free_space.erase(space);
        if (region_size > write_block_size)
            free_space[region_pos + write_block_size] = region_size - write_block_size;

        free_bytes -= write_block_size;

        STXXL_VERBOSE_WBTL("wbtl:nextwb  p" << FMT_A_S(region_pos, write_block_size) << " F    f" << FMT_A_C(free_bytes, free_space.size()));
        return region_pos;
    }

    STXXL_THROW_ERRNO(io_error, "OutOfSpace, probably fragmented");
}

void wbtl_file::check_corruption(offset_type region_pos, offset_type region_size,
                                 sortseq::iterator pred, sortseq::iterator succ)
{
    if (pred != free_space.end())
    {
        if (pred->first <= region_pos && pred->first + pred->second > region_pos)
        {
            STXXL_THROW(bad_ext_alloc, "Error: double deallocation of external memory " <<
                        "System info: P " << pred->first << " " << pred->second << " " << region_pos);
        }
    }
    if (succ != free_space.end())
    {
        if (region_pos <= succ->first && region_pos + region_size > succ->first)
        {
            STXXL_THROW(bad_ext_alloc, "Error: double deallocation of external memory "
                        << "System info: S " << region_pos << " " << region_size << " " << succ->first);
        }
    }
}

const char* wbtl_file::io_type() const
{
    return "wbtl";
}

STXXL_END_NAMESPACE

#endif  // #if STXXL_HAVE_WBTL_FILE
// vim: et:ts=4:sw=4
