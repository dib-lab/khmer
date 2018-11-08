/***************************************************************************
 *  lib/io/mem_file.cpp
 *
 *  Part of the STXXL. See http://stxxl.sourceforge.net
 *
 *  Copyright (C) 2008 Andreas Beckmann <beckmann@cs.uni-frankfurt.de>
 *  Copyright (C) 2013-2014 Timo Bingmann <tb@panthema.net>
 *
 *  Distributed under the Boost Software License, Version 1.0.
 *  (See accompanying file LICENSE_1_0.txt or copy at
 *  http://www.boost.org/LICENSE_1_0.txt)
 **************************************************************************/

#include <cstring>
#include <limits>
#include <cassert>

#include <stxxl/bits/io/mem_file.h>
#include <stxxl/bits/io/iostats.h>

STXXL_BEGIN_NAMESPACE

void mem_file::serve(void* buffer, offset_type offset, size_type bytes,
                     request::request_type type)
{
    scoped_mutex_lock lock(m_mutex);

    if (type == request::READ)
    {
        stats::scoped_read_timer read_timer(bytes);
        memcpy(buffer, m_ptr + offset, bytes);
    }
    else
    {
        stats::scoped_write_timer write_timer(bytes);
        memcpy(m_ptr + offset, buffer, bytes);
    }
}

const char* mem_file::io_type() const
{
    return "memory";
}

mem_file::~mem_file()
{
    free(m_ptr);
    m_ptr = NULL;
}

void mem_file::lock()
{
    // nothing to do
}

file::offset_type mem_file::size()
{
    return m_size;
}

void mem_file::set_size(offset_type newsize)
{
    scoped_mutex_lock lock(m_mutex);
    assert(newsize <= std::numeric_limits<offset_type>::max());

    m_ptr = (char*)realloc(m_ptr, (size_t)newsize);
    m_size = newsize;
}

void mem_file::discard(offset_type offset, offset_type size)
{
    scoped_mutex_lock lock(m_mutex);
#ifndef STXXL_MEMFILE_DONT_CLEAR_FREED_MEMORY
    // overwrite the freed region with uninitialized memory
    STXXL_VERBOSE("discard at " << offset << " len " << size);
    void* uninitialized = malloc(STXXL_BLOCK_ALIGN);
    while (size >= STXXL_BLOCK_ALIGN) {
        memcpy(m_ptr + offset, uninitialized, STXXL_BLOCK_ALIGN);
        offset += STXXL_BLOCK_ALIGN;
        size -= STXXL_BLOCK_ALIGN;
    }
    assert(size <= std::numeric_limits<offset_type>::max());
    if (size > 0)
        memcpy(m_ptr + offset, uninitialized, (size_t)size);
    free(uninitialized);
#else
    STXXL_UNUSED(offset);
    STXXL_UNUSED(size);
#endif
}

STXXL_END_NAMESPACE
