/***************************************************************************
 *  lib/io/fileperblock_file.cpp
 *
 *  Part of the STXXL. See http://stxxl.sourceforge.net
 *
 *  Copyright (C) 2008, 2009 Johannes Singler <singler@ira.uka.de>
 *  Copyright (C) 2008 Andreas Beckmann <beckmann@cs.uni-frankfurt.de>
 *
 *  Distributed under the Boost Software License, Version 1.0.
 *  (See accompanying file LICENSE_1_0.txt or copy at
 *  http://www.boost.org/LICENSE_1_0.txt)
 **************************************************************************/

#include <cassert>
#include <cstdio>
#include <iomanip>
#include <sstream>
#include <string>

#include <stxxl/bits/common/aligned_alloc.h>
#include <stxxl/bits/common/counting_ptr.h>
#include <stxxl/bits/common/error_handling.h>
#include <stxxl/bits/common/exceptions.h>
#include <stxxl/bits/config.h>
#include <stxxl/bits/io/boostfd_file.h>
#include <stxxl/bits/io/completion_handler.h>
#include <stxxl/bits/io/disk_queued_file.h>
#include <stxxl/bits/io/file.h>
#include <stxxl/bits/io/fileperblock_file.h>
#include <stxxl/bits/io/mmap_file.h>
#include <stxxl/bits/io/boostfd_file.h>
#include <stxxl/bits/io/wincall_file.h>
#include <stxxl/bits/io/request.h>
#include <stxxl/bits/io/serving_request.h>
#include <stxxl/bits/io/syscall_file.h>
#include <stxxl/bits/io/wincall_file.h>
#include <stxxl/bits/namespace.h>
#include <stxxl/bits/unused.h>
#include <stxxl/bits/verbose.h>
#include "ufs_platform.h"

STXXL_BEGIN_NAMESPACE

template <class base_file_type>
fileperblock_file<base_file_type>::fileperblock_file(
    const std::string& filename_prefix,
    int mode,
    int queue_id,
    int allocator_id,
    unsigned int device_id)
    : file(device_id),
      disk_queued_file(queue_id, allocator_id),
      filename_prefix(filename_prefix),
      mode(mode),
      current_size(0),
      lock_file_created(false),
      lock_file(filename_prefix + "_fpb_lock", mode, queue_id)
{ }

template <class base_file_type>
fileperblock_file<base_file_type>::~fileperblock_file()
{
    if (lock_file_created)
    {
        if (::remove((filename_prefix + "_fpb_lock").c_str()) != 0)
            STXXL_ERRMSG("remove() error on path=" << filename_prefix << "_fpb_lock error=" << strerror(errno));
    }
}

template <class base_file_type>
std::string fileperblock_file<base_file_type>::filename_for_block(offset_type offset)
{
    std::ostringstream name;
    //enough for 1 billion blocks
    name << filename_prefix << "_fpb_" << std::setw(20) << std::setfill('0') << offset;
    return name.str();
}

template <class base_file_type>
void fileperblock_file<base_file_type>::serve(void* buffer, offset_type offset,
                                              size_type bytes, request::request_type type)
{
    base_file_type base_file(filename_for_block(offset), mode, get_queue_id());
    base_file.set_size(bytes);
    base_file.serve(buffer, 0, bytes, type);
}

template <class base_file_type>
void fileperblock_file<base_file_type>::lock()
{
    if (!lock_file_created)
    {
        //create lock file and fill it with one page, an empty file cannot be locked
        const int page_size = STXXL_BLOCK_ALIGN;
        void* one_page = aligned_alloc<STXXL_BLOCK_ALIGN>(page_size);
#if STXXL_WITH_VALGRIND
        memset(one_page, 0, page_size);
#endif
        lock_file.set_size(page_size);
        request_ptr r = lock_file.awrite(one_page, 0, page_size);
        r->wait();
        aligned_dealloc<STXXL_BLOCK_ALIGN>(one_page);
        lock_file_created = true;
    }
    lock_file.lock();
}

template <class base_file_type>
void fileperblock_file<base_file_type>::discard(offset_type offset, offset_type length)
{
    STXXL_UNUSED(length);
#ifdef STXXL_FILEPERBLOCK_NO_DELETE
    if (::truncate(filename_for_block(offset).c_str(), 0) != 0)
        STXXL_ERRMSG("truncate() error on path=" << filename_for_block(offset) << " error=" << strerror(errno));
#else
    if (::remove(filename_for_block(offset).c_str()) != 0)
        STXXL_ERRMSG("remove() error on path=" << filename_for_block(offset) << " error=" << strerror(errno));
#endif

    STXXL_VERBOSE2("discard " << offset << " + " << length);
}

template <class base_file_type>
void fileperblock_file<base_file_type>::export_files(offset_type offset, offset_type length, std::string filename)
{
    std::string original(filename_for_block(offset));
    filename.insert(0, original.substr(0, original.find_last_of("/") + 1));
    if (::remove(filename.c_str()) != 0)
        STXXL_ERRMSG("remove() error on path=" << filename << " error=" << strerror(errno));

    if (::rename(original.c_str(), filename.c_str()) != 0)
        STXXL_ERRMSG("rename() error on path=" << filename << " to=" << original << " error=" << strerror(errno));

#if !STXXL_WINDOWS
    //TODO: implement on Windows
    if (::truncate(filename.c_str(), length) != 0) {
        STXXL_THROW_ERRNO(io_error, "Error doing truncate()");
    }
#else
    STXXL_UNUSED(length);
#endif
}

template <class base_file_type>
const char* fileperblock_file<base_file_type>::io_type() const
{
    return "fileperblock";
}

////////////////////////////////////////////////////////////////////////////

template class fileperblock_file<syscall_file>;

#if STXXL_HAVE_MMAP_FILE
template class fileperblock_file<mmap_file>;
#endif

#if STXXL_HAVE_WINCALL_FILE
template class fileperblock_file<wincall_file>;
#endif

#if STXXL_HAVE_BOOSTFD_FILE
template class fileperblock_file<boostfd_file>;
#endif

STXXL_END_NAMESPACE
