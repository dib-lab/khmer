/***************************************************************************
 *  include/stxxl/bits/io/mem_file.h
 *
 *  Part of the STXXL. See http://stxxl.sourceforge.net
 *
 *  Copyright (C) 2008 Andreas Beckmann <beckmann@cs.uni-frankfurt.de>
 *  Copyright (C) 2009 Johannes Singler <singler@ira.uka.de>
 *  Copyright (C) 2014 Timo Bingmann <tb@panthema.net>
 *
 *  Distributed under the Boost Software License, Version 1.0.
 *  (See accompanying file LICENSE_1_0.txt or copy at
 *  http://www.boost.org/LICENSE_1_0.txt)
 **************************************************************************/

#ifndef STXXL_IO_MEM_FILE_HEADER
#define STXXL_IO_MEM_FILE_HEADER

#include <stxxl/bits/io/disk_queued_file.h>
#include <stxxl/bits/io/request.h>

STXXL_BEGIN_NAMESPACE

//! \addtogroup fileimpl
//! \{

//! Implementation of file based on new[] and memcpy.
class mem_file : public disk_queued_file
{
    //! pointer to memory area of "file"
    char* m_ptr;

    //! size of memory area
    offset_type m_size;

    //! sequentialize function calls
    mutex m_mutex;

public:
    //! constructs file object.
    mem_file(
        int queue_id = DEFAULT_QUEUE,
        int allocator_id = NO_ALLOCATOR,
        unsigned int device_id = DEFAULT_DEVICE_ID)
        : file(device_id),
          disk_queued_file(queue_id, allocator_id),
          m_ptr(NULL), m_size(0)
    { }
    void serve(void* buffer, offset_type offset, size_type bytes,
               request::request_type type);
    ~mem_file();
    offset_type size();
    void set_size(offset_type newsize);
    void lock();
    void discard(offset_type offset, offset_type size);
    const char * io_type() const;
};

//! \}

STXXL_END_NAMESPACE

#endif // !STXXL_IO_MEM_FILE_HEADER
