/***************************************************************************
 *  include/stxxl/bits/io/linuxaio_file.h
 *
 *  Part of the STXXL. See http://stxxl.sourceforge.net
 *
 *  Copyright (C) 2011 Johannes Singler <singler@kit.edu>
 *  Copyright (C) 2014 Timo Bingmann <tb@panthema.net>
 *
 *  Distributed under the Boost Software License, Version 1.0.
 *  (See accompanying file LICENSE_1_0.txt or copy at
 *  http://www.boost.org/LICENSE_1_0.txt)
 **************************************************************************/

#ifndef STXXL_IO_LINUXAIO_FILE_HEADER
#define STXXL_IO_LINUXAIO_FILE_HEADER

#include <stxxl/bits/config.h>

#if STXXL_HAVE_LINUXAIO_FILE

#include <stxxl/bits/io/ufs_file_base.h>
#include <stxxl/bits/io/disk_queued_file.h>
#include <stxxl/bits/io/linuxaio_queue.h>

STXXL_BEGIN_NAMESPACE

class linuxaio_queue;

//! \addtogroup fileimpl
//! \{

//! Implementation of \c file based on the Linux kernel interface for
//! asynchronous I/O
class linuxaio_file : public ufs_file_base, public disk_queued_file
{
    friend class linuxaio_request;

private:
    int desired_queue_length;

public:
    //! Constructs file object
    //! \param filename path of file
    //! \param mode open mode, see \c stxxl::file::open_modes
    //! \param queue_id disk queue identifier
    //! \param allocator_id linked disk_allocator
    //! \param device_id physical device identifier
    //! \param desired_queue_length queue length requested from kernel
    linuxaio_file(
        const std::string& filename, int mode,
        int queue_id = DEFAULT_LINUXAIO_QUEUE,
        int allocator_id = NO_ALLOCATOR,
        unsigned int device_id = DEFAULT_DEVICE_ID,
        int desired_queue_length = 0)
        : file(device_id),
          ufs_file_base(filename, mode),
          disk_queued_file(queue_id, allocator_id),
          desired_queue_length(desired_queue_length)
    { }

    void serve(void* buffer, offset_type offset, size_type bytes,
               request::request_type type);
    request_ptr aread(void* buffer, offset_type pos, size_type bytes,
                      const completion_handler& on_cmpl = completion_handler());
    request_ptr awrite(void* buffer, offset_type pos, size_type bytes,
                       const completion_handler& on_cmpl = completion_handler());
    const char * io_type() const;

    int get_desired_queue_length() const
    {
        return desired_queue_length;
    }
};

//! \}

STXXL_END_NAMESPACE

#endif // #if STXXL_HAVE_LINUXAIO_FILE

#endif // !STXXL_IO_LINUXAIO_FILE_HEADER
// vim: et:ts=4:sw=4
