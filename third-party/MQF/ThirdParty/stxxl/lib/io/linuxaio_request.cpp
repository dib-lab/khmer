/***************************************************************************
 *  lib/io/linuxaio_request.cpp
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

#include <stxxl/bits/io/linuxaio_request.h>

#if STXXL_HAVE_LINUXAIO_FILE

#include <stxxl/bits/io/disk_queues.h>
#include <stxxl/bits/verbose.h>
#include <stxxl/bits/common/error_handling.h>

#include <unistd.h>
#include <sys/syscall.h>

STXXL_BEGIN_NAMESPACE

void linuxaio_request::completed(bool posted, bool canceled)
{
    STXXL_VERBOSE_LINUXAIO("linuxaio_request[" << this << "] completed(" <<
                           posted << "," << canceled << ")");

    if (!canceled)
    {
        if (m_type == READ)
            stats::get_instance()->read_finished();
        else
            stats::get_instance()->write_finished();
    }
    else if (posted)
    {
        if (m_type == READ)
            stats::get_instance()->read_canceled(m_bytes);
        else
            stats::get_instance()->write_canceled(m_bytes);
    }
    request_with_state::completed(canceled);
}

void linuxaio_request::fill_control_block()
{
    linuxaio_file* af = dynamic_cast<linuxaio_file*>(m_file);

    memset(&cb, 0, sizeof(cb));
    // indirection, so the I/O system retains a counting_ptr reference
    cb.aio_data = reinterpret_cast<__u64>(new request_ptr(this));
    cb.aio_fildes = af->file_des;
    cb.aio_lio_opcode = (m_type == READ) ? IOCB_CMD_PREAD : IOCB_CMD_PWRITE;
    cb.aio_reqprio = 0;
    cb.aio_buf = static_cast<__u64>((unsigned long)(m_buffer));
    cb.aio_nbytes = m_bytes;
    cb.aio_offset = m_offset;
}

//! Submits an I/O request to the OS
//! \returns false if submission fails
bool linuxaio_request::post()
{
    STXXL_VERBOSE_LINUXAIO("linuxaio_request[" << this << "] post()");

    fill_control_block();
    iocb* cb_pointer = &cb;
    // io_submit might considerable time, so we have to remember the current
    // time before the call.
    double now = timestamp();
    linuxaio_queue* queue = dynamic_cast<linuxaio_queue*>(
        disk_queues::get_instance()->get_queue(m_file->get_queue_id())
        );
    long success = syscall(SYS_io_submit, queue->get_io_context(), 1, &cb_pointer);
    if (success == 1)
    {
        if (m_type == READ)
            stats::get_instance()->read_started(m_bytes, now);
        else
            stats::get_instance()->write_started(m_bytes, now);
    }
    else if (success == -1 && errno != EAGAIN)
        STXXL_THROW_ERRNO(io_error, "linuxaio_request::post"
                          " io_submit()");

    return success == 1;
}

//! Cancel the request
//!
//! Routine is called by user, as part of the request interface.
bool linuxaio_request::cancel()
{
    STXXL_VERBOSE_LINUXAIO("linuxaio_request[" << this << "] cancel()");

    if (!m_file) return false;

    request_ptr req(this);
    linuxaio_queue* queue = dynamic_cast<linuxaio_queue*>(
        disk_queues::get_instance()->get_queue(m_file->get_queue_id())
        );
    return queue->cancel_request(req);
}

//! Cancel already posted request
bool linuxaio_request::cancel_aio()
{
    STXXL_VERBOSE_LINUXAIO("linuxaio_request[" << this << "] cancel_aio()");

    if (!m_file) return false;

    io_event event;
    linuxaio_queue* queue = dynamic_cast<linuxaio_queue*>(
        disk_queues::get_instance()->get_queue(m_file->get_queue_id())
        );
    long result = syscall(SYS_io_cancel, queue->get_io_context(), &cb, &event);
    if (result == 0)    //successfully canceled
        queue->handle_events(&event, 1, true);
    return result == 0;
}

STXXL_END_NAMESPACE

#endif // #if STXXL_HAVE_LINUXAIO_FILE
// vim: et:ts=4:sw=4
