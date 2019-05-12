/***************************************************************************
 *  include/stxxl/bits/io/linuxaio_queue.h
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

#ifndef STXXL_IO_LINUXAIO_QUEUE_HEADER
#define STXXL_IO_LINUXAIO_QUEUE_HEADER

#include <stxxl/bits/io/linuxaio_file.h>

#if STXXL_HAVE_LINUXAIO_FILE

#include <linux/aio_abi.h>
#include <list>

#include <stxxl/bits/io/request_queue_impl_worker.h>
#include <stxxl/bits/common/mutex.h>

STXXL_BEGIN_NAMESPACE

//! \addtogroup reqlayer
//! \{

//! Queue for linuxaio_file(s)
//!
//! Only one queue exists in a program, i.e. it is a singleton.
class linuxaio_queue : public request_queue_impl_worker
{
    friend class linuxaio_request;

    typedef linuxaio_queue self_type;

private:
    //! OS context
    aio_context_t context;

    //! storing linuxaio_request* would drop ownership
    typedef std::list<request_ptr> queue_type;

    // "waiting" request have submitted to this queue, but not yet to the OS,
    // those are "posted"
    mutex waiting_mtx, posted_mtx;
    queue_type waiting_requests;

    //! max number of OS requests
    int max_events;
    //! number of requests in waitings_requests
    semaphore num_waiting_requests, num_free_events, num_posted_requests;

    // two threads, one for posting, one for waiting
    thread_type post_thread, wait_thread;
    state<thread_state> post_thread_state, wait_thread_state;

    // Why do we need two threads, one for posting, and one for waiting?  Is
    // one not enough?
    // 1. User call cannot io_submit directly, since this tends to take
    //    considerable time sometimes
    // 2. A single thread cannot wait for the user program to post requests
    //    and the OS to produce I/O completion events at the same time
    //    (IOCB_CMD_NOOP does not seem to help here either)

    static const priority_op _priority_op = WRITE;

    static void * post_async(void* arg);   // thread start callback
    static void * wait_async(void* arg);   // thread start callback
    void post_requests();
    void handle_events(io_event* events, long num_events, bool canceled);
    void wait_requests();
    void suspend();

    // needed by linuxaio_request
    aio_context_t get_io_context() { return context; }

public:
    //! Construct queue. Requests max number of requests simultaneously
    //! submitted to disk, 0 means as many as possible
    linuxaio_queue(int desired_queue_length = 0);

    void add_request(request_ptr& req);
    bool cancel_request(request_ptr& req);
    void complete_request(request_ptr& req);
    ~linuxaio_queue();
};

//! \}

STXXL_END_NAMESPACE

#endif // #if STXXL_HAVE_LINUXAIO_FILE

#endif // !STXXL_IO_LINUXAIO_QUEUE_HEADER
// vim: et:ts=4:sw=4
