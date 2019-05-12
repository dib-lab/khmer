/***************************************************************************
 *  lib/io/linuxaio_queue.cpp
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

#include <stxxl/bits/io/linuxaio_queue.h>

#if STXXL_HAVE_LINUXAIO_FILE

#include <unistd.h>
#include <sys/syscall.h>

#include <stxxl/bits/verbose.h>
#include <stxxl/bits/mng/block_manager.h>
#include <stxxl/bits/common/error_handling.h>
#include <stxxl/bits/io/linuxaio_request.h>
#include <stxxl/bits/io/linuxaio_queue.h>

#include <algorithm>

#ifndef STXXL_CHECK_FOR_PENDING_REQUESTS_ON_SUBMISSION
#define STXXL_CHECK_FOR_PENDING_REQUESTS_ON_SUBMISSION 1
#endif

STXXL_BEGIN_NAMESPACE

linuxaio_queue::linuxaio_queue(int desired_queue_length)
    : num_waiting_requests(0), num_free_events(0), num_posted_requests(0),
      post_thread_state(NOT_RUNNING), wait_thread_state(NOT_RUNNING)
{
    if (desired_queue_length == 0) {
        // default value, 64 entries per queue (i.e. usually per disk) should
        // be enough
        max_events = 64;
    }
    else
        max_events = desired_queue_length;

    // negotiate maximum number of simultaneous events with the OS
    context = 0;
    long result;
    while ((result = syscall(SYS_io_setup, max_events, &context)) == -1 &&
           errno == EAGAIN && max_events > 1)
    {
        max_events <<= 1;               // try with half as many events
    }
    if (result != 0) {
        STXXL_THROW_ERRNO(io_error, "linuxaio_queue::linuxaio_queue"
                          " io_setup() nr_events=" << max_events);
    }

    for (int e = 0; e < max_events; ++e)
        num_free_events++;  // cannot set semaphore to value directly

    STXXL_MSG("Set up an linuxaio queue with " << max_events << " entries.");

    start_thread(post_async, static_cast<void*>(this), post_thread, post_thread_state);
    start_thread(wait_async, static_cast<void*>(this), wait_thread, wait_thread_state);
}

linuxaio_queue::~linuxaio_queue()
{
    stop_thread(post_thread, post_thread_state, num_waiting_requests);
    stop_thread(wait_thread, wait_thread_state, num_posted_requests);
    syscall(SYS_io_destroy, context);
}

void linuxaio_queue::add_request(request_ptr& req)
{
    if (req.empty())
        STXXL_THROW_INVALID_ARGUMENT("Empty request submitted to disk_queue.");
    if (post_thread_state() != RUNNING)
        STXXL_ERRMSG("Request submitted to stopped queue.");
    if (!dynamic_cast<linuxaio_request*>(req.get()))
        STXXL_ERRMSG("Non-LinuxAIO request submitted to LinuxAIO queue.");

    scoped_mutex_lock lock(waiting_mtx);

    waiting_requests.push_back(req);
    num_waiting_requests++;
}

bool linuxaio_queue::cancel_request(request_ptr& req)
{
    if (req.empty())
        STXXL_THROW_INVALID_ARGUMENT("Empty request canceled disk_queue.");
    if (post_thread_state() != RUNNING)
        STXXL_ERRMSG("Request canceled in stopped queue.");
    if (!dynamic_cast<linuxaio_request*>(req.get()))
        STXXL_ERRMSG("Non-LinuxAIO request submitted to LinuxAIO queue.");

    queue_type::iterator pos;
    {
        scoped_mutex_lock lock(waiting_mtx);

        pos = std::find(waiting_requests.begin(), waiting_requests.end(),
                        req _STXXL_FORCE_SEQUENTIAL);
        if (pos != waiting_requests.end())
        {
            waiting_requests.erase(pos);

            // polymorphic_downcast to linuxaio_request,
            // request is canceled, but was not yet posted.
            dynamic_cast<linuxaio_request*>(req.get())->completed(false, true);

            num_waiting_requests--; // will never block
            return true;
        }
    }

    scoped_mutex_lock lock(posted_mtx);

    return false;
}

// internal routines, run by the posting thread
void linuxaio_queue::post_requests()
{
    request_ptr req;
    io_event* events = new io_event[max_events];

    for ( ; ; ) // as long as thread is running
    {
        // might block until next request or message comes in
        int num_currently_waiting_requests = num_waiting_requests--;

        // terminate if termination has been requested
        if (post_thread_state() == TERMINATING && num_currently_waiting_requests == 0)
            break;

        scoped_mutex_lock lock(waiting_mtx);
        if (!waiting_requests.empty())
        {
            req = waiting_requests.front();
            waiting_requests.pop_front();
            lock.unlock();

            num_free_events--; // might block because too many requests are posted

            // polymorphic_downcast
            while (!dynamic_cast<linuxaio_request*>(req.get())->post())
            {
                // post failed, so first handle events to make queues (more)
                // empty, then try again.

                // wait for at least one event to complete, no time limit
                long num_events = syscall(SYS_io_getevents, context, 1, max_events, events, NULL);
                if (num_events < 0) {
                    STXXL_THROW_ERRNO(io_error, "linuxaio_queue::post_requests"
                                      " io_getevents() nr_events=" << num_events);
                }

                handle_events(events, num_events, false);
            }

            // request is finally posted
            num_posted_requests++;
        }
        else
        {
            lock.unlock();

            // num_waiting_requests-- was premature, compensate for that
            num_waiting_requests++;
        }
    }

    delete[] events;
}

void linuxaio_queue::handle_events(io_event* events, long num_events, bool canceled)
{
    for (int e = 0; e < num_events; ++e)
    {
        // unsigned_type is as long as a pointer, and like this, we avoid an icpc warning
        request_ptr* r = reinterpret_cast<request_ptr*>(static_cast<unsigned_type>(events[e].data));
        r->get()->completed(canceled);
        delete r;              // release auto_ptr reference
        num_free_events++;
        num_posted_requests--; // will never block
    }
}

// internal routines, run by the waiting thread
void linuxaio_queue::wait_requests()
{
    request_ptr req;
    io_event* events = new io_event[max_events];

    for ( ; ; ) // as long as thread is running
    {
        // might block until next request is posted or message comes in
        int num_currently_posted_requests = num_posted_requests--;

        // terminate if termination has been requested
        if (wait_thread_state() == TERMINATING && num_currently_posted_requests == 0)
            break;

        // wait for at least one of them to finish
        long num_events;
        while (1) {
            num_events = syscall(SYS_io_getevents, context, 1, max_events, events, NULL);
            if (num_events < 0) {
                if (errno == EINTR) {
                    // io_getevents may return prematurely in case a signal is received
                    continue;
                }

                STXXL_THROW_ERRNO(io_error, "linuxaio_queue::wait_requests"
                                  " io_getevents() nr_events=" << max_events);
            }
            break;
        }

        num_posted_requests++; // compensate for the one eaten prematurely above

        handle_events(events, num_events, false);
    }

    delete[] events;
}

void* linuxaio_queue::post_async(void* arg)
{
    (static_cast<linuxaio_queue*>(arg))->post_requests();

    self_type* pthis = static_cast<self_type*>(arg);
    pthis->post_thread_state.set_to(TERMINATED);

#if STXXL_STD_THREADS && STXXL_MSVC >= 1700
    // Workaround for deadlock bug in Visual C++ Runtime 2012 and 2013, see
    // request_queue_impl_worker.cpp. -tb
    ExitThread(NULL);
#else
    return NULL;
#endif
}

void* linuxaio_queue::wait_async(void* arg)
{
    (static_cast<linuxaio_queue*>(arg))->wait_requests();

    self_type* pthis = static_cast<self_type*>(arg);
    pthis->wait_thread_state.set_to(TERMINATED);

#if STXXL_STD_THREADS && STXXL_MSVC >= 1700
    // Workaround for deadlock bug in Visual C++ Runtime 2012 and 2013, see
    // request_queue_impl_worker.cpp. -tb
    ExitThread(NULL);
#else
    return NULL;
#endif
}

STXXL_END_NAMESPACE

#endif // #if STXXL_HAVE_LINUXAIO_FILE
// vim: et:ts=4:sw=4
