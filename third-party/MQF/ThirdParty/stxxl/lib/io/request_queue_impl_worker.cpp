/***************************************************************************
 *  lib/io/request_queue_impl_worker.cpp
 *
 *  Part of the STXXL. See http://stxxl.sourceforge.net
 *
 *  Copyright (C) 2002-2005 Roman Dementiev <dementiev@mpi-sb.mpg.de>
 *  Copyright (C) 2008, 2009 Andreas Beckmann <beckmann@cs.uni-frankfurt.de>
 *  Copyright (C) 2009 Johannes Singler <singler@ira.uka.de>
 *  Copyright (C) 2013 Timo Bingmann <tb@panthema.net>
 *
 *  Distributed under the Boost Software License, Version 1.0.
 *  (See accompanying file LICENSE_1_0.txt or copy at
 *  http://www.boost.org/LICENSE_1_0.txt)
 **************************************************************************/

#include <stxxl/bits/common/error_handling.h>
#include <stxxl/bits/common/semaphore.h>
#include <stxxl/bits/common/state.h>
#include <stxxl/bits/config.h>
#include <stxxl/bits/io/request_queue_impl_worker.h>
#include <stxxl/bits/namespace.h>

#include <cassert>
#include <cstddef>

#if STXXL_BOOST_THREADS
 #include <boost/bind.hpp>
#endif
#if STXXL_STD_THREADS && STXXL_MSVC >= 1700
 #include <windows.h>
#endif

STXXL_BEGIN_NAMESPACE

void request_queue_impl_worker::start_thread(void* (*worker)(void*), void* arg, thread_type& t, state<thread_state>& s)
{
    assert(s() == NOT_RUNNING);
#if STXXL_STD_THREADS
    t = new std::thread(worker, arg);
#elif STXXL_BOOST_THREADS
    t = new boost::thread(boost::bind(worker, arg));
#else
    STXXL_CHECK_PTHREAD_CALL(pthread_create(&t, NULL, worker, arg));
#endif
    s.set_to(RUNNING);
}

void request_queue_impl_worker::stop_thread(thread_type& t, state<thread_state>& s, semaphore& sem)
{
    assert(s() == RUNNING);
    s.set_to(TERMINATING);
    sem++;
#if STXXL_STD_THREADS
#if STXXL_MSVC >= 1700
    // In the Visual C++ Runtime 2012 and 2013, there is a deadlock bug, which
    // occurs when threads are joined after main() exits. Apparently, Microsoft
    // thinks this is not a big issue. It has not been fixed in VC++RT 2013.
    // https://connect.microsoft.com/VisualStudio/feedback/details/747145
    //
    // All STXXL threads are created by singletons, which are global variables
    // that are deleted after main() exits. The fix applied here it to use
    // std::thread::native_handle() and access the WINAPI to terminate the
    // thread directly (after it finished handling its i/o requests).

    WaitForSingleObject(t->native_handle(), INFINITE);
    CloseHandle(t->native_handle());
#else
    t->join();
    delete t;
#endif
    t = NULL;
#elif STXXL_BOOST_THREADS
    t->join();
    delete t;
    t = NULL;
#else
    STXXL_CHECK_PTHREAD_CALL(pthread_join(t, NULL));
#endif
    assert(s() == TERMINATED);
    s.set_to(NOT_RUNNING);
}

STXXL_END_NAMESPACE
// vim: et:ts=4:sw=4
