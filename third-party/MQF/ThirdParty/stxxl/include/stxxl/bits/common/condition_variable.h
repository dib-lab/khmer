/***************************************************************************
 *  include/stxxl/bits/common/condition_variable.h
 *
 *  Part of the STXXL. See http://stxxl.sourceforge.net
 *
 *  Copyright (C) 2013 Timo Bingmann <tb@panthema.net>
 *
 *  Distributed under the Boost Software License, Version 1.0.
 *  (See accompanying file LICENSE_1_0.txt or copy at
 *  http://www.boost.org/LICENSE_1_0.txt)
 **************************************************************************/

#ifndef STXXL_COMMON_CONDITION_VARIABLE_HEADER
#define STXXL_COMMON_CONDITION_VARIABLE_HEADER

#include <stxxl/bits/config.h>
#include <stxxl/bits/namespace.h>
#include <stxxl/bits/common/mutex.h>

#if STXXL_STD_THREADS
 #include <condition_variable>
#elif STXXL_BOOST_THREADS
 #include <boost/thread/condition.hpp>
#elif STXXL_POSIX_THREADS
 #include <pthread.h>
 #include <cerrno>

 #include <stxxl/bits/noncopyable.h>
 #include <stxxl/bits/common/error_handling.h>
#else
 #error "Thread implementation not detected."
#endif

STXXL_BEGIN_NAMESPACE

#if STXXL_STD_THREADS

typedef std::condition_variable condition_variable;

#elif STXXL_BOOST_THREADS

typedef boost::condition condition_variable;

#elif STXXL_POSIX_THREADS

class condition_variable : private noncopyable
{
    //! pthread handle to condition
    pthread_cond_t cond;

public:
    //! initialize condition variable
    condition_variable()
    {
        STXXL_CHECK_PTHREAD_CALL(pthread_cond_init(&cond, NULL));
    }
    //! destroy condition variable
    ~condition_variable() noexcept(false)
    {
        STXXL_CHECK_PTHREAD_CALL(pthread_cond_destroy(&cond));
    }
    //! notify one waiting thread
    void notify_one()
    {
        STXXL_CHECK_PTHREAD_CALL(pthread_cond_signal(&cond));
    }
    //! notify all waiting threads
    void notify_all()
    {
        STXXL_CHECK_PTHREAD_CALL(pthread_cond_broadcast(&cond));
    }
    //! wait for a signal on the condition variable
    void wait(scoped_mutex_lock& lock)
    {
        STXXL_CHECK_PTHREAD_CALL(pthread_cond_wait(&cond, &lock.native_handle()));
    }
};

#endif

STXXL_END_NAMESPACE

#endif // !STXXL_COMMON_CONDITION_VARIABLE_HEADER
