/***************************************************************************
 *  include/stxxl/bits/common/mutex.h
 *
 *  Part of the STXXL. See http://stxxl.sourceforge.net
 *
 *  Copyright (C) 2002 Roman Dementiev <dementiev@mpi-sb.mpg.de>
 *  Copyright (C) 2008 Andreas Beckmann <beckmann@cs.uni-frankfurt.de>
 *  Copyright (C) 2013 Timo Bingmann <tb@panthema.net>
 *
 *  Distributed under the Boost Software License, Version 1.0.
 *  (See accompanying file LICENSE_1_0.txt or copy at
 *  http://www.boost.org/LICENSE_1_0.txt)
 **************************************************************************/

#ifndef STXXL_COMMON_MUTEX_HEADER
#define STXXL_COMMON_MUTEX_HEADER

#include <stxxl/bits/config.h>
#include <stxxl/bits/namespace.h>

#if STXXL_STD_THREADS && STXXL_WINDOWS && STXXL_MSVC >= 1700
#include <atomic>
#endif

#if STXXL_STD_THREADS
 #include <mutex>
#elif STXXL_BOOST_THREADS
 #include <boost/thread/mutex.hpp>
#elif STXXL_POSIX_THREADS
 #include <pthread.h>

 #include <stxxl/bits/noncopyable.h>
 #include <stxxl/bits/common/error_handling.h>
#else
 #error "Thread implementation not detected."
#endif

STXXL_BEGIN_NAMESPACE

#if STXXL_STD_THREADS

typedef std::mutex mutex;

#elif STXXL_BOOST_THREADS

typedef boost::mutex mutex;

#elif STXXL_POSIX_THREADS

class mutex : private noncopyable
{
    //! mutex handle
    pthread_mutex_t m_mutex;

public:
    //! construct unlocked mutex
    mutex()
    {
        STXXL_CHECK_PTHREAD_CALL(pthread_mutex_init(&m_mutex, NULL));
    }
    //! destroy mutex handle
    ~mutex() noexcept(false)
    {
        // try simple delete first
        int res = pthread_mutex_destroy(&m_mutex);
        if (res == 0) return;

        // try to lock and unlock mutex
        res = pthread_mutex_trylock(&m_mutex);

        if (res == 0 || res == EBUSY) {
            STXXL_CHECK_PTHREAD_CALL(pthread_mutex_unlock(&m_mutex));
        }
        else {
            STXXL_THROW_ERRNO2(resource_error, "pthread_mutex_trylock() failed", res);
        }

        STXXL_CHECK_PTHREAD_CALL(pthread_mutex_destroy(&m_mutex));
    }
    //! lock mutex, may block
    void lock()
    {
        STXXL_CHECK_PTHREAD_CALL(pthread_mutex_lock(&m_mutex));
    }
    //! unlock mutex
    void unlock()
    {
        STXXL_CHECK_PTHREAD_CALL(pthread_mutex_unlock(&m_mutex));
    }
    //! return platform specific handle
    pthread_mutex_t & native_handle()
    {
        return m_mutex;
    }
};

#endif // STXXL_POSIX_THREADS

#if STXXL_STD_THREADS && STXXL_WINDOWS && STXXL_MSVC >= 1700

class spin_lock;
typedef spin_lock fastmutex;

#else

typedef mutex fastmutex;

#endif

#if STXXL_STD_THREADS

typedef std::unique_lock<std::mutex> scoped_mutex_lock;
typedef std::unique_lock<fastmutex> scoped_fast_mutex_lock;

#elif STXXL_BOOST_THREADS

typedef boost::mutex::scoped_lock scoped_mutex_lock;
typedef boost::mutex::scoped_lock scoped_fast_mutex_lock;

#else

//! Aquire a lock that's valid until the end of scope.
class scoped_mutex_lock
{
    //! mutex reference
    mutex& m_mutex;

    //! marker if already unlocked by this thread (needs no synchronization)
    bool is_locked;

public:
    //! lock mutex
    scoped_mutex_lock(mutex& m)
        : m_mutex(m), is_locked(true)
    {
        m_mutex.lock();
    }
    //! unlock mutex hold when object goes out of scope.
    ~scoped_mutex_lock()
    {
        unlock();
    }
    //! unlock mutex hold prematurely
    void unlock()
    {
        if (is_locked) {
            is_locked = false;
            m_mutex.unlock();
        }
    }
    //! return platform specific handle
    pthread_mutex_t & native_handle()
    {
        return m_mutex.native_handle();
    }
};

typedef scoped_mutex_lock scoped_fast_mutex_lock;

#endif

#if STXXL_STD_THREADS && STXXL_WINDOWS && STXXL_MSVC >= 1700

class spin_lock
{
public:
#if STXXL_MSVC < 1800
    spin_lock()
    {
        lck.clear(std::memory_order_release);
    }
#else
    spin_lock()
    { }
#endif

    void lock()
    {
        while (lck.test_and_set(std::memory_order_acquire))
        { }
    }

    void unlock()
    {
        lck.clear(std::memory_order_release);
    }

private:
#if STXXL_MSVC >= 1800
    std::atomic_flag lck = ATOMIC_FLAG_INIT;
    spin_lock(const spin_lock&) = delete;
    spin_lock& operator = (const spin_lock&) = delete;
#else
    std::atomic_flag lck;
    spin_lock(const spin_lock&);
    spin_lock& operator = (const spin_lock&);
#endif
};

#endif
STXXL_END_NAMESPACE

#endif // !STXXL_COMMON_MUTEX_HEADER
