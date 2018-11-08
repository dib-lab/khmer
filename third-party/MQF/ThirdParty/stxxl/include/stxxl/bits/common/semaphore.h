/***************************************************************************
 *  include/stxxl/bits/common/semaphore.h
 *
 *  Part of the STXXL. See http://stxxl.sourceforge.net
 *
 *  Copyright (C) 2002 Roman Dementiev <dementiev@mpi-sb.mpg.de>
 *  Copyright (C) 2013 Timo Bingmann <tb@panthema.net>
 *
 *  Distributed under the Boost Software License, Version 1.0.
 *  (See accompanying file LICENSE_1_0.txt or copy at
 *  http://www.boost.org/LICENSE_1_0.txt)
 **************************************************************************/

#ifndef STXXL_COMMON_SEMAPHORE_HEADER
#define STXXL_COMMON_SEMAPHORE_HEADER

#include <stxxl/bits/noncopyable.h>
#include <stxxl/bits/common/mutex.h>
#include <stxxl/bits/common/condition_variable.h>

STXXL_BEGIN_NAMESPACE

class semaphore : private noncopyable
{
    //! value of the semaphore
    int v;

    //! mutex for condition variable
    mutex m_mutex;

    //! condition variable
    condition_variable m_cond;

public:
    //! construct semaphore
    semaphore(int init_value = 1)
        : v(init_value)
    { }
    //! function increments the semaphore and signals any threads that are
    //! blocked waiting a change in the semaphore
    int operator ++ (int)
    {
        scoped_mutex_lock lock(m_mutex);
        int res = ++v;
        lock.unlock();
        m_cond.notify_one();
        return res;
    }
    //! function decrements the semaphore and blocks if the semaphore is <= 0
    //! until another thread signals a change
    int operator -- (int)
    {
        scoped_mutex_lock lock(m_mutex);
        while (v <= 0)
            m_cond.wait(lock);

        return --v;
    }
    //! function does NOT block but simply decrements the semaphore should not
    //! be used instead of down -- only for programs where multiple threads
    //! must up on a semaphore before another thread can go down, i.e., allows
    //! programmer to set the semaphore to a negative value prior to using it
    //! for synchronization.
    int decrement()
    {
        scoped_mutex_lock lock(m_mutex);
        return --v;
    }
#if 0
    //! function returns the value of the semaphore at the time the
    //! critical section is accessed.  obviously the value is not guaranteed
    //! after the function unlocks the critical section.
    int get_value()
    {
        scoped_mutex_lock lock(m_mutex);
        return v;
    }
#endif
};

STXXL_END_NAMESPACE

#endif // !STXXL_COMMON_SEMAPHORE_HEADER
