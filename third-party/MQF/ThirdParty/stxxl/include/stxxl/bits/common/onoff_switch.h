/***************************************************************************
 *  include/stxxl/bits/common/onoff_switch.h
 *
 *  Kind of binary semaphore: initially OFF, then multiple waiters can attach
 *  to the switch, which get notified one-by-one when switched ON.
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

#ifndef STXXL_COMMON_ONOFF_SWITCH_HEADER
#define STXXL_COMMON_ONOFF_SWITCH_HEADER

#include <stxxl/bits/config.h>

#include <stxxl/bits/noncopyable.h>
#include <stxxl/bits/common/mutex.h>
#include <stxxl/bits/common/condition_variable.h>

STXXL_BEGIN_NAMESPACE

class onoff_switch : private noncopyable
{
    //! mutex for condition variable
    mutex m_mutex;

    //! condition variable
    condition_variable m_cond;

    //! the switch's state
    bool m_on;

public:
    //! construct switch
    onoff_switch(bool flag = false)
        : m_on(flag)
    { }
    //! turn switch ON and notify one waiter
    void on()
    {
        scoped_mutex_lock lock(m_mutex);
        m_on = true;
        lock.unlock();
        m_cond.notify_one();
    }
    //! turn switch OFF and notify one waiter
    void off()
    {
        scoped_mutex_lock lock(m_mutex);
        m_on = false;
        lock.unlock();
        m_cond.notify_one();
    }
    //! wait for switch to turn ON
    void wait_for_on()
    {
        scoped_mutex_lock lock(m_mutex);
        if (!m_on)
            m_cond.wait(lock);
    }
    //! wait for switch to turn OFF
    void wait_for_off()
    {
        scoped_mutex_lock lock(m_mutex);
        if (m_on)
            m_cond.wait(lock);
    }
    //! return true if switch is ON
    bool is_on()
    {
        scoped_mutex_lock lock(m_mutex);
        return m_on;
    }
};

STXXL_END_NAMESPACE

#endif // !STXXL_COMMON_ONOFF_SWITCH_HEADER
