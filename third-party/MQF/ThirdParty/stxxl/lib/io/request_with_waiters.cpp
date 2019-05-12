/***************************************************************************
 *  lib/io/request_with_waiters.cpp
 *
 *  Part of the STXXL. See http://stxxl.sourceforge.net
 *
 *  Copyright (C) 2002 Roman Dementiev <dementiev@mpi-sb.mpg.de>
 *  Copyright (C) 2008 Andreas Beckmann <beckmann@cs.uni-frankfurt.de>
 *
 *  Distributed under the Boost Software License, Version 1.0.
 *  (See accompanying file LICENSE_1_0.txt or copy at
 *  http://www.boost.org/LICENSE_1_0.txt)
 **************************************************************************/

#include <stxxl/bits/common/mutex.h>
#include <stxxl/bits/common/onoff_switch.h>
#include <stxxl/bits/io/request_with_waiters.h>
#include <stxxl/bits/parallel.h>

#include <algorithm>
#include <functional>

STXXL_BEGIN_NAMESPACE

bool request_with_waiters::add_waiter(onoff_switch* sw)
{
    // this lock needs to be obtained before poll(), otherwise a race
    // condition might occur: the state might change and notify_waiters()
    // could be called between poll() and insert() resulting in waiter sw
    // never being notified
    scoped_mutex_lock lock(m_waiters_mutex);

    if (poll())                     // request already finished
    {
        return true;
    }

    m_waiters.insert(sw);

    return false;
}

void request_with_waiters::delete_waiter(onoff_switch* sw)
{
    scoped_mutex_lock lock(m_waiters_mutex);
    m_waiters.erase(sw);
}

void request_with_waiters::notify_waiters()
{
    scoped_mutex_lock lock(m_waiters_mutex);
    std::for_each(m_waiters.begin(),
                  m_waiters.end(),
                  std::mem_fun(&onoff_switch::on)
                  _STXXL_FORCE_SEQUENTIAL);
}

size_t request_with_waiters::num_waiters()
{
    scoped_mutex_lock lock(m_waiters_mutex);
    return m_waiters.size();
}

STXXL_END_NAMESPACE
// vim: et:ts=4:sw=4
