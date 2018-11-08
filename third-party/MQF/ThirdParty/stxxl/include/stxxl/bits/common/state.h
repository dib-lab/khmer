/***************************************************************************
 *  include/stxxl/bits/common/state.h
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

#ifndef STXXL_COMMON_STATE_HEADER
#define STXXL_COMMON_STATE_HEADER

#include <stxxl/bits/noncopyable.h>
#include <stxxl/bits/common/mutex.h>
#include <stxxl/bits/common/condition_variable.h>

STXXL_BEGIN_NAMESPACE

template <typename ValueType = int>
class state : private noncopyable
{
    typedef ValueType value_type;

    //! mutex for condition variable
    mutex m_mutex;

    //! condition variable
    condition_variable m_cond;

    //! current state
    value_type m_state;

public:
    state(const value_type& s)
        : m_state(s)
    { }

    void set_to(const value_type& new_state)
    {
        scoped_mutex_lock lock(m_mutex);
        m_state = new_state;
        lock.unlock();
        m_cond.notify_all();
    }

    void wait_for(const value_type& needed_state)
    {
        scoped_mutex_lock lock(m_mutex);
        while (needed_state != m_state)
            m_cond.wait(lock);
    }

    value_type operator () ()
    {
        scoped_mutex_lock lock(m_mutex);
        return m_state;
    }
};

STXXL_END_NAMESPACE

#endif // !STXXL_COMMON_STATE_HEADER
