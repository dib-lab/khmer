/***************************************************************************
 *  include/stxxl/bits/io/request_with_waiters.h
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

#ifndef STXXL_IO_REQUEST_WITH_WAITERS_HEADER
#define STXXL_IO_REQUEST_WITH_WAITERS_HEADER

#include <set>

#include <stxxl/bits/common/mutex.h>
#include <stxxl/bits/common/onoff_switch.h>
#include <stxxl/bits/io/request.h>
#include <stxxl/bits/namespace.h>

STXXL_BEGIN_NAMESPACE

//! \addtogroup reqlayer
//! \{

//! Request that is aware of threads waiting for it to complete.
class request_with_waiters : public request
{
    mutex m_waiters_mutex;
    std::set<onoff_switch*> m_waiters;

protected:
    bool add_waiter(onoff_switch* sw);
    void delete_waiter(onoff_switch* sw);
    void notify_waiters();

    //! returns number of waiters
    size_t num_waiters();

public:
    request_with_waiters(
        const completion_handler& on_cmpl,
        file* f,
        void* buf,
        offset_type off,
        size_type b,
        request_type t)
        : request(on_cmpl, f, buf, off, b, t)
    { }
};

//! \}

STXXL_END_NAMESPACE

#endif // !STXXL_IO_REQUEST_WITH_WAITERS_HEADER
// vim: et:ts=4:sw=4
