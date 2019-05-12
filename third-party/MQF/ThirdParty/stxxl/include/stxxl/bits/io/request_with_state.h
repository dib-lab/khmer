/***************************************************************************
 *  include/stxxl/bits/io/request_with_state.h
 *
 *  Part of the STXXL. See http://stxxl.sourceforge.net
 *
 *  Copyright (C) 2008 Andreas Beckmann <beckmann@cs.uni-frankfurt.de>
 *  Copyright (C) 2009 Johannes Singler <singler@ira.uka.de>
 *
 *  Distributed under the Boost Software License, Version 1.0.
 *  (See accompanying file LICENSE_1_0.txt or copy at
 *  http://www.boost.org/LICENSE_1_0.txt)
 **************************************************************************/

#ifndef STXXL_IO_REQUEST_WITH_STATE_HEADER
#define STXXL_IO_REQUEST_WITH_STATE_HEADER

#include <stxxl/bits/common/state.h>
#include <stxxl/bits/io/request.h>
#include <stxxl/bits/io/request_with_waiters.h>
#include <stxxl/bits/namespace.h>

STXXL_BEGIN_NAMESPACE

//! \addtogroup reqlayer
//! \{

//! Request with completion state.
class request_with_state : public request_with_waiters
{
protected:
    //! states of request
    //! OP - operating, DONE - request served, READY2DIE - can be destroyed
    enum request_state { OP = 0, DONE = 1, READY2DIE = 2 };

    state<request_state> m_state;

protected:
    request_with_state(
        const completion_handler& on_cmpl,
        file* f,
        void* buf,
        offset_type off,
        size_type b,
        request_type t)
        : request_with_waiters(on_cmpl, f, buf, off, b, t),
          m_state(OP)
    { }

public:
    virtual ~request_with_state();
    void wait(bool measure_time = true);
    bool poll();
    bool cancel();

protected:
    void completed(bool canceled);
};

//! \}

STXXL_END_NAMESPACE

#endif // !STXXL_IO_REQUEST_WITH_STATE_HEADER
// vim: et:ts=4:sw=4
