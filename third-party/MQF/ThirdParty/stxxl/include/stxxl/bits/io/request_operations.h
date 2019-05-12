/***************************************************************************
 *  include/stxxl/bits/io/request_operations.h
 *
 *  Part of the STXXL. See http://stxxl.sourceforge.net
 *
 *  Copyright (C) 2002 Roman Dementiev <dementiev@mpi-sb.mpg.de>
 *  Copyright (C) 2008, 2009 Andreas Beckmann <beckmann@cs.uni-frankfurt.de>
 *  Copyright (C) 2009 Johannes Singler <singler@ira.uka.de>
 *
 *  Distributed under the Boost Software License, Version 1.0.
 *  (See accompanying file LICENSE_1_0.txt or copy at
 *  http://www.boost.org/LICENSE_1_0.txt)
 **************************************************************************/

#ifndef STXXL_IO_REQUEST_OPERATIONS_HEADER
#define STXXL_IO_REQUEST_OPERATIONS_HEADER

#include <stxxl/bits/namespace.h>
#include <stxxl/bits/io/request.h>
#include <stxxl/bits/io/iostats.h>
#include <stxxl/bits/common/onoff_switch.h>

STXXL_BEGIN_NAMESPACE

//! \addtogroup reqlayer
//! \{

//! Collection of functions to track statuses of a number of requests.

//! Suspends calling thread until \b all given requests are completed.
//! \param reqs_begin begin of request sequence to wait for
//! \param reqs_end end of request sequence to wait for
template <class RequestIterator>
void wait_all(RequestIterator reqs_begin, RequestIterator reqs_end)
{
    for ( ; reqs_begin != reqs_end; ++reqs_begin)
        (request_ptr(*reqs_begin))->wait();
}

//! Suspends calling thread until \b all given requests are completed.
//! \param req_array array of request_ptr objects
//! \param count size of req_array
inline void wait_all(request_ptr req_array[], size_t count)
{
    wait_all(req_array, req_array + count);
}

//! Cancel requests.
//! The specified requests are canceled unless already being processed.
//! However, cancelation cannot be guaranteed.
//! Cancelled requests must still be waited for in order to ensure correct
//! operation.
//! \param reqs_begin begin of request sequence
//! \param reqs_end end of request sequence
//! \return number of request canceled
template <class RequestIterator>
typename std::iterator_traits<RequestIterator>::difference_type
cancel_all(RequestIterator reqs_begin, RequestIterator reqs_end)
{
    typename std::iterator_traits<RequestIterator>::difference_type num_canceled = 0;
    while (reqs_begin != reqs_end)
    {
        if ((request_ptr(*reqs_begin))->cancel())
            ++num_canceled;
        ++reqs_begin;
    }
    return num_canceled;
}

//! Polls requests.
//! \param reqs_begin begin of request sequence to poll
//! \param reqs_end end of request sequence to poll
//! \return \c true if any of requests is completed, then index contains valid value, otherwise \c false
template <class RequestIterator>
RequestIterator poll_any(RequestIterator reqs_begin, RequestIterator reqs_end)
{
    while (reqs_begin != reqs_end)
    {
        if ((request_ptr(*reqs_begin))->poll())
            return reqs_begin;

        ++reqs_begin;
    }
    return reqs_end;
}

//! Polls requests.
//! \param req_array array of request_ptr objects
//! \param count size of req_array
//! \param index contains index of the \b first completed request if any
//! \return \c true if any of requests is completed, then index contains valid value, otherwise \c false
inline bool poll_any(request_ptr req_array[], size_t count, size_t& index)
{
    request_ptr* res = poll_any(req_array, req_array + count);
    index = res - req_array;
    return res != (req_array + count);
}

//! Suspends calling thread until \b any of requests is completed.
//! \param reqs_begin begin of request sequence to wait for
//! \param reqs_end end of request sequence to wait for
//! \return index in req_array pointing to the \b first completed request
template <class RequestIterator>
RequestIterator wait_any(RequestIterator reqs_begin, RequestIterator reqs_end)
{
    stats::scoped_wait_timer wait_timer(stats::WAIT_OP_ANY);

    onoff_switch sw;

    RequestIterator cur = reqs_begin, result = reqs_end;

    for ( ; cur != reqs_end; cur++)
    {
        if ((request_ptr(*cur))->add_waiter(&sw))
        {
            // request is already done, no waiter was added to the request
            result = cur;

            if (cur != reqs_begin)
            {
                while (--cur != reqs_begin)
                    (request_ptr(*cur))->delete_waiter(&sw);

                (request_ptr(*cur))->delete_waiter(&sw);
            }

            (request_ptr(*result))->check_errors();

            return result;
        }
    }

    sw.wait_for_on();

    for (cur = reqs_begin; cur != reqs_end; cur++)
    {
        (request_ptr(*cur))->delete_waiter(&sw);
        if (result == reqs_end && (request_ptr(*cur))->poll())
            result = cur;
    }

    return result;
}

//! Suspends calling thread until \b any of requests is completed.
//! \param req_array array of \c request_ptr objects
//! \param count size of req_array
//! \return index in req_array pointing to the \b first completed request
inline size_t wait_any(request_ptr req_array[], size_t count)
{
    return wait_any(req_array, req_array + count) - req_array;
}

//! \}

STXXL_END_NAMESPACE

#endif // !STXXL_IO_REQUEST_OPERATIONS_HEADER
// vim: et:ts=4:sw=4
