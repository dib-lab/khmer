/***************************************************************************
 *  include/stxxl/bits/io/completion_handler.h
 *
 *  Loki-style completion handler (functors)
 *
 *  Part of the STXXL. See http://stxxl.sourceforge.net
 *
 *  Copyright (C) 2003 Roman Dementiev <dementiev@mpi-sb.mpg.de>
 *  Copyright (C) 2008 Andreas Beckmann <beckmann@cs.uni-frankfurt.de>
 *
 *  Distributed under the Boost Software License, Version 1.0.
 *  (See accompanying file LICENSE_1_0.txt or copy at
 *  http://www.boost.org/LICENSE_1_0.txt)
 **************************************************************************/

#ifndef STXXL_IO_COMPLETION_HANDLER_HEADER
#define STXXL_IO_COMPLETION_HANDLER_HEADER

#include <stxxl/bits/namespace.h>
#include <stxxl/bits/compat/unique_ptr.h>
#include <cstdlib>

STXXL_BEGIN_NAMESPACE

class request;

class completion_handler_impl
{
public:
    virtual void operator () (request*) = 0;
    virtual completion_handler_impl * clone() const = 0;
    virtual ~completion_handler_impl() { }
};

template <typename HandlerType>
class completion_handler1 : public completion_handler_impl
{
private:
    HandlerType m_handler;

public:
    completion_handler1(const HandlerType& handler)
        : m_handler(handler)
    { }
    completion_handler1 * clone() const
    {
        return new completion_handler1(*this);
    }
    void operator () (request* req)
    {
        m_handler(req);
    }
};

//! Completion handler class (Loki-style).
//!
//! In some situations one needs to execute some actions after completion of an
//! I/O request. In these cases one can use an I/O completion handler - a
//! function object that can be passed as a parameter to asynchronous I/O calls
//! \c stxxl::file::aread and \c stxxl::file::awrite .
class completion_handler
{
    compat_unique_ptr<completion_handler_impl>::result m_ptr;

public:
    //! Construct default, no operation completion handler.
    completion_handler()
        : m_ptr(static_cast<completion_handler_impl*>(NULL))
    { }

    //! Copy constructor.
    completion_handler(const completion_handler& obj)
        : m_ptr(obj.m_ptr.get() ? obj.m_ptr.get()->clone() : NULL)
    { }

    //! Construct a completion handler which calls some function.
    template <typename HandlerType>
    completion_handler(const HandlerType& handler)
        : m_ptr(new completion_handler1<HandlerType>(handler))
    { }

    //! Assignment operator
    completion_handler& operator = (const completion_handler& obj)
    {
        m_ptr.reset(obj.m_ptr.get() ? obj.m_ptr.get()->clone() : NULL);
        return *this;
    }

    //! Call the enclosed completion handler.
    void operator () (request* req)
    {
        if (m_ptr.get())
            (*m_ptr)(req);
    }
};

STXXL_END_NAMESPACE

#endif // !STXXL_IO_COMPLETION_HANDLER_HEADER
// vim: et:ts=4:sw=4
