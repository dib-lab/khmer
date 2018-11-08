/***************************************************************************
 *  include/stxxl/bits/io/request.h
 *
 *  Part of the STXXL. See http://stxxl.sourceforge.net
 *
 *  Copyright (C) 2002 Roman Dementiev <dementiev@mpi-sb.mpg.de>
 *  Copyright (C) 2008 Andreas Beckmann <beckmann@cs.uni-frankfurt.de>
 *  Copyright (C) 2013-2014 Timo Bingmann <tb@panthema.net>
 *
 *  Distributed under the Boost Software License, Version 1.0.
 *  (See accompanying file LICENSE_1_0.txt or copy at
 *  http://www.boost.org/LICENSE_1_0.txt)
 **************************************************************************/

#ifndef STXXL_IO_REQUEST_HEADER
#define STXXL_IO_REQUEST_HEADER

#include <cassert>

#include <stxxl/bits/namespace.h>
#include <stxxl/bits/io/request_interface.h>
#include <stxxl/bits/common/counting_ptr.h>
#include <stxxl/bits/common/exceptions.h>
#include <stxxl/bits/io/completion_handler.h>
#include <stxxl/bits/compat/unique_ptr.h>
#include <stxxl/bits/verbose.h>

STXXL_BEGIN_NAMESPACE

//! \addtogroup reqlayer
//! \{

#define STXXL_BLOCK_ALIGN 4096

class file;

//! Request object encapsulating basic properties like file and offset.
class request : virtual public request_interface, public atomic_counted_object
{
    friend class linuxaio_queue;

protected:
    completion_handler m_on_complete;
    compat_unique_ptr<stxxl::io_error>::result m_error;

protected:
    file* m_file;
    void* m_buffer;
    offset_type m_offset;
    size_type m_bytes;
    request_type m_type;

public:
    request(const completion_handler& on_compl,
            file* file,
            void* buffer,
            offset_type offset,
            size_type bytes,
            request_type type);

    virtual ~request() noexcept(false);

    file * get_file() const { return m_file; }
    void * get_buffer() const { return m_buffer; }
    offset_type get_offset() const { return m_offset; }
    size_type get_size() const { return m_bytes; }
    request_type get_type() const { return m_type; }

    void check_alignment() const;

    std::ostream & print(std::ostream& out) const;

    //! Inform the request object that an error occurred during the I/O
    //! execution.
    void error_occured(const char* msg)
    {
        m_error.reset(new stxxl::io_error(msg));
    }

    //! Inform the request object that an error occurred during the I/O
    //! execution.
    void error_occured(const std::string& msg)
    {
        m_error.reset(new stxxl::io_error(msg));
    }

    //! Rises an exception if there were error with the I/O.
    void check_errors()
    {
        if (m_error.get())
            throw *(m_error.get());
    }

    virtual const char * io_type() const;

protected:
    void check_nref(bool after = false)
    {
        if (get_reference_count() < 2)
            check_nref_failed(after);
    }

private:
    void check_nref_failed(bool after);
};

inline std::ostream& operator << (std::ostream& out, const request& req)
{
    return req.print(out);
}

//! A reference counting pointer for \c request.
typedef counting_ptr<request> request_ptr;

//! \}

STXXL_END_NAMESPACE

#endif // !STXXL_IO_REQUEST_HEADER
// vim: et:ts=4:sw=4
