/***************************************************************************
 *  include/stxxl/bits/io/linuxaio_request.h
 *
 *  Part of the STXXL. See http://stxxl.sourceforge.net
 *
 *  Copyright (C) 2011 Johannes Singler <singler@kit.edu>
 *  Copyright (C) 2014 Timo Bingmann <tb@panthema.net>
 *
 *  Distributed under the Boost Software License, Version 1.0.
 *  (See accompanying file LICENSE_1_0.txt or copy at
 *  http://www.boost.org/LICENSE_1_0.txt)
 **************************************************************************/

#ifndef STXXL_IO_LINUXAIO_REQUEST_HEADER
#define STXXL_IO_LINUXAIO_REQUEST_HEADER

#include <stxxl/bits/io/linuxaio_file.h>

#if STXXL_HAVE_LINUXAIO_FILE

#include <linux/aio_abi.h>
#include <stxxl/bits/io/request_with_state.h>

#define STXXL_VERBOSE_LINUXAIO(msg) STXXL_VERBOSE2(msg)

STXXL_BEGIN_NAMESPACE

//! \addtogroup reqlayer
//! \{

//! Request for an linuxaio_file.
class linuxaio_request : public request_with_state
{
    template <class base_file_type>
    friend class fileperblock_file;

    //! control block of async request
    iocb cb;

    void fill_control_block();

public:
    linuxaio_request(
        const completion_handler& on_cmpl,
        file* file,
        void* buffer,
        offset_type offset,
        size_type bytes,
        request_type type)
        : request_with_state(on_cmpl, file, buffer, offset, bytes, type)
    {
        assert(dynamic_cast<linuxaio_file*>(file));
        STXXL_VERBOSE_LINUXAIO("linuxaio_request[" << this << "]" <<
                               " linuxaio_request" <<
                               "(file=" << file << " buffer=" << buffer <<
                               " offset=" << offset << " bytes=" << bytes <<
                               " type=" << type << ")");
    }

    bool post();
    bool cancel();
    bool cancel_aio();
    void completed(bool posted, bool canceled);
    void completed(bool canceled) { completed(true, canceled); }
};

//! \}

STXXL_END_NAMESPACE

#endif // #if STXXL_HAVE_LINUXAIO_FILE

#endif // !STXXL_IO_LINUXAIO_REQUEST_HEADER
// vim: et:ts=4:sw=4
