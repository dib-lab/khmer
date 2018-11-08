/***************************************************************************
 *  lib/io/syscall_file.cpp
 *
 *  Part of the STXXL. See http://stxxl.sourceforge.net
 *
 *  Copyright (C) 2002 Roman Dementiev <dementiev@mpi-sb.mpg.de>
 *  Copyright (C) 2008, 2010 Andreas Beckmann <beckmann@cs.uni-frankfurt.de>
 *  Copyright (C) 2010 Johannes Singler <singler@kit.edu>
 *
 *  Distributed under the Boost Software License, Version 1.0.
 *  (See accompanying file LICENSE_1_0.txt or copy at
 *  http://www.boost.org/LICENSE_1_0.txt)
 **************************************************************************/

#include <stxxl/bits/common/error_handling.h>
#include <stxxl/bits/common/mutex.h>
#include <stxxl/bits/config.h>
#include <stxxl/bits/io/iostats.h>
#include <stxxl/bits/io/request.h>
#include <stxxl/bits/io/request_interface.h>
#include <stxxl/bits/io/syscall_file.h>
#include "ufs_platform.h"

STXXL_BEGIN_NAMESPACE

void syscall_file::serve(void* buffer, offset_type offset, size_type bytes,
                         request::request_type type)
{
    scoped_mutex_lock fd_lock(fd_mutex);

    char* cbuffer = static_cast<char*>(buffer);

    stats::scoped_read_write_timer read_write_timer(bytes, type == request::WRITE);

    while (bytes > 0)
    {
        off_t rc = ::lseek(file_des, offset, SEEK_SET);
        if (rc < 0)
        {
            STXXL_THROW_ERRNO
                (io_error,
                " this=" << this <<
                " call=::lseek(fd,offset,SEEK_SET)" <<
                " path=" << filename <<
                " fd=" << file_des <<
                " offset=" << offset <<
                " buffer=" << cbuffer <<
                " bytes=" << bytes <<
                " type=" << ((type == request::READ) ? "READ" : "WRITE") <<
                " rc=" << rc);
        }

        if (type == request::READ)
        {
#if STXXL_MSVC
            assert(bytes <= std::numeric_limits<unsigned int>::max());
            if ((rc = ::read(file_des, cbuffer, (unsigned int)bytes)) <= 0)
#else
            if ((rc = ::read(file_des, cbuffer, bytes)) <= 0)
#endif
            {
                STXXL_THROW_ERRNO
                    (io_error,
                    " this=" << this <<
                    " call=::read(fd,buffer,bytes)" <<
                    " path=" << filename <<
                    " fd=" << file_des <<
                    " offset=" << offset <<
                    " buffer=" << buffer <<
                    " bytes=" << bytes <<
                    " type=" << "READ" <<
                    " rc=" << rc);
            }
            bytes = (size_type)(bytes - rc);
            offset += rc;
            cbuffer += rc;

            if (bytes > 0 && offset == this->_size())
            {
                // read request extends past end-of-file
                // fill reminder with zeroes
                memset(cbuffer, 0, bytes);
                bytes = 0;
            }
        }
        else
        {
#if STXXL_MSVC
            assert(bytes <= std::numeric_limits<unsigned int>::max());
            if ((rc = ::write(file_des, cbuffer, (unsigned int)bytes)) <= 0)
#else
            if ((rc = ::write(file_des, cbuffer, bytes)) <= 0)
#endif
            {
                STXXL_THROW_ERRNO
                    (io_error,
                    " this=" << this <<
                    " call=::write(fd,buffer,bytes)" <<
                    " path=" << filename <<
                    " fd=" << file_des <<
                    " offset=" << offset <<
                    " buffer=" << buffer <<
                    " bytes=" << bytes <<
                    " type=" << "WRITE" <<
                    " rc=" << rc);
            }
            bytes = (size_type)(bytes - rc);
            offset += rc;
            cbuffer += rc;
        }
    }
}

const char* syscall_file::io_type() const
{
    return "syscall";
}

STXXL_END_NAMESPACE
// vim: et:ts=4:sw=4
