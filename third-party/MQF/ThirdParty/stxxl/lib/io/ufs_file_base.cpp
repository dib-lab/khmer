/***************************************************************************
 *  lib/io/ufs_file_base.cpp
 *
 *  Part of the STXXL. See http://stxxl.sourceforge.net
 *
 *  Copyright (C) 2002, 2005, 2008 Roman Dementiev <dementiev@mpi-sb.mpg.de>
 *  Copyright (C) 2008 Ilja Andronov <sni4ok@yandex.ru>
 *  Copyright (C) 2008-2010 Andreas Beckmann <beckmann@cs.uni-frankfurt.de>
 *  Copyright (C) 2009 Johannes Singler <singler@ira.uka.de>
 *  Copyright (C) 2013 Timo Bingmann <tb@panthema.net>
 *
 *  Distributed under the Boost Software License, Version 1.0.
 *  (See accompanying file LICENSE_1_0.txt or copy at
 *  http://www.boost.org/LICENSE_1_0.txt)
 **************************************************************************/

#include <stxxl/bits/common/error_handling.h>
#include <stxxl/bits/common/exceptions.h>
#include <stxxl/bits/config.h>
#include <stxxl/bits/io/file.h>
#include <stxxl/bits/io/ufs_file_base.h>
#include <stxxl/bits/verbose.h>
#include "ufs_platform.h"

STXXL_BEGIN_NAMESPACE

const char* ufs_file_base::io_type() const
{
    return "ufs_base";
}

ufs_file_base::ufs_file_base(
    const std::string& filename,
    int mode)
    : file_des(-1), m_mode(mode), filename(filename)
{
    int flags = 0;

    if (mode & RDONLY)
    {
        flags |= O_RDONLY;
    }

    if (mode & WRONLY)
    {
        flags |= O_WRONLY;
    }

    if (mode & RDWR)
    {
        flags |= O_RDWR;
    }

    if (mode & CREAT)
    {
        flags |= O_CREAT;
    }

    if (mode & TRUNC)
    {
        flags |= O_TRUNC;
    }

    if ((mode & DIRECT) || (mode & REQUIRE_DIRECT))
    {
#ifdef __APPLE__
        // no additional open flags are required for Mac OS X
#elif !STXXL_DIRECT_IO_OFF
        flags |= O_DIRECT;
#else
        if (mode & REQUIRE_DIRECT) {
            STXXL_ERRMSG("Error: open()ing " << filename << " with DIRECT mode required, but the system does not support it.");
            file_des = -1;
            return;
        }
        else {
            STXXL_MSG("Warning: open()ing " << filename << " without DIRECT mode, as the system does not support it.");
        }
#endif
    }

    if (mode & SYNC)
    {
        flags |= O_RSYNC;
        flags |= O_DSYNC;
        flags |= O_SYNC;
    }

#if STXXL_WINDOWS
    flags |= O_BINARY;                     // the default in MS is TEXT mode
#endif

#if STXXL_WINDOWS || defined(__MINGW32__)
    const int perms = S_IREAD | S_IWRITE;
#else
    const int perms = S_IREAD | S_IWRITE | S_IRGRP | S_IWGRP;
#endif

    if ((file_des = ::open(filename.c_str(), flags, perms)) >= 0)
    {
        _after_open();
        return;
    }

#if !STXXL_DIRECT_IO_OFF
    if ((mode & DIRECT) && !(mode & REQUIRE_DIRECT) && errno == EINVAL)
    {
        STXXL_MSG("open() error on path=" << filename << " flags=" << flags << ", retrying without O_DIRECT.");

        flags &= ~O_DIRECT;
        m_mode &= ~DIRECT;

        if ((file_des = ::open(filename.c_str(), flags, perms)) >= 0)
        {
            _after_open();
            return;
        }
    }
#endif

    STXXL_THROW_ERRNO(io_error, "open() rc=" << file_des << " path=" << filename << " flags=" << flags);
}

ufs_file_base::~ufs_file_base()
{
    close();
}

void ufs_file_base::_after_open()
{
    // stat file type
#if STXXL_WINDOWS || defined(__MINGW32__)
    struct _stat64 st;
    STXXL_THROW_ERRNO_NE_0(::_fstat64(file_des, &st), io_error,
                           "_fstat64() path=" << filename << " fd=" << file_des);
#else
    struct stat st;
    STXXL_THROW_ERRNO_NE_0(::fstat(file_des, &st), io_error,
                           "fstat() path=" << filename << " fd=" << file_des);
#endif
    m_is_device = S_ISBLK(st.st_mode) ? true : false;

#ifdef __APPLE__
    if (m_mode & REQUIRE_DIRECT) {
        STXXL_THROW_ERRNO_NE_0(fcntl(file_des, F_NOCACHE, 1), io_error,
                               "fcntl() path=" << filename << " fd=" << file_des);
        STXXL_THROW_ERRNO_NE_0(fcntl(file_des, F_RDAHEAD, 0), io_error,
                               "fcntl() path=" << filename << " fd=" << file_des);
    }
    else if (m_mode & DIRECT) {
        if (fcntl(file_des, F_NOCACHE, 1) != 0) {
            STXXL_MSG("fcntl(fd,F_NOCACHE,1) failed on path=" << filename <<
                      " fd=" << file_des << " : " << strerror(errno));
        }
        if (fcntl(file_des, F_RDAHEAD, 0) != 0) {
            STXXL_MSG("fcntl(fd,F_RDAHEAD,0) failed on path=" << filename <<
                      " fd=" << file_des << " : " << strerror(errno));
        }
    }
#endif

    // successfully opened file descriptor
    if (!(m_mode & NO_LOCK))
        lock();
}

void ufs_file_base::close()
{
    scoped_mutex_lock fd_lock(fd_mutex);

    if (file_des == -1)
        return;

    if (::close(file_des) < 0)
        STXXL_THROW_ERRNO(io_error, "close() fd=" << file_des);

    file_des = -1;
}

void ufs_file_base::lock()
{
#if STXXL_WINDOWS || defined(__MINGW32__)
    // not yet implemented
#else
    scoped_mutex_lock fd_lock(fd_mutex);
    struct flock lock_struct;
    lock_struct.l_type = (short)(m_mode & RDONLY ? F_RDLCK : F_RDLCK | F_WRLCK);
    lock_struct.l_whence = SEEK_SET;
    lock_struct.l_start = 0;
    lock_struct.l_len = 0; // lock all bytes
    if ((::fcntl(file_des, F_SETLK, &lock_struct)) < 0)
        STXXL_THROW_ERRNO(io_error, "fcntl(,F_SETLK,) path=" << filename << " fd=" << file_des);
#endif
}

file::offset_type ufs_file_base::_size()
{
    // We use lseek SEEK_END to find the file size. This works for raw devices
    // (where stat() returns zero), and we need not reset the position because
    // serve() always lseek()s before read/write.

    off_t rc = ::lseek(file_des, 0, SEEK_END);
    if (rc < 0)
        STXXL_THROW_ERRNO(io_error, "lseek(fd,0,SEEK_END) path=" << filename << " fd=" << file_des);

    // return value is already the total size
    return rc;
}

file::offset_type ufs_file_base::size()
{
    scoped_mutex_lock fd_lock(fd_mutex);
    return _size();
}

void ufs_file_base::set_size(offset_type newsize)
{
    scoped_mutex_lock fd_lock(fd_mutex);
    return _set_size(newsize);
}

void ufs_file_base::_set_size(offset_type newsize)
{
    offset_type cur_size = _size();

    if (!(m_mode & RDONLY) && !m_is_device)
    {
#if STXXL_WINDOWS || defined(__MINGW32__)
        HANDLE hfile = (HANDLE)::_get_osfhandle(file_des);
        STXXL_THROW_ERRNO_NE_0((hfile == INVALID_HANDLE_VALUE), io_error,
                               "_get_osfhandle() path=" << filename << " fd=" << file_des);

        LARGE_INTEGER desired_pos;
        desired_pos.QuadPart = newsize;

        if (!SetFilePointerEx(hfile, desired_pos, NULL, FILE_BEGIN))
            STXXL_THROW_WIN_LASTERROR(io_error,
                                      "SetFilePointerEx in ufs_file_base::set_size(..) oldsize=" << cur_size <<
                                      " newsize=" << newsize << " ");

        if (!SetEndOfFile(hfile))
            STXXL_THROW_WIN_LASTERROR(io_error,
                                      "SetEndOfFile oldsize=" << cur_size <<
                                      " newsize=" << newsize << " ");
#else
        STXXL_THROW_ERRNO_NE_0(::ftruncate(file_des, newsize), io_error,
                               "ftruncate() path=" << filename << " fd=" << file_des);
#endif
    }

#if !STXXL_WINDOWS
    if (newsize > cur_size)
        STXXL_THROW_ERRNO_LT_0(::lseek(file_des, newsize - 1, SEEK_SET), io_error,
                               "lseek() path=" << filename << " fd=" << file_des << " pos=" << newsize - 1);
#endif
}

void ufs_file_base::close_remove()
{
    close();

    if (m_is_device) {
        STXXL_ERRMSG("remove() path=" << filename << " skipped as file is device node");
        return;
    }

    if (::remove(filename.c_str()) != 0)
        STXXL_ERRMSG("remove() error on path=" << filename << " error=" << strerror(errno));
}

void ufs_file_base::unlink()
{
    if (m_is_device) {
        STXXL_ERRMSG("unlink() path=" << filename << " skipped as file is device node");
        return;
    }

    if (::unlink(filename.c_str()) != 0)
        STXXL_THROW_ERRNO(io_error, "unlink() path=" << filename << " fd=" << file_des);
}

bool ufs_file_base::is_device() const
{
    return m_is_device;
}

STXXL_END_NAMESPACE
// vim: et:ts=4:sw=4
