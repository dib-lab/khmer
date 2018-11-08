/***************************************************************************
 *  lib/io/boostfd_file.cpp
 *
 *  Part of the STXXL. See http://stxxl.sourceforge.net
 *
 *  Copyright (C) 2006 Roman Dementiev <dementiev@ira.uka.de>
 *  Copyright (C) 2009, 2010 Johannes Singler <singler@kit.edu>
 *  Copyright (C) 2008, 2010 Andreas Beckmann <beckmann@cs.uni-frankfurt.de>
 *  Copyright (C) 2014 Timo Bingmann <tb@panthema.net>
 *
 *  Distributed under the Boost Software License, Version 1.0.
 *  (See accompanying file LICENSE_1_0.txt or copy at
 *  http://www.boost.org/LICENSE_1_0.txt)
 **************************************************************************/

#include <stxxl/bits/io/boostfd_file.h>

#if STXXL_HAVE_BOOSTFD_FILE

#include <stxxl/bits/io/iostats.h>
#include <stxxl/bits/common/error_handling.h>

#include <boost/filesystem/operations.hpp>
#include <boost/filesystem/fstream.hpp>
#include <boost/version.hpp>

STXXL_BEGIN_NAMESPACE

void boostfd_file::serve(void* buffer, offset_type offset, size_type bytes,
                         request::request_type type)
{
    scoped_mutex_lock fd_lock(m_fd_mutex);

    try
    {
        m_file_des.seek(offset, BOOST_IOS::beg);
    }
    catch (const std::exception& ex)
    {
        STXXL_THROW_ERRNO
            (io_error,
            "Error doing seek() in boostfd_request::serve()" <<
            " offset=" << offset <<
            " this=" << this <<
            " buffer=" << buffer <<
            " bytes=" << bytes <<
            " type=" << ((type == request::READ) ? "READ" : "WRITE") <<
            " : " << ex.what());
    }

    stats::scoped_read_write_timer read_write_timer(bytes, type == request::WRITE);

    if (type == request::READ)
    {
        try
        {
            std::streamsize rc = m_file_des.read((char*)buffer, bytes);
            if (rc != std::streamsize(bytes)) {
                STXXL_THROW_ERRNO(io_error, " partial read: " << rc << " missing " << (bytes - rc) << " out of " << bytes << " bytes");
            }
        }
        catch (const std::exception& ex)
        {
            STXXL_THROW_ERRNO
                (io_error,
                "Error doing read() in boostfd_request::serve()" <<
                " offset=" << offset <<
                " this=" << this <<
                " buffer=" << buffer <<
                " bytes=" << bytes <<
                " type=" << ((type == request::READ) ? "READ" : "WRITE") <<
                " : " << ex.what());
        }
    }
    else
    {
        try
        {
            std::streamsize rc = m_file_des.write((char*)buffer, bytes);
            if (rc != std::streamsize(bytes)) {
                STXXL_THROW_ERRNO(io_error, " partial write: " << rc << " missing " << (bytes - rc) << " out of " << bytes << " bytes");
            }
        }
        catch (const std::exception& ex)
        {
            STXXL_THROW_ERRNO
                (io_error,
                "Error doing write() in boostfd_request::serve()" <<
                " offset=" << offset <<
                " this=" << this <<
                " buffer=" << buffer <<
                " bytes=" << bytes <<
                " type=" << ((type == request::READ) ? "READ" : "WRITE") <<
                " : " << ex.what());
        }
    }
}

const char* boostfd_file::io_type() const
{
    return "boostfd";
}

boostfd_file::boostfd_file(
    const std::string& filename,
    int mode,
    int queue_id, int allocator_id, unsigned int device_id)
    : file(device_id),
      disk_queued_file(queue_id, allocator_id),
      m_mode(mode)
{
    BOOST_IOS::openmode boostfd_mode =
        (mode & RDWR) ? (BOOST_IOS::out | BOOST_IOS::in) :
        (mode & WRONLY) ? (BOOST_IOS::out) :
        (mode & RDONLY) ? (BOOST_IOS::in) :
        BOOST_IOS::in;

#if defined(BOOST_FILESYSTEM_VERSION) && (BOOST_FILESYSTEM_VERSION >= 3)
    const boost::filesystem::path fspath(filename);
#else
    const boost::filesystem::path fspath(filename,
                                         boost::filesystem::native);
#endif

    if (mode & TRUNC)
    {
        if (boost::filesystem::exists(fspath))
        {
            boost::filesystem::remove(fspath);
            boost::filesystem::ofstream f(fspath);
            f.close();
            assert(boost::filesystem::exists(fspath));
        }
    }

    if (mode & CREAT)
    {
        // need to be emulated:
        if (!boost::filesystem::exists(fspath))
        {
            boost::filesystem::ofstream f(fspath);
            f.close();
            assert(boost::filesystem::exists(fspath));
        }
    }

    if (mode & DIRECT)
    {
        // direct mode not supported in Boost
        STXXL_MSG("Warning: open()ing " << filename << " without DIRECT mode, boostfd does not support it.");
    }

    if (mode & REQUIRE_DIRECT)
    {
        // direct mode not supported in Boost
        STXXL_ERRMSG("Error: open()ing " << filename << " with REQUIRE_DIRECT mode, but boostfd does not support it.");
        return;
    }

    if (mode & SYNC)
    {
        // ???
    }

#if (BOOST_VERSION >= 104100)
    m_file_des.open(filename, boostfd_mode);      // also compiles with earlier Boost versions, but differs semantically
#else
    m_file_des.open(filename, boostfd_mode, boostfd_mode);
#endif
}

boostfd_file::~boostfd_file()
{
    scoped_mutex_lock fd_lock(m_fd_mutex);
    m_file_des.close();
}

inline file::offset_type boostfd_file::_size()
{
    return m_file_des.seek(0, BOOST_IOS::end);
}

file::offset_type boostfd_file::size()
{
    scoped_mutex_lock fd_lock(m_fd_mutex);
    return _size();
}

void boostfd_file::set_size(offset_type newsize)
{
    scoped_mutex_lock fd_lock(m_fd_mutex);
    offset_type oldsize = _size();
    m_file_des.seek(newsize, BOOST_IOS::beg);
    m_file_des.seek(0, BOOST_IOS::beg); // not important ?
    STXXL_ASSERT(_size() >= oldsize);
}

void boostfd_file::lock()
{
    // FIXME: is there no locking possible/needed/... for boostfd?
}

STXXL_END_NAMESPACE

#endif  // #if STXXL_HAVE_BOOSTFD_FILE
// vim: et:ts=4:sw=4
