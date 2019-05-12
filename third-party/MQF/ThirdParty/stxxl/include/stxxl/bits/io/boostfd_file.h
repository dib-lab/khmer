/***************************************************************************
 *  include/stxxl/bits/io/boostfd_file.h
 *
 *  File implementation based on boost::iostreams::file_decriptor
 *
 *  Part of the STXXL. See http://stxxl.sourceforge.net
 *
 *  Copyright (C) 2006 Roman Dementiev <dementiev@ira.uka.de>
 *  Copyright (C) 2008 Andreas Beckmann <beckmann@cs.uni-frankfurt.de>
 *  Copyright (C) 2009 Johannes Singler <singler@ira.uka.de>
 *  Copyright (C) 2014 Timo Bingmann <tb@panthema.net>
 *
 *  Distributed under the Boost Software License, Version 1.0.
 *  (See accompanying file LICENSE_1_0.txt or copy at
 *  http://www.boost.org/LICENSE_1_0.txt)
 **************************************************************************/

#ifndef STXXL_IO_BOOSTFD_FILE_HEADER
#define STXXL_IO_BOOSTFD_FILE_HEADER

#include <stxxl/bits/config.h>

#ifndef STXXL_HAVE_BOOSTFD_FILE
#if STXXL_BOOST_CONFIG // if boost is available
 #define STXXL_HAVE_BOOSTFD_FILE 1
#else
 #define STXXL_HAVE_BOOSTFD_FILE 0
#endif
#endif

#if STXXL_HAVE_BOOSTFD_FILE

#include <stxxl/bits/io/disk_queued_file.h>
#include <stxxl/bits/io/request.h>

#include <boost/iostreams/device/file_descriptor.hpp>

STXXL_BEGIN_NAMESPACE

//! \addtogroup fileimpl
//! \{

//! Implementation based on boost::iostreams::file_decriptor.
class boostfd_file : public disk_queued_file
{
    typedef boost::iostreams::file_descriptor fd_type;

protected:
    //! sequentialize function calls involving m_file_des
    mutex m_fd_mutex;
    fd_type m_file_des;
    int m_mode;
    offset_type _size();

public:
    boostfd_file(
        const std::string& filename, int mode,
        int queue_id = DEFAULT_QUEUE,
        int allocator_id = NO_ALLOCATOR,
        unsigned int device_id = DEFAULT_DEVICE_ID);
    ~boostfd_file();
    offset_type size();
    void set_size(offset_type newsize);
    void lock();
    void serve(void* buffer, offset_type offset, size_type bytes,
               request::request_type type);
    const char * io_type() const;
};

//! \}

STXXL_END_NAMESPACE

#endif // #if STXXL_HAVE_BOOSTFD_FILE

#endif // !STXXL_IO_BOOSTFD_FILE_HEADER
