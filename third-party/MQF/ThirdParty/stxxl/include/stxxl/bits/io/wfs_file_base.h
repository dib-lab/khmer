/***************************************************************************
 *  include/stxxl/bits/io/wfs_file_base.h
 *
 *  Windows file system file base
 *
 *  Part of the STXXL. See http://stxxl.sourceforge.net
 *
 *  Copyright (C) 2005 Roman Dementiev <dementiev@ira.uka.de>
 *  Copyright (C) 2008, 2010 Andreas Beckmann <beckmann@cs.uni-frankfurt.de>
 *  Copyright (C) 2009, 2010 Johannes Singler <singler@kit.edu>
 *
 *  Distributed under the Boost Software License, Version 1.0.
 *  (See accompanying file LICENSE_1_0.txt or copy at
 *  http://www.boost.org/LICENSE_1_0.txt)
 **************************************************************************/

#ifndef STXXL_IO_WFS_FILE_BASE_HEADER
#define STXXL_IO_WFS_FILE_BASE_HEADER

#include <stxxl/bits/config.h>

#if STXXL_WINDOWS

#include <stxxl/bits/io/file.h>
#include <stxxl/bits/io/request.h>

STXXL_BEGIN_NAMESPACE

//! \addtogroup fileimpl
//! \{

//! Base for Windows file system implementations.
class wfs_file_base : public virtual file
{
protected:
    typedef void* HANDLE;

    mutex fd_mutex;        // sequentialize function calls involving file_des
    HANDLE file_des;       // file descriptor
    int mode_;             // open mode
    const std::string filename;
    offset_type bytes_per_sector;
    bool locked;
    wfs_file_base(const std::string& filename, int mode);
    offset_type _size();
    void close();

public:
    ~wfs_file_base();
    offset_type size();
    void set_size(offset_type newsize);
    void lock();
    const char * io_type() const;
    void close_remove();
};

//! \}

STXXL_END_NAMESPACE

#endif // STXXL_WINDOWS

#endif // !STXXL_IO_WFS_FILE_BASE_HEADER
