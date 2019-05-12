/***************************************************************************
 *  include/stxxl/bits/io/ufs_file_base.h
 *
 *  UNIX file system file base
 *
 *  Part of the STXXL. See http://stxxl.sourceforge.net
 *
 *  Copyright (C) 2002 Roman Dementiev <dementiev@mpi-sb.mpg.de>
 *  Copyright (C) 2008 Andreas Beckmann <beckmann@cs.uni-frankfurt.de>
 *  Copyright (C) 2009 Johannes Singler <singler@ira.uka.de>
 *  Copyright (C) 2013 Timo Bingmann <tb@panthema.net>
 *
 *  Distributed under the Boost Software License, Version 1.0.
 *  (See accompanying file LICENSE_1_0.txt or copy at
 *  http://www.boost.org/LICENSE_1_0.txt)
 **************************************************************************/

#ifndef STXXL_IO_UFS_FILE_BASE_HEADER
#define STXXL_IO_UFS_FILE_BASE_HEADER

#include <stxxl/bits/common/mutex.h>
#include <stxxl/bits/io/file.h>
#include <stxxl/bits/namespace.h>

#include <string>

STXXL_BEGIN_NAMESPACE

//! \addtogroup fileimpl
//! \{

//! Base for UNIX file system implementations.
class ufs_file_base : public virtual file
{
protected:
    mutex fd_mutex;        // sequentialize function calls involving file_des
    int file_des;          // file descriptor
    int m_mode;            // open mode
    const std::string filename;
    bool m_is_device;      //!< is special device node
    ufs_file_base(const std::string& filename, int mode);
    void _after_open();
    offset_type _size();
    void _set_size(offset_type newsize);
    void close();

public:
    ~ufs_file_base();
    offset_type size();
    void set_size(offset_type newsize);
    void lock();
    const char * io_type() const;
    void close_remove();
    //! unlink file without closing it.
    void unlink();
    //! return true if file is special device node
    bool is_device() const;
};

//! \}

STXXL_END_NAMESPACE

#endif // !STXXL_IO_UFS_FILE_BASE_HEADER
