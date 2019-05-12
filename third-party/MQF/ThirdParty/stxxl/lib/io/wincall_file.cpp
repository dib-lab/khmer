/***************************************************************************
 *  lib/io/wincall_file.cpp
 *
 *  Part of the STXXL. See http://stxxl.sourceforge.net
 *
 *  Copyright (C) 2005-2006 Roman Dementiev <dementiev@ira.uka.de>
 *  Copyright (C) 2008-2010 Andreas Beckmann <beckmann@cs.uni-frankfurt.de>
 *
 *  Distributed under the Boost Software License, Version 1.0.
 *  (See accompanying file LICENSE_1_0.txt or copy at
 *  http://www.boost.org/LICENSE_1_0.txt)
 **************************************************************************/

#include <stxxl/bits/io/wincall_file.h>

#if STXXL_HAVE_WINCALL_FILE

#include <stxxl/bits/io/iostats.h>
#include <stxxl/bits/common/error_handling.h>

#ifndef NOMINMAX
  #define NOMINMAX
#endif
#include <windows.h>

STXXL_BEGIN_NAMESPACE

void wincall_file::serve(void* buffer, offset_type offset, size_type bytes,
                         request::request_type type)
{
    scoped_mutex_lock fd_lock(fd_mutex);

    if (bytes > 32 * 1024 * 1024) {
        STXXL_ERRMSG("Using a block size larger than 32 MiB may not work with the " << io_type() << " filetype");
    }

    HANDLE handle = file_des;
    LARGE_INTEGER desired_pos;
    desired_pos.QuadPart = offset;
    if (!SetFilePointerEx(handle, desired_pos, NULL, FILE_BEGIN))
    {
        STXXL_THROW_WIN_LASTERROR(io_error,
                                  "SetFilePointerEx in wincall_request::serve()" <<
                                  " offset=" << offset <<
                                  " this=" << this <<
                                  " buffer=" << buffer <<
                                  " bytes=" << bytes <<
                                  " type=" << ((type == request::READ) ? "READ" : "WRITE"));
    }
    else
    {
        stats::scoped_read_write_timer read_write_timer(bytes, type == request::WRITE);

        if (type == request::READ)
        {
            DWORD NumberOfBytesRead = 0;
            assert(bytes <= std::numeric_limits<DWORD>::max());
            if (!ReadFile(handle, buffer, (DWORD)bytes, &NumberOfBytesRead, NULL))
            {
                STXXL_THROW_WIN_LASTERROR(io_error,
                                          "ReadFile" <<
                                          " this=" << this <<
                                          " offset=" << offset <<
                                          " buffer=" << buffer <<
                                          " bytes=" << bytes <<
                                          " type=" << ((type == request::READ) ? "READ" : "WRITE") <<
                                          " NumberOfBytesRead= " << NumberOfBytesRead);
            }
            else if (NumberOfBytesRead != bytes) {
                STXXL_THROW_WIN_LASTERROR(io_error, " partial read: missing " << (bytes - NumberOfBytesRead) << " out of " << bytes << " bytes");
            }
        }
        else
        {
            DWORD NumberOfBytesWritten = 0;
            assert(bytes <= std::numeric_limits<DWORD>::max());
            if (!WriteFile(handle, buffer, (DWORD)bytes, &NumberOfBytesWritten, NULL))
            {
                STXXL_THROW_WIN_LASTERROR(io_error,
                                          "WriteFile" <<
                                          " this=" << this <<
                                          " offset=" << offset <<
                                          " buffer=" << buffer <<
                                          " bytes=" << bytes <<
                                          " type=" << ((type == request::READ) ? "READ" : "WRITE") <<
                                          " NumberOfBytesWritten= " << NumberOfBytesWritten);
            }
            else if (NumberOfBytesWritten != bytes) {
                STXXL_THROW_WIN_LASTERROR(io_error, " partial write: missing " << (bytes - NumberOfBytesWritten) << " out of " << bytes << " bytes");
            }
        }
    }
}

const char* wincall_file::io_type() const
{
    return "wincall";
}

STXXL_END_NAMESPACE

#endif  // #if STXXL_HAVE_WINCALL_FILE
// vim: et:ts=4:sw=4
