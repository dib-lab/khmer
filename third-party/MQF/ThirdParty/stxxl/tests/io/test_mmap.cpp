/***************************************************************************
 *  tests/io/test_mmap.cpp
 *
 *  Part of the STXXL. See http://stxxl.sourceforge.net
 *
 *  Copyright (C) 2007 Roman Dementiev <dementiev@ira.uka.de>
 *  Copyright (C) 2013 Timo Bingmann <tb@panthema.net>
 *
 *  Distributed under the Boost Software License, Version 1.0.
 *  (See accompanying file LICENSE_1_0.txt or copy at
 *  http://www.boost.org/LICENSE_1_0.txt)
 **************************************************************************/

#include <stxxl.h>

struct my_handler
{
    void operator () (stxxl::request* ptr)
    {
        STXXL_MSG("Request completed: " << ptr);
    }
};

void testIO()
{
    const int size = 1024 * 384;
    char* buffer = static_cast<char*>(stxxl::aligned_alloc<STXXL_BLOCK_ALIGN>(size));
    memset(buffer, 0, size);
#if STXXL_WINDOWS
    const char* paths[2] = { "data1", "data2" };
#else
    const char* paths[2] = { "/var/tmp/data1", "/var/tmp/data2" };
    stxxl::mmap_file file1(paths[0], stxxl::file::CREAT | stxxl::file::RDWR, 0);
    file1.set_size(size * 1024);
#endif

    stxxl::syscall_file file2(paths[1], stxxl::file::CREAT | stxxl::file::RDWR, 1);

    stxxl::request_ptr req[16];
    unsigned i = 0;
    for ( ; i < 16; i++)
        req[i] = file2.awrite(buffer, i * size, size, my_handler());

    stxxl::wait_all(req, 16);

    stxxl::aligned_dealloc<STXXL_BLOCK_ALIGN>(buffer);

#if !STXXL_WINDOWS
    file1.close_remove();
#endif
    file2.close_remove();
}

void testIOException()
{
    stxxl::file::unlink("TestFile");
    // try to open non-existing files
    STXXL_CHECK_THROW(stxxl::mmap_file file1("TestFile", stxxl::file::RDWR, 0), stxxl::io_error);
    STXXL_CHECK_THROW(stxxl::syscall_file file1("TestFile", stxxl::file::RDWR, 0), stxxl::io_error);
}

int main()
{
    testIO();
    testIOException();
}
