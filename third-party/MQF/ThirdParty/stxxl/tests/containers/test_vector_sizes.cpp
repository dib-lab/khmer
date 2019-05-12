/***************************************************************************
 *  tests/containers/test_vector_sizes.cpp
 *
 *  Part of the STXXL. See http://stxxl.sourceforge.net
 *
 *  Copyright (C) 2010 Andreas Beckmann <beckmann@cs.uni-frankfurt.de>
 *
 *  Distributed under the Boost Software License, Version 1.0.
 *  (See accompanying file LICENSE_1_0.txt or copy at
 *  http://www.boost.org/LICENSE_1_0.txt)
 **************************************************************************/

#include <stxxl/io>
#include <stxxl/vector>

typedef int my_type;
typedef stxxl::VECTOR_GENERATOR<my_type>::result vector_type;
typedef vector_type::block_type block_type;

void test_write(const char* fn, const char* ft, stxxl::unsigned_type sz, my_type ofs)
{
    stxxl::file* f = stxxl::create_file(ft, fn, stxxl::file::CREAT | stxxl::file::DIRECT | stxxl::file::RDWR);
    {
        vector_type v(f);
        v.resize(sz);
        STXXL_MSG("writing " << v.size() << " elements");
        for (stxxl::unsigned_type i = 0; i < v.size(); ++i)
            v[i] = ofs + (int)i;
    }
    delete f;
}

template <typename Vector>
void test_rdwr(const char* fn, const char* ft, stxxl::unsigned_type sz, my_type ofs)
{
    stxxl::file* f = stxxl::create_file(ft, fn, stxxl::file::DIRECT | stxxl::file::RDWR);
    {
        Vector v(f);
        STXXL_MSG("reading " << v.size() << " elements (RDWR)");
        STXXL_CHECK(v.size() == sz);
        for (stxxl::unsigned_type i = 0; i < v.size(); ++i)
            STXXL_CHECK(v[i] == ofs + my_type(i));
    }
    delete f;
}

template <typename Vector>
void test_rdonly(const char* fn, const char* ft, stxxl::unsigned_type sz, my_type ofs)
{
    stxxl::file* f = stxxl::create_file(ft, fn, stxxl::file::DIRECT | stxxl::file::RDONLY);
    {
        Vector v(f);
        STXXL_MSG("reading " << v.size() << " elements (RDONLY)");
        STXXL_CHECK(v.size() == sz);
        for (stxxl::unsigned_type i = 0; i < v.size(); ++i)
            STXXL_CHECK(v[i] == ofs + my_type(i));
    }
    delete f;
}

void test(const char* fn, const char* ft, stxxl::unsigned_type sz, my_type ofs)
{
    test_write(fn, ft, sz, ofs);
    test_rdwr<const vector_type>(fn, ft, sz, ofs);
    test_rdwr<vector_type>(fn, ft, sz, ofs);

    // 2013-tb: there is a bug with read-only vectors on mmap backed files:
    // copying from mmapped area will fail for invalid ranges at the end,
    // whereas a usual read() will just stop short at the end. The read-only
    // vector however will always read the last block in full, thus causing a
    // segfault with mmap files. FIXME
    if (strcmp(ft, "mmap") == 0) return;

    test_rdonly<const vector_type>(fn, ft, sz, ofs);
    //-tb: vector always writes data! FIXME
    //test_rdonly<vector_type>(fn, ft, sz, ofs);
}

int main(int argc, char** argv)
{
    if (argc < 2)
    {
        std::cout << "Usage: " << argv[0] << " file [filetype]" << std::endl;
        return -1;
    }

    stxxl::config::get_instance();

    const char* fn = argv[1];
    const char* ft = (argc >= 3) ? argv[2] : "syscall";

    stxxl::unsigned_type start_elements = 42 * block_type::size;

    STXXL_MSG("using " << ft << " file");

    // multiple of block size
    STXXL_MSG("running test with " << start_elements << " items");
    test(fn, ft, start_elements, 100000000);

    // multiple of page size, but not block size
    STXXL_MSG("running test with " << start_elements << " + 4096 items");
    test(fn, ft, start_elements + 4096, 200000000);

    // multiple of neither block size nor page size
    STXXL_MSG("running test with " << start_elements << " + 4096 + 23 items");
    test(fn, ft, start_elements + 4096 + 23, 300000000);

    // truncate 1 byte
    {
        stxxl::syscall_file f(fn, stxxl::file::DIRECT | stxxl::file::RDWR);
        STXXL_MSG("file size is " << f.size() << " bytes");
        f.set_size(f.size() - 1);
        STXXL_MSG("truncated to " << f.size() << " bytes");
    }

    // will truncate after the last complete element
    test_rdwr<vector_type>(fn, ft, start_elements + 4096 + 23 - 1, 300000000);

    // truncate 1 more byte
    {
        stxxl::syscall_file f(fn, stxxl::file::DIRECT | stxxl::file::RDWR);
        STXXL_MSG("file size is " << f.size() << " bytes");
        f.set_size(f.size() - 1);
        STXXL_MSG("truncated to " << f.size() << " bytes");
    }

    // will not truncate
    //-tb: vector already writes data! TODO
    //test_rdonly<vector_type>(fn, ft, start_elements + 4096 + 23 - 2, 300000000);

    // check final size
    {
        stxxl::syscall_file f(fn, stxxl::file::DIRECT | stxxl::file::RDWR);
        STXXL_MSG("file size is " << f.size() << " bytes");
        STXXL_CHECK(f.size() == (start_elements + 4096 + 23 - 1) * sizeof(my_type) - 1);
    }

    {
        stxxl::syscall_file f(fn, stxxl::file::DIRECT | stxxl::file::RDWR);
        f.close_remove();
    }
}
// vim: et:ts=4:sw=4
