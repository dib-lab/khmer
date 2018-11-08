/***************************************************************************
 *  tests/io/test_sim_disk.cpp
 *
 *  Part of the STXXL. See http://stxxl.sourceforge.net
 *
 *  Copyright (C) 2002 Roman Dementiev <dementiev@mpi-sb.mpg.de>
 *
 *  Distributed under the Boost Software License, Version 1.0.
 *  (See accompanying file LICENSE_1_0.txt or copy at
 *  http://www.boost.org/LICENSE_1_0.txt)
 **************************************************************************/

// Tests sim_disk file implementation
// must be run on in-memory swap partition !

#include <cmath>
#include <stxxl/io>
#include <stxxl/random>
#include <stxxl/aligned_alloc>

using stxxl::file;
using stxxl::timestamp;

int main()
{
    const stxxl::int64 disk_size = stxxl::int64(1024 * 1024) * 1024 * 40;
    std::cout << sizeof(void*) << std::endl;
    const int block_size = 4 * 1024 * 1024;
    char* buffer = static_cast<char*>(stxxl::aligned_alloc<STXXL_BLOCK_ALIGN>(block_size));
    memset(buffer, 0, block_size);
    const char* paths[2] = { "/tmp/data1", "/tmp/data2" };
    stxxl::sim_disk_file file1(paths[0], file::CREAT | file::RDWR | file::DIRECT, 0);
    file1.set_size(disk_size);

    stxxl::sim_disk_file file2(paths[1], file::CREAT | file::RDWR | file::DIRECT, 1);
    file2.set_size(disk_size);

    unsigned i = 0;

    stxxl::int64 pos = 0;

    stxxl::request_ptr req;

    STXXL_MSG("Estimated time:" << block_size / stxxl::simdisk_geometry::s_average_speed);
    STXXL_MSG("Sequential write");

    for (i = 0; i < 40; i++)
    {
        double begin = timestamp();
        req = file1.awrite(buffer, pos, block_size);
        req->wait();
        double end = timestamp();

        STXXL_MSG("Pos: " << pos << " block_size:" << block_size << " time:" << (end - begin));
        pos += 1024 * 1024 * 1024;
    }

    double sum = 0.;
    double sum2 = 0.;
    STXXL_MSG("Random write");
    const unsigned int times = 80;
    stxxl::random_number<> rnd;
    for (i = 0; i < times; i++)
    {
        pos = (stxxl::int64)rnd(disk_size / block_size) * block_size;
        double begin = timestamp();
        req = file1.awrite(buffer, pos, block_size);
        req->wait();
        double diff = timestamp() - begin;

        sum += diff;
        sum2 += diff * diff;

        STXXL_MSG("Pos: " << pos << " block_size:" << block_size << " time:" << (diff));
    }

    sum = sum / double(times);
    sum2 = sum2 / double(times);
    STXXL_CHECK(sum2 - sum * sum >= 0.0);
    double err = sqrt(sum2 - sum * sum);
    STXXL_MSG("Standard Deviation: " << err << " s, " << 100. * (err / sum) << " %");

    stxxl::aligned_dealloc<STXXL_BLOCK_ALIGN>(buffer);

    file1.close_remove();
    file2.close_remove();

    return 0;
}
