/***************************************************************************
 *  tools/extras/iobench_scatter_in_place.cpp
 *
 *  Part of the STXXL. See http://stxxl.sourceforge.net
 *
 *  Copyright (C) 2009 Andreas Beckmann <beckmann@cs.uni-frankfurt.de>
 *
 *  Distributed under the Boost Software License, Version 1.0.
 *  (See accompanying file LICENSE_1_0.txt or copy at
 *  http://www.boost.org/LICENSE_1_0.txt)
 **************************************************************************/

#include <iomanip>
#include <vector>
#include <cstdio>

#include <stxxl/io>
#include <stxxl/aligned_alloc>
#include <stxxl/timer>

using stxxl::request_ptr;
using stxxl::file;
using stxxl::timer;
using stxxl::uint64;
using stxxl::unsigned_type;

#ifndef BLOCK_ALIGN
 #define BLOCK_ALIGN  4096
#endif

#define MB (1024 * 1024)
#define GB (1024 * 1024 * 1024)

void usage(const char* argv0)
{
    std::cout << "Usage: " << argv0 << " num_blocks blocks_per_round block_size file" << std::endl;
    std::cout << "    'block_size' in bytes" << std::endl;
    std::cout << "    'file' is split into 'num_blocks' files of size 'block_size'," << std::endl;
    std::cout << "    reading chunks of 'blocks_per_round' blocks starting from end-of-file" << std::endl;
    std::cout << "    and truncating the input file after each chunk was read," << std::endl;
    std::cout << "    before writing the chunk to new files" << std::endl;
    exit(-1);
}

// returns throughput in MiB/s
inline double throughput(stxxl::uint64 bytes, double seconds)
{
    if (seconds == 0.0)
        return 0.0;
    return (double)bytes / (1024 * 1024) / seconds;
}

int main(int argc, char* argv[])
{
    if (argc < 5)
        usage(argv[0]);

    unsigned_type num_blocks = (unsigned_type)stxxl::atouint64(argv[1]);
    unsigned_type blocks_per_round = (unsigned_type)stxxl::atouint64(argv[2]);
    unsigned_type block_size = (unsigned_type)stxxl::atouint64(argv[3]);
    const char* filebase = argv[4];

    unsigned_type num_rounds = stxxl::div_ceil(num_blocks, blocks_per_round);

    std::cout << "# Splitting '" << filebase << "' into "
              << num_rounds * blocks_per_round << " blocks of size "
              << block_size << ", reading chunks of "
              << blocks_per_round << " blocks" << std::endl;

    char* buffer = (char*)stxxl::aligned_alloc<BLOCK_ALIGN>(block_size * blocks_per_round);
    double totaltimeread = 0, totaltimewrite = 0;
    stxxl::int64 totalsizeread = 0, totalsizewrite = 0;
    double totaltimereadchunk = 0.0, totaltimewritechunk = 0.0;
    stxxl::int64 totalsizereadchunk = 0, totalsizewritechunk = 0;

    typedef stxxl::syscall_file file_type;

    file_type input_file(filebase, file::RDWR | file::DIRECT, 0);

    timer t_total(true);
    try {
        for (stxxl::unsigned_type r = num_rounds; r-- > 0; )
        {
            // read a chunk of blocks_per_round blocks
            timer t_read(true);
            for (stxxl::unsigned_type i = blocks_per_round; i-- > 0; )
            {
                const uint64 offset = (r * blocks_per_round + i) * block_size;
                timer t_op(true);
                // read a block
                {
                    input_file.aread(buffer + i * block_size, offset, block_size)->wait();
                }
                t_op.stop();
                totalsizeread += block_size;
                totaltimeread += t_op.seconds();
                if (blocks_per_round > 1) {
                    std::cout << "Offset         " << std::setw(8) << offset / MB << " MiB: " << std::fixed;
                    std::cout << std::setw(8) << std::setprecision(3) << throughput(block_size, t_op.seconds()) << " MiB/s read";
                    std::cout << std::endl;
                }
            }

            // truncate
            input_file.set_size(r * blocks_per_round * block_size);

            t_read.stop();
            totalsizereadchunk += blocks_per_round * block_size;
            totaltimereadchunk += t_read.seconds();

            // write the chunk
            timer t_write(true);
            for (stxxl::unsigned_type i = blocks_per_round; i-- > 0; )
            {
                const uint64 offset = (r * blocks_per_round + i) * block_size;
                timer t_op(true);
                // write a block
                {
                    char cfn[4096]; // PATH_MAX
                    snprintf(cfn, sizeof(cfn), "%s_%012llX", filebase, offset);
                    file_type chunk_file(cfn, file::CREAT | file::RDWR | file::DIRECT, 0);
                    chunk_file.awrite(buffer + i * block_size, 0, block_size)->wait();
                }
                t_op.stop();
                totalsizewrite += block_size;
                totaltimewrite += t_op.seconds();
                if (blocks_per_round > 1) {
                    std::cout << "Offset         " << std::setw(8) << offset / MB << " MiB: " << std::fixed;
                    std::cout << std::setw(8) << std::setprecision(3) << "" << "             ";
                    std::cout << std::setw(8) << std::setprecision(3) << throughput(block_size, t_op.seconds()) << " MiB/s write";
                    std::cout << std::endl;
                }
            }
            t_write.stop();
            totalsizewritechunk += blocks_per_round * block_size;
            totaltimewritechunk += t_write.seconds();

            const uint64 offset = r * blocks_per_round * block_size;
            std::cout << "Input offset   " << std::setw(8) << offset / MB << " MiB: " << std::fixed;
            std::cout << std::setw(8) << std::setprecision(3) << throughput(block_size * blocks_per_round, t_read.seconds()) << " MiB/s read, ";
            std::cout << std::setw(8) << std::setprecision(3) << throughput(block_size * blocks_per_round, t_write.seconds()) << " MiB/s write";
            std::cout << std::endl;
        }
    }
    catch (const std::exception& ex)
    {
        std::cout << std::endl;
        STXXL_ERRMSG(ex.what());
    }
    t_total.stop();

    const int ndisks = 1;

    std::cout << "=============================================================================================" << std::endl;
    std::cout << "# Average over " << std::setw(8) << stxxl::STXXL_MAX(totalsizewrite, totalsizeread) / MB << " MiB: ";
    std::cout << std::setw(8) << std::setprecision(3) << (throughput(totalsizeread, totaltimeread)) << " MiB/s read, ";
    std::cout << std::setw(8) << std::setprecision(3) << (throughput(totalsizewrite, totaltimewrite)) << " MiB/s write" << std::endl;
    if (totaltimeread != 0.0)
        std::cout << "# Read time    " << std::setw(8) << std::setprecision(3) << totaltimeread << " s" << std::endl;
    if (totaltimereadchunk != 0.0)
        std::cout << "# ChRd/trnk ti " << std::setw(8) << std::setprecision(3) << totaltimereadchunk << " s" << std::endl;
    if (totaltimewrite != 0.0)
        std::cout << "# Write time   " << std::setw(8) << std::setprecision(3) << totaltimewrite << " s" << std::endl;
    if (totaltimewritechunk != 0.0)
        std::cout << "# ChWrite time " << std::setw(8) << std::setprecision(3) << totaltimewritechunk << " s" << std::endl;
    std::cout << "# Non-I/O time " << std::setw(8) << std::setprecision(3) << (t_total.seconds() - totaltimewrite - totaltimeread) << " s, average throughput "
              << std::setw(8) << std::setprecision(3) << (throughput(totalsizewrite + totalsizeread, t_total.seconds() - totaltimewrite - totaltimeread) * ndisks) << " MiB/s"
              << std::endl;
    std::cout << "# Total time   " << std::setw(8) << std::setprecision(3) << t_total.seconds() << " s, average throughput "
              << std::setw(8) << std::setprecision(3) << (throughput(totalsizewrite + totalsizeread, t_total.seconds()) * ndisks) << " MiB/s"
              << std::endl;

    stxxl::aligned_dealloc<BLOCK_ALIGN>(buffer);

    return 0;
}

// vim: et:ts=4:sw=4
