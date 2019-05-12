/***************************************************************************
 *  tools/benchmark_disks.cpp
 *
 *  Part of the STXXL. See http://stxxl.sourceforge.net
 *
 *  Copyright (C) 2009 Johannes Singler <singler@ira.uka.de>
 *  Copyright (C) 2013 Timo Bingmann <tb@panthema.net>
 *
 *  Distributed under the Boost Software License, Version 1.0.
 *  (See accompanying file LICENSE_1_0.txt or copy at
 *  http://www.boost.org/LICENSE_1_0.txt)
 **************************************************************************/

/*
  This programm will benchmark the disks configured via .stxxl disk
  configuration files. The block manager is used to read and write blocks using
  the different allocation strategies.
*/

/*
   example gnuplot command for the output of this program:
   (x-axis: offset in GiB, y-axis: bandwidth in MiB/s)

   plot \
        "disk.log" using ($2/1024):($7) w l title "read", \
        "disk.log" using ($2/1024):($4)  w l title "write"
 */

#include <iomanip>
#include <vector>

#include <stxxl/io>
#include <stxxl/mng>
#include <stxxl/bits/common/cmdline.h>

using stxxl::timestamp;
using stxxl::unsigned_type;
using stxxl::uint64;

#ifdef BLOCK_ALIGN
 #undef BLOCK_ALIGN
#endif

#define BLOCK_ALIGN  4096

#define POLL_DELAY 1000

#define CHECK_AFTER_READ 0

#define KiB (1024)
#define MiB (1024 * 1024)

template <unsigned_type RawBlockSize, typename AllocStrategy>
int benchmark_disks_blocksize_alloc(uint64 length, uint64 start_offset, uint64 batch_size,
                                    std::string optrw)
{
    uint64 endpos = start_offset + length;

    if (length == 0)
        endpos = std::numeric_limits<uint64>::max();

    bool do_read = (optrw.find('r') != std::string::npos);
    bool do_write = (optrw.find('w') != std::string::npos);

    // initialize disk configuration
    stxxl::block_manager::get_instance();

    // construct block type

    const unsigned_type raw_block_size = RawBlockSize;
    const unsigned_type block_size = raw_block_size / sizeof(int);

    typedef stxxl::typed_block<raw_block_size, unsigned> block_type;
    typedef stxxl::BID<raw_block_size> BID_type;

    if (batch_size == 0)
        batch_size = stxxl::config::get_instance()->disks_number();

    // calculate total bytes processed in a batch
    batch_size = raw_block_size * batch_size;

    unsigned_type num_blocks_per_batch = (unsigned_type)stxxl::div_ceil(batch_size, raw_block_size);
    batch_size = num_blocks_per_batch * raw_block_size;

    block_type* buffer = new block_type[num_blocks_per_batch];
    stxxl::request_ptr* reqs = new stxxl::request_ptr[num_blocks_per_batch];
    std::vector<BID_type> blocks;
    double totaltimeread = 0, totaltimewrite = 0;
    uint64 totalsizeread = 0, totalsizewrite = 0;

    std::cout << "# Batch size: "
              << stxxl::add_IEC_binary_multiplier(batch_size, "B") << " ("
              << num_blocks_per_batch << " blocks of "
              << stxxl::add_IEC_binary_multiplier(raw_block_size, "B") << ")"
              << " using " << AllocStrategy().name()
              << std::endl;

    // touch data, so it is actually allcoated
    for (unsigned j = 0; j < num_blocks_per_batch; ++j)
        for (unsigned i = 0; i < block_size; ++i)
            buffer[j][i] = (unsigned)(j * block_size + i);

    try {
        AllocStrategy alloc;
        uint64 current_batch_size;

        for (uint64 offset = 0; offset < endpos; offset += current_batch_size)
        {
            current_batch_size = std::min<uint64>(batch_size, endpos - offset);
#if CHECK_AFTER_READ
            const uint64 current_batch_size_int = current_batch_size / sizeof(int);
#endif
            const unsigned_type current_num_blocks_per_batch = (unsigned_type)stxxl::div_ceil(current_batch_size, raw_block_size);

            unsigned_type num_total_blocks = blocks.size();
            blocks.resize(num_total_blocks + current_num_blocks_per_batch);
            stxxl::block_manager::get_instance()->new_blocks(alloc, blocks.begin() + num_total_blocks, blocks.end());

            if (offset < start_offset)
                continue;

            std::cout << "Offset    " << std::setw(7) << offset / MiB << " MiB: " << std::fixed;

            double begin = timestamp(), end, elapsed;

            if (do_write)
            {
                for (unsigned j = 0; j < current_num_blocks_per_batch; j++)
                    reqs[j] = buffer[j].write(blocks[num_total_blocks + j]);

                wait_all(reqs, current_num_blocks_per_batch);

                end = timestamp();
                elapsed = end - begin;
                totalsizewrite += current_batch_size;
                totaltimewrite += elapsed;
            }
            else
                elapsed = 0.0;

            std::cout << std::setw(5) << std::setprecision(1) << (double(current_batch_size) / MiB / elapsed) << " MiB/s write, ";

            begin = timestamp();

            if (do_read)
            {
                for (unsigned j = 0; j < current_num_blocks_per_batch; j++)
                    reqs[j] = buffer[j].read(blocks[num_total_blocks + j]);

                wait_all(reqs, current_num_blocks_per_batch);

                end = timestamp();
                elapsed = end - begin;
                totalsizeread += current_batch_size;
                totaltimeread += elapsed;
            }
            else
                elapsed = 0.0;

            std::cout << std::setw(5) << std::setprecision(1) << (double(current_batch_size) / MiB / elapsed) << " MiB/s read" << std::endl;

#if CHECK_AFTER_READ
            for (unsigned j = 0; j < current_num_blocks_per_batch; j++)
            {
                for (unsigned i = 0; i < block_size; i++)
                {
                    if (buffer[j][i] != j * block_size + i)
                    {
                        int ibuf = i / current_batch_size_int;
                        int pos = i % current_batch_size_int;

                        std::cout << "Error on disk " << ibuf << " position " << std::hex << std::setw(8) << offset + pos * sizeof(int)
                                  << "  got: " << std::hex << std::setw(8) << buffer[j][i] << " wanted: " << std::hex << std::setw(8) << (j * block_size + i)
                                  << std::dec << std::endl;

                        i = (ibuf + 1) * current_batch_size_int; // jump to next
                    }
                }
            }
#endif
        }
    }
    catch (const std::exception& ex)
    {
        std::cout << std::endl;
        STXXL_ERRMSG(ex.what());
    }

    std::cout << "=============================================================================================" << std::endl;
    std::cout << "# Average over " << std::setw(7) << totalsizewrite / MiB << " MiB: ";
    std::cout << std::setw(5) << std::setprecision(1) << (double(totalsizewrite) / MiB / totaltimewrite) << " MiB/s write, ";
    std::cout << std::setw(5) << std::setprecision(1) << (double(totalsizeread) / MiB / totaltimeread) << " MiB/s read" << std::endl;

    delete[] reqs;
    delete[] buffer;

    return 0;
}

template <typename AllocStrategy>
int benchmark_disks_alloc(uint64 length, uint64 offset, uint64 batch_size,
                          unsigned_type block_size, std::string optrw)
{
#define run(bs) benchmark_disks_blocksize_alloc<bs, AllocStrategy>(length, offset, batch_size, optrw)
    if (block_size == 4 * KiB)
        run(4 * KiB);
    else if (block_size == 8 * KiB)
        run(8 * KiB);
    else if (block_size == 16 * KiB)
        run(16 * KiB);
    else if (block_size == 32 * KiB)
        run(32 * KiB);
    else if (block_size == 64 * KiB)
        run(64 * KiB);
    else if (block_size == 128 * KiB)
        run(128 * KiB);
    else if (block_size == 256 * KiB)
        run(256 * KiB);
    else if (block_size == 512 * KiB)
        run(512 * KiB);
    else if (block_size == 1 * MiB)
        run(1 * MiB);
    else if (block_size == 2 * MiB)
        run(2 * MiB);
    else if (block_size == 4 * MiB)
        run(4 * MiB);
    else if (block_size == 8 * MiB)
        run(8 * MiB);
    else if (block_size == 16 * MiB)
        run(16 * MiB);
    else if (block_size == 32 * MiB)
        run(32 * MiB);
    else if (block_size == 64 * MiB)
        run(64 * MiB);
    else if (block_size == 128 * MiB)
        run(128 * MiB);
    else
        std::cerr << "Unsupported block_size " << block_size << "." << std::endl
                  << "Available are only powers of two from 4 KiB to 128 MiB. "
                  << "You must use 'ki' instead of 'k'." << std::endl;
#undef run

    return 0;
}

int benchmark_disks(int argc, char* argv[])
{
    // parse command line

    stxxl::cmdline_parser cp;

    uint64 length = 0, offset = 0;
    unsigned int batch_size = 0;
    unsigned_type block_size = 8 * MiB;
    std::string optrw = "rw", allocstr;

    cp.add_param_bytes("size", length,
                       "Amount of data to write/read from disks (e.g. 10GiB)");
    cp.add_opt_param_string(
        "r|w", optrw,
        "Only read or write blocks (default: both write and read)");
    cp.add_opt_param_string(
        "alloc", allocstr,
        "Block allocation strategy: RC, SR, FR, striping. (default: RC)");

    cp.add_uint('b', "batch", batch_size,
                "Number of blocks written/read in one batch (default: D * B)");
    cp.add_bytes('B', "block_size", block_size,
                 "Size of blocks written in one syscall. (default: B = 8MiB)");
    cp.add_bytes('o', "offset", offset,
                 "Starting offset of operation range. (default: 0)");

    cp.set_description(
        "This program will benchmark the disks configured by the standard "
        ".stxxl disk configuration files mechanism. Blocks of 8 MiB are "
        "written and/or read in sequence using the block manager. The batch "
        "size describes how many blocks are written/read in one batch. The "
        "are taken from block_manager using given the specified allocation "
        "strategy. If size == 0, then writing/reading operation are done "
        "until an error occurs. ");

    if (!cp.process(argc, argv))
        return -1;

    if (allocstr.size())
    {
        if (allocstr == "RC")
            return benchmark_disks_alloc<stxxl::RC>(
                length, offset, batch_size, block_size, optrw);
        if (allocstr == "SR")
            return benchmark_disks_alloc<stxxl::SR>(
                length, offset, batch_size, block_size, optrw);
        if (allocstr == "FR")
            return benchmark_disks_alloc<stxxl::FR>(
                length, offset, batch_size, block_size, optrw);
        if (allocstr == "striping")
            return benchmark_disks_alloc<stxxl::striping>(
                length, offset, batch_size, block_size, optrw);

        std::cout << "Unknown allocation strategy '" << allocstr << "'" << std::endl;
        cp.print_usage();
        return -1;
    }

    return benchmark_disks_alloc<STXXL_DEFAULT_ALLOC_STRATEGY>(
        length, offset, batch_size, block_size, optrw);
}
