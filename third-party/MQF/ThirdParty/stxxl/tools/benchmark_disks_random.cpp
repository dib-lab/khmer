/***************************************************************************
 *  tools/benchmark_disks_random.cpp
 *
 *  Part of the STXXL. See http://stxxl.sourceforge.net
 *
 *  Copyright (C) 2009 Johannes Singler <singler@ira.uka.de>
 *  Copyright (C) 2009 Andreas Beckmann <beckmann@cs.uni-frankfurt.de>
 *  Copyright (C) 2013 Timo Bingmann <tb@panthema.net>
 *
 *  Distributed under the Boost Software License, Version 1.0.
 *  (See accompanying file LICENSE_1_0.txt or copy at
 *  http://www.boost.org/LICENSE_1_0.txt)
 **************************************************************************/

/*
   example gnuplot command for the output of this program:
   (x-axis: offset in GiB, y-axis: bandwidth in MiB/s)

   plot \
        "disk.log" using ($2/1024):($7) w l title "read", \
        "disk.log" using ($2/1024):($4)  w l title "write"
 */

#include <iomanip>
#include <vector>
#include <ctime>

#include <stxxl/io>
#include <stxxl/mng>
#include <stxxl/cmdline>

using stxxl::request_ptr;
using stxxl::timestamp;

#define KiB (1024)
#define MiB (1024 * 1024)

struct print_number
{
    int n;

    print_number(int n) : n(n) { }

    void operator () (stxxl::request_ptr)
    {
        //std::cout << n << " " << std::flush;
    }
};

template <unsigned BlockSize, typename AllocStrategy>
void run_test(stxxl::int64 span, stxxl::int64 worksize, bool do_init, bool do_read, bool do_write)
{
    const unsigned raw_block_size = BlockSize;

    typedef stxxl::typed_block<raw_block_size, unsigned> block_type;
    typedef stxxl::BID<raw_block_size> BID_type;

    stxxl::unsigned_type num_blocks =
        (stxxl::unsigned_type)stxxl::div_ceil(worksize, raw_block_size);
    stxxl::unsigned_type num_blocks_in_span =
        (stxxl::unsigned_type)stxxl::div_ceil(span, raw_block_size);

    num_blocks = stxxl::STXXL_MIN(num_blocks, num_blocks_in_span);
    if (num_blocks == 0) num_blocks = num_blocks_in_span;

    worksize = num_blocks * raw_block_size;

    block_type* buffer = new block_type;
    request_ptr* reqs = new request_ptr[num_blocks_in_span];
    std::vector<BID_type> blocks;

    //touch data, so it is actually allocated
    for (unsigned i = 0; i < block_type::size; ++i)
        (*buffer)[i] = i;

    try {
        AllocStrategy alloc;

        blocks.resize(num_blocks_in_span);
        stxxl::block_manager::get_instance()->new_blocks(alloc, blocks.begin(), blocks.end());

        std::cout << "# Span size: "
                  << stxxl::add_IEC_binary_multiplier(span, "B") << " ("
                  << num_blocks_in_span << " blocks of "
                  << stxxl::add_IEC_binary_multiplier(raw_block_size, "B") << ")" << std::endl;

        std::cout << "# Work size: "
                  << stxxl::add_IEC_binary_multiplier(worksize, "B") << " ("
                  << num_blocks << " blocks of "
                  << stxxl::add_IEC_binary_multiplier(raw_block_size, "B") << ")" << std::endl;

        double begin, end, elapsed;

        if (do_init)
        {
            begin = timestamp();
            std::cout << "First fill up space by writing sequentially..." << std::endl;
            for (unsigned j = 0; j < num_blocks_in_span; j++)
                reqs[j] = buffer->write(blocks[j]);
            wait_all(reqs, num_blocks_in_span);
            end = timestamp();
            elapsed = end - begin;
            std::cout << "Written "
                      << std::setw(12) << num_blocks_in_span << " blocks in " << std::fixed << std::setw(9) << std::setprecision(2) << elapsed << " seconds: "
                      << std::setw(9) << std::setprecision(1) << (double(num_blocks_in_span) / elapsed) << " blocks/s "
                      << std::setw(7) << std::setprecision(1) << (double(num_blocks_in_span * raw_block_size) / MiB / elapsed) << " MiB/s write " << std::endl;
        }

        std::cout << "Random block access..." << std::endl;

        srand((unsigned int)time(NULL));
        std::random_shuffle(blocks.begin(), blocks.end());

        begin = timestamp();
        if (do_read)
        {
            for (unsigned j = 0; j < num_blocks; j++)
                reqs[j] = buffer->read(blocks[j], print_number(j));
            wait_all(reqs, num_blocks);

            end = timestamp();
            elapsed = end - begin;
            std::cout << "Read    " << num_blocks << " blocks in " << std::fixed << std::setw(5) << std::setprecision(2) << elapsed << " seconds: "
                      << std::setw(5) << std::setprecision(1) << (double(num_blocks) / elapsed) << " blocks/s "
                      << std::setw(5) << std::setprecision(1) << (double(num_blocks * raw_block_size) / MiB / elapsed) << " MiB/s read" << std::endl;
        }

        std::random_shuffle(blocks.begin(), blocks.end());

        begin = timestamp();
        if (do_write)
        {
            for (unsigned j = 0; j < num_blocks; j++)
                reqs[j] = buffer->write(blocks[j], print_number(j));
            wait_all(reqs, num_blocks);

            end = timestamp();
            elapsed = end - begin;
            std::cout << "Written " << num_blocks << " blocks in " << std::fixed << std::setw(5) << std::setprecision(2) << elapsed << " seconds: "
                      << std::setw(5) << std::setprecision(1) << (double(num_blocks) / elapsed) << " blocks/s "
                      << std::setw(5) << std::setprecision(1) << (double(num_blocks * raw_block_size) / MiB / elapsed) << " MiB/s write " << std::endl;
        }
    }
    catch (const std::exception& ex)
    {
        std::cout << std::endl;
        STXXL_ERRMSG(ex.what());
    }

    delete[] reqs;
    delete buffer;

    stxxl::block_manager::get_instance()->delete_blocks(blocks.begin(), blocks.end());
}

template <typename AllocStrategy>
int benchmark_disks_random_alloc(stxxl::uint64 span, stxxl::uint64 block_size, stxxl::uint64 worksize,
                                 const std::string& optirw)
{
    bool do_init = (optirw.find('i') != std::string::npos);
    bool do_read = (optirw.find('r') != std::string::npos);
    bool do_write = (optirw.find('w') != std::string::npos);

#define run(bs) run_test<bs, AllocStrategy>(span, worksize, do_init, do_read, do_write)
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
                  << "Available are only powers of two from 4 KiB to 128 MiB. You must use 'ki' instead of 'k'." << std::endl;
#undef run

    return 0;
}

int benchmark_disks_random(int argc, char* argv[])
{
    // parse command line

    stxxl::cmdline_parser cp;

    stxxl::uint64 span, block_size = 8 * MiB, worksize = 0;
    std::string optirw = "irw", allocstr;

    cp.add_param_bytes(
        "span", span,
        "Span of external memory to write/read to (e.g. 10GiB).");
    cp.add_opt_param_bytes(
        "block_size", block_size,
        "Size of blocks to randomly write/read (default: 8MiB).");
    cp.add_opt_param_bytes(
        "size", worksize,
        "Amount of data to operate on (e.g. 2GiB), default: whole span.");
    cp.add_opt_param_string(
        "i|r|w", optirw,
        "Operations: [i]nitialize, [r]ead, and/or [w]rite (default: all).");
    cp.add_opt_param_string(
        "alloc", allocstr,
        "Block allocation strategy: RC, SR, FR, striping (default: RC).");

    cp.set_description(
        "This program will benchmark _random_ block access on the disks "
        "configured by the standard .stxxl disk configuration files mechanism. "
        "Available block sizes are power of two from 4 KiB to 128 MiB. "
        "A set of three operations can be performed: sequential initialization, "
        "random reading and random writing."
        );

    if (!cp.process(argc, argv))
        return -1;

#define run_alloc(alloc) benchmark_disks_random_alloc<alloc>(span, block_size, worksize, optirw)
    if (allocstr.size())
    {
        if (allocstr == "RC")
            return run_alloc(stxxl::RC);
        if (allocstr == "SR")
            return run_alloc(stxxl::SR);
        if (allocstr == "FR")
            return run_alloc(stxxl::FR);
        if (allocstr == "striping")
            return run_alloc(stxxl::striping);

        std::cout << "Unknown allocation strategy '" << allocstr << "'" << std::endl;
        cp.print_usage();
        return -1;
    }

    return run_alloc(STXXL_DEFAULT_ALLOC_STRATEGY);
#undef run_alloc
}

// vim: et:ts=4:sw=4
