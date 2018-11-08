/***************************************************************************
 *  tools/benchmark_files.cpp
 *
 *  Part of the STXXL. See http://stxxl.sourceforge.net
 *
 *  Copyright (C) 2002-2003 Roman Dementiev <dementiev@mpi-sb.mpg.de>
 *  Copyright (C) 2007-2011 Andreas Beckmann <beckmann@cs.uni-frankfurt.de>
 *  Copyright (C) 2009 Johannes Singler <singler@ira.uka.de>
 *  Copyright (C) 2013 Timo Bingmann <tb@panthema.net>
 *
 *  Distributed under the Boost Software License, Version 1.0.
 *  (See accompanying file LICENSE_1_0.txt or copy at
 *  http://www.boost.org/LICENSE_1_0.txt)
 **************************************************************************/

/*
   example gnuplot command for the output of this program:
   (x-axis: file offset in GiB, y-axis: bandwidth in MiB/s)

   plot \
        "file.log" using ($3/1024):($14) w l title "read", \
        "file.log" using ($3/1024):($7)  w l title "write"
 */

#include <iomanip>
#include <cstring>
#include <vector>

#include <stxxl/io>
#include <stxxl/aligned_alloc>
#include <stxxl/timer>
#include <stxxl/bits/version.h>
#include <stxxl/bits/common/cmdline.h>

using stxxl::request_ptr;
using stxxl::file;
using stxxl::timestamp;
using stxxl::unsigned_type;
using stxxl::uint64;

#ifdef BLOCK_ALIGN
 #undef BLOCK_ALIGN
#endif

#define BLOCK_ALIGN  4096

#define POLL_DELAY 1000

#if STXXL_WINDOWS
const char* default_file_type = "wincall";
#else
const char* default_file_type = "syscall";
#endif

#ifdef WATCH_TIMES
void watch_times(request_ptr reqs[], unsigned n, double* out)
{
    bool* finished = new bool[n];
    unsigned count = 0;
    for (unsigned i = 0; i < n; i++)
        finished[i] = false;

    while (count != n)
    {
        usleep(POLL_DELAY);
        unsigned i = 0;
        for (i = 0; i < n; i++)
        {
            if (!finished[i])
                if (reqs[i]->poll())
                {
                    finished[i] = true;
                    out[i] = timestamp();
                    count++;
                }
        }
    }
    delete[] finished;
}

void out_stat(double start, double end, double* times, unsigned n, const std::vector<std::string>& names)
{
    for (unsigned i = 0; i < n; i++)
    {
        std::cout << i << " " << names[i] << " took " <<
            100. * (times[i] - start) / (end - start) << " %" << std::endl;
    }
}
#endif

#define MB (1024 * 1024)

// returns throughput in MiB/s
static inline double throughput(stxxl::int64 bytes, double seconds)
{
    if (seconds == 0.0)
        return 0.0;
    return (double)bytes / (1024 * 1024) / seconds;
}

int benchmark_files(int argc, char* argv[])
{
    uint64 offset = 0, length = 0;

    bool no_direct_io = false;
    bool sync_io = false;
    bool resize_after_open = false;
    std::string file_type = default_file_type;
    unsigned_type block_size = 0;
    unsigned int batch_size = 1;
    std::string opstr = "wv";
    unsigned pattern = 0;

    std::vector<std::string> files_arr;

    stxxl::cmdline_parser cp;

    cp.add_param_bytes("length", length,
                       "Length to write in file.");

    cp.add_param_stringlist("filename", files_arr,
                            "File path to run benchmark on.");

    cp.add_bytes('o', "offset", offset,
                 "Starting offset to write in file.");

    cp.add_flag(0, "no-direct", no_direct_io,
                "open files without O_DIRECT");

    cp.add_flag(0, "sync", sync_io,
                "open files with O_SYNC|O_DSYNC|O_RSYNC");

    cp.add_flag(0, "resize", resize_after_open,
                "resize the file size after opening, "
                "needed e.g. for creating mmap files");

    cp.add_bytes(0, "block_size", block_size,
                 "block size for operations (default 8 MiB)");

    cp.add_uint(0, "batch_size", batch_size,
                "increase (default 1) to submit several I/Os at once "
                "and report average rate");

    cp.add_string('f', "file-type", file_type,
                  "Method to open file (syscall|mmap|wincall|boostfd|...) "
                  "default: " + file_type);

    cp.add_string('p', "operations", opstr,
                  "[w]rite pattern, [r]ead without verification, "
                  "read and [v]erify pattern (default: 'wv')");

    cp.add_uint(0, "pattern", pattern,
                "32-bit pattern to write (default: block index)");

    cp.set_description(
        "Open a file using one of STXXL's file abstractions and perform "
        "write/read/verify tests on the file. "
        "Block sizes and batch size can be adjusted via command line. "
        "If length == 0 , then operation will continue till end of space "
        "(please ignore the write error). "
        "Memory consumption: block_size * batch_size * num_files");

    if (!cp.process(argc, argv))
        return -1;

    uint64 endpos = offset + length;

    if (block_size == 0)
        block_size = 8 * MB;

    if (batch_size == 0)
        batch_size = 1;

    bool do_read = false, do_write = false, do_verify = false;

    // deprecated, use --no-direct instead
    if (opstr.find("nd") != std::string::npos || opstr.find("ND") != std::string::npos) {
        no_direct_io = true;
    }

    if (opstr.find('r') != std::string::npos || opstr.find('R') != std::string::npos) {
        do_read = true;
    }
    if (opstr.find('v') != std::string::npos || opstr.find('V') != std::string::npos) {
        do_verify = true;
    }
    if (opstr.find('w') != std::string::npos || opstr.find('W') != std::string::npos) {
        do_write = true;
    }

    const char* myself = strrchr(argv[0], '/');
    if (!myself || !*(++myself))
        myself = argv[0];
    std::cout << "# " << myself << " " << stxxl::get_version_string_long();
#if STXXL_DIRECT_IO_OFF
    std::cout << " STXXL_DIRECT_IO_OFF";
#endif
    std::cout << std::endl;

    for (size_t ii = 0; ii < files_arr.size(); ii++)
    {
        std::cout << "# Add file: " << files_arr[ii] << std::endl;
    }

    const size_t nfiles = files_arr.size();
    bool verify_failed = false;

    const unsigned_type step_size = block_size * batch_size;
    const unsigned_type block_size_int = block_size / sizeof(int);
    const uint64 step_size_int = step_size / sizeof(int);

    unsigned* buffer = (unsigned*)stxxl::aligned_alloc<BLOCK_ALIGN>(step_size * nfiles);
    file** files = new file*[nfiles];
    request_ptr* reqs = new request_ptr[nfiles * batch_size];

#ifdef WATCH_TIMES
    double* r_finish_times = new double[nfiles];
    double* w_finish_times = new double[nfiles];
#endif

    double totaltimeread = 0, totaltimewrite = 0;
    stxxl::int64 totalsizeread = 0, totalsizewrite = 0;

    // fill buffer with pattern
    for (unsigned i = 0; i < nfiles * step_size_int; i++)
        buffer[i] = (pattern ? pattern : i);

    // open files
    for (unsigned i = 0; i < nfiles; i++)
    {
        int openmode = file::CREAT | file::RDWR;
        if (!no_direct_io) {
            openmode |= file::DIRECT;
        }
        if (sync_io) {
            openmode |= file::SYNC;
        }

        files[i] = stxxl::create_file(file_type, files_arr[i], openmode, i);
        if (resize_after_open)
            files[i]->set_size(endpos);
    }

    std::cout << "# Step size: "
              << step_size << " bytes per file ("
              << batch_size << " block" << (batch_size == 1 ? "" : "s") << " of "
              << block_size << " bytes)"
              << " file_type=" << file_type
              << " O_DIRECT=" << (no_direct_io ? "no" : "yes")
              << " O_SYNC=" << (sync_io ? "yes" : "no")
              << std::endl;

    stxxl::timer t_total(true);
    try {
        while (offset + uint64(step_size) <= endpos || length == 0)
        {
            const uint64 current_step_size = (length == 0) ? stxxl::int64(step_size) : std::min<stxxl::int64>(step_size, endpos - offset);
            const uint64 current_step_size_int = current_step_size / sizeof(int);
            const unsigned_type current_num_blocks = (unsigned_type)stxxl::div_ceil(current_step_size, block_size);

            std::cout << "File offset    " << std::setw(8) << offset / MB << " MiB: " << std::fixed;

            double begin = timestamp(), end = begin, elapsed;

            if (do_write)
            {
                // write block number (512 byte blocks) into each block at position 42 * sizeof(unsigned)
                for (uint64 j = 42, b = offset >> 9; j < current_step_size_int; j += 512 / sizeof(unsigned), ++b)
                {
                    for (unsigned i = 0; i < nfiles; i++)
                        buffer[current_step_size_int * i + j] = (unsigned int)b;
                }

                for (unsigned i = 0; i < nfiles; i++)
                {
                    for (unsigned_type j = 0; j < current_num_blocks; j++)
                        reqs[i * current_num_blocks + j] =
                            files[i]->awrite(buffer + current_step_size_int * i + j * block_size_int,
                                             offset + j * block_size,
                                             (unsigned_type)block_size);
                }

 #ifdef WATCH_TIMES
                watch_times(reqs, nfiles, w_finish_times);
 #else
                wait_all(reqs, nfiles * current_num_blocks);
 #endif

                end = timestamp();
                elapsed = end - begin;
                totalsizewrite += current_step_size;
                totaltimewrite += elapsed;
            }
            else {
                elapsed = 0.0;
            }

#if 0
            std::cout << "# WRITE\nFiles: " << nfiles
                      << " \nElapsed time: " << end - begin
                      << " \nThroughput: " << int(double(current_step_size * nfiles) / MB / (end - begin))
                      << " MiB/s \nPer one file:"
                      << int(double(current_step_size) / MB / (end - begin)) << " MiB/s"
                      << std::endl;
#endif

 #ifdef WATCH_TIMES
            out_stat(begin, end, w_finish_times, nfiles, files_arr);
 #endif
            std::cout << std::setw(2) << nfiles << " * "
                      << std::setw(8) << std::setprecision(3)
                      << (throughput(current_step_size, elapsed)) << " = "
                      << std::setw(8) << std::setprecision(3)
                      << (throughput(current_step_size, elapsed) * (double)nfiles) << " MiB/s write,";

            begin = end = timestamp();

            if (do_read || do_verify)
            {
                for (unsigned i = 0; i < nfiles; i++)
                {
                    for (unsigned j = 0; j < current_num_blocks; j++)
                        reqs[i * current_num_blocks + j] =
                            files[i]->aread(buffer + current_step_size_int * i + j * block_size_int,
                                            offset + j * block_size,
                                            (unsigned_type)block_size);
                }

 #ifdef WATCH_TIMES
                watch_times(reqs, nfiles, r_finish_times);
 #else
                wait_all(reqs, nfiles * current_num_blocks);
 #endif

                end = timestamp();
                elapsed = end - begin;
                totalsizeread += current_step_size;
                totaltimeread += elapsed;
            }
            else {
                elapsed = 0.0;
            }

#if 0
            std::cout << "# READ\nFiles: " << nfiles
                      << " \nElapsed time: " << end - begin
                      << " \nThroughput: " << int(double(current_step_size * nfiles) / MB / (end - begin))
                      << " MiB/s \nPer one file:"
                      << int(double(current_step_size) / MB / (end - begin)) << " MiB/s"
                      << std::endl;
#endif

            std::cout << std::setw(2) << nfiles << " * "
                      << std::setw(8) << std::setprecision(3)
                      << (throughput(current_step_size, elapsed)) << " = "
                      << std::setw(8) << std::setprecision(3)
                      << (throughput(current_step_size, elapsed) * (double)nfiles) << " MiB/s read";

#ifdef WATCH_TIMES
            out_stat(begin, end, r_finish_times, nfiles, files_arr);
#endif

            if (do_verify)
            {
                for (unsigned d = 0; d < nfiles; ++d)
                {
                    for (unsigned s = 0; s < (current_step_size >> 9); ++s) {
                        uint64 i = d * current_step_size_int + s * (512 / sizeof(unsigned)) + 42;
                        uint64 b = (offset >> 9) + s;
                        if (buffer[i] != b)
                        {
                            verify_failed = true;
                            std::cout << "Error on file " << d << " sector " << std::hex << std::setw(8) << b
                                      << " got: " << std::hex << std::setw(8) << buffer[i] << " wanted: " << std::hex << std::setw(8) << b
                                      << std::dec << std::endl;
                        }
                        buffer[i] = (pattern ? pattern : (unsigned int)i);
                    }
                }

                for (uint64 i = 0; i < nfiles * current_step_size_int; i++)
                {
                    if (buffer[i] != (pattern ? pattern : i))
                    {
                        stxxl::int64 ibuf = i / current_step_size_int;
                        uint64 pos = i % current_step_size_int;

                        std::cout << std::endl
                                  << "Error on file " << ibuf << " position " << std::hex << std::setw(8) << offset + pos * sizeof(int)
                                  << "  got: " << std::hex << std::setw(8) << buffer[i] << " wanted: " << std::hex << std::setw(8) << i
                                  << std::dec << std::endl;

                        i = (ibuf + 1) * current_step_size_int; // jump to next

                        verify_failed = true;
                    }
                }
            }
            std::cout << std::endl;

            offset += current_step_size;
        }
    }
    catch (const std::exception& ex)
    {
        std::cout << std::endl;
        STXXL_ERRMSG(ex.what());
    }
    t_total.stop();

    std::cout << "=============================================================================================" << std::endl;
    // the following line of output is parsed by misc/filebench-avgplot.sh
    std::cout << "# Average over " << std::setw(8) << stxxl::STXXL_MAX(totalsizewrite, totalsizeread) / MB << " MiB: ";
    std::cout << std::setw(2) << nfiles << " * "
              << std::setw(8) << std::setprecision(3)
              << (throughput(totalsizewrite, totaltimewrite)) << " = "
              << std::setw(8) << std::setprecision(3)
              << (throughput(totalsizewrite, totaltimewrite) * (double)nfiles) << " MiB/s write,";

    std::cout << std::setw(2) << nfiles << " * "
              << std::setw(8) << std::setprecision(3)
              << (throughput(totalsizeread, totaltimeread)) << " = "
              << std::setw(8) << std::setprecision(3)
              << (throughput(totalsizeread, totaltimeread) * (double)nfiles) << " MiB/s read"
              << std::endl;

    if (totaltimewrite != 0.0)
        std::cout << "# Write time   " << std::setw(8) << std::setprecision(3) << totaltimewrite << " s" << std::endl;
    if (totaltimeread != 0.0)
        std::cout << "# Read time    " << std::setw(8) << std::setprecision(3) << totaltimeread << " s" << std::endl;

    std::cout << "# Non-I/O time " << std::setw(8) << std::setprecision(3)
              << (t_total.seconds() - totaltimewrite - totaltimeread) << " s, average throughput "
              << std::setw(8) << std::setprecision(3)
              << (throughput(totalsizewrite + totalsizeread, t_total.seconds() - totaltimewrite - totaltimeread) * (double)nfiles) << " MiB/s"
              << std::endl;

    std::cout << "# Total time   " << std::setw(8) << std::setprecision(3) << t_total.seconds() << " s, average throughput "
              << std::setw(8) << std::setprecision(3)
              << (throughput(totalsizewrite + totalsizeread, t_total.seconds()) * (double)nfiles) << " MiB/s"
              << std::endl;

    if (do_verify)
    {
        std::cout << "# Verify: " << (verify_failed ? "FAILED." : "all okay.") << std::endl;
    }

#ifdef WATCH_TIMES
    delete[] r_finish_times;
    delete[] w_finish_times;
#endif
    delete[] reqs;
    for (unsigned i = 0; i < nfiles; i++)
        delete files[i];
    delete[] files;
    stxxl::aligned_dealloc<BLOCK_ALIGN>(buffer);

    return (verify_failed ? 1 : 0);
}

// vim: et:ts=4:sw=4
