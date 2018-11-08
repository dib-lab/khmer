/***************************************************************************
 *  lib/io/iostats.cpp
 *
 *  Part of the STXXL. See http://stxxl.sourceforge.net
 *
 *  Copyright (C) 2002-2004 Roman Dementiev <dementiev@mpi-sb.mpg.de>
 *  Copyright (C) 2008-2010 Andreas Beckmann <beckmann@cs.uni-frankfurt.de>
 *  Copyright (C) 2009, 2010 Johannes Singler <singler@kit.edu>
 *
 *  Distributed under the Boost Software License, Version 1.0.
 *  (See accompanying file LICENSE_1_0.txt or copy at
 *  http://www.boost.org/LICENSE_1_0.txt)
 **************************************************************************/

#include <stxxl/bits/io/iostats.h>
#include <stxxl/bits/common/log.h>
#include <stxxl/bits/verbose.h>
#include <stxxl/bits/common/mutex.h>
#include <stxxl/bits/common/timer.h>
#include <stxxl/bits/common/types.h>
#include <stxxl/bits/namespace.h>

#include <string>
#include <sstream>
#include <iomanip>

STXXL_BEGIN_NAMESPACE

stats::stats()
    : reads(0),
      writes(0),
      volume_read(0),
      volume_written(0),
      c_reads(0),
      c_writes(0),
      c_volume_read(0),
      c_volume_written(0),
      t_reads(0.0),
      t_writes(0.0),
      p_reads(0.0),
      p_writes(0.0),
      p_begin_read(0.0),
      p_begin_write(0.0),
      p_ios(0.0),
      p_begin_io(0.0),
      t_waits(0.0),
      p_waits(0.0),
      p_begin_wait(0.0),
      t_wait_read(0.0),
      p_wait_read(0.0),
      p_begin_wait_read(0.0),
      t_wait_write(0.0),
      p_wait_write(0.0),
      p_begin_wait_write(0.0),
      acc_reads(0), acc_writes(0),
      acc_ios(0),
      acc_waits(0),
      acc_wait_read(0), acc_wait_write(0),
      last_reset(timestamp())
{ }

#ifndef STXXL_IO_STATS_RESET_FORBIDDEN
void stats::reset()
{
    {
        scoped_mutex_lock ReadLock(read_mutex);

        //assert(acc_reads == 0);
        if (acc_reads)
            STXXL_ERRMSG("Warning: " << acc_reads <<
                         " read(s) not yet finished");

        reads = 0;
        volume_read = 0;
        c_reads = 0;
        c_volume_read = 0;
        t_reads = 0;
        p_reads = 0.0;
    }
    {
        scoped_mutex_lock WriteLock(write_mutex);

        //assert(acc_writes == 0);
        if (acc_writes)
            STXXL_ERRMSG("Warning: " << acc_writes <<
                         " write(s) not yet finished");

        writes = 0;
        volume_written = 0;
        c_writes = 0;
        c_volume_written = 0;
        t_writes = 0.0;
        p_writes = 0.0;
    }
    {
        scoped_mutex_lock IOLock(io_mutex);

        //assert(acc_ios == 0);
        if (acc_ios)
            STXXL_ERRMSG("Warning: " << acc_ios <<
                         " io(s) not yet finished");

        p_ios = 0.0;
    }
    {
        scoped_mutex_lock WaitLock(wait_mutex);

        //assert(acc_waits == 0);
        if (acc_waits)
            STXXL_ERRMSG("Warning: " << acc_waits <<
                         " wait(s) not yet finished");

        t_waits = 0.0;
        p_waits = 0.0;
        t_wait_read = 0.0;
        p_wait_read = 0.0;
        t_wait_write = 0.0;
        p_wait_write = 0.0;
    }

    last_reset = timestamp();
}
#endif

#if STXXL_IO_STATS
void stats::write_started(unsigned_type size_, double now)
{
    if (now == 0.0)
        now = timestamp();
    {
        scoped_mutex_lock WriteLock(write_mutex);

        ++writes;
        volume_written += size_;
        double diff = now - p_begin_write;
        t_writes += double(acc_writes) * diff;
        p_begin_write = now;
        p_writes += (acc_writes++) ? diff : 0.0;
    }
    {
        scoped_mutex_lock IOLock(io_mutex);

        double diff = now - p_begin_io;
        p_ios += (acc_ios++) ? diff : 0.0;
        p_begin_io = now;
    }
}

void stats::write_canceled(unsigned_type size_)
{
    {
        scoped_mutex_lock WriteLock(write_mutex);

        --writes;
        volume_written -= size_;
    }
    write_finished();
}

void stats::write_finished()
{
    double now = timestamp();
    {
        scoped_mutex_lock WriteLock(write_mutex);

        double diff = now - p_begin_write;
        t_writes += double(acc_writes) * diff;
        p_begin_write = now;
        p_writes += (acc_writes--) ? diff : 0.0;
    }
    {
        scoped_mutex_lock IOLock(io_mutex);

        double diff = now - p_begin_io;
        p_ios += (acc_ios--) ? diff : 0.0;
        p_begin_io = now;
    }
}

void stats::write_cached(unsigned_type size_)
{
    scoped_mutex_lock WriteLock(write_mutex);

    ++c_writes;
    c_volume_written += size_;
}

void stats::read_started(unsigned_type size_, double now)
{
    if (now == 0.0)
        now = timestamp();
    {
        scoped_mutex_lock ReadLock(read_mutex);

        ++reads;
        volume_read += size_;
        double diff = now - p_begin_read;
        t_reads += double(acc_reads) * diff;
        p_begin_read = now;
        p_reads += (acc_reads++) ? diff : 0.0;
    }
    {
        scoped_mutex_lock IOLock(io_mutex);

        double diff = now - p_begin_io;
        p_ios += (acc_ios++) ? diff : 0.0;
        p_begin_io = now;
    }
}

void stats::read_canceled(unsigned_type size_)
{
    {
        scoped_mutex_lock ReadLock(read_mutex);

        --reads;
        volume_read -= size_;
    }
    read_finished();
}

void stats::read_finished()
{
    double now = timestamp();
    {
        scoped_mutex_lock ReadLock(read_mutex);

        double diff = now - p_begin_read;
        t_reads += double(acc_reads) * diff;
        p_begin_read = now;
        p_reads += (acc_reads--) ? diff : 0.0;
    }
    {
        scoped_mutex_lock IOLock(io_mutex);

        double diff = now - p_begin_io;
        p_ios += (acc_ios--) ? diff : 0.0;
        p_begin_io = now;
    }
}

void stats::read_cached(unsigned_type size_)
{
    scoped_mutex_lock ReadLock(read_mutex);

    ++c_reads;
    c_volume_read += size_;
}
#endif

#ifndef STXXL_DO_NOT_COUNT_WAIT_TIME
void stats::wait_started(wait_op_type wait_op)
{
    double now = timestamp();
    {
        scoped_mutex_lock WaitLock(wait_mutex);

        double diff = now - p_begin_wait;
        t_waits += double(acc_waits) * diff;
        p_begin_wait = now;
        p_waits += (acc_waits++) ? diff : 0.0;

        if (wait_op == WAIT_OP_READ) {
            diff = now - p_begin_wait_read;
            t_wait_read += double(acc_wait_read) * diff;
            p_begin_wait_read = now;
            p_wait_read += (acc_wait_read++) ? diff : 0.0;
        }
        else /* if (wait_op == WAIT_OP_WRITE) */ {
            // wait_any() is only used from write_pool and buffered_writer, so account WAIT_OP_ANY for WAIT_OP_WRITE, too
            diff = now - p_begin_wait_write;
            t_wait_write += double(acc_wait_write) * diff;
            p_begin_wait_write = now;
            p_wait_write += (acc_wait_write++) ? diff : 0.0;
        }
    }
}

void stats::wait_finished(wait_op_type wait_op)
{
    double now = timestamp();
    {
        scoped_mutex_lock WaitLock(wait_mutex);

        double diff = now - p_begin_wait;
        t_waits += double(acc_waits) * diff;
        p_begin_wait = now;
        p_waits += (acc_waits--) ? diff : 0.0;

        if (wait_op == WAIT_OP_READ) {
            double diff2 = now - p_begin_wait_read;
            t_wait_read += double(acc_wait_read) * diff2;
            p_begin_wait_read = now;
            p_wait_read += (acc_wait_read--) ? diff2 : 0.0;
        }
        else /* if (wait_op == WAIT_OP_WRITE) */ {
            double diff2 = now - p_begin_wait_write;
            t_wait_write += double(acc_wait_write) * diff2;
            p_begin_wait_write = now;
            p_wait_write += (acc_wait_write--) ? diff2 : 0.0;
        }
#ifdef STXXL_WAIT_LOG_ENABLED
        std::ofstream* waitlog = stxxl::logger::get_instance()->waitlog_stream();
        if (waitlog)
            *waitlog << (now - last_reset) << "\t"
                     << ((wait_op == WAIT_OP_READ) ? diff : 0.0) << "\t"
                     << ((wait_op != WAIT_OP_READ) ? diff : 0.0) << "\t"
                     << t_wait_read << "\t" << t_wait_write << std::endl << std::flush;
#endif
    }
}
#endif

void stats::_reset_io_wait_time()
{
#ifndef STXXL_DO_NOT_COUNT_WAIT_TIME
    {
        scoped_mutex_lock WaitLock(wait_mutex);

        //assert(acc_waits == 0);
        if (acc_waits)
            STXXL_ERRMSG("Warning: " << acc_waits <<
                         " wait(s) not yet finished");

        t_waits = 0.0;
        p_waits = 0.0;
    }
#endif
}

std::string format_with_SI_IEC_unit_multiplier(uint64 number, const char* unit, int multiplier)
{
    // may not overflow, std::numeric_limits<uint64>::max() == 16 EB
    static const char* endings[] = { "", "k", "M", "G", "T", "P", "E" };
    static const char* binary_endings[] = { "", "Ki", "Mi", "Gi", "Ti", "Pi", "Ei" };
    std::ostringstream out;
    out << number << ' ';
    int scale = 0;
    double number_d = (double)number;
    double multiplier_d = multiplier;
    while (number_d >= multiplier_d)
    {
        number_d /= multiplier_d;
        ++scale;
    }
    if (scale > 0)
        out << '(' << std::fixed << std::setprecision(3) << number_d << ' '
            << (multiplier == 1024 ? binary_endings[scale] : endings[scale])
            << (unit ? unit : "") << ") ";
    else if (unit && *unit)
        out << unit << ' ';
    return out.str();
}

std::ostream& operator << (std::ostream& o, const stats_data& s)
{
#define hr add_IEC_binary_multiplier
    o << "STXXL I/O statistics" << std::endl;
#if STXXL_IO_STATS
    o << " total number of reads                      : " << hr(s.get_reads()) << std::endl;
    o << " average block size (read)                  : "
      << hr(s.get_reads() ? s.get_read_volume() / s.get_reads() : 0, "B") << std::endl;
    o << " number of bytes read from disks            : " << hr(s.get_read_volume(), "B") << std::endl;
    o << " time spent in serving all read requests    : " << s.get_read_time() << " s"
      << " @ " << ((double)s.get_read_volume() / 1048576.0 / s.get_read_time()) << " MiB/s"
      << std::endl;
    o << " time spent in reading (parallel read time) : " << s.get_pread_time() << " s"
      << " @ " << ((double)s.get_read_volume() / 1048576.0 / s.get_pread_time()) << " MiB/s"
      << std::endl;
    if (s.get_cached_reads()) {
        o << " total number of cached reads               : " << hr(s.get_cached_reads()) << std::endl;
        o << " average block size (cached read)           : " << hr(s.get_cached_read_volume() / s.get_cached_reads(), "B") << std::endl;
        o << " number of bytes read from cache            : " << hr(s.get_cached_read_volume(), "B") << std::endl;
    }
    if (s.get_cached_writes()) {
        o << " total number of cached writes              : " << hr(s.get_cached_writes()) << std::endl;
        o << " average block size (cached write)          : " << hr(s.get_cached_written_volume() / s.get_cached_writes(), "B") << std::endl;
        o << " number of bytes written to cache           : " << hr(s.get_cached_written_volume(), "B") << std::endl;
    }
    o << " total number of writes                     : " << hr(s.get_writes()) << std::endl;
    o << " average block size (write)                 : "
      << hr(s.get_writes() ? s.get_written_volume() / s.get_writes() : 0, "B") << std::endl;
    o << " number of bytes written to disks           : " << hr(s.get_written_volume(), "B") << std::endl;
    o << " time spent in serving all write requests   : " << s.get_write_time() << " s"
      << " @ " << ((double)s.get_written_volume() / 1048576.0 / s.get_write_time()) << " MiB/s"
      << std::endl;
    o << " time spent in writing (parallel write time): " << s.get_pwrite_time() << " s"
      << " @ " << ((double)s.get_written_volume() / 1048576.0 / s.get_pwrite_time()) << " MiB/s"
      << std::endl;
    o << " time spent in I/O (parallel I/O time)      : " << s.get_pio_time() << " s"
      << " @ " << ((double)(s.get_read_volume() + s.get_written_volume()) / 1048576.0 / s.get_pio_time()) << " MiB/s"
      << std::endl;
#else
    o << " n/a" << std::endl;
#endif
#ifndef STXXL_DO_NOT_COUNT_WAIT_TIME
    o << " I/O wait time                              : " << s.get_io_wait_time() << " s" << std::endl;
    if (s.get_wait_read_time() != 0.0)
        o << " I/O wait4read time                         : " << s.get_wait_read_time() << " s" << std::endl;
    if (s.get_wait_write_time() != 0.0)
        o << " I/O wait4write time                        : " << s.get_wait_write_time() << " s" << std::endl;
#endif
    o << " Time since the last reset                  : " << s.get_elapsed_time() << " s" << std::endl;
    return o;
#undef hr
}

STXXL_END_NAMESPACE
// vim: et:ts=4:sw=4
