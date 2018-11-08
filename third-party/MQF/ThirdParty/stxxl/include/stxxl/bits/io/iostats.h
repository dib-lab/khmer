/***************************************************************************
 *  include/stxxl/bits/io/iostats.h
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

#ifndef STXXL_IO_IOSTATS_HEADER
#define STXXL_IO_IOSTATS_HEADER

#ifndef STXXL_IO_STATS
 #define STXXL_IO_STATS 1
#endif

#include <stxxl/bits/namespace.h>
#include <stxxl/bits/deprecated.h>
#include <stxxl/bits/common/mutex.h>
#include <stxxl/bits/common/timer.h>
#include <stxxl/bits/common/types.h>
#include <stxxl/bits/common/utils.h>
#include <stxxl/bits/unused.h>
#include <stxxl/bits/singleton.h>

#include <iostream>
#include <string>

STXXL_BEGIN_NAMESPACE

//! \addtogroup iolayer
//!
//! \{

//! Collects various I/O statistics.
//! \remarks is a singleton
class stats : public singleton<stats>
{
    friend class singleton<stats>;

    unsigned reads, writes;                     // number of operations
    int64 volume_read, volume_written;          // number of bytes read/written
    unsigned c_reads, c_writes;                 // number of cached operations
    int64 c_volume_read, c_volume_written;      // number of bytes read/written from/to cache
    double t_reads, t_writes;                   // seconds spent in operations
    double p_reads, p_writes;                   // seconds spent in parallel operations
    double p_begin_read, p_begin_write;         // start time of parallel operation
    double p_ios;                               // seconds spent in all parallel I/O operations (read and write)
    double p_begin_io;
    double t_waits, p_waits;                    // seconds spent waiting for completion of I/O operations
    double p_begin_wait;
    double t_wait_read, p_wait_read;
    double p_begin_wait_read;
    double t_wait_write, p_wait_write;
    double p_begin_wait_write;
    int acc_reads, acc_writes;                  // number of requests, participating in parallel operation
    int acc_ios;
    int acc_waits;
    int acc_wait_read, acc_wait_write;
    double last_reset;
    mutex read_mutex, write_mutex, io_mutex, wait_mutex;

    stats();

public:
    enum wait_op_type {
        WAIT_OP_ANY,
        WAIT_OP_READ,
        WAIT_OP_WRITE
    };

    class scoped_read_write_timer
    {
        typedef unsigned_type size_type;

        bool is_write;
#if STXXL_IO_STATS
        bool running;
#endif

    public:
        scoped_read_write_timer(size_type size, bool is_write = false)
            : is_write(is_write)
#if STXXL_IO_STATS
              , running(false)
#endif
        {
            start(size);
        }

        ~scoped_read_write_timer()
        {
            stop();
        }

        void start(size_type size)
        {
#if STXXL_IO_STATS
            if (!running) {
                running = true;
                if (is_write)
                    stats::get_instance()->write_started(size);
                else
                    stats::get_instance()->read_started(size);
            }
#else
            STXXL_UNUSED(size);
#endif
        }

        void stop()
        {
#if STXXL_IO_STATS
            if (running) {
                if (is_write)
                    stats::get_instance()->write_finished();
                else
                    stats::get_instance()->read_finished();
                running = false;
            }
#endif
        }
    };

    class scoped_write_timer
    {
        typedef unsigned_type size_type;

#if STXXL_IO_STATS
        bool running;
#endif

    public:
        scoped_write_timer(size_type size)
#if STXXL_IO_STATS
            : running(false)
#endif
        {
            start(size);
        }

        ~scoped_write_timer()
        {
            stop();
        }

        void start(size_type size)
        {
#if STXXL_IO_STATS
            if (!running) {
                running = true;
                stats::get_instance()->write_started(size);
            }
#else
            STXXL_UNUSED(size);
#endif
        }

        void stop()
        {
#if STXXL_IO_STATS
            if (running) {
                stats::get_instance()->write_finished();
                running = false;
            }
#endif
        }
    };

    class scoped_read_timer
    {
        typedef unsigned_type size_type;

#if STXXL_IO_STATS
        bool running;
#endif

    public:
        scoped_read_timer(size_type size)
#if STXXL_IO_STATS
            : running(false)
#endif
        {
            start(size);
        }

        ~scoped_read_timer()
        {
            stop();
        }

        void start(size_type size)
        {
#if STXXL_IO_STATS
            if (!running) {
                running = true;
                stats::get_instance()->read_started(size);
            }
#else
            STXXL_UNUSED(size);
#endif
        }

        void stop()
        {
#if STXXL_IO_STATS
            if (running) {
                stats::get_instance()->read_finished();
                running = false;
            }
#endif
        }
    };

    class scoped_wait_timer
    {
#ifndef STXXL_DO_NOT_COUNT_WAIT_TIME
        bool running;
        wait_op_type wait_op;
#endif

    public:
        scoped_wait_timer(wait_op_type wait_op, bool measure_time = true)
#ifndef STXXL_DO_NOT_COUNT_WAIT_TIME
            : running(false), wait_op(wait_op)
#endif
        {
            if (measure_time)
                start();
        }

        ~scoped_wait_timer()
        {
            stop();
        }

        void start()
        {
#ifndef STXXL_DO_NOT_COUNT_WAIT_TIME
            if (!running) {
                running = true;
                stats::get_instance()->wait_started(wait_op);
            }
#endif
        }

        void stop()
        {
#ifndef STXXL_DO_NOT_COUNT_WAIT_TIME
            if (running) {
                stats::get_instance()->wait_finished(wait_op);
                running = false;
            }
#endif
        }
    };

public:
    //! Returns total number of reads.
    //! \return total number of reads
    unsigned get_reads() const
    {
        return reads;
    }

    //! Returns total number of writes.
    //! \return total number of writes
    unsigned get_writes() const
    {
        return writes;
    }

    //! Returns number of bytes read from disks.
    //! \return number of bytes read
    int64 get_read_volume() const
    {
        return volume_read;
    }

    //! Returns number of bytes written to the disks.
    //! \return number of bytes written
    int64 get_written_volume() const
    {
        return volume_written;
    }

    //! Returns total number of reads served from cache.
    //! \return total number of cached reads
    unsigned get_cached_reads() const
    {
        return c_reads;
    }

    //! Returns total number of cached writes.
    //! \return total number of cached writes
    unsigned get_cached_writes() const
    {
        return c_writes;
    }

    //! Returns number of bytes read from cache.
    //! \return number of bytes read from cache
    int64 get_cached_read_volume() const
    {
        return c_volume_read;
    }

    //! Returns number of bytes written to the cache.
    //! \return number of bytes written to cache
    int64 get_cached_written_volume() const
    {
        return c_volume_written;
    }

    //! Time that would be spent in read syscalls if all parallel reads were serialized.
    //! \return seconds spent in reading
    double get_read_time() const
    {
        return t_reads;
    }

    //! Time that would be spent in write syscalls if all parallel writes were serialized.
    //! \return seconds spent in writing
    double get_write_time() const
    {
        return t_writes;
    }

    //! Period of time when at least one I/O thread was executing a read.
    //! \return seconds spent in reading
    double get_pread_time() const
    {
        return p_reads;
    }

    //! Period of time when at least one I/O thread was executing a write.
    //! \return seconds spent in writing
    double get_pwrite_time() const
    {
        return p_writes;
    }

    //! Period of time when at least one I/O thread was executing a read or a write.
    //! \return seconds spent in I/O
    double get_pio_time() const
    {
        return p_ios;
    }

    //! I/O wait time counter.
    //! \return number of seconds spent in I/O waiting functions \link
    //! request::wait request::wait \endlink, \c wait_any and \c wait_all
    double get_io_wait_time() const
    {
        return t_waits;
    }

    double get_wait_read_time() const
    {
        return t_wait_read;
    }

    double get_wait_write_time() const
    {
        return t_wait_write;
    }

    //! Return time of the last reset.
    //! \return seconds passed from the last reset()
    double get_last_reset_time() const
    {
        return last_reset;
    }

#ifndef STXXL_IO_STATS_RESET_FORBIDDEN
    //! Resets I/O time counters (including I/O wait counter).
    STXXL_DEPRECATED(void reset());
#endif

    //! Resets I/O wait time counter.
    STXXL_DEPRECATED(void _reset_io_wait_time());

    // for library use
    void write_started(unsigned_type size_, double now = 0.0);
    void write_canceled(unsigned_type size_);
    void write_finished();
    void write_cached(unsigned_type size_);
    void read_started(unsigned_type size_, double now = 0.0);
    void read_canceled(unsigned_type size_);
    void read_finished();
    void read_cached(unsigned_type size_);
    void wait_started(wait_op_type wait_op);
    void wait_finished(wait_op_type wait_op);
};

#if !STXXL_IO_STATS
inline void stats::write_started(unsigned_type size_, double now)
{
    STXXL_UNUSED(size_);
    STXXL_UNUSED(now);
}
inline void stats::write_cached(unsigned_type size_)
{
    STXXL_UNUSED(size_);
}
inline void stats::write_finished() { }
inline void stats::read_started(unsigned_type size_, double now)
{
    STXXL_UNUSED(size_);
    STXXL_UNUSED(now);
}
inline void stats::read_cached(unsigned_type size_)
{
    STXXL_UNUSED(size_);
}
inline void stats::read_finished() { }
#endif
#ifdef STXXL_DO_NOT_COUNT_WAIT_TIME
inline void stats::wait_started(wait_op_type) { }
inline void stats::wait_finished(wait_op_type) { }
#endif

class stats_data
{
    //! number of operations
    unsigned reads, writes;
    //! number of bytes read/written
    int64 volume_read, volume_written;
    //! number of cached operations
    unsigned c_reads, c_writes;
    //! number of bytes read/written from/to cache
    int64 c_volume_read, c_volume_written;
    //! seconds spent in operations
    double t_reads, t_writes;
    //! seconds spent in parallel operations
    double p_reads, p_writes;
    //! seconds spent in all parallel I/O operations (read and write)
    double p_ios;
    //! seconds spent waiting for completion of I/O operations
    double t_wait;
    double t_wait_read, t_wait_write;
    double elapsed;

public:
    stats_data()
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
          p_ios(0.0),
          t_wait(0.0),
          t_wait_read(0.0),
          t_wait_write(0.0),
          elapsed(0.0)
    { }

    stats_data(const stats& s)
        : reads(s.get_reads()),
          writes(s.get_writes()),
          volume_read(s.get_read_volume()),
          volume_written(s.get_written_volume()),
          c_reads(s.get_cached_reads()),
          c_writes(s.get_cached_writes()),
          c_volume_read(s.get_cached_read_volume()),
          c_volume_written(s.get_cached_written_volume()),
          t_reads(s.get_read_time()),
          t_writes(s.get_write_time()),
          p_reads(s.get_pread_time()),
          p_writes(s.get_pwrite_time()),
          p_ios(s.get_pio_time()),
          t_wait(s.get_io_wait_time()),
          t_wait_read(s.get_wait_read_time()),
          t_wait_write(s.get_wait_write_time()),
          elapsed(timestamp() - s.get_last_reset_time())
    { }

    stats_data operator + (const stats_data& a) const
    {
        stats_data s;
        s.reads = reads + a.reads;
        s.writes = writes + a.writes;
        s.volume_read = volume_read + a.volume_read;
        s.volume_written = volume_written + a.volume_written;
        s.c_reads = c_reads + a.c_reads;
        s.c_writes = c_writes + a.c_writes;
        s.c_volume_read = c_volume_read + a.c_volume_read;
        s.c_volume_written = c_volume_written + a.c_volume_written;
        s.t_reads = t_reads + a.t_reads;
        s.t_writes = t_writes + a.t_writes;
        s.p_reads = p_reads + a.p_reads;
        s.p_writes = p_writes + a.p_writes;
        s.p_ios = p_ios + a.p_ios;
        s.t_wait = t_wait + a.t_wait;
        s.t_wait_read = t_wait_read + a.t_wait_read;
        s.t_wait_write = t_wait_write + a.t_wait_write;
        s.elapsed = elapsed + a.elapsed;
        return s;
    }

    stats_data operator - (const stats_data& a) const
    {
        stats_data s;
        s.reads = reads - a.reads;
        s.writes = writes - a.writes;
        s.volume_read = volume_read - a.volume_read;
        s.volume_written = volume_written - a.volume_written;
        s.c_reads = c_reads - a.c_reads;
        s.c_writes = c_writes - a.c_writes;
        s.c_volume_read = c_volume_read - a.c_volume_read;
        s.c_volume_written = c_volume_written - a.c_volume_written;
        s.t_reads = t_reads - a.t_reads;
        s.t_writes = t_writes - a.t_writes;
        s.p_reads = p_reads - a.p_reads;
        s.p_writes = p_writes - a.p_writes;
        s.p_ios = p_ios - a.p_ios;
        s.t_wait = t_wait - a.t_wait;
        s.t_wait_read = t_wait_read - a.t_wait_read;
        s.t_wait_write = t_wait_write - a.t_wait_write;
        s.elapsed = elapsed - a.elapsed;
        return s;
    }

    unsigned get_reads() const
    {
        return reads;
    }

    unsigned get_writes() const
    {
        return writes;
    }

    int64 get_read_volume() const
    {
        return volume_read;
    }

    int64 get_written_volume() const
    {
        return volume_written;
    }

    unsigned get_cached_reads() const
    {
        return c_reads;
    }

    unsigned get_cached_writes() const
    {
        return c_writes;
    }

    int64 get_cached_read_volume() const
    {
        return c_volume_read;
    }

    int64 get_cached_written_volume() const
    {
        return c_volume_written;
    }

    double get_read_time() const
    {
        return t_reads;
    }

    double get_write_time() const
    {
        return t_writes;
    }

    double get_pread_time() const
    {
        return p_reads;
    }

    double get_pwrite_time() const
    {
        return p_writes;
    }

    double get_pio_time() const
    {
        return p_ios;
    }

    double get_elapsed_time() const
    {
        return elapsed;
    }

    double get_io_wait_time() const
    {
        return t_wait;
    }

    double get_wait_read_time() const
    {
        return t_wait_read;
    }

    double get_wait_write_time() const
    {
        return t_wait_write;
    }
};

std::ostream& operator << (std::ostream& o, const stats_data& s);

inline std::ostream& operator << (std::ostream& o, const stats& s)
{
    o << stxxl::stats_data(s);
    return o;
}

std::string format_with_SI_IEC_unit_multiplier(uint64 number, const char* unit = "", int multiplier = 1000);

inline std::string add_IEC_binary_multiplier(uint64 number, const char* unit = "")
{
    return format_with_SI_IEC_unit_multiplier(number, unit, 1024);
}

inline std::string add_SI_multiplier(uint64 number, const char* unit = "")
{
    return format_with_SI_IEC_unit_multiplier(number, unit, 1000);
}

//! \}

STXXL_END_NAMESPACE

#endif // !STXXL_IO_IOSTATS_HEADER
// vim: et:ts=4:sw=4
