/***************************************************************************
 *  include/stxxl/bits/common/timer.h
 *
 *  Part of the STXXL. See http://stxxl.sourceforge.net
 *
 *  Copyright (C) 2002, 2005 Roman Dementiev <dementiev@mpi-sb.mpg.de>
 *  Copyright (C) 2007-2009 Andreas Beckmann <beckmann@cs.uni-frankfurt.de>
 *  Copyright (C) 2008 Johannes Singler <singler@ira.uka.de>
 *  Copyright (C) 2013-2014 Timo Bingmann <tb@panthema.net>
 *
 *  Distributed under the Boost Software License, Version 1.0.
 *  (See accompanying file LICENSE_1_0.txt or copy at
 *  http://www.boost.org/LICENSE_1_0.txt)
 **************************************************************************/

#ifndef STXXL_COMMON_TIMER_HEADER
#define STXXL_COMMON_TIMER_HEADER

#include <stxxl/bits/config.h>
#include <stxxl/bits/namespace.h>
#include <stxxl/bits/verbose.h>
#include <stxxl/bits/common/utils.h>

#if STXXL_BOOST_TIMESTAMP
  #include <boost/date_time/posix_time/posix_time.hpp>
  #include <cmath>
#elif STXXL_WINDOWS
  #ifndef NOMINMAX
    #define NOMINMAX
  #endif
  #include <windows.h>
#else
  #include <ctime>
  #include <sys/time.h>
#endif

STXXL_BEGIN_NAMESPACE

//! \addtogroup support
//! \{

//! Returns number of seconds since the epoch, high resolution.
inline double
timestamp()
{
#if STXXL_BOOST_TIMESTAMP
    boost::posix_time::ptime MyTime = boost::posix_time::microsec_clock::local_time();
    boost::posix_time::time_duration Duration =
        MyTime - boost::posix_time::time_from_string("1970-01-01 00:00:00.000");
    double sec = double(Duration.hours()) * 3600. +
                 double(Duration.minutes()) * 60. +
                 double(Duration.seconds()) +
                 double(Duration.fractional_seconds()) / (pow(10., Duration.num_fractional_digits()));
    return sec;
#elif STXXL_WINDOWS
    return GetTickCount() / 1000.0;
#else
    struct timeval tp;
    gettimeofday(&tp, NULL);
    return double(tp.tv_sec) + double(tp.tv_usec) / 1000000.;
#endif
}

/*!
 * Class timer is a simple stop watch timer. It uses the timestamp() function
 * to get the current time when start() is called. Then, after some processing,
 * the function stop() functions can be called, or seconds() and other
 * accessors can be called directly.
 */
class timer
{
    //! boolean whether the stopwatch timer is currently running
    bool running;

    //! total accumulated time in seconds.
    double accumulated;

    //! last start time of the stopwatch
    double last_clock;

    //! return current timestamp
    static inline double timestamp()
    {
        return stxxl::timestamp();
    }

public:
    //! boolean indicating that this class does real timing
    static const bool is_real = true;

    //! initialize and optionally immediately start the timer
    inline timer(bool start_immediately = false)
        : running(false), accumulated(0), last_clock(0)
    {
        if (start_immediately) start();
    }

    //! start timer
    inline void start()
    {
        running = true;
        last_clock = timestamp();
    }

    //! stop timer
    inline void stop()
    {
        running = false;
        accumulated += timestamp() - last_clock;
    }

    //! return accumulated time
    inline void reset()
    {
        accumulated = 0.;
        last_clock = timestamp();
    }

    //! return currently accumulated time in milliseconds
    inline double mseconds() const
    {
        if (running)
            return (accumulated + timestamp() - last_clock) * 1000.;

        return (accumulated * 1000.);
    }

    //! return currently accumulated time in microseconds
    inline double useconds() const
    {
        if (running)
            return (accumulated + timestamp() - last_clock) * 1000000.;

        return (accumulated * 1000000.);
    }

    //! return currently accumulated time in seconds (as double)
    inline double seconds() const
    {
        if (running)
            return (accumulated + timestamp() - last_clock);

        return (accumulated);
    }

    //! accumulate elapsed time from another timer
    inline timer& operator += (const timer& tm)
    {
#if STXXL_PARALLEL
#pragma omp atomic
#endif
        accumulated += tm.seconds();
        return *this;
    }

    //! direct <<-operator for ostream. Can be used for printing with std::cout.
    friend std::ostream& operator << (std::ostream& os, const timer& t)
    {
        return os << t.seconds() << 's';
    }
};

/*!
 * Class fake_timer is a drop-in replacement for timer, which does
 * nothing. Using the fake class, timers can quickly be disabled in release
 * builds, but still be available for debugging session.
 *
 * \see timer
 */
class fake_timer
{
public:
    //! boolean indicating that this class does NOT do real timing
    static const bool is_real = false;

    //! initialize and optionally immediately start the timer
    fake_timer(bool = false)
    { }

    //! start timer
    void start()
    { }

    //! stop timer
    void stop()
    { }

    //! return accumulated time
    void reset()
    { }

    //! return currently accumulated time in milliseconds
    double mseconds() const
    {
        return std::numeric_limits<double>::quiet_NaN();
    }

    //! return currently accumulated time in microseconds
    double useconds() const
    {
        return std::numeric_limits<double>::quiet_NaN();
    }

    //! return currently accumulated time in seconds (as double)
    double seconds() const
    {
        return std::numeric_limits<double>::quiet_NaN();
    }

    //! accumulate elapsed time from another timer
    inline fake_timer& operator += (const fake_timer&)
    {
        return *this;
    }

    //! direct <<-operator for ostream. Can be used for printing with std::cout.
    friend std::ostream& operator << (std::ostream& os, const fake_timer& t)
    {
        return os << t.seconds() << 's';
    }
};

/*!
 * Simple scoped timer, which takes a text message and prints the duration
 * until the scope is destroyed.
 */
class scoped_print_timer
{
protected:
    //! message
    std::string m_message;

    //! bytes processed
    uint64 m_bytes;

    //! timer
    stxxl::timer m_timer;

public:
    //! save message and start timer
    scoped_print_timer(const std::string& message, const uint64 bytes = 0)
        : m_message(message),
          m_bytes(bytes),
          m_timer(true)
    {
        STXXL_MSG("Starting " << message);
    }

    //! on destruction: tell the time
    ~scoped_print_timer()
    {
        if (m_bytes == 0) {
            STXXL_MSG("Finished "
                      << m_message
                      << " after " << m_timer.seconds() << " seconds");
        }
        else {
            double bps = (double)m_bytes / m_timer.seconds();

            STXXL_MSG("Finished "
                      << m_message
                      << " after " << m_timer.seconds() << " seconds. "
                      << "Processed " << format_IEC_size(m_bytes) << "B"
                      << " @ " << format_IEC_size((uint64)bps) << "B/s");
        }
    }

    //! constant access to enclosed timer
    const stxxl::timer & timer() const
    {
        return m_timer;
    }
};

//! \}

STXXL_END_NAMESPACE

#endif // !STXXL_COMMON_TIMER_HEADER
// vim: et:ts=4:sw=4
