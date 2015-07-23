//
// This file is part of khmer, https://github.com/dib-lab/khmer/, and is
// Copyright (C) Michigan State University, 2009-2015. It is licensed under
// the three-clause BSD license; see LICENSE.
// Contact: khmer-project@idyll.org
//

#ifndef PERF_METRICS_HH
#define PERF_METRICS_HH


#include <cstring>
#include <ctime>


#include "khmer.hh"

namespace khmer
{

#ifdef WITH_INTERNAL_METRICS
struct InvalidPerformanceMetricsKey : public khmer_exception {
};


struct IPerformanceMetrics {

    IPerformanceMetrics( );
    virtual ~IPerformanceMetrics( );

    inline void	    start_timers( )
    {
#if defined (__linux__)
        clock_gettime( CLOCK_REALTIME, &_temp_clock_start );
        clock_gettime( CLOCK_THREAD_CPUTIME_ID, &_temp_cpu_start );
// TODO: Create proper stopwatches for MacOS X.
#else
        memset( &_temp_clock_start, 0, sizeof( timespec ) );
        memset( &_temp_cpu_start, 0, sizeof( timespec ) );
#endif
    }
    inline void	    stop_timers( )
    {
#if defined (__linux__)
        clock_gettime( CLOCK_THREAD_CPUTIME_ID, &_temp_cpu_stop );
        clock_gettime( CLOCK_REALTIME, &_temp_clock_stop );
// TODO: Create proper stopwatches for MacOS X.
#else
        memset( &_temp_cpu_stop, 0, sizeof( timespec ) );
        memset( &_temp_clock_stop, 0, sizeof( timespec ) );
#endif
    }
    virtual void    accumulate_timer_deltas( uint32_t metrics_key )	= 0;

    // TODO: Add a printing or log file feature.

protected:

    timespec	_temp_cpu_start;
    timespec	_temp_cpu_stop;
    timespec	_temp_clock_start;
    timespec	_temp_clock_stop;

    uint64_t const  _timespec_diff_in_nsecs(
        timespec const &start, timespec const &stop
    );

};

#endif // WITH_INTERNAL_METRICS

} // namespace khmer
#endif // PERF_METRICS_HH

// vim: set ft=cpp sts=4 sw=4 tw=79:
