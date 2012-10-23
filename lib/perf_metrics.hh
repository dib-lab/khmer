#ifndef PERF_METRICS_HH
#define PERF_METRICS_HH


#if (__cplusplus >= 201103L)
#   include <cstdint>
#else
extern "C"
{
#   include <stdint.h>
}
#endif

#include <exception>


namespace khmer
{


struct InvalidPerformanceMetricsKey : public std:: exception
{ };


struct IPerformanceMetrics
{

	    IPerformanceMetrics( );
    virtual ~IPerformanceMetrics( );
    
    inline void	    start_timers( )
    {
	clock_gettime( CLOCK_REALTIME, &_temp_clock_start );
	clock_gettime( CLOCK_THREAD_CPUTIME_ID, &_temp_cpu_start );
    }
    inline void	    stop_timers( )
    {
	clock_gettime( CLOCK_THREAD_CPUTIME_ID, &_temp_cpu_stop );
	clock_gettime( CLOCK_REALTIME, &_temp_clock_stop );
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

} // namespace khmer

#endif // PERF_METRICS_HH

// vim: set ft=cpp sts=4 sw=4 tw=79:
