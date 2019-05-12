/***************************************************************************
 *  include/stxxl/bits/parallel/timing.h
 *
 *  Provides a simple tool to do performance debugging, also in parallel code.
 *  Extracted from MCSTL - http://algo2.iti.uni-karlsruhe.de/singler/mcstl/
 *
 *  Part of the STXXL. See http://stxxl.sourceforge.net
 *
 *  Copyright (C) 2006 Johannes Singler <singler@ira.uka.de>
 *  Copyright (C) 2014 Timo Bingmann <tb@panthema.net>
 *
 *  Distributed under the Boost Software License, Version 1.0.
 *  (See accompanying file LICENSE_1_0.txt or copy at
 *  http://www.boost.org/LICENSE_1_0.txt)
 **************************************************************************/

#ifndef STXXL_PARALLEL_TIMING_HEADER
#define STXXL_PARALLEL_TIMING_HEADER

#include <cstdio>
#include <cstring>
#include <cassert>

#include <stxxl/bits/config.h>
#include <stxxl/bits/parallel/tags.h>

#if STXXL_PARALLEL
  #include <omp.h>
#endif

STXXL_BEGIN_NAMESPACE

namespace parallel {

/** Type of of point in time, used for the Timing classes. */
typedef double point_in_time;

template <typename tag, typename must_be_int = int>
class Timing;

#if STXXL_PARALLEL

/** A class that provides simple run time measurements, also for parallel code.
 *  \param tag If active_tag, then the measurements are actually done.
 *  Otherwise, no code at all is emitted by the compiler. */
template <typename must_be_int>
class Timing<active_tag, must_be_int>
{
private:
    static const int max_points_in_time = 100;
    point_in_time points_in_time[max_points_in_time];
    point_in_time active, last_start;
    int pos;
    char* str;
    const char* tags[max_points_in_time];

public:
    Timing()
    {
        str = NULL;
        pos = 0;
        active = 0.0;
        last_start = -1.0;
    }

    ~Timing()
    {
        delete[] str;
    }

    /** Take a running time measurement.
     *  \param  tag  Optional description that will be output again with the timings.
     *  It should describe the operation before the tic(). To time a series of \c n
     *  operations, there should be \c n+1 calls to tic(), and one call to print(). */
    inline void tic(const char* tag = NULL)
    {
        points_in_time[pos] = omp_get_wtime();
        tags[pos] = tag;
        pos++;
    }

    /** Start the running time measurement.
     *
     *  Should be paired with stop(). */
    inline void start()
    {
        assert(last_start == -1.0);
        last_start = omp_get_wtime();
    }

    /** Stop the running time measurement.
     *
     *  Should be paired with start(). */
    inline void stop()
    {
        assert(last_start != -1.0);
        active += (omp_get_wtime() - last_start);
        last_start = -1.0;
    }

    /** Reset running time accumulation. */
    inline void reset()
    {
        active = 0.0;
        last_start = -1.0;
    }

    /** Accumulate the time between all pairs of start() and stop() so far */
    inline point_in_time active_time()
    {
        return active;
    }

    /** Total time between first and last tic() */
    inline point_in_time total_time()
    {
        return (points_in_time[pos - 1] - points_in_time[0]) * 1000.0;
    }

private:
    /** Construct string to print out, presenting the timings. */
    const char * c_str()
    {
        //avoid stream library here, to avoid cyclic dependencies in header files

        char tmp[1000];

        if (!str)
            str = new char[pos * 200];
        else
            str[0] = '\0';

        sprintf(str, "t %2d      T[ms]", omp_get_thread_num());
        strcat(str, "\n");

        for (int i = 0; i < pos; )
        {
            point_in_time last = points_in_time[i];
            i++;
            if (i == pos)
                break;
            if (tags[i] == NULL)
                sprintf(tmp, "%2d:     ", i - 1);
            else
                sprintf(tmp, "%20s:     ", tags[i]);
            strcat(str, tmp);

            sprintf(tmp, "%7.2f     ", (points_in_time[i] - last) * 1000.0);
            strcat(str, tmp);
            strcat(str, "\n");
        }

        return str;
    }

public:
    /** Print the running times between the tic()s. */
    void print()
    {
        printf("print\n");

#pragma omp barrier

#pragma omp master
        printf("\n\n");

#pragma omp critical
        printf("%s\n", c_str());
    }
};

#endif // STXXL_PARALLEL

/** A class that provides simple run time measurements, also for parallel code.
 *  \param tag If active_tag, then the measurements are actually done,
 *  otherwise, no code at all is emitted by the compiler. */
template <typename must_be_int>
class Timing<inactive_tag, must_be_int>
{
private:
    static const char* empty_string;

public:
    inline void tic(const char* /*tag*/ = NULL) { }
    inline void start() { }
    inline void stop() { }
    inline void reset() { }
    inline point_in_time active_time() { return -1.0; }
    inline point_in_time total_time() { return -1.0; }
    inline const char * c_str() { return empty_string; }
    inline void print() { }
};

template <typename must_be_int>
const char* Timing<inactive_tag, must_be_int>::empty_string = "";

} // namespace parallel

STXXL_END_NAMESPACE

#endif // !STXXL_PARALLEL_TIMING_HEADER
