//
// This file is part of khmer, http://github.com/ged-lab/khmer/, and is
// Copyright (C) Michigan State University, 2009-2013. It is licensed under
// the three-clause BSD license; see doc/LICENSE.txt. Contact: ctb@msu.edu
//

#include "perf_metrics.hh"

namespace khmer
{


IPerformanceMetrics::
IPerformanceMetrics( )
{ }


IPerformanceMetrics::
~IPerformanceMetrics( )
{ }


uint64_t const
IPerformanceMetrics::
_timespec_diff_in_nsecs( timespec const &start, timespec const &stop )
{
    return
	    ((stop.tv_sec * 1000000000U) + (uint64_t)stop.tv_nsec)
	-   ((start.tv_sec * 1000000000U) + (uint64_t)start.tv_nsec);
}

} // namespace khmer

// vim: set ft=cpp sts=4 sw=4 tw=79:
