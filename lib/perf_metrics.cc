//
// This file is part of khmer, https://github.com/dib-lab/khmer/, and is
// Copyright (C) Michigan State University, 2009-2015. It is licensed under
// the three-clause BSD license; see LICENSE.
// Contact: khmer-project@idyll.org
//

#include "perf_metrics.hh"

namespace khmer
{

#ifdef WITH_INTERNAL_METRICS
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
#endif
} // namespace khmer

// vim: set ft=cpp sts=4 sw=4 tw=79:
