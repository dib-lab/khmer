/***************************************************************************
 *  include/stxxl/bits/parallel/compiletime_settings.h
 *
 *  Defines on options concerning debugging and performance, at compile-time.
 *  Extracted from MCSTL - http://algo2.iti.uni-karlsruhe.de/singler/mcstl/
 *
 *  Part of the STXXL. See http://stxxl.sourceforge.net
 *
 *  Copyright (C) 2007 Johannes Singler <singler@ira.uka.de>
 *  Copyright (C) 2014 Timo Bingmann <tb@panthema.net>
 *
 *  Distributed under the Boost Software License, Version 1.0.
 *  (See accompanying file LICENSE_1_0.txt or copy at
 *  http://www.boost.org/LICENSE_1_0.txt)
 **************************************************************************/

#ifndef STXXL_PARALLEL_COMPILETIME_SETTINGS_HEADER
#define STXXL_PARALLEL_COMPILETIME_SETTINGS_HEADER

#include <stxxl/bits/namespace.h>
#include <stxxl/bits/config.h>
#include <stxxl/bits/parallel/settings.h>
#include <cstdio>

STXXL_BEGIN_NAMESPACE

namespace parallel {

/** STXXL_PARALLEL_PCALL Macro to produce log message when entering a
 *  function.
 *  \param n Input size.
 *  \see STXXL_VERBOSE_LEVEL */
#if (STXXL_VERBOSE_LEVEL <= 0)
#define STXXL_PARALLEL_PCALL(n)
#endif

#if (STXXL_VERBOSE_LEVEL >= 1)
#if STXXL_PARALLEL
#define STXXL_PARALLEL_PCALL(n)                        \
    STXXL_MSG("   " << __FUNCTION__ << ":\n"           \
              "iam = " << omp_get_thread_num() << ", " \
              "n = " << (n) << ", "                    \
              "num_threads = " << SETTINGS::num_threads);
#else
#define STXXL_PARALLEL_PCALL(n)              \
    STXXL_MSG("   " << __FUNCTION__ << ":\n" \
              "iam = single-threaded, "      \
              "n = " << (n) << ", "          \
              "num_threads = " << SETTINGS::num_threads);
#endif
#endif

/** First copy the data, sort it locally, and merge it back (0); or copy it
 * back after everyting is done (1).
 *
 *  Recommendation: 0 */
#define STXXL_MULTIWAY_MERGESORT_COPY_LAST 0

} // namespace parallel

STXXL_END_NAMESPACE

#endif // !STXXL_PARALLEL_COMPILETIME_SETTINGS_HEADER
