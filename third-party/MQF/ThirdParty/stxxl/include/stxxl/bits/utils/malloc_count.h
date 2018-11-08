/***************************************************************************
 *  include/stxxl/bits/utils/malloc_count.h
 *
 *  Header containing prototypes of user-callable functions to retrieve
 *  run-time information about malloc()/free() allocation. Used with
 *  malloc_count.c, a light-weight malloc()-based heap allocation counter based
 *  on overriding malloc/free/realloc() functions.
 *
 *  See http://panthema.net/2013/malloc_count/ for a tutorial.
 *
 *  Part of the STXXL. See http://stxxl.sourceforge.net
 *
 *  Copyright (C) 2013-2015 Timo Bingmann <tb@panthema.net>
 *
 *  Distributed under the Boost Software License, Version 1.0.
 *  (See accompanying file LICENSE_1_0.txt or copy at
 *  http://www.boost.org/LICENSE_1_0.txt)
 **************************************************************************/

#ifndef STXXL_UTILS_MALLOC_COUNT_HEADER
#define STXXL_UTILS_MALLOC_COUNT_HEADER

#include <cstdlib>

#ifdef __cplusplus
extern "C" { /* for inclusion from C++ */
#endif

/* returns the currently allocated amount of memory */
extern size_t malloc_count_current(void);

/* returns the current peak memory allocation */
extern size_t malloc_count_peak(void);

/* resets the peak memory allocation to current */
extern void malloc_count_reset_peak(void);

/* typedef of callback function */
typedef void (* malloc_count_callback_type)(void* cookie, size_t current);

/* supply malloc_count with a callback function that is invoked on each change
 * of the current allocation. The callback function must not use
 * malloc()/realloc()/free() or it will go into an endless recursive loop! */
extern void malloc_count_set_callback(malloc_count_callback_type cb,
                                      void* cookie);

/* user function which prints current and peak allocation to stderr */
extern void malloc_count_print_status(void);

#ifdef __cplusplus
} /* extern "C" */
#endif

#endif // !STXXL_UTILS_MALLOC_COUNT_HEADER
