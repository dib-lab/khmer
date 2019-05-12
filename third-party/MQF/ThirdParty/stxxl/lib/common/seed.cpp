/***************************************************************************
 *  lib/common/seed.cpp
 *
 *  Part of the STXXL. See http://stxxl.sourceforge.net
 *
 *  Copyright (C) 2007, 2008 Andreas Beckmann <beckmann@cs.uni-frankfurt.de>
 *  Copyright (C) 2008 Roman Dementiev <dementiev@ira.uka.de>
 *
 *  Distributed under the Boost Software License, Version 1.0.
 *  (See accompanying file LICENSE_1_0.txt or copy at
 *  http://www.boost.org/LICENSE_1_0.txt)
 **************************************************************************/

#include <cassert>
#include <stxxl/bits/config.h>
#include <stxxl/bits/common/mutex.h>
#include <stxxl/bits/namespace.h>

#if STXXL_WINDOWS
  #ifndef NOMINMAX
    #define NOMINMAX
  #endif
  #include <io.h>
  #include <windows.h>
#else
  #include <unistd.h>
  #include <sys/time.h>
#endif

STXXL_BEGIN_NAMESPACE

inline unsigned initial_seed();

struct seed_generator_t {
    unsigned seed;
    fastmutex mtx;

    seed_generator_t(unsigned s) : seed(s)
    { }
};

seed_generator_t & seed_generator()
{
    static seed_generator_t sg(initial_seed());
    return sg;
}

inline unsigned initial_seed()
{
#ifndef NDEBUG
    static bool initialized = false;
    assert(!initialized); // this should only be called once!

    initialized = true;
#endif // NDEBUG
#if STXXL_WINDOWS
    // GetTickCount():  ms since system start
    return GetTickCount() ^ GetCurrentProcessId();
#else
    struct timeval tv;
    gettimeofday(&tv, 0);

    return (unsigned)(tv.tv_sec ^ tv.tv_usec ^ (getpid() << 16));
#endif
}

void set_seed(unsigned seed)
{
    scoped_fast_mutex_lock Lock(seed_generator().mtx);
    seed_generator().seed = seed;
}

unsigned get_next_seed()
{
    scoped_fast_mutex_lock Lock(seed_generator().mtx);
    return seed_generator().seed++;
}

STXXL_END_NAMESPACE
