/***************************************************************************
 *  tools/mlock.cpp
 *
 *  Allocate some memory and mlock() it to consume physical memory.
 *  Needs to run as root to block more than 64 KiB in default settings.
 *
 *  Part of the STXXL. See http://stxxl.sourceforge.net
 *
 *  Copyright (C) 2010 Andreas Beckmann <beckmann@cs.uni-frankfurt.de>
 *  Copyright (C) 2013 Timo Bingmann <tb@panthema.net>
 *
 *  Distributed under the Boost Software License, Version 1.0.
 *  (See accompanying file LICENSE_1_0.txt or copy at
 *  http://www.boost.org/LICENSE_1_0.txt)
 **************************************************************************/

#include <stxxl/bits/config.h>
#include <stxxl/bits/verbose.h>

#if STXXL_HAVE_MLOCK_PROTO

#include <cerrno>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <sys/mman.h>
#include <unistd.h>

#include <stxxl/cmdline>

int do_mlock(int argc, char* argv[])
{
    // parse command line
    stxxl::cmdline_parser cp;

    cp.set_description(
        "Allocate some memory and mlock() it to consume physical memory. "
        "Needs to run as root to block more than 64 KiB in default settings.");
    cp.set_author("Andreas Beckmann <beckmann@cs.uni-frankfurt.de>");

    stxxl::unsigned_type M;
    cp.add_param_bytes("size", M,
                       "Amount of memory to allocate (e.g. 4GiB)");

    if (!cp.process(argc, argv))
        return -1;

    // allocate and fill
    char* c = (char*)malloc(M);
    memset(c, 42, M);

    if (mlock(c, M) == 0)
    {
        std::cout << "mlock(" << (void*)c << ", " << M << ") successful, press Ctrl-C to exit." << std::endl;
        while (1)
            sleep(86400);
    }
    else {
        std::cerr << "mlock(" << (void*)c << ", " << M << ") failed: " << strerror(errno) << std::endl;
        return 1;
    }
}

#else // !STXXL_HAVE_MLOCK_PROTO

int do_mlock(int, char*[])
{
    STXXL_MSG("Sorry, mlock() is not supported on this platform.");
    return -1;
}

#endif // STXXL_HAVE_MLOCK_PROTO

// vim: et:ts=4:sw=4
