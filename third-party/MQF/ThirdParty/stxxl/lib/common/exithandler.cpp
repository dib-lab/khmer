/***************************************************************************
 *  lib/common/exithandler.cpp
 *
 *  Part of the STXXL. See http://stxxl.sourceforge.net
 *
 *  Copyright (C) 2009 Andreas Beckmann <beckmann@cs.uni-frankfurt.de>
 *
 *  Distributed under the Boost Software License, Version 1.0.
 *  (See accompanying file LICENSE_1_0.txt or copy at
 *  http://www.boost.org/LICENSE_1_0.txt)
 **************************************************************************/

#include <stxxl/bits/common/exithandler.h>
#include <stxxl/bits/namespace.h>

// 1. do nothing for default handler
// 2. #define STXXL_NON_DEFAULT_EXIT_HANDLER for a handler that does not use atexit()
// 3. #define STXXL_EXTERNAL_EXIT_HANDLER to provide your own implementation

#ifndef STXXL_EXTERNAL_EXIT_HANDLER
#ifndef STXXL_NON_DEFAULT_EXIT_HANDLER

#include <cstdlib>

STXXL_BEGIN_NAMESPACE

// default exit handler
int register_exit_handler(void (* function)(void))
{
    return atexit(function);
}

// default exit handler
void run_exit_handlers()
{
    // nothing to do
}

STXXL_END_NAMESPACE

#else // STXXL_NON_DEFAULT_EXIT_HANDLER

#include <vector>
#include <stxxl/bits/common/mutex.h>

STXXL_BEGIN_NAMESPACE

mutex exit_handler_mutex;
std::vector<void (*)(void)> exit_handlers;

int register_exit_handler(void (* function)(void))
{
    scoped_mutex_lock lock(exit_handler_mutex);
    exit_handlers.push_back(function);
    return 0;
}

// default exit handler
void run_exit_handlers()
{
    scoped_mutex_lock lock(exit_handler_mutex);
    while (!exit_handlers.empty()) {
        (*(exit_handlers.back()))();
        exit_handlers.pop_back();
    }
}

STXXL_END_NAMESPACE

#endif // STXXL_NON_DEFAULT_EXIT_HANDLER
#endif // STXXL_EXTERNAL_EXIT_HANDLER

// vim: et:ts=4:sw=4
