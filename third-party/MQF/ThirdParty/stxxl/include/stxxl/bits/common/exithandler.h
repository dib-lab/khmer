/***************************************************************************
 *  include/stxxl/bits/common/exithandler.h
 *
 *  Part of the STXXL. See http://stxxl.sourceforge.net
 *
 *  Copyright (C) 2009 Andreas Beckmann <beckmann@cs.uni-frankfurt.de>
 *
 *  Distributed under the Boost Software License, Version 1.0.
 *  (See accompanying file LICENSE_1_0.txt or copy at
 *  http://www.boost.org/LICENSE_1_0.txt)
 **************************************************************************/

#ifndef STXXL_COMMON_EXITHANDLER_HEADER
#define STXXL_COMMON_EXITHANDLER_HEADER

#include <stxxl/bits/namespace.h>

STXXL_BEGIN_NAMESPACE

// There are several possibilities for the exit handlers.  To use the default
// implementation (which uses atexit()), nothing special has to be done.
//
// To work around problems with atexit() being used in a dll you may #define
// STXXL_NON_DEFAULT_EXIT_HANDLER at library compilation time.  In this case
// the library/application should call stxxl::run_exit_handlers() during
// shutdown.
//
// To provide your own exit handler implementation, #define
// STXXL_EXTERNAL_EXIT_HANDLER and implement stxxl::register_exit_handler(void
// (*)(void)) and stxxl::run_exit_handlers() in your application.

int register_exit_handler(void (* function)(void));
void run_exit_handlers();

STXXL_END_NAMESPACE

#endif // !STXXL_COMMON_EXITHANDLER_HEADER
// vim: et:ts=4:sw=4
