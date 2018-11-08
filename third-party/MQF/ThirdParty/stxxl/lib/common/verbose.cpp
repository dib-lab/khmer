/***************************************************************************
 *  lib/common/verbose.cpp
 *
 *  Part of the STXXL. See http://stxxl.sourceforge.net
 *
 *  Copyright (C) 2009, 2010 Andreas Beckmann <beckmann@cs.uni-frankfurt.de>
 *
 *  Distributed under the Boost Software License, Version 1.0.
 *  (See accompanying file LICENSE_1_0.txt or copy at
 *  http://www.boost.org/LICENSE_1_0.txt)
 **************************************************************************/

#include <iostream>
#include <cstdio>
#include <cmath>
#include <stxxl/bits/verbose.h>
#include <stxxl/bits/common/log.h>
#include <stxxl/bits/common/timer.h>
#include <stxxl/bits/msvc_compatibility.h>

#ifndef STXXL_THREAD_ID
# if STXXL_STD_THREADS || STXXL_BOOST_THREADS
#  define STXXL_THREAD_ID (-1)
# else
#  define STXXL_THREAD_ID pthread_self()
# endif
#endif

STXXL_BEGIN_NAMESPACE

static const double program_start_time_stamp = timestamp();

void print_msg(const char* label, const std::string& msg, unsigned flags)
{
    std::string s;
#ifdef STXXL_PRINT_TIMESTAMP_ALWAYS
    const bool timestamp_always = true;
#else
    const bool timestamp_always = false;
#endif
    if (timestamp_always || (flags & _STXXL_PRNT_TIMESTAMP)) {
        double t = timestamp() - program_start_time_stamp;
        char tstr[23]; /* "[364:23:59:59.999999] " */
        snprintf(tstr, sizeof(tstr), "[%d.%02d:%02d:%02d.%06d] ",
                 int(t / (24 * 60 * 60)),
                 int(t / (60 * 60)) % 24,
                 int(t / 60) % 60, int(t) % 60,
                 int((t - floor(t)) * 1000000));
        s += tstr;
    }
    if (label) {
        s += '[';
        s += label;
        s += "] ";
    }
    if (flags & _STXXL_PRNT_THREAD_ID) {
        char tstr[32];
        snprintf(tstr, sizeof(tstr), "[T%ld] ", long(STXXL_THREAD_ID));
        s += tstr;
    }
    s += msg;
    if (flags & _STXXL_PRNT_ADDNEWLINE)
        s += '\n';
    if (flags & _STXXL_PRNT_COUT)
        std::cout << s << std::flush;
    if (flags & _STXXL_PRNT_CERR)
        std::cerr << s << std::flush;
    logger* logger_instance = logger::get_instance();
    if (flags & _STXXL_PRNT_LOG)
        logger_instance->log_stream() << s << std::flush;
    if (flags & _STXXL_PRNT_ERRLOG)
        logger_instance->errlog_stream() << s << std::flush;
}

STXXL_END_NAMESPACE

// vim: et:ts=4:sw=4
