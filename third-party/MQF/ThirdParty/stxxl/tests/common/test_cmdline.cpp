/***************************************************************************
 *  tests/common/test_cmdline.cpp
 *
 *  Part of the STXXL. See http://stxxl.sourceforge.net
 *
 *  Copyright (C) 2013 Timo Bingmann <tb@panthema.net>
 *
 *  Distributed under the Boost Software License, Version 1.0.
 *  (See accompanying file LICENSE_1_0.txt or copy at
 *  http://www.boost.org/LICENSE_1_0.txt)
 **************************************************************************/

#include <stxxl/bits/common/cmdline.h>
#include <stxxl/bits/verbose.h>

#include <sstream>

void test1()
{
    int a_int = 0;
    std::string a_str;

    stxxl::cmdline_parser cp;
    cp.add_int('i', "int", "<N>", a_int, "an integer");
    cp.add_string('f', "filename", "<F>", a_str, "a filename");

    cp.set_description("Command Line Parser Test");
    cp.set_author("Timo Bingmann <tb@panthema.net>");

    // good command line
    const char* cmdline1[] =
    { "test", "-i", "42", "-f", "somefile", NULL };

    std::ostringstream os1;
    STXXL_CHECK(cp.process(5, cmdline1, os1));

    STXXL_CHECK(a_int == 42);
    STXXL_CHECK(a_str == "somefile");

    // bad command line
    const char* cmdline2[] =
    { "test", "-i", "dd", "-f", "somefile", NULL };

    std::ostringstream os2;
    STXXL_CHECK(!cp.process(5, cmdline2, os2));
}

int main(int, char**)
{
    test1();
    return 0;
}
