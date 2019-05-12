/***************************************************************************
 *  examples/common/cmdline.cpp
 *
 *  Part of the STXXL. See http://stxxl.sourceforge.net
 *
 *  Copyright (C) 2013 Timo Bingmann <tb@panthema.net>
 *
 *  Distributed under the Boost Software License, Version 1.0.
 *  (See accompanying file LICENSE_1_0.txt or copy at
 *  http://www.boost.org/LICENSE_1_0.txt)
 **************************************************************************/

// [example]
#include <stxxl/cmdline>

int main(int argc, char* argv[])
{
    stxxl::cmdline_parser cp;

    // add description and author
    cp.set_description("This may some day be a useful program, which solves "
                       "many serious problems of the real world and achives "
                       "global peace.");
    cp.set_author("Timo Bingmann <tb@panthema.net>");

    // add an unsigned integer option --rounds <N>
    unsigned int rounds = 0;
    cp.add_uint('r', "rounds", "N", rounds,
                "Run N rounds of the experiment.");

    // add a byte size argument which the user can enter like '1gi'
    stxxl::uint64 a_size = 0;
    cp.add_bytes('s', "size", a_size,
                 "Number of bytes to process.");

    // add a required parameter
    std::string a_filename;
    cp.add_param_string("filename", a_filename,
                        "A filename to process");

    // process command line
    if (!cp.process(argc, argv))
        return -1; // some error occurred and help was always written to user.

    std::cout << "Command line parsed okay." << std::endl;

    // output for debugging
    cp.print_result();

    // do something useful
}
// [example]
