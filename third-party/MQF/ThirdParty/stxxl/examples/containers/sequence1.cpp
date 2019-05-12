/***************************************************************************
 *  examples/containers/sequence1.cpp
 *
 *  Part of the STXXL. See http://stxxl.sourceforge.net
 *
 *  Copyright (C) 2013 Daniel Feist <daniel.feist@student.kit.edu>
 *
 *  Distributed under the Boost Software License, Version 1.0.
 *  (See accompanying file LICENSE_1_0.txt or copy at
 *  http://www.boost.org/LICENSE_1_0.txt)
 **************************************************************************/

//! [example]
#include <stxxl/sequence>
#include <iostream>

int main()
{
    typedef stxxl::sequence<int> sequence_type;
    sequence_type my_sequence;

    for (int i = 0; i < 100; ++i)
    {
        my_sequence.push_back(i * i);
    }

    sequence_type::stream forward_stream = my_sequence.get_stream();

    while (!forward_stream.empty())
    {
        std::cout << *forward_stream << " ";
        my_sequence.pop_back();
        ++forward_stream;
    }

    return 0;
}
//! [example]
