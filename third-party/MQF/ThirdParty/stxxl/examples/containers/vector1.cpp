/***************************************************************************
 *  examples/containers/vector1.cpp
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
#include <stxxl/vector>
#include <iostream>

int main()
{
    typedef stxxl::VECTOR_GENERATOR<int>::result vector;
    vector my_vector;

    for (int i = 0; i < 1024 * 1024; i++)
    {
        my_vector.push_back(i + 2);
    }

    std::cout << my_vector[99] << std::endl;
    my_vector[100] = 0;

    while (!my_vector.empty())
    {
        my_vector.pop_back();
    }

    return 0;
}
//! [example]
