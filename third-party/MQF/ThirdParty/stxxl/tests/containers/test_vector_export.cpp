/***************************************************************************
 *  tests/containers/test_vector_export.cpp
 *
 *  Part of the STXXL. See http://stxxl.sourceforge.net
 *
 *  Copyright (C) 2008 Johannes Singler <singler@ira.uka.de>
 *
 *  Distributed under the Boost Software License, Version 1.0.
 *  (See accompanying file LICENSE_1_0.txt or copy at
 *  http://www.boost.org/LICENSE_1_0.txt)
 **************************************************************************/

//! \example containers/test_vector_export.cpp
//! This is an example of use of \c stxxl::vector::export_files

#include <iostream>
#include <algorithm>
#include <stxxl/vector>
#include <stxxl/scan>

typedef stxxl::int64 int64;

int main()
{
    // use non-randomized striping to avoid side effects on random generator
    typedef stxxl::VECTOR_GENERATOR<int64, 2, 2, (2* 1024* 1024), stxxl::striping>::result vector_type;
    vector_type v(int64(64 * 1024 * 1024) / sizeof(int64));

    stxxl::random_number32 rnd;
    int offset = rnd();

    STXXL_MSG("write " << v.size() << " elements");

    stxxl::ran32State = 0xdeadbeef;
    vector_type::size_type i;

    // fill the vector with increasing sequence of integer numbers
    for (i = 0; i < v.size(); ++i)
    {
        v[i] = i + offset;
        STXXL_CHECK(v[i] == int64(i + offset));
    }

    v.flush();

    STXXL_MSG("export files");
    v.export_files("exported_");

    return 0;
}
