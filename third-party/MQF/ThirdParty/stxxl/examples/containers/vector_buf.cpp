/***************************************************************************
 *  examples/containers/vector_buf.cpp
 *
 *  Part of the STXXL. See http://stxxl.sourceforge.net
 *
 *  Copyright (C) 2013 Timo Bingmann <tb@panthema.net>
 *
 *  Distributed under the Boost Software License, Version 1.0.
 *  (See accompanying file LICENSE_1_0.txt or copy at
 *  http://www.boost.org/LICENSE_1_0.txt)
 **************************************************************************/

#include <iostream>
#include <stxxl/vector>
#include <stxxl/bits/config.h>

using stxxl::uint64;

void test_vector_element(uint64 size)
{
    stxxl::scoped_print_timer tm("vector element access", 2 * size * sizeof(uint64));

//! [element]
    typedef stxxl::VECTOR_GENERATOR<uint64>::result vector_type;

    vector_type vec(size);

    for (uint64 i = 0; i < vec.size(); ++i)
        vec[i] = (i % 1024);

    uint64 sum = 0;
    for (uint64 i = 0; i < vec.size(); ++i)
        sum += vec[i];
//! [element]

    std::cout << "sum: " << sum << std::endl;
    STXXL_CHECK(sum == size / 1024 * (1024 * 1023 / 2));
}

void test_vector_iterator(uint64 size)
{
    stxxl::scoped_print_timer tm("vector iterator access", 2 * size * sizeof(uint64));

//! [iterator]
    typedef stxxl::VECTOR_GENERATOR<uint64>::result vector_type;

    vector_type vec(size);

    uint64 i = 0;
    for (vector_type::iterator it = vec.begin(); it != vec.end(); ++it, ++i)
        *it = (i % 1024);

    uint64 sum = 0;
    for (vector_type::const_iterator it = vec.begin(); it != vec.end(); ++it)
        sum += *it;
//! [iterator]

    std::cout << "sum: " << sum << std::endl;
    STXXL_CHECK(sum == size / 1024 * (1024 * 1023 / 2));
}

void test_vector_buffered(uint64 size)
{
    stxxl::scoped_print_timer tm("vector buffered access", 2 * size * sizeof(uint64));

//! [buffered]
    typedef stxxl::VECTOR_GENERATOR<uint64>::result vector_type;

    vector_type vec(size);

    // write using vector_bufwriter
    vector_type::bufwriter_type writer(vec);

    for (uint64 i = 0; i < vec.size(); ++i)
        writer << (i % 1024);

    // required to flush out the last block (or destruct the bufwriter)
    writer.finish();

    // now read using vector_bufreader
    uint64 sum = 0;

    for (vector_type::bufreader_type reader(vec); !reader.empty(); ++reader)
    {
        sum += *reader;
    }
//! [buffered]

    std::cout << "sum: " << sum << std::endl;
    STXXL_CHECK(sum == size / 1024 * (1024 * 1023 / 2));
}

#if STXXL_HAVE_CXX11_RANGE_FOR_LOOP
void test_vector_cxx11(uint64 size)
{
    stxxl::scoped_print_timer tm("vector C++11 loop access", 2 * size * sizeof(uint64));

    typedef stxxl::VECTOR_GENERATOR<uint64>::result vector_type;

    vector_type vec(size);

    {
        vector_type::bufwriter_type writer(vec);

        for (uint64 i = 0; i < vec.size(); ++i)
            writer << (i % 1024);
    }

//! [cxx11]
    // now using vector_bufreader adaptor to C++11 for loop
    uint64 sum = 0;

    for (auto it : vector_type::bufreader_type(vec))
    {
        sum += it;
    }
//! [cxx11]

    std::cout << "sum: " << sum << std::endl;
    STXXL_CHECK(sum == size / 1024 * (1024 * 1023 / 2));
}
#endif

int main(int argc, char* argv[])
{
    int multi = (argc >= 2 ? atoi(argv[1]) : 64);
    const uint64 size = multi * 1024 * uint64(1024) / sizeof(uint64);

    stxxl::block_manager::get_instance();

    test_vector_element(size);
    test_vector_iterator(size);
    test_vector_buffered(size);

#if STXXL_HAVE_CXX11_RANGE_FOR_LOOP
    test_vector_cxx11(size);
#endif

    return 0;
}
