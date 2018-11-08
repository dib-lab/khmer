/***************************************************************************
 *  tests/containers/test_vector_buf.cpp
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
#include <algorithm>
#include <stxxl/vector>
#include <stxxl/scan>
#include <stxxl/stream>

using stxxl::uint64;

struct my_type  // 24 bytes, not a power of 2 intentionally
{
    uint64 key;
    uint64 load0;
    uint64 load1;

    my_type(uint64 i = 0)
        : key(i),
          load0(i + 1),
          load1(1 + 42)
    { }

    bool operator == (const my_type& b) const
    {
        return (key == b.key) && (load0 == b.load0) && (load1 == b.load1);
    }
};

//! Verify contents of the vector
template <typename VectorType>
void check_vector(const VectorType& v)
{
    typedef typename VectorType::value_type value_type;

    for (uint64 i = 0; i < v.size(); ++i)
    {
        STXXL_CHECK(v[i] == value_type(i));
    }
}

//! Stream object generating lots of ValueTypes
template <typename ValueType>
class MyStream
{
    uint64 i;

public:
    typedef ValueType value_type;

    MyStream()
        : i(0)
    { }

    value_type operator * () const
    {
        return value_type(i);
    }

    MyStream& operator ++ ()
    {
        ++i;
        return *this;
    }

    bool empty() const
    {
        return false;
    }
};

template <typename ValueType>
void test_vector_buf(uint64 size)
{
    typedef typename stxxl::VECTOR_GENERATOR<ValueType>::result vector_type;

    typedef typename vector_type::iterator vector_iterator_type;

    {   // fill vector using element access
        stxxl::scoped_print_timer tm("element access");

        vector_type vec(size);

        for (uint64 i = 0; i < size; ++i)
            vec[i] = ValueType(i);

        check_vector(vec);
    }
    {   // fill vector using iterator access
        stxxl::scoped_print_timer tm("iterator access");

        vector_type vec(size);

        vector_iterator_type vi = vec.begin();

        for (uint64 i = 0; i < size; ++i, ++vi)
            *vi = ValueType(i);

        check_vector(vec);
    }
    {   // fill vector using vector_bufwriter
        stxxl::scoped_print_timer tm("vector_bufwriter");

        vector_type vec(size);

        typename vector_type::bufwriter_type writer(vec.begin());

        for (uint64 i = 0; i < size; ++i)
            writer << ValueType(i);

        writer.finish();

        check_vector(vec);
    }
    {   // fill empty vector using vector_bufwriter
        stxxl::scoped_print_timer tm("empty vector_bufwriter");

        vector_type vec;

        typename vector_type::bufwriter_type writer(vec);

        for (uint64 i = 0; i < size; ++i)
            writer << ValueType(i);

        writer.finish();

        check_vector(vec);
    }

    vector_type vec(size);

    {   // fill vector using materialize
        stxxl::scoped_print_timer tm("materialize");

        MyStream<ValueType> stream;
        stxxl::stream::materialize(stream, vec.begin(), vec.end());

        check_vector(vec);
    }
    {   // read vector using vector_bufreader
        stxxl::scoped_print_timer tm("vector_bufreader");

        const vector_type& cvec = vec;

        typename vector_type::bufreader_type reader(cvec.begin(), cvec.end());

        for (uint64 i = 0; i < size; ++i)
        {
            STXXL_CHECK(!reader.empty());
            STXXL_CHECK(reader.size() == size - i);

            ValueType pv = *reader;

            ValueType v;
            reader >> v;

            STXXL_CHECK(v == ValueType(i));
            STXXL_CHECK(pv == v);
            STXXL_CHECK(reader.size() == size - i - 1);
        }

        STXXL_CHECK(reader.empty());

        // rewind reader and read again
        reader.rewind();

        for (uint64 i = 0; i < size; ++i)
        {
            STXXL_CHECK(!reader.empty());
            STXXL_CHECK(reader.size() == size - i);

            ValueType pv = *reader;

            ValueType v;
            reader >> v;

            STXXL_CHECK(v == ValueType(i));
            STXXL_CHECK(pv == v);
            STXXL_CHECK(reader.size() == size - i - 1);
        }

        STXXL_CHECK(reader.empty());
    }
    {   // read vector using vector_bufreader_reverse
        stxxl::scoped_print_timer tm("vector_bufreader_reverse");

        const vector_type& cvec = vec;

        typename vector_type::bufreader_reverse_type reader(cvec.begin(), cvec.end());

        for (uint64 i = 0; i < size; ++i)
        {
            STXXL_CHECK(!reader.empty());
            STXXL_CHECK(reader.size() == size - i);

            ValueType pv = *reader;

            ValueType v;
            reader >> v;

            STXXL_CHECK(v == ValueType(size - i - 1));
            STXXL_CHECK(pv == v);
            STXXL_CHECK(reader.size() == size - i - 1);
        }

        STXXL_CHECK(reader.empty());

        // rewind reader and read again
        reader.rewind();

        for (uint64 i = 0; i < size; ++i)
        {
            STXXL_CHECK(!reader.empty());
            STXXL_CHECK(reader.size() == size - i);

            ValueType pv = *reader;

            ValueType v;
            reader >> v;

            STXXL_CHECK(v == ValueType(size - i - 1));
            STXXL_CHECK(pv == v);
            STXXL_CHECK(reader.size() == size - i - 1);
        }

        STXXL_CHECK(reader.empty());
    }
#if STXXL_HAVE_CXX11_RANGE_FOR_LOOP
    {   // read vector using C++11 for loop construct
        stxxl::scoped_print_timer tm("C++11 for loop");

        uint64 i = 0;

        for (auto it : vec)
        {
            STXXL_CHECK(it == ValueType(i));
            ++i;
        }

        STXXL_CHECK(i == vec.size());
    }
    {   // read vector using C++11 for loop construct
        stxxl::scoped_print_timer tm("C++11 bufreader for loop");
        typedef typename vector_type::bufreader_type bufreader_type;

        uint64 i = 0;

        for (auto it : bufreader_type(vec))
        {
            STXXL_CHECK(it == ValueType(i));
            ++i;
        }

        STXXL_CHECK(i == vec.size());
    }
#endif
}

int main(int argc, char* argv[])
{
    int size = (argc > 1) ? atoi(argv[1]) : 16;

    STXXL_MSG("Testing stxxl::vector<int> with even size");
    test_vector_buf<int>(size * 1024 * 1024);

    STXXL_MSG("Testing stxxl::vector<int> with odd size");
    test_vector_buf<int>(size * 1024 * 1024 + 501 + 42);

    STXXL_MSG("Testing stxxl::vector<uint64>");
    test_vector_buf<uint64>(size * 1024 * 1024 + 501 + 42);

    STXXL_MSG("Testing stxxl::vector<my_type>");
    test_vector_buf<my_type>(size * 1024 * 1024);

    return 0;
}
