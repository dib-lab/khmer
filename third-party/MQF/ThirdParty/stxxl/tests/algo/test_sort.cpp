/***************************************************************************
 *  tests/algo/test_sort.cpp
 *
 *  Part of the STXXL. See http://stxxl.sourceforge.net
 *
 *  Copyright (C) 2002 Roman Dementiev <dementiev@mpi-sb.mpg.de>
 *  Copyright (C) 2010 Andreas Beckmann <beckmann@cs.uni-frankfurt.de>
 *
 *  Distributed under the Boost Software License, Version 1.0.
 *  (See accompanying file LICENSE_1_0.txt or copy at
 *  http://www.boost.org/LICENSE_1_0.txt)
 **************************************************************************/

//! \example algo/test_sort.cpp
//! This is an example of how to use \c stxxl::sort() algorithm

#include <stxxl/mng>
#include <stxxl/sort>
#include <stxxl/vector>

#define RECORD_SIZE 8

struct my_type
{
    typedef unsigned key_type;

    key_type m_key;
    char m_data[RECORD_SIZE - sizeof(key_type)];
    key_type key() const
    {
        return m_key;
    }

    my_type() { }
    my_type(key_type k) : m_key(k)
    {
#if STXXL_WITH_VALGRIND
        memset(m_data, 0, sizeof(m_data));
#endif
    }

    static my_type min_value()
    {
        return my_type(std::numeric_limits<key_type>::min());
    }
    static my_type max_value()
    {
        return my_type(std::numeric_limits<key_type>::max());
    }

    ~my_type() { }
};

std::ostream& operator << (std::ostream& o, const my_type& obj)
{
    o << obj.m_key;
    return o;
}

bool operator < (const my_type& a, const my_type& b)
{
    return a.key() < b.key();
}

bool operator == (const my_type& a, const my_type& b)
{
    return a.key() == b.key();
}

bool operator != (const my_type& a, const my_type& b)
{
    return a.key() != b.key();
}

struct cmp : public std::less<my_type>
{
    my_type min_value() const
    {
        return my_type::min_value();
    }
    my_type max_value() const
    {
        return my_type::max_value();
    }
};

int main()
{
#if STXXL_PARALLEL_MULTIWAY_MERGE
    STXXL_MSG("STXXL_PARALLEL_MULTIWAY_MERGE");
#endif
    unsigned memory_to_use = 128 * 1024 * 1024;
    typedef stxxl::vector<my_type> vector_type;

    {
        // test small vector that can be sorted internally
        vector_type v(3);
        v[0] = 42;
        v[1] = 0;
        v[2] = 23;
        STXXL_MSG("small vector unsorted " << v[0] << " " << v[1] << " " << v[2]);
        //stxxl::sort(v.begin(), v.end(), cmp(), memory_to_use);
        stxxl::stl_in_memory_sort(v.begin(), v.end(), cmp());
        STXXL_MSG("small vector sorted   " << v[0] << " " << v[1] << " " << v[2]);
        STXXL_CHECK(stxxl::is_sorted(v.begin(), v.end(), cmp()));
    }

    const stxxl::int64 n_records =
        stxxl::int64(384) * stxxl::int64(1024 * 1024) / sizeof(my_type);
    vector_type v(n_records);

    stxxl::random_number32 rnd;
    STXXL_MSG("Filling vector..., input size = " << v.size() << " elements (" << ((v.size() * sizeof(my_type)) >> 20) << " MiB)");
    for (vector_type::size_type i = 0; i < v.size(); i++)
        v[i].m_key = 1 + (rnd() % 0xfffffff);

    STXXL_MSG("Checking order...");
    STXXL_CHECK(!stxxl::is_sorted(v.begin(), v.end(), cmp()));

    STXXL_MSG("Sorting (using " << (memory_to_use >> 20) << " MiB of memory)...");
    stxxl::sort(v.begin(), v.end(), cmp(), memory_to_use);

    STXXL_MSG("Checking order...");
    STXXL_CHECK(stxxl::is_sorted(v.begin(), v.end(), cmp()));

    STXXL_MSG("Done, output size=" << v.size());

    return 0;
}

// vim: et:ts=4:sw=4
