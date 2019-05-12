/***************************************************************************
 *  tests/algo/test_ksort.cpp
 *
 *  Part of the STXXL. See http://stxxl.sourceforge.net
 *
 *  Copyright (C) 2002 Roman Dementiev <dementiev@mpi-sb.mpg.de>
 *
 *  Distributed under the Boost Software License, Version 1.0.
 *  (See accompanying file LICENSE_1_0.txt or copy at
 *  http://www.boost.org/LICENSE_1_0.txt)
 **************************************************************************/

//! \example algo/test_ksort.cpp
//! This is an example of how to use \c stxxl::ksort() algorithm

#include <stxxl/mng>
#include <stxxl/ksort>
#include <stxxl/vector>

struct my_type
{
    typedef stxxl::uint64 key_type1;

    key_type1 m_key;
    key_type1 m_key_copy;
    char m_data[32 - 2 * sizeof(key_type1)];

    key_type1 key() const
    {
        return m_key;
    }

    my_type() : m_key(0), m_key_copy(0) { }
    my_type(key_type1 k) : m_key(k), m_key_copy(k) { }

    my_type min_value() const { return my_type(std::numeric_limits<key_type1>::min()); }
    my_type max_value() const { return my_type(std::numeric_limits<key_type1>::max()); }
};

std::ostream& operator << (std::ostream& o, const my_type& obj)
{
    o << obj.m_key << " " << obj.m_key_copy;
    return o;
}

struct get_key
{
    typedef my_type::key_type1 key_type;
    my_type dummy;
    key_type operator () (const my_type& obj) const
    {
        return obj.m_key;
    }
    my_type min_value() const
    {
        return my_type((dummy.min_value)());
    }
    my_type max_value() const
    {
        return my_type((dummy.max_value)());
    }
};

bool operator < (const my_type& a, const my_type& b)
{
    return a.key() < b.key();
}

int main()
{
    STXXL_MSG("Check config...");
    try {
        stxxl::block_manager::get_instance();
    }
    catch (std::exception& e)
    {
        STXXL_MSG("Exception: " << e.what());
        abort();
    }
#if STXXL_PARALLEL_MULTIWAY_MERGE
    STXXL_MSG("STXXL_PARALLEL_MULTIWAY_MERGE");
#endif
    unsigned memory_to_use = 32 * 1024 * 1024;
    typedef stxxl::VECTOR_GENERATOR<my_type, 4, 4>::result vector_type;
    const stxxl::int64 n_records = 3 * 32 * stxxl::int64(1024 * 1024) / sizeof(my_type);
    vector_type v(n_records);

    stxxl::random_number32 rnd;
    STXXL_MSG("Filling vector... ");
    for (vector_type::size_type i = 0; i < v.size(); i++)
    {
        v[i].m_key = rnd() + 1;
        v[i].m_key_copy = v[i].m_key;
    }

    STXXL_MSG("Checking order...");
    STXXL_CHECK(!stxxl::is_sorted(v.begin(), v.end()));

    STXXL_MSG("Sorting...");
    stxxl::ksort(v.begin(), v.end(), get_key(), memory_to_use);
    //stxxl::ksort(v.begin(),v.end(),memory_to_use);

    STXXL_MSG("Checking order...");
    STXXL_CHECK(stxxl::is_sorted(v.begin(), v.end()));
    STXXL_MSG("Checking content...");
    my_type prev;
    for (vector_type::size_type i = 0; i < v.size(); i++)
    {
        if (v[i].m_key != v[i].m_key_copy)
        {
            STXXL_MSG("Bug at position " << i);
            abort();
        }
        if (i > 0 && prev.m_key == v[i].m_key)
        {
            STXXL_MSG("Duplicate at position " << i << " key=" << v[i].m_key);
            //abort();
        }
        prev = v[i];
    }
    STXXL_MSG("OK");

    return 0;
}
