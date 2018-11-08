/***************************************************************************
 *  examples/algo/sort_file.cpp
 *
 *  Part of the STXXL. See http://stxxl.sourceforge.net
 *
 *  Copyright (C) 2002-2003 Roman Dementiev <dementiev@mpi-sb.mpg.de>
 *  Copyright (C) 2009 Andreas Beckmann <beckmann@cs.uni-frankfurt.de>
 *
 *  Distributed under the Boost Software License, Version 1.0.
 *  (See accompanying file LICENSE_1_0.txt or copy at
 *  http://www.boost.org/LICENSE_1_0.txt)
 **************************************************************************/

//! \example algo/sort_file.cpp
//! This example imports a file into an \c stxxl::vector without copying its
//! content and then sorts it using stxxl::sort / stxxl::ksort / ...

#include <stxxl/io>
#include <stxxl/mng>
#include <stxxl/ksort>
#include <stxxl/sort>
#include <stxxl/stable_ksort>
#include <stxxl/vector>

struct my_type
{
    typedef unsigned key_type;

    key_type m_key;
    char m_data[128 - sizeof(key_type)];

    key_type key() const
    {
        return m_key;
    }

    my_type() { }
    my_type(key_type k) : m_key(k) { }

    static my_type min_value()
    {
        return my_type(std::numeric_limits<key_type>::min());
    }
    static my_type max_value()
    {
        return my_type(std::numeric_limits<key_type>::max());
    }
};

inline bool operator < (const my_type& a, const my_type& b)
{
    return a.key() < b.key();
}

inline bool operator == (const my_type& a, const my_type& b)
{
    return a.key() == b.key();
}

struct Cmp
{
    typedef my_type first_argument_type;
    typedef my_type second_argument_type;
    typedef bool result_type;
    bool operator () (const my_type& a, const my_type& b) const
    {
        return a < b;
    }
    static my_type min_value()
    {
        return my_type::min_value();
    }
    static my_type max_value()
    {
        return my_type::max_value();
    }
};

std::ostream& operator << (std::ostream& o, const my_type& obj)
{
    o << obj.key();
    return o;
}

int main(int argc, char** argv)
{
    if (argc < 3)
    {
        std::cout << "Usage: " << argv[0] << " action file" << std::endl;
        std::cout << "       where action is one of generate, sort, ksort, stable_sort, stable_ksort" << std::endl;
        return -1;
    }

    const stxxl::unsigned_type block_size = sizeof(my_type) * 4096;

    if (strcmp(argv[1], "generate") == 0) {
        const my_type::key_type num_elements = 1 * 1024 * 1024;
        const stxxl::unsigned_type records_in_block = block_size / sizeof(my_type);
        stxxl::syscall_file f(argv[2], stxxl::file::CREAT | stxxl::file::RDWR);
        my_type* array = (my_type*)stxxl::aligned_alloc<STXXL_BLOCK_ALIGN>(block_size);
        memset(array, 0, block_size);

        my_type::key_type cur_key = num_elements;
        for (unsigned i = 0; i < num_elements / records_in_block; i++)
        {
            for (unsigned j = 0; j < records_in_block; j++)
                array[j].m_key = cur_key--;

            stxxl::request_ptr req = f.awrite((void*)array, stxxl::int64(i) * block_size, block_size);
            req->wait();
        }
        stxxl::aligned_dealloc<STXXL_BLOCK_ALIGN>(array);
    }
    else {
#if STXXL_PARALLEL_MULTIWAY_MERGE
        STXXL_MSG("STXXL_PARALLEL_MULTIWAY_MERGE");
#endif
        stxxl::syscall_file f(argv[2], stxxl::file::DIRECT | stxxl::file::RDWR);
        unsigned memory_to_use = 50 * 1024 * 1024;
        typedef stxxl::vector<my_type, 1, stxxl::lru_pager<8>, block_size> vector_type;
        vector_type v(&f);

        /*
        STXXL_MSG("Printing...");
        for(stxxl::int64 i=0; i < v.size(); i++)
            STXXL_MSG(v[i].key());
         */

        STXXL_MSG("Checking order...");
        STXXL_MSG((stxxl::is_sorted(v.begin(), v.end()) ? "OK" : "WRONG"));

        STXXL_MSG("Sorting...");
        if (strcmp(argv[1], "sort") == 0) {
            stxxl::sort(v.begin(), v.end(), Cmp(), memory_to_use);
#if 0       // stable_sort is not yet implemented
        }
        else if (strcmp(argv[1], "stable_sort") == 0) {
            stxxl::stable_sort(v.begin(), v.end(), memory_to_use);
#endif
        }
        else if (strcmp(argv[1], "ksort") == 0) {
            stxxl::ksort(v.begin(), v.end(), memory_to_use);
        }
        else if (strcmp(argv[1], "stable_ksort") == 0) {
            stxxl::stable_ksort(v.begin(), v.end(), memory_to_use);
        }
        else {
            STXXL_MSG("Not implemented: " << argv[1]);
        }

        STXXL_MSG("Checking order...");
        STXXL_MSG((stxxl::is_sorted(v.begin(), v.end()) ? "OK" : "WRONG"));
    }

    return 0;
}

// vim: et:ts=4:sw=4
