/***************************************************************************
 *  tests/containers/hash_map/test_hash_map_iterators.cpp
 *
 *  Part of the STXXL. See http://stxxl.sourceforge.net
 *
 *  Copyright (C) 2007 Markus Westphal <marwes@users.sourceforge.net>
 *  Copyright (C) 2014 Timo Bingmann <tb@panthema.net>
 *
 *  Distributed under the Boost Software License, Version 1.0.
 *  (See accompanying file LICENSE_1_0.txt or copy at
 *  http://www.boost.org/LICENSE_1_0.txt)
 **************************************************************************/

#include <iostream>

#include <stxxl.h>
#include <stxxl/bits/common/seed.h>
#include <stxxl/bits/common/rand.h>
#include <stxxl/bits/containers/hash_map/hash_map.h>
#include <stxxl/bits/compat/hash_map.h>

using stxxl::unsigned_type;

struct rand_pairs
{
    stxxl::random_number32& rand_;

    rand_pairs(stxxl::random_number32& rand)
        : rand_(rand)
    { }

    std::pair<int, int> operator () ()
    {
        int v = rand_();
        return std::pair<int, int>(v, v);
    }
};

struct hash_int
{
    size_t operator () (int key) const
    {
        // a simple integer hash function
        return (size_t)(key * 2654435761u);
    }
};

struct cmp : public std::less<int>
{
    int min_value() const { return std::numeric_limits<int>::min(); }
    int max_value() const { return std::numeric_limits<int>::max(); }
};

////////////////////////////////////////////////////////////////////////////////
void cmp_with_internal_map()
{
    typedef std::pair<int, int> value_type;
    const unsigned_type value_size = sizeof(value_type);

    const unsigned_type n_values = 15000;
    const unsigned_type n_tests = 7500;

    // make sure all changes will be buffered
    const unsigned_type buffer_size = 5 * n_values * (value_size + sizeof(int*));
    const unsigned_type mem_to_sort = 32 * 1024 * 1024;

    const unsigned_type subblock_raw_size = 4 * 1024;
    const unsigned_type block_size = 4;

    typedef stxxl::hash_map::hash_map<int, int, hash_int, cmp,
                                      subblock_raw_size, block_size> hash_map;
    typedef hash_map::const_iterator const_iterator;

    typedef stxxl::compat_hash_map<int, int>::result int_hash_map;

    stxxl::stats_data stats_begin = *stxxl::stats::get_instance();

    hash_map map;
    map.max_buffer_size(buffer_size);
    const hash_map& cmap = map;
    int_hash_map int_map;

    // generate random values
    stxxl::random_number32 rand32;
    std::vector<value_type> values1(n_values);
    std::vector<value_type> values2(n_values);
    std::vector<value_type> values3(n_values);
    std::generate(values1.begin(), values1.end(), rand_pairs(rand32) _STXXL_FORCE_SEQUENTIAL);
    std::generate(values2.begin(), values2.end(), rand_pairs(rand32) _STXXL_FORCE_SEQUENTIAL);
    std::generate(values3.begin(), values3.end(), rand_pairs(rand32) _STXXL_FORCE_SEQUENTIAL);

    // --- initial import: create a nice mix of externally (values1) and
    // --- internally (values2) stored values
    std::cout << "Initial import...";

    map.insert(values1.begin(), values1.end(), mem_to_sort);
    int_map.insert(values1.begin(), values1.end());

    std::vector<value_type>::iterator val_it = values2.begin();
    for ( ; val_it != values2.end(); ++val_it) {
        map.insert_oblivious(*val_it);
        int_map.insert(*val_it);
    }

    // --- erase and overwrite some external values
    std::random_shuffle(values1.begin(), values1.end());
    val_it = values1.begin();
    for ( ; val_it != values1.begin() + n_tests; ++val_it) {
        map.erase_oblivious(val_it->first);
        int_map.erase(val_it->first);
    }
    for ( ; val_it != values1.begin() + 2 * n_tests; ++val_it) {
        map.insert_oblivious(*val_it);
        int_map.insert(*val_it);
    }

    // --- scan and compare with internal memory hash-map
    std::cout << "Compare with internal-memory map...";
    STXXL_CHECK(int_map.size() == map.size());
    const_iterator cit = cmap.begin();
    for ( ; cit != cmap.end(); ++cit) {
        int key = (*cit).first;
        STXXL_CHECK(int_map.find(key) != int_map.end());
    }
    std::cout << "passed" << std::endl;

    // --- another bulk insert
    std::cout << "Compare with internal-memory map after another bulk-insert...";
    map.insert(values3.begin(), values3.end(), mem_to_sort);
    int_map.insert(values3.begin(), values3.end());
    STXXL_CHECK(map.size() == map.size());
    cit = cmap.begin();
    for ( ; cit != cmap.end(); ++cit) {
        int key = (*cit).first;
        STXXL_CHECK(int_map.find(key) != int_map.end());
    }
    std::cout << "passed" << std::endl;
    STXXL_MSG(stxxl::stats_data(*stxxl::stats::get_instance()) - stats_begin);
}

////////////////////////////////////////////////////////////////////////////////
void basic_iterator_test()
{
    typedef std::pair<int, int> value_type;
    const unsigned_type value_size = sizeof(value_type);

    const unsigned_type n_values = 15000;
    const unsigned_type n_tests = 7500;

    // make sure all changes will be buffered
    const unsigned_type buffer_size = 5 * n_values * (value_size + sizeof(int*));

    const unsigned_type mem_to_sort = 32 * 1024 * 1024;

    const unsigned_type subblock_raw_size = 4 * 1024;
    const unsigned_type block_size = 4;

    typedef stxxl::hash_map::hash_map<int, int, hash_int, cmp,
                                      subblock_raw_size, block_size> hash_map;
    typedef hash_map::iterator iterator;
    typedef hash_map::const_iterator const_iterator;

    stxxl::stats_data stats_begin = *stxxl::stats::get_instance();

    hash_map map;
    map.max_buffer_size(buffer_size);
    const hash_map& cmap = map;

    // generate random values
    stxxl::random_number32 rand32;

    std::vector<value_type> values1(n_values);
    std::vector<value_type> values2(n_values);
    std::generate(values1.begin(), values1.end(), rand_pairs(rand32) _STXXL_FORCE_SEQUENTIAL);
    std::generate(values2.begin(), values2.end(), rand_pairs(rand32) _STXXL_FORCE_SEQUENTIAL);

    // --- initial import: create a nice mix of externally (values1) and
    // --- internally (values2) stored values
    std::cout << "Initial import...";

    STXXL_CHECK(map.begin() == map.end());
    map.insert(values1.begin(), values1.end(), mem_to_sort);
    for (std::vector<value_type>::iterator val_it = values2.begin();
         val_it != values2.end(); ++val_it)
        map.insert_oblivious(*val_it);
    STXXL_CHECK(map.begin() != map.end());
    STXXL_CHECK(map.size() == 2 * n_values);
    std::cout << "passed" << std::endl;

    // --- actual testing begins: modfiy random values via iterator
    std::cout << "Lookup and modify...";
    std::random_shuffle(values1.begin(), values1.end());
    std::random_shuffle(values2.begin(), values2.end());
    for (unsigned_type i = 0; i < n_tests; ++i) {
        iterator it1 = map.find(values1[i].first);
        iterator it2 = map.find(values2[i].first);
        STXXL_CHECK(it1 != map.end());
        STXXL_CHECK(it2 != map.end());
        (*it1).second++;
        (*it2).second++;
    }
    // check again
    for (unsigned_type i = 0; i < n_tests; ++i) {
        const_iterator cit1 = cmap.find(values1[i].first);
        const_iterator cit2 = cmap.find(values2[i].first);
        STXXL_CHECK(cit1 != map.end());
        STXXL_CHECK(cit2 != map.end());
        value_type value1 = *cit1;
        value_type value2 = *cit2;
        STXXL_CHECK(value1.second == value1.first + 1);
        STXXL_CHECK(value2.second == value2.first + 1);
    }
    std::cout << "passed" << std::endl;

    // --- scan and modify
    std::cout << "Scan and modify...";
    {
        for (iterator it = map.begin(); it != map.end(); ++it)
            (*it).second = (*it).first + 1;

        for (const_iterator cit = cmap.begin(); cit != cmap.end(); ++cit) {
            STXXL_CHECK((*cit).second == (*cit).first + 1);
        }
    }
    std::cout << "passed" << std::endl;

    // --- interator-value altered by insert_oblivious
    std::cout << "Iterator-value altered by insert_oblivious...";
    std::random_shuffle(values1.begin(), values1.end());
    std::random_shuffle(values2.begin(), values2.end());
    for (unsigned_type i = 0; i < n_tests; i++) {
        int key1 = values1[i].first;
        int key2 = values2[i].first;
        const_iterator cit1 = cmap.find(key1);
        STXXL_CHECK(cit1 != cmap.end());
        const_iterator cit2 = cmap.find(key2);
        STXXL_CHECK(cit2 != cmap.end());

        map.insert_oblivious(value_type(key1, key1 + 3));
        map.insert_oblivious(value_type(key2, key2 + 3));

        STXXL_CHECK((*cit1).second == key1 + 3);
        STXXL_CHECK((*cit2).second == key2 + 3);
    }
    std::cout << "passed" << std::endl;

    // --- iterator-value altered by other iterator
    std::cout << "Iterator-value altered by other iterator...";
    std::random_shuffle(values1.begin(), values1.end());
    std::random_shuffle(values2.begin(), values2.end());
    for (unsigned_type i = 0; i < n_tests; i++) {
        const_iterator cit1 = cmap.find(values1[i].first);
        STXXL_CHECK(cit1 != cmap.end());
        const_iterator cit2 = cmap.find(values2[i].first);
        STXXL_CHECK(cit2 != cmap.end());
        iterator it1 = map.find(values1[i].first);
        STXXL_CHECK(it1 != map.end());
        iterator it2 = map.find(values2[i].first);
        STXXL_CHECK(it2 != map.end());

        (*it1).second = (*it1).first + 5;
        (*it2).second = (*it2).first + 5;
        STXXL_CHECK((*cit1).second == (*cit1).first + 5);
        STXXL_CHECK((*cit2).second == (*cit2).first + 5);
    }
    std::cout << "passed" << std::endl;

    // --- erase by iterator
    std::cout << "Erase by iterator...";
    std::random_shuffle(values1.begin(), values1.end());
    std::random_shuffle(values2.begin(), values2.end());
    for (unsigned_type i = 0; i < n_tests; i++) {
        const_iterator cit1 = cmap.find(values1[i].first);
        STXXL_CHECK(cit1 != cmap.end());
        const_iterator cit2 = cmap.find(values2[i].first);
        STXXL_CHECK(cit2 != cmap.end());
        iterator it1 = map.find(values1[i].first);
        STXXL_CHECK(it1 != map.end());
        iterator it2 = map.find(values2[i].first);
        STXXL_CHECK(it2 != map.end());

        map.erase(it1);
        map.erase(it2);
        STXXL_CHECK(cit1 == cmap.end());
        STXXL_CHECK(cit2 == cmap.end());
    }
    std::cout << "passed" << std::endl;

    // --- erase by value (key)
    std::cout << "Erase by key...";
    for (unsigned_type i = 0; i < n_tests; i++) {
        const_iterator cit1 = cmap.find(values1[i + n_tests].first);
        STXXL_CHECK(cit1 != cmap.end());
        const_iterator cit2 = cmap.find(values2[i + n_tests].first);
        STXXL_CHECK(cit2 != cmap.end());

        map.erase_oblivious(values1[i + n_tests].first);
        map.erase_oblivious(values2[i + n_tests].first);
        STXXL_CHECK(cit1 == cmap.end());
        STXXL_CHECK(cit2 == cmap.end());
    }
    std::cout << "passed" << std::endl;

    STXXL_MSG(stxxl::stats_data(*stxxl::stats::get_instance()) - stats_begin);
}

////////////////////////////////////////////////////////////////////////////////
void more_iterator_test()
{
    typedef std::pair<int, int> value_type;
    const unsigned_type value_size = sizeof(value_type);

    const unsigned_type n_values = 15000;

    // make sure all changes will be buffered
    const unsigned_type buffer_size = 5 * n_values * (value_size + sizeof(int*));
    const unsigned_type mem_to_sort = 32 * 1024 * 1024;

    const unsigned_type subblock_raw_size = 4 * 1024;
    const unsigned_type block_size = 4;

    typedef stxxl::hash_map::hash_map<int, int, hash_int, cmp,
                                      subblock_raw_size, block_size> hash_map;
    typedef hash_map::const_iterator const_iterator;

    stxxl::stats_data stats_begin = *stxxl::stats::get_instance();

    hash_map map;
    map.max_buffer_size(buffer_size);
    const hash_map& cmap = map;

    // generate random values
    stxxl::random_number32 rand32;
    std::vector<value_type> values1(n_values);
    std::vector<value_type> values2(n_values);
    std::generate(values1.begin(), values1.end(), rand_pairs(rand32) _STXXL_FORCE_SEQUENTIAL);
    std::generate(values2.begin(), values2.end(), rand_pairs(rand32) _STXXL_FORCE_SEQUENTIAL);

    // --- initial import
    map.insert(values1.begin(), values1.end(), mem_to_sort);
    for (std::vector<value_type>::iterator val_it = values2.begin();
         val_it != values2.end(); ++val_it)
        map.insert_oblivious(*val_it);

    // --- store some iterators, rebuild and check
    std::cout << "Rebuild test...";
    std::random_shuffle(values1.begin(), values1.end());
    std::random_shuffle(values2.begin(), values2.end());
    {
        const_iterator cit1 = cmap.find(values1[17].first);
        const_iterator cit2 = cmap.find(values2[19].first);
        *cit1;
        *cit2;
        map.rehash();
        STXXL_CHECK(map.size() == 2 * n_values);
        STXXL_CHECK((*cit1).first == values1[17].first);
        STXXL_CHECK((*cit2).first == values2[19].first);
    }
    std::cout << "passed" << std::endl;

    // --- unusual cases while scanning
    std::cout << "Another scan-test...";
    {
        const_iterator cit1 = cmap.find(values1[n_values / 2].first);
        const_iterator cit2 = cit1;
        ++cit1;
        int key1 = (*cit1).first;
        ++cit1;
        int key2 = (*cit1).first;
        map.erase_oblivious(key1);
        map.insert_oblivious(value_type(key2, key2 + 2));

        STXXL_CHECK((*cit1).second == key2 + 2);
        ++cit2;
        STXXL_CHECK(cit1 == cit2);
    }
    std::cout << "passed" << std::endl;

    STXXL_MSG(stxxl::stats_data(*stxxl::stats::get_instance()) - stats_begin);
}

////////////////////////////////////////////////////////////////////////////////

int main()
{
    cmp_with_internal_map();
    basic_iterator_test();
    more_iterator_test();
    return 0;
}
