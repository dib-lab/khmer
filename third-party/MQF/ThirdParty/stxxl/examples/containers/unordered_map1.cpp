/***************************************************************************
 *  examples/containers/unordered_map1.cpp
 *
 *  Part of the STXXL. See http://stxxl.sourceforge.net
 *
 *  Copyright (C) 2013 Daniel Feist <daniel.feist@student.kit.edu>
 *  Copyright (C) 2014 Timo Bingmann <tb@panthema.net>
 *
 *  Distributed under the Boost Software License, Version 1.0.
 *  (See accompanying file LICENSE_1_0.txt or copy at
 *  http://www.boost.org/LICENSE_1_0.txt)
 **************************************************************************/

//! [example]
#include <stxxl/unordered_map>
#include <iostream>

//! [hash]
struct HashFunctor
{
    size_t operator () (int key) const
    {
        // a simple integer hash function
        return (size_t)(key * 2654435761u);
    }
};
//! [hash]

//! [comparator]
struct CompareLess
{
    bool operator () (const int& a, const int& b) const
    { return a < b; }

    static int min_value() { return std::numeric_limits<int>::min(); }
    static int max_value() { return std::numeric_limits<int>::max(); }
};
//! [comparator]

int main()
{
//! [construction]
#define SUB_BLOCK_SIZE 8192
#define SUB_BLOCKS_PER_BLOCK 256

    // template parameter <KeyType, MappedType, HashType, CompareType, SubBlockSize, SubBlocksPerBlock>
    typedef stxxl::unordered_map<
            int, char, HashFunctor, CompareLess, SUB_BLOCK_SIZE, SUB_BLOCKS_PER_BLOCK
            > unordered_map_type;

    // constructor: use defaults for all parameters
    unordered_map_type my_map;
//! [construction]

    // insert some items and delete one
    my_map.insert(std::make_pair(1, 'a'));
    my_map.insert(std::make_pair(2, 'b'));
    my_map.insert(std::make_pair(3, 'c'));
    my_map.insert(std::make_pair(4, 'd'));

    my_map.erase(3);

    // iterate over all items in the unordered_map
    unordered_map_type::iterator iter;

    std::cout << "my_map contains:\n";
    for (iter = my_map.begin(); iter != my_map.end(); ++iter)
    {
        std::cout << iter->first << " => " << iter->second << std::endl;
    }

    // direct operator[] access to items
    std::cout << "my_map[2] = " << my_map[2] << std::endl;

    // efficient bulk-insert into hash map by sorting by hash keys
    std::vector<unordered_map_type::value_type> value_array;

    for (int i = 0; i < 128; ++i)
        value_array.push_back(std::make_pair(i, (char)i));

    my_map.insert(value_array.begin(), value_array.end(), 8 * 1024 * 1024);

    // check results of insertion
    std::cout << "my_map[42] = " << my_map[42] << std::endl;
    std::cout << "my_map.size() = " << my_map.size() << std::endl;

    return 0;
}
//! [example]
