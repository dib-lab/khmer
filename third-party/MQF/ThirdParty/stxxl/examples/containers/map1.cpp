/***************************************************************************
 *  examples/containers/map1.cpp
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
#include <stxxl/map>
#include <iostream>

#define DATA_NODE_BLOCK_SIZE (4096)
#define DATA_LEAF_BLOCK_SIZE (4096)

//! [comparator]
struct CompareGreater
{
    bool operator () (const int& a, const int& b) const
    { return a > b; }

    static int max_value()
    { return std::numeric_limits<int>::min(); }
};
//! [comparator]

int main()
{
    // template parameter <KeyType, DataType, CompareType, RawNodeSize, RawLeafSize, PDAllocStrategy (optional)>
    typedef stxxl::map<int, char, CompareGreater, DATA_NODE_BLOCK_SIZE, DATA_LEAF_BLOCK_SIZE> map_type;

    // Constructor map(node_cache_size_in_bytes, leaf_cache_size_in_bytes)
    map_type my_map((map_type::node_block_type::raw_size)*3, (map_type::leaf_block_type::raw_size)*3);

    my_map.insert(std::pair<int, char>(1, 'a'));
    my_map.insert(std::pair<int, char>(2, 'b'));
    my_map.insert(std::pair<int, char>(3, 'c'));
    my_map.insert(std::pair<int, char>(4, 'd'));

    my_map.erase(3);

    map_type::iterator iter;

    std::cout << "my_map contains:\n";
    for (iter = my_map.begin(); iter != my_map.end(); ++iter)
    {
        std::cout << iter->first << " => " << iter->second << std::endl;
    }

    map_type::iterator iter_low, iter_up;

    iter_low = my_map.lower_bound(1); // iter_low points to (1,a) in this case
    iter_up = my_map.upper_bound(3);  // iter_up points to (2,b) in this case

    std::cout << "lower bound " << iter_low->second << ", upper bound " << iter_up->second << std::endl;

    return 0;
}
//! [example]
