/***************************************************************************
 *  tests/containers/btree/test_btree_const_scan.cpp
 *
 *  Part of the STXXL. See http://stxxl.sourceforge.net
 *
 *  Copyright (C) 2006 Roman Dementiev <dementiev@ira.uka.de>
 *
 *  Distributed under the Boost Software License, Version 1.0.
 *  (See accompanying file LICENSE_1_0.txt or copy at
 *  http://www.boost.org/LICENSE_1_0.txt)
 **************************************************************************/

#include <iostream>

#include <stxxl/bits/containers/btree/btree.h>
#include <stxxl/timer>

struct comp_type : public std::less<int>
{
    static int max_value()
    {
        return std::numeric_limits<int>::max();
    }
};

#define NODE_BLOCK_SIZE 4096
#define LEAF_BLOCK_SIZE 128 * 1024

struct my_type
{
    stxxl::uint64 data;
    char filler[24];
};

std::ostream& operator << (std::ostream& o, const my_type& obj)
{
    o << " " << obj.data;
    return o;
}

std::ostream& operator << (std::ostream& o, const std::pair<int, double>& obj)
{
    o << obj.first << " " << obj.second;
    return o;
}

typedef stxxl::btree::btree<int, my_type, comp_type,
                            NODE_BLOCK_SIZE, LEAF_BLOCK_SIZE, stxxl::SR> btree_type;

// forced instantiation
template class stxxl::btree::btree<int, my_type, comp_type,
                                   NODE_BLOCK_SIZE, LEAF_BLOCK_SIZE, stxxl::SR>;

#define node_cache_size (25 * 1024 * 1024)
#define leaf_cache_size (6 * LEAF_BLOCK_SIZE)

stxxl::uint64 checksum = 0;

void NC(btree_type& BTree)
{
    stxxl::uint64 sum = 0;
    stxxl::timer Timer1;
    Timer1.start();
    btree_type::iterator it = BTree.begin(), end = BTree.end();
    for ( ; it != end; ++it)
        sum += it->second.data;

    Timer1.stop();
    STXXL_MSG("Scanning with non const iterator: " << Timer1.mseconds() << " msec");
    STXXL_CHECK(sum == checksum);
}

void C(btree_type& BTree)
{
    stxxl::uint64 sum = 0;
    stxxl::timer Timer1;
    Timer1.start();
    btree_type::const_iterator it = BTree.begin(), end = BTree.end();
    for ( ; it != end; ++it)
        sum += it->second.data;

    Timer1.stop();
    STXXL_MSG("Scanning with const iterator: " << Timer1.mseconds() << " msec");
    STXXL_CHECK(sum == checksum);
}

int main(int argc, char* argv[])
{
    if (argc < 2)
    {
        STXXL_MSG("Usage: " << argv[0] << " #ins");
        return -1;
    }

    const unsigned nins = atoi(argv[1]);

    STXXL_MSG("Data set size  : " << nins * sizeof(std::pair<int, my_type>) << " bytes");
    STXXL_MSG("Node cache size: " << node_cache_size << " bytes");
    STXXL_MSG("Leaf cache size: " << leaf_cache_size << " bytes");

    //stxxl::random_number32 rnd;

    std::vector<std::pair<int, my_type> > Data(nins);

    for (unsigned int i = 0; i < nins; ++i)
    {
        Data[i].first = i;
        Data[i].second.data = i;
        checksum += i;
    }
    {
        btree_type BTree1(Data.begin(), Data.end(), comp_type(), node_cache_size, leaf_cache_size, true);
        btree_type BTree2(Data.begin(), Data.end(), comp_type(), node_cache_size, leaf_cache_size, true);

        //STXXL_MSG(*stxxl::stats::get_instance());

        C(BTree1);

        //STXXL_MSG(*stxxl::stats::get_instance());

        NC(BTree2);

        //STXXL_MSG(*stxxl::stats::get_instance());
    }

    {
        btree_type BTree1(Data.begin(), Data.end(), comp_type(), node_cache_size, leaf_cache_size, true);
        btree_type BTree2(Data.begin(), Data.end(), comp_type(), node_cache_size, leaf_cache_size, true);

        STXXL_MSG("Disabling prefetching");
        BTree1.disable_prefetching();
        BTree2.disable_prefetching();

        //STXXL_MSG(*stxxl::stats::get_instance());

        C(BTree1);

        //STXXL_MSG(*stxxl::stats::get_instance());

        NC(BTree2);

        //STXXL_MSG(*stxxl::stats::get_instance());
    }
    STXXL_MSG("All tests passed successfully");

    return 0;
}
