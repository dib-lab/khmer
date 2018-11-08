/***************************************************************************
 *  tests/containers/test_map_random.cpp
 *
 *  Part of the STXXL. See http://stxxl.sourceforge.net
 *
 *  Copyright (C) 2004, 2005 Thomas Nowak <t.nowak@imail.de>
 *  Copyright (C) 2005, 2006 Roman Dementiev <dementiev@ira.uka.de>
 *
 *  Distributed under the Boost Software License, Version 1.0.
 *  (See accompanying file LICENSE_1_0.txt or copy at
 *  http://www.boost.org/LICENSE_1_0.txt)
 **************************************************************************/

//! \file containers/test_map_random.cpp
//! File for testing functionality of stxxl::map.

//! \example containers/test_map_random.cpp
//! This is an example of use of \c stxxl::map container.

#include <stxxl/map>
#include "map_test_handlers.h"

typedef int key_type;
typedef int data_type;

struct cmp2 : public std::less<int>
{
    static int max_value()
    {
        return std::numeric_limits<int>::max();
    }
};

#define DATA_NODE_BLOCK_SIZE (4096)
#define DATA_LEAF_BLOCK_SIZE (4096)

typedef std::map<key_type, data_type, cmp2> std_map_type;
typedef stxxl::map<key_type, data_type, cmp2,
                   DATA_NODE_BLOCK_SIZE, DATA_LEAF_BLOCK_SIZE> xxl_map_type;

#define PERCENT_CLEAR 1
#define PERCENT_ERASE_BULK 9
#define PERCENT_ERASE_KEY 90
#define PERCENT_ERASE_ITERATOR 100
#define PERCENT_INSERT_PAIR 100
#define PERCENT_INSERT_BULK 100
#define PERCENT_SIZING 100
#define PERCENT_LOWER 100
#define PERCENT_UPPER 200
#define PERCENT_FIND 100
#define PERCENT_ITERATOR 100

//#define MAX_KEY 1000
#define MAX_KEY 10000

//#define MAX_STEP 0x0001000

#define NODE_BLOCK_SIZE xxl_map_type::node_block_type::raw_size
#define LEAF_BLOCK_SIZE xxl_map_type::leaf_block_type::raw_size
#define NODE_MELEMENTS xxl_map_type::node_block_type::size
#define LEAF_MELEMENTS xxl_map_type::leaf_block_type::size

int main(int argc, char* argv[])
{
    typedef std::vector<std::pair<key_type, data_type> > vector_type;

    STXXL_MSG("Node block size: " << NODE_BLOCK_SIZE << " bytes");
    STXXL_MSG("Leaf block size: " << LEAF_BLOCK_SIZE << " bytes");
    STXXL_MSG("Node max elements: " << NODE_MELEMENTS);
    STXXL_MSG("Leaf max elements: " << LEAF_MELEMENTS);

    stxxl::random_number32 rnd;
    //stxxl::ran32State = 1141225706;
    STXXL_MSG("Init random seed: " << stxxl::ran32State);

    int a = (PERCENT_CLEAR +
             PERCENT_SIZING +
             PERCENT_ERASE_BULK +
             PERCENT_ERASE_KEY +
             PERCENT_ERASE_ITERATOR +
             PERCENT_INSERT_PAIR +
             PERCENT_INSERT_BULK +
             PERCENT_LOWER +
             PERCENT_UPPER +
             PERCENT_FIND +
             PERCENT_ITERATOR);

    STXXL_CHECK(a == 1000);

    if (argc < 2)
    {
        STXXL_MSG("Usage: " << argv[0] << " STEP ");
        STXXL_MSG("Note, that STEP must be > 1000");
        return -1;
    }
    stxxl::uint64 MAX_STEP = atoi(argv[1]);
    STXXL_CHECK(MAX_STEP > 1000);
    std_map_type stdmap;
    xxl_map_type xxlmap(NODE_BLOCK_SIZE * 4, LEAF_BLOCK_SIZE * 3);

    for (stxxl::uint64 i = 0; i < MAX_STEP; i++)
    {
        // ***************************************************
        // A random number is created to determine which kind
        // of operation we will be called.
        // ***************************************************

        long step = rnd() % 1000;
        int percent = 0;

        if (i % (MAX_STEP / 100) == 0)
        {
            STXXL_MSG("Step=" << i << " (" << (unsigned)stdmap.size() << ")");
        }

        // *********************************************************
        // The clear function will be called
        // *********************************************************
        if (step < (percent += PERCENT_CLEAR))
        {
            if ((unsigned)rand() % 1000 < stdmap.size())
            {
                stdmap.clear();
                xxlmap.clear();

                STXXL_CHECK(stdmap.empty());
                STXXL_CHECK(xxlmap.empty());
            }
        }
        // *********************************************************
        // The size function will be called
        // *********************************************************
        else if (step < (percent += PERCENT_SIZING))
        {
            std_map_type::size_type size1 = stdmap.size();
            xxl_map_type::size_type size2 = xxlmap.size();

            STXXL_CHECK(size1 == size2);
        }
        // *********************************************************
        // The erase range function will be called
        // *********************************************************
        else if (step < (percent += PERCENT_ERASE_BULK))
        {
            key_type key1 = rand() % MAX_KEY;
            key_type key2 = rand() % MAX_KEY;

            if (key1 > key2)
            {
                std::swap(key1, key2);
            }

            stdmap.erase(stdmap.lower_bound(key1), stdmap.upper_bound(key2));
            xxlmap.erase(xxlmap.lower_bound(key1), xxlmap.upper_bound(key2));

            STXXL_CHECK(stdmap.size() == xxlmap.size());

            STXXL_CHECK(stdmap.lower_bound(key1) == stdmap.end() ||
                        stdmap.lower_bound(key1) == stdmap.upper_bound(key2));
            STXXL_CHECK(xxlmap.lower_bound(key1) == xxlmap.end() ||
                        xxlmap.lower_bound(key1) == xxlmap.upper_bound(key2));
        }
        // *********************************************************
        // The erase a key function will be called
        // *********************************************************
        else if (step < (percent += PERCENT_ERASE_KEY))
        {
            key_type key = rnd() % MAX_KEY;

            stdmap.erase(key);
            xxlmap.erase(key);

            STXXL_CHECK(stxxl::not_there(stdmap, key));
            STXXL_CHECK(stxxl::not_there(xxlmap, key));
        }
        // *********************************************************
        // The erase function will be called
        // *********************************************************
        else if (step < (percent += PERCENT_ERASE_ITERATOR))
        {
            key_type key = rnd() % MAX_KEY;

            std_map_type::iterator stditer = stdmap.find(key);
            xxl_map_type::iterator xxliter = xxlmap.find(key);

            STXXL_CHECK(stxxl::is_end(stdmap, stditer) == is_end(xxlmap, xxliter));

            if (stditer != stdmap.end())
                stdmap.erase(stditer);

            if (xxliter != xxlmap.end())
                xxlmap.erase(xxliter);

            STXXL_CHECK(stxxl::not_there(stdmap, key));
            STXXL_CHECK(stxxl::not_there(xxlmap, key));
        }
        // *********************************************************
        // The insert function will be called
        // *********************************************************
        else if (step < (percent += PERCENT_INSERT_PAIR))
        {
            key_type key = rnd() % MAX_KEY;
            stdmap.insert(std::pair<key_type, data_type>(key, 2 * key));
            xxlmap.insert(std::pair<key_type, data_type>(key, 2 * key));

            STXXL_CHECK(stxxl::there(stdmap, key, 2 * key));
            STXXL_CHECK(stxxl::there(xxlmap, key, 2 * key));
        }
        // *********************************************************
        // The bulk insert function will be called
        // *********************************************************
        else if (step < (percent += PERCENT_INSERT_BULK))
        {
            unsigned lower = rnd() % MAX_KEY;
            unsigned upper = rnd() % MAX_KEY;
            if (lower > upper)
                std::swap(lower, upper);

            vector_type v2(upper - lower);
            for (unsigned j = 0; j < (unsigned)(upper - lower); j++)
            {
                v2[j].first = lower + j;
                v2[j].second = 2 * v2[j].first;
            }

            stdmap.insert(v2.begin(), v2.end());
            xxlmap.insert(v2.begin(), v2.end());

            for (unsigned i = lower; i < upper; i++)
                STXXL_CHECK(stxxl::there(stdmap, i, 2 * i));

            for (unsigned i = lower; i < upper; i++)
                STXXL_CHECK(stxxl::there(xxlmap, i, 2 * i));
        }
        // *********************************************************
        // The lower_bound function will be called
        // *********************************************************
        else if (step < (percent += PERCENT_LOWER))
        {
            key_type key1 = rand() % MAX_KEY;
            key_type key2 = rand() % MAX_KEY;
            if (key1 > key2)
            {
                std::swap(key1, key2);
            }

            while (key1 < key2)
            {
                std_map_type::iterator stditer = stdmap.lower_bound(key1);
                xxl_map_type::iterator xxliter = xxlmap.lower_bound(key1);

                STXXL_CHECK(stxxl::is_end(stdmap, stditer) == is_end(xxlmap, xxliter));
                if (!stxxl::is_end(stdmap, stditer)) {
                    STXXL_CHECK(stxxl::is_same(*(stditer), *(xxliter)));
                }

                key1++;
            }
        }
        // *********************************************************
        // The upper_bound function will be called
        // *********************************************************
        else if (step < (percent += PERCENT_UPPER))
        {
            key_type key1 = rand() % MAX_KEY;
            key_type key2 = rand() % MAX_KEY;
            if (key1 > key2)
            {
                std::swap(key1, key2);
            }

            while (key1 < key2)
            {
                std_map_type::iterator stditer = stdmap.upper_bound(key1);
                xxl_map_type::iterator xxliter = xxlmap.upper_bound(key1);

                STXXL_CHECK(stxxl::is_end(stdmap, stditer) == is_end(xxlmap, xxliter));
                if (!stxxl::is_end(stdmap, stditer)) {
                    STXXL_CHECK(stxxl::is_same(*(stditer), *(xxliter)));
                }

                key1++;
            }
        }
        // *********************************************************
        // The find function will be called
        // *********************************************************
        else if (step < (percent += PERCENT_FIND))
        {
            key_type key1 = rand() % MAX_KEY;
            key_type key2 = rand() % MAX_KEY;
            if (key1 > key2)
            {
                std::swap(key1, key2);
            }

            while (key1 < key2)
            {
                std_map_type::iterator stditer = stdmap.find(key1);
                xxl_map_type::iterator xxliter = xxlmap.find(key1);

                STXXL_CHECK(stxxl::is_end(stdmap, stditer) == stxxl::is_end(xxlmap, xxliter));
                if (!stxxl::is_end(stdmap, stditer)) {
                    STXXL_CHECK(stxxl::is_same(*(stditer), *(xxliter)));
                }

                key1++;
            }
        }
        // *********************************************************
        // The iterate functions will be called
        // *********************************************************
        else if (step < (percent += PERCENT_ITERATOR))
        {
            std_map_type::const_iterator siter1 = stdmap.begin();
            xxl_map_type::const_iterator xiter1 = xxlmap.begin();

            std_map_type::const_iterator siter2 = siter1;
            xxl_map_type::const_iterator xiter2 = xiter1;

            while (siter1 != stdmap.end())
            {
                STXXL_CHECK(xiter1 != xxlmap.end());
                STXXL_CHECK(stxxl::is_same(*(siter1++), *(xiter1++)));
                if (siter1 != stdmap.end()) {
                    STXXL_CHECK(!stxxl::is_same(*siter1, *siter2));
                }
                if (xiter1 != xxlmap.end()) {
                    STXXL_CHECK(!stxxl::is_same(*xiter1, *xiter2));
                }
            }
            STXXL_CHECK(xiter1 == xxlmap.end());
            STXXL_CHECK(siter2 == stdmap.begin());
            STXXL_CHECK(xiter2 == xxlmap.begin());
        }
    }
    return 0;
}
