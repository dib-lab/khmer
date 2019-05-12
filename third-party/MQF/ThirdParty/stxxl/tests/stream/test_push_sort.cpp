/***************************************************************************
 *  tests/stream/test_push_sort.cpp
 *
 *  Part of the STXXL. See http://stxxl.sourceforge.net
 *
 *  Copyright (C) 2004 Roman Dementiev <dementiev@mpi-sb.mpg.de>
 *  Copyright (C) 2009, 2010 Johannes Singler <singler@kit.edu>
 *
 *  Distributed under the Boost Software License, Version 1.0.
 *  (See accompanying file LICENSE_1_0.txt or copy at
 *  http://www.boost.org/LICENSE_1_0.txt)
 **************************************************************************/

//! \example stream/test_push_sort.cpp
//! This is an example of how to use some basic algorithms from
//! stream package. This example shows how to create
//! \c sorted_runs data structure
//! using \c stream::use_push specialization of \c stream::runs_creator class

#include <limits>
#include <stxxl/stream>

const unsigned long long megabyte = 1024 * 1024;
const unsigned int block_size = 1 * megabyte;

typedef unsigned value_type;

struct Cmp : public std::binary_function<value_type, value_type, bool>
{
    typedef unsigned value_type;
    bool operator () (const value_type& a, const value_type& b) const
    {
        return a < b;
    }
    value_type min_value()
    {
        return std::numeric_limits<value_type>::min();
    }
    value_type max_value()
    {
        return std::numeric_limits<value_type>::max();
    }
};

// special parameter type
typedef stxxl::stream::use_push<value_type> InputType;
typedef stxxl::stream::runs_creator<InputType, Cmp, 4096, stxxl::RC> CreateRunsAlg;
typedef CreateRunsAlg::sorted_runs_type SortedRunsType;

// forced instantiation
template class stxxl::stream::runs_merger<SortedRunsType, Cmp>;

int main()
{
#if STXXL_PARALLEL_MULTIWAY_MERGE
    STXXL_MSG("STXXL_PARALLEL_MULTIWAY_MERGE");
#endif

    unsigned input_size = (50 * megabyte / sizeof(value_type));

    Cmp c;
    CreateRunsAlg SortedRuns(c, 10 * megabyte);
    value_type checksum_before(0);

    stxxl::random_number32 rnd;

    for (unsigned cnt = input_size; cnt > 0; --cnt)
    {
        const value_type element = rnd();
        checksum_before += element;
        SortedRuns.push(element);               // push into the sorter
    }

    SortedRunsType Runs = SortedRuns.result();  // get sorted_runs data structure
    STXXL_CHECK(stxxl::stream::check_sorted_runs(Runs, Cmp()));

    // merge the runs
    stxxl::stream::runs_merger<SortedRunsType, Cmp> merger(Runs, Cmp(), 10 * megabyte);
    stxxl::vector<value_type, 4, stxxl::lru_pager<8>, block_size, STXXL_DEFAULT_ALLOC_STRATEGY> array;
    STXXL_MSG(input_size << " " << Runs->elements);
    STXXL_MSG("checksum before: " << checksum_before);
    value_type checksum_after(0);
    for (unsigned i = 0; i < input_size; ++i)
    {
        checksum_after += *merger;
        array.push_back(*merger);
        ++merger;
    }
    STXXL_MSG("checksum after:  " << checksum_after);
    STXXL_CHECK(stxxl::is_sorted(array.begin(), array.end(), Cmp()));
    STXXL_CHECK(checksum_before == checksum_after);
    STXXL_CHECK(merger.empty());

    return 0;
}
