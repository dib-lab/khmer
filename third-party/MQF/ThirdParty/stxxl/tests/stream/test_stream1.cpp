/***************************************************************************
 *  tests/stream/test_stream1.cpp
 *
 *  Part of the STXXL. See http://stxxl.sourceforge.net
 *
 *  Copyright (C) 2003 Roman Dementiev <dementiev@mpi-sb.mpg.de>
 *
 *  Distributed under the Boost Software License, Version 1.0.
 *  (See accompanying file LICENSE_1_0.txt or copy at
 *  http://www.boost.org/LICENSE_1_0.txt)
 **************************************************************************/

#include <limits>
#include <stxxl/stream>

struct Input
{
    typedef unsigned value_type;
    value_type value;
    value_type rnd_value;
    stxxl::random_number32 rnd;
    value_type crc;
    Input(value_type init) : value(init)
    {
        rnd_value = rnd();
        crc = rnd_value;
    }
    bool empty() const
    {
        return value == 1;
    }
    Input& operator ++ ()
    {
        --value;
        rnd_value = rnd();
        if (!empty())
            crc += rnd_value;

        return *this;
    }
    const value_type& operator * () const
    {
        return rnd_value;
    }
};

struct Cmp : std::binary_function<unsigned, unsigned, bool>
{
    typedef unsigned value_type;
    bool operator () (const value_type& a, const value_type& b) const
    {
        return a < b;
    }
    value_type min_value() const
    {
        return std::numeric_limits<value_type>::min();
    }
    value_type max_value() const
    {
        return std::numeric_limits<value_type>::max();
    }
};

#define MULT (1000)

int main()
{
    typedef stxxl::stream::runs_creator<Input, Cmp, 4096* MULT, stxxl::RC> CreateRunsAlg;
    typedef CreateRunsAlg::sorted_runs_type SortedRunsType;

    stxxl::stats* s = stxxl::stats::get_instance();

    std::cout << *s;

    STXXL_MSG("Size of block type " << sizeof(CreateRunsAlg::block_type));
    unsigned size = MULT * 1024 * 128 / (unsigned)(sizeof(Input::value_type) * 2);
    Input in(size + 1);
    CreateRunsAlg SortedRuns(in, Cmp(), 1024 * 128 * MULT);
    SortedRunsType Runs = SortedRuns.result();
    STXXL_CHECK(stxxl::stream::check_sorted_runs(Runs, Cmp()));
    // merge the runs
    stxxl::stream::runs_merger<SortedRunsType, Cmp> merger(Runs, Cmp(), MULT * 1024 * 128);
    stxxl::vector<Input::value_type> array;
    STXXL_MSG(size << " " << Runs->elements);
    STXXL_MSG("CRC: " << in.crc);
    Input::value_type crc(0);
    for (unsigned i = 0; i < size; ++i)
    {
        crc += *merger;
        array.push_back(*merger);
        ++merger;
    }
    STXXL_MSG("CRC: " << crc);
    STXXL_CHECK(stxxl::is_sorted(array.begin(), array.end(), Cmp()));
    STXXL_CHECK(merger.empty());

    std::cout << *s;

    return 0;
}
