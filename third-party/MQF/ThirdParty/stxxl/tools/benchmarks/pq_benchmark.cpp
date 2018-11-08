/***************************************************************************
 *  tools/benchmarks/pq_benchmark.cpp
 *
 *  Part of the STXXL. See http://stxxl.sourceforge.net
 *
 *  Copyright (C) 2006 Roman Dementiev <dementiev@ira.uka.de>
 *  Copyright (C) 2009 Andreas Beckmann <beckmann@cs.uni-frankfurt.de>
 *
 *  Distributed under the Boost Software License, Version 1.0.
 *  (See accompanying file LICENSE_1_0.txt or copy at
 *  http://www.boost.org/LICENSE_1_0.txt)
 **************************************************************************/

//! \example containers/pq_benchmark.cpp
//! This is a benchmark mentioned in the paper
//! R. Dementiev, L. Kettner, P. Sanders "STXXL: standard template library for XXL data sets"
//! Software: Practice and Experience
//! Volume 38, Issue 6, Pages 589-637, May 2008
//! DOI: 10.1002/spe.844

#include <limits>
#include <stxxl/priority_queue>
#include <stxxl/stats>
#include <stxxl/timer>

#define TOTAL_PQ_MEM_SIZE    (768 * 1024 * 1024)

#define PQ_MEM_SIZE                                     (512 * 1024 * 1024)

#define PREFETCH_POOL_SIZE                      ((TOTAL_PQ_MEM_SIZE - PQ_MEM_SIZE) / 2)
#define WRITE_POOL_SIZE                                         (PREFETCH_POOL_SIZE)

#define MAX_ELEMENTS (2000 * 1024 * 1024)

struct my_record
{
    int key;
    int data;
    my_record() : key(0), data(0) { }
    my_record(int k, int d) : key(k), data(d) { }
};

std::ostream& operator << (std::ostream& o, const my_record& obj)
{
    o << obj.key << " " << obj.data;
    return o;
}

bool operator == (const my_record& a, const my_record& b)
{
    return a.key == b.key;
}

bool operator != (const my_record& a, const my_record& b)
{
    return a.key != b.key;
}

bool operator < (const my_record& a, const my_record& b)
{
    return a.key < b.key;
}

bool operator > (const my_record& a, const my_record& b)
{
    return a.key > b.key;
}

struct comp_type : std::binary_function<my_record, my_record, bool>
{
    bool operator () (const my_record& a, const my_record& b) const
    {
        return a > b;
    }
    static my_record min_value()
    {
        return my_record(std::numeric_limits<int>::max(), 0);
    }
};

typedef stxxl::PRIORITY_QUEUE_GENERATOR<my_record, comp_type,
                                        PQ_MEM_SIZE, MAX_ELEMENTS / (1024 / 8)>::result pq_type;

typedef pq_type::block_type block_type;

#define    BLOCK_SIZE block_type::raw_size

#if 1
unsigned ran32State = 0xdeadbeef;
inline int myrand()
{
    return ((int)((ran32State = 1664525 * ran32State + 1013904223) >> 1)) - 1;
}
#else // a longer pseudo random sequence
long long unsigned ran32State = 0xdeadbeef;
inline long long unsigned myrand()
{
    return (ran32State = (ran32State * 0x5DEECE66DULL + 0xBULL) & 0xFFFFFFFFFFFFULL);
}
#endif

void run_stxxl_insert_all_delete_all(stxxl::uint64 ops)
{
    pq_type PQ(PREFETCH_POOL_SIZE, WRITE_POOL_SIZE);

    stxxl::uint64 i;

    my_record cur;

    stxxl::stats_data stats_begin(*stxxl::stats::get_instance());

    stxxl::timer Timer;
    Timer.start();

    for (i = 0; i < ops; ++i)
    {
        cur.key = myrand();
        PQ.push(cur);
    }

    Timer.stop();

    STXXL_MSG("Records in PQ: " << PQ.size());
    if (i != PQ.size())
    {
        STXXL_MSG("Size does not match");
        abort();
    }

    STXXL_MSG("Insertions elapsed time: " << (Timer.mseconds() / 1000.) <<
              " seconds : " << (double(ops) / (Timer.mseconds() / 1000.)) << " key/data pairs per sec");

    std::cout << stxxl::stats_data(*stxxl::stats::get_instance()) - stats_begin;

    ////////////////////////////////////////////////

    stats_begin = *stxxl::stats::get_instance();
    Timer.reset();
    Timer.start();

    for (i = 0; i < ops; ++i)
    {
        PQ.pop();
    }

    Timer.stop();

    STXXL_MSG("Records in PQ: " << PQ.size());
    if (!PQ.empty())
    {
        STXXL_MSG("PQ must be empty");
        abort();
    }

    STXXL_MSG("Deletions elapsed time: " << (Timer.mseconds() / 1000.) <<
              " seconds : " << (double(ops) / (Timer.mseconds() / 1000.)) << " key/data pairs per sec");

    std::cout << stxxl::stats_data(*stxxl::stats::get_instance()) - stats_begin;
}

void run_stxxl_intermixed(stxxl::uint64 ops)
{
    pq_type PQ(PREFETCH_POOL_SIZE, WRITE_POOL_SIZE);

    stxxl::uint64 i;

    my_record cur;

    stxxl::stats_data stats_begin(*stxxl::stats::get_instance());

    stxxl::timer Timer;
    Timer.start();

    for (i = 0; i < ops; ++i)
    {
        cur.key = myrand();
        PQ.push(cur);
    }

    Timer.stop();

    STXXL_MSG("Records in PQ: " << PQ.size());
    if (i != PQ.size())
    {
        STXXL_MSG("Size does not match");
        abort();
    }

    STXXL_MSG("Insertions elapsed time: " << (Timer.mseconds() / 1000.) <<
              " seconds : " << (double(ops) / (Timer.mseconds() / 1000.)) << " key/data pairs per sec");

    std::cout << stxxl::stats_data(*stxxl::stats::get_instance()) - stats_begin;

    ////////////////////////////////////////////////

    stats_begin = *stxxl::stats::get_instance();
    Timer.reset();
    Timer.start();

    for (i = 0; i < ops; ++i)
    {
        int o = myrand() % 3;
        if (o == 0)
        {
            cur.key = myrand();
            PQ.push(cur);
        }
        else
        {
            assert(!PQ.empty());
            PQ.pop();
        }
    }

    Timer.stop();

    STXXL_MSG("Records in PQ: " << PQ.size());

    STXXL_MSG("Deletions/Insertion elapsed time: " << (Timer.mseconds() / 1000.) <<
              " seconds : " << (double(ops) / (Timer.mseconds() / 1000.)) << " key/data pairs per sec");

    std::cout << stxxl::stats_data(*stxxl::stats::get_instance()) - stats_begin;
}

int main(int argc, char* argv[])
{
    STXXL_MSG("stxxl::pq lock size: " << BLOCK_SIZE << " bytes");

#if STXXL_DIRECT_IO_OFF
    STXXL_MSG("STXXL_DIRECT_IO_OFF is defined");
#else
    STXXL_MSG("STXXL_DIRECT_IO_OFF is NOT defined");
#endif

    if (argc < 3)
    {
        STXXL_MSG("Usage: " << argv[0] << " version #ops");
        STXXL_MSG("\t version = 1: insert-all-delete-all stxxl pq");
        STXXL_MSG("\t version = 2: intermixed insert/delete stxxl pq");
        return -1;
    }

    int version = atoi(argv[1]);
    stxxl::uint64 ops = stxxl::atouint64(argv[2]);

    STXXL_MSG("Running version      : " << version);
    STXXL_MSG("Operations to perform: " << ops);

    switch (version)
    {
    case 1:
        run_stxxl_insert_all_delete_all(ops);
        break;
    case 2:
        run_stxxl_intermixed(ops);
        break;
    default:
        STXXL_MSG("Unsupported version " << version);
    }
}
