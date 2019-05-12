/***************************************************************************
 *  tools/benchmark_pqueue.cpp
 *
 *  Part of the STXXL. See http://stxxl.sourceforge.net
 *
 *  Copyright (C) 2003 Roman Dementiev <dementiev@mpi-sb.mpg.de>
 *  Copyright (C) 2009 Andreas Beckmann <beckmann@cs.uni-frankfurt.de>
 *  Copyright (C) 2013 Timo Bingmann <tb@panthema.net>
 *
 *  Distributed under the Boost Software License, Version 1.0.
 *  (See accompanying file LICENSE_1_0.txt or copy at
 *  http://www.boost.org/LICENSE_1_0.txt)
 **************************************************************************/

static const char* description =
    "Benchmark the priority queue implementation using a sequence of "
    "operations. The PQ contains pairs of 32- or 64-bit integers, or a "
    "24 byte struct. The operation sequence is either a simple fill/delete "
    "cycle or fill/intermixed inserts/deletes. Because the memory parameters "
    "of the PQ must be set a compile-time, the benchmark provides only "
    "three PQ sizes: for 256 MiB, 1 GiB and 8 GiB of RAM, with the maximum "
    "number of items set accordingly.";

#include <limits>
#include <iomanip>
#include <stxxl/priority_queue>
#include <stxxl/timer>
#include <stxxl/random>
#include <stxxl/cmdline>
#include <stxxl/bits/common/tuple.h>

using stxxl::uint32;
using stxxl::uint64;
using stxxl::internal_size_type;

#define MiB (1024 * 1024)
#define PRINTMOD (16 * MiB)

// *** Integer Pair Types

typedef stxxl::tuple<uint32, uint32> uint32_pair_type;

typedef stxxl::tuple<uint64, uint64> uint64_pair_type;

// *** Larger Structure Type

#define MY_TYPE_SIZE 24

struct my_type : public uint32_pair_type
{
    typedef uint32 key_type;

    char data[MY_TYPE_SIZE - sizeof(uint32_pair_type)];
    my_type() { }
    my_type(const key_type& k1, const key_type& k2)
        : uint32_pair_type(k1, k2)
    {
#if STXXL_WITH_VALGRIND
        memset(data, 0, sizeof(data));
#endif
    }

    static my_type max_value()
    {
        return my_type(std::numeric_limits<key_type>::max(),
                       std::numeric_limits<key_type>::max());
    }
};

template <typename ValueType>
struct my_cmp : public std::binary_function<ValueType, ValueType, bool>
{
    bool operator () (const ValueType& a, const ValueType& b) const
    {
        // PQ is a max priority queue, thus compare greater
        return a.first > b.first;
    }

    ValueType min_value() const
    {
        return ValueType::max_value();
    }
};

static inline void progress(const char* text, uint64 i, uint64 nelements)
{
    if ((i % PRINTMOD) == 0)
        STXXL_MSG(text << " " << i << " ("
                       << std::setprecision(5)
                       << ((double)i * 100.0 / (double)nelements) << " %)");
}

template <typename PQType>
void run_pqueue_insert_delete(uint64 nelements, internal_size_type mem_for_pools)
{
    typedef typename PQType::value_type ValueType;

    // construct priority queue
    PQType pq(mem_for_pools / 2, mem_for_pools / 2);

    pq.dump_sizes();

    STXXL_MSG("Internal memory consumption of the priority queue: " << pq.mem_cons() << " B");
    stxxl::stats_data stats_begin(*stxxl::stats::get_instance());

    {
        stxxl::scoped_print_timer timer("Filling PQ", nelements * sizeof(ValueType));

        for (stxxl::uint64 i = 0; i < nelements; i++)
        {
            progress("Inserting element", i, nelements);

            pq.push(ValueType((int)(nelements - i), 0));
        }
    }

    STXXL_CHECK(pq.size() == nelements);

    STXXL_MSG("Internal memory consumption of the priority queue: " << pq.mem_cons() << " B");

    pq.dump_sizes();

    std::cout << stxxl::stats_data(*stxxl::stats::get_instance()) - stats_begin;
    stats_begin = *stxxl::stats::get_instance();

    {
        stxxl::scoped_print_timer timer("Reading PQ", nelements * sizeof(ValueType));

        for (stxxl::uint64 i = 0; i < nelements; ++i)
        {
            STXXL_CHECK(!pq.empty());
            STXXL_CHECK(pq.top().first == i + 1);

            pq.pop();

            progress("Popped element", i, nelements);
        }
    }

    STXXL_MSG("Internal memory consumption of the priority queue: " << pq.mem_cons() << " B");
    std::cout << stxxl::stats_data(*stxxl::stats::get_instance()) - stats_begin;
}

template <typename PQType>
void run_pqueue_insert_intermixed(uint64 nelements, internal_size_type mem_for_pools)
{
    typedef typename PQType::value_type ValueType;

    // construct priority queue
    PQType pq(mem_for_pools / 2, mem_for_pools / 2);

    pq.dump_sizes();

    STXXL_MSG("Internal memory consumption of the priority queue: " << pq.mem_cons() << " B");
    stxxl::stats_data stats_begin(*stxxl::stats::get_instance());

    {
        stxxl::scoped_print_timer timer("Filling PQ", nelements * sizeof(ValueType));

        for (stxxl::uint64 i = 0; i < nelements; i++)
        {
            progress("Inserting element", i, nelements);

            pq.push(ValueType((int)(nelements - i), 0));
        }
    }

    STXXL_CHECK(pq.size() == nelements);

    STXXL_MSG("Internal memory consumption of the priority queue: " << pq.mem_cons() << " B");

    pq.dump_sizes();

    std::cout << stxxl::stats_data(*stxxl::stats::get_instance()) - stats_begin;
    stats_begin = *stxxl::stats::get_instance();

    stxxl::random_number32 rand;

    {
        stxxl::scoped_print_timer timer("Intermixed Insert/Delete", nelements * sizeof(ValueType));

        for (stxxl::uint64 i = 0; i < nelements; ++i)
        {
            int o = rand() % 3;
            if (o == 0)
            {
                pq.push(ValueType((int)(nelements - i), 0));
            }
            else
            {
                STXXL_CHECK(!pq.empty());
                pq.pop();
            }

            progress("Intermixed element", i, nelements);
        }
    }

    STXXL_MSG("Internal memory consumption of the priority queue: " << pq.mem_cons() << " B");

    pq.dump_sizes();

    std::cout << stxxl::stats_data(*stxxl::stats::get_instance()) - stats_begin;
}

template <typename ValueType,
          internal_size_type mib_for_queue, internal_size_type mib_for_pools,
          uint64 maxvolume>
int do_benchmark_pqueue(uint64 volume, int opseq)
{
    const internal_size_type mem_for_queue = mib_for_queue * MiB;
    const internal_size_type mem_for_pools = mib_for_pools * MiB;

    typedef typename stxxl::PRIORITY_QUEUE_GENERATOR<
            ValueType, my_cmp<ValueType>,
            mem_for_queue,
            maxvolume* MiB / sizeof(ValueType)> gen;

    typedef typename gen::result pq_type;

    STXXL_MSG("Given PQ parameters: " << mib_for_queue << " MiB for queue, "
                                      << mib_for_pools << " MiB for pools, " << maxvolume << " GiB maximum volume.");

    STXXL_MSG("Selected PQ parameters:");
    STXXL_MSG("element size: " << sizeof(ValueType));
    STXXL_MSG("block size: " << pq_type::BlockSize);
    STXXL_MSG("insertion buffer size (N): " << pq_type::N << " items ("
                                            << pq_type::N * sizeof(ValueType) << " B)");
    STXXL_MSG("delete buffer size: " << pq_type::delete_buffer_size);
    STXXL_MSG("maximal arity for internal mergers (AI): " << pq_type::IntKMAX);
    STXXL_MSG("maximal arity for external mergers (AE): " << pq_type::ExtKMAX);
    STXXL_MSG("internal groups: " << pq_type::num_int_groups);
    STXXL_MSG("external groups: " << pq_type::num_ext_groups);
    STXXL_MSG("X : " << gen::X);

    if (volume == 0) volume = 2 * (mem_for_queue + mem_for_pools);

    stxxl::uint64 nelements = volume / sizeof(ValueType);

    STXXL_MSG("Number of elements: " << nelements);

    if (opseq == 0)
    {
        run_pqueue_insert_delete<pq_type>(nelements, mem_for_pools);
        run_pqueue_insert_intermixed<pq_type>(nelements, mem_for_pools);
    }
    else if (opseq == 1)
        run_pqueue_insert_delete<pq_type>(nelements, mem_for_pools);
    else if (opseq == 2)
        run_pqueue_insert_intermixed<pq_type>(nelements, mem_for_pools);
    else
        STXXL_ERRMSG("Invalid operation sequence.");

    return 1;
}

template <typename ValueType>
int do_benchmark_pqueue_config(unsigned pqconfig, uint64 size, unsigned opseq)
{
    if (pqconfig == 0)
    {
        do_benchmark_pqueue_config<ValueType>(1, size, opseq);
        do_benchmark_pqueue_config<ValueType>(2, size, opseq);
        do_benchmark_pqueue_config<ValueType>(3, size, opseq);
        return 1;
    }
    else if (pqconfig == 1)
        return do_benchmark_pqueue<ValueType, 128, 128, 16>(size, opseq);
    else if (pqconfig == 2)
        return do_benchmark_pqueue<ValueType, 512, 512, 64>(size, opseq);
#if __x86_64__ || __LP64__ || (__WORDSIZE == 64)
    else if (pqconfig == 3)
        return do_benchmark_pqueue<ValueType, 4096, 4096, 512>(size, opseq);
#endif
    else
        return 0;
}

int do_benchmark_pqueue_type(unsigned type, unsigned pqconfig, uint64 size, unsigned opseq)
{
    if (type == 0)
    {
        do_benchmark_pqueue_type(1, pqconfig, size, opseq);
        do_benchmark_pqueue_type(2, pqconfig, size, opseq);
        do_benchmark_pqueue_type(3, pqconfig, size, opseq);
        return 1;
    }
    else if (type == 1)
        return do_benchmark_pqueue_config<uint32_pair_type>(pqconfig, size, opseq);
    else if (type == 2)
        return do_benchmark_pqueue_config<uint64_pair_type>(pqconfig, size, opseq);
    else if (type == 3)
        return do_benchmark_pqueue_config<my_type>(pqconfig, size, opseq);
    else
        return 0;
}

int benchmark_pqueue(int argc, char* argv[])
{
    // parse command line
    stxxl::cmdline_parser cp;

    cp.set_description(description);

    uint64 size = 0;
    cp.add_opt_param_bytes("size", size,
                           "Amount of data to insert (e.g. 1GiB)");

    unsigned type = 2;
    cp.add_uint('t', "type", type,
                "Value type of tested priority queue:\n"
                " 1 = pair of uint32,\n"
                " 2 = pair of uint64 (default),\n"
                " 3 = 24 byte struct\n"
                " 0 = all of the above");

    unsigned pqconfig = 2;
    cp.add_uint('p', "pq", pqconfig,
                "Priority queue configuration to test:\n"
                "1 = small (256 MiB RAM, 4 GiB elements)\n"
                "2 = medium (1 GiB RAM, 16 GiB elements) (default)\n"
#if __x86_64__ || __LP64__ || (__WORDSIZE == 64)
                "3 = big (8 GiB RAM, 64 GiB elements)\n"
#endif
                "0 = all of the above");

    unsigned opseq = 1;
    cp.add_uint('o', "opseq", opseq,
                "Operation sequence to perform:\n"
                " 1 = insert all, delete all (default)\n"
                " 2 = insert all, intermixed insert/delete\n"
                " 0 = all of the above");

    if (!cp.process(argc, argv))
        return -1;

    stxxl::config::get_instance();

    if (!do_benchmark_pqueue_type(type, pqconfig, size, opseq))
    {
        STXXL_ERRMSG("Invalid (type,pqconfig) combination.");
    }

    return 0;
}
