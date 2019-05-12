/***************************************************************************
 *  tools/benchmarks/monotonic_pq.cpp
 *
 *  Part of the STXXL. See http://stxxl.sourceforge.net
 *
 *  Copyright (C) 2003 Roman Dementiev <dementiev@mpi-sb.mpg.de>
 *  Copyright (C) 2007, 2009 Johannes Singler <singler@ira.uka.de>
 *  Copyright (C) 2008, 2009 Andreas Beckmann <beckmann@cs.uni-frankfurt.de>
 *
 *  Distributed under the Boost Software License, Version 1.0.
 *  (See accompanying file LICENSE_1_0.txt or copy at
 *  http://www.boost.org/LICENSE_1_0.txt)
 **************************************************************************/

#include <queue>
#include <limits>

#define STXXL_PARALLEL_PQ_MULTIWAY_MERGE_EXTERNAL 1
#define STXXL_PARALLEL_PQ_MULTIWAY_MERGE_INTERNAL 1
#define STXXL_PARALLEL_PQ_MULTIWAY_MERGE_DELETE_BUFFER 1

#define TINY_PQ 0
#define MANUAL_PQ 0

#define SIDE_PQ 1       // compare with second, in-memory PQ (needs a lot of memory)

#include <stxxl/priority_queue>
#include <stxxl/stats>
#include <stxxl/timer>

const stxxl::unsigned_type mega = 1024 * 1024;

#define RECORD_SIZE 16
#define LOAD 0

typedef stxxl::uint64 my_key_type;

#define MAGIC 123

struct my_type
{
    typedef my_key_type key_type;

    key_type key;
#if LOAD
    key_type load;
    char data[RECORD_SIZE - 2 * sizeof(key_type)];
#else
    char data[RECORD_SIZE - sizeof(key_type)];
#endif

    my_type() { }
    my_type(key_type k) : key(k) { }
#if LOAD
    my_type(key_type k, key_type l) : key(k), load(l) { }
#endif

    void operator = (const key_type& k) { key = k; }
#if LOAD
    void operator = (const my_type& mt)
    {
        key = mt.key;
        load = mt.load;
    }
    bool operator == (const my_type& mt) { return (key == mt.key) && (load = mt.load); }
#else
    void operator = (const my_type& mt) { key = mt.key; }
    bool operator == (const my_type& mt) { return key == mt.key; }
#endif
};

std::ostream& operator << (std::ostream& o, const my_type& obj)
{
    o << obj.key;
#if LOAD
    o << "/" << obj.load;
#endif
    return o;
}

//STXXL priority queue is a _maximum_ PQ. "Greater" comparator makes this a "minimum" PQ again.

struct my_cmp /*: public std::binary_function<my_type, my_type, bool>*/ // greater
{
    typedef my_type first_argument_type;
    typedef my_type second_argument_type;
    typedef bool result_type;

    bool operator () (const my_type& a, const my_type& b) const
    {
        return a.key > b.key;
    }

    my_type min_value() const
    {
#if LOAD
        return my_type(std::numeric_limits<my_type::key_type>::max(), MAGIC);
#else
        return my_type(std::numeric_limits<my_type::key_type>::max());
#endif
    }
    my_type max_value() const
    {
#if LOAD
        return my_type(std::numeric_limits<my_type::key_type>::min(), MAGIC);
#else
        return my_type(std::numeric_limits<my_type::key_type>::min());
#endif
    }
};

int main(int argc, char* argv[])
{
    if (argc < 3)
    {
        std::cout << "Usage: " << argv[0] << " [n in MiB]"
            #if defined(STXXL_PARALLEL)
            << " [p threads]"
            #endif
            << std::endl;
        return -1;
    }

    STXXL_MSG("----------------------------------------");

    stxxl::config::get_instance();
    std::string Flags = std::string("")
#if STXXL_CHECK_ORDER_IN_SORTS
                        + " STXXL_CHECK_ORDER_IN_SORTS"
#endif
#ifdef NDEBUG
                        + " NDEBUG"
#endif
#if TINY_PQ
                        + " TINY_PQ"
#endif
#if MANUAL_PQ
                        + " MANUAL_PQ"
#endif
#if SIDE_PQ
                        + " SIDE_PQ"
#endif
#if STXXL_PARALLEL_PQ_MULTIWAY_MERGE_INTERNAL
                        + " STXXL_PARALLEL_PQ_MULTIWAY_MERGE_INTERNAL"
#endif
#if STXXL_PARALLEL_PQ_MULTIWAY_MERGE_EXTERNAL
                        + " STXXL_PARALLEL_PQ_MULTIWAY_MERGE_EXTERNAL"
#endif
#if STXXL_PARALLEL_PQ_MULTIWAY_MERGE_DELETE_BUFFER
                        + " STXXL_PARALLEL_PQ_MULTIWAY_MERGE_DELETE_BUFFER"
#endif
    ;
    STXXL_MSG("Flags:" << Flags);

    unsigned long megabytes = atoi(argv[1]);
#if defined(STXXL_PARALLEL_MODE)
    int num_threads = atoi(argv[2]);
    STXXL_MSG("Threads: " << num_threads);

    omp_set_num_threads(num_threads);
    __gnu_parallel::_Settings parallel_settings(__gnu_parallel::_Settings::get());
    parallel_settings.sort_algorithm = __gnu_parallel::QS_BALANCED;
    parallel_settings.sort_splitting = __gnu_parallel::SAMPLING;
    parallel_settings.sort_minimal_n = 1000;
    parallel_settings.sort_mwms_oversampling = 10;

    parallel_settings.merge_splitting = __gnu_parallel::SAMPLING;
    parallel_settings.merge_minimal_n = 1000;
    parallel_settings.merge_oversampling = 10;

    parallel_settings.multiway_merge_algorithm = __gnu_parallel::LOSER_TREE;
    parallel_settings.multiway_merge_splitting = __gnu_parallel::EXACT;
    parallel_settings.multiway_merge_oversampling = 10;
    parallel_settings.multiway_merge_minimal_n = 1000;
    parallel_settings.multiway_merge_minimal_k = 2;
    __gnu_parallel::_Settings::set(parallel_settings);
#endif

    const stxxl::unsigned_type mem_for_queue = 512 * mega;
    const stxxl::unsigned_type mem_for_pools = 512 * mega;

#if TINY_PQ
    stxxl::STXXL_UNUSED(mem_for_queue);
    const unsigned BufferSize1 = 32;               // equalize procedure call overheads etc.
    const unsigned N = (1 << 9) / sizeof(my_type); // minimal sequence length
    const unsigned IntKMAX = 8;                    // maximal arity for internal mergersq
    const unsigned IntLevels = 2;                  // number of internal levels
    const unsigned BlockSize = (4 * mega);
    const unsigned ExtKMAX = 8;                    // maximal arity for external mergers
    const unsigned ExtLevels = 2;                  // number of external levels
    typedef stxxl::priority_queue<
            stxxl::priority_queue_config<
                my_type,
                my_cmp,
                BufferSize1,
                N,
                IntKMAX,
                IntLevels,
                BlockSize,
                ExtKMAX,
                ExtLevels
                >
            > pq_type;
#elif MANUAL_PQ
    stxxl::STXXL_UNUSED(mem_for_queue);
    const unsigned BufferSize1 = 32;                    // equalize procedure call overheads etc.
    const unsigned N = (1 << 20) / sizeof(my_type);     // minimal sequence length
    const unsigned IntKMAX = 16;                        // maximal arity for internal mergersq
    const unsigned IntLevels = 2;                       // number of internal levels
    const unsigned BlockSize = (4 * mega);
    const unsigned ExtKMAX = 32;                        // maximal arity for external mergers
    const unsigned ExtLevels = 2;                       // number of external levels
    typedef stxxl::priority_queue<
            stxxl::priority_queue_config<
                my_type,
                my_cmp,
                BufferSize1,
                N,
                IntKMAX,
                IntLevels,
                BlockSize,
                ExtKMAX,
                ExtLevels
                >
            > pq_type;
#else
    const stxxl::uint64 volume = stxxl::uint64(200000) * mega;     // in bytes
    typedef stxxl::PRIORITY_QUEUE_GENERATOR<my_type, my_cmp, mem_for_queue, volume / sizeof(my_type) / 1024 + 1> gen;
    typedef gen::result pq_type;
//         BufferSize1 = Config::BufferSize1,
//         N = Config::N,
//         IntKMAX = Config::IntKMAX,
//         IntLevels = Config::IntLevels,
//         ExtLevels = Config::ExtLevels,
//         Levels = Config::IntLevels + Config::ExtLevels,
//         BlockSize = Config::BlockSize,
//         ExtKMAX = Config::ExtKMAX

/*  STXXL_MSG ( "Blocks fitting into internal memory m: "<<gen::m );
  STXXL_MSG ( "X : "<<gen::X );  //maximum number of internal elements //X = B * (settings::k - m) / settings::E,
  STXXL_MSG ( "Expected internal memory consumption: "<< (gen::EConsumption / 1048576) << " MiB");*/
#endif
    STXXL_MSG("Internal arity: " << pq_type::IntKMAX);
    STXXL_MSG("N : " << pq_type::N); //X / (AI * AI)
    STXXL_MSG("External arity: " << pq_type::ExtKMAX);
    STXXL_MSG("Block size B: " << pq_type::BlockSize);
    //EConsumption = X * settings::E + settings::B * AE + ((MaxS_ / X) / AE) * settings::B * 1024

    STXXL_MSG("Data type size: " << sizeof(my_type));
    STXXL_MSG("");

    stxxl::stats_data sd_start(*stxxl::stats::get_instance());
    stxxl::timer Timer;
    Timer.start();

    pq_type p(mem_for_pools / 2, mem_for_pools / 2);
    stxxl::int64 nelements = stxxl::int64(megabytes * mega / sizeof(my_type)), i;

    STXXL_MSG("Internal memory consumption of the priority queue: " << p.mem_cons() << " B");
    STXXL_MSG("Peak number of elements (n): " << nelements);
    STXXL_MSG("Max number of elements to contain: " << (stxxl::uint64(pq_type::N) * pq_type::IntKMAX * pq_type::IntKMAX * pq_type::ExtKMAX * pq_type::ExtKMAX));
    srand(5);
    my_cmp cmp;
    my_key_type r, sum_input = 0, sum_output = 0;
    my_type least(0), last_least(0);

    const my_key_type modulo = 0x10000000;

#if SIDE_PQ
    std::priority_queue<my_type, std::vector<my_type>, my_cmp> side_pq;
#endif

    my_type side_pq_least;

    STXXL_MSG("op-sequence(monotonic pq): ( push, pop, push ) * n");
    for (i = 0; i < nelements; ++i)
    {
        if ((i % mega) == 0)
            STXXL_MSG(
                std::fixed << std::setprecision(2) << std::setw(5)
                           << (100.0 * (double)i / (double)nelements) << "% "
                           << "Inserting element " << i << " top() == " << least.key << " @ "
                           << std::setprecision(3) << Timer.seconds() << " s"
                           << std::setprecision(6) << std::resetiosflags(std::ios_base::floatfield));

        //monotone priority queue
        r = least.key + rand() % modulo;
        sum_input += r;
        p.push(my_type(r));
#if SIDE_PQ
        side_pq.push(my_type(r));
#endif

        least = p.top();
        sum_output += least.key;
        p.pop();
#if SIDE_PQ
        side_pq_least = side_pq.top();
        side_pq.pop();
        if (!(side_pq_least == least))
            STXXL_MSG("Wrong result at  " << i << "  " << side_pq_least.key << " != " << least.key);
#endif

        if (cmp(last_least, least))
        {
            STXXL_MSG("Wrong order at  " << i << "  " << last_least.key << " > " << least.key);
        }
        else
            last_least = least;

        r = least.key + rand() % modulo;
        sum_input += r;
        p.push(my_type(r));
#if SIDE_PQ
        side_pq.push(my_type(r));
#endif
    }
    Timer.stop();
    STXXL_MSG("Time spent for filling: " << Timer.seconds() << " s");

    STXXL_MSG("Internal memory consumption of the priority queue: " << p.mem_cons() << " B");
    stxxl::stats_data sd_middle(*stxxl::stats::get_instance());
    std::cout << sd_middle - sd_start;
    Timer.reset();
    Timer.start();

    STXXL_MSG("op-sequence(monotonic pq): ( pop, push, pop ) * n");
    for (i = 0; i < (nelements); ++i)
    {
        assert(!p.empty());

        least = p.top();
        sum_output += least.key;
        p.pop();
#if SIDE_PQ
        side_pq_least = side_pq.top();
        side_pq.pop();
        if (!(side_pq_least == least))
        {
            STXXL_VERBOSE1("" << side_pq_least << " != " << least);
        }
#endif
        if (cmp(last_least, least))
        {
            STXXL_MSG("Wrong result at " << i << "  " << last_least.key << " > " << least.key);
        }
        else
            last_least = least;

        r = least.key + rand() % modulo;
        sum_input += r;
        p.push(my_type(r));
#if SIDE_PQ
        side_pq.push(my_type(r));
#endif

        least = p.top();
        sum_output += least.key;
        p.pop();
#if SIDE_PQ
        side_pq_least = side_pq.top();
        side_pq.pop();
        if (!(side_pq_least == least))
        {
            STXXL_VERBOSE1("" << side_pq_least << " != " << least);
        }
#endif
        if (cmp(last_least, least))
        {
            STXXL_MSG("Wrong result at " << i << "  " << last_least.key << " > " << least.key);
        }
        else
            last_least = least;

        if ((i % mega) == 0)
            STXXL_MSG(
                std::fixed << std::setprecision(2) << std::setw(5)
                           << (100.0 * (double)i / (double)nelements) << "% "
                           << "Popped element " << i << " == " << least.key << " @ "
                           << std::setprecision(3) << Timer.seconds() << " s"
                           << std::setprecision(6) << std::resetiosflags(std::ios_base::floatfield));
    }
    STXXL_MSG("Last element " << i << " popped");
    Timer.stop();

    if (sum_input != sum_output)
        STXXL_MSG("WRONG sum! " << sum_input << " - " << sum_output << " = " << (sum_output - sum_input) << " / " << (sum_input - sum_output));

    STXXL_MSG("Time spent for removing elements: " << Timer.seconds() << " s");
    STXXL_MSG("Internal memory consumption of the priority queue: " << p.mem_cons() << " B");
    std::cout << stxxl::stats_data(*stxxl::stats::get_instance()) - sd_middle;
    std::cout << *stxxl::stats::get_instance();

    assert(sum_input == sum_output);
}
