/***************************************************************************
 *  tests/containers/ppq/benchmark_pqs.cpp
 *
 *  Part of the STXXL. See http://stxxl.sourceforge.net
 *
 *  Copyright (C) 2014-2015 Thomas Keh <thomas.keh@student.kit.edu>
 *  Copyright (C) 2014-2015 Timo Bingmann <tb@panthema.net>
 *
 *  Distributed under the Boost Software License, Version 1.0.
 *  (See accompanying file LICENSE_1_0.txt or copy at
 *  http://www.boost.org/LICENSE_1_0.txt)
 **************************************************************************/

static const char* description =
    "Benchmark different priority queue implementations using a sequence "
    "of operations. The PQ contains of 24 byte values. "
    "The operation sequence is either a simple fill/delete "
    "cycle or fill/intermixed inserts/deletes.";

#include <cstdlib>
#include <ctime>
#include <iomanip>
#include <limits>
#include <queue>

#include <stxxl/bits/common/tuple.h>
#include <stxxl/bits/containers/priority_queue.h>
#include <stxxl/bits/containers/parallel_priority_queue.h>
#include <stxxl/bits/verbose.h>
#include <stxxl/bits/unused.h>
#include <stxxl/cmdline>
#include <stxxl/random>
#include <stxxl/sorter>
#include <stxxl/timer>

#if STXXL_PARALLEL
  #include <omp.h>
#endif

using stxxl::int64;
using stxxl::uint64;

using stxxl::scoped_print_timer;
using stxxl::random_number32_r;

static const int64 KiB = 1024L;
static const int64 MiB = 1024L * KiB;
static const int64 GiB = 1024L * MiB;

// process indicator
static const bool g_progress = false;
static const int64 g_printmod = 16 * MiB;

// constant values for STXXL PQ
const int64 _RAM = 16 * GiB;
const int64 _volume = 1024 * GiB;
const unsigned _value_size = 8;
const uint64 _num_elements = _volume / _value_size;
const int64 _mem_for_queue = _RAM / 2;
const int64 mem_for_prefetch_pool = _RAM / 4;
const int64 mem_for_write_pool = _RAM / 4;

std::string g_testset;

uint64 RAM = _RAM;
uint64 volume = _volume;
uint64 num_elements = _num_elements;

uint64 single_heap_ram = 1 * MiB;
uint64 extract_buffer_ram = 0;
unsigned num_write_buffers = 14;
double num_prefetchers = 2.0;
unsigned num_read_blocks = 1;
int64 block_size = STXXL_DEFAULT_BLOCK_SIZE(value_type);

#if STXXL_PARALLEL
const unsigned g_max_threads = omp_get_max_threads();
#else
const unsigned g_max_threads = 1;
#endif

static inline int get_thread_num()
{
#if STXXL_PARALLEL
    return omp_get_thread_num();
#else
    return 0;
#endif
}

int64 ram_write_buffers;
uint64 bulk_size = 123456;

int64 num_arrays = 0;
int64 max_merge_buffer_ram_stat = 0;

uint64 value_universe_size = 100000000;

unsigned int g_seed = 12345;

random_number32_r g_rand(g_seed);

/*
 * The value type and corresponding benchmark functions.
 */

template <typename KeyType, int ValueSize>
struct my_type
{
    enum { value_size = ValueSize };
    typedef KeyType key_type;

    union {
        key_type key;
        char padding[value_size];
    };

    my_type()
    { }

    my_type(const key_type& _key)
        : key(_key)
    {
#ifdef STXXL_VALGRIND_AVOID_UNINITIALIZED_WRITE_ERRORS
        memset(padding + sizeof(key_type), 0, sizeof(padding) - sizeof(key_type));
#endif
    }

    friend std::ostream& operator << (std::ostream& o, const my_type& obj)
    {
        o << obj.key;
        return o;
    }
};

template <typename ValueType>
struct my_type_cmp_smaller
    : public std::binary_function<ValueType, ValueType, bool>
{
    typedef ValueType value_type;
    typedef typename value_type::key_type key_type;
    bool operator () (const value_type& a, const value_type& b) const
    {
        return a.key < b.key;
    }
    value_type min_value() const
    {
        return value_type(std::numeric_limits<key_type>::min());
    }
    value_type max_value() const
    {
        return value_type(std::numeric_limits<key_type>::max());
    }
};

template <typename ValueType>
struct my_type_cmp_greater
    : public std::binary_function<ValueType, ValueType, bool>
{
    typedef ValueType value_type;
    typedef typename value_type::key_type key_type;
    bool operator () (const value_type& a, const value_type& b) const
    {
        return a.key > b.key;
    }
    value_type min_value() const
    {
        return value_type(std::numeric_limits<key_type>::max());
    }
};

template <typename ValueType>
struct less_min_max : public std::binary_function<ValueType, ValueType, bool>
{
    bool operator () (const ValueType& a, const ValueType& b) const
    {
        return a < b;
    }
    ValueType min_value() const
    {
        return std::numeric_limits<ValueType>::min();
    }
    ValueType max_value() const
    {
        return std::numeric_limits<ValueType>::max();
    }
};

/*
 * my_type specializations used in benchmarks
 */

typedef my_type<uint64, sizeof(uint64)> my8_type;
typedef my_type<uint64, 24> my24_type;

/*
 * Progress messages
 */

static inline void progress(const char* text, int64 i, int64 nelements)
{
    if (!g_progress) return;

    if ((i % g_printmod) == 0) {
        STXXL_MSG(text << " " << i << " (" << std::setprecision(5) <<
                  (static_cast<double>(i) * 100. / static_cast<double>(nelements)) << " %)");
    }
}

/*!
 * Abstract Base Container class with derived classes for sorters and PQs.
 */

class PQBase
{
protected:
    std::string m_name, m_shortname;

    bool m_have_warned;

public:
    PQBase(const std::string& name, const std::string& shortname)
        : m_name(name), m_shortname(shortname),
          m_have_warned(false)
    { }

    //! Return a longer text name
    const std::string & name() const
    {
        return m_name;
    }

    //! Return a short id
    const std::string & shortname() const
    {
        return m_shortname;
    }

    //! Operation called after filling the container.
    void fill_end()
    { }

    //! Check whether the container supports thread-safe bulk pushes
    bool supports_threaded_bulks() const
    {
        return false;
    }
};

// *** Base Wrapper containing many Stubs

class PQBaseDefaults : public PQBase
{
public:
    PQBaseDefaults(const std::string& name, const std::string& shortname)
        : PQBase(name, shortname)
    { }

    void bulk_push_begin(size_t /* this_bulk_size */)
    {
        if (!m_have_warned) {
            std::cout << name() << ": bulk_push_begin not supported\n";
            m_have_warned = true;
        }
    }

    void bulk_push_end()
    {
        if (!m_have_warned) {
            std::cout << name() << ": bulk_push_end not supported\n";
        }
    }

    template <typename value_type>
    void limit_begin(const value_type& /* v */, size_t /* bulk_size */)
    {
        if (!m_have_warned) {
            std::cout << name() << ": limit_begin not supported\n";
            m_have_warned = true;
        }
    }

    void limit_end()
    {
        if (!m_have_warned) {
            std::cout << name() << ": limit_end not supported\n";
        }
    }
};

// *** Wrapper for STXXL PQ

template <typename ValueType>
class CStxxlPQ
    : public PQBaseDefaults,
      public stxxl::PRIORITY_QUEUE_GENERATOR<
          ValueType,
          my_type_cmp_greater<ValueType>,
          _mem_for_queue,
          _num_elements / 1024>::result
{
public:
    typedef ValueType value_type;

    typedef typename stxxl::PRIORITY_QUEUE_GENERATOR<
            value_type,
            my_type_cmp_greater<value_type>,
            _mem_for_queue,
            _num_elements / 1024>::result pq_type;

    CStxxlPQ() :
#if STXXL_PARALLEL
        PQBaseDefaults("STXXL PQ (parallelized)", "spq"),
#else
        PQBaseDefaults("STXXL PQ (sequential)", "spqs"),
#endif
        pq_type(mem_for_prefetch_pool, mem_for_write_pool)
    { }

    //! Simple single top() and pop() operation.
    value_type top_pop()
    {
        value_type item = pq_type::top();
        pq_type::pop();
        return item;
    }

    //! Emulate bulk_push() using plain push().
    void bulk_push(const value_type& v, int /* thread_num */)
    {
        return pq_type::push(v);
    }

    //! Emulate bulk_pop() using iterated pop().
    void bulk_pop(std::vector<value_type>& v, size_t max_size)
    {
        v.clear();
        v.reserve(max_size);

        while (!pq_type::empty() && max_size > 0)
        {
            v.push_back(top_pop());
            --max_size;
        }
    }

    //! Emulate limit_push() using plain push().
    void limit_push(const value_type& v, int /* thread_num */)
    {
        pq_type::push(v);
    }

    //! Emulate limit_top() using top().
    const value_type & limit_top()
    {
        return pq_type::top();
    }

    value_type limit_top_pop()
    {
        return top_pop();
    }
};

// *** Wrapper for STXXL Parallel PQ

template <typename ValueType>
class CStxxlParallePQ
    : public PQBase,
      public stxxl::parallel_priority_queue<ValueType, my_type_cmp_greater<ValueType> >
{
public:
    typedef ValueType value_type;

    typedef stxxl::parallel_priority_queue<ValueType, my_type_cmp_greater<ValueType> > pq_type;

    CStxxlParallePQ()
        : PQBase("Parallel PQ", "ppq"),
          pq_type(my_type_cmp_greater<ValueType>(),
                  RAM, num_prefetchers, num_write_buffers,
                  g_max_threads, single_heap_ram, extract_buffer_ram)
    { }

    //! Check whether the container supports thread-safe bulk pushes
    bool supports_threaded_bulks() const
    {
        return true;
    }

    value_type top_pop()
    {
        value_type top = pq_type::top();
        pq_type::pop();
        return top;
    }

    value_type limit_top_pop()
    {
        value_type top = pq_type::limit_top();
        pq_type::limit_pop();
        return top;
    }
};

// *** Wrapper for STL main memory PQ

template <typename ValueType>
class CStdPriorityQueue
    : public PQBaseDefaults,
      public std::priority_queue<
          ValueType, std::vector<ValueType>,
          my_type_cmp_greater<ValueType>
          >
{
public:
    typedef ValueType value_type;

    typedef std::priority_queue<
            ValueType, std::vector<ValueType>,
            my_type_cmp_greater<ValueType>
            > pq_type;

    CStdPriorityQueue()
        : PQBaseDefaults("STL PQ", "stl"), pq_type()
    { }

    value_type top_pop()
    {
        value_type item = pq_type::top();
        pq_type::pop();
        return item;
    }

    void bulk_push(const value_type& v, int /* thread_num */)
    {
        return pq_type::push(v);
    }

    void bulk_pop(std::vector<value_type>& v, size_t max_size)
    {
        v.clear();
        v.reserve(max_size);

        while (!this->empty() && max_size > 0)
        {
            v.push_back(top_pop());
            --max_size;
        }
    }

    void limit_push(const value_type& v, int /* thread_num */)
    {
        pq_type::push(v);
    }

    const value_type & limit_top()
    {
        return pq_type::top();
    }

    value_type limit_top_pop()
    {
        return top_pop();
    }
};

// *** Wrapper for STL Sorter

template <typename ValueType>
class CStxxlSorter
    : public PQBaseDefaults,
      public stxxl::sorter<ValueType, my_type_cmp_smaller<ValueType> >
{
public:
    typedef ValueType value_type;

    // note that we use SMALLER here:
    typedef stxxl::sorter<ValueType, my_type_cmp_smaller<ValueType> > pq_type;

    CStxxlSorter()
        : PQBaseDefaults("STXXL Sorter", "sorter"),
          pq_type(my_type_cmp_smaller<ValueType>(), RAM)
    { }

    void fill_end()
    {
        pq_type::sort(); // switch state
    }

    const value_type & top()
    {
        return pq_type::operator * ();
    }

    void pop()
    {
        pq_type::operator ++ ();
    }

    value_type top_pop()
    {
        value_type item = top();
        pop();
        return item;
    }

    void bulk_push(const value_type& v, int /* thread_num */)
    {
        return pq_type::push(v);
    }

    void bulk_pop(std::vector<value_type>& v, size_t max_size)
    {
        v.clear();
        v.reserve(max_size);

        while (!this->empty() && max_size > 0)
        {
            v.push_back(top_pop());
            --max_size;
        }
    }

    void limit_push(const value_type& v, int /* thread_num */)
    {
        pq_type::push(v);
    }

    const value_type & limit_top()
    {
        return top();
    }

    value_type limit_top_pop()
    {
        return top_pop();
    }
};

/*
 * Timer and Statistics
 */

class scoped_stats : public stxxl::scoped_print_timer
{
protected:
    //! short name of tested PQ
    std::string m_shortname;

    //! operation on PQ being timed
    std::string m_op;

    //! stxxl stats at start of operation
    stxxl::stats_data m_stxxl_stats;

    //! size of values in container
    size_t m_value_size;

public:
    template <typename ContainerType>
    scoped_stats(ContainerType& c, const std::string& op,
                 const std::string& message,
                 unsigned int rwcycles = 1)
        : stxxl::scoped_print_timer(message,
                                    rwcycles * num_elements *
                                    ContainerType::value_type::value_size),
          m_shortname(c.shortname()),
          m_op(op),
          m_stxxl_stats(*stxxl::stats::get_instance()),
          m_value_size(ContainerType::value_type::value_size)
    { }

    ~scoped_stats()
    {
        stxxl::stats_data s =
            stxxl::stats_data(*stxxl::stats::get_instance()) - m_stxxl_stats;

        const char* hostname = getenv("HOSTNAME");

        std::cout << "RESULT"
                  << " host=" << hostname
                  << " pqname=" << m_shortname
                  << " op=" << m_op
                  << " value_size=" << m_value_size
                  << " num_elements=" << num_elements
                  << " volume=" << num_elements * m_value_size
                  << " time=" << timer().seconds()
                  << " read_count=" << s.get_reads()
                  << " read_vol=" << s.get_read_volume()
                  << " read_time=" << s.get_read_time()
                  << " pread_time=" << s.get_pread_time()
                  << " write_count=" << s.get_writes()
                  << " write_vol=" << s.get_written_volume()
                  << " write_time=" << s.get_write_time()
                  << " pwrite_time=" << s.get_pwrite_time()
                  << " wait_io_time=" << s.get_io_wait_time()
                  << " wait_read_time=" << s.get_wait_read_time()
                  << " wait_write_time=" << s.get_wait_write_time()
                  << " ram=" << RAM
                  << " block_size=" << block_size
                  << " bulk_size=" << bulk_size
                  << " num_threads=" << g_max_threads
                  << " single_heap_ram=" << single_heap_ram
                  << " universe_size=" << value_universe_size
                  << std::endl;
    }
};

/*
 * Random Generator spaced onto different cache lines
 */

class random_number32_r_bigger : public stxxl::random_number32_r
{
protected:
    // padding to put random number state onto different cache lines
    char padding[32];
};

/*
 * Benchmark Functions
 */

template <typename ContainerType>
void do_push_asc(ContainerType& c)
{
    typedef typename ContainerType::value_type value_type;

    scoped_stats stats(c, g_testset + "|push-asc",
                       "Filling " + c.name() + " sequentially");

    for (uint64 i = 0; i < num_elements; i++)
    {
        progress("Inserting element", i, num_elements);
        c.push(value_type(num_elements - i - 1));
    }

    c.fill_end();
}

template <typename ContainerType>
void do_push_rand(ContainerType& c, unsigned int seed)
{
    typedef typename ContainerType::value_type value_type;

    scoped_stats stats(c, g_testset + "|push-rand",
                       "Filling " + c.name() + " randomly");

    random_number32_r datarng(seed);

    for (uint64 i = 0; i < num_elements; i++)
    {
        int64 k = datarng() % value_universe_size;
        progress("Inserting element", i, num_elements);
        c.push(value_type(k));
    }

    c.fill_end();
}

template <typename ContainerType>
void do_pop(ContainerType& c)
{
    scoped_stats stats(c, g_testset + "|pop",
                       "Reading " + c.name());

    STXXL_CHECK_EQUAL(c.size(), num_elements);

    for (uint64 i = 0; i < num_elements; ++i)
    {
        STXXL_CHECK(!c.empty());
        STXXL_CHECK_EQUAL(c.size(), num_elements - i);
        c.pop();
        progress("Popped element", i, num_elements);
    }
}

template <typename ContainerType>
void do_pop_asc_check(ContainerType& c)
{
    typedef typename ContainerType::value_type value_type;

    scoped_stats stats(c, g_testset + "|pop-check",
                       "Reading " + c.name() + " and checking order");

    STXXL_CHECK_EQUAL(c.size(), num_elements);

    for (uint64 i = 0; i < num_elements; ++i)
    {
        STXXL_CHECK(!c.empty());
        STXXL_CHECK_EQUAL(c.size(), num_elements - i);

        value_type top = c.top_pop();

        STXXL_CHECK_EQUAL(top.key, i);
        progress("Popped element", i, num_elements);
    }
}

template <typename ContainerType>
void do_pop_rand_check(ContainerType& c, unsigned int seed)
{
    typedef typename ContainerType::value_type value_type;

    stxxl::sorter<uint64, less_min_max<uint64> >
    sorted_vals(less_min_max<uint64>(), RAM / 2);

    {
        scoped_print_timer timer("Filling sorter for comparison");

        random_number32_r datarng(seed);

        for (uint64 i = 0; i < num_elements; ++i)
        {
            uint64 k = datarng() % value_universe_size;
            sorted_vals.push(k);
        }
    }
    sorted_vals.sort();

    {
        scoped_stats stats(c, g_testset + "|pop-check",
                           "Reading " + c.name() + " and check order");

        for (uint64 i = 0; i < num_elements; ++i)
        {
            STXXL_CHECK(!c.empty());
            STXXL_CHECK_EQUAL(c.size(), num_elements - i);
            STXXL_CHECK_EQUAL(c.size(), sorted_vals.size());

            value_type top = c.top_pop();

            STXXL_CHECK_EQUAL(top.key, *sorted_vals);
            progress("Popped element", i, num_elements);
            ++sorted_vals;
        }
    }
}

template <typename ContainerType>
void do_bulk_push_asc(ContainerType& c, bool parallel)
{
    typedef typename ContainerType::value_type value_type;

    std::string parallel_str = parallel ? " in parallel" : " sequentially";

    scoped_stats stats(c, g_testset + "|bulk-push-asc",
                       "Filling " + c.name() + " with bulks" + parallel_str);

    if (parallel && !c.supports_threaded_bulks()) {
        STXXL_ERRMSG("Container " + c.name() +
                     " does not support parallel bulks operations");
        abort();
    }

    for (uint64 i = 0; i < num_elements / bulk_size; ++i)
    {
        c.bulk_push_begin(bulk_size);

#if STXXL_PARALLEL
        #pragma omp parallel if(parallel) num_threads(g_max_threads)
#endif
        {
            const unsigned thread_id = parallel ? get_thread_num() : 0;

#if STXXL_PARALLEL
            #pragma omp for schedule(static)
#endif
            for (int64 j = 0; j < (int64)bulk_size; ++j)
            {
                c.bulk_push(value_type(num_elements - 1 - i * bulk_size - j), thread_id);
            }
        }

        c.bulk_push_end();

        progress("Inserting element", i * bulk_size, num_elements);
    }

    c.bulk_push_begin(num_elements % bulk_size);

#if STXXL_PARALLEL
#pragma omp parallel if(parallel)
#endif
    {
        const unsigned thread_id = parallel ? get_thread_num() : 0;

#if STXXL_PARALLEL
        #pragma omp for schedule(static)
#endif
        for (int64 j = (num_elements / bulk_size) * bulk_size;
             j < (int64)num_elements; j++)
        {
            c.bulk_push(value_type(num_elements - 1 - j), thread_id);
        }
    }

    c.bulk_push_end();

    progress("Inserting element", num_elements, num_elements);

    STXXL_CHECK_EQUAL(c.size(), num_elements);

    c.fill_end();
}

template <typename ContainerType>
void do_bulk_push_rand(ContainerType& c,
                       unsigned int _seed, bool parallel)
{
    typedef typename ContainerType::value_type value_type;

    std::string parallel_str = parallel ? " in parallel" : " sequentially";

    scoped_stats stats(c, g_testset + "|bulk-push-rand",
                       "Filling " + c.name() + " with bulks" + parallel_str);

    if (parallel && !c.supports_threaded_bulks()) {
        STXXL_ERRMSG("Container " + c.name() +
                     " does not support parallel bulks operations");
        abort();
    }

    // independent random generator (also: in different cache lines)
    std::vector<random_number32_r_bigger> datarng(g_max_threads);
    for (unsigned int t = 0; t < g_max_threads; ++t)
        datarng[t].set_seed(t * _seed);

    for (uint64 i = 0; i < num_elements / bulk_size; ++i)
    {
        c.bulk_push_begin(bulk_size);

#if STXXL_PARALLEL
#pragma omp parallel if(parallel) num_threads(g_max_threads)
#endif
        {
            const unsigned thread_id = parallel ? get_thread_num() : 0;

#if STXXL_PARALLEL
            #pragma omp for schedule(static)
#endif
            for (int64 j = 0; j < (int64)bulk_size; ++j)
            {
                int64 k = datarng[thread_id]() % value_universe_size;
                c.bulk_push(value_type(k), thread_id);
            }
        }

        c.bulk_push_end();

        progress("Inserting element", i * bulk_size, num_elements);
    }

    STXXL_CHECK_EQUAL(c.size(), num_elements - (num_elements % bulk_size));

    int64 bulk_remain = num_elements % bulk_size;

    c.bulk_push_begin(bulk_remain);

#if STXXL_PARALLEL
#pragma omp parallel if(parallel)
#endif
    {
        const unsigned thread_id = parallel ? get_thread_num() : 0;

#if STXXL_PARALLEL
        #pragma omp for schedule(static)
#endif
        for (int64 j = 0; j < bulk_remain; ++j)
        {
            int64 k = datarng[thread_id]() % value_universe_size;
            c.bulk_push(value_type(k), thread_id);
        }
    }

    c.bulk_push_end();

    progress("Inserting element", num_elements, num_elements);

    STXXL_CHECK_EQUAL(c.size(), num_elements);

    c.fill_end();
}

template <typename ContainerType>
void do_bulk_pop(ContainerType& c)
{
    typedef typename ContainerType::value_type value_type;

    scoped_stats stats(c, g_testset + "|bulk-pop",
                       "Reading " + c.name() + " in bulks");

    STXXL_CHECK_EQUAL(c.size(), num_elements);

    std::vector<value_type> work;

    for (uint64 i = 0; i < num_elements; i += work.size())
    {
        STXXL_CHECK(!c.empty());
        STXXL_CHECK_EQUAL(c.size(), num_elements - i);

        c.bulk_pop(work, bulk_size);

        STXXL_CHECK(work.size() <= std::min(bulk_size, num_elements - i));
        progress("Popped element", i, num_elements);
    }

    STXXL_CHECK(c.empty());
}

template <typename ContainerType>
void do_bulk_pop_check_asc(ContainerType& c)
{
    typedef typename ContainerType::value_type value_type;

    scoped_stats stats(c, g_testset + "|bulk-pop-check",
                       "Reading " + c.name() + " in bulks");

    STXXL_CHECK_EQUAL(c.size(), num_elements);

    std::vector<value_type> work;

    for (uint64 i = 0; i < num_elements; i += work.size())
    {
        STXXL_CHECK(!c.empty());
        STXXL_CHECK_EQUAL(c.size(), num_elements - i);

        c.bulk_pop(work, bulk_size);

        STXXL_CHECK(work.size() <= std::min(bulk_size, num_elements - i));

        for (uint64 j = 0; j < work.size(); ++j)
            STXXL_CHECK_EQUAL(work[j].key, i + j);

        progress("Popped element", i, num_elements);
    }

    STXXL_CHECK(c.empty());
}

template <typename ContainerType>
void do_bulk_pop_check_rand(ContainerType& c,
                            unsigned int _seed, bool parallel)
{
    typedef typename ContainerType::value_type value_type;

    std::string parallel_str = parallel ? " in parallel" : " sequentially";

    stxxl::sorter<uint64, less_min_max<uint64> >
    sorted_vals(less_min_max<uint64>(), RAM / 2);

    {
        scoped_print_timer timer("Filling sorter for comparison",
                                 num_elements * value_type::value_size);

        // independent random generator (also: in different cache lines)
        std::vector<random_number32_r_bigger> datarng(g_max_threads);
        for (unsigned int t = 0; t < g_max_threads; ++t)
            datarng[t].set_seed(t * _seed);

        for (uint64 i = 0; i < num_elements / bulk_size; ++i)
        {
#if STXXL_PARALLEL
#pragma omp parallel if(parallel) num_threads(g_max_threads)
#endif
            {
                const unsigned thread_id = parallel ? get_thread_num() : 0;

#if STXXL_PARALLEL
                #pragma omp for schedule(static)
#endif
                for (int64 j = 0; j < (int64)bulk_size; ++j)
                {
                    uint64 k = datarng[thread_id]() % value_universe_size;
#if STXXL_PARALLEL
                    #pragma omp critical
#endif
                    sorted_vals.push(k);
                }
            }

            progress("Inserting element", i * bulk_size, num_elements);
        }

        STXXL_CHECK_EQUAL(sorted_vals.size(),
                          num_elements - (num_elements % bulk_size));

        int64 bulk_remain = num_elements % bulk_size;

#if STXXL_PARALLEL
#pragma omp parallel if(parallel)
#endif
        {
            const unsigned thread_id = parallel ? get_thread_num() : 0;

#if STXXL_PARALLEL
            #pragma omp for schedule(static)
#endif
            for (int64 j = 0; j < bulk_remain; ++j)
            {
                uint64 k = datarng[thread_id]() % value_universe_size;
#if STXXL_PARALLEL
                #pragma omp critical
#endif
                sorted_vals.push(k);
            }
        }

        progress("Inserting element", num_elements, num_elements);
    }

    sorted_vals.sort();
    STXXL_CHECK_EQUAL(c.size(), sorted_vals.size());

    {
        scoped_stats stats(c, g_testset + "|bulk-pop-check",
                           "Reading " + c.name() + " in bulks");

        STXXL_CHECK_EQUAL(c.size(), num_elements);

        std::vector<value_type> work;

        for (uint64 i = 0; i < num_elements; i += bulk_size)
        {
            STXXL_CHECK(!c.empty());
            STXXL_CHECK_EQUAL(c.size(), num_elements - i);

            c.bulk_pop(work, bulk_size);

            STXXL_CHECK_EQUAL(work.size(), std::min(bulk_size, num_elements - i));

            for (uint64 j = 0; j < work.size(); ++j)
            {
                STXXL_CHECK_EQUAL(work[j].key, *sorted_vals);
                ++sorted_vals;
            }

            progress("Popped element", i, num_elements);
        }
    }

    STXXL_CHECK(c.empty());
    STXXL_CHECK(sorted_vals.empty());
}

#if 0

template <typename ContainerType>
void do_rand_intermixed(ContainerType& c, unsigned int seed, bool filled)
{
    int64 num_inserts = filled ? num_elements : 0;
    int64 num_deletes = 0;

    scoped_stats stats(c, "intermix-rand",
                       c.name() + ": Intermixed rand insert and delete",
                       2);

    random_number32_r datarng(seed);

    for (int64 i = 0; i < 2 * num_elements; ++i)
    {
        unsigned r = datarng() % 2;

        if (num_deletes < num_inserts &&
            (r > 0 || num_inserts >= (filled ? 2 * num_elements : num_elements)))
        {
            STXXL_CHECK(!c.empty());

            c.pop();

            progress("Deleting / inserting element", i, 2 * num_elements);
            ++num_deletes;
            STXXL_CHECK_EQUAL(c.size(), num_inserts - num_deletes);
        }
        else
        {
            progress("Inserting/deleting element", i, 2 * num_elements);

            int64 k = datarng() % value_universe_size;
            c.push(value_type(k));

            num_inserts++;
            STXXL_CHECK_EQUAL(c.size(), num_inserts - num_deletes);
        }
    }
}

template <typename ContainerType>
void do_bulk_rand_intermixed(ContainerType& c,
                             unsigned int _seed, bool filled,
                             bool parallel)
{
    std::string parallel_str = parallel ? " in parallel" : "";

    int64 num_inserts = filled ? num_elements : 0;
    int64 num_deletes = 0;

    scoped_stats stats(c, "intermix-bulk-rand",
                       c.name() + ": Intermixed parallel rand bulk insert"
                       + parallel_str + " and delete",
                       (filled ? 3 : 2));

    for (int64 i = 0; i < (filled ? 3 : 2) * num_elements; ++i)
    {
        int64 r = g_rand() % ((filled ? 2 : 1) * bulk_size);

        // probability for an extract is bulk_size times larger than for a
        // bulk_insert.  when already filled with num_elements elements:
        // probability for an extract is 2*bulk_size times larger than for
        // a bulk_insert.
        if (num_deletes < num_inserts &&
            num_deletes < (filled ? 2 : 1) * num_elements &&
            (r > 0 || num_inserts >= (filled ? 2 : 1) * num_elements))
        {
            STXXL_CHECK(!c.empty());

            c.top_pop();

            progress("Deleting / inserting element", i, (filled ? 3 : 2) * num_elements);
            ++num_deletes;
            STXXL_CHECK_EQUAL(c.size(), num_inserts - num_deletes);
        }
        else {
            int64 this_bulk_size =
                (num_inserts + bulk_size > (filled ? 2 : 1) * num_elements)
                ? num_elements % bulk_size
                : bulk_size;

            c.bulk_push_begin(this_bulk_size);

#if STXXL_PARALLEL
#pragma omp parallel if(parallel)
#endif
            {
                const unsigned int thread_num =
                    parallel ? get_thread_num() : g_rand() % g_max_threads;

                unsigned int seed = thread_num * static_cast<unsigned>(i) * _seed;
                random_number32_r datarng(seed);

#if STXXL_PARALLEL
                #pragma omp for
#endif

                for (int64 j = 0; j < this_bulk_size; j++)
                {
                    progress("Inserting / deleting element",
                             i + j, (filled ? 3 : 2) * num_elements);

                    int64 k = datarng() % value_universe_size;
                    c.bulk_push(value_type(k), thread_num);
                }
            }

            c.bulk_push_end();

            num_inserts += this_bulk_size;
            i += this_bulk_size - 1;
            STXXL_CHECK_EQUAL(c.size(), num_inserts - num_deletes);
        }
    }
}

template <typename ContainerType>
void do_bulk_intermixed_check(ContainerType& c, bool parallel)
{
    std::string parallel_str = parallel ? " in parallel" : "";

    int64 num_inserts = 0, num_deletes = 0, num_failures = 0;
    std::vector<bool> extracted_values(num_elements, false);

    {
        scoped_stats stats(c, "intermix-bulk-check",
                           c.name() + ": Intermixed parallel bulk insert"
                           + parallel_str + " and delete",
                           2);

        for (int64 i = 0; i < 2 * num_elements; ++i)
        {
            int64 r = g_rand() % bulk_size;

            if (num_deletes < num_inserts && num_deletes < num_elements &&
                (r > 0 || num_inserts >= num_elements))
            {
                STXXL_CHECK(!c.empty());

                value_type top = c.top_pop();

                if (top.key < 1 || top.key > num_elements) {
                    std::cout << top.key << " has never been inserted." << std::endl;
                    ++num_failures;
                }
                else if (extracted_values[top.key - 1]) {
                    std::cout << top.key << " has already been extracted." << std::endl;
                    ++num_failures;
                }
                else {
                    extracted_values[top.key - 1] = true;
                }

                progress("Deleting / inserting element", i, 2 * num_elements);
                ++num_deletes;
                STXXL_CHECK_EQUAL(c.size(), num_inserts - num_deletes);
            }
            else
            {
                int64 this_bulk_size = (num_inserts + bulk_size > num_elements)
                                       ? num_elements % bulk_size
                                       : bulk_size;

                c.bulk_push_begin(this_bulk_size);
#if STXXL_PARALLEL
#pragma omp parallel if(parallel)
#endif
                {
                    const unsigned int thread_num =
                        parallel ? get_thread_num() : g_rand() % g_max_threads;
#if STXXL_PARALLEL
                    #pragma omp for
#endif

                    for (int64 j = 0; j < this_bulk_size; j++)
                    {
                        int64 k = num_elements - num_inserts - j;
                        progress("Inserting/deleting element", i + j, 2 * num_elements);

                        c.bulk_push(value_type(k), thread_num);
                    }
                }

                c.bulk_push_end();

                num_inserts += this_bulk_size;
                i += this_bulk_size - 1;
                STXXL_CHECK_EQUAL(c.size(), num_inserts - num_deletes);
            }
        }
    }
    {
        scoped_print_timer stats("Check if all values have been extracted");
        for (int64 i = 0; i < num_elements; ++i) {
            if (!extracted_values[i]) {
                std::cout << i + 1 << " has never been extracted." << std::endl;
                ++num_failures;
            }
        }
    }
    std::cout << num_failures << " failures." << std::endl;
}

#endif

// *****************************************************************************
// Bulk Limit and Pop/Push Operations

template <typename ContainerType>
void do_bulk_prefill(ContainerType& c, bool parallel,
                     typename ContainerType::value_type::key_type& windex)
{
    typedef typename ContainerType::value_type value_type;

    // fill PQ with num_elements items

    scoped_stats stats(c, g_testset + "|prefill", c.name() + ": prefill", 1);

    c.bulk_push_begin(num_elements / 2);

#if STXXL_PARALLEL
#pragma omp parallel if(parallel)
#endif
    {
        const unsigned int thread_num =
            parallel ? get_thread_num() : 0;
#if STXXL_PARALLEL
#pragma omp for
#endif
        for (int64 i = 0; i < (int64)num_elements; ++i)
            c.bulk_push(value_type(i), thread_num);
    }

    c.bulk_push_end();

    windex += num_elements;

    STXXL_CHECK_EQUAL(c.size(), num_elements);
}

template <typename ContainerType>
void do_bulk_postempty(ContainerType& c, bool /* parallel */,
                       typename ContainerType::value_type::key_type& rindex)
{
    typedef typename ContainerType::value_type value_type;

    // extract remaining items

    scoped_stats stats(c, g_testset + "|postempty",
                       c.name() + ": post empty", 1);

    STXXL_CHECK_EQUAL(c.size(), num_elements);

    std::vector<value_type> work;
    while (!c.empty())
    {
        c.bulk_pop(work, std::min<size_t>(bulk_size, c.size()));

        for (int64 j = 0; j < (int64)work.size(); ++j)
        {
            STXXL_CHECK_EQUAL(work[j].key, rindex);
            ++rindex;
        }
    }

    STXXL_CHECK(c.empty());
}

template <bool RandomizeBulkSize, typename ContainerType>
void do_bulk_limit(ContainerType& c, bool parallel)
{
    typedef typename ContainerType::value_type value_type;
    typedef typename value_type::key_type key_type;

    key_type windex = 0; // continuous insertion index
    key_type rindex = 0; // continuous pop index

    do_bulk_prefill(c, parallel, windex);

    // extract bulk_size items, and reinsert them with higher indexes
    {
        scoped_stats stats(c, g_testset + "|cycles",
                           c.name() + ": running bulk-limit cycles", 1);

        size_t cycles = 0;

        for (uint64 elem = 0; elem < num_elements; ++cycles)
        {
            uint64 this_bulk_size = RandomizeBulkSize
                                    ? g_rand() % bulk_size
                                    : bulk_size;

            this_bulk_size = std::min(num_elements, this_bulk_size);
            elem += this_bulk_size * 2;

            if (0) // simple test procedure
            {
                for (uint64 i = 0; i < this_bulk_size; ++i)
                {
                    value_type top = c.top_pop();
                    STXXL_CHECK_EQUAL(top.key, rindex);
                    ++rindex;

                    c.push(value_type(windex++));
                }
            }
            else // bulk-limit procedure
            {
                value_type limit = value_type(rindex + this_bulk_size);
                c.limit_begin(limit, this_bulk_size);

                unsigned iterations = 0;

                while (c.limit_top().key < limit.key)
                {
                    value_type top = c.limit_top_pop();
                    STXXL_CHECK_EQUAL(top.key, rindex);
                    ++rindex;

                    c.limit_push(value_type(windex++), 0);
                    ++iterations;
                }

                STXXL_CHECK_EQUAL(c.limit_top().key, rindex);

                c.limit_end();
                STXXL_CHECK_EQUAL(iterations, this_bulk_size);
            }
        }

        std::cout << "bulk-limit cycles: " << cycles << "\n";
    }

    STXXL_CHECK_EQUAL(c.size(), (size_t)windex - rindex);
    STXXL_CHECK_EQUAL(c.size(), num_elements);

    do_bulk_postempty(c, parallel, rindex);
}

template <bool RandomizeBulkSize, typename ContainerType>
void do_bulk_pop_push(ContainerType& c, bool parallel)
{
    typedef typename ContainerType::value_type value_type;
    typedef typename value_type::key_type key_type;

    key_type windex = 0; // continuous insertion index
    key_type rindex = 0; // continuous pop index

    do_bulk_prefill(c, parallel, windex);

    // extract bulk_size items, and reinsert them with higher indexes
    {
        scoped_stats stats(c, g_testset + "|cycles",
                           c.name() + ": running bulk-pop-push cycles", 1);

        size_t cycles = 0;
        std::vector<value_type> work;

        for (uint64 elem = 0; elem < num_elements; ++cycles)
        {
            uint64 this_bulk_size = RandomizeBulkSize
                                    ? g_rand() % bulk_size
                                    : bulk_size;

            this_bulk_size = std::min(num_elements, this_bulk_size);
            elem += this_bulk_size * 2;

            c.bulk_pop(work, this_bulk_size);
            c.bulk_push_begin(this_bulk_size);

            if (!parallel)
            {
                for (uint64 i = 0; i < work.size(); ++i)
                {
                    STXXL_CHECK_EQUAL(work[i].key, rindex + i);
                    c.bulk_push(value_type(windex + i), 0);
                }
            }
            else // (parallel)
            {
#if STXXL_PARALLEL
#pragma omp parallel
#endif
                {
                    const unsigned int thread_num = get_thread_num();
#if STXXL_PARALLEL
#pragma omp for schedule(static)
#endif
                    for (int64 i = 0; i < (int64)work.size(); ++i)
                    {
                        STXXL_CHECK_EQUAL(work[i].key, rindex + i);
                        c.bulk_push(value_type(windex + i), thread_num);
                    }
                }
            }
            c.bulk_push_end();

            rindex += work.size();
            windex += work.size();
        }

        std::cout << "bulk-pop-push cycles: " << cycles << "\n";
    }

    STXXL_CHECK_EQUAL(c.size(), (size_t)windex - rindex);
    STXXL_CHECK_EQUAL(c.size(), num_elements);

    do_bulk_postempty(c, parallel, rindex);
}

/*
 * Dijkstra Graph SSP Algorithm with different PQs on a random graph.
 */

#if 0
template <typename ContainerType>
class dijkstra_graph
{
public:
    typedef int64 size_type;

protected:
    size_type n, m;
    std::vector<size_type> nodes;
    std::vector<bool> node_visited;
    std::vector<size_type> edge_targets;
    std::vector<size_type> edge_lengths;
    std::vector<std::vector<std::pair<size_type, size_type> > > temp_edges;

public:
    //! construct a random graph for applying Dijkstra's SSP algorithm
    dijkstra_graph(size_type n, size_type m)
        : n(n),
          m(m),
          nodes(n + 1),
          node_visited(n, false),
          edge_targets(m),
          edge_lengths(m),
          temp_edges(n)
    {
        unsigned int seed = 12345;
        random_number32_r datarng(seed);

        for (size_type i = 0; i < m; ++i) {
            size_type source = datarng() % n;
            size_type target = datarng() % n;
            size_type length = datarng();
            temp_edges[source].push_back(std::make_pair(target, length));
        }

        size_type offset = 0;

        for (size_type i = 0; i < n; ++i) {
            nodes[i] = offset;

            for (size_type j = 0; j < temp_edges[i].size(); ++j) {
                edge_targets[offset + j] = temp_edges[i][j].key;
                edge_lengths[offset + j] = temp_edges[i][j].second;
            }

            offset += temp_edges[i].size();
        }

        nodes[n] = offset;
        temp_edges.clear();
    }

    //! run using container c
    void run_random_dijkstra(ContainerType& c)
    {
        unsigned int seed = 678910;
        random_number32_r datarng(seed);

        size_type source = datarng() % n;
        c.push(value_type(source, 0));

        while (!c.empty())
        {
            value_type top = c.top_pop();

            size_type target = top.key;
            size_type distance = top.second;

            if (node_visited[target]) {
                continue;
            }
            else {
                node_visited[target] = true;
            }

            size_type edges_begin = nodes[target];
            size_type edges_end = nodes[target + 1];

            size_type num_edges = edges_end - edges_begin;
            c.bulk_push_begin(num_edges);

#if STXXL_PARALLEL
#pragma omp parallel
#endif
            {
                const unsigned thread_num = get_thread_num();
#if STXXL_PARALLEL
                #pragma omp for
#endif
                for (size_type i = edges_begin; i < edges_end; ++i)
                {
                    c.bulk_push(value_type(edge_targets[i], distance + edge_lengths[i]),
                                thread_num);
                }
            }

            c.bulk_push_end();
        }
    }
};

template <typename ContainerType>
void do_dijkstra(ContainerType& c)
{
    int64 num_nodes = num_elements / value_size;
    int64 num_edges = 100 * num_nodes;

    if (c.name() == "sorter") {
        STXXL_ERRMSG("Dijkstra benchmark does not work with the stxxl::sorter.");
        abort();
    }

    std::cout << "Generating Graph" << std::endl;

    dijkstra_graph<ContainerType> dg(num_nodes, num_edges);

    scoped_stats stats(c, "dijkstra", "Running random dijkstra");
    dg.run_random_dijkstra(c);
}
#endif

/*
 * Print statistics and parameters.
 */

void print_params()
{
    STXXL_MEMDUMP(RAM);
    STXXL_MEMDUMP(volume);
    STXXL_MEMDUMP(block_size);
    STXXL_MEMDUMP(single_heap_ram);
    STXXL_VARDUMP(bulk_size);
    STXXL_MEMDUMP(extract_buffer_ram);
    STXXL_VARDUMP(num_elements);
    STXXL_VARDUMP(g_max_threads);
    STXXL_VARDUMP(num_write_buffers);
    STXXL_VARDUMP(num_prefetchers);
    STXXL_VARDUMP(num_read_blocks);
}

/*
 * Run the selected benchmark procedure.
 */

bool opt_bulk = false;
bool opt_intermixed = false;
bool opt_parallel = false;

const char* benchmark_desc =
    "Select benchmark test to run:\n"
    "  push-asc           push ascending sequence\n"
    "  push-rand          push random sequence\n"
    "  push-asc-pop       push-asc and then pop all\n"
    "  push-rand-pop      push-rand and then pop all\n"
    "  push-asc-popcheck  push-asc and then pop and verify all\n"
    "  push-rand-popcheck push-rand and then pop and verify all\n"
    "  bulk-pop-push      intermixed bulk-pop and bulk-push\n"
    "  bulk-limit         intermixed bulk-limit sequence\n"
    "  rbulk-pop-push     intermixed bulk-pop and bulk-push, random bulks\n"
    "  rbulk-limit        intermixed bulk-limit sequence, random bulks\n"
    "  dijkstra           a Dijkstra SSSP test\n";

template <typename ContainerType>
void run_benchmark(ContainerType& c, const std::string& testset)
{
    std::string testdesc = "Running benchmark " + testset + " on " + c.name();

    typedef typename ContainerType::value_type value_type;

    g_testset = testset;
    // calculate_depending_parameters
    num_elements = volume / sizeof(value_type);

    if (testset == "push-asc") {
        scoped_stats stats(c, testset, testdesc, 1);
        if (!opt_bulk)
            do_push_asc(c);
        else
            do_bulk_push_asc(c, opt_parallel);
    }
    else if (testset == "push-rand") {
        scoped_stats stats(c, testset, testdesc, 1);
        if (!opt_bulk)
            do_push_rand(c, g_seed);
        else
            do_bulk_push_rand(c, g_seed, opt_parallel);
    }
    else if (testset == "push-asc-pop") {
        scoped_stats stats(c, testset, testdesc, 2);
        if (!opt_bulk) {
            do_push_asc(c);
            do_pop(c);
        }
        else {
            do_bulk_push_asc(c, opt_parallel);
            do_bulk_pop(c);
        }
    }
    else if (testset == "push-rand-pop") {
        scoped_stats stats(c, testset, testdesc, 2);
        if (!opt_bulk) {
            do_push_rand(c, g_seed);
            do_pop(c);
        }
        else {
            do_bulk_push_rand(c, g_seed, opt_parallel);
            do_bulk_pop(c);
        }
    }
    else if (testset == "push-asc-popcheck") {
        scoped_stats stats(c, testset, testdesc, 2);
        if (!opt_bulk) {
            do_push_asc(c);
            do_pop_asc_check(c);
        }
        else {
            do_bulk_push_asc(c, opt_parallel);
            do_bulk_pop_check_asc(c);
        }
    }
    else if (testset == "push-rand-popcheck") {
        scoped_stats stats(c, testset, testdesc, 4);
        if (!opt_bulk) {
            do_push_rand(c, g_seed);
            do_pop_rand_check(c, g_seed);
        }
        else {
            do_bulk_push_rand(c, g_seed, opt_parallel);
            do_bulk_pop_check_rand(c, g_seed, opt_parallel);
        }
    }
    else if (testset == "bulk-limit") {
        scoped_stats stats(c, testset, testdesc, 4);
        do_bulk_limit<false>(c, opt_parallel);
    }
    else if (testset == "rbulk-limit") {
        scoped_stats stats(c, testset, testdesc, 4);
        do_bulk_limit<true>(c, opt_parallel);
    }
    else if (testset == "bulk-pop-push") {
        scoped_stats stats(c, testset, testdesc, 4);
        do_bulk_pop_push<false>(c, opt_parallel);
    }
    else if (testset == "rbulk-pop-push") {
        scoped_stats stats(c, testset, testdesc, 4);
        do_bulk_pop_push<true>(c, opt_parallel);
    }
    else if (testset == "dijkstra") {
        //do_dijkstra(c);
    }
    else {
        STXXL_ERRMSG("Unknown benchmark test '" << testset << "'");
    }
}

/*
 * The main procedure.
 */

int main(int argc, char* argv[])
{
    /*
     * Parse flags
     */

    std::string opt_queue = "ppq";
    std::string opt_benchmark = "undefined";

    stxxl::cmdline_parser cp;
    cp.set_description(description);

    cp.add_param_string(
        "pqtype", opt_queue,
        "Select priority queue to test:\n"
        "  ppq - stxxl::parallel_priority_queue (default)\n"
        "  spq - stxxl::priority_queue\n"
        "  sorter - stxxl::sorter\n"
        "  stl - std::priority_queue");

    cp.add_flag('p', "parallel", opt_parallel, "Insert bulks in parallel");
    cp.add_flag('b', "bulk", opt_bulk, "Use bulk insert");
    cp.add_flag('i', "intermixed", opt_intermixed, "Intermixed insert/delete");

    cp.add_param_string("benchmark", opt_benchmark, benchmark_desc);

    cp.add_bytes('v', "volume", volume,
                 "Amount of data to insert in Bytes (not recommended for STXXL PQ)");

    cp.add_bytes('m', "ram", RAM,
                 "Amount of main memory in Bytes (not possible for STXXL PQ)");

    cp.add_bytes('k', "bulk_size", (uint64&)bulk_size,
                 "Number of elements per bulk");

    cp.add_bytes('u', "universe_size", value_universe_size,
                 "Range for random values");

    cp.add_uint('w', "write_buffers", num_write_buffers,
                "Number of buffered blocks when writing to external memory");

    cp.add_bytes('x', "extract_buffer", extract_buffer_ram,
                 "Maximum memory consumption of the extract buffer in Bytes");

    cp.add_bytes('y', "heap_ram", single_heap_ram,
                 "Size of a single insertion heap in Bytes");

    cp.add_double('z', "prefetchers", num_prefetchers,
                  "Number of prefetched blocks for each external array");

    if (!cp.process(argc, argv))
        return -1;

    ram_write_buffers = num_write_buffers * block_size;

    print_params();

#if !STXXL_PARALLEL
    if (opt_parallel) {
        STXXL_MSG("STXXL is compiled without STXXL_PARALLEL flag. "
                  "Parallel processing cannot be used.");
        return EXIT_FAILURE;
    }
#else
    if (!omp_get_nested()) {
        omp_set_nested(1);
        if (!omp_get_nested()) {
            STXXL_ERRMSG("Could not enable OpenMP's nested parallelism, "
                         "however, the PPQ requires this OpenMP feature.");
            abort();
        }
    }

    unsigned thread_initialize = 0;
    #pragma omp parallel
    {
        thread_initialize++;
    }
#endif

    stxxl::stats_data stats_begin(*stxxl::stats::get_instance());

    if (opt_queue == "ppq")
    {
        typedef CStxxlParallePQ<my8_type> ppq_type;
        ppq_type* ppq = new ppq_type();
        run_benchmark(*ppq, opt_benchmark);
        delete ppq;
    }
    else if (opt_queue == "spq")
    {
        typedef CStxxlPQ<my8_type> spq_type;
        spq_type* spq = new spq_type();

        STXXL_DEBUG1("PQ parameters:" <<
                     " total_memory=" << spq->mem_cons() <<
                     " delete_buffer_size=" << spq->delete_buffer_size <<
                     " N=" << spq->N <<
                     " IntKMAX=" << spq->IntKMAX <<
                     " num_int_groups=" << spq->num_int_groups <<
                     " num_ext_groups=" << spq->num_ext_groups <<
                     " total_num_groups=" << spq->total_num_groups <<
                     " BlockSize=" << spq->BlockSize <<
                     " ExtKMAX=" << spq->ExtKMAX);

        run_benchmark(*spq, opt_benchmark);
        delete spq;
    }
    else if (opt_queue == "sorter")
    {
        CStxxlSorter<my8_type> sorter;
        run_benchmark(sorter, opt_benchmark);
    }
    else if (opt_queue == "stl")
    {
        CStdPriorityQueue<my8_type> stlpq;
        run_benchmark(stlpq, opt_benchmark);
    }
#if 1
    else if (opt_queue == "ppq24")
    {
        typedef CStxxlParallePQ<my24_type> ppq_type;
        ppq_type* ppq = new ppq_type();
        run_benchmark(*ppq, opt_benchmark);
        delete ppq;
    }
    else if (opt_queue == "spq24")
    {
        typedef CStxxlPQ<my24_type> spq_type;
        spq_type* spq = new spq_type();
        run_benchmark(*spq, opt_benchmark);
        delete spq;
    }
    else if (opt_queue == "sorter24")
    {
        CStxxlSorter<my24_type> sorter;
        run_benchmark(sorter, opt_benchmark);
    }
    else if (opt_queue == "stl24")
    {
        CStdPriorityQueue<my24_type> stlpq;
        run_benchmark(stlpq, opt_benchmark);
    }
#endif
    else {
        STXXL_MSG("Please choose a queue container type. Use -h for help.");
        return EXIT_FAILURE;
    }

    std::cout << stxxl::stats_data(*stxxl::stats::get_instance()) - stats_begin;

    return EXIT_SUCCESS;
}
