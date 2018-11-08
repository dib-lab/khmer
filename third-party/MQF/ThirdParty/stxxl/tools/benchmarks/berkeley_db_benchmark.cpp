/***************************************************************************
 *  tools/benchmarks/berkeley_db_benchmark.cpp
 *
 *  Part of the STXXL. See http://stxxl.sourceforge.net
 *
 *  Copyright (C) 2006 Roman Dementiev <dementiev@ira.uka.de>
 *  Copyright (C) 2009, 2010 Andreas Beckmann <beckmann@cs.uni-frankfurt.de>
 *
 *  Distributed under the Boost Software License, Version 1.0.
 *  (See accompanying file LICENSE_1_0.txt or copy at
 *  http://www.boost.org/LICENSE_1_0.txt)
 **************************************************************************/

//! \example containers/berkeley_db_benchmark.cpp
//! This is a benchmark mentioned in the paper
//! R. Dementiev, L. Kettner, P. Sanders "STXXL: standard template library for XXL data sets"
//! Software: Practice and Experience
//! Volume 38, Issue 6, Pages 589-637, May 2008
//! DOI: 10.1002/spe.844

#include <stxxl/vector>
#include <stxxl/map>
#include <stxxl/timer>
#include <stxxl/stream>

///// BDB header ////////////
#include <db_cxx.h>

///////// TPIE headers //////////
#include "app_config.h"
#include <portability.h>
#include <ami_btree.h>
////////////////////////////

#define BDB_BULK_SCAN

#define KEY_SIZE                8
#define DATA_SIZE               32

#define NODE_BLOCK_SIZE         (32 * 1024)
#define LEAF_BLOCK_SIZE         (32 * 1024)

#define LEAF_BLOCK_SIZE         (32 * 1024)

#define TOTAL_CACHE_SIZE        (750 * 1024 * 1024)
//#define TOTAL_CACHE_SIZE      (150 * 1024 * 1024)

//#define NODE_CACHE_SIZE       (1 * (TOTAL_CACHE_SIZE / 40))
//#define LEAF_CACHE_SIZE       (39 * (TOTAL_CACHE_SIZE / 40))

#define NODE_CACHE_SIZE         (64 * 1024 * 1024)
#define LEAF_CACHE_SIZE         (TOTAL_CACHE_SIZE - NODE_CACHE_SIZE)

#define SORTER_MEM              (TOTAL_CACHE_SIZE - 1024 * 1024 * 2 * 4)

#define SCAN_LIMIT(x)   (x)

//#define BDB_FILE "/data3/bdb_file"
#define BDB_FILE "/var/tmp/bdb_file"

// BDB settings
u_int32_t pagesize = LEAF_BLOCK_SIZE;
u_int bulkbufsize = 4 * 1024 * 1024;
u_int logbufsize = 8 * 1024 * 1024;
u_int cachesize = (TOTAL_CACHE_SIZE < 500 * 1024 * 1024) ? (4 * (TOTAL_CACHE_SIZE / 5)) : (TOTAL_CACHE_SIZE - 100 * 1024 * 1024);
u_int datasize = DATA_SIZE;
u_int keysize = KEY_SIZE;
u_int numitems = 0;

const char* letters = "abcdefghijklmnopqrstuvwxuz";

struct my_key
{
    char keybuf[KEY_SIZE];
};

std::ostream& operator << (std::ostream& o, const my_key& obj)
{
    for (int i = 0; i < KEY_SIZE; ++i)
        o << obj.keybuf[i];

    return o;
}

bool operator == (const my_key& a, const my_key& b)
{
    return strncmp(a.keybuf, b.keybuf, KEY_SIZE) == 0;
}

bool operator != (const my_key& a, const my_key& b)
{
    return strncmp(a.keybuf, b.keybuf, KEY_SIZE) != 0;
}

bool operator < (const my_key& a, const my_key& b)
{
    return strncmp(a.keybuf, b.keybuf, KEY_SIZE) < 0;
}

bool operator > (const my_key& a, const my_key& b)
{
    return strncmp(a.keybuf, b.keybuf, KEY_SIZE) > 0;
}

bool operator <= (const my_key& a, const my_key& b)
{
    return strncmp(a.keybuf, b.keybuf, KEY_SIZE) <= 0;
}

bool operator >= (const my_key& a, const my_key& b)
{
    return strncmp(a.keybuf, b.keybuf, KEY_SIZE) >= 0;
}

struct my_data
{
    char databuf[DATA_SIZE];
};

std::ostream& operator << (std::ostream& o, const my_data& obj)
{
    o << "DATA(size=" << sizeof(obj) << ") ";

    return o;
}

my_key min_key, max_key;

struct comp_type : std::binary_function<my_key, my_key, bool>
{
    bool operator () (const my_key& a, const my_key& b) const
    {
        return strncmp(a.keybuf, b.keybuf, KEY_SIZE) < 0;
    }
    static my_key max_value()
    {
        return max_key;
    }
    static my_key min_value()
    {
        return min_key;
    }
};

/// TPIE  declarations
// Key type.
typedef my_key bkey_t;

// Element type for the btree.
struct el_t {
    bkey_t key_;
    my_data data_;
    el_t(bkey_t k) : key_(k) { }
    el_t() { }
};
struct key_from_el {
    bkey_t operator () (const el_t& v) const
    {
        return v.key_;
    }
};

// Temporary distinction btw UN*X and WIN, since there are some
// problems with the MMAP collection implementation.
#ifdef _WIN32
typedef AMI_btree<bkey_t, el_t, less<bkey_t>, key_from_el, BTE_COLLECTION_UFS> u_btree_t;
#else
typedef AMI_btree<bkey_t, el_t, less<bkey_t>, key_from_el> u_btree_t;
#endif

void init()
{
    memset(max_key.keybuf, std::numeric_limits<unsigned char>::max(), KEY_SIZE);
    memset(min_key.keybuf, std::numeric_limits<unsigned char>::min(), KEY_SIZE);
}

typedef stxxl::map<my_key, my_data, comp_type, NODE_BLOCK_SIZE, LEAF_BLOCK_SIZE> map_type;

#define REAL_NODE_BLOCK_SIZE map_type::node_block_type::raw_size
#define REAL_LEAF_BLOCK_SIZE map_type::leaf_block_type::raw_size
#define REAL_NODE_MELEMENTS map_type::node_block_type::size
#define REAL_LEAF_MELEMENTS map_type::leaf_block_type::size

typedef stxxl::VECTOR_GENERATOR<std::pair<my_key, my_data>, 1, 1>::result vector_type;
//typedef stxxl::vector<std::pair<my_key,my_data>,1,stxxl::lru_pager<1>,512*1024>  vector_type;

//#define KEYPOS        (i % KEY_SIZE)
//#define VALUE         (myrand() % 26)

#if 0
unsigned ran32State = 0xdeadbeef;
inline unsigned myrand()
{
    return (ran32State = 1664525 * ran32State + 1013904223);
}
inline void rand_key(stxxl::int64 pos, my_key& Key)
{
    for (int i = 0; i < KEY_SIZE; ++i)
        Key.keybuf[i] = letters[myrand() % 26];
}
#else // a longer pseudo random sequence
long long unsigned ran32State = 0xdeadbeef;
inline long long unsigned myrand()
{
    return (ran32State = (ran32State * 0x5DEECE66DULL + 0xBULL) & 0xFFFFFFFFFFFFULL);
}
inline void rand_key(stxxl::int64 /*pos*/, my_key& Key)
{
    long long unsigned r = myrand();
    for (int i = 0; i < KEY_SIZE; ++i)
    {
        Key.keybuf[i] = letters[r % 26];
        r >>= 5;
    }
}
#endif

void run_bdb_btree(stxxl::int64 ops)
{
    const char* filename = BDB_FILE;

    my_key key1_storage;
    my_data data1_storage;

    unlink(filename);

    memset(key1_storage.keybuf, 'a', KEY_SIZE);
    memset(data1_storage.databuf, 'b', DATA_SIZE);

    Db db(NULL, 0);             // Instantiate the Db object

    try {
        db.set_errfile(stderr);
        db.set_pagesize(pagesize);
        db.set_cachesize(0, cachesize, 1);

        // Open the database
        db.open(NULL,           // Transaction pointer
                filename,       // Database file name
                NULL,           // Optional logical database name
                DB_BTREE,       // Database access method
                DB_CREATE,      // Open flags
                0);             // File mode (using defaults)

        // here we start with the tests
        Dbt key1(key1_storage.keybuf, KEY_SIZE);
        Dbt data1(data1_storage.databuf, DATA_SIZE);

        stxxl::timer Timer;
        stxxl::int64 n_inserts = ops, n_locates = ops, n_range_queries = ops, n_deletes = ops;
        stxxl::int64 i;
        //comp_type cmp_;

        ran32State = 0xdeadbeef;

        DB_BTREE_STAT* dbstat;

        db.stat(NULL, &dbstat, 0);
        STXXL_MSG("Records in map: " << dbstat->bt_ndata);

        db.get_env()->memp_stat(NULL, NULL, DB_STAT_CLEAR);

        Timer.start();

        for (i = 0; i < n_inserts; ++i)
        {
            //key1_storage.keybuf[KEYPOS] = letters[VALUE];
            rand_key(i, key1_storage);
            db.put(NULL, &key1, &data1, DB_NOOVERWRITE);
        }

        Timer.stop();
        db.stat(NULL, &dbstat, 0);
        STXXL_MSG("Records in map: " << dbstat->bt_ndata);
        STXXL_MSG("Insertions elapsed time: " << (Timer.mseconds() / 1000.) <<
                  " seconds : " << (double(n_inserts) / (Timer.mseconds() / 1000.)) << " key/data pairs per sec");

        db.get_env()->memp_stat_print(DB_STAT_CLEAR);

        /////////////////////////////////////////
        Timer.reset();
        Timer.start();

        Dbc* cursorp;
        db.cursor(NULL, &cursorp, 0);

        for (i = 0; i < n_locates; ++i)
        {
            //key1_storage.keybuf[KEYPOS] = letters[VALUE];
            rand_key(i, key1_storage);
            Dbt keyx(key1_storage.keybuf, KEY_SIZE);
            Dbt datax(data1_storage.databuf, DATA_SIZE);

            cursorp->get(&keyx, &datax, DB_SET_RANGE);
        }

        Timer.stop();
        STXXL_MSG("Locates elapsed time: " << (Timer.mseconds() / 1000.) <<
                  " seconds : " << (double(ops) / (Timer.mseconds() / 1000.)) << " key/data pairs per sec");

        db.get_env()->memp_stat_print(DB_STAT_CLEAR);

        ////////////////////////////////////
        Timer.reset();

        Timer.start();

        stxxl::int64 n_scanned = 0;

        for (i = 0; i < n_range_queries; ++i)
        {
            //key1_storage.keybuf[KEYPOS] = letters[VALUE];
            rand_key(i, key1_storage);
            my_key last_key = key1_storage;
            //key1_storage.keybuf[KEYPOS] = letters[VALUE];
            rand_key(i, key1_storage);
            if (last_key < key1_storage)
                std::swap(last_key, key1_storage);

            Dbt keyx(key1_storage.keybuf, KEY_SIZE);
            Dbt datax(data1_storage.databuf, DATA_SIZE);

            if (cursorp->get(&keyx, &datax, DB_SET_RANGE) == DB_NOTFOUND)
                continue;

            while (*((my_key*)keyx.get_data()) <= last_key)
            {
                ++n_scanned;
                if (cursorp->get(&keyx, &datax, DB_NEXT) == DB_NOTFOUND)
                    break;
            }

            if (n_scanned >= 10 * n_range_queries)
            {
                ++i;
                break;
            }
        }

        n_range_queries = i;

        Timer.stop();
        if (cursorp != NULL)
            cursorp->close();

        STXXL_MSG("Range query elapsed time: " << (Timer.mseconds() / 1000.) <<
                  " seconds : " << (double(n_scanned) / (Timer.mseconds() / 1000.)) <<
                  " key/data pairs per sec, #queries " << n_range_queries << " #scanned elements: " << n_scanned);

        db.get_env()->memp_stat_print(DB_STAT_CLEAR);

        //////////////////////////////////////

        ran32State = 0xdeadbeef;
        memset(key1_storage.keybuf, 'a', KEY_SIZE);

        Timer.reset();
        Timer.start();

        for (i = 0; i < n_deletes; ++i)
        {
            //key1_storage.keybuf[KEYPOS] = letters[VALUE];
            rand_key(i, key1_storage);
            Dbt keyx(key1_storage.keybuf, KEY_SIZE);
            db.del(NULL, &keyx, 0);
        }

        Timer.stop();
        db.stat(NULL, &dbstat, 0);
        STXXL_MSG("Records in map: " << dbstat->bt_ndata);
        STXXL_MSG("Erase elapsed time: " << (Timer.mseconds() / 1000.) <<
                  " seconds : " << (double(ops) / (Timer.mseconds() / 1000.)) << " key/data pairs per sec");

        db.get_env()->memp_stat_print(DB_STAT_CLEAR);

        db.close(0);
    }
    catch (DbException& e)
    {
        STXXL_ERRMSG("DbException happened");
    }
    catch (std::exception& e)
    {
        STXXL_ERRMSG("std::exception happened");
    }

    unlink(filename);
}

void run_stxxl_map(stxxl::int64 ops)
{
    map_type Map(NODE_CACHE_SIZE, LEAF_CACHE_SIZE);
    Map.enable_prefetching();
    stxxl::stats* Stats = stxxl::stats::get_instance();

    std::pair<my_key, my_data> element;

    memset(element.first.keybuf, 'a', KEY_SIZE);
    memset(element.second.databuf, 'b', DATA_SIZE);

    stxxl::timer Timer;
    stxxl::int64 n_inserts = ops, n_locates = ops, n_range_queries = ops, n_deletes = ops;
    stxxl::int64 i;
    //comp_type cmp_;

    ran32State = 0xdeadbeef;

    //stxxl::random_number32 myrand;

    STXXL_MSG("Records in map: " << Map.size());

    stxxl::stats_data stats_begin(*Stats);
    Timer.start();

    for (i = 0; i < n_inserts; ++i)
    {
        //element.first.keybuf[KEYPOS] = letters[VALUE];
        rand_key(i, element.first);
        Map.insert(element);
    }

    Timer.stop();

    STXXL_MSG("Records in map: " << Map.size());
    STXXL_MSG("Insertions elapsed time: " << (Timer.mseconds() / 1000.) <<
              " seconds : " << (double(n_inserts) / (Timer.mseconds() / 1000.)) << " key/data pairs per sec");

    std::cout << (stxxl::stats_data(*Stats) - stats_begin);

    ////////////////////////////////////////////////

    const map_type& CMap(Map);      // const map reference

    stats_begin = stxxl::stats_data(*Stats);
    Timer.reset();
    Timer.start();

    for (i = 0; i < n_locates; ++i)
    {
        //element.first.keybuf[KEYPOS] = letters[VALUE];
        rand_key(i, element.first);
        CMap.lower_bound(element.first);
    }

    Timer.stop();
    STXXL_MSG("Locates elapsed time: " << (Timer.mseconds() / 1000.) <<
              " seconds : " << (double(n_locates) / (Timer.mseconds() / 1000.)) << " key/data pairs per sec");

    std::cout << (stxxl::stats_data(*Stats) - stats_begin);

    ////////////////////////////////////

    stats_begin = stxxl::stats_data(*Stats);
    Timer.reset();
    Timer.start();

    stxxl::int64 n_scanned = 0; //, skipped_qieries = 0;

    for (i = 0; i < n_range_queries; ++i)
    {
        //element.first.keybuf[KEYPOS] = letters[VALUE];
        rand_key(i, element.first);
        my_key begin_key = element.first;
        map_type::const_iterator begin = CMap.lower_bound(element.first);
        //element.first.keybuf[KEYPOS] = letters[VALUE];
        rand_key(i, element.first);
        map_type::const_iterator beyond = CMap.lower_bound(element.first);
        if (element.first < begin_key)
            std::swap(begin, beyond);

        while (begin != beyond)
        {
            my_data tmp = begin->second;
            stxxl::STXXL_UNUSED(tmp);
            ++n_scanned;
            ++begin;
        }
        if (n_scanned >= 10 * n_range_queries)
        {
            ++i;
            break;
        }
    }

    n_range_queries = i;

    Timer.stop();
    STXXL_MSG("Range query elapsed time: " << (Timer.mseconds() / 1000.) <<
              " seconds : " << (double(n_scanned) / (Timer.mseconds() / 1000.)) <<
              " key/data pairs per sec, #queries " << n_range_queries << " #scanned elements: " << n_scanned);

    std::cout << (stxxl::stats_data(*Stats) - stats_begin);

    //////////////////////////////////////
    ran32State = 0xdeadbeef;
    memset(element.first.keybuf, 'a', KEY_SIZE);
    memset(element.second.databuf, 'b', DATA_SIZE);

    stats_begin = stxxl::stats_data(*Stats);
    Timer.reset();
    Timer.start();

    for (i = 0; i < n_deletes; ++i)
    {
        //element.first.keybuf[KEYPOS] = letters[VALUE];
        rand_key(i, element.first);
        Map.erase(element.first);
    }

    Timer.stop();
    STXXL_MSG("Records in map: " << Map.size());
    STXXL_MSG("Erase elapsed time: " << (Timer.mseconds() / 1000.) <<
              " seconds : " << (double(n_deletes) / (Timer.mseconds() / 1000.)) << " key/data pairs per sec");

    std::cout << (stxxl::stats_data(*Stats) - stats_begin);
}

class rand_key_gen
{
    stxxl::int64 counter;
    my_key& current;
    stxxl::random_number32 myrand;
    rand_key_gen();

public:
    typedef my_key value_type;

    rand_key_gen(stxxl::int64 el, my_key& cur)
        : counter(el), current(cur)
    {
        //const stxxl::int64  & i = counter;
        //current.keybuf[KEYPOS] = letters[VALUE];
        rand_key(counter, current);
    }
    const my_key& operator * () { return current; }
    const my_key* operator -> () { return &current; }

    rand_key_gen& operator ++ ()
    {
        --counter;
        //const stxxl::int64  & i = counter;
        //current.keybuf[KEYPOS] = letters[VALUE];
        rand_key(counter, current);
        return *this;
    }
    bool empty() const { return counter == 0; }
};

template <class InputType>
class key2pair
{
    InputType& in;
    std::pair<my_key, my_data> current;
    key2pair();

public:
    typedef std::pair<my_key, my_data> value_type;

    key2pair(InputType& in_) : in(in_)
    {
        if (!in.empty())
            current.first = *in;
    }
    const value_type& operator * () { return current; }
    const value_type* operator -> () { return &current; }

    key2pair& operator ++ ()
    {
        ++in;
        if (!empty())
            current.first = *in;

        return *this;
    }
    bool empty() const { return in.empty(); }
};

void run_stxxl_map_big(stxxl::int64 n, unsigned ops)
{
    stxxl::stats* Stats = stxxl::stats::get_instance();

    std::pair<my_key, my_data> element;

    memset(element.first.keybuf, 'a', KEY_SIZE);
    memset(element.second.databuf, 'b', DATA_SIZE);

    stxxl::timer Timer;
    stxxl::int64 n_inserts = ops, n_locates = ops, n_range_queries = ops, n_deletes = ops;
    stxxl::int64 i;
    //comp_type cmp_;

    ran32State = 0xdeadbeef;

    //stxxl::random_number32 myrand;

    Timer.start();

    vector_type SortedSeq(n);
    const vector_type& CSortedSeq(SortedSeq);
    {
        rand_key_gen Gen(n, element.first);
        typedef stxxl::stream::sort<rand_key_gen, comp_type> sorter_type;
        sorter_type Sorter(Gen, comp_type(), SORTER_MEM);
        typedef key2pair<sorter_type> key2pair_type;
        key2pair_type Key2Pair(Sorter);
        stxxl::stream::materialize(Key2Pair, SortedSeq.begin());
    }

    Timer.stop();

    STXXL_MSG("Finished sorting input. Elapsed time: " <<
              (Timer.mseconds() / 1000.) << " seconds.");

    stxxl::stats_data stats_begin(*Stats);
    Timer.reset();
    Timer.start();

    // bulk construction
    map_type Map(CSortedSeq.begin(),
                 CSortedSeq.end(),
                 NODE_CACHE_SIZE, LEAF_CACHE_SIZE, true);

    Timer.stop();

    STXXL_MSG("Records in map: " << Map.size());
    STXXL_MSG("Construction elapsed time: " << (Timer.mseconds() / 1000.) <<
              " seconds : " << (double(n) / (Timer.mseconds() / 1000.)) << " key/data pairs per sec");

    Map.print_statistics(cout);
    Map.reset_statistics();
    std::cout << (stxxl::stats_data(*Stats) - stats_begin);

    ////////////////////////////////////////

    Map.disable_prefetching();

    stats_begin = stxxl::stats_data(*Stats);
    Timer.reset();
    Timer.start();

    for (i = 0; i < n_inserts; ++i)
    {
        //element.first.keybuf[KEYPOS] = letters[VALUE];
        rand_key(i, element.first);
        Map.insert(element);
    }

    Timer.stop();

    STXXL_MSG("Records in map: " << Map.size());
    STXXL_MSG("Insertions elapsed time: " << (Timer.mseconds() / 1000.) <<
              " seconds : " << (double(n_inserts) / (Timer.mseconds() / 1000.)) << " key/data pairs per sec");

    Map.print_statistics(cout);
    Map.reset_statistics();
    std::cout << (stxxl::stats_data(*Stats) - stats_begin);

    ////////////////////////////////////

    const map_type& CMap(Map);      // const map reference

    Timer.reset();
    Timer.start();

    for (i = 0; i < n_locates; ++i)
    {
        //element.first.keybuf[KEYPOS] = letters[VALUE];
        rand_key(i, element.first);
        CMap.lower_bound(element.first);
    }

    Timer.stop();
    STXXL_MSG("Locates elapsed time: " << (Timer.mseconds() / 1000.) <<
              " seconds : " << (double(ops) / (Timer.mseconds() / 1000.)) << " key/data pairs per sec");

    Map.print_statistics(cout);
    Map.reset_statistics();
    std::cout << (stxxl::stats_data(*Stats) - stats_begin);

    ////////////////////////////////////

    Map.enable_prefetching();

    stats_begin = stxxl::stats_data(*Stats);
    Timer.reset();
    Timer.start();

    stxxl::int64 n_scanned = 0; //, skipped_qieries = 0;

    for (i = 0; i < n_range_queries; ++i)
    {
        //element.first.keybuf[KEYPOS] = letters[VALUE];
        rand_key(i, element.first);
        my_key begin_key = element.first;
        map_type::const_iterator begin = CMap.lower_bound(element.first);
        //element.first.keybuf[KEYPOS] = letters[VALUE];
        rand_key(i, element.first);
        map_type::const_iterator beyond = CMap.lower_bound(element.first);
        if (element.first < begin_key)
            std::swap(begin, beyond);

        /*
           STXXL_MSG("Looking     "<<element.first<<" scanned: "<<n_scanned);

           if(beyond==CMap.end())
           {
                STXXL_MSG("Upper bound "<<"END");
           }
           else
           {
                STXXL_MSG("Upper bound "<<((element.first>begin_key)?element.first:begin_key));
           }*/

        while (begin != beyond)
        {
            ++n_scanned;
            ++begin;
        }

        if (n_scanned >= SCAN_LIMIT(n))
        {
            ++i;
            break;
        }
    }

    n_range_queries = i;

    Timer.stop();
    STXXL_MSG("Range query elapsed time: " << (Timer.mseconds() / 1000.) <<
              " seconds : " << (double(n_scanned) / (Timer.mseconds() / 1000.)) <<
              " key/data pairs per sec, #queries " << n_range_queries << " #scanned elements: " << n_scanned);

    Map.print_statistics(cout);
    Map.reset_statistics();
    std::cout << (stxxl::stats_data(*Stats) - stats_begin);

    //////////////////////////////////////

    ran32State = 0xdeadbeef;
    memset(element.first.keybuf, 'a', KEY_SIZE);
    memset(element.second.databuf, 'b', DATA_SIZE);

    Map.disable_prefetching();

    stats_begin = stxxl::stats_data(*Stats);
    Timer.reset();
    Timer.start();

    for (i = n_deletes; i > 0; --i)
    {
        //element.first.keybuf[KEYPOS] = letters[VALUE];
        rand_key(i, element.first);
        Map.erase(element.first);
    }

    Timer.stop();
    STXXL_MSG("Records in map: " << Map.size());
    STXXL_MSG("Erase elapsed time: " << (Timer.mseconds() / 1000.) <<
              " seconds : " << (double(ops) / (Timer.mseconds() / 1000.)) << " key/data pairs per sec");

    Map.print_statistics(cout);
    Map.reset_statistics();
    std::cout << (stxxl::stats_data(*Stats) - stats_begin);
}

/////////////////////////////////////////////////////////////////////////

typedef AMI_STREAM<el_t> stream_t;

char dummy;

class MyFilter
{
public:
    bool operator () (const el_t& v) const
    {
        dummy += v.key_.keybuf[0];         // touch element
        return true;
    }
};

void run_tpie_btree_big(stxxl::int64 n, unsigned ops)
{
    el_t element;

    memset(element.key_.keybuf, 'a', KEY_SIZE);
    memset(element.data_.databuf, 'b', DATA_SIZE);

    // Log debugging info from the application, but not from the library.
    tpie_log_init(TPIE_LOG_APP_DEBUG);
    MM_manager.set_memory_limit(TOTAL_CACHE_SIZE);
    MM_manager.enforce_memory_limit();

    stream_t* is = new stream_t;

    stxxl::timer Timer;
    stxxl::int64 n_inserts = ops, n_locates = ops, n_range_queries = ops, n_deletes = ops;
    stxxl::int64 i;
    //comp_type cmp_;

    ran32State = 0xdeadbeef;

    //stxxl::random_number32 myrand;

    Timer.start();

    {
        rand_key_gen Gen(n, element.key_);
        typedef stxxl::stream::sort<rand_key_gen, comp_type> sorter_type;
        sorter_type Sorter(Gen, comp_type(), SORTER_MEM);
        //typedef key2pair<sorter_type> key2pair_type;
        //key2pair_type Key2Pair(Sorter);
        while (!Sorter.empty())
        {
            is->write_item(*Sorter);
            ++Sorter;
        }
    }

    Timer.stop();

    STXXL_MSG("Finished sorting input. Elapsed time: " <<
              (Timer.mseconds() / 1000.) << " seconds.");

    Timer.reset();
    Timer.start();

    // bulk construction
    u_btree_t* u_btree;
    AMI_btree_params params;
    params.node_block_factor = NODE_BLOCK_SIZE / 4096;
    params.leaf_block_factor = LEAF_BLOCK_SIZE / 4096;
    params.leaf_cache_size = LEAF_CACHE_SIZE / LEAF_BLOCK_SIZE;
    params.node_cache_size = NODE_CACHE_SIZE / NODE_BLOCK_SIZE;

    u_btree = new u_btree_t(params);

    using std::cout;
    using std::cerr;

    if (!u_btree->is_valid()) {
        cerr << "Error initializing btree. Aborting.\n";
        delete u_btree;
        exit(1);
    }

    if (u_btree->load_sorted(is, 1.0, 1.0) != AMI_ERROR_NO_ERROR)
        cerr << "Error during bulk loading.\n";

    Timer.stop();

    STXXL_MSG("Records in map: " << u_btree->size());
    STXXL_MSG("Construction elapsed time: " << (Timer.mseconds() / 1000.) <<
              " seconds : " << (double(n) / (Timer.mseconds() / 1000.)) << " key/data pairs per sec");

    ////////////////////////////////////////
    Timer.reset();

    Timer.start();

    for (i = 0; i < n_inserts; ++i)
    {
        //element.first.keybuf[KEYPOS] = letters[VALUE];
        rand_key(i, element.key_);
        u_btree->insert(element);
    }

    Timer.stop();

    STXXL_MSG("Records in map: " << u_btree->size());
    STXXL_MSG("Insertions elapsed time: " << (Timer.mseconds() / 1000.) <<
              " seconds : " << (double(n_inserts) / (Timer.mseconds() / 1000.)) << " key/data pairs per sec");

    ////////////////////////////////////////////////
    Timer.reset();

    Timer.start();

    el_t result;
    for (i = 0; i < n_locates; ++i)
    {
        //element.first.keybuf[KEYPOS] = letters[VALUE];
        rand_key(i, element.key_);
        u_btree->succ(element.key_, result);
    }

    Timer.stop();
    STXXL_MSG("Locates elapsed time: " << (Timer.mseconds() / 1000.) <<
              " seconds : " << (double(ops) / (Timer.mseconds() / 1000.)) << " key/data pairs per sec");

    ////////////////////////////////////
    Timer.reset();

    Timer.start();

    stxxl::int64 n_scanned = 0; //, skipped_qieries = 0;
    MyFilter filter;

    for (i = 0; i < n_range_queries; ++i)
    {
        rand_key(i, element.key_);
        my_key begin_key = element.key_;
        rand_key(i, element.key_);
        if (element.key_ < begin_key)
            n_scanned += u_btree->range_query(element.key_, begin_key, NULL, filter);

        else
            n_scanned += u_btree->range_query(begin_key, element.key_, NULL, filter);

        if (n_scanned >= SCAN_LIMIT(n))
        {
            ++i;
            break;
        }
    }

    n_range_queries = i;

    Timer.stop();
    STXXL_MSG("Range query elapsed time: " << (Timer.mseconds() / 1000.) <<
              " seconds : " << (double(n_scanned) / (Timer.mseconds() / 1000.)) <<
              " key/data pairs per sec, #queries " << n_range_queries << " #scanned elements: " << n_scanned);

    //////////////////////////////////////
    ran32State = 0xdeadbeef;
    memset(element.key_.keybuf, 'a', KEY_SIZE);
    memset(element.data_.databuf, 'b', DATA_SIZE);

    Timer.reset();

    Timer.start();

    for (i = n_deletes; i > 0; --i)
    {
        //element.first.keybuf[KEYPOS] = letters[VALUE];
        rand_key(i, element.key_);
        u_btree->erase(element.key_);
    }

    Timer.stop();
    STXXL_MSG("Records in map: " << u_btree->size());
    STXXL_MSG("Erase elapsed time: " << (Timer.mseconds() / 1000.) <<
              " seconds : " << (double(ops) / (Timer.mseconds() / 1000.)) << " key/data pairs per sec");
    delete u_btree;
    delete is;
}

void run_bdb_btree_big(stxxl::int64 n, unsigned ops)
{
    const char* filename = BDB_FILE;

    my_key key1_storage;
    my_data data1_storage;

#ifdef BDB_BULK_SCAN
    int* bulk_buffer = new int[logbufsize / sizeof(int)];
#endif

    unlink(filename);

    memset(key1_storage.keybuf, 'a', KEY_SIZE);
    memset(data1_storage.databuf, 'b', DATA_SIZE);

    Db db(NULL, 0);                   // Instantiate the Db object

    try {
        // here we start with the tests
        Dbt key1(key1_storage.keybuf, KEY_SIZE);
        Dbt data1(data1_storage.databuf, DATA_SIZE);

        stxxl::timer Timer;
        stxxl::int64 n_inserts = ops, n_locates = ops, n_range_queries = ops, n_deletes = ops;
        stxxl::int64 i;
        //comp_type cmp_;

        ran32State = 0xdeadbeef;

        Timer.start();

        vector_type SortedSeq(n);
        //const vector_type & CSortedSeq(SortedSeq);
        {
            rand_key_gen Gen(n, key1_storage);
            typedef stxxl::stream::sort<rand_key_gen, comp_type> sorter_type;
            sorter_type Sorter(Gen, comp_type(), SORTER_MEM);
            typedef key2pair<sorter_type> key2pair_type;
            key2pair_type Key2Pair(Sorter);
            stxxl::stream::materialize(Key2Pair, SortedSeq.begin());
        }

        Timer.stop();

        STXXL_MSG("Finished sorting input. Elapsed time: " << (Timer.mseconds() / 1000.) << " seconds.");

        db.set_errfile(stderr);
        db.set_pagesize(pagesize);
        db.set_cachesize(0, cachesize, 10);

        STXXL_MSG("BDB cache size set.");

        // Open the database
        db.open(NULL,           // Transaction pointer
                filename,       // Database file name
                NULL,           // Optional logical database name
                DB_BTREE,       // Database access method
                DB_CREATE,      // Open flags
                0);             // File mode (using defaults)

        db.get_env()->memp_stat(NULL, NULL, DB_STAT_CLEAR);

        Timer.reset();
        Timer.start();

        // DBD does not have bulk construction
        // however inserting in sorted order might help
        // to improve performance
        vector_type::const_iterator cit = SortedSeq.begin();
        for (i = 0; i < n; ++i, ++cit)
        {
            key1_storage = cit->first;
            db.put(NULL, &key1, &data1, DB_NOOVERWRITE);
        }

        Timer.stop();

        DB_BTREE_STAT* dbstat;
        db.stat(NULL, &dbstat, 0);
        STXXL_MSG("Records in map: " << dbstat->bt_ndata);
        STXXL_MSG("Construction elapsed time: " << (Timer.mseconds() / 1000.) <<
                  " seconds : " << (double(n) / (Timer.mseconds() / 1000.)) << " key/data pairs per sec");

        db.stat_print(0);
        db.get_env()->memp_stat_print(DB_STAT_CLEAR);
        ////////////////////////////////////////

        Timer.reset();
        Timer.start();

        for (i = 0; i < n_inserts; ++i)
        {
            //key1_storage.keybuf[KEYPOS] = letters[VALUE];
            rand_key(i, key1_storage);
            db.put(NULL, &key1, &data1, DB_NOOVERWRITE);
        }

        Timer.stop();
        db.stat(NULL, &dbstat, 0);
        STXXL_MSG("Records in map: " << dbstat->bt_ndata);
        STXXL_MSG("Insertions elapsed time: " << (Timer.mseconds() / 1000.) <<
                  " seconds : " << (double(n_inserts) / (Timer.mseconds() / 1000.)) << " key/data pairs per sec");

        db.stat_print(0);
        db.get_env()->memp_stat_print(DB_STAT_CLEAR);

        /////////////////////////////////////////
        Timer.reset();
        Timer.start();

        Dbc* cursorp;
        db.cursor(NULL, &cursorp, 0);

        for (i = 0; i < n_locates; ++i)
        {
            //key1_storage.keybuf[KEYPOS] = letters[VALUE];
            rand_key(i, key1_storage);
            Dbt keyx(key1_storage.keybuf, KEY_SIZE);
            Dbt datax(data1_storage.databuf, DATA_SIZE);

            cursorp->get(&keyx, &datax, DB_SET_RANGE);
        }

        Timer.stop();
        STXXL_MSG("Locates elapsed time: " << (Timer.mseconds() / 1000.) <<
                  " seconds : " << (double(ops) / (Timer.mseconds() / 1000.)) << " key/data pairs per sec");

        db.stat_print(0);
        db.get_env()->memp_stat_print(DB_STAT_CLEAR);

        ////////////////////////////////////
        Timer.reset();

        Timer.start();

        stxxl::int64 n_scanned = 0;

        for (i = 0; i < n_range_queries; ++i)
        {
            //key1_storage.keybuf[KEYPOS] = letters[VALUE];
            rand_key(i, key1_storage);
            my_key last_key = key1_storage;
            //key1_storage.keybuf[KEYPOS] = letters[VALUE];
            rand_key(i, key1_storage);
            if (last_key < key1_storage)
                std::swap(last_key, key1_storage);

            //STXXL_MSG("Looking     "<<key1_storage<<" scanned: "<<n_scanned);
            //STXXL_MSG("Upper bound "<<last_key);

            Dbt keyx(key1_storage.keybuf, KEY_SIZE);
#ifdef BDB_BULK_SCAN
            Dbt datax(bulk_buffer, DATA_SIZE);
            datax.set_ulen(logbufsize);
            datax.set_flags(DB_DBT_USERMEM);
#else
            Dbt datax(data1_storage.databuf, DATA_SIZE);
#endif

#ifdef BDB_BULK_SCAN
            if (cursorp->get(&keyx, &datax, DB_SET_RANGE | DB_MULTIPLE_KEY) == DB_NOTFOUND)
                continue;

            do
            {
                DbMultipleKeyDataIterator BulkIterator(datax);
                Dbt key1, data1;
                while (BulkIterator.next(key1, data1) &&
                       *((my_key*)key1.get_data()) <= last_key)
                {
                    ++n_scanned;
                    //STXXL_MSG("Result      "<<*((my_key *)key1.get_data()));
                }
                if (cursorp->get(&keyx, &datax, DB_NEXT | DB_MULTIPLE_KEY) == DB_NOTFOUND)
                    break;

                if (*((my_key*)keyx.get_data()) > last_key)
                {
                    break;
                }
            } while (1);

#else
            if (cursorp->get(&keyx, &datax, DB_SET_RANGE) == DB_NOTFOUND)
                continue;

            while (*((my_key*)keyx.get_data()) <= last_key)
            {
                ++n_scanned;
                if (cursorp->get(&keyx, &datax, DB_NEXT) == DB_NOTFOUND)
                    break;
            }
#endif

            if (n_scanned >= SCAN_LIMIT(n))
            {
                ++i;
                break;
            }
        }

        n_range_queries = i;

        Timer.stop();
        if (cursorp != NULL)
            cursorp->close();

        STXXL_MSG("Range query elapsed time: " << (Timer.mseconds() / 1000.) <<
                  " seconds : " << (double(n_scanned) / (Timer.mseconds() / 1000.)) <<
                  " key/data pairs per sec, #queries " << n_range_queries << " #scanned elements: " << n_scanned);

        db.stat_print(0);
        db.get_env()->memp_stat_print(DB_STAT_CLEAR);

        //////////////////////////////////////

        ran32State = 0xdeadbeef;
        memset(key1_storage.keybuf, 'a', KEY_SIZE);

        Timer.reset();
        Timer.start();

        for (i = 0; i < n_deletes; ++i)
        {
            //key1_storage.keybuf[KEYPOS] = letters[VALUE];
            rand_key(i, key1_storage);
            Dbt keyx(key1_storage.keybuf, KEY_SIZE);
            db.del(NULL, &keyx, 0);
        }

        Timer.stop();
        db.stat(NULL, &dbstat, 0);
        STXXL_MSG("Records in map: " << dbstat->bt_ndata);
        STXXL_MSG("Erase elapsed time: " << (Timer.mseconds() / 1000.) <<
                  " seconds : " << (double(ops) / (Timer.mseconds() / 1000.)) << " key/data pairs per sec");

        db.stat_print(0);
        db.get_env()->memp_stat_print(DB_STAT_CLEAR);

        db.close(0);
    }
    catch (DbException& e)
    {
        STXXL_ERRMSG("DbException happened");
    }
    catch (std::exception& e)
    {
        STXXL_ERRMSG("std::exception happened");
    }

    unlink(filename);

#ifdef BDB_BULK_SCAN
    delete[]  bulk_buffer;
#endif
}

int main(int argc, char* argv[])
{
    STXXL_MSG("stxxl::map Real Node block size: " << REAL_NODE_BLOCK_SIZE << " bytes");
    STXXL_MSG("stxxl::map Real Leaf block size: " << REAL_LEAF_BLOCK_SIZE << " bytes");
    STXXL_MSG("stxxl::map Node max elements   : " << REAL_NODE_MELEMENTS);
    STXXL_MSG("stxxl::map Leaf max elements   : " << REAL_LEAF_MELEMENTS);
#if STXXL_DIRECT_IO_OFF
    STXXL_MSG("STXXL_DIRECT_IO_OFF is defined");
#else
    STXXL_MSG("STXXL_DIRECT_IO_OFF is NOT defined");
#endif

    if (argc < 3)
    {
        STXXL_MSG("Usage: " << argv[0] << " version #ops");
        STXXL_MSG("\t version = 1: test stxxl map");
        STXXL_MSG("\t version = 2: test Berkeley DB btree");
        STXXL_MSG("\t version = 3: big test stxxl map");
        STXXL_MSG("\t version = 4: big test Berkeley DB btree");
        STXXL_MSG("\t version = 5: big test TPIE btree");
        return -1;
    }

    init();

    int version = atoi(argv[1]);
    stxxl::uint64 ops = stxxl::atouint64(argv[2]);

    STXXL_MSG("Running version      : " << version);
    STXXL_MSG("Operations to perform: " << ops);
    STXXL_MSG("Btree cache size     : " << TOTAL_CACHE_SIZE << " bytes");
    STXXL_MSG("Leaf block size      : " << LEAF_BLOCK_SIZE << " bytes");

    switch (version)
    {
    case 1:
        run_stxxl_map(ops);
        break;
    case 2:
        run_bdb_btree(ops);
        break;
    case 3:
        run_stxxl_map_big(ops, 100000);
        break;
    case 4:
        run_bdb_btree_big(ops, 100000);
        break;
    case 5:
        run_tpie_btree_big(ops, 100000);
        break;
    default:
        STXXL_MSG("Unsupported version " << version);
    }
}
