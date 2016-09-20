#ifndef BASE_FILTER_HH
#define BASE_FILTER_HH

#include <stddef.h>
#include <stdint.h>
#include <string.h>
#include <fstream>
#include <iostream>
#include <list>
#include <map>
#include <queue>
#include <queue>
#include <set>
#include <string>
#include <vector>

#include "khmer.hh"
#include "khmer_exception.hh"
#include "kmer_hash.hh"
#include "read_parsers.hh"
#include "traversal.hh"

namespace khmer
{

template <typename HashFunctorType, typename HashType>
class BaseFilter
{
protected:

    HashFunctorType _hash;

    uint64_t _n_unique;
    uint64_t _occupied_bins;

    std::vector<uint64_t> _tablesizes;
    size_t _n_tables;
    Byte ** _counts;

    explicit BaseFilter(std::vector<uint64_t>& tablesizes,
                        HashFunctorType hash_function)
        : hash_function(hash_function), _tablesizes(tablesizes),
          _n_unique_kmers(0), _occupied_bins(0)
    {
    }

    explicit BaseFilter(uint64_t single_tablesize,
                        HashFunctorType hash_function)
        : hash_function(hash_function), _n_unique_kmers(0), _occupied_bins(0)
    {
        _tablesizes.push_back(single_tablesize);
    }

    explicit BaseFilter(const BaseFilter&);
    BaseFilter& operator=(const BaseFilter&);

    virtual ~BaseFilter( ) {}

    virtual void _allocate_counters() = 0;

public:

    HashFunctorType get_hash_function() {
        return _hash;
    }

    // Accessors for protected/private table info members
    std::vector<uint64_t> get_tablesizes() const
    {
        return _tablesizes;
    }

    const size_t n_tables() const
    {
        return _n_tables;
    }

    // Writing to the tables outside of defined methods has undefined behavior!
    // As such, this should only be used to return read-only interfaces
    Byte ** get_raw_tables()
    {
        return _counts;
    }

    // count number of occupied bins
    virtual const uint64_t n_occupied() const
    {
        return _occupied_bins;
    }

    virtual const uint64_t n_unique() const
    {
        return _n_unique;
    }

    virtual BoundedCounterType test_and_set_bits(const char * str) = 0;
    virtual BoudnedCounterType test_and_set_bits(HashType elem) = 0;

    virtual void count(const char * str) = 0;
    virtual void count(HashType elem) = 0;

    // get the count for the given k-mer.
    virtual const BoundedCounterType get_count(const char * str) const = 0;
    virtual const BoundedCounterType get_count(HashType elem) const = 0;
};


template <typename HashFunctorType, typename HashType>
class BloomFilter: public khmer::BaseFilter<HashFunctorType, HashType>
{
protected:

    virtual void _allocate_counters()
    {
        _n_tables = _tablesizes.size();

        _counts = new Byte*[_n_tables];

        for (size_t i = 0; i < _n_tables; i++) {
            uint64_t tablesize = _tablesizes[i];
            uint64_t tablebytes = tablesize / 8 + 1;

            _counts[i] = new Byte[tablebytes];
            memset(_counts[i], 0, tablebytes);
        }
    }

public:
    BloomFilter(std::vector<HashType>& tablesizes,
                HashFunctorType hash_function)
        : khmer::BaseFilter<HashFunctorType, HashType>(tablesizes, hash_function)
    {
        _allocate_counters();
    }

    ~BloomFilter()
    {
        if (_counts) {
            for (size_t i = 0; i < _n_tables; i++) {
                delete[] _counts[i];
                _counts[i] = NULL;
            }
            delete[] _counts;
            _counts = NULL;

            _n_tables = 0;
        }

    }

    inline
    virtual
    BoundedCounterType
    test_and_set_bits(const char * str)
    {
        return test_and_set_bits(_hash(str));
    }

    // Get and set the hashbits for the given kmer hash.
    // Generally, it is better to keep tests and mutations separate,
    // but, in the interests of efficiency and thread safety,
    // tests and mutations are being blended here against conventional
    // software engineering wisdom.
    inline
    virtual
    BoundedCounterType
    test_and_set_bits( HashType elem )
    {
        bool is_new = false;

        for (size_t i = 0; i < _n_tables; i++) {
            HashIntoType bin = elem % _tablesizes[i];
            HashIntoType byte = bin / 8;
            unsigned char bit = (unsigned char)(1 << (bin % 8));

            unsigned char bits_orig = __sync_fetch_and_or( *(_counts + i) +
                                      byte, bit );
            if (!(bits_orig & bit)) {
                if (i == 0) {
                    __sync_add_and_fetch( &_occupied_bins, 1 );
                }
                is_new = true;
            }
        } // iteration over hashtables

        if (is_new) {
            __sync_add_and_fetch( &_n_unique, 1 );
            return 1; // kmer not seen before
        }

        return 0; // kmer already seen
    } // test_and_set_bits

    virtual void count(const char * str)
    {
        count(_hash(str));
    }

    virtual void count(HashType elem)
    {
        test_and_set_bits(elem);
    }

    // get the count for the given k-mer.
    virtual const BoundedCounterType get_count(const char * str) const
    {
        return get_count(_hash(str));
    }

    // get the count for the given k-mer hash.
    virtual const BoundedCounterType get_count(HashType elem) const
    {
        for (size_t i = 0; i < _n_tables; i++) {
            uint64_t bin = elem % _tablesizes[i];
            uint64_t byte = bin / 8;
            unsigned char bit = bin % 8;

            if (!(_counts[i][byte] & (1 << bit))) {
                return 0;
            }
        }
        return 1;
    }

    void update_from(const BloomFilter<HashFunctorType, HashType> &other) {
        if (_tablesizes != other._tablesizes) {
            throw khmer_exception("both BloomFilters must have same table sizes");
        }
        Byte tmp = 0;
        for (unsigned int table_num = 0; table_num < _n_tables; table_num++) {
            Byte * me = _counts[table_num];
            Byte * ot = other._counts[table_num];
            uint64_t tablesize = _tablesizes[table_num];
            uint64_t tablebytes = tablesize / 8 + 1;

            for (uint64_t index = 0; index < tablebytes; index++) {
                // Bloom filters can be unioned with bitwise OR.
                // First, get the new value
                tmp = me[index] | ot[index];
                if (table_num == 0) {
                    // We'd like for the merged filter to have an accurate
                    // count of occupied bins.  First, observe that
                    // HammingDistance(x,y) is equivalent to
                    // HammingWeight(x^y).  Then, observe that the number
                    // of additional occupied bins from the update is the
                    // hamming distance between the original bin and the
                    // OR'd bin. Thus, we can use the builtin popcountll
                    // function, which calls a hardware instruction for
                    // hamming weight, with the original and merged bin,
                    // to find the number of additional occupied bins.
                    _occupied_bins += __builtin_popcountll(me[index] ^ tmp);
                }
                me[index] = tmp;
            }
        }
    }
};


template <typename HashFunctorType, typename HashType>
class CountingFilter: public khmer::BaseFilter<HashFunctorType, HashType>
{

protected:

    bool _use_bigcount;		// keep track of counts > Bloom filter hash count threshold?
    uint32_t _bigcount_spin_lock;

    unsigned int    _max_count;
    unsigned int    _max_bigcount;

    virtual void _allocate_counters()
    {
        _n_tables = _tablesizes.size();

        _counts = new Byte*[_n_tables];
        for (size_t i = 0; i < _n_tables; i++) {
            _counts[i] = new Byte[_tablesizes[i]];
            memset(_counts[i], 0, _tablesizes[i]);
        }
    }

public:

    std::map<HashType, BoundedCounterType> _bigcounts;

    CountingHash( uint64_t single_tablesize,
                  HashFunctorType hash_function) :
        khmer::BaseFilter<HashFunctorType, HashType>(single_tablesize, hash_function), 
        _use_bigcount(false),
        _bigcount_spin_lock(false),
        _max_count( MAX_KCOUNT ),
        _max_bigcount( MAX_BIGCOUNT )
    {
        _allocate_counters();
    }

    CountingHash(std::vector<uint64_t>& tablesizes,
                 HashFunctorType hash_function) :
        khmer::BaseFilter<HashFunctorType, HashType>(tablesizes, hash_function),
        _use_bigcount(false),
        _bigcount_spin_lock(false),
        _max_count( MAX_KCOUNT ),
        _max_bigcount( MAX_BIGCOUNT )
    {
        _allocate_counters();
    }

    virtual ~CountingHash()
    {
        if (_counts) {
            for (size_t i = 0; i < _n_tables; i++) {
                if (_counts[i]) {
                    delete[] _counts[i];
                    _counts[i] = NULL;
                }
            }

            delete[] _counts;
            _counts = NULL;

            _n_tables = 0;
        }
    }

    virtual BoundedCounterType test_and_set_bits(const char * str)
    {
        BoundedCounterType x = get_count(str); // @CTB just hash it, yo.
        count(str);
        return !x;
    }

    virtual BoundedCounterType test_and_set_bits(HashType elem)
    {
        BoundedCounterType x = get_count(elem);
        count(elem);
        return !x;
    }

    virtual void count(const char * str)
    {
        count(_hash(str));
    }

    virtual void count(HashType elem)
    {
        bool is_new = false;
        unsigned int  n_full	  = 0;

        for (unsigned int i = 0; i < _n_tables; i++) {
            const uint64_t bin = elem % _tablesizes[i];
            Byte current_count = _counts[ i ][ bin ];
            if (!is_new) {
                if (current_count == 0) {
                    is_new = true;
                    if (i == 0) {
                        __sync_add_and_fetch(&_occupied_bins, 1);
                    }
                }
            }
            // NOTE: Technically, multiple threads can cause the bin to spill
            //	 over max_count a little, if they all read it as less than
            //	 max_count before any of them increment it.
            //	 However, do we actually care if there is a little
            //	 bit of slop here? It can always be trimmed off later, if
            //	 that would help with stats.
            if ( _max_count > current_count ) {
                __sync_add_and_fetch( *(_counts + i) + bin, 1 );
            } else {
                n_full++;
            }
        } // for each table

        if (n_full == _n_tables && _use_bigcount) {
            while (!__sync_bool_compare_and_swap( &_bigcount_spin_lock, 0, 1 ));
            if (_bigcounts[elem] == 0) {
                _bigcounts[elem] = _max_count + 1;
            } else {
                if (_bigcounts[elem] < _max_bigcount) {
                    _bigcounts[elem] += 1;
                }
            }
            __sync_bool_compare_and_swap( &_bigcount_spin_lock, 1, 0 );
        }

        if (is_new) {
            __sync_add_and_fetch(&_n_unique, 1);
        }

    } // count

    // Alias for the `count function`
    virtual void add(const char * str)
    {
        count(str);
    }

    // Alias for the `count function`
    virtual void add(HashType elem)
    {
        count(elem);
    }

    // get the count for the given k-mer.
    virtual const BoundedCounterType get_count(const char * str) const
    {
        return get_count(_hash(str));
    }

    // get the count for the given k-mer hash.
    virtual const BoundedCounterType get_count(HashType elem) const
    {
        unsigned int	  max_count	= _max_count;
        BoundedCounterType  min_count	= max_count;
        for (unsigned int i = 0; i < _n_tables; i++) {
            BoundedCounterType the_count = _counts[i][elem % _tablesizes[i]];
            if (the_count < min_count) {
                min_count = the_count;
            }
        }
        if (min_count == max_count && _use_bigcount) {
            KmerCountMap::const_iterator it = _bigcounts.find(elem);
            if (it != _bigcounts.end()) {
                min_count = it->second;
            }
        }
        return min_count;
    }


    void set_use_bigcount(bool b)
    {
        _use_bigcount = b;
    }

    bool get_use_bigcount()
    {
        return _use_bigcount;
    }

};

}
