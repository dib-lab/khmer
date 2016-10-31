#ifndef BITSTORAGE_HH
#define BITSTORAGE_HH

namespace khmer
{

typedef std::map<HashIntoType, BoundedCounterType> KmerCountMap;
class Hashbits;

class Storage {
public:
    bool _use_bigcount;

    Storage() : _use_bigcount(false) { } ;
    virtual ~Storage() { }
    virtual std::vector<uint64_t> get_tablesizes() const = 0;
    virtual const size_t n_tables() const = 0;
    virtual void save(std::string) = 0;
    virtual void load(std::string) = 0;
    virtual const uint64_t n_occupied() const = 0;
    virtual const uint64_t n_unique_kmers() const = 0;
    virtual BoundedCounterType test_and_set_bits( HashIntoType khash ) = 0;
    virtual void count(HashIntoType khash) = 0;
    virtual const BoundedCounterType get_count(HashIntoType khash) const = 0;
    virtual Byte ** get_raw_tables() = 0;

    void set_use_bigcount(bool b);
    bool get_use_bigcount();
};


class BitStorage : public Storage
{
friend class Hashbits;
protected:
    std::vector<uint64_t> _tablesizes;
    size_t _n_tables;
    uint64_t _occupied_bins;
    uint64_t _n_unique_kmers;
    Byte ** _counts;

    BitStorage(std::vector<uint64_t>& tablesizes) :
        _tablesizes(tablesizes) 
    {
        _occupied_bins = 0;
        _n_unique_kmers = 0;

        _allocate_counters();
    }
    ~BitStorage() {
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
    
    void _allocate_counters()
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

    // Accessors for protected/private table info members
    std::vector<uint64_t> get_tablesizes() const
    {
        return _tablesizes;
    }

    const size_t n_tables() const
    {
        return _n_tables;
    }

    void save(std::string) { ; }
    void load(std::string) { ; }

    // count number of occupied bins
    const uint64_t n_occupied() const
    {
        return _occupied_bins;
    }

    const uint64_t n_unique_kmers() const
    {
        return _n_unique_kmers;
    }

    // Get and set the hashbits for the given kmer hash.
    // Generally, it is better to keep tests and mutations separate,
    // but, in the interests of efficiency and thread safety,
    // tests and mutations are being blended here against conventional
    // software engineering wisdom.
    inline
    BoundedCounterType
    test_and_set_bits( HashIntoType khash )
    {
        bool is_new_kmer = false;

        for (size_t i = 0; i < _n_tables; i++) {
            uint64_t bin = khash % _tablesizes[i];
            uint64_t byte = bin / 8;
            unsigned char bit = (unsigned char)(1 << (bin % 8));

            unsigned char bits_orig = __sync_fetch_and_or( *(_counts + i) +
                                      byte, bit );
            if (!(bits_orig & bit)) {
                if (i == 0) {
                    __sync_add_and_fetch( &_occupied_bins, 1 );
                }
                is_new_kmer = true;
            }
        } // iteration over hashtables

        if (is_new_kmer) {
            __sync_add_and_fetch( &_n_unique_kmers, 1 );
            return 1; // kmer not seen before
        }

        return 0; // kmer already seen
    } // test_and_set_bits

    void count(HashIntoType khash)
    {
        test_and_set_bits(khash);
    }

    // get the count for the given k-mer hash.
    const BoundedCounterType get_count(HashIntoType khash) const
    {
        for (size_t i = 0; i < _n_tables; i++) {
            uint64_t bin = khash % _tablesizes[i];
            uint64_t byte = bin / 8;
            unsigned char bit = bin % 8;

            if (!(_counts[i][byte] & (1 << bit))) {
                return 0;
            }
        }
        return 1;
    }

    // Writing to the tables outside of defined methods has undefined behavior!
    // As such, this should only be used to return read-only interfaces
    Byte ** get_raw_tables()
    {
        return _counts;
    }

    // void update_from(const Hashbits &other) { ; }
};


class CountingHash;

class ByteStorage : public Storage {
friend class CountingHash;
protected:
    unsigned int    _max_count;
    unsigned int    _max_bigcount;

    uint32_t _bigcount_spin_lock;
    std::vector<uint64_t> _tablesizes;
    size_t _n_tables;
    uint64_t _n_unique_kmers;
    uint64_t _occupied_bins;

    Byte ** _counts;

    void _allocate_counters()
    {
        _n_tables = _tablesizes.size();

        _counts = new Byte*[_n_tables];
        for (size_t i = 0; i < _n_tables; i++) {
            _counts[i] = new Byte[_tablesizes[i]];
            memset(_counts[i], 0, _tablesizes[i]);
        }
    }
public:
    KmerCountMap _bigcounts;

    ByteStorage(uint64_t single_tablesize ) :
        _max_count(MAX_KCOUNT),
        _bigcount_spin_lock(false), _n_unique_kmers(0), _occupied_bins(0)
    {
        _tablesizes.push_back(single_tablesize);

        _allocate_counters();
    }

    ByteStorage(std::vector<uint64_t>& tablesizes ) :
        _max_count(MAX_KCOUNT),
        _bigcount_spin_lock(false), _tablesizes(tablesizes),
        _n_unique_kmers(0), _occupied_bins(0)
    {

        _allocate_counters();
    }

    ~ByteStorage()
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

    // Writing to the tables outside of defined methods has undefined behavior!
    // As such, this should only be used to return read-only interfaces
    Byte ** get_raw_tables()
    {
        return _counts;
    }

    BoundedCounterType test_and_set_bits(HashIntoType khash)
    {
        BoundedCounterType x = get_count(khash);
        count(khash);
        return !x;
    }

    std::vector<uint64_t> get_tablesizes() const
    {
        return _tablesizes;
    }

    const uint64_t n_unique_kmers() const
    {
        return _n_unique_kmers;
    }

    void save(std::string) { ; }
    void load(std::string) { ; }

    const size_t n_tables() const
    {
        return _n_tables;
    }

    // count number of occupied bins
    const uint64_t n_occupied() const
    {
        return _occupied_bins;
    }

    void count(HashIntoType khash)
    {
        bool is_new_kmer = false;
        unsigned int  n_full	  = 0;

        for (unsigned int i = 0; i < _n_tables; i++) {
            const uint64_t bin = khash % _tablesizes[i];
            Byte current_count = _counts[ i ][ bin ];
            if (!is_new_kmer) {
                if (current_count == 0) {
                    is_new_kmer = true;
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
            if (_bigcounts[khash] == 0) {
                _bigcounts[khash] = _max_count + 1;
            } else {
                if (_bigcounts[khash] < _max_bigcount) {
                    _bigcounts[khash] += 1;
                }
            }
            __sync_bool_compare_and_swap( &_bigcount_spin_lock, 1, 0 );
        }

        if (is_new_kmer) {
            __sync_add_and_fetch(&_n_unique_kmers, 1);
        }

    } // count

    // get the count for the given k-mer hash.
    const BoundedCounterType get_count(HashIntoType khash) const
    {
        unsigned int	  max_count	= _max_count;
        BoundedCounterType  min_count	= max_count;
        for (unsigned int i = 0; i < _n_tables; i++) {
            BoundedCounterType the_count = _counts[i][khash % _tablesizes[i]];
            if (the_count < min_count) {
                min_count = the_count;
            }
        }
        if (min_count == max_count && _use_bigcount) {
            KmerCountMap::const_iterator it = _bigcounts.find(khash);
            if (it != _bigcounts.end()) {
                min_count = it->second;
            }
        }
        return min_count;
    }
};
}

#endif // BITSTORAGE_HH
