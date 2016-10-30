#ifndef BYTESTORAGE_HH
#define BYTESTORAGE_HH

namespace khmer
{

class CountingHash;

class ByteStorage : public Storage {
friend class CountingHash;
protected:
    unsigned int    _max_count;
    unsigned int    _max_bigcount;

    bool _use_bigcount;		// keep track of counts > Bloom filter hash count threshold?
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
        _max_count(MAX_KCOUNT), _use_bigcount(false),
        _bigcount_spin_lock(false), _n_unique_kmers(0), _occupied_bins(0)
    {
        _tablesizes.push_back(single_tablesize);

        _allocate_counters();
    }

    ByteStorage(std::vector<uint64_t>& tablesizes ) :
        _max_count(MAX_KCOUNT), _use_bigcount(false),
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
        std::cout << "XXX " << khash << std::endl;
        bool is_new_kmer = false;
        unsigned int  n_full	  = 0;

        for (unsigned int i = 0; i < _n_tables; i++) {
            const uint64_t bin = khash % _tablesizes[i];
            std::cout << "XXX2 " << bin << std::endl;
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
                std::cout << "XXX3 " << (int) _counts[i][bin] << "\n";
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
        std::cout << "YYY " << khash << "\n";
        unsigned int	  max_count	= _max_count;
        BoundedCounterType  min_count	= max_count;
        for (unsigned int i = 0; i < _n_tables; i++) {
            std::cout << "XXX2 " << khash % _tablesizes[i] << std::endl;
            BoundedCounterType the_count = _counts[i][khash % _tablesizes[i]];
            std::cout << "XXX3 " << (int) the_count << std::endl;
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
#endif // BYTESTORAGE_HH
