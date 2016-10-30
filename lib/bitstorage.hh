#ifndef BITSTORAGE_HH
#define BITSTORAGE_HH

namespace khmer
{

typedef std::map<HashIntoType, BoundedCounterType> KmerCountMap;
class Hashbits;

class Storage {
public:
    bool _use_bigcount;

    virtual ~Storage();
    virtual std::vector<uint64_t> get_tablesizes() const = 0;
    virtual const size_t n_tables() const = 0;
    virtual void save(std::string) = 0;
    virtual void load(std::string) = 0;
    virtual const uint64_t n_occupied() const = 0;
    virtual const uint64_t n_unique_kmers() const = 0;
    virtual BoundedCounterType     test_and_set_bits( HashIntoType khash ) = 0;
    virtual void count(HashIntoType khash) = 0;
    virtual const BoundedCounterType get_count(HashIntoType khash) const = 0;
    virtual Byte ** get_raw_tables() = 0;

    void set_use_bigcount(bool b)
    {
        _use_bigcount = b;
    }
    bool get_use_bigcount()
    {
        return _use_bigcount;
    }

};


class BitStorage : Storage
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
}

#endif // BITSTORAGE_HH
