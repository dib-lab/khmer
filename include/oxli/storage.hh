/*
This file is part of khmer, https://github.com/dib-lab/khmer/, and is
Copyright (C) 2016, The Regents of the University of California.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are
met:

    * Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.

    * Redistributions in binary form must reproduce the above
      copyright notice, this list of conditions and the following
      disclaimer in the documentation and/or other materials provided
      with the distribution.

    * Neither the name of the University of California nor the names
      of its contributors may be used to endorse or promote products
      derived from this software without specific prior written
      permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
"AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
LICENSE (END)

Contact: khmer-project@idyll.org
*/

#ifndef STORAGE_HH
#define STORAGE_HH

#include <cassert>
#include <array>
#include <mutex>
#include <unordered_map>
using MuxGuard = std::lock_guard<std::mutex>;

#include "gqf.h"

namespace oxli {
typedef std::unordered_map<HashIntoType, BoundedCounterType> KmerCountMap;

//
// base Storage class for hashtable-related storage of information in memory.
//

class Storage
{
protected:
    bool _supports_bigcount;
    bool _use_bigcount;

public:
    Storage() : _supports_bigcount(false), _use_bigcount(false) { } ;
    virtual ~Storage() { }
    virtual std::vector<uint64_t> get_tablesizes() const = 0;
    virtual const size_t n_tables() const = 0;
    virtual void save(std::string, WordLength) = 0;
    virtual void load(std::string, WordLength&) = 0;
    virtual const uint64_t n_occupied() const = 0;
    virtual const uint64_t n_unique_kmers() const = 0;
    virtual BoundedCounterType test_and_set_bits( HashIntoType khash ) = 0;
    virtual bool add(HashIntoType khash) = 0;
    virtual const BoundedCounterType get_count(HashIntoType khash) const = 0;
    virtual Byte ** get_raw_tables() = 0;

    void set_use_bigcount(bool b);
    bool get_use_bigcount();
};


/*
 * \class BitStorage
 *
 * \brief A Bloom filter implementation.
 *
 * BitStorage is used to track presence/absence of k-mers by Hashtable
 * and derived classes.  It contains 'n_tables' different tables of
 * bitsizes specified in 'tablesizes' (so 1/8 for bytesizes).
 *
 * Like other Storage classes, BitStorage manages setting the bits and
 * tracking statistics, as well as save/load, and not much else.
 *
 */

class BitStorage : public Storage
{
protected:
    std::vector<uint64_t> _tablesizes;
    size_t _n_tables;
    uint64_t _occupied_bins;
    uint64_t _n_unique_kmers;
    Byte ** _counts;

public:
    BitStorage(std::vector<uint64_t>& tablesizes) :
        _tablesizes(tablesizes)
    {
        _occupied_bins = 0;
        _n_unique_kmers = 0;

        _allocate_counters();
    }
    ~BitStorage()
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

    void save(std::string, WordLength ksize);
    void load(std::string, WordLength& ksize);

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

    inline bool add(HashIntoType khash)
    {
        return test_and_set_bits(khash);
    }

    // get the count for the given k-mer hash.
    inline const BoundedCounterType get_count(HashIntoType khash) const
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

    void update_from(const BitStorage&);
    double similarity(const BitStorage&);
    double containment(const BitStorage&);
};


/*
 * \class NibbleStorage
 *
 * \brief A A CountMin sketch implementation using 4bit counters.
 *
 * NibbleStorage is used to track counts of k-mers by Hashtable
 * and derived classes.  It contains 'n_tables' different tables of
 * 'tablesizes' entries. It allocates half a byte per table entry.
 *
 * Like other Storage classes, NibbleStorage manages setting the bits and
 * tracking statistics, as well as save/load, and not much else.
 *
 */
class NibbleStorage : public Storage
{
protected:
    // table size is measured in number of entries in the table, not in bytes
    std::vector<uint64_t> _tablesizes;
    size_t _n_tables;
    uint64_t _occupied_bins;
    uint64_t _n_unique_kmers;
    std::array<std::mutex, 32> mutexes;
    static constexpr uint8_t _max_count{15};
    Byte ** _counts;

    // Compute index into the table, this retrieves the correct byte
    // which you then need to select the correct nibble from
    uint64_t _table_index(const HashIntoType k, const uint64_t tablesize) const
    {
        return (k % tablesize) / 2;
    }
    // Compute which half of the byte to use for this hash value
    uint8_t _mask(const HashIntoType k, const uint64_t tablesize) const
    {
        return (k%tablesize)%2 ? 15 : 240;
    }
    // Compute which half of the byte to use for this hash value
    uint8_t _shift(const HashIntoType k, const uint64_t tablesize) const
    {
        return (k%tablesize)%2 ? 0 : 4;
    }

public:
    NibbleStorage(std::vector<uint64_t>& tablesizes) :
        _tablesizes{tablesizes},
        _occupied_bins{0}, _n_unique_kmers{0}
    {
        // to allow more than 32 tables increase the size of mutex pool
        assert(_n_tables <= 32);
        _allocate_counters();
    }

    ~NibbleStorage()
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

    void _allocate_counters()
    {
        _n_tables = _tablesizes.size();

        _counts = new Byte*[_n_tables];

        for (size_t i = 0; i < _n_tables; i++) {
            const uint64_t tablesize = _tablesizes[i];
            const uint64_t tablebytes = tablesize / 2 + 1;

            _counts[i] = new Byte[tablebytes];
            memset(_counts[i], 0, tablebytes);
        }
    }


    BoundedCounterType test_and_set_bits(HashIntoType khash)
    {
        BoundedCounterType x = get_count(khash);
        add(khash);
        return !x;
    }

    bool add(HashIntoType khash)
    {
        bool is_new_kmer = false;

        for (unsigned int i = 0; i < _n_tables; i++) {
            MuxGuard g(mutexes[i]);
            Byte* const table(_counts[i]);
            const uint64_t idx = _table_index(khash, _tablesizes[i]);
            const uint8_t mask = _mask(khash, _tablesizes[i]);
            const uint8_t shift = _shift(khash, _tablesizes[i]);
            const uint8_t current_count = (table[idx] & mask) >> shift;

            if (!is_new_kmer) {
                if (current_count == 0) {
                    is_new_kmer = true;

                    // track occupied bins in the first table only, as proxy
                    // for all.
                    if (i == 0) {
                        __sync_add_and_fetch(&_occupied_bins, 1);
                    }
                }
            }
            // if we have reached the maximum count stop incrementing the
            // counter. This avoids overflowing it.
            if (current_count == _max_count) {
                continue;
            }

            // increase count, no checking for overflow
            const uint8_t new_count = (current_count + 1) << shift;
            table[idx] = (table[idx] & ~mask) | (new_count & mask);
        }

        if (is_new_kmer) {
            __sync_add_and_fetch(&_n_unique_kmers, 1);
        }

        return is_new_kmer;
    }

    // get the count for the given k-mer hash.
    const BoundedCounterType get_count(HashIntoType khash) const
    {
        uint8_t min_count = _max_count; // bound count by maximum

        // get the minimum count across all tables
        for (unsigned int i = 0; i < _n_tables; i++) {
            const Byte* table(_counts[i]);
            const uint64_t idx = _table_index(khash, _tablesizes[i]);
            const uint8_t mask = _mask(khash, _tablesizes[i]);
            const uint8_t shift = _shift(khash, _tablesizes[i]);
            const uint8_t the_count = (table[idx] & mask) >> shift;

            if (the_count < min_count) {
                min_count = the_count;
            }
        }
        return min_count;
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
    const uint64_t n_unique_kmers() const
    {
        return _n_unique_kmers;
    }
    const uint64_t n_occupied() const
    {
        return _occupied_bins;
    }
    void save(std::string outfilename, WordLength ksize);
    void load(std::string infilename, WordLength& ksize);

    Byte ** get_raw_tables()
    {
        return _counts;
    }
};


/*
 * \class QFStorage
 *
 * \brief A Quotient Filter storage
 */
 class QFStorage : public Storage {
protected:
  QF cf;

public:
  QFStorage(int size) {
    // size is the power of two to specify the number of slots in
    // the filter (2**size). Third argument sets the number of bits used
    // in the key (current value of size+8 is copied from the CQF example)
    // Final argument is the number of bits allocated for the value, which
    // we do not use.
    qf_init(&cf, (1ULL << size), size+8, 0);
  }

  ~QFStorage() { qf_destroy(&cf); }

  BoundedCounterType test_and_set_bits(HashIntoType khash) {
    BoundedCounterType x = get_count(khash);
    add(khash);
    return !x;
  }

  //
  bool add(HashIntoType khash) {
      bool is_new = get_count(khash) == 0;
      qf_insert(&cf, khash % cf.range, 0, 1);
      return is_new;
  }

  // get the count for the given k-mer hash.
  const BoundedCounterType get_count(HashIntoType khash) const {
    return qf_count_key_value(&cf, khash % cf.range, 0);
  }

  // Accessors for protected/private table info members
  // xnslots is larger than nslots. It includes some extra slots to deal
  // with some details of how the counting is implemented
  std::vector<uint64_t> get_tablesizes() const { return {cf.xnslots}; }
  const size_t n_tables() const { return 1; }
  const uint64_t n_unique_kmers() const { return cf.ndistinct_elts; }
  const uint64_t n_occupied() const { return cf.noccupied_slots; }
  void save(std::string outfilename, WordLength ksize);
  void load(std::string infilename, WordLength &ksize);

  Byte **get_raw_tables() { return nullptr; }
};


/*
 * \class ByteStorage
 *
 * \brief A CountMin sketch implementation.
 *
 * ByteStorage is used to track counts of k-mers by Hashtable
 * and derived classes.  It contains 'n_tables' different tables of
 * bytesizes specified in 'tablesizes'.
 *
 * Like other Storage classes, ByteStorage manages setting the bits and
 * tracking statistics, as well as save/load, and not much else.
 *
 */

class ByteStorageFile;
class ByteStorageFileReader;
class ByteStorageFileWriter;
class ByteStorageGzFileReader;
class ByteStorageGzFileWriter;

class ByteStorage : public Storage
{
    friend class ByteStorageFile;
    friend class ByteStorageFileReader;
    friend class ByteStorageFileWriter;
    friend class ByteStorageGzFileReader;
    friend class ByteStorageGzFileWriter;
    friend class CountGraph;
protected:
    unsigned int    _max_count;
    unsigned int    _max_bigcount;

    uint32_t _bigcount_spin_lock;
    std::vector<uint64_t> _tablesizes;
    size_t _n_tables;
    uint64_t _n_unique_kmers;
    uint64_t _occupied_bins;

    Byte ** _counts;

    // initialize counts with empty hashtables.
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

    // constructor: create an empty CountMin sketch.
    ByteStorage(std::vector<uint64_t>& tablesizes ) :
        _max_count(MAX_KCOUNT), _max_bigcount(MAX_BIGCOUNT),
        _bigcount_spin_lock(false), _tablesizes(tablesizes),
        _n_unique_kmers(0), _occupied_bins(0)
    {
        _supports_bigcount = true;
        _allocate_counters();
    }

    // destructor: clear out the memory.
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

    std::vector<uint64_t> get_tablesizes() const
    {
        return _tablesizes;
    }

    const uint64_t n_unique_kmers() const
    {
        return _n_unique_kmers;
    }
    const size_t n_tables() const
    {
        return _n_tables;
    }
    const uint64_t n_occupied() const
    {
        return _occupied_bins;
    }

    void save(std::string, WordLength);
    void load(std::string, WordLength&);

    inline BoundedCounterType test_and_set_bits(HashIntoType khash)
    {
        BoundedCounterType x = get_count(khash);
        add(khash);
        return !x;
    }

    inline bool add(HashIntoType khash)
    {
        bool is_new_kmer = false;
        unsigned int  n_full	  = 0;

        // add one to each entry in each table.
        for (unsigned int i = 0; i < _n_tables; i++) {
            const uint64_t bin = khash % _tablesizes[i];
            Byte current_count = _counts[ i ][ bin ];

            if (!is_new_kmer) {
                if (current_count == 0) {
                    is_new_kmer = true;

                    // track occupied bins in the first table only, as proxy
                    // for all.
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

        // if all tables are full for this position, then add in bigcounts.
        if (n_full == _n_tables && _use_bigcount) {
            while (!__sync_bool_compare_and_swap(&_bigcount_spin_lock, 0, 1));
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

        return is_new_kmer;
    }

    // get the count for the given k-mer hash.
    inline const BoundedCounterType get_count(HashIntoType khash) const
    {
        unsigned int	  max_count	= _max_count;
        BoundedCounterType  min_count	= max_count; // bound count by max.

        // first, get the min count across all tables (standard CMS).
        for (unsigned int i = 0; i < _n_tables; i++) {
            BoundedCounterType the_count = _counts[i][khash % _tablesizes[i]];
            if (the_count < min_count) {
                min_count = the_count;
            }
        }

        // if the count is saturated, check in the bigcount structure to
        // see if we've accumulated more counts.
        if (min_count == max_count && _use_bigcount) {
            KmerCountMap::const_iterator it = _bigcounts.find(khash);
            if (it != _bigcounts.end()) {
                min_count = it->second;
            }
        }
        return min_count;
    }
    // Get direct access to the counts.
    //
    // Note:
    // Writing to the tables outside of defined methods has undefined behavior!
    // As such, this should only be used to return read-only interfaces
    Byte ** get_raw_tables()
    {
        return _counts;
    }

};

// Helper classes for saving ByteStorage objs to disk & loading them.

class ByteStorageFile
{
public:
    static void load(const std::string &infilename,
                     WordLength &ksize,
                     ByteStorage &store);
    static void save(const std::string &outfilename,
                     const WordLength ksize,
                     const ByteStorage &store);
};

class ByteStorageFileReader : public ByteStorageFile
{
public:
    ByteStorageFileReader(const std::string &infilename,
                          WordLength &ksize,
                          ByteStorage &store);
};

class ByteStorageGzFileReader : public ByteStorageFile
{
public:
    ByteStorageGzFileReader(const std::string &infilename,
                            WordLength &ksize,
                            ByteStorage &store);
};


class ByteStorageFileWriter : public ByteStorageFile
{
public:
    ByteStorageFileWriter(const std::string &outfilename,
                          const WordLength ksize,
                          const ByteStorage &store);
};

class ByteStorageGzFileWriter : public ByteStorageFile
{
public:
    ByteStorageGzFileWriter(const std::string &outfilename,
                            const WordLength ksize,
                            const ByteStorage &store);
};
}

#endif // STORAGE_HH
