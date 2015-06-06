//
// This file is part of khmer, https://github.com/dib-lab/khmer/, and is
// Copyright (C) Michigan State University, 2009-2015. It is licensed under
// the three-clause BSD license; see LICENSE.
// Contact: khmer-project@idyll.org
//

#ifndef COUNTING_HH
#define COUNTING_HH

#include "hashtable.hh"
#include <vector>

namespace khmer
{
typedef std::map<HashIntoType, BoundedCounterType> KmerCountMap;

class CountingHashIntersect;
class CountingHashFile;
class CountingHashFileReader;
class CountingHashFileWriter;
class CountingHashGzFileReader;
class CountingHashGzFileWriter;

class CountingHash : public khmer::Hashtable
{
    friend class CountingHashIntersect;
    friend class CountingHashFile;
    friend class CountingHashFileReader;
    friend class CountingHashFileWriter;
    friend class CountingHashGzFileReader;
    friend class CountingHashGzFileWriter;

protected:
    bool _use_bigcount;		// keep track of counts > Bloom filter hash count threshold?
    uint32_t _bigcount_spin_lock;
    std::vector<HashIntoType> _tablesizes;
    size_t _n_tables;
    HashIntoType _n_unique_kmers;

    Byte ** _counts;

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
    KmerCountMap _bigcounts;

    CountingHash( WordLength ksize, HashIntoType single_tablesize ) :
        khmer::Hashtable(ksize), _use_bigcount(false),
        _bigcount_spin_lock(false), _n_unique_kmers(0)
    {
        _tablesizes.push_back(single_tablesize);

        _allocate_counters();
    }

    CountingHash( WordLength ksize, std::vector<HashIntoType>& tablesizes ) :
        khmer::Hashtable(ksize), _use_bigcount(false),
        _bigcount_spin_lock(false), _tablesizes(tablesizes), _n_unique_kmers(0)
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

    // Writing to the tables outside of defined methods has undefined behavior!
    // As such, this should only be used to return read-only interfaces
    Byte ** get_raw_tables()
    {
        return _counts;
    }

    virtual BoundedCounterType test_and_set_bits(const char * kmer)
    {
        BoundedCounterType x = get_count(kmer); // @CTB just hash it, yo.
        count(kmer);
        return !x;
    }

    virtual BoundedCounterType test_and_set_bits(HashIntoType khash)
    {
        BoundedCounterType x = get_count(khash);
        count(khash);
        return !x;
    }

    std::vector<HashIntoType> get_tablesizes() const
    {
        return _tablesizes;
    }

    virtual const HashIntoType n_unique_kmers() const;

    void set_use_bigcount(bool b)
    {
        _use_bigcount = b;
    }
    bool get_use_bigcount()
    {
        return _use_bigcount;
    }

    virtual void save(std::string);
    virtual void load(std::string);

    // accessors to get table info
    const HashIntoType n_entries() const
    {
        return _tablesizes[0];
    }

    const size_t n_tables() const
    {
        return _n_tables;
    }

    // count number of occupied bins
    virtual const HashIntoType n_occupied(HashIntoType start=0,
                                          HashIntoType stop=0) const
    {
        HashIntoType n = 0;
        if (stop == 0) {
            stop = _tablesizes[0];
        }
        for (HashIntoType i = start; i < stop; i++) {
            if (_counts[0][i % _tablesizes[0]]) {
                n++;
            }
        }
        return n;
    }

    virtual void count(const char * kmer)
    {
        HashIntoType hash = _hash(kmer, _ksize);
        count(hash);
    }

    virtual void count(HashIntoType khash)
    {
        bool is_new_kmer = true;
        unsigned int  n_full	  = 0;

        for (unsigned int i = 0; i < _n_tables; i++) {
            const HashIntoType bin = khash % _tablesizes[i];
            Byte current_count = _counts[ i ][ bin ];
            if (is_new_kmer && current_count != 0) {
                is_new_kmer = false;
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

    // get the count for the given k-mer.
    virtual const BoundedCounterType get_count(const char * kmer) const
    {
        HashIntoType hash = _hash(kmer, _ksize);
        return get_count(hash);
    }

    // get the count for the given k-mer hash.
    virtual const BoundedCounterType get_count(HashIntoType khash) const
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

    void output_fasta_kmer_pos_freq(const std::string &inputfile,
                                    const std::string &outputfile);

    BoundedCounterType get_min_count(const std::string &s);

    BoundedCounterType get_max_count(const std::string &s);

    void get_kadian_count(const std::string &s,
                          BoundedCounterType &kadian,
                          unsigned int nk = 1);

    HashIntoType * abundance_distribution(read_parsers::IParser * parser,
                                          Hashbits * tracking);
    HashIntoType * abundance_distribution(std::string filename,
                                          Hashbits * tracking);

    HashIntoType * fasta_count_kmers_by_position(const std::string &inputfile,
            const unsigned int max_read_len,
            BoundedCounterType limit_by_count=0,
            CallbackFn callback = NULL,
            void * callback_data = NULL);

    void fasta_dump_kmers_by_abundance(const std::string &inputfile,
                                       BoundedCounterType limit_by_count,
                                       CallbackFn callback = NULL,
                                       void * callback_data = NULL);

    unsigned long trim_on_abundance(std::string seq,
                                    BoundedCounterType min_abund) const;
    unsigned long trim_below_abundance(std::string seq,
                                       BoundedCounterType max_abund) const;
    std::vector<unsigned int> find_spectral_error_positions(std::string seq,
            BoundedCounterType min_abund) const;

    void collect_high_abundance_kmers(const std::string &infilename,
                                      unsigned int lower_count,
                                      unsigned int upper_count,
                                      SeenSet& kmers);
};


class CountingHashFile
{
public:
    static void load(const std::string &infilename, CountingHash &ht);
    static void save(const std::string &outfilename, const CountingHash &ht);
};

class CountingHashFileReader : public CountingHashFile
{
public:
    CountingHashFileReader(const std::string &infilename, CountingHash &ht);
};

class CountingHashGzFileReader : public CountingHashFile
{
public:
    CountingHashGzFileReader(const std::string &infilename, CountingHash &ht);
};


class CountingHashFileWriter : public CountingHashFile
{
public:
    CountingHashFileWriter(const std::string &outfilename, const CountingHash &ht);
};

class CountingHashGzFileWriter : public CountingHashFile
{
public:
    CountingHashGzFileWriter(const std::string &outfilename,
                             const CountingHash &ht);
};
};

#endif // COUNTING_HH

// vim: set sts=2 sw=2:
