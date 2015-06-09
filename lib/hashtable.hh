//
// This file is part of khmer, https://github.com/dib-lab/khmer/, and is
// Copyright (C) Michigan State University, 2009-2015. It is licensed under
// the three-clause BSD license; see LICENSE.
// Contact: khmer-project@idyll.org
//

#ifndef HASHTABLE_HH
#define HASHTABLE_HH


#include <vector>
#include <iostream>
#include <list>
#include <queue>

#include <fstream>
#include <string>
#include <set>
#include <map>
#include <queue>

#include "khmer.hh"
#include "khmer_exception.hh"
#include "read_parsers.hh"
#include "subset.hh"
#include "kmer_hash.hh"

#define MAX_KEEPER_SIZE int(1e6)

#define next_f(kmer_f, ch) ((((kmer_f) << 2) & bitmask) | (twobit_repr(ch)))
#define next_r(kmer_r, ch) (((kmer_r) >> 2) | (twobit_comp(ch) << rc_left_shift))

#define prev_f(kmer_f, ch) ((kmer_f) >> 2 | twobit_repr(ch) << rc_left_shift)
#define prev_r(kmer_r, ch) ((((kmer_r) << 2) & bitmask) | (twobit_comp(ch)))

#define set_contains(s, e) ((s).find(e) != (s).end())

#define CALLBACK_PERIOD 100000

namespace khmer
{
#ifdef WITH_INTERNAL_METRICS
struct HashTablePerformanceMetrics : public IPerformanceMetrics {

    enum {
        MKEY_TIME_NORM_READ,
        MKEY_TIME_HASH_KMER,
        MKEY_TIME_UPDATE_TALLIES
    };

    uint64_t	clock_nsecs_norm_read;
    uint64_t	cpu_nsecs_norm_read;
    uint64_t	clock_nsecs_hash_kmer;
    uint64_t	cpu_nsecs_hash_kmer;
    uint64_t	clock_nsecs_update_tallies;
    uint64_t	cpu_nsecs_update_tallies;

    HashTablePerformanceMetrics( );
    virtual ~HashTablePerformanceMetrics( );

    virtual void	accumulate_timer_deltas( uint32_t metrics_key );

};
#endif

//
// Sequence iterator class, test.  Not really a C++ iterator yet.
//

class KMerIterator
{
protected:
    const char * _seq;
    const unsigned char _ksize;

    HashIntoType _kmer_f, _kmer_r;
    HashIntoType bitmask;
    unsigned int _nbits_sub_1;
    unsigned int index;
    size_t length;
    bool initialized;
public:
    KMerIterator(const char * seq, unsigned char k) : _seq(seq), _ksize(k)
    {
        bitmask = 0;
        for (unsigned char i = 0; i < _ksize; i++) {
            bitmask = (bitmask << 2) | 3;
        }
        _nbits_sub_1 = (_ksize*2 - 2);

        index = _ksize - 1;
        length = strlen(seq);
        _kmer_f = 0;
        _kmer_r = 0;

        initialized = false;
    }

    HashIntoType first(HashIntoType& f, HashIntoType& r)
    {
        HashIntoType x;
        x = _hash(_seq, _ksize, _kmer_f, _kmer_r);

        f = _kmer_f;
        r = _kmer_r;

        index = _ksize;

        return x;
    }

    HashIntoType next(HashIntoType& f, HashIntoType& r)
    {
        if (done()) {
            throw khmer_exception();
        }

        if (!initialized) {
            initialized = true;
            return first(f, r);
        }

        unsigned char ch = _seq[index];
        index++;
        if (!(index <= length)) {
            throw khmer_exception();
        }

        // left-shift the previous hash over
        _kmer_f = _kmer_f << 2;

        // 'or' in the current nt
        _kmer_f |= twobit_repr(ch);

        // mask off the 2 bits we shifted over.
        _kmer_f &= bitmask;

        // now handle reverse complement
        _kmer_r = _kmer_r >> 2;
        _kmer_r |= (twobit_comp(ch) << _nbits_sub_1);

        f = _kmer_f;
        r = _kmer_r;

        return uniqify_rc(_kmer_f, _kmer_r);
    }

    HashIntoType first()
    {
        return first(_kmer_f, _kmer_r);
    }
    HashIntoType next()
    {
        return next(_kmer_f, _kmer_r);
    }

    bool done()
    {
        return index >= length;
    }

    unsigned int get_start_pos() const
    {
        return index - _ksize;
    }

    unsigned int get_end_pos() const
    {
        return index;
    }
}; // class KMerIterator

class Hashtable  		// Base class implementation of a Bloom ht.
{
    friend class SubsetPartition;
    friend class LabelHash;
protected:
    unsigned int _tag_density;

    unsigned int    _max_count;
    unsigned int    _max_bigcount;

    WordLength	    _ksize;
    HashIntoType    bitmask;
    unsigned int    _nbits_sub_1;

    Hashtable( WordLength ksize )
        : _max_count( MAX_KCOUNT ),
          _max_bigcount( MAX_BIGCOUNT ),
          _ksize( ksize )
    {
        _tag_density = DEFAULT_TAG_DENSITY;
        if (!(_tag_density % 2 == 0)) {
            throw khmer_exception();
        }
        partition = new SubsetPartition(this);
        _init_bitstuff();
        _all_tags_spin_lock = 0;

    }

    virtual ~Hashtable( )
    {
        delete partition;
    }

    void _init_bitstuff()
    {
        bitmask = 0;
        for (unsigned int i = 0; i < _ksize; i++) {
            bitmask = (bitmask << 2) | 3;
        }
        _nbits_sub_1 = (_ksize*2 - 2);
    }

    HashIntoType _next_hash(char ch, HashIntoType &h, HashIntoType &r) const
    {
        // left-shift the previous hash over
        h = h << 2;

        // 'or' in the current nt
        h |= twobit_repr(ch);

        // mask off the 2 bits we shifted over.
        h &= bitmask;

        // now handle reverse complement
        r = r >> 2;
        r |= (twobit_comp(ch) << _nbits_sub_1);

        return uniqify_rc(h, r);
    }

    void _clear_all_partitions()
    {
        if (partition != NULL) {
            partition->_clear_all_partitions();
        }
    }

    uint32_t _all_tags_spin_lock;

    NONCOPYABLE(Hashtable);

public:
    SubsetPartition * partition;
    SeenSet all_tags;
    SeenSet stop_tags;
    SeenSet repart_small_tags;

    // accessor to get 'k'
    const WordLength ksize() const
    {
        return _ksize;
    }

    virtual void count(const char * kmer) = 0;
    virtual void count(HashIntoType khash) = 0;

    // get the count for the given k-mer.
    virtual const BoundedCounterType get_count(const char * kmer) const = 0;
    virtual const BoundedCounterType get_count(HashIntoType khash) const = 0;

    virtual void save(std::string) = 0;
    virtual void load(std::string) = 0;

    // count every k-mer in the string.
    unsigned int consume_string(const std::string &s);

    // checks each read for non-ACGT characters
    bool check_and_normalize_read(std::string &read) const;

    // check each read for non-ACGT characters, and then consume it.
    unsigned int check_and_process_read(std::string &read,
                                        bool &is_valid);

    // Count every k-mer in a FASTA or FASTQ file.
    // Note: Yes, the name 'consume_fasta' is a bit misleading,
    //	     but the FASTA format is effectively a subset of the FASTQ format
    //	     and the FASTA portion is what we care about in this case.
    void consume_fasta(
        std::string const   &filename,
        unsigned int	    &total_reads,
        unsigned long long  &n_consumed
    );
    // Count every k-mer from a stream of FASTA or FASTQ reads,
    // using the supplied parser.
    void consume_fasta(
        read_parsers:: IParser *	    parser,
        unsigned int	    &total_reads,
        unsigned long long  &n_consumed
    );

    bool median_at_least(const std::string &s,
                         unsigned int cutoff);

    void get_median_count(const std::string &s,
                          BoundedCounterType &median,
                          float &average,
                          float &stddev);
    virtual const HashIntoType n_unique_kmers() const = 0;

    // partitioning stuff
    void _validate_pmap()
    {
        if (partition) {
            partition->_validate_pmap();
        }
    }

    virtual void save_tagset(std::string);
    virtual void load_tagset(std::string, bool clear_tags=true);

    // for debugging/testing purposes only!
    void _set_tag_density(unsigned int d)
    {
        if (!(d % 2 == 0) || !all_tags.empty()) { // must be even and tags must exist
            throw khmer_exception();
        }
        _tag_density = d;
    }

    unsigned int _get_tag_density() const
    {
        return _tag_density;
    }

    void add_tag(HashIntoType tag)
    {
        all_tags.insert(tag);
    }
    void add_stop_tag(HashIntoType tag)
    {
        stop_tags.insert(tag);
    }

    // Partitioning stuff.

    size_t n_tags() const
    {
        return all_tags.size();
    }

    void divide_tags_into_subsets(unsigned int subset_size, SeenSet& divvy);

    void add_kmer_to_tags(HashIntoType kmer)
    {
        all_tags.insert(kmer);
    }

    void clear_tags()
    {
        all_tags.clear();
    }

    // Count every k-mer in a FASTA or FASTQ file.
    // Tag certain ones on the connectivity graph.
    void consume_fasta_and_tag(
        std::string const	  &filename,
        unsigned int	  &total_reads,
        unsigned long long  &n_consumed
    );

    // Count every k-mer from a stream of FASTA or FASTQ reads,
    // using the supplied parser.
    // Tag certain ones on the connectivity graph.
    void consume_fasta_and_tag(
        read_parsers:: IParser *	    parser,
        unsigned int	    &total_reads,
        unsigned long long  &n_consumed
    );

    void consume_sequence_and_tag(const std::string& seq,
                                  unsigned long long& n_consumed,
                                  SeenSet * new_tags = 0);


    void consume_fasta_and_tag_with_stoptags(const std::string &filename,
            unsigned int &total_reads,
            unsigned long long &n_consumed);
    void consume_fasta_and_traverse(const std::string &filename,
                                    unsigned int distance,
                                    unsigned int big_threshold,
                                    unsigned int transfer_threshold,
                                    CountingHash &counting);

    void consume_partitioned_fasta(const std::string &filename,
                                   unsigned int &total_reads,
                                   unsigned long long &n_consumed);

    virtual BoundedCounterType test_and_set_bits(const char * kmer) = 0;
    virtual BoundedCounterType test_and_set_bits(HashIntoType khash) = 0;

    virtual std::vector<HashIntoType> get_tablesizes() const = 0;
    virtual const size_t n_tables() const = 0;
    virtual const HashIntoType n_occupied(HashIntoType start=0,
                                          HashIntoType stop=0) const = 0;
    virtual const HashIntoType n_entries() const = 0;

    void filter_if_present(const std::string &infilename,
                           const std::string &outputfilename);

    unsigned int count_kmers_within_radius(HashIntoType kmer_f,
                                           HashIntoType kmer_r,
                                           unsigned int radius,
                                           unsigned int max_count,
                                           const SeenSet * seen=0) const;

    size_t trim_on_stoptags(std::string sequence) const;

    void traverse_from_tags(unsigned int distance,
                            unsigned int threshold,
                            unsigned int num_high_todo,
                            CountingHash &counting);

    unsigned int traverse_from_kmer(HashIntoType start,
                                    unsigned int radius,
                                    SeenSet &keeper) const;

    unsigned int count_and_transfer_to_stoptags(SeenSet &keeper,
            unsigned int threshold,
            CountingHash &counting);

    virtual void print_tagset(std::string);
    virtual void print_stop_tags(std::string);
    virtual void save_stop_tags(std::string);
    void load_stop_tags(std::string filename, bool clear_tags=true);

    void identify_stop_tags_by_position(std::string sequence,
                                        std::vector<unsigned int> &posns)
    const;

    void extract_unique_paths(std::string seq,
                              unsigned int min_length,
                              float min_unique_f,
                              std::vector<std::string> &results);

    void calc_connected_graph_size(const char * kmer,
                                   unsigned long long& count,
                                   SeenSet& keeper,
                                   const unsigned long long threshold=0,
                                   bool break_on_circum=false) const
    {
        HashIntoType r, f;
        _hash(kmer, _ksize, f, r);
        calc_connected_graph_size(f, r, count, keeper, threshold, break_on_circum);
    }

    void calc_connected_graph_size(const HashIntoType kmer_f,
                                   const HashIntoType kmer_r,
                                   unsigned long long& count,
                                   SeenSet& keeper,
                                   const unsigned long long threshold=0,
                                   bool break_on_circum=false) const;

    typedef void (*kmer_cb)(const char * k, unsigned int n_reads, void *data);


    unsigned int kmer_degree(HashIntoType kmer_f, HashIntoType kmer_r) const;
    unsigned int kmer_degree(const char * kmer_s) const
    {
        HashIntoType kmer_f, kmer_r;
        _hash(kmer_s, _ksize, kmer_f, kmer_r);

        return kmer_degree(kmer_f, kmer_r);
    }

    // return all k-mer substrings, on the forward strand.
    void get_kmers(const std::string &s, std::vector<std::string> &kmers)
    const;

    // return hash values for all k-mer substrings
    void get_kmer_hashes(const std::string &s,
                         std::vector<HashIntoType> &kmers) const;

    // return counts of all k-mers in this string.
    void get_kmer_counts(const std::string &s,
                         std::vector<BoundedCounterType> &counts) const;
};
};



#define ACQUIRE_ALL_TAGS_SPIN_LOCK \
  while (!__sync_bool_compare_and_swap( &_all_tags_spin_lock, 0, 1 ));

#define RELEASE_ALL_TAGS_SPIN_LOCK \
  __sync_bool_compare_and_swap( &_all_tags_spin_lock, 1, 0 );

#endif // HASHTABLE_HH
