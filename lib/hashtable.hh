/*
This file is part of khmer, https://github.com/dib-lab/khmer/, and is
Copyright (C) 2010-2015, Michigan State University.
Copyright (C) 2015, The Regents of the University of California.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are
met:

    * Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.

    * Redistributions in binary form must reproduce the above
      copyright notice, this list of conditions and the following
      disclaimer in the documentation and/or other materials provided
      with the distribution.

    * Neither the name of the Michigan State University nor the names
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
#ifndef HASHTABLE_HH
#define HASHTABLE_HH


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
#include "subset.hh"

namespace khmer
{
class CountingHash;
class Hashtable;

namespace read_parsers
{
struct IParser;
}  // namespace read_parsers
}  // namespace khmer

#define MAX_KEEPER_SIZE int(1e6)

#define next_f(kmer_f, ch) ((((kmer_f) << 2) & bitmask) | (twobit_repr(ch)))
#define next_r(kmer_r, ch) (((kmer_r) >> 2) | (twobit_comp(ch) << rc_left_shift))

#define prev_f(kmer_f, ch) ((kmer_f) >> 2 | twobit_repr(ch) << rc_left_shift)
#define prev_r(kmer_r, ch) ((((kmer_r) << 2) & bitmask) | (twobit_comp(ch)))

#define set_contains(s, e) ((s).find(e) != (s).end())

#define CALLBACK_PERIOD 100000

namespace khmer
{

class Hashtable: public
    KmerFactory  		// Base class implementation of a Bloom ht.
{
    friend class SubsetPartition;
    friend class LabelHash;
    friend class Traverser;

protected:
    unsigned int _tag_density;

    unsigned int    _max_count;
    unsigned int    _max_bigcount;

    //WordLength	    _ksize;
    HashIntoType    bitmask;
    unsigned int    _nbits_sub_1;

    explicit Hashtable( WordLength ksize )
        : KmerFactory( ksize ),
          _max_count( MAX_KCOUNT ),
          _max_bigcount( MAX_BIGCOUNT )
    {
        _tag_density = DEFAULT_TAG_DENSITY;
        if (!(_tag_density % 2 == 0)) {
            throw khmer_exception();
        }
        _init_bitstuff();
        partition = new SubsetPartition(this);
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

    explicit Hashtable(const Hashtable&);
    Hashtable& operator=(const Hashtable&);

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

    // number of unique k-mers
    virtual const HashIntoType n_unique_kmers() const = 0;

    // count number of occupied bins
    virtual const HashIntoType n_occupied() const = 0;

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

    void filter_if_present(const std::string &infilename,
                           const std::string &outputfilename);

    size_t trim_on_stoptags(std::string sequence) const;

    void traverse_from_tags(unsigned int distance,
                            unsigned int threshold,
                            unsigned int num_high_todo,
                            CountingHash &counting);

    unsigned int traverse_from_kmer(Kmer start,
                                    unsigned int radius,
                                    KmerSet &keeper,
                                    unsigned int max_count = MAX_KEEPER_SIZE)
    const;

    unsigned int count_and_transfer_to_stoptags(KmerSet &keeper,
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

    void calc_connected_graph_size(Kmer node,
                                   unsigned long long& count,
                                   KmerSet& keeper,
                                   const unsigned long long threshold=0,
                                   bool break_on_circum=false) const;

    typedef void (*kmer_cb)(const char * k, unsigned int n_reads, void *data);


    unsigned int kmer_degree(HashIntoType kmer_f, HashIntoType kmer_r);
    unsigned int kmer_degree(const char * kmer_s);

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
}



#define ACQUIRE_ALL_TAGS_SPIN_LOCK \
  while (!__sync_bool_compare_and_swap( &_all_tags_spin_lock, 0, 1 ));

#define RELEASE_ALL_TAGS_SPIN_LOCK \
  __sync_bool_compare_and_swap( &_all_tags_spin_lock, 1, 0 );

#endif // HASHTABLE_HH
