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
#ifndef HASHGRAPH_HH
#define HASHGRAPH_HH

#include <stddef.h>
#include <stdint.h>
#include <string.h>
#include <fstream>
#include <iostream>
#include <list>
#include <queue>
#include <set>
#include <string>
#include <vector>
#include <memory>

#include "hashtable.hh"
#include "traversal.hh"
#include "subset.hh"

namespace oxli
{
    namespace read_parsers
    {
        template<typename SeqIO> class ReadParser;
        class FastxReader;
    }
}

#define MAX_KEEPER_SIZE int(1e6)

#define next_f(kmer_f, ch) ((((kmer_f) << 2) & bitmask) | (twobit_repr(ch)))
#define next_r(kmer_r, ch) (((kmer_r) >> 2) | (twobit_comp(ch) << rc_left_shift))

#define prev_f(kmer_f, ch) ((kmer_f) >> 2 | twobit_repr(ch) << rc_left_shift)
#define prev_r(kmer_r, ch) ((((kmer_r) << 2) & bitmask) | (twobit_comp(ch)))

#define set_contains(s, e) ((s).find(e) != (s).end())

#define CALLBACK_PERIOD 100000

namespace oxli
{
//
// Hashgraph: Extension of Hashtable to support graph operations.
//

class Hashgraph: public Hashtable
{

    friend class SubsetPartition;
    friend class LabelHash;
    friend class Traverser;

protected:
    unsigned int _tag_density;

    explicit Hashgraph(WordLength ksize, Storage * s)
        : Hashtable(ksize, s)
    {
        _tag_density = DEFAULT_TAG_DENSITY;
        if (!(_tag_density % 2 == 0)) {
            throw oxli_exception();
        }
        partition = make_shared<SubsetPartition>(this);
        _all_tags_spin_lock = 0;
    }


    // empty the partition structure
    void _clear_all_partitions()
    {
        if (partition != NULL) {
            partition->_clear_all_partitions();
        }
    }

    uint32_t _all_tags_spin_lock;
public:
    // default master partitioning
    shared_ptr<SubsetPartition> partition;

    // tags for sparse graph implementation
    SeenSet all_tags;

    // tags at which to stop traversal
    SeenSet stop_tags;

    // tags used in repartitioning
    SeenSet repart_small_tags;

    // set the minimum density of tagging.
    void _set_tag_density(unsigned int d)
    {
        // must be odd; can't be set if tags exist.
        if (!(d % 2 == 0) || !all_tags.empty()) {
            throw oxli_exception();
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

    bool has_tag(HashIntoType tag) const
    {
        return set_contains(all_tags, tag);
    }

    bool has_stop_tag(HashIntoType stop_tag) const
    {
        return set_contains(stop_tags, stop_tag);
    }

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

    // Consume reads & build sparse graph.
    template<typename SeqIO>
    void consume_seqfile_and_tag(
        std::string const &filename,
        unsigned int &total_reads,
        unsigned long long &n_consumed
    );

    // Count every k-mer from a stream of FASTA or FASTQ reads,
    // using the supplied parser.
    // Tag certain ones on the connectivity graph.
    template<typename SeqIO>
    void consume_seqfile_and_tag(
        read_parsers::ReadParserPtr<SeqIO>& parser,
        unsigned int &total_reads,
        unsigned long long &n_consumed
    );

    // consume a string & add sparse graph nodes.
    void consume_sequence_and_tag(const std::string& seq,
                                  unsigned long long& n_consumed,
                                  SeenSet * new_tags = 0);

    // get the tags present in this sequence.
    void get_tags_for_sequence(const std::string& seq,
                               SeenSet& tags) const;

    // consume an already-partitioned file & load in the partition IDs
    template<typename SeqIO>
    void consume_partitioned_fasta(const std::string &filename,
                                   unsigned int &total_reads,
                                   unsigned long long &n_consumed);

    // trim the given sequence on stoptags
    size_t trim_on_stoptags(std::string sequence) const;

    // @@
    unsigned int traverse_from_kmer(Kmer start,
                                    unsigned int radius,
                                    KmerSet &keeper,
                                    unsigned int max_count = MAX_KEEPER_SIZE)
    const;

    // print, save, and load the set of tags.
    void print_tagset(std::string);
    void save_tagset(std::string);
    void load_tagset(std::string, bool clear_tags=true);

    // print, save and load the set of stop tags.
    void print_stop_tags(std::string);
    void save_stop_tags(std::string);
    void load_stop_tags(std::string filename, bool clear_tags=true);

    // @@
    void extract_unique_paths(std::string seq,
                              unsigned int min_length,
                              float min_unique_f,
                              std::vector<std::string> &results);

    // @@
    void calc_connected_graph_size(Kmer node,
                                   unsigned long long& count,
                                   KmerSet& keeper,
                                   const unsigned long long threshold=0,
                                   bool break_on_circum=false) const;

    // Calculate the graph degree of the given k-mer.
    unsigned int kmer_degree(HashIntoType kmer_f, HashIntoType kmer_r);
    unsigned int kmer_degree(const char * kmer_s);

    // Find all nodes with a degree > 2.
    void find_high_degree_nodes(const char * sequence,
                                SeenSet& high_degree_nodes) const;

    // Find the maximal linear path (nodes degree <= 2)
    unsigned int traverse_linear_path(const Kmer start_kmer,
                                      SeenSet &adjacencies,
                                      SeenSet &nodes, Hashtable& bf,
                                      SeenSet &high_degree_nodes) const;

    //
    // for debugging/testing purposes only!
    //

    // check partition map validity.
    void _validate_pmap()
    {
        if (partition) {
            partition->_validate_pmap();
        }
    }

};

// Hashgraph-derived class with ByteStorage.
class Countgraph : public oxli::Hashgraph
{
public:
    explicit Countgraph(WordLength ksize, std::vector<uint64_t> sizes)
        : Hashgraph(ksize, new ByteStorage(sizes)) { } ;
};

// Hashgraph-derived class with NibbleStorage.
class SmallCountgraph : public oxli::Hashgraph
{
public:
    explicit SmallCountgraph(WordLength ksize, std::vector<uint64_t> sizes)
        : Hashgraph(ksize, new NibbleStorage(sizes)) { } ;
};

// Hashgraph-derived class with BitStorage.
class Nodegraph : public Hashgraph
{
public:
    explicit Nodegraph(WordLength ksize, std::vector<uint64_t> sizes)
        : Hashgraph(ksize, new BitStorage(sizes)) { } ;

    void update_from(const Nodegraph &other);
    double similarity(const Nodegraph &other);
    double containment(const Nodegraph &other);
};

}

#define ACQUIRE_ALL_TAGS_SPIN_LOCK \
  while (!__sync_bool_compare_and_swap( &_all_tags_spin_lock, 0, 1 ));

#define RELEASE_ALL_TAGS_SPIN_LOCK \
  __sync_bool_compare_and_swap( &_all_tags_spin_lock, 1, 0 );

#endif // HASHGRAPH_HH
