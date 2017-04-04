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
#ifndef SUBSET_HH
#define SUBSET_HH

#include <stddef.h>
#include <queue>
#include <string>

#include "oxli.hh"

namespace oxli
{
class Countgraph;
class Hashgraph;

struct pre_partition_info {
    HashIntoType kmer;
    SeenSet tagged_kmers;

    explicit pre_partition_info(HashIntoType _kmer) : kmer(_kmer) {};
};

class SubsetPartition
{
    friend class Hashgraph;
protected:
    unsigned int next_partition_id;
    Hashgraph * _ht;
    PartitionMap partition_map;
    ReversePartitionMap reverse_pmap;

    void _clear_all_partitions();

    PartitionID * _merge_two_partitions(PartitionID *orig_pp,
                                        PartitionID *new_pp);
    PartitionID * _join_partitions_by_tags(const SeenSet& tagged_kmers,
                                           const HashIntoType kmer);

public:
    explicit SubsetPartition(Hashgraph * ht);

    ~SubsetPartition()
    {
        _clear_all_partitions();
    }

    PartitionID assign_partition_id(HashIntoType kmer, SeenSet& tagged_kmers);

    void set_partition_id(HashIntoType kmer, PartitionID p);
    void set_partition_id(std::string kmer_s, PartitionID p);
    PartitionID join_partitions(PartitionID orig, PartitionID join);
    PartitionID get_partition_id(std::string kmer_s);
    PartitionID get_partition_id(HashIntoType kmer);

    PartitionID * get_new_partition()
    {
        PartitionID* pp = new PartitionID(next_partition_id);
        next_partition_id++;
        return pp;
    }

    void merge(SubsetPartition *);
    void merge_from_disk(std::string);
    void _merge_from_disk_consolidate(PartitionPtrMap&);

    void save_partitionmap(std::string outfile);
    void load_partitionmap(std::string infile);
    void _validate_pmap();

    void find_all_tags(Kmer start_kmer,
                       SeenSet& tagged_kmers,
                       const SeenSet& all_tags,
                       bool break_on_stop_tags=false,
                       bool stop_big_traversals=false);

    unsigned int sweep_for_tags(const std::string& seq,
                                SeenSet& tagged_kmers,
                                const SeenSet& all_tags,
                                unsigned int range,
                                bool break_on_stop_tags,
                                bool stop_big_traversals);

    void find_all_tags_truncate_on_abundance(Kmer start_kmer,
            SeenSet& tagged_kmers,
            const SeenSet& all_tags,
            BoundedCounterType min_count,
            BoundedCounterType max_count,
            bool break_on_stop_tags=false,
            bool stop_big_traversals=false);

    void do_partition(HashIntoType first_kmer,
                      HashIntoType last_kmer,
                      bool break_on_stop_tags=false,
                      bool stop_big_traversals=false,
                      CallbackFn callback=0,
                      void * callback_data=0);

    void do_partition_with_abundance(HashIntoType first_kmer,
                                     HashIntoType last_kmer,
                                     BoundedCounterType min_count,
                                     BoundedCounterType max_count,
                                     bool break_on_stop_tags=false,
                                     bool stop_big_traversals=false,
                                     CallbackFn callback=0,
                                     void * callback_data=0);

    void count_partitions(size_t& n_partitions,
                          size_t& n_unassigned);

    size_t output_partitioned_file(const std::string &infilename,
                                   const std::string &outputfilename,
                                   bool output_unassigned=false,
                                   CallbackFn callback=0,
                                   void * callback_data=0);

    void partition_sizes(PartitionCountMap &cm,
                         unsigned int& n_unassigned) const;

    void partition_size_distribution(PartitionCountDistribution &d,
                                     unsigned int& n_unassigned) const;

    void partition_average_coverages(PartitionCountMap &cm,
                                     Countgraph * ht) const;

    unsigned long long repartition_largest_partition(unsigned int, unsigned int,
            unsigned int, Countgraph&);

    void repartition_a_partition(const SeenSet& partition_tags);
    void _clear_partition(PartitionID, SeenSet& partition_tags);

    void _merge_other(HashIntoType tag,
                      PartitionID other_partition,
                      PartitionPtrMap& diskp_to_pp);

    void report_on_partitions();
};
}

#endif // SUBSET_HH
