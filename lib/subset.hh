//
// This file is part of khmer, https://github.com/dib-lab/khmer/, and is
// Copyright (C) Michigan State University, 2009-2015. It is licensed under
// the three-clause BSD license; see LICENSE.
// Contact: khmer-project@idyll.org
//

#ifndef SUBSET_HH
#define SUBSET_HH

#include "khmer.hh"

namespace khmer
{
class CountingHash;
class Hashtable;
class Hashbits;

struct pre_partition_info {
    HashIntoType kmer;
    SeenSet tagged_kmers;

    pre_partition_info(HashIntoType _kmer) : kmer(_kmer) {};
};

class SubsetPartition
{
    friend class Hashtable;
protected:
    unsigned int next_partition_id;
    Hashtable * _ht;
    PartitionMap partition_map;
    ReversePartitionMap reverse_pmap;

    void _clear_all_partitions();

    PartitionID * _merge_two_partitions(PartitionID *orig_pp,
                                        PartitionID *new_pp);
    PartitionID * _join_partitions_by_tags(const SeenSet& tagged_kmers,
                                           const HashIntoType kmer);

public:
    SubsetPartition(Hashtable * ht) : next_partition_id(2), _ht(ht)
    {
        ;
    };

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

    void queue_neighbors(HashIntoType kmer_f,
                         HashIntoType kmer_r,
                         unsigned int breadth,
                         SeenSet& traversed_kmers,
                         NodeQueue& node_q,
                         std::queue<unsigned int>& breadth_q);

    void find_all_tags(HashIntoType kmer_f, HashIntoType kmer_r,
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

    void find_all_tags_truncate_on_abundance(HashIntoType kmer_f,
            HashIntoType kmer_r,
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

    unsigned int find_unpart(const std::string &infilename,
                             bool traverse,
                             bool stop_big_traversals,
                             CallbackFn callback=0,
                             void * callback_data=0);

    bool is_single_partition(std::string sequence);

    void join_partitions_by_path(std::string sequence);

    void partition_sizes(PartitionCountMap &cm,
                         unsigned int& n_unassigned) const;

    void partition_size_distribution(PartitionCountDistribution &d,
                                     unsigned int& n_unassigned) const;

    void partition_average_coverages(PartitionCountMap &cm,
                                     CountingHash * ht) const;

    unsigned long long repartition_largest_partition(unsigned int, unsigned int,
            unsigned int, CountingHash&);

    void repartition_a_partition(const SeenSet& partition_tags);
    void _clear_partition(PartitionID, SeenSet& partition_tags);

    void _merge_other(HashIntoType tag,
                      PartitionID other_partition,
                      PartitionPtrMap& diskp_to_pp);

    void report_on_partitions();

    void compare_to_partition(PartitionID, SubsetPartition *, PartitionID,
                              unsigned int &n_only1,
                              unsigned int &n_only2,
                              unsigned int &n_shared);
};
}

#endif // SUBSET_HH
