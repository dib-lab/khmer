#ifndef SUBSET_HH
#define SUBSET_HH

#include "hashtable.hh"

namespace khmer {
  class CountingHash;
  class Hashbits;

  class pre_partition_info {
  public:
    HashIntoType kmer;
    SeenSet tagged_kmers;

    pre_partition_info(HashIntoType _kmer) : kmer(_kmer) {};
  };

  class SubsetPartition {
    friend class Hashbits;
  protected:
    unsigned int next_partition_id;
    Hashbits * _ht;
    PartitionMap partition_map;
    ReversePartitionMap reverse_pmap;

    void _clear_all_partitions();

    PartitionID * _merge_two_partitions(PartitionID *orig_pp,
					PartitionID *new_pp);
    PartitionID * _join_partitions_by_tags(const SeenSet& tagged_kmers,
					   const HashIntoType kmer);

  public:
    SubsetPartition(Hashbits * ht) : next_partition_id(2), _ht(ht) {
      ;
    };

    ~SubsetPartition() { _clear_all_partitions(); }

    PartitionID assign_partition_id(HashIntoType kmer, SeenSet& tagged_kmers);

    void set_partition_id(HashIntoType kmer, PartitionID p);
    void set_partition_id(std::string kmer_s, PartitionID p);
    PartitionID join_partitions(PartitionID orig, PartitionID join);
    PartitionID get_partition_id(std::string kmer_s);
    PartitionID get_partition_id(HashIntoType kmer);

    PartitionID * get_new_partition() {
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

    void find_all_tags(HashIntoType kmer_f, HashIntoType kmer_r,
		       SeenSet& tagged_kmers,
		       const SeenSet& all_tags,
		       bool break_on_stop_tags=false);

    void do_partition(HashIntoType first_kmer,
		      HashIntoType last_kmer,
		      bool break_on_stop_tags=false,
		      CallbackFn callback=0,
		      void * callback_data=0);

    void count_partitions(unsigned int& n_partitions,
			  unsigned int& n_unassigned);

    unsigned int output_partitioned_file(const std::string infilename,
					 const std::string outputfilename,
					 bool output_unassigned=false,
					 CallbackFn callback=0,
					 void * callback_data=0);

    bool is_single_partition(std::string sequence);

    void join_partitions_by_path(std::string sequence);

    void partition_size_distribution(PartitionCountDistribution &d,
				    unsigned int& n_unassigned) const;

    unsigned int repartition_largest_partition(unsigned int, unsigned int,
					       unsigned int, CountingHash&);

    void repartition_a_partition(const SeenSet& partition_tags);
    void _clear_partition(PartitionID, SeenSet& partition_tags);

    void _merge_other(HashIntoType tag,
		      PartitionID other_partition,
		      PartitionPtrMap& diskp_to_pp);
  };
}

#endif // SUBSET_HH
