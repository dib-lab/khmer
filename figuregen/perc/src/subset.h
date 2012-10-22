#ifndef SUBSET_HH
#define SUBSET_HH

#include "hashtable.h"

namespace khmer {
  class Hashbits;

  class SubsetPartition {
    friend class Hashbits;
  protected:
    unsigned int next_partition_id;
    Hashbits * _ht;
    PartitionMap partition_map;
    ReversePartitionMap reverse_pmap;

    void _clear_partitions();

    void _add_partition_ptr(PartitionID *orig_pp, PartitionID *new_pp);
    PartitionID * _add_partition_ptr2(PartitionID *orig_pp, PartitionID *new_pp);
    PartitionID * _reassign_partition_ids(SeenSet& tagged_kmers,
					  const HashIntoType kmer);

  public:
    SubsetPartition(Hashbits * ht) : next_partition_id(2), _ht(ht) {
      ;
    };

    ~SubsetPartition() { _clear_partitions(); }

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
		       SeenSet& tagged_kmers);

    void do_partition(HashIntoType first_kmer,
		      HashIntoType last_kmer,
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
  };
}

#endif // SUBSET_HH
