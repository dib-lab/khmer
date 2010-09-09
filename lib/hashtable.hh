#ifndef HASHTABLE_HH
#define HASHTABLE_HH

#include <fstream>
#include <string>
#include <set>
#include <map>
#include <queue>

#include "khmer.hh"
#include "storage.hh"

#define SURRENDER_PARTITION 1

namespace khmer {
  typedef unsigned int PartitionID;
  typedef std::set<HashIntoType> SeenSet;
  typedef std::set<PartitionID> PartitionSet;
  typedef std::map<HashIntoType, PartitionID*> PartitionMap;
  typedef std::map<PartitionID, PartitionID*> PartitionPtrMap;
  typedef std::map<PartitionID, SeenSet*> PartitionsToTagsMap;
  typedef std::set<PartitionID *> PartitionPtrSet;
  typedef std::map<PartitionID, PartitionPtrSet*> ReversePartitionMap;
  typedef std::queue<HashIntoType> NodeQueue;
  typedef std::map<PartitionID, PartitionID*> PartitionToPartitionPMap;
  typedef std::map<HashIntoType, unsigned int> TagCountMap;
  typedef std::map<PartitionID, unsigned int> PartitionCountMap;

  class Hashtable;

  class SubsetPartition {
    friend class Hashtable;
  protected:
    unsigned int next_partition_id;
    Hashtable * _ht;
    PartitionMap partition_map;
    ReversePartitionMap reverse_pmap;

    void _clear_partitions();

    void _add_partition_ptr(PartitionID *orig_pp, PartitionID *new_pp);
    PartitionID * _reassign_partition_ids(SeenSet& tagged_kmers,
					const HashIntoType kmer_f);

    bool _is_tagged_kmer(const HashIntoType kmer_f,
			 const HashIntoType kmer_r,
			 HashIntoType& tagged_kmer);

    bool _do_continue(const HashIntoType kmer,
		      const SeenSet& keeper);

  public:
    SubsetPartition(Hashtable * ht) : next_partition_id(2), _ht(ht) {
      PartitionPtrSet * s = new PartitionPtrSet();
      PartitionID * p = new PartitionID(SURRENDER_PARTITION);
      s->insert(p);
      reverse_pmap[SURRENDER_PARTITION] = s;
    };

    ~SubsetPartition() { _clear_partitions(); }

    PartitionID assign_partition_id(HashIntoType kmer_f,
				    SeenSet& tagged_kmers,
				    bool surrender);

    void set_partition_id(HashIntoType kmer_f, PartitionID p);
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

    void save_partitionmap(std::string outfile);
    void load_partitionmap(std::string infile);
    void _validate_pmap();

    void find_all_tags(HashIntoType kmer_f, HashIntoType kmer_r,
		       SeenSet& tagged_kmers,
		       bool& surrender, bool do_initial_check);

    void do_partition(HashIntoType first_kmer,
		      HashIntoType last_kmer,
		      CallbackFn callback=0,
		      void * callback_data=0);

    void count_partitions(unsigned int& n_partitions,
			  unsigned int& n_unassigned,
			  unsigned int& n_surrendered);

    unsigned int output_partitioned_file(const std::string infilename,
					 const std::string outputfilename,
					 bool output_unassigned=false,
					 CallbackFn callback=0,
					 void * callback_data=0);

    void maxify_partition_size(TagCountMap& tag_map);
    void filter_against_tags(TagCountMap& tag_map);
  };

  class Hashtable {
    friend class SubsetPartition;
  protected:
    WordLength _ksize;
    HashIntoType _tablesize;
    HashIntoType _tablebytes;
    HashIntoType bitmask;

    BoundedCounterType * _counts;

    SeenSet all_tags;

    void _allocate_counters() {
      _tablebytes = _tablesize / 8 + 1;
      _counts = new BoundedCounterType[_tablebytes];
      memset(_counts, 0, _tablebytes * sizeof(BoundedCounterType));
    }

    void _clear_partitions() {
      if (partition != NULL) {
	partition->_clear_partitions();
      }
    }

  public:
    SubsetPartition * partition;

    void _validate_pmap() {
      if (partition) { partition->_validate_pmap(); }
    }

    Hashtable(WordLength ksize, HashIntoType tablesize) :
      _ksize(ksize), _tablesize(tablesize) {
      partition = new SubsetPartition(this);

      bitmask = 0;
      for (unsigned int i = 0; i < _ksize; i++) {
	bitmask = (bitmask << 2) | 3;
      }
      _allocate_counters();
    }

    ~Hashtable() {
      if (_counts) { delete _counts; _counts = NULL; }
      _clear_partitions();
    }

    void save(std::string);
    void load(std::string);
    void save_tagset(std::string);
    void load_tagset(std::string);

    // accessor to get 'k'
    const WordLength ksize() const { return _ksize; }

    // accessors to get table info
    const HashIntoType n_entries() const { return _tablesize; }

    // count number of occupied bins
    const HashIntoType n_occupied(HashIntoType start=0,
				  HashIntoType stop=0) const {
      HashIntoType n = 0;
      if (stop == 0) { stop = _tablesize; }
      for (HashIntoType i = start; i < stop; i++) {
	unsigned int byte = i / 8;
	unsigned char bit = i % 8;
	if (_counts[byte] & (1 << bit)) {
	  n++;
	}
      }
      return n;
    }

    void count(const char * kmer) {
      HashIntoType bin = _hash(kmer, _ksize) % _tablesize;
      unsigned int byte = bin / 8;
      unsigned char bit = bin % 8;

      _counts[byte] |= (1 << bit);
    }

    void count(HashIntoType khash) {
      HashIntoType bin = khash % _tablesize;
      unsigned int byte = bin / 8;
      unsigned char bit = bin % 8;

      _counts[byte] |= (1 << bit);
    }

    // get the count for the given k-mer.
    const BoundedCounterType get_count(const char * kmer) const {
      HashIntoType bin = _hash(kmer, _ksize) % _tablesize;
      unsigned int byte = bin / 8;
      unsigned char bit = bin % 8;
      
      if (_counts[byte] & (1 << bit)) {
	return 1;
      }
      return 0;
    }

    // get the count for the given k-mer hash.
    const BoundedCounterType get_count(HashIntoType khash) const {
      HashIntoType bin = khash % _tablesize;
      unsigned int byte = bin / 8;
      unsigned char bit = bin % 8;
      
      if (_counts[byte] & (1 << bit)) {
	return 1;
      }
      return 0;
    }

    // count every k-mer in the string.
    unsigned int consume_string(const std::string &s,
				HashIntoType lower_bound = 0,
				HashIntoType upper_bound = 0);

    // checks each read for non-ACGT characters
    bool check_read(const std::string &read);

    // check each read for non-ACGT characters, and then consume it.
    unsigned int check_and_process_read(const std::string &read,
					bool &is_valid,
					HashIntoType lower_bound = 0,
					HashIntoType upper_bound = 0);

    // count every k-mer in the FASTA file.
    void consume_fasta(const std::string &filename,
		       unsigned int &total_reads,
		       unsigned long long &n_consumed,
		       HashIntoType lower_bound = 0,
		       HashIntoType upper_bound = 0,
		       ReadMaskTable ** readmask = NULL,
		       bool update_readmask = true,
		       CallbackFn callback = NULL,
		       void * callback_data = NULL);

    MinMaxTable * fasta_file_to_minmax(const std::string &inputfile,
				       unsigned int total_reads,
				       ReadMaskTable * readmask = NULL,
				       CallbackFn callback = NULL,
				       void * callback_data = NULL);

    ReadMaskTable * filter_fasta_file_any(MinMaxTable &minmax,
					  BoundedCounterType threshold,
					  ReadMaskTable * readmask = NULL,
					  CallbackFn callback = NULL,
					  void * callback_data = NULL);

    ReadMaskTable * filter_fasta_file_all(MinMaxTable &minmax,
					  BoundedCounterType threshold,
					  ReadMaskTable * readmask = NULL,
					  CallbackFn callback = NULL,
					  void * callback_data = NULL);

    ReadMaskTable * filter_fasta_file_limit_n(const std::string &readsfile,
                                              MinMaxTable &minmax,
                                              BoundedCounterType threshold,
                                              BoundedCounterType n, 
                                              ReadMaskTable * old_readmask = NULL,
                                              CallbackFn callback = NULL,
                                              void * callback_data = NULL);

    ReadMaskTable * filter_fasta_file_run(const std::string &inputfile,
					  unsigned int total_reads,
					  BoundedCounterType threshold,
					  unsigned int runlength,
					  ReadMaskTable * old_readmask = NULL,
					  CallbackFn callback = NULL,
					  void * callback_data = NULL);

    void output_fasta_kmer_pos_freq(const std::string &inputfile,
                                    const std::string &outputfile);

    BoundedCounterType get_min_count(const std::string &s,
				     HashIntoType lower_bound = 0,
				     HashIntoType upper_bound = 0);
				     
    BoundedCounterType get_max_count(const std::string &s,
				     HashIntoType lower_bound = 0,
				     HashIntoType upper_bound = 0);

    HashIntoType * abundance_distribution() const;

    HashIntoType * fasta_count_kmers_by_position(const std::string &inputfile,
					 const unsigned int max_read_len,
					 ReadMaskTable * old_readmask = NULL,
					 BoundedCounterType limit_by_count=0,
						 CallbackFn callback = NULL,
						 void * callback_data = NULL);

    void fasta_dump_kmers_by_abundance(const std::string &inputfile,
				       ReadMaskTable * readmask,
				       BoundedCounterType limit_by_count,
				       CallbackFn callback = NULL,
				       void * callback_data = NULL);

    void trim_graphs(const std::string infilename,
		     const std::string outfilename,
		     unsigned int min_size,
		     CallbackFn callback = NULL,
		     void * callback_data = NULL);

    HashIntoType * graphsize_distribution(const unsigned int &max_size);

    ReadMaskTable * filter_file_connected(const std::string &est,
                                          const std::string &readsfile,
                                          unsigned int total_reads);

    void calc_connected_graph_size(const char * kmer,
				   unsigned long long& count,
				   SeenSet& keeper,
				   const unsigned long long threshold=0) const{
      HashIntoType r, f;
      _hash(kmer, _ksize, f, r);
      calc_connected_graph_size(f, r, count, keeper, threshold);
    }

    void calc_connected_graph_size(const HashIntoType kmer_f,
				   const HashIntoType kmer_r,
				   unsigned long long& count,
				   SeenSet& keeper,
				   const unsigned long long threshold=0) const;

    typedef void (*kmer_cb)(const char * k, unsigned int n_reads, void *data);

    void dump_kmers_and_counts(kmer_cb cb_fn = NULL, void * data = NULL) const {
      for (HashIntoType i = 0; i < _tablesize; i++) {
	BoundedCounterType count = _counts[i] & 127;
	if (_counts[i]) {
	  if (cb_fn) {
	    cb_fn(_revhash(i, _ksize).c_str(), count, data);
	  } else{
	    std::cout << _revhash(i, _ksize) << " " << count << std::endl;
	  }
	}
      }
    }

    // Partitioning stuff.

    unsigned int n_tags() const { return all_tags.size(); }

    void divide_tags_into_subsets(unsigned int subset_size, SeenSet& divvy);

    void add_kmer_to_tags(HashIntoType kmer) {
      all_tags.insert(kmer);
    }

    void clear_tags() { all_tags.clear(); }

    void consume_fasta_and_tag(const std::string &filename,
			       unsigned int &total_reads,
			       unsigned long long &n_consumed,
			       CallbackFn callback = 0,
			       void * callback_data = 0);

    void consume_partitioned_fasta(const std::string &filename,
				   unsigned int &total_reads,
				   unsigned long long &n_consumed,
				   CallbackFn callback = 0,
				   void * callback_data = 0);

    void do_truncated_partition(const std::string infilename,
				CallbackFn callback=0,
				void * callback_data=0);
    void do_threaded_partition(const std::string infilename,
			       CallbackFn callback=0,
			       void * callback_data=0);

    void tags_to_map(TagCountMap& tag_map);
    void discard_tags(TagCountMap& tag_map, unsigned int threshold);
  };

  class HashtableIntersect {
  protected:
    khmer::Hashtable * _kh1;
    khmer::Hashtable * _kh2;

  public:
    HashtableIntersect(WordLength ksize,
		       HashIntoType tablesize1, HashIntoType tablesize2)
    {
      _kh1 = new Hashtable(ksize, tablesize1);
      _kh2 = new Hashtable(ksize, tablesize2);
    }

    ~HashtableIntersect()
    {
      delete _kh1;
      delete _kh2;
    }

    // count every k-mer in the string.
    void consume_string(const std::string &s)
    {
      _kh1->consume_string(s);
      _kh2->consume_string(s);
    }

    BoundedCounterType get_min_count(const std::string &s)
    {
      BoundedCounterType kh1Min = _kh1->get_min_count(s);
      BoundedCounterType kh2Min = _kh2->get_min_count(s);

      if (kh1Min < kh2Min) {
        return kh1Min;
      } else {
        return kh2Min;
      }
    }

    BoundedCounterType get_max_count(const std::string &s)
    {
      BoundedCounterType kh1Max = _kh1->get_max_count(s);
      BoundedCounterType kh2Max = _kh2->get_max_count(s);

      if (kh1Max > kh2Max) {
        return kh1Max;
      } else {
        return kh2Max;
      }
    }
  };
};

#endif // HASHTABLE_HH
