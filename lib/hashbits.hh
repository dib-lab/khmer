#ifndef HASHBITS_HH
#define HASHBITS_HH

#include "hashtable.hh"
#include "subset.hh"

#define N_TABLES 2
#define PARTITION_ALL_TAG_DEPTH 500

namespace khmer {
  extern HashIntoType primes[];

  class Hashbits : public Hashtable {
    friend class SubsetPartition;
  protected:
    BoundedCounterType * _counts[N_TABLES];
    HashIntoType _tablebytes;
    SeenSet all_tags;

    virtual void _allocate_counters() {
      for (unsigned int i = 0; i < N_TABLES; i++) {
	_tablesize = primes[i];
	_tablebytes = _tablesize / 8 + 1;

	_counts[i] = new BoundedCounterType[_tablebytes];
	memset(_counts[i], 0, _tablebytes * sizeof(BoundedCounterType));
      }
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

    Hashbits(WordLength ksize, HashIntoType tablesize) :
      Hashtable(ksize, tablesize) {
      partition = new SubsetPartition(this);

      _allocate_counters();
    }

    ~Hashbits() {
      if (_counts[0]) {
	for (unsigned int i = 0; i < N_TABLES; i++) {
	  delete _counts[i];
	  _counts[i] = NULL;
	}
      }
      _clear_partitions();
    }

    virtual void save(std::string);
    virtual void load(std::string);
    virtual void save_tagset(std::string);
    virtual void load_tagset(std::string);

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

    // count number of occupied bins
    virtual const HashIntoType n_occupied(HashIntoType start=0,
				  HashIntoType stop=0) const {
      HashIntoType n = 0;
      if (stop == 0) { stop = _tablesize; }
      for (HashIntoType i = start; i < stop; i++) {
	unsigned int byte = i / 8;
	unsigned char bit = i % 8;
	if (_counts[0][byte] & (1 << bit)) {
	  n++;
	}
      }
      return n;
    }

    virtual void count(const char * kmer) {
      HashIntoType hash = _hash(kmer, _ksize);

      for (unsigned int i = 0; i < N_TABLES; i++) {
	HashIntoType bin = hash % primes[i];
	unsigned int byte = bin / 8;
	unsigned char bit = bin % 8;

	_counts[i][byte] |= (1 << bit);
      }
    }

    virtual void count(HashIntoType khash) {
      for (unsigned int i = 0; i < N_TABLES; i++) {
	HashIntoType bin = khash % primes[i];
	unsigned int byte = bin / 8;
	unsigned char bit = bin % 8;

	_counts[i][byte] |= (1 << bit);
      }
    }

    // get the count for the given k-mer.
    virtual const BoundedCounterType get_count(const char * kmer) const {
      HashIntoType hash = _hash(kmer, _ksize);

      for (unsigned int i = 0; i < N_TABLES; i++) {
	HashIntoType bin = hash % primes[i];
	unsigned int byte = bin / 8;
	unsigned char bit = bin % 8;
      
	if (!(_counts[i][byte] & (1 << bit))) {
	  return 0;
	}
      }
      return 1;
    }

    // get the count for the given k-mer hash.
    virtual const BoundedCounterType get_count(HashIntoType khash) const {
      for (unsigned int i = 0; i < N_TABLES; i++) {
	HashIntoType bin = khash % primes[i];
	unsigned int byte = bin / 8;
	unsigned char bit = bin % 8;
      
	if (!(_counts[i][byte] & (1 << bit))) {
	  return 0;
	}
      }
      return 1;
    }
  };

};

#endif // HASHBITS_HH
