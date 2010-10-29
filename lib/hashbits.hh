#ifndef HASHBITS_HH
#define HASHBITS_HH

#include <vector>
#include "hashtable.hh"
#include "subset.hh"

namespace khmer {
  class Hashbits : public Hashtable {
    friend class SubsetPartition;
  protected:
    std::vector<HashIntoType> _tablesizes;
    unsigned int n_tables;
    unsigned int _tag_density;
    HashIntoType occupied_bins;
    HashIntoType n_unique_kmers;

    BoundedCounterType ** _counts;
    SeenSet all_tags;

    virtual void _allocate_counters() {
      n_tables = _tablesizes.size();

      HashIntoType tablebytes;
      HashIntoType tablesize;

      _counts = new BoundedCounterType*[n_tables];

      for (unsigned int i = 0; i < n_tables; i++) {
	tablesize = _tablesizes[i];
	tablebytes = tablesize / 8 + 1;

	_counts[i] = new BoundedCounterType[tablebytes];
	memset(_counts[i], 0, tablebytes * sizeof(BoundedCounterType));
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

    Hashbits(WordLength ksize, std::vector<HashIntoType>& tablesizes) :
      Hashtable(ksize, 1), _tablesizes(tablesizes) {
      _tablesize = 0;
      _tag_density = TAG_DENSITY;
      partition = new SubsetPartition(this);
      occupied_bins = 0;
      n_unique_kmers = 0;

      _allocate_counters();
    }

    ~Hashbits() {
      if (_counts) {
	for (unsigned int i = 0; i < n_tables; i++) {
	  delete _counts[i];
	  _counts[i] = NULL;
	}
      }
      n_tables = 0;

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

    void connectivity_distribution(const std::string infilename,
				   HashIntoType dist[9],
				   CallbackFn callback=0,
				   void * callback_data=0);

    void tags_to_map(TagCountMap& tag_map);
    void discard_tags(TagCountMap& tag_map, unsigned int threshold);

    // count number of occupied bins
    virtual const HashIntoType n_occupied(HashIntoType start=0,
				  HashIntoType stop=0) const {
      return occupied_bins/n_tables;
    }
      
    virtual const HashIntoType n_kmers(HashIntoType start=0,
                  HashIntoType stop=0) const {
      return n_unique_kmers;
    }

    virtual void count(const char * kmer) {
      HashIntoType hash = _hash(kmer, _ksize);
      int flag = 1; // if the kmer appears in any hashtable
      for (unsigned int i = 0; i < n_tables; i++) {
	HashIntoType bin = hash % _tablesizes[i];
	unsigned int byte = bin / 8;
	unsigned char bit = bin % 8;
    if (!( _counts[i][byte] & (1<<bit))) {
        occupied_bins += 1;
        flag = 0; // change the value
    }

	_counts[i][byte] |= (1 << bit);
      }
      if (flag == 0) {
          n_unique_kmers +=1;
      }
    }

    virtual void count(HashIntoType khash) {
      int flag = 1;
      for (unsigned int i = 0; i < n_tables; i++) {
	HashIntoType bin = khash % _tablesizes[i];
	unsigned int byte = bin / 8;
	unsigned char bit = bin % 8;
    if (!( _counts[i][byte] & (1<<bit))) {
        occupied_bins += 1;
        flag = 0;
    }
          
	_counts[i][byte] |= (1 << bit);
      }
        if (flag == 0) {
            n_unique_kmers +=1;
        }
    }

    // get the count for the given k-mer.
    virtual const BoundedCounterType get_count(const char * kmer) const {
      HashIntoType hash = _hash(kmer, _ksize);

      for (unsigned int i = 0; i < n_tables; i++) {
	HashIntoType bin = hash % _tablesizes[i];
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
      for (unsigned int i = 0; i < n_tables; i++) {
	HashIntoType bin = khash % _tablesizes[i];
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
