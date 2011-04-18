#ifndef HASHBITS_HH
#define HASHBITS_HH

#include <vector>
#include "hashtable.hh"
#include "subset.hh"
#include "counting.hh"

#define next_f(kmer_f, ch) ((((kmer_f) << 2) & bitmask) | twobit_repr(ch))
#define next_r(kmer_r, ch) (((kmer_r) >> 2) | (twobit_comp(ch) << rc_left_shift))

#define prev_f(kmer_f, ch) ((kmer_f) >> 2 | twobit_repr(ch) << rc_left_shift)
#define prev_r(kmer_r, ch) (((kmer_r) << 2) & bitmask | twobit_comp(ch));

namespace khmer {
  class Hashbits : public khmer::Hashtable {
    friend class SubsetPartition;
  protected:
    std::vector<HashIntoType> _tablesizes;
    unsigned int _n_tables;
    unsigned int _tag_density;
    HashIntoType _occupied_bins;
    HashIntoType _n_unique_kmers;

    Byte ** _counts;
    SeenSet all_tags;
    SeenSet stop_tags;

    virtual void _allocate_counters() {
      _n_tables = _tablesizes.size();

      HashIntoType tablebytes;
      HashIntoType tablesize;

      _counts = new Byte*[_n_tables];

      for (unsigned int i = 0; i < _n_tables; i++) {
	tablesize = _tablesizes[i];
	tablebytes = tablesize / 8 + 1;

	_counts[i] = new Byte[tablebytes];
	memset(_counts[i], 0, tablebytes);
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
      khmer::Hashtable(ksize), _tablesizes(tablesizes) {
      _tag_density = DEFAULT_TAG_DENSITY;
      assert(_tag_density % 2 == 0);
      partition = new SubsetPartition(this);
      _occupied_bins = 0;
      _n_unique_kmers = 0;

      _allocate_counters();
    }

    ~Hashbits() {
      if (_counts) {
	for (unsigned int i = 0; i < _n_tables; i++) {
	  delete _counts[i];
	  _counts[i] = NULL;
	}
	delete _counts;
	_counts = NULL;

	_n_tables = 0;
      }

      _clear_partitions();
    }

    virtual void save(std::string);
    virtual void load(std::string);
    virtual void save_tagset(std::string);
    virtual void load_tagset(std::string, bool clear_tags=true);

    // for debugging/testing purposes only!
    void _set_tag_density(unsigned int d) {
      assert(d % 2 == 0);	// must be even
      assert(all_tags.size() == 0); // no tags exist!
      _tag_density = d;
    }

    unsigned int _get_tag_density() const {
      return _tag_density;
    }

    void add_tag(HashIntoType tag) { all_tags.insert(tag); }
    void add_stop_tag(HashIntoType tag) { stop_tags.insert(tag); }

    void calc_connected_graph_size(const char * kmer,
				   unsigned long long& count,
				   SeenSet& keeper,
				   const unsigned long long threshold=0,
				   bool break_on_circum=false) const{
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

    unsigned int kmer_degree(HashIntoType kmer_f, HashIntoType kmer_r) const;
    unsigned int kmer_degree(const char * kmer_s) const {
      HashIntoType kmer_f, kmer_r;
      _hash(kmer_s, _ksize, kmer_f, kmer_r);

      return kmer_degree(kmer_f, kmer_r);
    }

    // count number of occupied bins
    virtual const HashIntoType n_occupied(HashIntoType start=0,
				  HashIntoType stop=0) const {
      // @@ CTB need to be able to *save* this...
      return _occupied_bins/_n_tables;
    }
      
    virtual const HashIntoType n_kmers(HashIntoType start=0,
                  HashIntoType stop=0) const {
      return _n_unique_kmers;	// @@ CTB need to be able to *save* this...
    }

    virtual void count(const char * kmer) {
      HashIntoType hash = _hash(kmer, _ksize);
      count(hash);
    }

    virtual void count(HashIntoType khash) {
      bool is_new_kmer = false;

      for (unsigned int i = 0; i < _n_tables; i++) {
	HashIntoType bin = khash % _tablesizes[i];
	HashIntoType byte = bin / 8;
	unsigned char bit = bin % 8;
	if (!( _counts[i][byte] & (1<<bit))) {
	  _occupied_bins += 1;
	  is_new_kmer = true;
	}
	_counts[i][byte] |= (1 << bit);
      }
      if (is_new_kmer) {
	_n_unique_kmers +=1;
      }
    }

    // get the count for the given k-mer.
    virtual const BoundedCounterType get_count(const char * kmer) const {
      HashIntoType hash = _hash(kmer, _ksize);
      return get_count(hash);
    }

    // get the count for the given k-mer hash.
    virtual const BoundedCounterType get_count(HashIntoType khash) const {
      for (unsigned int i = 0; i < _n_tables; i++) {
	HashIntoType bin = khash % _tablesizes[i];
	HashIntoType byte = bin / 8;
	unsigned char bit = bin % 8;
      
	if (!(_counts[i][byte] & (1 << bit))) {
	  return 0;
	}
      }
      return 1;
    }

    void filter_if_present(const std::string infilename,
			   const std::string outputfilename,
			   CallbackFn callback=0,
			   void * callback_data=0);

    unsigned int count_kmers_within_radius(HashIntoType kmer_f,
					   HashIntoType kmer_r,
					   unsigned int radius,
					   unsigned int max_count,
					   const SeenSet * seen=0) const;
    unsigned int count_kmers_within_depth(HashIntoType kmer_f,
					  HashIntoType kmer_r,
					  unsigned int depth,
					  unsigned int max_count,
					  SeenSet * seen) const;

    unsigned int find_radius_for_volume(HashIntoType kmer_f,
					HashIntoType kmer_r,
					unsigned int max_count,
					unsigned int max_radius) const;

    unsigned int count_kmers_on_radius(HashIntoType kmer_f,
				       HashIntoType kmer_r,
				       unsigned int radius,
				       unsigned int max_volume) const;

    unsigned int trim_on_degree(std::string sequence, unsigned int max_degree)
      const;
    unsigned int trim_on_sodd(std::string sequence, unsigned int max_degree)
      const;

    unsigned int trim_on_density_explosion(std::string sequence, unsigned int radius, unsigned int max_volume)
      const;

    void load_stop_tags(std::string filename, bool clear_tags=true);
    unsigned int trim_on_stoptags(std::string sequence) const;

    void traverse_from_tags(unsigned int distance,
			    unsigned int frequency,
			    CountingHash &counting) const;
    unsigned int _traverse_from_tag(HashIntoType start,
			    unsigned int radius,
			    CountingHash &counting) const;
    void hitraverse_to_stoptags(CountingHash &counting,
				unsigned int cutoff);
  };
};

#endif // HASHBITS_HH
