#ifndef HASHBITS_HH
#define HASHBITS_HH

#include <vector>
#include "hashtable.hh"
#include "subset.hh"

#define next_f(kmer_f, ch) ((((kmer_f) << 2) & bitmask) | (twobit_repr(ch)))
#define next_r(kmer_r, ch) (((kmer_r) >> 2) | (twobit_comp(ch) << rc_left_shift))

#define prev_f(kmer_f, ch) ((kmer_f) >> 2 | twobit_repr(ch) << rc_left_shift)
#define prev_r(kmer_r, ch) ((((kmer_r) << 2) & bitmask) | (twobit_comp(ch)))

#define set_contains(s, e) ((s).find(e) != (s).end())

namespace khmer {
  class CountingHash;

  class Hashbits : public khmer::Hashtable {
    friend class SubsetPartition;
  protected:
    std::vector<HashIntoType> _tablesizes;
    unsigned int _n_tables;
    unsigned int _tag_density;
    HashIntoType _occupied_bins;
    HashIntoType _n_unique_kmers;
	HashIntoType _n_overlap_kmers;
    Byte ** _counts;

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
            
    void _clear_all_partitions() {
      if (partition != NULL) {
	partition->_clear_all_partitions();
      }
    }

    uint32_t _all_tags_spin_lock;

  public:
    SubsetPartition * partition;
    SeenSet all_tags;
    SeenSet stop_tags;
    SeenSet repart_small_tags;

    void _validate_pmap() {
      if (partition) { partition->_validate_pmap(); }
    }

    Hashbits(WordLength ksize, std::vector<HashIntoType>& tablesizes)
    : khmer::Hashtable(ksize),
      _tablesizes(tablesizes),
      _all_tags_spin_lock( 0 )
    {
      _tag_density = DEFAULT_TAG_DENSITY;
      assert(_tag_density % 2 == 0);
      partition = new SubsetPartition(this);
      _occupied_bins = 0;
      _n_unique_kmers = 0;
      _n_overlap_kmers = 0;

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

      _clear_all_partitions();
    }

    std::vector<HashIntoType> get_tablesizes() const {
      return _tablesizes;
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

    // Count every k-mer in a FASTA or FASTQ file.
    // Tag certain ones on the connectivity graph.
    void consume_fasta_and_tag(
      std::string const	  &filename,
      unsigned int	  &total_reads,
      unsigned long long  &n_consumed,
      CallbackFn	  callback	  = NULL,
      void *		  callback_data	  = NULL
    );
    // Count every k-mer from a stream of FASTA or FASTQ reads, 
    // using the supplied parser.
    // Tag certain ones on the connectivity graph.
    void consume_fasta_and_tag(
	read_parsers:: IParser *	    parser,
	unsigned int	    &total_reads,
	unsigned long long  &n_consumed,
	CallbackFn	    callback	    = NULL,
	void *		    callback_data   = NULL
    );

    void consume_sequence_and_tag(const std::string& seq,
				  unsigned long long& n_consumed,
				  SeenSet * new_tags = 0);


    void consume_fasta_and_tag_with_stoptags(const std::string &filename,
					     unsigned int &total_reads,
					     unsigned long long &n_consumed,
					     CallbackFn callback = 0,
					     void * callback_data = 0);

    void consume_fasta_and_traverse(const std::string &filename,
				    unsigned int distance,
				    unsigned int big_threshold,
				    unsigned int transfer_threshold,
				    CountingHash &counting);

    void consume_partitioned_fasta(const std::string &filename,
				   unsigned int &total_reads,
				   unsigned long long &n_consumed,
				   CallbackFn callback = 0,
				   void * callback_data = 0);

    // for overlap k-mer counting
    void consume_fasta_overlap(const std::string &filename,HashIntoType curve[2][100],
                              khmer::Hashbits &ht2,
			      unsigned int &total_reads,
			      unsigned long long &n_consumed,
			      HashIntoType lower_bound,
			      HashIntoType upper_bound,
			      CallbackFn callback,
			      void * callback_data);



    // just for overlap k-mer counting!
    unsigned int check_and_process_read_overlap(std::string &read,
					    bool &is_valid,HashIntoType lower_bound,
                                            HashIntoType upper_bound,
                                            khmer::Hashbits &ht2);
    // for overlap k-mer counting!
    unsigned int consume_string_overlap(const std::string &s,
				       HashIntoType lower_bound,
				       HashIntoType upper_bound,khmer::Hashbits &ht2);




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

    // Get and set the hashbits for the given kmer.
    inline
    virtual
    const
    BoundedCounterType
    test_and_set_bits(const char * kmer)
    {
      HashIntoType hash = _hash(kmer, _ksize);
      return test_and_set_bits(hash);
    }

    // Get and set the hashbits for the given kmer hash.
    // Generally, it is better to keep tests and mutations separate, 
    // but, in the interests of efficiency and thread safety, 
    // tests and mutations are being blended here against conventional 
    // software engineering wisdom.
    inline
    virtual
    const
    bool
    test_and_set_bits( HashIntoType khash ) 
    {
      bool is_new_kmer = false;

      for (unsigned int i = 0; i < _n_tables; i++)
      {
        HashIntoType bin = khash % _tablesizes[i];
	HashIntoType byte = bin / 8;
	unsigned char bit = (unsigned char)(1 << (bin % 8));

	unsigned char bits_orig = __sync_fetch_and_or( *(_counts + i) + byte, bit );
	if (!(bits_orig & bit))
	{
	  __sync_add_and_fetch( &_occupied_bins, 1 );
	  is_new_kmer = true;
	}
      } // iteration over hashtables

      if (is_new_kmer)
      {
	__sync_add_and_fetch( &_n_unique_kmers, 1 );
	return true; // kmer not seen before
      }

      return false; // kmer already seen
    } // test_and_set_bits

    virtual const HashIntoType n_overlap_kmers(HashIntoType start=0,
                  HashIntoType stop=0) const {
      return _n_overlap_kmers;	// @@ CTB need to be able to *save* this...
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

	virtual bool check_overlap(HashIntoType khash, Hashbits &ht2) {

	  for (unsigned int i = 0; i < ht2._n_tables; i++) {
		HashIntoType bin = khash % ht2._tablesizes[i];
		HashIntoType byte = bin / 8;
		unsigned char bit = bin % 8;
		if (!( ht2._counts[i][byte] & (1<<bit))) {
		  return false;
	}
      }
	  return true;
	  }

    virtual void count_overlap(const char * kmer, Hashbits &ht2) {
      HashIntoType hash = _hash(kmer, _ksize);
      count_overlap(hash,ht2);
    }

    virtual void count_overlap(HashIntoType khash, Hashbits &ht2) {
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
	if (check_overlap(khash,ht2)){
		_n_overlap_kmers +=1;
      }
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

    unsigned int trim_on_stoptags(std::string sequence) const;

    void traverse_from_tags(unsigned int distance,
			    unsigned int threshold,
			    unsigned int num_high_todo,
			    CountingHash &counting);

    unsigned int traverse_from_kmer(HashIntoType start,
				    unsigned int radius,
				    SeenSet &keeper) const;

    unsigned int count_and_transfer_to_stoptags(SeenSet &keeper,
						unsigned int threshold,
						CountingHash &counting);

    void traverse_from_reads(std::string filename,
			     unsigned int radius,
			     unsigned int big_threshold,
			     unsigned int transfer_threshold,
			     CountingHash &counting);

    void hitraverse_to_stoptags(std::string filename,
				CountingHash &counting,
				unsigned int cutoff);

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
  };
};

#include "counting.hh"

#endif // HASHBITS_HH

// vim: set sts=2 sw=2:
