#ifndef HASHTABLE_HH
#define HASHTABLE_HH

#include <fstream>
#include <string>
#include <set>
#include <map>

#include "khmer.hh"
#include "storage.hh"

namespace khmer {
  typedef std::set<HashIntoType> SeenSet;
  typedef std::map<HashIntoType, unsigned int> PartitionMap;
  typedef std::map<unsigned int, SeenSet*> ReversePartitionMap;

  class Hashtable {
  protected:
    const WordLength _ksize;
    const HashIntoType _tablesize;
    HashIntoType bitmask;

    BoundedCounterType * _counts;

    void (*_writelock_acquire)(void * data);
    void (*_writelock_release)(void * data);
    void * _writelock_data;

    void writelock_acquire() {
      if (_writelock_acquire) {
	_writelock_acquire(_writelock_data);
      }
    }

    void writelock_release() {
      if (_writelock_release) {
	_writelock_release(_writelock_data);
      }
    }

    void _allocate_counters() {
      _counts = new BoundedCounterType[_tablesize];
      memset(_counts, 0, _tablesize * sizeof(BoundedCounterType));
    }

  public:
    Hashtable(WordLength ksize, HashIntoType tablesize) :
      _ksize(ksize), _tablesize(tablesize) {
      bitmask = 0;
      for (unsigned int i = 0; i < _ksize; i++) {
	bitmask = (bitmask << 2) | 3;
      }
      _allocate_counters();

      _writelock_acquire = NULL;
      _writelock_release = NULL;
      _writelock_data = NULL;
    }

    ~Hashtable() {
      if (_counts) { delete _counts; _counts = NULL; }
    }

#if 0
    // setter to set the writelock functions.
    void set_writelock_functions(void (*acquire)(void *),
				 void (*release)(void *),
				 void * data) {
      _writelock_acquire = acquire;
      _writelock_release = release;
      _writelock_data = data;
    }
#endif //0

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
	if (_counts[i]) {
	  n++;
	}
      }
      return n;
    }

    void count(const char * kmer) {
      HashIntoType bin = _hash(kmer, _ksize) % _tablesize;
      if (_counts[bin] == MAX_COUNT) { return; }
      writelock_acquire();
      try {
	_counts[bin]++;
      } catch (...) {
	writelock_release();
	throw;
      }
      writelock_release();
    }

    void count(HashIntoType khash) {
      HashIntoType bin = khash % _tablesize;
      if (_counts[bin] == MAX_COUNT) { return; }
      writelock_acquire();
      try {
	_counts[bin]++;
      } catch (...) {
	writelock_release();
	throw;
      }
      writelock_release();
    }

    // get the count for the given k-mer.
    const BoundedCounterType get_count(const char * kmer) const {
      HashIntoType bin = _hash(kmer, _ksize) % _tablesize;
      return _counts[bin];
    }

    // get the count for the given k-mer hash.
    const BoundedCounterType get_count(HashIntoType khash) const {
      HashIntoType bin = khash % _tablesize;
      return _counts[bin];
    }

    // count every k-mer in the string.
    unsigned int consume_string(const std::string &s,
				HashIntoType lower_bound = 0,
				HashIntoType upper_bound = 0);

    // checks each read for non-ACGT characters
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

    void partition_set_id(const HashIntoType kmer_f,
			  const HashIntoType kmer_r,
			  SeenSet& keeper,
			  const unsigned int partition_id,
			  PartitionMap& partition_map);

    void do_partition(const std::string infilename,
				 CallbackFn callback,
				 void * callback_data);
    void do_partition2(const std::string infilename,
				 CallbackFn callback,
				 void * callback_data);

    void partition_find_id(const HashIntoType kmer_f,
			   const HashIntoType kmer_r,
			   SeenSet& keeper,
			   SeenSet& tagged_kmers,
			   PartitionMap& partition_map,
			   ReversePartitionMap& rev_pmap,
			   bool& done);
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
