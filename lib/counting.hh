#ifndef COUNTING_HH
#define COUNTING_HH

#include <vector>
#include "hashtable.hh"

namespace khmer {
  class CountingHashIntersect;

  class CountingHash : public khmer::Hashtable {
    friend class CountingHashIntersect;

  protected:
    HashIntoType _tablesize;

    BoundedCounterType * _counts;

    virtual void _allocate_counters() {
      _counts = new BoundedCounterType[_tablesize];
      memset(_counts, 0, _tablesize * sizeof(BoundedCounterType));
    }

  public:
    CountingHash(WordLength ksize, HashIntoType tablesize) :
      khmer::Hashtable(ksize), _tablesize(tablesize) {

      _allocate_counters();
    }

    virtual ~CountingHash() {
      if (_counts) {
	delete _counts;
	_counts = NULL;
      }
    }

    virtual void save(std::string);
    virtual void load(std::string);

    // accessors to get table info
    const HashIntoType n_entries() const { return _tablesize; }

    // count number of occupied bins
    virtual const HashIntoType n_occupied(HashIntoType start=0,
					  HashIntoType stop=0) const {
      HashIntoType n = 0;
      if (stop == 0) { stop = _tablesize; }
      for (HashIntoType i = start; i < stop; i++) {
	if (_counts[i % _tablesize]) {
	  n++;
	}
      }
      return n;
    }

    virtual void count(const char * kmer) {
      HashIntoType hash = _hash(kmer, _ksize);
      HashIntoType bin = hash % _tablesize;

      if (_counts[bin] < MAX_COUNT) {
	_counts[bin] += 1;
      }
    }

    virtual void count(HashIntoType khash) {
      HashIntoType bin = khash % _tablesize;

      if (_counts[bin] < MAX_COUNT) {
	_counts[bin] += 1;
      }
    }

    // get the count for the given k-mer.
    virtual const BoundedCounterType get_count(const char * kmer) const {
      HashIntoType hash = _hash(kmer, _ksize);

      HashIntoType bin = hash % _tablesize;
      return _counts[bin];
    }

    // get the count for the given k-mer hash.
    virtual const BoundedCounterType get_count(HashIntoType khash) const {
      HashIntoType bin = khash % _tablesize;

      return _counts[bin];
    }

    //

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
  };

  class CountingHashIntersect {
  protected:
    khmer::CountingHash * _kh1;
    khmer::CountingHash * _kh2;

  public:
    CountingHashIntersect(WordLength ksize,
		       HashIntoType tablesize1, HashIntoType tablesize2)
    {
      _kh1 = new CountingHash(ksize, tablesize1);
      _kh2 = new CountingHash(ksize, tablesize2);
    }

    ~CountingHashIntersect()
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

#endif // COUNTING_HH
