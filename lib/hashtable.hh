#ifndef HASHTABLE_HH
#define HASHTABLE_HH

#include <fstream>
#include <string>

#include "khmer.hh"
#include "storage.hh"

namespace khmer {
  class Hashtable {
  protected:
    const WordLength _ksize;
    const HashIntoType _tablesize;

    BoundedCounterType * _counts;

    void _allocate_counters() {
      _counts = new BoundedCounterType[_tablesize];
      memset(_counts, 0, _tablesize * sizeof(BoundedCounterType));
    }

  public:
    Hashtable(WordLength ksize, HashIntoType tablesize) :
      _ksize(ksize), _tablesize(tablesize) {
      _allocate_counters();
    }

    ~Hashtable() {
      if (_counts) { delete _counts; _counts = NULL; }
    }

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
      _counts[bin]++;
    }

    void count(HashIntoType khash) {
      HashIntoType bin = khash % _tablesize;
      if (_counts[bin] == MAX_COUNT) { return; }
      _counts[bin]++;
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

    ReadMaskTable * filter_fasta_file_any(const std::string &inputfile,
					  MinMaxTable &minmax,
					  BoundedCounterType threshold,
					  ReadMaskTable * readmask = NULL,
					  CallbackFn callback = NULL,
					  void * callback_data = NULL);

    ReadMaskTable * filter_fasta_file_all(const std::string &inputfile,
					  MinMaxTable &minmax,
					  BoundedCounterType threshold,
					  ReadMaskTable * readmask = NULL,
					  CallbackFn callback = NULL,
					  void * callback_data = NULL);

    ReadMaskTable * filter_fasta_file_run(const std::string &inputfile,
					  unsigned int total_reads,
					  BoundedCounterType threshold,
					  unsigned int runlength,
					  ReadMaskTable * old_readmask = NULL,
					  CallbackFn callback = NULL,
					  void * callback_data = NULL);

    BoundedCounterType get_min_count(const std::string &s,
				     HashIntoType lower_bound = 0,
				     HashIntoType upper_bound = 0);
				     
    BoundedCounterType get_max_count(const std::string &s,
				     HashIntoType lower_bound = 0,
				     HashIntoType upper_bound = 0);

    HashIntoType * abundance_distribution() const;
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
