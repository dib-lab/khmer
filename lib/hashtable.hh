#ifndef HASHTABLE_HH
#define HASHTABLE_HH

#include <iostream>
#include <list>
#include <queue>

#include <fstream>
#include <string>
#include <set>
#include <map>
#include <queue>

#include "khmer.hh"
#include "storage.hh"

#define CALLBACK_PERIOD 10000

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
  typedef std::map<unsigned long long, unsigned long long> PartitionCountDistribution;

  //
  // Sequence iterator class, test.  Not really a C++ iterator yet.
  //

  class KMerIterator {
  protected:
    const char * _seq;
    const unsigned char _ksize;
    
    HashIntoType _kmer_f, _kmer_r;
    HashIntoType bitmask;
    unsigned int _nbits_sub_1;
    unsigned int index, length;
    bool initialized;
  public:
    KMerIterator(const char * seq, unsigned char k) : _seq(seq), _ksize(k) {
      bitmask = 0;
      for (unsigned int i = 0; i < _ksize; i++) {
	bitmask = (bitmask << 2) | 3;
      }
      _nbits_sub_1 = (_ksize*2 - 2);

      index = _ksize - 1;
      length = strlen(seq);
      initialized = false;
    }

    HashIntoType first(HashIntoType& f, HashIntoType& r) {
      HashIntoType x;
      x = _hash(_seq, _ksize, _kmer_f, _kmer_r);

      f = _kmer_f;
      r = _kmer_r;

      index = _ksize;

      return x;
    }

    HashIntoType next(HashIntoType& f, HashIntoType& r) {
      if (done()) {
	throw std::exception();
      }

      if (!initialized) {
	initialized = true;
	return first(f, r);
      }

      unsigned char ch = _seq[index];
      index++;
      assert(index <= length);

      // left-shift the previous hash over
      _kmer_f = _kmer_f << 2;

      // 'or' in the current nt
      _kmer_f |= twobit_repr(ch);

      // mask off the 2 bits we shifted over.
      _kmer_f &= bitmask;

      // now handle reverse complement
      _kmer_r = _kmer_r >> 2;
      _kmer_r |= (twobit_comp(ch) << _nbits_sub_1);

      f = _kmer_f;
      r = _kmer_r;

      return uniqify_rc(_kmer_f, _kmer_r);
    }

    HashIntoType first() { return first(_kmer_f, _kmer_r); }
    HashIntoType next() { return next(_kmer_f, _kmer_r); }

    bool done() { return index >= length; }
  };

  class ReadMaskTable;

  class Hashtable {		// Base class implementation of a Bloom ht.
  protected:
    WordLength _ksize;
    HashIntoType bitmask;
    unsigned int _nbits_sub_1;

    Hashtable(WordLength ksize) : _ksize(ksize) {
      _init_bitstuff();
    }

    virtual ~Hashtable() {}

    void _init_bitstuff() {
      bitmask = 0;
      for (unsigned int i = 0; i < _ksize; i++) {
	bitmask = (bitmask << 2) | 3;
      }
      _nbits_sub_1 = (_ksize*2 - 2);
    }

    HashIntoType _next_hash(char ch, HashIntoType &h, HashIntoType &r) const {
      // left-shift the previous hash over
      h = h << 2;

      // 'or' in the current nt
      h |= twobit_repr(ch);

      // mask off the 2 bits we shifted over.
      h &= bitmask;

      // now handle reverse complement
      r = r >> 2;
      r |= (twobit_comp(ch) << _nbits_sub_1);

      return uniqify_rc(h, r);
    }

  public:

    // accessor to get 'k'
    const WordLength ksize() const { return _ksize; }

    virtual void count(const char * kmer) = 0;
    virtual void count(HashIntoType khash) = 0;

    // get the count for the given k-mer.
    virtual const BoundedCounterType get_count(const char * kmer) const = 0;
    virtual const BoundedCounterType get_count(HashIntoType khash) const = 0;

    virtual void save(std::string) = 0;
    virtual void load(std::string) = 0;

    // count every k-mer in the string.
    unsigned int consume_string(const std::string &s,
				HashIntoType lower_bound = 0,
				HashIntoType upper_bound = 0);

    // checks each read for non-ACGT characters
    bool check_and_normalize_read(std::string &read) const;

    // check each read for non-ACGT characters, and then consume it.
    unsigned int check_and_process_read(std::string &read,
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
  };

			  

};

#endif // HASHTABLE_HH
