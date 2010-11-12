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

  class Hashtable {		// Base class implementation of a Bloom ht.
  protected:
    WordLength _ksize;
    HashIntoType bitmask;

    Hashtable(WordLength ksize) : _ksize(ksize) {
      bitmask = 0;
      for (unsigned int i = 0; i < _ksize; i++) {
	bitmask = (bitmask << 2) | 3;
      }
    }

    virtual ~Hashtable() {}

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
  };

};

#endif // HASHTABLE_HH
