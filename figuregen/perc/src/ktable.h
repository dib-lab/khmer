#ifndef KTABLE_HH
#define KTABLE_HH

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <string>
#include <assert.h>

#include "khmer.h"

// test validity
#define is_valid_dna(ch) ((toupper(ch)) == 'A' || (toupper(ch)) == 'C' || \
			  (toupper(ch)) == 'G' || (toupper(ch)) == 'T')

// bit representation of A/T/C/G.
#define twobit_repr(ch) ((toupper(ch)) == 'A' ? 0LL : \
                         (toupper(ch)) == 'T' ? 1LL : \
                         (toupper(ch)) == 'C' ? 2LL : 3LL)

#define revtwobit_repr(n) ((n) == 0 ? 'A' : \
                           (n) == 1 ? 'T' : \
                           (n) == 2 ? 'C' : 'G')

#define twobit_comp(ch) ((toupper(ch)) == 'A' ? 1LL : \
                         (toupper(ch)) == 'T' ? 0LL : \
                         (toupper(ch)) == 'C' ? 3LL : 2LL)

// choose wisely between forward and rev comp.
#if !NO_UNIQUE_RC
#define uniqify_rc(f, r) ((f) < (r) ? (f) : (r))
#else
#define uniqify_rc(f,r)(f)
#endif


namespace khmer {
  // two-way hash functions.
  HashIntoType _hash(const char * kmer, const WordLength k);
  HashIntoType _hash(const char * kmer, const WordLength k,
		     HashIntoType& h, HashIntoType& r);
  HashIntoType _hash_forward(const char * kmer, WordLength k);

  std::string _revhash(HashIntoType hash, WordLength k);

  //
  // KTable class: keep track of k-mer prevalences.
  //
  // the main (so far only...) class in khmer.
  //

  class KTable {
  protected:
    WordLength _ksize;	// 'k'
    HashIntoType _max_hash;	// 4**k

    ExactCounterType * _counts;	// counts table.

    // allocate the counts table.
    void _allocate_counters() {
      // allocate.
      _counts = new ExactCounterType[n_entries()];
      memset(_counts, 0, n_entries() * sizeof(ExactCounterType));
    }
  public:

    // Constructor: initialize stuff.
    KTable(WordLength ksize) : _ksize(ksize) {
      _max_hash = (unsigned int) pow(4, _ksize) - 1;
      _allocate_counters();
    }

    // destructor: free.
    ~KTable() {
      delete _counts; _counts = NULL;
    }

    // accessor to get 'k'
    const WordLength ksize() const { return _ksize; }

    // accessors to get table info
    const HashIntoType max_hash() const { return _max_hash; }
    const HashIntoType n_entries() const { return _max_hash + 1; }

    // add the given k-mer into the counts table.
    void count(const char * kmer) {
      _counts[_hash(kmer, _ksize)]++;
    }

    // get the count for the given k-mer.
    const ExactCounterType get_count(const char * kmer) const {
      return _counts[_hash(kmer, _ksize)];
    }

    // get the count for the given k-mer hash.
    const ExactCounterType get_count(HashIntoType i) const {
      return _counts[i];
    }

    // set the count for the given k-mer.
    void set_count(const char * kmer, ExactCounterType c) {
      _counts[_hash(kmer, _ksize)] = c;
    }

    // set the count for the given k-mer hash.
    void set_count(HashIntoType i, ExactCounterType c) {
      assert(i <= max_hash());
      _counts[i] = c;
    }

    // count every k-mer in the string.
    void consume_string(const std::string &s);

    // reset the table.
    void clear() {
      delete _counts;
      _allocate_counters();
    }

    // set operations
    void update(const KTable &other);
    KTable * intersect(const KTable &other) const;
  };
};

#endif // KTABLE_HH
