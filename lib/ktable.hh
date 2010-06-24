#ifndef KTABLE_HH
#define KTABLE_HH

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <assert.h>

// bit representation of A/T/C/G.
#define twobit_repr(ch) ((toupper(ch)) == 'A' ? 0 : \
                         (toupper(ch)) == 'T' ? 1 : \
                         (toupper(ch)) == 'C' ? 2 : 3)

#define revtwobit_repr(n) ((n) == 0 ? 'A' : \
                           (n) == 1 ? 'T' : \
                           (n) == 2 ? 'C' : 'G')

#define twobit_comp(ch) ((toupper(ch)) == 'A' ? 1 : \
                         (toupper(ch)) == 'T' ? 0 : \
                         (toupper(ch)) == 'C' ? 3 : 2)


namespace khmer {
  typedef long long CounterType;

  // two-way hash functions.
  unsigned long long int _hash(const char * kmer, unsigned int k);
  unsigned long long int _hash(const char * kmer, unsigned int k,
                               unsigned long long int * h, unsigned long long int * r);
  std::string _revhash(unsigned int hash, unsigned int k);

  //
  // KTable class: keep track of k-mer prevalences.
  //
  // the main (so far only...) class in khmer.
  //

  class KTable {
  protected:
    const unsigned int _ksize;	// 'k'
    unsigned long long int _max_hash;	// 4**k

    CounterType * _counts;	// counts table.

    // allocate the counts table.
    void _allocate_counters() {
      // allocate.
      _counts = new CounterType[n_entries()];
      memset(_counts, 0, n_entries() * sizeof(CounterType));
    }
  public:

    // Constructor: initialize stuff.
    KTable(unsigned int ksize) : _ksize(ksize) {
      _max_hash = (unsigned int) pow(4, _ksize) - 1;
      _allocate_counters();
    }

    // destructor: free.
    ~KTable() {
      delete _counts; _counts = NULL;
    }

    // accessor to get 'k'
    const unsigned int ksize() const { return _ksize; }

    // accessors to get table info
    const unsigned long long int max_hash() const { return _max_hash; }
    const unsigned long long int n_entries() const { return _max_hash + 1; }

    // add the given k-mer into the counts table.
    void count(const char * kmer) {
      _counts[_hash(kmer, _ksize)]++;
    }

    // get the count for the given k-mer.
    const unsigned long long int get_count(const char * kmer) const {
      return _counts[_hash(kmer, _ksize)];
    }

    // get the count for the given k-mer hash.
    const unsigned long long int get_count(unsigned long long int i) const {
      return _counts[i];
    }

    // set the count for the given k-mer.
    void set_count(const char * kmer, unsigned int c) {
      _counts[_hash(kmer, _ksize)] = c;
    }

    // set the count for the given k-mer hash.
    void set_count(unsigned int i, unsigned int c) {
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
}

#endif // KTABLE_HH
