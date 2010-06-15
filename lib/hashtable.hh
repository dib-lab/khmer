#ifndef HASHTABLE_HH
#define HASHTABLE_HH

#include "khmer.hh"

namespace khmer {
  typedef unsigned char HashcountType;
  class Hashtable {
  protected:
    const unsigned int _tablesize;
    const unsigned int _ksize;

    HashcountType * _counts;

    void _allocate_counters() {
      _counts = new HashcountType[_tablesize];
      memset(_counts, 0, _tablesize * sizeof(HashcountType));
    }

  public:
    Hashtable(unsigned int ksize, unsigned int tablesize) :
      _ksize(ksize), _tablesize(tablesize) {
      _allocate_counters();
    }

    ~Hashtable() {
      delete _counts; _counts = NULL;
    }

    void count(const char * kmer) {
      unsigned int bin = _hash(kmer, _ksize) % _tablesize;
      _counts[bin]++;
    }

    void count(unsigned int khash) {
      unsigned int bin = khash % _tablesize;
      _counts[bin]++;
    }

    // get the count for the given k-mer.
    const unsigned int get_count(const char * kmer) const {
      unsigned int bin = _hash(kmer, _ksize) % _tablesize;
      return _counts[bin];
    }

    // get the count for the given k-mer hash.
    const unsigned int get_count(unsigned int khash) const {
      unsigned int bin = khash % _tablesize;
      return _counts[bin];
    }

    // count every k-mer in the string.
    void consume_string(const std::string &s);
  };
};

#endif // HASHTABLE_HH
