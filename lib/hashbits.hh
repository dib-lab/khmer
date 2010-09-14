#ifndef HASHBITS_HH
#define HASHBITS_HH

#include "hashtable.hh"

namespace khmer {
  class Hashbits : public Hashtable {
  protected:
    BoundedCounterType * _counts[N_TABLES];
    HashIntoType _tablebytes;

    virtual void _allocate_counters() {
      for (unsigned int i = 0; i < N_TABLES; i++) {
	_tablesize = primes[i];
	_tablebytes = _tablesize / 8 + 1;

	_counts[i] = new BoundedCounterType[_tablebytes];
	memset(_counts[i], 0, _tablebytes * sizeof(BoundedCounterType));
      }
    }
  public:
    Hashbits(WordLength ksize, HashIntoType tablesize) :
      Hashtable(ksize, tablesize) { }

    ~Hashbits() {
      if (_counts[0]) {
	for (unsigned int i = 0; i < N_TABLES; i++) {
	  delete _counts[i];
	  _counts[i] = NULL;
	}
      }
      _clear_partitions();
    }

    virtual void save(std::string);
    virtual void load(std::string);

    // count number of occupied bins
    virtual const HashIntoType n_occupied(HashIntoType start=0,
				  HashIntoType stop=0) const {
      HashIntoType n = 0;
      if (stop == 0) { stop = _tablesize; }
      for (HashIntoType i = start; i < stop; i++) {
	unsigned int byte = i / 8;
	unsigned char bit = i % 8;
	if (_counts[0][byte] & (1 << bit)) {
	  n++;
	}
      }
      return n;
    }

    virtual void count(const char * kmer) {
      HashIntoType hash = _hash(kmer, _ksize);

      for (unsigned int i = 0; i < N_TABLES; i++) {
	HashIntoType bin = hash % primes[i];
	unsigned int byte = bin / 8;
	unsigned char bit = bin % 8;

	_counts[i][byte] |= (1 << bit);
      }
    }

    virtual void count(HashIntoType khash) {
      for (unsigned int i = 0; i < N_TABLES; i++) {
	HashIntoType bin = khash % primes[i];
	unsigned int byte = bin / 8;
	unsigned char bit = bin % 8;

	_counts[i][byte] |= (1 << bit);
      }
    }

    // get the count for the given k-mer.
    virtual const BoundedCounterType get_count(const char * kmer) const {
      HashIntoType hash = _hash(kmer, _ksize);

      for (unsigned int i = 0; i < N_TABLES; i++) {
	HashIntoType bin = hash % primes[i];
	unsigned int byte = bin / 8;
	unsigned char bit = bin % 8;
      
	if (!(_counts[i][byte] & (1 << bit))) {
	  return 0;
	}
      }
      return 1;
    }

    // get the count for the given k-mer hash.
    virtual const BoundedCounterType get_count(HashIntoType khash) const {
      for (unsigned int i = 0; i < N_TABLES; i++) {
	HashIntoType bin = khash % primes[i];
	unsigned int byte = bin / 8;
	unsigned char bit = bin % 8;
      
	if (!(_counts[i][byte] & (1 << bit))) {
	  return 0;
	}
      }
      return 1;
    }
  };

};

#endif // HASHBITS_HH
