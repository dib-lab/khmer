#ifndef HASHBITS_HH
#define HASHBITS_HH

#include "hashtable.hh"

namespace khmer {
  class Hashbits : public Hashtable {
  protected:
    // BoundedCounterType * _counts[N_TABLES];

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
  };

};

#endif // HASHBITS_HH
