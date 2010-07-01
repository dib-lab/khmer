#ifndef STORAGE_HH
#define STORAGE_HH

#include <fstream>
#include "khmer.hh"
#include "hashtable.hh"

namespace khmer {
  class ReadMaskTable {
  protected:
    const unsigned int _tablesize;
    BoundedCounterType * _mask;

    void _allocate() {
      _mask = new BoundedCounterType[_tablesize];
      for (unsigned int i = 0; i < _tablesize; i++) {
	_mask[i] = 1;
      }
    }

  public:
    ReadMaskTable(unsigned int tablesize) :
      _tablesize(tablesize), _mask(NULL) {
      _allocate();
    }

    ~ReadMaskTable() {
      delete _mask;
    }

    const bool get(unsigned int index) {
      assert(index < _tablesize);
      return _mask[index];
    }

    void set(unsigned int index, bool keep) {
      assert(index < _tablesize);
      _mask[index] = keep ? 1 : 0;
    }
  };
}

#endif // STORAGE_HH
