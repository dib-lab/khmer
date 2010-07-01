#ifndef STORAGE_HH
#define STORAGE_HH

#include <fstream>
#include "khmer.hh"
#include "hashtable.hh"

namespace khmer {
  class ReadMaskTable {
  protected:
    const unsigned int _tablesize;
    unsigned char * _mask;

    void _allocate() {
      _mask = new unsigned char[_tablesize];
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
      if (_mask) { delete _mask; _mask = NULL; }
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

  // ** //

  typedef struct {
    BoundedCounterType min_val;
    BoundedCounterType max_val;
  } MinMaxValue;;

  class MinMaxTable {
  protected:
    const unsigned int _tablesize;
    MinMaxValue * _table;

    void _allocate() {
      _table = new MinMaxValue[_tablesize];
      memset(_table, 0, sizeof(MinMaxValue) * _tablesize);
    }

  public:
    MinMaxTable(unsigned int tablesize) :
      _tablesize(tablesize), _table(NULL) {
      _allocate();
    }

    ~MinMaxTable() {
      if (_table) { delete _table; _table = NULL; }
    }

    const BoundedCounterType get_min(unsigned int index) {
      assert(index < _tablesize);
      return _table[index].min_val;
    }

    const BoundedCounterType get_max(unsigned int index) {
      assert(index < _tablesize);
      return _table[index].max_val;
    }

    BoundedCounterType add_max(unsigned int index, unsigned int val) {
      assert(index < _tablesize);

      if (val > MAX_COUNT) { val = MAX_COUNT; }

      MinMaxValue * mmv = &_table[index];

      if (val > mmv->max_val) {
	mmv->max_val = (unsigned char) val;
      }

      return mmv->max_val;
    }

    BoundedCounterType add_min(unsigned int index, unsigned int val) {
      assert(index < _tablesize);

      if (val > MAX_COUNT) { val = MAX_COUNT; }

      MinMaxValue * mmv = &_table[index];

      if (mmv->min_val == 0) {
	mmv->min_val = (unsigned char) val;
      } else if (val < mmv->min_val) {
	mmv->min_val = (unsigned char) val;
      }

      return mmv->min_val;
    }

    void clear(unsigned int index) {
      assert(index < _tablesize);

      MinMaxValue * mmv = &_table[index];
      mmv->min_val = mmv->max_val = 0;
    }
  };

}

#endif // STORAGE_HH
