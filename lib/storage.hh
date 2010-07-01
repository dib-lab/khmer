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

    void _allocate(bool initialize=true) {
      _mask = new unsigned char[_tablesize];
      if (initialize) {
	for (unsigned int i = 0; i < _tablesize; i++) {
	  _mask[i] = 1;
	}
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
      if (index >= _tablesize) { return false; } // @CTB throw?

      return _mask[index];
    }

    void set(unsigned int index, bool keep) {
      if (index >= _tablesize) { return; } // @CTB throw?

      _mask[index] = keep ? 1 : 0;
    }

    void merge(ReadMaskTable &other) {
      if (this ->_tablesize != other._tablesize) { return; } // @CTB throw?

      for (unsigned int i = 0; i < _tablesize; i++) {
	_mask[i] = _mask[i] && other._mask[i] ? 1 : 0;
      }
    }

    void save(const std::string &outputfile) {
      std::ofstream outfile;
      outfile.open(outputfile.c_str(), std::ofstream::binary);

      outfile.write((const char *) &_tablesize, sizeof(_tablesize));
      outfile.write((const char *) _mask, _tablesize);
      outfile.close();
    }

    void load(const std::string &inputfile) {
      std::ifstream infile;
      infile.open(inputfile.c_str(), std::ifstream::binary);

      infile.read((char *) &_tablesize, sizeof(_tablesize));
      _allocate(false);
      infile.read((char *) _mask, _tablesize);

      infile.close();
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

      if (val == 0) {
	;
      } else if (mmv->min_val == 0) {
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

    void merge(MinMaxTable &other) {
      assert(this->_tablesize == other._tablesize);
      for (unsigned int i = 0; i < _tablesize; i++) {
	this->add_min(i, other.get_min(i));
	this->add_max(i, other.get_max(i));
      }
    }
  };

}

#endif // STORAGE_HH
