//
// This file is part of khmer, http://github.com/ged-lab/khmer/, and is
// Copyright (C) Michigan State University, 2009-2013. It is licensed under
// the three-clause BSD license; see doc/LICENSE.txt. Contact: ctb@msu.edu
//
#ifndef STORAGE_HH
#define STORAGE_HH

#include <fstream>
#include "khmer.h"
#include "ktable.h"
#include "hashtable.h"

namespace khmer {
  class ReadMaskTable;
  unsigned int output_filtered_fasta_file(const std::string &inputfile,
					  const std::string &outputfile,
					  ReadMaskTable * readmask,
					  CallbackFn callback = NULL,
					  void * callback_data = NULL);

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

    const unsigned int get_tablesize() const {
      return _tablesize;
    }

    const unsigned int n_kept() const {
      unsigned int n = 0;
      for (unsigned int i = 0; i < _tablesize; i++) {
	if (_mask[i]) { n++; }
      }
      return n;
    }

    const bool get(unsigned int index) const {
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

    void invert() {
      for (unsigned int i = 0; i < _tablesize; i++) {
	_mask[i] = !_mask[i];
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

    unsigned int filter_fasta_file(const std::string &inputfile,
				   const std::string &outputfile,
				   CallbackFn callback = NULL,
				   void * callback_data = NULL) {
      return output_filtered_fasta_file(inputfile, outputfile, this,
					callback, callback_data);
    }
  };

  // ** //

  typedef struct {
    BoundedCounterType min_val;
    BoundedCounterType max_val;
  } MinMaxValue;

  class MinMaxTable {
  protected:
    const unsigned int _tablesize;
    MinMaxValue * _table;

    void _allocate(bool initialize=true) {
      _table = new MinMaxValue[_tablesize];
      if (initialize) {
	memset(_table, 0, sizeof(MinMaxValue) * _tablesize);
      }
    }

  public:
    MinMaxTable(unsigned int tablesize) :
      _tablesize(tablesize), _table(NULL) {
      _allocate();
    }

    ~MinMaxTable() {
      if (_table) { delete _table; _table = NULL; }
    }

    const unsigned int get_tablesize() const {
      return _tablesize;
    }

    const BoundedCounterType get_min(unsigned int index) const {
      assert(index < _tablesize);
      return _table[index].min_val;
    }

    const BoundedCounterType get_max(unsigned int index) const {
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

    void save(const std::string &outputfile) const {
      std::ofstream outfile;
      outfile.open(outputfile.c_str(), std::ofstream::binary);

      outfile.write((const char *) &_tablesize, sizeof(_tablesize));
      outfile.write((const char *) _table, _tablesize * sizeof(MinMaxValue));
      outfile.close();
    }

    void load(const std::string &inputfile) {
      std::ifstream infile;
      infile.open(inputfile.c_str(), std::ifstream::binary);

      infile.read((char *) &_tablesize, sizeof(_tablesize));
      _allocate(false);
      infile.read((char *) _table, _tablesize * sizeof(MinMaxValue));

      infile.close();
    }
  };

}

#endif // STORAGE_HH
