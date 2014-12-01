//
// This file is part of khmer, http://github.com/ged-lab/khmer/, and is
// Copyright (C) Michigan State University, 2009-2013. It is licensed under
// the three-clause BSD license; see doc/LICENSE.txt.
// Contact: khmer-project@idyll.org
//

//
// A module for Python that exports khmer C++ library functions.
//

#ifndef KHMERMODULE_HH
#define KHMERMODULE_HH

#include <Python.h>
#include "hashtable.hh"

#include "khmer.hh"
#include "kmer_hash.hh"
#include "hashbits.hh"
#include "counting.hh"
#include "read_aligner.hh"
#include "labelhash.hh"

#include "async/_khmerasyncmodule.hh"
#include "khmer_exception.hh"

using namespace khmer;

class _khmer_exception
{
private:
    std::string _message;
public:
    _khmer_exception(std::string message) : _message(message) { };
    inline const std::string get_message() const
    {
        return _message;
    };
};

class _khmer_signal : public _khmer_exception
{
public:
    _khmer_signal(std::string message) : _khmer_exception(message) { };
};

typedef pre_partition_info _pre_partition_info;

extern "C" {

namespace khmer {
namespace python {

typedef struct {
    PyObject_HEAD
    //! Pointer to the low-level genomic read object.
    read_parsers:: Read *   read;
} Read_Object;

PyAPI_DATA(PyTypeObject) Read_Type;

typedef struct {
    PyObject_HEAD
    //! Pointer to the low-level parser object.
    read_parsers:: IParser *  parser;
} ReadParser_Object;

PyAPI_DATA(PyTypeObject) ReadParser_Type;

typedef struct {
    PyObject_HEAD
    //! Pointer to Python parser object for reference counting purposes.
    PyObject *  parent;
    //! Persistent value of pair mode across invocations.
    int pair_mode;
} ReadPairIterator_Object;

PyAPI_DATA(PyTypeObject) ReadPairIterator_Type;

} // namespace python
} // namespace khmer

typedef struct {
    PyObject_HEAD
    CountingHash * counting;
} khmer_KCountingHashObject;

PyAPI_DATA(PyTypeObject) khmer_KCountingHashType;

typedef struct {
    PyObject_HEAD
    SubsetPartition * subset;
} khmer_KSubsetPartitionObject;

PyAPI_DATA(PyTypeObject) khmer_KSubsetPartitionType;

typedef struct {
    PyObject_HEAD
    Hashbits * hashbits;
} khmer_KHashbitsObject;

PyAPI_DATA(PyTypeObject) khmer_KHashbitsType;

typedef struct {
    PyObject_HEAD
    ReadAligner * aligner;
} khmer_ReadAlignerObject;

PyAPI_DATA(PyTypeObject) khmer_ReadAlignerType;

typedef struct {
    khmer_KHashbitsObject khashbits;
    LabelHash * labelhash;
} khmer_KLabelHashObject;

PyAPI_DATA(PyTypeObject) KLabelHashType;

}

#endif  // KHMERMODULE_HH

// vim: set ft=cpp sts=4 sw=4 tw=79:
