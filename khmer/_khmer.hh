/*
This file is part of khmer, https://github.com/dib-lab/khmer/, and is
Copyright (C) 2010-2015, Michigan State University.
Copyright (C) 2015-2016, The Regents of the University of California.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are
met:

    * Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.

    * Redistributions in binary form must reproduce the above
      copyright notice, this list of conditions and the following
      disclaimer in the documentation and/or other materials provided
      with the distribution.

    * Neither the name of the Michigan State University nor the names
      of its contributors may be used to endorse or promote products
      derived from this software without specific prior written
      permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
"AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
LICENSE (END)

Contact: khmer-project@idyll.org
*/

//
// A module for Python that exports khmer C++ library functions.
//
#ifndef _KHMER_HH
#define _KHMER_HH

#include <Python.h>

#include <iostream>

#include "khmer.hh"
#include "kmer_hash.hh"
#include "hashtable.hh"
#include "hashgraph.hh"
#include "assembler.hh"
#include "read_aligner.hh"
#include "labelhash.hh"
#include "khmer_exception.hh"
#include "hllcounter.hh"


//
// Python 2/3 compatibility: PyInt and PyLong
//

#if (PY_MAJOR_VERSION >= 3)
#define PyInt_Check(arg) PyLong_Check(arg)
#define PyInt_AsLong(arg) PyLong_AsLong(arg)
#define PyInt_FromLong(arg) PyLong_FromLong(arg)
#define Py_TPFLAGS_HAVE_ITER 0
#endif

//
// Python 2/3 compatibility: PyBytes and PyString
// https://docs.python.org/2/howto/cporting.html#str-unicode-unification
//

#include "bytesobject.h"

//
// Python 2/3 compatibility: Module initialization
// http://python3porting.com/cextensions.html#module-initialization
//

#if PY_MAJOR_VERSION >= 3
#define MOD_ERROR_VAL NULL
#define MOD_SUCCESS_VAL(val) val
#define MOD_INIT(name) PyMODINIT_FUNC PyInit_##name(void)
#define MOD_DEF(ob, name, doc, methods) \
          static struct PyModuleDef moduledef = { \
            PyModuleDef_HEAD_INIT, name, doc, -1, methods, }; \
          ob = PyModule_Create(&moduledef);
#else
#define MOD_ERROR_VAL
#define MOD_SUCCESS_VAL(val)
#define MOD_INIT(name) void init##name(void)
#define MOD_DEF(ob, name, doc, methods) \
          ob = Py_InitModule3(name, methods, doc);
#endif

namespace khmer {
namespace python {


typedef struct {
    PyObject_HEAD
    //! Pointer to the low-level genomic read object.
    read_parsers:: Read *   read;
} khmer_Read_Object;


typedef struct {
    PyObject_HEAD
    //! Pointer to the low-level parser object.
    read_parsers:: IParser *  parser;
} khmer_ReadParser_Object;


typedef struct {
    PyObject_HEAD
    //! Pointer to Python parser object for reference counting purposes.
    PyObject *  parent;
    //! Persistent value of pair mode across invocations.
    int pair_mode;
} khmer_ReadPairIterator_Object;

} //python

typedef struct {
    PyObject_HEAD
    pre_partition_info *   PrePartitionInfo;
} khmer_PrePartitionInfo_Object;


typedef struct {
    PyObject_HEAD
    SeenSet * hashes;
    WordLength ksize;
} khmer_HashSet_Object;


typedef struct {
    PyObject_HEAD
    khmer_HashSet_Object * parent;
    SeenSet::iterator * it;
} _HashSet_iterobj;


typedef struct {
    PyObject_HEAD
    Hashtable * hashtable;
} khmer_KHashtable_Object;


typedef struct {
    khmer_KHashtable_Object khashtable;
    Hashgraph * hashgraph;
} khmer_KHashgraph_Object;


typedef struct {
    khmer_KHashgraph_Object khashgraph;
    Nodegraph * nodegraph;
} khmer_KNodegraph_Object;


typedef struct {
    khmer_KHashgraph_Object khashgraph;
    Countgraph * countgraph;
} khmer_KCountgraph_Object;

typedef struct {
    khmer_KHashtable_Object khashtable;
    Counttable * counttable;
} khmer_KCounttable_Object;


typedef struct {
    khmer_KHashtable_Object khashtable;
    Nodetable * nodetable;
} khmer_KNodetable_Object;


typedef struct {
    PyObject_HEAD
    SubsetPartition * subset;
} khmer_KSubsetPartition_Object;


typedef struct {
    PyObject_HEAD
    ReadAligner * aligner;
} khmer_ReadAligner_Object;


typedef struct {
    PyObject_HEAD
    LabelHash * labelhash;
} khmer_KGraphLabels_Object;


typedef struct {
    PyObject_HEAD
    HLLCounter * hllcounter;
} khmer_KHLLCounter_Object;


typedef struct {
    PyObject_HEAD
    LinearAssembler * assembler;
} khmer_KLinearAssembler_Object;


typedef struct {
    PyObject_HEAD
    SimpleLabeledAssembler * assembler;
} khmer_KSimpleLabeledAssembler_Object;


typedef struct {
    PyObject_HEAD
    JunctionCountAssembler * assembler;
} khmer_KJunctionCountAssembler_Object;


static bool convert_Pytablesizes_to_vector(PyListObject * sizes_list_o,
                                           std::vector<uint64_t>& sizes)
{
    Py_ssize_t sizes_list_o_length = PyList_GET_SIZE(sizes_list_o);
    if (sizes_list_o_length < 1) {
        PyErr_SetString(PyExc_ValueError,
                        "tablesizes needs to be one or more numbers");
        return false;
    }
    for (Py_ssize_t i = 0; i < sizes_list_o_length; i++) {
        PyObject * size_o = PyList_GET_ITEM(sizes_list_o, i);
        if (PyLong_Check(size_o)) {
            sizes.push_back(PyLong_AsUnsignedLongLong(size_o));
        } else if (PyInt_Check(size_o)) {
            sizes.push_back(PyInt_AsLong(size_o));
        } else if (PyFloat_Check(size_o)) {
            sizes.push_back(PyFloat_AS_DOUBLE(size_o));
        } else {
            PyErr_SetString(PyExc_TypeError,
                            "2nd argument must be a list of ints, longs, or floats");
            return false;
        }
    }
    return true;
}
}

#endif
