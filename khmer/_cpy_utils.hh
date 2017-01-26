#ifndef _CPY_UTILS_HH
#define _CPY_UTILS_HH


#include <Python.h>
#include <vector>
#include "khmer.hh"
#include "kmer_hash.hh"
#include "hashtable.hh"
#include "khmer_exception.hh"

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


namespace khmer {


bool 
convert_HashIntoType_to_PyObject(const HashIntoType &hashval, 
                                 PyObject **value);


bool 
convert_PyLong_to_HashIntoType(PyObject * value, HashIntoType &hashval);


bool 
convert_PyObject_to_HashIntoType(PyObject * value, 
                                 HashIntoType& hashval, 
                                 WordLength ksize);


bool 
ht_convert_PyObject_to_HashIntoType(PyObject * value, 
                                    HashIntoType& hashval, 
                                    const Hashtable * ht);


bool 
ht_convert_PyObject_to_Kmer(PyObject * value, 
                            Kmer& kmer, 
                            const Hashtable * ht);


bool convert_Pytablesizes_to_vector(PyListObject * sizes_list_o,
                                           std::vector<uint64_t>& sizes);

}

#endif
