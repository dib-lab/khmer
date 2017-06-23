#ifndef _CPY_UTILS_HH
#define _CPY_UTILS_HH


#include <Python.h>
#include <vector>
#include "oxli/oxli.hh"
#include "oxli/kmer_hash.hh"
#include "oxli/hashtable.hh"
#include "oxli/oxli_exception.hh"

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
convert_HashIntoType_to_PyObject(const oxli::HashIntoType &hashval, 
                                 PyObject **value);


bool 
convert_PyLong_to_HashIntoType(PyObject * value,
                               oxli::HashIntoType &hashval);


bool 
convert_PyObject_to_HashIntoType(PyObject * value, 
                                 oxli::HashIntoType& hashval, 
                                 oxli::WordLength ksize);


bool 
ht_convert_PyObject_to_HashIntoType(PyObject * value, 
                                    oxli::HashIntoType& hashval, 
                                    const oxli::Hashtable * ht);


bool 
ht_convert_PyObject_to_Kmer(PyObject * value, 
                            oxli::Kmer& kmer, 
                            const oxli::Hashtable * ht);


bool convert_Pytablesizes_to_vector(PyListObject * sizes_list_o,
                                    std::vector<uint64_t>& sizes);

}

#endif
