#ifndef _MINHASH_HH
#define _MINHASH_HH

//
// Python 2/3 compatibility: PyInt and PyLong
//

// Must be first.
#include <Python.h>

#if (PY_MAJOR_VERSION >= 3)
#define PyInt_Check(arg) PyLong_Check(arg)
#define PyInt_AsLong(arg) PyLong_AsLong(arg)
#define PyInt_FromLong(arg) PyLong_FromLong(arg)
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

#include "kmer_min_hash.hh"

typedef struct {
    PyObject_HEAD
    khmer::KmerMinHash * mh;
} MinHash_Object;

khmer::KmerMinHash * extract_KmerMinHash(PyObject * mh_obj)
{
    //if (!PyObject_TypeCheck(mh_obj, &MinHash_Type)) {
    //  return NULL;
    //}
    MinHash_Object * obj = (MinHash_Object *) mh_obj;
    return obj->mh;
}

typedef struct {
    PyObject_HEAD
    khmer::NeighborhoodMinHash * nbhd_mh;
} NeighborhoodMinHash_Object;

khmer::NeighborhoodMinHash * extract_NeighborhoodMinHash(PyObject * nbhd_obj)
{
    //if (!PyObject_TypeCheck(nbhd_obj, &NeighborhoodMinHash_Type)) {
    //  return NULL;
    //}
    NeighborhoodMinHash_Object * obj = (NeighborhoodMinHash_Object *) nbhd_obj;
    return obj->nbhd_mh;
}

typedef struct {
    PyObject_HEAD
    khmer::CombinedMinHash * combined_mh;
} CombinedMinHash_Object;

khmer::CombinedMinHash * extract_CombinedMinHash(PyObject * combined_obj)
{
    //if (!PyObject_TypeCheck(combined_obj, &CombinedMinHash_Type)) {
    //  return NULL;
    //}
    CombinedMinHash_Object * obj = (CombinedMinHash_Object *) combined_obj;
    return obj->combined_mh;
}


#endif // _MINHASH_HH
