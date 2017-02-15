
#ifndef _CPY_HASHSET_HH
#define _CPY_HASHSET_HH

#include <Python.h>
#include "khmer.hh"
#include "_cpy_utils.hh"

namespace khmer {

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

extern PyMethodDef khmer_HashSet_methods[]; 

extern PySequenceMethods khmer_HashSet_seqmethods[];

extern PyTypeObject khmer_HashSet_Type;

extern PyTypeObject _HashSet_iter_Type;


khmer_HashSet_Object * create_HashSet_Object(SeenSet * h, WordLength k);

void
khmer_HashSet_dealloc(khmer_HashSet_Object * obj);

PyObject* khmer_HashSet_new(PyTypeObject * type, PyObject * args,
                                   PyObject * kwds);

void _HashSet_iter_dealloc(_HashSet_iterobj * obj);


PyObject * _HashSet_iter(PyObject * self);

PyObject * _HashSet_iternext(PyObject * self);

PyObject * khmer_HashSet_iter(PyObject * self);


int khmer_HashSet_len(khmer_HashSet_Object * o);

PyObject * khmer_HashSet_concat(khmer_HashSet_Object * o,
                                       khmer_HashSet_Object * o2);

PyObject * khmer_HashSet_concat_inplace(khmer_HashSet_Object * o,
        khmer_HashSet_Object * o2);

int khmer_HashSet_contains(khmer_HashSet_Object * o, PyObject * val);

PyObject *
hashset_add(khmer_HashSet_Object * me, PyObject * args);

PyObject *
hashset_remove(khmer_HashSet_Object * me, PyObject * args);

PyObject *
hashset_update(khmer_HashSet_Object * me, PyObject * args);
}

#endif
