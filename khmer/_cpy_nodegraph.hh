#ifndef _CPY_NODEGRAPH_HH
#define _CPY_NODEGRAPH_HH

#include <Python.h>
#include "_cpy_utils.hh"
#include "_cpy_hashgraph.hh"
#include "hashgraph.hh"

namespace khmer {

typedef struct {
    khmer_KHashgraph_Object khashgraph;
    Nodegraph * nodegraph;
} khmer_KNodegraph_Object;

void khmer_nodegraph_dealloc(khmer_KNodegraph_Object * obj);

// __new__ for nodegraph; necessary for proper subclassing
// This will essentially do what the old factory function did. Unlike many __new__
// methods, we take our arguments here, because there's no "uninitialized" nodegraph
// object; we have to have k and the table sizes before creating the new objects
PyObject* khmer_nodegraph_new(PyTypeObject * type, PyObject * args,
                                    PyObject * kwds);

extern PyTypeObject khmer_KNodegraph_Type
CPYCHECKER_TYPE_OBJECT_FOR_TYPEDEF("khmer_KNodegraph_Object");

extern PyMethodDef khmer_nodegraph_methods[];

PyObject *
nodegraph_update(khmer_KNodegraph_Object * me, PyObject * args);


PyObject *
nodegraph_get_raw_tables(khmer_KNodegraph_Object * self, PyObject * args);



}

#endif
