#ifndef _CPY_COUNTGRAPH_HH
#define _CPY_COUNTGRAPH_HH

#include <Python.h>
#include "_cpy_utils.hh"
#include "_cpy_hashgraph.hh"

namespace khmer {

typedef struct {
    khmer_KHashgraph_Object khashgraph;
    Countgraph * countgraph;
} khmer_KCountgraph_Object;

extern PyTypeObject khmer_KCountgraph_Type
CPYCHECKER_TYPE_OBJECT_FOR_TYPEDEF("khmer_KCountgraph_Object");

extern PyMethodDef khmer_countgraph_methods[];

PyObject* khmer_countgraph_new(PyTypeObject * type, PyObject * args,
                                      PyObject * kwds);

void khmer_countgraph_dealloc(khmer_KCountgraph_Object * obj);

PyObject *
count_get_raw_tables(khmer_KCountgraph_Object * self, PyObject * args);

PyObject *
count_do_subset_partition_with_abundance(khmer_KCountgraph_Object * me,
        PyObject * args);

}

#endif
