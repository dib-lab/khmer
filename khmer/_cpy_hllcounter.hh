#ifndef _CPY_HLLCOUNTER_HH
#define _CPY_HLLCOUNTER_HH

#include <Python.h>
#include "_cpy_utils.hh"
#include "hllcounter.hh"

namespace khmer {


typedef struct {
    PyObject_HEAD
    HLLCounter * hllcounter;
} khmer_KHLLCounter_Object;


//
// KHLLCounter object
//

extern PyMethodDef khmer_hllcounter_methods[];

extern PyGetSetDef khmer_hllcounter_getseters[];

extern PyTypeObject khmer_KHLLCounter_Type;




PyObject* khmer_hllcounter_new(PyTypeObject * type, PyObject * args,
                                      PyObject * kwds);



void khmer_hllcounter_dealloc(khmer_KHLLCounter_Object * obj);



PyObject *
hllcounter_add(khmer_KHLLCounter_Object * me, PyObject * args);



PyObject *
hllcounter_estimate_cardinality(khmer_KHLLCounter_Object * me, PyObject * args);



PyObject *
hllcounter_consume_string(khmer_KHLLCounter_Object * me, PyObject * args);


PyObject * hllcounter_consume_fasta(khmer_KHLLCounter_Object * me,
        PyObject * args, PyObject * kwds);

PyObject * hllcounter_merge(khmer_KHLLCounter_Object * me,
                                   PyObject * args);


PyObject *
hllcounter_get_erate(khmer_KHLLCounter_Object * me);



PyObject *
hllcounter_get_ksize(khmer_KHLLCounter_Object * me);



int
hllcounter_set_ksize(khmer_KHLLCounter_Object * me, PyObject *value,
                     void *closure);


int
hllcounter_set_erate(khmer_KHLLCounter_Object * me, PyObject *value,
                     void *closure);


PyObject *
hllcounter_getalpha(khmer_KHLLCounter_Object * me);


PyObject *
hllcounter_getcounters(khmer_KHLLCounter_Object * me);


PyObject * hllcounter_merge(khmer_KHLLCounter_Object * me,
                                   PyObject * args);



}

#endif
