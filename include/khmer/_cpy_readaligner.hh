#ifndef _CPY_READALIGNER_HH
#define _CPY_READALIGNER_HH

#include <Python.h>
#include "read_aligner.hh"

namespace khmer {

typedef struct {
    PyObject_HEAD
    ReadAligner * aligner;
} khmer_ReadAligner_Object;

extern PyTypeObject khmer_ReadAlignerType;

extern PyMethodDef khmer_ReadAligner_methods[];


PyObject* khmer_ReadAligner_new(PyTypeObject *type, PyObject * args,
                                       PyObject *kwds);


void khmer_readaligner_dealloc(khmer_ReadAligner_Object* obj);


PyObject * readaligner_align(khmer_ReadAligner_Object * me,
                                    PyObject * args);


PyObject * readaligner_align_forward(khmer_ReadAligner_Object * me,
        PyObject * args);


PyObject* khmer_ReadAligner_get_scoring_matrix(
    khmer_ReadAligner_Object * me, PyObject * args);


PyObject* khmer_ReadAligner_get_transition_probabilities(
    khmer_ReadAligner_Object * me, PyObject * args);

}

#endif
