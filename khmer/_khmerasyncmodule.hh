//
// This file is part of khmer, http://github.com/ged-lab/khmer/, and is
// Copyright (C) Michigan State University, 2009-2013. It is licensed under
// the three-clause BSD license; see doc/LICENSE.txt.
// Contact: khmer-project@idyll.org
//

//
// A module for Python that exports khmer C++ library functions.
//

// Must be first.
#include <Python.h>
#include "_khmermodule.hh"

#include "khmer_async.hh"

using namespace khmer;


typedef struct {
    PyObject_HEAD
    AsyncSequenceProcessor * async_sp;
} khmer_AsyncSequenceProcessorObject;

static void khmer_asyncseqproc_dealloc(PyObject *);
static int khmer_asyncseqproc_init(khmer_AsyncSequenceProcessorObject * self, 
        PyObject * args, PyObject * kwds);
static PyObject * khmer_asyncseqproc_new(PyTypeObject * type, 
        PyObject * args, PyObject * kwds);
static PyObject * asyncseqproc_stop(PyObject * self, PyObject * args);
static PyObject * asyncseqproc_processed_iter(PyObject * self);
static PyObject * asyncseqproc_processed_iternext(PyObject * self);
static PyObject * asyncseqproc_n_parsed(PyObject * self, PyObject * args);
static PyObject * asyncseqproc_n_processed(PyObject * self, PyObject * args);


static PyMethodDef khmer_asyncseqproc_methods[] = {
    { "processed", (PyCFunction)asyncseqproc_processed_iter, METH_NOARGS, "Iterator over processed reads" },
    { "stop", asyncseqproc_stop, METH_VARARGS, "Stop processors, join threads" },
    { "n_processed", asyncseqproc_n_processed, METH_NOARGS, "Number of reads processed" },
    { "n_parsed", asyncseqproc_n_parsed, METH_NOARGS, "Number of reads parsed" },

    {NULL, NULL, 0, NULL}           /* sentinel */
};

static PyTypeObject khmer_AsyncSequenceProcessorType = {
    PyObject_HEAD_INIT(NULL)
    0,                       /* ob_size */
    "AsyncSequenceProcessor",            /* tp_name */
    sizeof(khmer_AsyncSequenceProcessorObject), /* tp_basicsize */
    0,                       /* tp_itemsize */
    (destructor)khmer_asyncseqproc_dealloc, /* tp_dealloc */
    0,                       /* tp_print */
    0,                       /* tp_getattr */
    0,                       /* tp_setattr */
    0,                       /* tp_compare */
    0,                       /* tp_repr */
    0,                       /* tp_as_number */
    0,                       /* tp_as_sequence */
    0,                       /* tp_as_mapping */
    0,                       /* tp_hash */
    0,                       /* tp_call */
    0,                       /* tp_str */
    0,                       /* tp_getattro */
    0,                       /* tp_setattro */
    0,                       /* tp_as_buffer */
    Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE | Py_TPFLAGS_HAVE_ITER,   /* tp_flags */
    0,                       /* tp_doc */
    0,                       /* tp_traverse */
    0,                       /* tp_clear */
    0,                       /* tp_richcompare */
    0,                       /* tp_weaklistoffset */
    PyObject_SelfIter,                       /* tp_iter */
    (iternextfunc) asyncseqproc_processed_iternext,                       /* tp_iternext */
    khmer_asyncseqproc_methods, /* tp_methods */
    0,                       /* tp_members */
    0,                       /* tp_getset */
    0,                       /* tp_base */
    0,                       /* tp_dict */
    0,                       /* tp_descr_get */
    0,                       /* tp_descr_set */
    0,                       /* tp_dictoffset */
    (initproc)khmer_asyncseqproc_init,   /* tp_init */
    0,                       /* tp_alloc */
};

static void khmer_asyncseqproc_dealloc(PyObject* obj);

static PyObject * khmer_asyncseqproc_new(PyTypeObject *type, 
        PyObject *args, PyObject *kwds);

static int khmer_asyncseqproc_init(khmer_AsyncSequenceProcessorObject * self, 
        PyObject *args, PyObject *kwds);
////////////////////
// AsyncDiginorm
////////////////////

typedef struct {
    khmer_AsyncSequenceProcessorObject async_sp;
    AsyncDiginorm * async_diginorm;
} khmer_AsyncDiginormObject;

static void khmer_asyncdiginorm_dealloc(PyObject *);
static int khmer_asyncdiginorm_init(khmer_AsyncDiginormObject * self, 
        PyObject * args, PyObject * kwds);
static PyObject * khmer_asyncdiginorm_new(PyTypeObject * type, 
        PyObject * args, PyObject * kwds);
static PyObject * asyncdiginorm_start(PyObject * self, PyObject * args);
static PyObject * asyncdiginorm_n_kept(PyObject * self, PyObject * args);

static PyMethodDef khmer_asyncdiginorm_methods[] = {
    { "start", asyncdiginorm_start, METH_VARARGS, "Start processing" },
    { "n_kept", asyncdiginorm_n_kept, METH_NOARGS, "Number of reads kept" },

    {NULL, NULL, 0, NULL}           /* sentinel */
};

static PyTypeObject khmer_AsyncDiginormType = {
    PyObject_HEAD_INIT(NULL)
    0,                       /* ob_size */
    "AsyncDiginorm",            /* tp_name */
    sizeof(khmer_AsyncDiginormObject), /* tp_basicsize */
    0,                       /* tp_itemsize */
    (destructor)khmer_asyncdiginorm_dealloc, /* tp_dealloc */
    0,                       /* tp_print */
    0,  /* khmer_labelhash_getattr, tp_getattr */
    0,                       /* tp_setattr */
    0,                       /* tp_compare */
    0,                       /* tp_repr */
    0,                       /* tp_as_number */
    0,                       /* tp_as_sequence */
    0,                       /* tp_as_mapping */
    0,                       /* tp_hash */
    0,                       /* tp_call */
    0,                       /* tp_str */
    0,                       /* tp_getattro */
    0,                       /* tp_setattro */
    0,                       /* tp_as_buffer */
    Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE | Py_TPFLAGS_HAVE_ITER,   /* tp_flags */
    0,                       /* tp_doc */
    0,                       /* tp_traverse */
    0,                       /* tp_clear */
    0,                       /* tp_richcompare */
    0,                       /* tp_weaklistoffset */
    PyObject_SelfIter,                       /* tp_iter */
    (iternextfunc) asyncseqproc_processed_iternext,                       /* tp_iternext */
    khmer_asyncseqproc_methods, /* tp_methods */
    0,                       /* tp_members */
    0,                       /* tp_getset */
    0,                       /* tp_base */
    0,                       /* tp_dict */
    0,                       /* tp_descr_get */
    0,                       /* tp_descr_set */
    0,                       /* tp_dictoffset */
    (initproc)khmer_asyncdiginorm_init,   /* tp_init */
    0,                       /* tp_alloc */
};

static void khmer_asyncdiginorm_dealloc(PyObject* obj);

static PyObject * khmer_asyncdiginorm_new(PyTypeObject *type, 
        PyObject *args, PyObject *kwds);

static int khmer_asyncdiginorm_init(khmer_AsyncDiginormObject * self,
        PyObject *args, PyObject *kwds);

