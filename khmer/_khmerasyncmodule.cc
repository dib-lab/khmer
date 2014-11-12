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

#include "_khmerasyncmodule.hh"

using namespace khmer;


extern "C" {
    void init_khmer_async();
}

////////////////////
//
// AsyncSequenceProcessor
//
////////////////////

static PyObject * asyncseqproc_stop(PyObject * self, PyObject * args)
{
    khmer_AsyncSequenceProcessorObject * me = (khmer_AsyncSequenceProcessorObject *) self;
    AsyncSequenceProcessor * async_sp = me->async_sp;
    
    async_sp->stop();

    Py_RETURN_NONE;
}

static PyObject * asyncseqproc_processed_iter(PyObject * self) {
    return PyObject_SelfIter(self);
}

static PyObject * asyncseqproc_processed_iternext(PyObject * self) {
    khmer_AsyncSequenceProcessorObject * me = (khmer_AsyncSequenceProcessorObject *) self;
    AsyncSequenceProcessor * async_sp = me->async_sp;
   
    khmer::read_parsers::Read * read_ptr = new Read();

    while(!(async_sp->pop(read_ptr))) {
        if(async_sp->iter_stop()) {
            //async_sp->lock_stdout();
            //std::cout << "\nITER STOP: kept " << async_diginorm->n_kept() 
            //    << " of " << async_diginorm->n_processed()
            //    << " processed (" << async_diginorm->n_parsed()
            //    << " parsed)" << std::endl;
            //async_diginorm->unlock_stdout();
            delete read_ptr;
            PyErr_SetNone(PyExc_StopIteration);
            return NULL;
        }
    }
    /*
    if (async_diginorm->n_kept() % 50000 == 0)
        std::cout << "ITER: (kept, processed, written): " 
            << async_diginorm->n_kept() << ", " 
            << async_diginorm->n_processed() << ", "
            << async_diginorm->n_written() << std::endl;
    */
    PyObject * read_obj = khmer::python::Read_Type.tp_alloc(&khmer::python::Read_Type, 1);
    ((khmer::python::Read_Object *) read_obj)->read = read_ptr;
    return read_obj;
}

static PyObject * asyncseqproc_n_parsed(PyObject * self, PyObject * args)
{
    khmer_AsyncSequenceProcessorObject * me = (khmer_AsyncSequenceProcessorObject *) self;
    AsyncSequenceProcessor * async_sp = me->async_sp;

    unsigned int n = async_sp->n_parsed();

    return PyLong_FromUnsignedLongLong(n);
}
static PyObject * asyncseqproc_n_processed(PyObject * self, PyObject * args)
{
    khmer_AsyncSequenceProcessorObject * me = (khmer_AsyncSequenceProcessorObject *) self;
    AsyncSequenceProcessor * async_sp = me->async_sp;

    unsigned int n = async_sp->n_processed();

    return PyLong_FromUnsignedLongLong(n);
}

static void khmer_asyncseqproc_dealloc(PyObject* obj)
{
    khmer_AsyncSequenceProcessorObject * self = (khmer_AsyncSequenceProcessorObject *) obj;

    self->async_sp->stop();
    //delete self->async_diginorm;
    self->async_sp = NULL;

    obj->ob_type->tp_free((PyObject*)self);
}

static PyObject * khmer_asyncseqproc_new(PyTypeObject *type, PyObject *args,
                                      PyObject *kwds)
{
    khmer_AsyncSequenceProcessorObject *self;
    self = (khmer_AsyncSequenceProcessorObject*)type->tp_alloc(type, 0);

    if (self != NULL) {
        khmer_KCountingHashObject * counting_o;

        if (!PyArg_ParseTuple(args, "O!", &khmer_KCountingHashType, &counting_o)) {
            return NULL;
        }
        
        //self->async_sp = new AsyncSequenceProcessor(counting_o->counting);
    }

    return (PyObject *) self;
}

static int khmer_asyncseqproc_init(khmer_AsyncSequenceProcessorObject * self, PyObject *args,
                                PyObject *kwds)
{
    return 0;
}

////////////////////
// AsyncDiginorm
////////////////////

static void khmer_asyncdiginorm_dealloc(PyObject *);
static int khmer_asyncdiginorm_init(khmer_AsyncDiginormObject * self, PyObject * args,
                                PyObject * kwds);
static PyObject * khmer_asyncdiginorm_new(PyTypeObject * type, PyObject * args,
                                PyObject * kwds);

static PyObject * asyncdiginorm_start(PyObject * self, PyObject * args)
{
    khmer_AsyncDiginormObject * me = (khmer_AsyncDiginormObject *) self;
    AsyncDiginorm * async_diginorm = me->async_diginorm;

    const char * filename;
    unsigned int cutoff = 20;
    unsigned int n_threads = 1;

    if (!PyArg_ParseTuple(args, "s|II", &filename, &cutoff, &n_threads)) {
        return NULL;
    }
    async_diginorm->start(filename, cutoff, n_threads);

    Py_RETURN_NONE;
}

static PyObject * asyncdiginorm_n_kept(PyObject * self, PyObject * args)
{
    khmer_AsyncDiginormObject * me = (khmer_AsyncDiginormObject *) self;
    AsyncDiginorm * async_diginorm = me->async_diginorm;

    unsigned int n = async_diginorm->n_kept();

    return PyLong_FromUnsignedLongLong(n);
}

static void khmer_asyncdiginorm_dealloc(PyObject* obj)
{
    khmer_AsyncDiginormObject * self = (khmer_AsyncDiginormObject *) obj;

    self->async_diginorm->stop();
    //delete self->async_diginorm;
    self->async_diginorm = NULL;

    obj->ob_type->tp_free((PyObject*)self);
}

static PyObject * khmer_asyncdiginorm_new(PyTypeObject *type, PyObject *args,
                                      PyObject *kwds)
{
    khmer_AsyncDiginormObject *self;
    self = (khmer_AsyncDiginormObject*)type->tp_alloc(type, 0);

    if (self != NULL) {
        khmer_KCountingHashObject * counting_o;

        if (!PyArg_ParseTuple(args, "O!", &khmer_KCountingHashType, &counting_o)) {
            return NULL;
        }
        
        self->async_diginorm = new AsyncDiginorm(counting_o->counting);
        self->async_sp.async_sp = (AsyncSequenceProcessor *)self->async_diginorm;
    }

    return (PyObject *) self;
}

static int khmer_asyncdiginorm_init(khmer_AsyncDiginormObject * self, PyObject *args,
                                PyObject *kwds)
{
    if (khmer_AsyncSequenceProcessorType.tp_init((PyObject *)self, args, kwds) < 0) {
        return -1;
    }
    return 0;
}


PyMODINIT_FUNC
init_khmer_async(void)
{
    using namespace python;

    khmer_AsyncSequenceProcessorType.tp_new = khmer_asyncseqproc_new;
    if (PyType_Ready(&khmer_AsyncSequenceProcessorType) < 0) {
        return;
    }

    khmer_AsyncDiginormType.tp_new = khmer_asyncdiginorm_new;
    khmer_AsyncDiginormType.tp_base = &khmer_AsyncSequenceProcessorType;
    if (PyType_Ready(&khmer_AsyncDiginormType) < 0) {
        return;
    }

    PyObject * m;
    m = Py_InitModule3( "_khmer_async", NULL,
                        "interface for the khmer_async module low-level extensions" );
    if (m == NULL) {
        return;
    }

    Py_INCREF(&khmer_AsyncSequenceProcessorType);
    PyModule_AddObject(m, "AsyncSequenceProcessor", 
        (PyObject *)&khmer_AsyncSequenceProcessorType);

    Py_INCREF(&khmer_AsyncDiginormType);
    PyModule_AddObject(m, "AsyncDiginorm", (PyObject *)&khmer_AsyncDiginormType);
}


