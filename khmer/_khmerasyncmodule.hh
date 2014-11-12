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

#ifndef KHMER_ASYNC_MODULE_HH
#define KHMER_ASYNC_MODULE_HH

#include <Python.h>

#include "khmer_async.hh"

using namespace khmer;


typedef struct {
    PyObject_HEAD
    AsyncSequenceProcessor * async_sp;
} khmer_AsyncSequenceProcessorObject;

PyAPI_DATA(PyTypeObject) khmer_AsyncSequenceProcessorType;

////////////////////
// AsyncDiginorm
////////////////////

typedef struct {
    khmer_AsyncSequenceProcessorObject async_sp;
    AsyncDiginorm * async_diginorm;
} khmer_AsyncDiginormObject;

PyAPI_DATA(PyTypeObject) khmer_AsyncDiginormType;

#endif
