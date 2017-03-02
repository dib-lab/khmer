/*
This file is part of khmer, https://github.com/dib-lab/khmer/, and is
Copyright (C) 2010-2015, Michigan State University.
Copyright (C) 2015-2016, The Regents of the University of California.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are
met:

    * Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.

    * Redistributions in binary form must reproduce the above
      copyright notice, this list of conditions and the following
      disclaimer in the documentation and/or other materials provided
      with the distribution.

    * Neither the name of the Michigan State University nor the names
      of its contributors may be used to endorse or promote products
      derived from this software without specific prior written
      permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
"AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
LICENSE (END)

Contact: khmer-project@idyll.org
*/

//
// A module for Python that exports khmer C++ library functions.
//

// Must be first.
#include <Python.h>

#include <iostream>

#include "khmer.hh"
#include "kmer_hash.hh"
#include "hashtable.hh"
#include "hashgraph.hh"
#include "assembler.hh"
#include "read_aligner.hh"
#include "labelhash.hh"
#include "khmer_exception.hh"
#include "hllcounter.hh"

using namespace khmer;
using namespace read_parsers;

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

using namespace khmer;

//
// Function necessary for Python loading:
//

extern "C" {
    MOD_INIT(_khmer);
}

/***********************************************************************/

// Convert a hash to a python long object.
static bool convert_HashIntoType_to_PyObject(const HashIntoType &hashval,
        PyObject **value)
{
    *value = PyLong_FromUnsignedLongLong(hashval);
    return true;
}


// Convert a python long to a hash
static bool convert_PyLong_to_HashIntoType(PyObject * value,
        HashIntoType &hashval)
{
    if (PyLong_Check(value)) {
        //(PyLongObject *)
        hashval = PyLong_AsUnsignedLongLong(value);
        return true;
    } else if (PyInt_Check(value)) {
        hashval = PyInt_AsLong(value);
        return true;
    } else {
        PyErr_SetString(PyExc_ValueError, "could not convert to hash");
        return false;
    }
}


// Take a Python object and (try to) convert it to a HashIntoType.
// Note: will set error condition and return false if cannot do.

static bool convert_PyObject_to_HashIntoType(PyObject * value,
        HashIntoType& hashval,
        WordLength ksize)
{
    if (PyInt_Check(value) || PyLong_Check(value)) {
        return convert_PyLong_to_HashIntoType(value, hashval);
    } else {
        PyErr_SetString(PyExc_ValueError,
                        "must use a hash");
        return false;
    }
}

// Take a Python object and (try to) convert it to a HashIntoType.
// Note: will set error condition and return false if cannot do.
// Further note: the main difference between this and
// ht_convert_PyObject_to_Kmer is that this will not pass HashIntoType
// numbers through the Kmer class, which means reverse complements
// will not be calculated.  There is a test in test_nodegraph.py
// that checks this.

static bool ht_convert_PyObject_to_HashIntoType(PyObject * value,
        HashIntoType& hashval,
        const Hashtable * ht)
{
    if (PyInt_Check(value) || PyLong_Check(value)) {
        return convert_PyLong_to_HashIntoType(value, hashval);
    } else if (PyUnicode_Check(value))  {
        PyObject* val_as_str = PyUnicode_AsEncodedString(value,
                               "utf-8", "strict");
        std::string s = PyBytes_AsString(val_as_str);
        if (strlen(s.c_str()) != ht->ksize()) {
            Py_DECREF(val_as_str);
            PyErr_SetString(PyExc_ValueError,
                            "k-mer length must equal the k-mer size");
            return false;
        }

        try {
            hashval = ht->hash_dna(s.c_str());
        } catch (khmer_exception &e) {
            PyErr_SetString(PyExc_ValueError, e.what());
            Py_DECREF(val_as_str);
            return false;
        }

        Py_DECREF(val_as_str);
        return true;

    } else if (PyBytes_Check(value)) {
        std::string s = PyBytes_AsString(value);
        if (strlen(s.c_str()) != ht->ksize()) {
            PyErr_SetString(PyExc_ValueError,
                            "k-mer length must equal the k-mer size");
            return false;
        }
        try {
            hashval = ht->hash_dna(s.c_str());
        } catch (khmer_exception &e) {
            PyErr_SetString(PyExc_ValueError, e.what());
            return false;
        }
        return true;
    } else {
        PyErr_SetString(PyExc_ValueError,
                        "k-mers must be either a hash or a string");
        return false;
    }
}

// Take a Python object and (try to) convert it to a khmer::Kmer.
// Note: will set error condition and return false if cannot do.

static bool ht_convert_PyObject_to_Kmer(PyObject * value,
                                        Kmer& kmer, const Hashtable * ht)
{
    if (PyInt_Check(value) || PyLong_Check(value)) {
        HashIntoType h;
        if (!convert_PyLong_to_HashIntoType(value, h)) {
            return false;
        }
        kmer.set_from_unique_hash(h, ht->ksize());
        return true;
    } else if (PyUnicode_Check(value))  {
        std::string s = PyBytes_AsString(PyUnicode_AsEncodedString(
                                             value, "utf-8", "strict"));
        if (strlen(s.c_str()) != ht->ksize()) {
            PyErr_SetString(PyExc_ValueError,
                            "k-mer length must equal the k-mer size");
            return false;
        }
        kmer = Kmer(s, ht->ksize());
        return true;

    } else if (PyBytes_Check(value)) {
        std::string s = PyBytes_AsString(value);
        if (strlen(s.c_str()) != ht->ksize()) {
            PyErr_SetString(PyExc_ValueError,
                            "k-mer length must equal the k-mer size");
            return false;
        }
        kmer = Kmer(s, ht->ksize());
        return true;
    } else {
        PyErr_SetString(PyExc_ValueError,
                        "k-mers must be either a hash or a string");
        return false;
    }
}


static bool convert_Pytablesizes_to_vector(PyListObject * sizes_list_o,
        std::vector<uint64_t>& sizes)
{
    Py_ssize_t sizes_list_o_length = PyList_GET_SIZE(sizes_list_o);
    if (sizes_list_o_length < 1) {
        PyErr_SetString(PyExc_ValueError,
                        "tablesizes needs to be one or more numbers");
        return false;
    }
    for (Py_ssize_t i = 0; i < sizes_list_o_length; i++) {
        PyObject * size_o = PyList_GET_ITEM(sizes_list_o, i);
        if (PyLong_Check(size_o)) {
            sizes.push_back(PyLong_AsUnsignedLongLong(size_o));
        } else if (PyInt_Check(size_o)) {
            sizes.push_back(PyInt_AsLong(size_o));
        } else if (PyFloat_Check(size_o)) {
            sizes.push_back(PyFloat_AS_DOUBLE(size_o));
        } else {
            PyErr_SetString(PyExc_TypeError,
                            "2nd argument must be a list of ints, longs, or floats");
            return false;
        }
    }
    return true;
}


static FastxParserPtr& _PyObject_to_khmer_ReadParser(PyObject * py_object);

/***********************************************************************/

//
// Read object -- name, sequence, and FASTQ stuff
//

namespace khmer
{

namespace python
{

typedef struct {
    PyObject_HEAD
    //! Pointer to the low-level genomic read object.
    read_parsers:: Read *   read;
} khmer_Read_Object;


static
PyObject*
khmer_Read_new(PyTypeObject * type, PyObject * args, PyObject * kwds)
{
    khmer_Read_Object * self;
    self = (khmer_Read_Object *)type->tp_alloc(type, 0);
    if (self != NULL) {
        try {
            self->read = new Read;
        } catch (std::bad_alloc &exc) {
            Py_DECREF(self);
            return PyErr_NoMemory();
        }
    }
    return (PyObject *)self;
}

static
int
khmer_Read_init(khmer_Read_Object *self, PyObject *args, PyObject *kwds)
{
    const char * name{};
    const char * description{};
    const char * sequence{};
    const char * quality{};
    char *kwlist[5] = {
        const_cast<char *>("name"), const_cast<char *>("sequence"),
        const_cast<char *>("quality"), const_cast<char *>("description"), NULL
    };

    if (!PyArg_ParseTupleAndKeywords(args, kwds, "ss|zz", kwlist,
                                     &name, &sequence, &quality, &description)) {
        return -1;
    }

    if (name != NULL) {
        self->read->name = name;
    }
    if (sequence != NULL) {
        self->read->sequence = sequence;
        self->read->set_clean_seq();
    }
    if (quality != NULL) {
        self->read->quality = quality;
    }
    if (description != NULL) {
        self->read->description = description;
    }
    return 0;
}

static
void
khmer_Read_dealloc(khmer_Read_Object * obj)
{
    delete obj->read;
    obj->read = NULL;
    Py_TYPE(obj)->tp_free((PyObject*)obj);
}


static Py_ssize_t
khmer_Read_len(khmer_Read_Object* obj)
{
    return obj->read->sequence.size();
}

static PySequenceMethods khmer_Read_sequence_methods = {
    (lenfunc)khmer_Read_len,                  /* sq_length */
};


static
PyObject *
Read_get_name(khmer_Read_Object * obj, void * closure )
{
    if (obj->read->name.size() > 0) {
        return PyUnicode_FromString(obj->read->name.c_str());
    } else {
        PyErr_SetString(PyExc_AttributeError,
                        "'Read' object has no attribute 'name'.");
        return NULL;
    }
}


static
PyObject *
Read_get_sequence(khmer_Read_Object * obj, void * closure)
{
    if (obj->read->sequence.size() > 0) {
        return PyUnicode_FromString(obj->read->sequence.c_str());
    } else {
        PyErr_SetString(PyExc_AttributeError,
                        "'Read' object has no attribute 'sequence'.");
        return NULL;
    }
}


static
PyObject *
Read_get_quality(khmer_Read_Object * obj, void * closure)
{
    if (obj->read->quality.size() > 0) {
        return PyUnicode_FromString(obj->read->quality.c_str());
    } else {
        PyErr_SetString(PyExc_AttributeError,
                        "'Read' object has no attribute 'quality'.");
        return NULL;
    }
}


static
PyObject *
Read_get_description(khmer_Read_Object * obj, void * closure)
{
    if (obj->read->description.size() > 0) {
        return PyUnicode_FromString(obj->read->description.c_str());
    } else {
        PyErr_SetString(PyExc_AttributeError,
                        "'Read' object has no attribute 'description'.");
        return NULL;
    }
}


static
PyObject *
Read_get_cleaned_seq(khmer_Read_Object * obj, void * closure)
{
    if (obj->read->cleaned_seq.size() > 0) {
        return PyUnicode_FromString(obj->read->cleaned_seq.c_str());
    } else if (obj->read->sequence.size() > 0) {
        obj->read->set_clean_seq();
        return PyUnicode_FromString(obj->read->cleaned_seq.c_str());
    } else {
        PyErr_SetString(PyExc_AttributeError,
                        "'Read' object has no attribute 'cleaned_seq'.");
        return NULL;
    }
}


// TODO? Implement setters.


static PyGetSetDef khmer_Read_accessors [ ] = {
    {
        (char *)"name",
        (getter)Read_get_name, (setter)NULL,
        (char *)"Read identifier.", NULL
    },
    {
        (char *)"sequence",
        (getter)Read_get_sequence, (setter)NULL,
        (char *)"Genomic sequence.", NULL
    },
    {
        (char *)"quality",
        (getter)Read_get_quality, (setter)NULL,
        (char *)"Quality scores.", NULL
    },
    {
        (char *)"description",
        (getter)Read_get_description, (setter)NULL,
        (char *)"Description.", NULL
    },
    {
        (char *)"cleaned_seq",
        (getter)Read_get_cleaned_seq, (setter)NULL,
        (char *)"Cleaned sequence.", NULL
    },

    { NULL, NULL, NULL, NULL, NULL } // sentinel
};


static PyTypeObject khmer_Read_Type = {
    PyVarObject_HEAD_INIT(NULL, 0)        /* init & ob_size */
    "khmer.Read",                         /* tp_name */
    sizeof(khmer_Read_Object),            /* tp_basicsize */
    0,                                    /* tp_itemsize */
    (destructor)khmer_Read_dealloc,       /* tp_dealloc */
    0,                                    /* tp_print */
    0,                                    /* tp_getattr */
    0,                                    /* tp_setattr */
    0,                                    /* tp_compare */
    0,                                    /* tp_repr */
    0,                                    /* tp_as_number */
    &khmer_Read_sequence_methods,         /* tp_as_sequence */
    0,                                    /* tp_as_mapping */
    0,                                    /* tp_hash */
    0,                                    /* tp_call */
    0,                                    /* tp_str */
    0,                                    /* tp_getattro */
    0,                                    /* tp_setattro */
    0,                                    /* tp_as_buffer */
    Py_TPFLAGS_DEFAULT,                   /* tp_flags */
    "A FASTQ record plus some metadata.", /* tp_doc */
    0,                                    /* tp_traverse */
    0,                                    /* tp_clear */
    0,                                    /* tp_richcompare */
    0,                                    /* tp_weaklistoffset */
    0,                                    /* tp_iter */
    0,                                    /* tp_iternext */
    0,                                    /* tp_methods */
    0,                                    /* tp_members */
    (PyGetSetDef *)khmer_Read_accessors,  /* tp_getset */
    0,                                    /* tp_base */
    0,                                    /* tp_dict */
    0,                                    /* tp_descr_get */
    0,                                    /* tp_descr_set */
    0,                                    /* tp_dictoffset */
    (initproc)khmer_Read_init,            /* tp_init */
    0,                                    /* tp_alloc */
    khmer_Read_new,                       /* tp_new */
};


/***********************************************************************/

//
// ReadParser object -- parse reads directly from streams
// ReadPairIterator -- return pairs of Read objects
//


typedef struct {
    PyObject_HEAD
    FastxParserPtr parser;
} khmer_ReadParser_Object;


typedef struct {
    PyObject_HEAD
    //! Pointer to Python parser object for reference counting purposes.
    PyObject * parent;
    //! Persistent value of pair mode across invocations.
    int pair_mode;
} khmer_ReadPairIterator_Object;


static
void
_ReadParser_dealloc(khmer_ReadParser_Object * obj)
{
    Py_TYPE(obj)->tp_free((PyObject*)obj);
}


static
void
khmer_ReadPairIterator_dealloc(khmer_ReadPairIterator_Object * obj)
{
    Py_DECREF(obj->parent);
    obj->parent = NULL;
    Py_TYPE(obj)->tp_free((PyObject*)obj);
}


static PyObject *
_ReadParser_new( PyTypeObject * subtype, PyObject * args, PyObject * kwds )
{
    const char *      ifile_name_CSTR;

    if (!PyArg_ParseTuple(args, "s", &ifile_name_CSTR )) {
        return NULL;
    }
    std:: string    ifile_name( ifile_name_CSTR );

    PyObject * self     = subtype->tp_alloc( subtype, 1 );
    if (self == NULL) {
        return NULL;
    }
    khmer_ReadParser_Object * myself  = (khmer_ReadParser_Object *)self;

    // Wrap the low-level parser object.
    try {
        myself->parser = get_parser<FastxReader>(ifile_name);
    } catch (khmer_file_exception &exc) {
        PyErr_SetString( PyExc_OSError, exc.what() );
        return NULL;
    }
    return self;
}


static
PyObject *
_ReadParser_iternext( PyObject * self )
{
    khmer_ReadParser_Object * myself  = (khmer_ReadParser_Object *)self;
    FastxParserPtr& parser = myself->parser;
    std::string exc_string;

    bool        stop_iteration  = false;
    const char *value_exception = NULL;
    const char *file_exception  = NULL;
    Read       *the_read_PTR    = NULL;
    try {
        the_read_PTR = new Read( );
    } catch (std::bad_alloc &exc) {
        return PyErr_NoMemory();
    }

    Py_BEGIN_ALLOW_THREADS
    stop_iteration = parser->is_complete( );
    if (!stop_iteration) {
        try {
            *the_read_PTR = parser->get_next_read();
        } catch (NoMoreReadsAvailable &exc) {
            stop_iteration = true;
        } catch (khmer_file_exception &exc) {
            exc_string = exc.what();
            file_exception = exc_string.c_str();
        } catch (khmer_value_exception &exc) {
            exc_string = exc.what();
            value_exception = exc_string.c_str();
        }
    }
    Py_END_ALLOW_THREADS

    // Note: Can simply return NULL instead of setting the StopIteration
    //       exception.
    if (stop_iteration) {
        delete the_read_PTR;
        return NULL;
    }

    if (file_exception != NULL) {
        delete the_read_PTR;
        PyErr_SetString(PyExc_OSError, file_exception);
        return NULL;
    }
    if (value_exception != NULL) {
        delete the_read_PTR;
        PyErr_SetString(PyExc_ValueError, value_exception);
        return NULL;
    }

    PyObject * the_read_OBJECT = khmer_Read_Type.tp_alloc( &khmer_Read_Type, 1 );
    ((khmer_Read_Object *)the_read_OBJECT)->read = the_read_PTR;
    return the_read_OBJECT;
}


static
PyObject *
_ReadPairIterator_iternext(khmer_ReadPairIterator_Object * myself)
{
    khmer_ReadParser_Object * parent = (khmer_ReadParser_Object*)myself->parent;
    FastxParserPtr& parser = parent->parser;
    uint8_t     pair_mode = myself->pair_mode;

    ReadPair    the_read_pair;
    bool        stop_iteration  = false;
    const char *value_exception = NULL;
    const char *file_exception  = NULL;
    std::string exc_string;

    Py_BEGIN_ALLOW_THREADS
    stop_iteration = parser->is_complete( );
    if (!stop_iteration) {
        try {
            the_read_pair = parser->get_next_read_pair(pair_mode);
        } catch (NoMoreReadsAvailable &exc) {
            stop_iteration = true;
        } catch (khmer_file_exception &exc) {
            exc_string = exc.what();
            file_exception = exc_string.c_str();
        } catch (khmer_value_exception &exc) {
            exc_string = exc.what();
            value_exception = exc_string.c_str();
        }
    }
    Py_END_ALLOW_THREADS

    // Note: Can return NULL instead of setting the StopIteration exception.
    if (stop_iteration) {
        return NULL;
    }
    if (file_exception != NULL) {
        PyErr_SetString(PyExc_OSError, file_exception);
        return NULL;
    }
    if (value_exception != NULL) {
        PyErr_SetString(PyExc_ValueError, value_exception);
        return NULL;
    }

    // Copy elements of 'ReadPair' object into Python tuple.
    // TODO? Replace dummy reads with 'None' object.
    PyObject * read_1_OBJECT = khmer_Read_Type.tp_alloc( &khmer_Read_Type, 1 );
    try {
        ((khmer_Read_Object *)read_1_OBJECT)->read = new Read( the_read_pair.first );
    } catch (std::bad_alloc &e) {
        return PyErr_NoMemory();
    }
    PyObject * read_2_OBJECT = khmer_Read_Type.tp_alloc( &khmer_Read_Type, 1 );
    try {
        ((khmer_Read_Object *)read_2_OBJECT)->read = new Read( the_read_pair.second );
    } catch (std::bad_alloc &e) {
        delete ((khmer_Read_Object *)read_1_OBJECT)->read;
        return PyErr_NoMemory();
    }
    PyObject * tup = PyTuple_Pack( 2, read_1_OBJECT, read_2_OBJECT );
    Py_XDECREF(read_1_OBJECT);
    Py_XDECREF(read_2_OBJECT);
    return tup;
}

static PyTypeObject khmer_ReadPairIterator_Type = {
    PyVarObject_HEAD_INIT(NULL, 0)              /* init & ob_size */
    "_khmer.ReadPairIterator",                   /* tp_name */
    sizeof(khmer_ReadPairIterator_Object),      /* tp_basicsize */
    0,                                          /* tp_itemsize */
    (destructor)khmer_ReadPairIterator_dealloc, /* tp_dealloc */
    0,                                          /* tp_print */
    0,                                          /* tp_getattr */
    0,                                          /* tp_setattr */
    0,                                          /* tp_compare */
    0,                                          /* tp_repr */
    0,                                          /* tp_as_number */
    0,                                          /* tp_as_sequence */
    0,                                          /* tp_as_mapping */
    0,                                          /* tp_hash */
    0,                                          /* tp_call */
    0,                                          /* tp_str */
    0,                                          /* tp_getattro */
    0,                                          /* tp_setattro */
    0,                                          /* tp_as_buffer */
    Py_TPFLAGS_DEFAULT,                         /* tp_flags */
    "Iterates over 'ReadParser' objects and returns read pairs.", /* tp_doc */
    0,                                          /* tp_traverse */
    0,                                          /* tp_clear */
    0,                                          /* tp_richcompare */
    0,                                          /* tp_weaklistoffset */
    PyObject_SelfIter,                          /* tp_iter */
    (iternextfunc)_ReadPairIterator_iternext,   /* tp_iternext */
};



static
PyObject *
ReadParser_iter_reads(PyObject * self, PyObject * args )
{
    return PyObject_SelfIter( self );
}

static
PyObject *
ReadParser_get_num_reads(khmer_ReadParser_Object * me)
{
    return PyLong_FromLong(me->parser->get_num_reads());
}

static
PyObject *
ReadParser_iter_read_pairs(PyObject * self, PyObject * args )
{
    int pair_mode = ReadParser<FastxReader>::PAIR_MODE_ERROR_ON_UNPAIRED;

    if (!PyArg_ParseTuple( args, "|i", &pair_mode )) {
        return NULL;
    }

    // Capture existing read parser.
    PyObject * obj = khmer_ReadPairIterator_Type.tp_alloc(
                         &khmer_ReadPairIterator_Type, 1
                     );
    if (obj == NULL) {
        return NULL;
    }
    khmer_ReadPairIterator_Object * rpi   = (khmer_ReadPairIterator_Object *)obj;
    rpi->parent             = self;
    rpi->pair_mode          = pair_mode;

    // Increment reference count on existing ReadParser object so that it
    // will not go away until all ReadPairIterator instances have gone away.
    Py_INCREF( self );

    return obj;
}


PyObject *
ReadParser_close(PyObject * self, PyObject * args)
{
    FastxParserPtr& rparser = _PyObject_to_khmer_ReadParser(self);
    rparser->close();

    Py_INCREF(Py_None);
    return Py_None;
}

static PyMethodDef _ReadParser_methods [ ] = {
    {
        "iter_reads",       (PyCFunction)ReadParser_iter_reads,
        METH_NOARGS,        "Iterates over reads."
    },
    {
        "iter_read_pairs",  (PyCFunction)ReadParser_iter_read_pairs,
        METH_VARARGS,       "Iterates over paired reads as pairs."
    },
    {
        "close",  (PyCFunction)ReadParser_close,
        METH_NOARGS,       "Close associated files."
    },
    { NULL, NULL, 0, NULL } // sentinel
};

static PyGetSetDef khmer_ReadParser_accessors[] = {
    {
        (char *)"num_reads",
        (getter)ReadParser_get_num_reads, NULL,
        (char *)"count of reads processed thus far.",
        NULL
    },
    {NULL, NULL, NULL, NULL, NULL} /* Sentinel */
};

static PyTypeObject khmer_ReadParser_Type
CPYCHECKER_TYPE_OBJECT_FOR_TYPEDEF("khmer_ReadParser_Object")
= {
    PyVarObject_HEAD_INIT(NULL, 0)             /* init & ob_size */
    "_khmer.ReadParser",                        /* tp_name */
    sizeof(khmer_ReadParser_Object),           /* tp_basicsize */
    0,                                         /* tp_itemsize */
    (destructor)_ReadParser_dealloc,           /* tp_dealloc */
    0,                                         /* tp_print */
    0,                                         /* tp_getattr */
    0,                                         /* tp_setattr */
    0,                                         /* tp_compare */
    0,                                         /* tp_repr */
    0,                                         /* tp_as_number */
    0,                                         /* tp_as_sequence */
    0,                                         /* tp_as_mapping */
    0,                                         /* tp_hash */
    0,                                         /* tp_call */
    0,                                         /* tp_str */
    0,                                         /* tp_getattro */
    0,                                         /* tp_setattro */
    0,                                         /* tp_as_buffer */
    Py_TPFLAGS_DEFAULT,                        /* tp_flags */
    "Parses streams from various file formats, " \
    "such as FASTA and FASTQ.",                /* tp_doc */
    0,                                         /* tp_traverse */
    0,                                         /* tp_clear */
    0,                                         /* tp_richcompare */
    0,                                         /* tp_weaklistoffset */
    PyObject_SelfIter,                         /* tp_iter */
    (iternextfunc)_ReadParser_iternext,        /* tp_iternext */
    _ReadParser_methods,                       /* tp_methods */
    0,                                         /* tp_members */
    khmer_ReadParser_accessors,                /* tp_getset */
    0,                                         /* tp_base */
    0,                                         /* tp_dict */
    0,                                         /* tp_descr_get */
    0,                                         /* tp_descr_set */
    0,                                         /* tp_dictoffset */
    0,                                         /* tp_init */
    0,                                         /* tp_alloc */
    _ReadParser_new,                           /* tp_new */
};

void _init_ReadParser_Type_constants()
{
    PyObject * cls_attrs_DICT = PyDict_New( );
    if (cls_attrs_DICT == NULL) {
        return;
    }

    // Place pair mode constants into class dictionary.
    int result;
    PyObject *value = PyLong_FromLong(
                          ReadParser<FastxReader>::PAIR_MODE_IGNORE_UNPAIRED);
    if (value == NULL) {
        Py_DECREF(cls_attrs_DICT);
        return;
    }
    result = PyDict_SetItemString(cls_attrs_DICT,
                                  "PAIR_MODE_IGNORE_UNPAIRED", value );
    Py_XDECREF(value);
    if (!result) {
        Py_DECREF(cls_attrs_DICT);
        return;
    }

    value = PyLong_FromLong(ReadParser<FastxReader>::PAIR_MODE_ERROR_ON_UNPAIRED);
    if (value == NULL) {
        Py_DECREF(cls_attrs_DICT);
        return;
    }
    result = PyDict_SetItemString(cls_attrs_DICT,
                                  "PAIR_MODE_ERROR_ON_UNPAIRED", value);
    Py_XDECREF(value);
    if (!result) {
        Py_DECREF(cls_attrs_DICT);
        return;
    }

    khmer_ReadParser_Type.tp_dict     = cls_attrs_DICT;
}

} // namespace python

} // namespace khmer


static FastxParserPtr& _PyObject_to_khmer_ReadParser(PyObject * py_object)
{
    // TODO: Add type-checking.

    return ((python:: khmer_ReadParser_Object *)py_object)->parser;
}

typedef struct {
    PyObject_HEAD
    pre_partition_info *   PrePartitionInfo;
} khmer_PrePartitionInfo_Object;

static
void
khmer_PrePartitionInfo_dealloc(khmer_PrePartitionInfo_Object * obj)
{
    delete obj->PrePartitionInfo;
    obj->PrePartitionInfo = NULL;
    Py_TYPE(obj)->tp_free((PyObject*)obj);
}

static PyTypeObject khmer_PrePartitionInfo_Type = {
    PyVarObject_HEAD_INIT(NULL, 0)        /* init & ob_size */
    "_khmer.PrePartitionInfo",            /* tp_name */
    sizeof(khmer_PrePartitionInfo_Object),/* tp_basicsize */
    0,                                    /* tp_itemsize */
    (destructor)khmer_PrePartitionInfo_dealloc,       /* tp_dealloc */
    0,                                    /* tp_print */
    0,                                    /* tp_getattr */
    0,                                    /* tp_setattr */
    0,                                    /* tp_compare */
    0,                                    /* tp_repr */
    0,                                    /* tp_as_number */
    0,                                    /* tp_as_sequence */
    0,                                    /* tp_as_mapping */
    0,                                    /* tp_hash */
    0,                                    /* tp_call */
    0,                                    /* tp_str */
    0,                                    /* tp_getattro */
    0,                                    /* tp_setattro */
    0,                                    /* tp_as_buffer */
    Py_TPFLAGS_DEFAULT,                   /* tp_flags */
    "Stores a k-kmer and a set of tagged seen k-mers.", /* tp_doc */
};


/***********************************************************************/
/***********************************************************************/

typedef struct {
    PyObject_HEAD
    SeenSet * hashes;
    WordLength ksize;
} khmer_HashSet_Object;

static khmer_HashSet_Object * create_HashSet_Object(SeenSet * h, WordLength k);

static
void
khmer_HashSet_dealloc(khmer_HashSet_Object * obj)
{
    delete obj->hashes;
    obj->hashes = NULL;
    obj->ksize = 0;
    Py_TYPE(obj)->tp_free((PyObject*)obj);
}

static PyObject* khmer_HashSet_new(PyTypeObject * type, PyObject * args,
                                   PyObject * kwds)
{
    khmer_HashSet_Object * self;

    self = (khmer_HashSet_Object *)type->tp_alloc(type, 0);

    if (self != NULL) {
        PyObject * list_o = NULL;
        WordLength k;
        if (!PyArg_ParseTuple(args, "b|O!", &k, &PyList_Type, &list_o)) {
            Py_DECREF(self);
            return NULL;
        }

        try {
            self->hashes = new SeenSet;
            self->ksize = k;
        } catch (std::bad_alloc &e) {
            Py_DECREF(self);
            return PyErr_NoMemory();
        }

        if (list_o) {
            Py_ssize_t size = PyList_Size(list_o);
            for (Py_ssize_t i = 0; i < size; i++) {
                PyObject * item = PyList_GET_ITEM(list_o, i);
                HashIntoType h;

                if (!convert_PyObject_to_HashIntoType(item, h, self->ksize)) {
                    return NULL;
                }
                self->hashes->insert(h);
            }
        }
    }
    return (PyObject *) self;
}

/***********************************************************************/

typedef struct {
    PyObject_HEAD
    khmer_HashSet_Object * parent;
    SeenSet::iterator * it;
} _HashSet_iterobj;

static
void
_HashSet_iter_dealloc(_HashSet_iterobj * obj)
{
    delete obj->it;
    obj->it = NULL;
    Py_DECREF(obj->parent);
    Py_TYPE(obj)->tp_free((PyObject*)obj);
}

static PyObject * _HashSet_iter(PyObject * self)
{
    Py_INCREF(self);
    return self;
}

static PyObject * _HashSet_iternext(PyObject * self)
{
    _HashSet_iterobj * iter_obj = (_HashSet_iterobj *) self;
    SeenSet * hashes = iter_obj->parent->hashes;
    if (*iter_obj->it != hashes->end()) {
        PyObject * ret = nullptr;
        convert_HashIntoType_to_PyObject((**iter_obj->it), &ret);
        (*(iter_obj->it))++;
        return ret;
    }

    PyErr_SetString(PyExc_StopIteration, "end of HashSet");
    return NULL;
}

static PyTypeObject _HashSet_iter_Type = {
    PyVarObject_HEAD_INIT(NULL, 0)        /* init & ob_size */
    "_khmer.HashSet_iter",                /* tp_name */
    sizeof(_HashSet_iterobj),             /* tp_basicsize */
    0,                                    /* tp_itemsize */
    (destructor)_HashSet_iter_dealloc,    /* tp_dealloc */
    0,                                    /* tp_print */
    0,                                    /* tp_getattr */
    0,                                    /* tp_setattr */
    0,                                    /* tp_compare */
    0,                                    /* tp_repr */
    0,                                    /* tp_as_number */
    0,                                    /* tp_as_sequence */
    0,                                    /* tp_as_mapping */
    0,                                    /* tp_hash */
    0,                                    /* tp_call */
    0,                                    /* tp_str */
    0,                                    /* tp_getattro */
    0,                                    /* tp_setattro */
    0,                                    /* tp_as_buffer */
    Py_TPFLAGS_DEFAULT | Py_TPFLAGS_HAVE_ITER, /* tp_flags */
    "iterator object for HashSet objects.", /* tp_doc */
    0,                                    /* tp_traverse */
    0,                                    /* tp_clear */
    0,                                    /* tp_richcompare */
    0,                                    /* tp_weaklistoffset */
    _HashSet_iter,                        /* tp_iter */
    _HashSet_iternext,                    /* tp_iternext */
    0,                                    /* tp_methods */
    0,                                    /* tp_members */
    0,                                    /* tp_getset */
    0,                                    /* tp_base */
    0,                                    /* tp_dict */
    0,                                    /* tp_descr_get */
    0,                                    /* tp_descr_set */
    0,                                    /* tp_dictoffset */
    0,                                    /* tp_init */
    0,                                    /* tp_alloc */
    0,                                    /* tp_new */
};

static PyObject * khmer_HashSet_iter(PyObject * self)
{
    khmer_HashSet_Object * me = (khmer_HashSet_Object *) self;
    _HashSet_iterobj * iter_obj = (_HashSet_iterobj *)
                                  _HashSet_iter_Type.tp_alloc(&_HashSet_iter_Type, 0);
    if (iter_obj != NULL) {
        Py_INCREF(me);
        iter_obj->parent = me;

        iter_obj->it = new SeenSet::iterator;
        *iter_obj->it = me->hashes->begin();
    }
    return (PyObject *) iter_obj;
}

static int khmer_HashSet_len(khmer_HashSet_Object * o)
{
    return (Py_ssize_t) o->hashes->size();
}

static PyObject * khmer_HashSet_concat(khmer_HashSet_Object * o,
                                       khmer_HashSet_Object * o2)
{
    if (o->ksize != o2->ksize) {
        PyErr_SetString(PyExc_ValueError,
                        "cannot add HashSets with different ksizes");
        return NULL;
    }
    khmer_HashSet_Object * no = create_HashSet_Object(new SeenSet,
                                o->ksize);
    no->hashes->insert(o->hashes->begin(), o->hashes->end());
    no->hashes->insert(o2->hashes->begin(), o2->hashes->end());

    return (PyObject *) no;
}

static PyObject * khmer_HashSet_concat_inplace(khmer_HashSet_Object * o,
        khmer_HashSet_Object * o2)
{
    if (o->ksize != o2->ksize) {
        PyErr_SetString(PyExc_ValueError,
                        "cannot add HashSets with different ksizes");
        return NULL;
    }
    o->hashes->insert(o2->hashes->begin(), o2->hashes->end());

    Py_INCREF(o);
    return (PyObject *) o;
}

static int khmer_HashSet_contains(khmer_HashSet_Object * o, PyObject * val)
{
    HashIntoType v;

    if (convert_PyObject_to_HashIntoType(val, v, 0)) {
        if (set_contains(*o->hashes, v)) {
            return 1;
        }
    }
    return 0;
}

static PyObject *
hashset_add(khmer_HashSet_Object * me, PyObject * args)
{
    PyObject * hash_obj;
    HashIntoType h;
    if (!PyArg_ParseTuple(args, "O", &hash_obj)) {
        return NULL;
    }

    if (!convert_PyObject_to_HashIntoType(hash_obj, h, 0)) {
        return NULL;
    }
    me->hashes->insert(h);

    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *
hashset_remove(khmer_HashSet_Object * me, PyObject * args)
{
    PyObject * hash_obj;
    HashIntoType h;
    if (!PyArg_ParseTuple(args, "O", &hash_obj)) {
        return NULL;
    }

    if (!convert_PyObject_to_HashIntoType(hash_obj, h, 0)) {
        return NULL;
    }
    SeenSet::iterator it = me->hashes->find(h);
    if (it == me->hashes->end()) {
        PyErr_SetString(PyExc_ValueError, "h not in list");
        return NULL;
    }
    me->hashes->erase(it);

    Py_INCREF(Py_None);
    return Py_None;
}

static PyObject *
hashset_update(khmer_HashSet_Object * me, PyObject * args)
{
    PyObject * obj;
    if (!PyArg_ParseTuple(args, "O", &obj)) {
        return NULL;
    }

    PyObject * iterator = PyObject_GetIter(obj);
    if (iterator == NULL) {
        return NULL;
    }
    PyObject * item = PyIter_Next(iterator);
    while(item) {
        HashIntoType h;

        if (!convert_PyObject_to_HashIntoType(item, h, 0)) {
            PyErr_SetString(PyExc_ValueError, "unknown item type for update");
            Py_DECREF(item);
            return NULL;
        }
        me->hashes->insert(h);

        Py_DECREF(item);
        item = PyIter_Next(iterator);
    }
    Py_DECREF(iterator);
    if (PyErr_Occurred()) {
        return NULL;
    }

    Py_INCREF(Py_None);
    return Py_None;
}

static PyMethodDef khmer_HashSet_methods[] = {
    {
        "add",
        (PyCFunction)hashset_add, METH_VARARGS,
        "Add element to the HashSet."
    },
    {
        "remove",
        (PyCFunction)hashset_remove, METH_VARARGS,
        "Remove an element from the HashSet."
    },
    {
        "update",
        (PyCFunction)hashset_update, METH_VARARGS,
        "Add a list of elements to the HashSet."
    },
    {NULL, NULL, 0, NULL}           /* sentinel */
};

static PySequenceMethods khmer_HashSet_seqmethods[] = {
    (lenfunc)khmer_HashSet_len, /* sq_length */
    (binaryfunc)khmer_HashSet_concat,      /* sq_concat */
    0,                          /* sq_repeat */
    0,                          /* sq_item */
    0,                          /* sq_slice */
    0,                          /* sq_ass_item */
    0,                          /* sq_ass_slice */
    (objobjproc)khmer_HashSet_contains, /* sq_contains */
    (binaryfunc)khmer_HashSet_concat_inplace,      /* sq_inplace_concat */
    0                           /* sq_inplace_repeat */
};

static PyTypeObject khmer_HashSet_Type = {
    PyVarObject_HEAD_INIT(NULL, 0)        /* init & ob_size */
    "_khmer.HashSet",                     /* tp_name */
    sizeof(khmer_HashSet_Object),         /* tp_basicsize */
    0,                                    /* tp_itemsize */
    (destructor)khmer_HashSet_dealloc,    /* tp_dealloc */
    0,                                    /* tp_print */
    0,                                    /* tp_getattr */
    0,                                    /* tp_setattr */
    0,                                    /* tp_compare */
    0,                                    /* tp_repr */
    0,                                    /* tp_as_number */
    khmer_HashSet_seqmethods,             /* tp_as_sequence */
    0,                                    /* tp_as_mapping */
    0,                                    /* tp_hash */
    0,                                    /* tp_call */
    0,                                    /* tp_str */
    0,                                    /* tp_getattro */
    0,                                    /* tp_setattro */
    0,                                    /* tp_as_buffer */
    Py_TPFLAGS_DEFAULT | Py_TPFLAGS_HAVE_ITER, /* tp_flags */
    "Stores a set of hashed k-mers.",     /* tp_doc */
    0,                                    /* tp_traverse */
    0,                                    /* tp_clear */
    0,                                    /* tp_richcompare */
    0,                                    /* tp_weaklistoffset */
    khmer_HashSet_iter,                   /* tp_iter */
    0,                                    /* tp_iternext */
    khmer_HashSet_methods,                /* tp_methods */
    0,                                    /* tp_members */
    0,                                    /* tp_getset */
    0,                                    /* tp_base */
    0,                                    /* tp_dict */
    0,                                    /* tp_descr_get */
    0,                                    /* tp_descr_set */
    0,                                    /* tp_dictoffset */
    0,                                    /* tp_init */
    0,                                    /* tp_alloc */
    khmer_HashSet_new,                    /* tp_new */
};

static khmer_HashSet_Object * create_HashSet_Object(SeenSet * h, WordLength k)
{
    khmer_HashSet_Object * self;

    self = (khmer_HashSet_Object *)
           khmer_HashSet_Type.tp_alloc(&khmer_HashSet_Type, 0);
    if (self != NULL) {
        self->hashes = h;
        self->ksize = k;
    }
    return self;
}

/***********************************************************************/

typedef struct {
    PyObject_HEAD
    Hashtable * hashtable;
} khmer_KHashtable_Object;

typedef struct {
    khmer_KHashtable_Object khashtable;
    Hashgraph * hashgraph;
} khmer_KHashgraph_Object;

typedef struct {
    PyObject_HEAD
    SubsetPartition * subset;
} khmer_KSubsetPartition_Object;

static void khmer_subset_dealloc(khmer_KSubsetPartition_Object * obj);

static PyTypeObject khmer_KSubsetPartition_Type = {
    PyVarObject_HEAD_INIT(NULL, 0)         /* init & ob_size */
    "_khmer.KSubsetPartition",              /* tp_name */
    sizeof(khmer_KSubsetPartition_Object), /* tp_basicsize */
    0,                                     /* tp_itemsize */
    (destructor)khmer_subset_dealloc,      /*tp_dealloc*/
    0,                                     /*tp_print*/
    0,                                     /*tp_getattr*/
    0,                                     /*tp_setattr*/
    0,                                     /*tp_compare*/
    0,                                     /*tp_repr*/
    0,                                     /*tp_as_number*/
    0,                                     /*tp_as_sequence*/
    0,                                     /*tp_as_mapping*/
    0,                                     /*tp_hash */
    0,                                     /*tp_call*/
    0,                                     /*tp_str*/
    0,                                     /*tp_getattro*/
    0,                                     /*tp_setattro*/
    0,                                     /*tp_as_buffer*/
    Py_TPFLAGS_DEFAULT,                    /*tp_flags*/
    "subset object",                       /* tp_doc */
};

typedef struct {
    khmer_KHashgraph_Object khashgraph;
    Nodegraph * nodegraph;
} khmer_KNodegraph_Object;

static void khmer_nodegraph_dealloc(khmer_KNodegraph_Object * obj);
static PyObject* khmer_nodegraph_new(PyTypeObject * type, PyObject * args,
                                     PyObject * kwds);

static PyTypeObject khmer_KNodegraph_Type
CPYCHECKER_TYPE_OBJECT_FOR_TYPEDEF("khmer_KNodegraph_Object")
= {
    PyVarObject_HEAD_INIT(NULL, 0) /* init & ob_size */
    "_khmer.Nodegraph",             /* tp_name */
    sizeof(khmer_KNodegraph_Object), /* tp_basicsize */
    0,                             /* tp_itemsize */
    (destructor)khmer_nodegraph_dealloc, /*tp_dealloc*/
    0,              /*tp_print*/
    0,              /*tp_getattr*/
    0,              /*tp_setattr*/
    0,              /*tp_compare*/
    0,              /*tp_repr*/
    0,              /*tp_as_number*/
    0,              /*tp_as_sequence*/
    0,              /*tp_as_mapping*/
    0,              /*tp_hash */
    0,              /*tp_call*/
    0,              /*tp_str*/
    0,              /*tp_getattro*/
    0,              /*tp_setattro*/
    0,              /*tp_as_buffer*/
    Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE,       /*tp_flags*/
    "nodegraph object",           /* tp_doc */
    0,                       /* tp_traverse */
    0,                       /* tp_clear */
    0,                       /* tp_richcompare */
    0,                       /* tp_weaklistoffset */
    0,                       /* tp_iter */
    0,                       /* tp_iternext */
    0,  /* tp_methods */
    0,                       /* tp_members */
    0,                       /* tp_getset */
    0,                       /* tp_base */
    0,                       /* tp_dict */
    0,                       /* tp_descr_get */
    0,                       /* tp_descr_set */
    0,                       /* tp_dictoffset */
    0,                       /* tp_init */
    0,                       /* tp_alloc */
    khmer_nodegraph_new,                  /* tp_new */
};


static
PyObject *
hashtable_ksize(khmer_KHashtable_Object * me, PyObject * args)
{
    Hashtable * hashtable = me->hashtable;

    if (!PyArg_ParseTuple(args, "")) {
        return NULL;
    }

    unsigned int k = hashtable->ksize();

    return PyLong_FromLong(k);
}

static
PyObject *
hashtable_hash(khmer_KHashtable_Object * me, PyObject * args)
{
    Hashtable * hashtable = me->hashtable;

    char * kmer;
    if (!PyArg_ParseTuple(args, "s", &kmer)) {
        return NULL;
    }

    if (strlen(kmer) != hashtable->ksize()) {
        PyErr_SetString(PyExc_ValueError,
                        "provided k-mer is wrong length");
        return NULL;
    }

    try {
        PyObject * hash = nullptr;
        const HashIntoType h(hashtable->hash_dna(kmer));
        convert_HashIntoType_to_PyObject(h, &hash);
        return hash;
    } catch (khmer_exception &e) {
        PyErr_SetString(PyExc_ValueError, e.what());
        return NULL;
    }
}

static
PyObject *
hashtable_reverse_hash(khmer_KHashtable_Object * me, PyObject * args)
{
    Hashtable * hashtable = me->hashtable;

    PyObject * val_o;
    HashIntoType val;
    if (!PyArg_ParseTuple(args, "O", &val_o)) {
        return NULL;
    }

    if (!ht_convert_PyObject_to_HashIntoType(val_o, val, hashtable)) {
        return NULL;
    }

    try {
        return PyUnicode_FromString(hashtable->unhash_dna(val).c_str());
    } catch (khmer_exception &e) {
        PyErr_SetString(PyExc_ValueError, e.what());
        return NULL;
    }
}

static
PyObject *
hashtable_n_occupied(khmer_KHashtable_Object * me, PyObject * args)
{
    Hashtable * hashtable = me->hashtable;

    if (!PyArg_ParseTuple(args, "")) {
        return NULL;
    }

    uint64_t n = hashtable->n_occupied();

    return PyLong_FromUnsignedLongLong(n);
}

static
PyObject *
hashtable_n_unique_kmers(khmer_KHashtable_Object * me, PyObject * args)
{
    Hashtable * hashtable = me->hashtable;

    uint64_t n = hashtable->n_unique_kmers();

    return PyLong_FromUnsignedLongLong(n);
}

static
PyObject *
hashtable_count(khmer_KHashtable_Object * me, PyObject * args)
{
    Hashtable * hashtable = me->hashtable;

    PyObject * v;
    if (!PyArg_ParseTuple(args, "O", &v)) {
        return NULL;
    }

    HashIntoType hashval;

    if (!ht_convert_PyObject_to_HashIntoType(v, hashval, hashtable)) {
        return NULL;
    }

    hashtable->count(hashval);

    return PyLong_FromLong(1);
}

static
PyObject *
hashtable_consume_seqfile(khmer_KHashtable_Object * me, PyObject * args)
{
    Hashtable * hashtable  = me->hashtable;

    const char * filename;

    if (!PyArg_ParseTuple(args, "s", &filename)) {
        return NULL;
    }

    // call the C++ function, and trap signals => Python
    unsigned long long  n_consumed    = 0;
    unsigned int          total_reads   = 0;
    try {
        hashtable->consume_seqfile<FastxReader>(filename, total_reads, n_consumed);
    } catch (khmer_file_exception &exc) {
        PyErr_SetString(PyExc_OSError, exc.what());
        return NULL;
    } catch (khmer_value_exception &exc) {
        PyErr_SetString(PyExc_ValueError, exc.what());
        return NULL;
    }

    return Py_BuildValue("IK", total_reads, n_consumed);
}

static
PyObject *
hashtable_consume_seqfile_with_reads_parser(khmer_KHashtable_Object * me,
        PyObject * args)
{
    Hashtable * hashtable = me->hashtable;

    PyObject * rparser_obj = NULL;

    if (!PyArg_ParseTuple(args, "O", &rparser_obj)) {
        return NULL;
    }

    FastxParserPtr& rparser = _PyObject_to_khmer_ReadParser( rparser_obj );

    // call the C++ function, and trap signals => Python
    unsigned long long  n_consumed      = 0;
    unsigned int        total_reads     = 0;
    const char         *value_exception = NULL;
    const char         *file_exception  = NULL;
    std::string exc_string;

    Py_BEGIN_ALLOW_THREADS
    try {
        hashtable->consume_seqfile<FastxReader>(rparser, total_reads, n_consumed);
    } catch (khmer_file_exception &exc) {
        exc_string = exc.what();
        file_exception = exc_string.c_str();
    } catch (khmer_value_exception &exc) {
        exc_string = exc.what();
        value_exception = exc_string.c_str();
    }
    Py_END_ALLOW_THREADS

    if (file_exception != NULL) {
        PyErr_SetString(PyExc_OSError, file_exception);
        return NULL;
    }
    if (value_exception != NULL) {
        PyErr_SetString(PyExc_ValueError, value_exception);
        return NULL;
    }

    return Py_BuildValue("IK", total_reads, n_consumed);
}

static
PyObject *
hashtable_consume(khmer_KHashtable_Object * me, PyObject * args)
{
    Hashtable * hashtable = me->hashtable;

    const char * long_str;

    if (!PyArg_ParseTuple(args, "s", &long_str)) {
        return NULL;
    }

    if (strlen(long_str) < hashtable->ksize()) {
        PyErr_SetString(PyExc_ValueError,
                        "string length must >= the hashtable k-mer size");
        return NULL;
    }

    unsigned int n_consumed;
    n_consumed = hashtable->consume_string(long_str);

    return PyLong_FromLong(n_consumed);
}

static
PyObject *
hashtable_get(khmer_KHashtable_Object * me, PyObject * args)
{
    Hashtable * hashtable = me->hashtable;

    PyObject * arg;

    if (!PyArg_ParseTuple(args, "O", &arg)) {
        return NULL;
    }

    HashIntoType hashval;

    if (!ht_convert_PyObject_to_HashIntoType(arg, hashval, hashtable)) {
        return NULL;
    }

    unsigned int count = hashtable->get_count(hashval);
    return PyLong_FromLong(count);
}

static
PyObject *
hashtable_set_use_bigcount(khmer_KHashtable_Object * me, PyObject * args)
{
    Hashtable * hashtable = me->hashtable;

    PyObject * x;
    if (!PyArg_ParseTuple(args, "O", &x)) {
        return NULL;
    }
    int setme = PyObject_IsTrue(x);
    if (setme < 0) {
        return NULL;
    }
    try {
        hashtable->set_use_bigcount((bool)setme);
    } catch (khmer_exception &e) {
        PyErr_SetString(PyExc_ValueError, e.what());
        return NULL;
    }

    Py_RETURN_NONE;
}

static
PyObject *
hashtable_get_use_bigcount(khmer_KHashtable_Object * me, PyObject * args)
{
    Hashtable * hashtable = me->hashtable;

    if (!PyArg_ParseTuple(args, "")) {
        return NULL;
    }

    bool val = hashtable->get_use_bigcount();

    return PyBool_FromLong((int)val);
}

static
PyObject *
hashtable_get_min_count(khmer_KHashtable_Object * me, PyObject * args)
{
    Hashtable * hashtable = me->hashtable;

    const char * long_str;

    if (!PyArg_ParseTuple(args, "s", &long_str)) {
        return NULL;
    }

    if (strlen(long_str) < hashtable->ksize()) {
        PyErr_SetString(PyExc_ValueError,
                        "string length must >= the hashtable k-mer size");
        return NULL;
    }

    BoundedCounterType c = hashtable->get_min_count(long_str);
    unsigned int N = c;

    return PyLong_FromLong(N);
}

static
PyObject *
hashtable_get_max_count(khmer_KHashtable_Object * me, PyObject * args)
{
    Hashtable * hashtable = me->hashtable;

    const char * long_str;

    if (!PyArg_ParseTuple(args, "s", &long_str)) {
        return NULL;
    }

    if (strlen(long_str) < hashtable->ksize()) {
        PyErr_SetString(PyExc_ValueError,
                        "string length must >= the hashtable k-mer size");
        return NULL;
    }

    BoundedCounterType c = hashtable->get_max_count(long_str);
    unsigned int N = c;

    return PyLong_FromLong(N);
}

// fwd decls to replace function definitions
static PyObject * hashtable_abundance_distribution_with_reads_parser(khmer_KHashtable_Object * me, PyObject * args);
static PyObject * hashtable_abundance_distribution(khmer_KHashtable_Object * me, PyObject * args);

static
PyObject *
hashtable_trim_on_abundance(khmer_KHashtable_Object * me, PyObject * args)
{
    Hashtable * hashtable = me->hashtable;

    const char * seq = NULL;
    unsigned int min_count_i = 0;

    if (!PyArg_ParseTuple(args, "sI", &seq, &min_count_i)) {
        return NULL;
    }

    unsigned long trim_at;
    Py_BEGIN_ALLOW_THREADS

    BoundedCounterType min_count = min_count_i;

    trim_at = hashtable->trim_on_abundance(seq, min_count);

    Py_END_ALLOW_THREADS;

    PyObject * trim_seq = PyUnicode_FromStringAndSize(seq, trim_at);
    if (trim_seq == NULL) {
        return NULL;
    }
    PyObject * ret = Py_BuildValue("Ok", trim_seq, trim_at);
    Py_DECREF(trim_seq);

    return ret;
}

static
PyObject *
hashtable_trim_below_abundance(khmer_KHashtable_Object * me, PyObject * args)
{
    Hashtable * hashtable = me->hashtable;

    const char * seq = NULL;
    BoundedCounterType max_count_i = 0;

    if (!PyArg_ParseTuple(args, "sH", &seq, &max_count_i)) {
        return NULL;
    }

    unsigned long trim_at;
    Py_BEGIN_ALLOW_THREADS

    BoundedCounterType max_count = max_count_i;

    trim_at = hashtable->trim_below_abundance(seq, max_count);

    Py_END_ALLOW_THREADS;

    PyObject * trim_seq = PyUnicode_FromStringAndSize(seq, trim_at);
    if (trim_seq == NULL) {
        return NULL;
    }
    PyObject * ret = Py_BuildValue("Ok", trim_seq, trim_at);
    Py_DECREF(trim_seq);

    return ret;
}

static
PyObject *
hashtable_find_spectral_error_positions(khmer_KHashtable_Object * me,
                                        PyObject * args)
{
    Hashtable * hashtable = me->hashtable;

    const char * seq = NULL;
    BoundedCounterType max_count = 0; // unsigned short int

    if (!PyArg_ParseTuple(args, "sH", &seq, &max_count)) {
        return NULL;
    }

    std::vector<unsigned int> posns;

    try {
        posns = hashtable->find_spectral_error_positions(seq, max_count);
    } catch (khmer_exception &e) {
        PyErr_SetString(PyExc_ValueError, e.what());
        return NULL;
    }

    Py_ssize_t posns_size = posns.size();

    PyObject * x = PyList_New(posns_size);
    if (x == NULL) {
        return NULL;
    }
    for (Py_ssize_t i = 0; i < posns_size; i++) {
        PyList_SET_ITEM(x, i, PyLong_FromLong(posns[i]));
    }

    return x;
}

static
PyObject *
hashtable_load(khmer_KHashtable_Object * me, PyObject * args)
{
    Hashtable * hashtable = me->hashtable;

    const char * filename = NULL;

    if (!PyArg_ParseTuple(args, "s", &filename)) {
        return NULL;
    }

    try {
        hashtable->load(filename);
    } catch (khmer_file_exception &e) {
        PyErr_SetString(PyExc_OSError, e.what());
        return NULL;
    }

    Py_RETURN_NONE;
}

static
PyObject *
hashtable_save(khmer_KHashtable_Object * me, PyObject * args)
{
    Hashtable * hashtable = me->hashtable;

    const char * filename = NULL;

    if (!PyArg_ParseTuple(args, "s", &filename)) {
        return NULL;
    }

    try {
        hashtable->save(filename);
    } catch (khmer_file_exception &e) {
        PyErr_SetString(PyExc_OSError, e.what());
        return NULL;
    }

    Py_RETURN_NONE;
}

static
PyObject *
hashtable_get_hashsizes(khmer_KHashtable_Object * me, PyObject * args)
{
    Hashtable * hashtable = me->hashtable;


    if (!PyArg_ParseTuple(args, "")) {
        return NULL;
    }

    std::vector<uint64_t> ts = hashtable->get_tablesizes();

    PyObject * x = PyList_New(ts.size());
    for (size_t i = 0; i < ts.size(); i++) {
        PyList_SET_ITEM(x, i, PyLong_FromUnsignedLongLong(ts[i]));
    }

    return x;
}

static
PyObject *
hashtable_get_median_count(khmer_KHashtable_Object * me, PyObject * args)
{
    Hashtable * hashtable = me->hashtable;

    const char * long_str;

    if (!PyArg_ParseTuple(args, "s", &long_str)) {
        return NULL;
    }

    if (strlen(long_str) < hashtable->ksize()) {
        PyErr_SetString(PyExc_ValueError,
                        "string length must >= the hashtable k-mer size");
        return NULL;
    }

    BoundedCounterType med = 0;
    float average = 0, stddev = 0;

    hashtable->get_median_count(long_str, med, average, stddev);

    return Py_BuildValue("iff", med, average, stddev);
}

static
PyObject *
hashtable_median_at_least(khmer_KHashtable_Object * me, PyObject * args)
{
    Hashtable * hashtable = me->hashtable;

    const char * long_str;
    unsigned int cutoff;

    if (!PyArg_ParseTuple(args, "sI", &long_str, &cutoff)) {
        return NULL;
    }

    if (strlen(long_str) < hashtable->ksize()) {
        PyErr_SetString(PyExc_ValueError,
                        "string length must >= the hashtable k-mer size");
        return NULL;
    }

    if (hashtable->median_at_least(long_str, cutoff)) {
        Py_RETURN_TRUE;
    }
    Py_RETURN_FALSE;

}

static
PyObject *
hashtable_get_kmers(khmer_KHashtable_Object * me, PyObject * args)
{
    Hashtable * hashtable = me->hashtable;
    const char * sequence;

    if (!PyArg_ParseTuple(args, "s", &sequence)) {
        return NULL;
    }

    std::vector<std::string> kmers;

    hashtable->get_kmers(sequence, kmers);

    PyObject * x = PyList_New(kmers.size());
    for (unsigned int i = 0; i < kmers.size(); i++) {
        PyObject * obj = PyUnicode_FromString(kmers[i].c_str());
        PyList_SET_ITEM(x, i, obj);
    }

    return x;
}

static
PyObject *
hashtable_get_kmer_counts(khmer_KHashtable_Object * me, PyObject * args)
{
    Hashtable * hashtable = me->hashtable;
    const char * sequence;

    if (!PyArg_ParseTuple(args, "s", &sequence)) {
        return NULL;
    }

    std::vector<BoundedCounterType> counts;
    try {
        hashtable->get_kmer_counts(sequence, counts);
    } catch (khmer_exception &e) {
        PyErr_SetString(PyExc_ValueError, e.what());
        return NULL;
    }

    PyObject * x = PyList_New(counts.size());
    for (unsigned int i = 0; i <counts.size(); i++) {
        PyObject * obj = PyInt_FromLong(counts[i]);
        PyList_SET_ITEM(x, i, obj);
    }

    return x;
}


static
PyObject *
hashtable_get_kmer_hashes(khmer_KHashtable_Object * me, PyObject * args)
{
    Hashtable * hashtable = me->hashtable;
    const char * sequence;

    if (!PyArg_ParseTuple(args, "s", &sequence)) {
        return NULL;
    }

    std::vector<HashIntoType> hashes;
    try {
        hashtable->get_kmer_hashes(sequence, hashes);
    } catch (khmer_exception &e) {
        PyErr_SetString(PyExc_ValueError, e.what());
        return NULL;
    }

    PyObject * x = PyList_New(hashes.size());
    for (unsigned int i = 0; i < hashes.size(); i++) {
        PyObject * obj = nullptr;
        convert_HashIntoType_to_PyObject(hashes[i], &obj);
        PyList_SET_ITEM(x, i, obj);
    }

    return x;
}


static
PyObject *
hashtable_get_kmer_hashes_as_hashset(khmer_KHashtable_Object * me,
                                     PyObject * args)
{
    Hashtable * hashtable = me->hashtable;
    const char * sequence;

    if (!PyArg_ParseTuple(args, "s", &sequence)) {
        return NULL;
    }

    SeenSet * hashes = new SeenSet;
    try {
        hashtable->get_kmer_hashes_as_hashset(sequence, *hashes);
    } catch (khmer_exception &e) {
        PyErr_SetString(PyExc_ValueError, e.what());
        return NULL;
    }

    PyObject * x = (PyObject *) create_HashSet_Object(hashes,
                   hashtable->ksize());

    return x;
}


static PyMethodDef khmer_hashtable_methods[] = {
    //
    // Basic methods
    //

    {
        "ksize",
        (PyCFunction)hashtable_ksize, METH_VARARGS,
        "Returns the k-mer size of this graph."
    },
    {
        "hash",
        (PyCFunction)hashtable_hash, METH_VARARGS,
        "Returns the hash of this k-mer. For Nodetables and Counttables, this "
        "function will fail if the supplied k-mer contains non-ACGT "
        "characters."
    },
    {
        "reverse_hash",
        (PyCFunction)hashtable_reverse_hash, METH_VARARGS,
        "Turns a k-mer hash back into a DNA k-mer, if possible."
    },
    {
        "hashsizes",
        (PyCFunction)hashtable_get_hashsizes, METH_VARARGS,
        "" },
    {
        "n_unique_kmers",
        (PyCFunction)hashtable_n_unique_kmers, METH_VARARGS,
        "Count the number of unique kmers in this graph."
    },
    {
        "n_occupied", (PyCFunction)hashtable_n_occupied, METH_VARARGS,
        "Count the number of occupied bins."
    },
    {
        "count",
        (PyCFunction)hashtable_count, METH_VARARGS,
        "Increment the count of this k-mer."
    },
    {
        "add",
        (PyCFunction)hashtable_count, METH_VARARGS,
        "Increment the count of this k-mer. (Synonym for 'count'.)"
    },
    {
        "consume",
        (PyCFunction)hashtable_consume, METH_VARARGS,
        "Increment the counts of all of the k-mers in the string."
    },
    {
        "consume_seqfile",
        (PyCFunction)hashtable_consume_seqfile, METH_VARARGS,
        "Increment the counts of all the k-mers in the sequences in the "
        "given file"
    },
    {
        "consume_seqfile_with_reads_parser",
        (PyCFunction)hashtable_consume_seqfile_with_reads_parser, METH_VARARGS,
        "Count all k-mers retrieved with this reads parser object."
    },
    {
        "get",
        (PyCFunction)hashtable_get, METH_VARARGS,
        "Retrieve the count for the given k-mer. For Nodetables and "
        "Counttables, this function will fail if the supplied k-mer contains "
        "non-ACGT characters."
    },
    {
        "load",
        (PyCFunction)hashtable_load, METH_VARARGS,
        "Load the graph from the specified file."
    },
    {
        "save",
        (PyCFunction)hashtable_save, METH_VARARGS,
        "Save the graph to the specified file."
    },
    {
        "get_kmers",
        (PyCFunction)hashtable_get_kmers, METH_VARARGS,
        "Generate an ordered list of all substrings of length k in the string."
    },
    {
        "get_kmer_hashes",
        (PyCFunction)hashtable_get_kmer_hashes, METH_VARARGS,
        "Retrieve an ordered list of all hashes of all k-mers in the string."
    },
    {
        "get_kmer_hashes_as_hashset",
        (PyCFunction)hashtable_get_kmer_hashes_as_hashset, METH_VARARGS,
        "Retrieve a HashSet containing all the k-mers in the string."
    },
    {
        "get_kmer_counts",
        (PyCFunction)hashtable_get_kmer_counts, METH_VARARGS,
        "Retrieve an ordered list of the counts of all k-mers in the string."
    },

    {
        "set_use_bigcount",
        (PyCFunction)hashtable_set_use_bigcount, METH_VARARGS,
        "Count past maximum binsize of hashtable (set to T/F)"
    },
    {
        "get_use_bigcount",
        (PyCFunction)hashtable_get_use_bigcount, METH_VARARGS,
        "Get value of bigcount flag (T/F)"
    },
    {
        "get_min_count",
        (PyCFunction)hashtable_get_min_count, METH_VARARGS,
        "Get the smallest count of all the k-mers in the string"
    },
    {
        "get_max_count",
        (PyCFunction)hashtable_get_max_count, METH_VARARGS,
        "Get the largest count of all the k-mers in the string"
    },
    {
        "trim_on_abundance",
        (PyCFunction)hashtable_trim_on_abundance, METH_VARARGS,
        "Trim string at first k-mer below the given abundance"
    },
    {
        "trim_below_abundance",
        (PyCFunction)hashtable_trim_below_abundance, METH_VARARGS,
        "Trim string at first k-mer above the given abundance"
    },
    {
        "find_spectral_error_positions",
        (PyCFunction)hashtable_find_spectral_error_positions, METH_VARARGS,
        "Identify positions of low-abundance k-mers"
    },
    {
        "abundance_distribution",
        (PyCFunction)hashtable_abundance_distribution, METH_VARARGS,
        "Calculate the k-mer abundance distribution of the given file"
    },
    {
        "abundance_distribution_with_reads_parser",
        (PyCFunction)hashtable_abundance_distribution_with_reads_parser,
        METH_VARARGS,
        "Calculate the k-mer abundance distribution for a reads parser handle"
    },
    {
        "get_median_count",
        (PyCFunction)hashtable_get_median_count, METH_VARARGS,
        "Get the median, average, and stddev of the k-mer counts in the string"
    },
    {
        "median_at_least",
        (PyCFunction)hashtable_median_at_least, METH_VARARGS,
        "Return true if the median is at least the given cutoff"
    },
    {NULL, NULL, 0, NULL}           /* sentinel */
};

static PyTypeObject khmer_KHashtable_Type
CPYCHECKER_TYPE_OBJECT_FOR_TYPEDEF("khmer_KHashtable_Object")
= {
    PyVarObject_HEAD_INIT(NULL, 0)       /* init & ob_size */
    "_khmer.KHashtable   ",              /*tp_name*/
    sizeof(khmer_KHashtable_Object),     /*tp_basicsize*/
    0,                                   /*tp_itemsize*/
    0,                                   /*tp_dealloc*/
    0,                                   /*tp_print*/
    0,                                   /*tp_getattr*/
    0,                                   /*tp_setattr*/
    0,                                   /*tp_compare*/
    0,                                   /*tp_repr*/
    0,                                   /*tp_as_number*/
    0,                                   /*tp_as_sequence*/
    0,                                   /*tp_as_mapping*/
    0,                                   /*tp_hash */
    0,                                   /*tp_call*/
    0,                                   /*tp_str*/
    0,                                   /*tp_getattro*/
    0,                                   /*tp_setattro*/
    0,                                   /*tp_as_buffer*/
    Py_TPFLAGS_DEFAULT,                  /*tp_flags*/
    "base hashtable object",             /* tp_doc */
    0,                                   /* tp_traverse */
    0,                                   /* tp_clear */
    0,                                   /* tp_richcompare */
    0,                                   /* tp_weaklistoffset */
    0,                                   /* tp_iter */
    0,                                   /* tp_iternext */
    khmer_hashtable_methods,             /* tp_methods */
    0,                                   /* tp_members */
    0,                                   /* tp_getset */
    0,                                   /* tp_base */
    0,                                   /* tp_dict */
    0,                                   /* tp_descr_get */
    0,                                   /* tp_descr_set */
    0,                                   /* tp_dictoffset */
    0,                                   /* tp_init */
    0,                                   /* tp_alloc */
    0,                                   /* tp_new */
};

#include "_cpy_nodetable.hh"
#include "_cpy_counttable.hh"
#include "_cpy_smallcounttable.hh"
#include "_cpy_hashgraph.hh"
#include "_cpy_smallcountgraph.hh"

//
// KCountgraph object
//

typedef struct {
    khmer_KHashgraph_Object khashgraph;
    Countgraph * countgraph;
} khmer_KCountgraph_Object;

typedef struct {
    PyObject_HEAD
    ReadAligner * aligner;
} khmer_ReadAligner_Object;

static void khmer_countgraph_dealloc(khmer_KCountgraph_Object * obj);

static
PyObject *
count_get_raw_tables(khmer_KCountgraph_Object * self, PyObject * args)
{
    Countgraph * countgraph = self->countgraph;

    khmer::Byte ** table_ptrs = countgraph->get_raw_tables();
    std::vector<uint64_t> sizes = countgraph->get_tablesizes();

    PyObject * raw_tables = PyList_New(sizes.size());
    for (unsigned int i=0; i<sizes.size(); ++i) {
        Py_buffer buffer;
        int res = PyBuffer_FillInfo(&buffer, NULL, table_ptrs[i], sizes[i], 0,
                                    PyBUF_FULL_RO);
        if (res == -1) {
            return NULL;
        }
        PyObject * buf = PyMemoryView_FromBuffer(&buffer);
        if(!PyMemoryView_Check(buf)) {
            return NULL;
        }
        PyList_SET_ITEM(raw_tables, i, buf);
    }

    return raw_tables;
}

static
PyObject *
count_do_subset_partition_with_abundance(khmer_KCountgraph_Object * me,
        PyObject * args)
{
    Countgraph * countgraph = me->countgraph;

    HashIntoType start_kmer = 0, end_kmer = 0;
    PyObject * break_on_stop_tags_o = NULL;
    PyObject * stop_big_traversals_o = NULL;
    BoundedCounterType min_count, max_count;

    if (!PyArg_ParseTuple(args, "HH|KKOO",
                          &min_count, &max_count,
                          &start_kmer, &end_kmer,
                          &break_on_stop_tags_o,
                          &stop_big_traversals_o)) {
        return NULL;
    }

    bool break_on_stop_tags = false;
    if (break_on_stop_tags_o && PyObject_IsTrue(break_on_stop_tags_o)) {
        break_on_stop_tags = true;
    }
    bool stop_big_traversals = false;
    if (stop_big_traversals_o && PyObject_IsTrue(stop_big_traversals_o)) {
        stop_big_traversals = true;
    }

    SubsetPartition * subset_p = NULL;
    try {
        Py_BEGIN_ALLOW_THREADS
        subset_p = new SubsetPartition(countgraph);
        subset_p->do_partition_with_abundance(start_kmer, end_kmer,
                                              min_count, max_count,
                                              break_on_stop_tags,
                                              stop_big_traversals);
        Py_END_ALLOW_THREADS
    } catch (std::bad_alloc &e) {
        return PyErr_NoMemory();
    }

    khmer_KSubsetPartition_Object * subset_obj = (khmer_KSubsetPartition_Object *)\
            PyObject_New(khmer_KSubsetPartition_Object, &khmer_KSubsetPartition_Type);

    if (subset_obj == NULL) {
        delete subset_p;
        return NULL;
    }

    subset_obj->subset = subset_p;

    return (PyObject *) subset_obj;
}

static PyMethodDef khmer_countgraph_methods[] = {
    {
        "get_raw_tables",
        (PyCFunction)count_get_raw_tables, METH_VARARGS,
        "Get a list of the raw storage tables as memoryview objects."
    },
    { "do_subset_partition_with_abundance", (PyCFunction)count_do_subset_partition_with_abundance, METH_VARARGS, "" },
    {NULL, NULL, 0, NULL}           /* sentinel */
};

static PyObject* khmer_countgraph_new(PyTypeObject * type, PyObject * args,
                                      PyObject * kwds);

static PyTypeObject khmer_KCountgraph_Type
CPYCHECKER_TYPE_OBJECT_FOR_TYPEDEF("khmer_KCountgraph_Object")
= {
    PyVarObject_HEAD_INIT(NULL, 0)       /* init & ob_size */
    "_khmer.Countgraph",                 /*tp_name*/
    sizeof(khmer_KCountgraph_Object),  /*tp_basicsize*/
    0,                                   /*tp_itemsize*/
    (destructor)khmer_countgraph_dealloc,  /*tp_dealloc*/
    0,                                   /*tp_print*/
    0,                                   /*tp_getattr*/
    0,                                   /*tp_setattr*/
    0,                                   /*tp_compare*/
    0,                                   /*tp_repr*/
    0,                                   /*tp_as_number*/
    0,                                   /*tp_as_sequence*/
    0,                                   /*tp_as_mapping*/
    0,                                   /*tp_hash */
    0,                                   /*tp_call*/
    0,                                   /*tp_str*/
    0,                                   /*tp_getattro*/
    0,                                   /*tp_setattro*/
    0,                                   /*tp_as_buffer*/
    Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE,                  /*tp_flags*/
    "countgraph hash object",              /* tp_doc */
    0,                                   /* tp_traverse */
    0,                                   /* tp_clear */
    0,                                   /* tp_richcompare */
    0,                                   /* tp_weaklistoffset */
    0,                                   /* tp_iter */
    0,                                   /* tp_iternext */
    khmer_countgraph_methods,              /* tp_methods */
    0,                                   /* tp_members */
    0,                                   /* tp_getset */
    0,                                   /* tp_base */
    0,                                   /* tp_dict */
    0,                                   /* tp_descr_get */
    0,                                   /* tp_descr_set */
    0,                                   /* tp_dictoffset */
    0,                                   /* tp_init */
    0,                                   /* tp_alloc */
    khmer_countgraph_new,                /* tp_new */
};

//
// khmer_countgraph_new
//

static PyObject* khmer_countgraph_new(PyTypeObject * type, PyObject * args,
                                      PyObject * kwds)
{
    khmer_KCountgraph_Object * self;

    self = (khmer_KCountgraph_Object *)type->tp_alloc(type, 0);

    if (self != NULL) {
        WordLength k = 0;
        PyListObject * sizes_list_o = NULL;

        if (!PyArg_ParseTuple(args, "bO!", &k, &PyList_Type, &sizes_list_o)) {
            Py_DECREF(self);
            return NULL;
        }

        std::vector<uint64_t> sizes;
        if (!convert_Pytablesizes_to_vector(sizes_list_o, sizes)) {
            Py_DECREF(self);
            return NULL;
        }

        try {
            self->countgraph = new Countgraph(k, sizes);
        } catch (std::bad_alloc &e) {
            Py_DECREF(self);
            return PyErr_NoMemory();
        }
        self->khashgraph.khashtable.hashtable =
            dynamic_cast<Hashtable*>(self->countgraph);
        self->khashgraph.hashgraph = dynamic_cast<Hashgraph*>(self->countgraph);
    }

    return (PyObject *) self;
}

static
PyObject *
nodegraph_update(khmer_KNodegraph_Object * me, PyObject * args)
{
    Nodegraph * nodegraph = me->nodegraph;
    Nodegraph * other;
    khmer_KNodegraph_Object * other_o;

    if (!PyArg_ParseTuple(args, "O!", &khmer_KNodegraph_Type, &other_o)) {
        return NULL;
    }

    other = other_o->nodegraph;

    try {
        nodegraph->update_from(*other);
    } catch (khmer_exception &e) {
        PyErr_SetString(PyExc_ValueError, e.what());
        return NULL;
    }

    Py_RETURN_NONE;
}

static
PyObject *
nodegraph_get_raw_tables(khmer_KNodegraph_Object * self, PyObject * args)
{
    Nodegraph * countgraph = self->nodegraph;

    khmer::Byte ** table_ptrs = countgraph->get_raw_tables();
    std::vector<uint64_t> sizes = countgraph->get_tablesizes();

    PyObject * raw_tables = PyList_New(sizes.size());
    for (unsigned int i=0; i<sizes.size(); ++i) {
        Py_buffer buffer;
        int res = PyBuffer_FillInfo(&buffer, NULL, table_ptrs[i], sizes[i], 0,
                                    PyBUF_FULL_RO);
        if (res == -1) {
            return NULL;
        }
        PyObject * buf = PyMemoryView_FromBuffer(&buffer);
        if(!PyMemoryView_Check(buf)) {
            return NULL;
        }
        PyList_SET_ITEM(raw_tables, i, buf);
    }

    return raw_tables;
}

static PyMethodDef khmer_nodegraph_methods[] = {
    {
        "update",
        (PyCFunction) nodegraph_update, METH_VARARGS,
        "a set update: update this nodegraph with all the entries from the other"
    },
    {
        "get_raw_tables",
        (PyCFunction) nodegraph_get_raw_tables, METH_VARARGS,
        "Get a list of the raw tables as memoryview objects"
    },
    {NULL, NULL, 0, NULL}           /* sentinel */
};

// __new__ for nodegraph; necessary for proper subclassing
// This will essentially do what the old factory function did. Unlike many __new__
// methods, we take our arguments here, because there's no "uninitialized" nodegraph
// object; we have to have k and the table sizes before creating the new objects
static PyObject* khmer_nodegraph_new(PyTypeObject * type, PyObject * args,
                                     PyObject * kwds)
{
    khmer_KNodegraph_Object * self;
    self = (khmer_KNodegraph_Object *)type->tp_alloc(type, 0);

    if (self != NULL) {
        WordLength k = 0;
        PyListObject* sizes_list_o = NULL;

        if (!PyArg_ParseTuple(args, "bO!", &k, &PyList_Type, &sizes_list_o)) {
            Py_DECREF(self);
            return NULL;
        }

        std::vector<uint64_t> sizes;
        if (!convert_Pytablesizes_to_vector(sizes_list_o, sizes)) {
            Py_DECREF(self);
            return NULL;
        }

        try {
            self->nodegraph = new Nodegraph(k, sizes);
        } catch (std::bad_alloc &e) {
            Py_DECREF(self);
            return PyErr_NoMemory();
        }
        self->khashgraph.khashtable.hashtable =
            dynamic_cast<Hashtable*>(self->nodegraph);
        self->khashgraph.hashgraph = dynamic_cast<Hashgraph*>(self->nodegraph);
    }
    return (PyObject *) self;
}

////////////////////////////////////////////////////////////////////////////

static
PyObject *
subset_count_partitions(khmer_KSubsetPartition_Object * me, PyObject * args)
{
    SubsetPartition * subset_p = me->subset;

    if (!PyArg_ParseTuple(args, "")) {
        return NULL;
    }

    size_t n_partitions = 0, n_unassigned = 0;
    subset_p->count_partitions(n_partitions, n_unassigned);

    return Py_BuildValue("nn", (Py_ssize_t) n_partitions,
                         (Py_ssize_t) n_unassigned);
}

static
PyObject *
subset_report_on_partitions(khmer_KSubsetPartition_Object * me, PyObject * args)
{
    SubsetPartition * subset_p = me->subset;

    if (!PyArg_ParseTuple(args, "")) {
        return NULL;
    }

    subset_p->report_on_partitions();

    Py_RETURN_NONE;
}

static
PyObject *
subset_partition_size_distribution(khmer_KSubsetPartition_Object * me,
                                   PyObject * args)
{
    SubsetPartition * subset_p = me->subset;

    if (!PyArg_ParseTuple(args, "")) {
        return NULL;
    }

    PartitionCountDistribution d;

    unsigned int n_unassigned = 0;
    subset_p->partition_size_distribution(d, n_unassigned);

    PyObject * x = PyList_New(d.size());
    if (x == NULL) {
        return NULL;
    }
    PartitionCountDistribution::iterator di;

    unsigned int i;
    for (i = 0, di = d.begin(); di != d.end(); ++di, i++) {
        PyObject * tup = Py_BuildValue("KK", di->first, di->second);
        if (tup != NULL) {
            PyList_SET_ITEM(x, i, tup);
        }
        Py_XDECREF(tup);
    }
    if (!(i == d.size())) {
        throw khmer_exception();
    }

    PyObject * ret = Py_BuildValue("OI", x, n_unassigned);
    Py_DECREF(x);
    return ret;
}

static
PyObject *
subset_partition_sizes(khmer_KSubsetPartition_Object * me, PyObject * args)
{
    SubsetPartition * subset_p = me->subset;

    unsigned int min_size = 0;

    if (!PyArg_ParseTuple(args, "|I", &min_size)) {
        return NULL;
    }

    PartitionCountMap cm;
    unsigned int n_unassigned = 0;
    subset_p->partition_sizes(cm, n_unassigned);

    unsigned int i = 0;
    PartitionCountMap::const_iterator mi;
    for (mi = cm.begin(); mi != cm.end(); ++mi) {
        if (mi->second >= min_size) {
            i++;
        }
    }

    PyObject * x = PyList_New(i);
    if (x == NULL) {
        return NULL;
    }

    // this should probably be a dict. @CTB
    for (i = 0, mi = cm.begin(); mi != cm.end(); ++mi) {
        if (mi->second >= min_size) {
            PyObject * tup = Py_BuildValue("II", mi->first, mi->second);
            if (tup != NULL) {
                PyList_SET_ITEM(x, i, tup);
            }
            i++;
        }
    }

    PyObject * ret = Py_BuildValue("OI", x, n_unassigned);
    Py_DECREF(x);

    return ret;
}

static
PyObject *
subset_partition_average_coverages(khmer_KSubsetPartition_Object * me,
                                   PyObject * args)
{
    SubsetPartition * subset_p = me->subset;

    khmer_KCountgraph_Object * countgraph_o;

    if (!PyArg_ParseTuple(args, "O!", &khmer_KCountgraph_Type, &countgraph_o)) {
        return NULL;
    }

    PartitionCountMap cm;
    subset_p->partition_average_coverages(cm, countgraph_o -> countgraph);

    unsigned int i;
    PartitionCountMap::iterator mi;

    PyObject * x = PyList_New(cm.size());
    if (x == NULL) {
        return NULL;
    }

    // this should probably be a dict. @CTB
    for (i = 0, mi = cm.begin(); mi != cm.end(); ++mi, i++) {
        PyObject * tup = Py_BuildValue("II", mi->first, mi->second);
        if (tup != NULL) {
            PyList_SET_ITEM(x, i, tup);
        }
    }

    return x;
}

static PyMethodDef khmer_subset_methods[] = {
    {
        "count_partitions",
        (PyCFunction)subset_count_partitions,
        METH_VARARGS,
        ""
    },
    {
        "report_on_partitions",
        (PyCFunction)subset_report_on_partitions,
        METH_VARARGS,
        ""
    },
    {
        "partition_size_distribution",
        (PyCFunction)subset_partition_size_distribution,
        METH_VARARGS,
        ""
    },
    {
        "partition_sizes",
        (PyCFunction)subset_partition_sizes,
        METH_VARARGS,
        ""
    },
    {
        "partition_average_coverages",
        (PyCFunction)subset_partition_average_coverages,
        METH_VARARGS,
        ""
    },
    {NULL, NULL, 0, NULL}           /* sentinel */
};

typedef struct {
    PyObject_HEAD
    LabelHash * labelhash;
} khmer_KGraphLabels_Object;

static PyObject * khmer_graphlabels_new(PyTypeObject * type, PyObject *args,
                                        PyObject *kwds);

static void khmer_graphlabels_dealloc(khmer_KGraphLabels_Object * obj)
{
    delete obj->labelhash;
    obj->labelhash = NULL;

    Py_TYPE(obj)->tp_free((PyObject*)obj);
}

static PyObject * khmer_graphlabels_new(PyTypeObject *type, PyObject *args,
                                        PyObject *kwds)
{
    khmer_KGraphLabels_Object *self;
    self = (khmer_KGraphLabels_Object*)type->tp_alloc(type, 0);

    if (self != NULL) {
        PyObject * hashgraph_o;
        khmer::Hashgraph * hashgraph = NULL; // @CTB

        if (!PyArg_ParseTuple(args, "O", &hashgraph_o)) {
            Py_DECREF(self);
            return NULL;
        }

        if (PyObject_TypeCheck(hashgraph_o, &khmer_KNodegraph_Type)) {
            khmer_KNodegraph_Object * kho = (khmer_KNodegraph_Object *) hashgraph_o;
            hashgraph = kho->nodegraph;
        } else if (PyObject_TypeCheck(hashgraph_o, &khmer_KCountgraph_Type)) {
            khmer_KCountgraph_Object * cho = (khmer_KCountgraph_Object *) hashgraph_o;
            hashgraph = cho->countgraph;
        } else {
            PyErr_SetString(PyExc_ValueError,
                            "graph object must be a NodeGraph or CountGraph");
            Py_DECREF(self);
            return NULL;
        }

        try {
            self->labelhash = new LabelHash(hashgraph);
        } catch (std::bad_alloc &e) {
            Py_DECREF(self);
            return PyErr_NoMemory();
        }
    }

    return (PyObject *) self;
}

static
PyObject *
labelhash_get_all_labels(khmer_KGraphLabels_Object * me, PyObject * args)
{
    LabelHash * hb = me->labelhash;

    PyObject * d = PyList_New(hb->all_labels.size());
    if (d == NULL) {
        return NULL;
    }
    LabelSet::iterator it;

    unsigned long long i = 0;
    for (it = hb->all_labels.begin(); it != hb->all_labels.end(); ++it) {
        PyObject * val = Py_BuildValue("K", *it);
        if (val != NULL) {
            PyList_SetItem(d, i, val);
        }
        i++;
    }

    return d;
}

static
PyObject *
labelhash_consume_seqfile_and_tag_with_labels(khmer_KGraphLabels_Object * me,
        PyObject * args)
{
    LabelHash * hb = me->labelhash;

    const char * filename;

    if (!PyArg_ParseTuple(args, "s", &filename)) {
        return NULL;
    }

    const char         *value_exception = NULL;
    const char         *file_exception  = NULL;
    unsigned long long  n_consumed      = 0;
    unsigned int        total_reads     = 0;
    std::string exc_string;

    //Py_BEGIN_ALLOW_THREADS
    try {
        hb->consume_seqfile_and_tag_with_labels<FastxReader>(filename, total_reads,
                n_consumed);
    } catch (khmer_file_exception &exc) {
        exc_string = exc.what();
        file_exception = exc_string.c_str();
    } catch (khmer_value_exception &exc) {
        exc_string = exc.what();
        value_exception = exc_string.c_str();
    }
    //Py_END_ALLOW_THREADS

    if (file_exception != NULL) {
        PyErr_SetString(PyExc_OSError, file_exception);
        return NULL;
    }
    if (value_exception != NULL) {
        PyErr_SetString(PyExc_ValueError, value_exception);
        return NULL;
    }

    return Py_BuildValue("IK", total_reads, n_consumed);
}

static
PyObject *
labelhash_consume_partitioned_fasta_and_tag_with_labels(
    khmer_KGraphLabels_Object * me, PyObject * args)
{
    LabelHash * labelhash = me->labelhash;

    const char * filename;

    if (!PyArg_ParseTuple(args, "s", &filename)) {
        return NULL;
    }

    // call the C++ function, and trap signals => Python

    unsigned long long  n_consumed  = 0;
    unsigned int        total_reads = 0;

    try {
        labelhash->consume_partitioned_fasta_and_tag_with_labels<FastxReader>(filename,
                total_reads, n_consumed);
    } catch (khmer_file_exception &exc) {
        PyErr_SetString(PyExc_OSError, exc.what());
        return NULL;
    } catch (khmer_value_exception &exc) {
        PyErr_SetString(PyExc_ValueError, exc.what());
        return NULL;
    }

    return Py_BuildValue("IK", total_reads, n_consumed);
}

static
PyObject *
labelhash_consume_sequence_and_tag_with_labels(khmer_KGraphLabels_Object * me,
        PyObject * args)
{
    LabelHash * hb = me->labelhash;
    const char * seq = NULL;
    unsigned long long c = 0;
    if (!PyArg_ParseTuple(args, "sK", &seq, &c)) {
        return NULL;
    }
    unsigned long long n_consumed = 0;

    hb->consume_sequence_and_tag_with_labels(seq, n_consumed, c);
    return Py_BuildValue("K", n_consumed);
}

static
PyObject *
labelhash_sweep_label_neighborhood(khmer_KGraphLabels_Object * me,
                                   PyObject * args)
{
    LabelHash * hb = me->labelhash;

    const char * seq = NULL;
    int r = 0;
    PyObject * break_on_stop_tags_o = NULL;
    PyObject * stop_big_traversals_o = NULL;

    if (!PyArg_ParseTuple(args, "s|iOO", &seq, &r,
                          &break_on_stop_tags_o,
                          &stop_big_traversals_o)) {
        return NULL;
    }

    unsigned int range = (2 * hb->graph->_get_tag_density()) + 1;
    if (r >= 0) {
        range = r;
    }

    bool break_on_stop_tags = false;
    if (break_on_stop_tags_o && PyObject_IsTrue(break_on_stop_tags_o)) {
        break_on_stop_tags = true;
    }
    bool stop_big_traversals = false;
    if (stop_big_traversals_o && PyObject_IsTrue(stop_big_traversals_o)) {
        stop_big_traversals = true;
    }

    if (strlen(seq) < hb->graph->ksize()) {
        PyErr_SetString(PyExc_ValueError,
                        "string length must >= the hashtable k-mer size");
        return NULL;
    }

    //std::pair<TagLabelPair::iterator, TagLabelPair::iterator> ret;
    LabelSet found_labels;

    //unsigned int num_traversed = 0;
    //Py_BEGIN_ALLOW_THREADS
    hb->sweep_label_neighborhood(seq, found_labels, range, break_on_stop_tags,
                                 stop_big_traversals);
    //Py_END_ALLOW_THREADS

    //printf("...%u kmers traversed\n", num_traversed);

    PyObject * x =  PyList_New(found_labels.size());
    LabelSet::const_iterator si;
    unsigned long long i = 0;
    for (si = found_labels.begin(); si != found_labels.end(); ++si) {
        PyList_SET_ITEM(x, i, Py_BuildValue("K", *si));
        i++;
    }

    return x;
}

// Similar to find_all_tags, but returns tags in a way actually usable by python
// need a tags_in_sequence iterator or function in c++ land for reuse in all
// these functions

static
PyObject *
labelhash_sweep_tag_neighborhood(khmer_KGraphLabels_Object * me,
                                 PyObject * args)
{
    LabelHash * labelhash = me->labelhash;

    const char * seq = NULL;
    int r = 0;
    PyObject * break_on_stop_tags_o = NULL;
    PyObject * stop_big_traversals_o = NULL;

    if (!PyArg_ParseTuple(args, "s|iOO", &seq, &r,
                          &break_on_stop_tags_o,
                          &stop_big_traversals_o)) {
        return NULL;
    }

    unsigned int range = (2 * labelhash->graph->_get_tag_density()) + 1;
    if (r >= 0) {
        range = r;
    }

    bool break_on_stop_tags = false;
    if (break_on_stop_tags_o && PyObject_IsTrue(break_on_stop_tags_o)) {
        break_on_stop_tags = true;
    }
    bool stop_big_traversals = false;
    if (stop_big_traversals_o && PyObject_IsTrue(stop_big_traversals_o)) {
        stop_big_traversals = true;
    }

    if (strlen(seq) < labelhash->graph->ksize()) {
        PyErr_SetString(PyExc_ValueError,
                        "string length must >= the hashtable k-mer size");
        return NULL;
    }

    SeenSet * tagged_kmers = new SeenSet;

    //Py_BEGIN_ALLOW_THREADS

    labelhash->graph->partition->sweep_for_tags(seq, *tagged_kmers,
            labelhash->graph->all_tags,
            range, break_on_stop_tags,
            stop_big_traversals);

    //Py_END_ALLOW_THREADS

    PyObject * x = (PyObject *) create_HashSet_Object(tagged_kmers,
                   labelhash->graph->ksize());
    return x;
}

static
PyObject *
labelhash_get_tag_labels(khmer_KGraphLabels_Object * me, PyObject * args)
{
    LabelHash * labelhash = me->labelhash;

    PyObject * tag_o;
    HashIntoType tag;

    if (!PyArg_ParseTuple(args, "O", &tag_o)) {
        return NULL;
    }
    if (!ht_convert_PyObject_to_HashIntoType(tag_o, tag,
            labelhash->graph)) {
        return NULL;
    }

    LabelSet labels;
    labelhash->get_tag_labels(tag, labels);

    PyObject * x =  PyList_New(labels.size());
    LabelSet::const_iterator si;
    unsigned long long i = 0;
    for (si = labels.begin(); si != labels.end(); ++si) {
        PyList_SET_ITEM(x, i, Py_BuildValue("K", *si));
        i++;
    }

    return x;
}

static
PyObject *
labelhash_n_labels(khmer_KGraphLabels_Object * me, PyObject * args)
{
    LabelHash * labelhash = me->labelhash;

    if (!PyArg_ParseTuple(args, "")) {
        return NULL;
    }

    return PyLong_FromSize_t(labelhash->n_labels());
}

static
PyObject *
labelhash_label_across_high_degree_nodes(khmer_KGraphLabels_Object * me,
        PyObject * args)
{
    LabelHash * labelhash = me->labelhash;

    const char * long_str;
    khmer_HashSet_Object * hdn_o = NULL;
    Label label;

    if (!PyArg_ParseTuple(args, "sO!K", &long_str,
                          &khmer_HashSet_Type, &hdn_o, &label)) {
        return NULL;
    }

    if (strlen(long_str) < labelhash->graph->ksize()) {
        Py_INCREF(Py_None);
        return Py_None;
    }

    labelhash->label_across_high_degree_nodes(long_str, *hdn_o->hashes, label);

    Py_INCREF(Py_None);
    return Py_None;
}

static
PyObject *
labelhash_assemble_labeled_path(khmer_KGraphLabels_Object * me,
                                PyObject * args)
{
    LabelHash* labelhash = me->labelhash;

    PyObject * val_o;
    khmer_KNodegraph_Object * nodegraph_o = NULL;
    Nodegraph * stop_bf = NULL;

    if (!PyArg_ParseTuple(args, "O|O!", &val_o,
                          &khmer_KNodegraph_Type, &nodegraph_o)) {
        return NULL;
    }

    Kmer start_kmer;
    if (!ht_convert_PyObject_to_Kmer(val_o, start_kmer, labelhash->graph)) {
        return NULL;
    }

    if (nodegraph_o) {
        stop_bf = nodegraph_o->nodegraph;
    }

    SimpleLabeledAssembler assembler(labelhash);
    std::vector<std::string> contigs = assembler.assemble(start_kmer, stop_bf);

    PyObject * ret = PyList_New(contigs.size());
    for (unsigned int i = 0; i < contigs.size(); i++) {
        PyList_SET_ITEM(ret, i, PyUnicode_FromString(contigs[i].c_str()));
    }

    return ret;
}

static
PyObject *
labelhash_save_labels_and_tags(khmer_KGraphLabels_Object * me, PyObject * args)
{
    const char * filename = NULL;
    LabelHash * labelhash = me->labelhash;

    if (!PyArg_ParseTuple(args, "s", &filename)) {
        return NULL;
    }

    try {
        labelhash->save_labels_and_tags(filename);
    } catch (khmer_file_exception &e) {
        PyErr_SetString(PyExc_OSError, e.what());
        return NULL;
    }

    Py_RETURN_NONE;
}

static
PyObject *
labelhash_load_labels_and_tags(khmer_KGraphLabels_Object * me, PyObject * args)
{
    const char * filename = NULL;
    LabelHash * labelhash = me->labelhash;

    if (!PyArg_ParseTuple(args, "s", &filename)) {
        return NULL;
    }

    try {
        labelhash->load_labels_and_tags(filename);
    } catch (khmer_file_exception &e) {
        PyErr_SetString(PyExc_OSError, e.what());
        return NULL;
    }

    Py_RETURN_NONE;
}

static PyMethodDef khmer_graphlabels_methods[] = {
    { "consume_seqfile_and_tag_with_labels", (PyCFunction)labelhash_consume_seqfile_and_tag_with_labels, METH_VARARGS, "" },
    { "sweep_label_neighborhood", (PyCFunction)labelhash_sweep_label_neighborhood, METH_VARARGS, "" },
    {"consume_partitioned_fasta_and_tag_with_labels", (PyCFunction)labelhash_consume_partitioned_fasta_and_tag_with_labels, METH_VARARGS, "" },
    {"sweep_tag_neighborhood", (PyCFunction)labelhash_sweep_tag_neighborhood, METH_VARARGS, "" },
    {"get_tag_labels", (PyCFunction)labelhash_get_tag_labels, METH_VARARGS, ""},
    {"consume_sequence_and_tag_with_labels", (PyCFunction)labelhash_consume_sequence_and_tag_with_labels, METH_VARARGS, "" },
    {"n_labels", (PyCFunction)labelhash_n_labels, METH_VARARGS, ""},
    {"get_all_labels", (PyCFunction)labelhash_get_all_labels, METH_VARARGS, "" },
    {
        "label_across_high_degree_nodes",
        (PyCFunction)labelhash_label_across_high_degree_nodes, METH_VARARGS,
        "Connect graph across high degree nodes using labels.",
    },
    {
        "assemble_labeled_path",
        (PyCFunction)labelhash_assemble_labeled_path, METH_VARARGS,
        "Assemble all paths, using labels to negotiate tricky bits."
    },
    { "save_labels_and_tags", (PyCFunction)labelhash_save_labels_and_tags, METH_VARARGS, "" },
    { "load_labels_and_tags", (PyCFunction)labelhash_load_labels_and_tags, METH_VARARGS, "" },    {NULL, NULL, 0, NULL}           /* sentinel */
};

static PyTypeObject khmer_KGraphLabels_Type = {
    PyVarObject_HEAD_INIT(NULL, 0)  /* init & ob_size */
    "_khmer.LabelHash",            /* tp_name */
    sizeof(khmer_KGraphLabels_Object), /* tp_basicsize */
    0,                       /* tp_itemsize */
    (destructor)khmer_graphlabels_dealloc, /* tp_dealloc */
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
    Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE,   /* tp_flags */
    0,                       /* tp_doc */
    0,                       /* tp_traverse */
    0,                       /* tp_clear */
    0,                       /* tp_richcompare */
    0,                       /* tp_weaklistoffset */
    0,                       /* tp_iter */
    0,                       /* tp_iternext */
    khmer_graphlabels_methods, /* tp_methods */
    0,                       /* tp_members */
    0,                       /* tp_getset */
    0,                       /* tp_base */
    0,                       /* tp_dict */
    0,                       /* tp_descr_get */
    0,                       /* tp_descr_set */
    0,                       /* tp_dictoffset */
    0,                       /* tp_init */
    0,                       /* tp_alloc */
    khmer_graphlabels_new,      /* tp_new */
};

static
PyObject *
hashgraph_repartition_largest_partition(khmer_KHashgraph_Object * me,
                                        PyObject * args)
{
    Hashgraph * hashgraph = me->hashgraph;
    khmer_KCountgraph_Object * countgraph_o = NULL;
    PyObject * subset_o = NULL;
    SubsetPartition * subset_p;
    unsigned int distance, threshold, frequency;

    if (!PyArg_ParseTuple(args, "OO!III",
                          &subset_o,
                          &khmer_KCountgraph_Type, &countgraph_o,
                          &distance, &threshold, &frequency)) {
        return NULL;
    }

    if (PyObject_TypeCheck(subset_o, &khmer_KSubsetPartition_Type)) {
        subset_p = ((khmer_KSubsetPartition_Object *) subset_o)->subset;
    } else {
        subset_p = hashgraph->partition;
    }

    Countgraph * countgraph = countgraph_o->countgraph;

    unsigned long next_largest;
    try {
        next_largest = subset_p->repartition_largest_partition(distance,
                       threshold, frequency, *countgraph);
    } catch (khmer_exception &e) {
        PyErr_SetString(PyExc_RuntimeError, e.what());
        return NULL;
    }

    return PyLong_FromLong(next_largest);
}

static PyObject * readaligner_align(khmer_ReadAligner_Object * me,
                                    PyObject * args)
{
    const char * read;

    if (!PyArg_ParseTuple(args, "s", &read)) {
        return NULL;
    }

    /*if (strlen(read) < (unsigned int)aligner->ksize()) {
        PyErr_SetString(PyExc_ValueError,
                        "string length must >= the hashtable k-mer size");
        return NULL;
    }*/

    Alignment * aln = me->aligner->Align(read);

    const char* alignment = aln->graph_alignment.c_str();
    const char* readAlignment = aln->read_alignment.c_str();
    PyObject * ret = Py_BuildValue("dssO", aln->score, alignment,
                                   readAlignment, (aln->truncated)? Py_True : Py_False);
    delete aln;

    return ret;
}

static PyObject * readaligner_align_forward(khmer_ReadAligner_Object * me,
        PyObject * args)
{
    ReadAligner * aligner = me->aligner;

    const char * read;

    if (!PyArg_ParseTuple(args, "s", &read)) {
        return NULL;
    }

    /*if (strlen(read) < (unsigned int)aligner->ksize()) {
        PyErr_SetString(PyExc_ValueError,
                        "string length must >= the hashtable k-mer size");
        return NULL;
    }*/

    Alignment * aln;
    aln = aligner->AlignForward(read);

    const char* alignment = aln->graph_alignment.c_str();
    const char* readAlignment = aln->read_alignment.c_str();
    PyObject * x = PyList_New(aln->covs.size());
    for (size_t i = 0; i < aln->covs.size(); i++ ) {
        PyList_SET_ITEM(x, i, PyLong_FromLong(aln->covs[i]));
    }

    PyObject * ret = Py_BuildValue("dssOO", aln->score, alignment,
                                   readAlignment,
                                   (aln->truncated)? Py_True : Py_False,
                                   x);
    delete aln;
    Py_DECREF(x);

    return ret;
}

static PyObject* khmer_ReadAligner_get_scoring_matrix(
    khmer_ReadAligner_Object * me, PyObject * args)
{

    if (!PyArg_ParseTuple(args, "")) {
        return NULL;
    }
    ScoringMatrix matrix = me->aligner->getScoringMatrix();

    return Py_BuildValue( "dddd", matrix.trusted_match, matrix.trusted_mismatch,
                          matrix.untrusted_match, matrix.untrusted_mismatch);
}

static PyObject* khmer_ReadAligner_get_transition_probabilities(
    khmer_ReadAligner_Object * me, PyObject * args)
{

    if (!PyArg_ParseTuple(args, "")) {
        return NULL;
    }
    ScoringMatrix matrix = me->aligner->getScoringMatrix();

    return Py_BuildValue( "(dddddd)(dddd)(dddd)(dddddd)(dddd)(dddd)",
                          matrix.tsc[0], matrix.tsc[1], matrix.tsc[2],
                          matrix.tsc[3], matrix.tsc[4], matrix.tsc[5],
                          matrix.tsc[6], matrix.tsc[7], matrix.tsc[8],
                          matrix.tsc[9], matrix.tsc[10], matrix.tsc[11],
                          matrix.tsc[12], matrix.tsc[13], matrix.tsc[14],
                          matrix.tsc[15], matrix.tsc[16], matrix.tsc[17],
                          matrix.tsc[18], matrix.tsc[19], matrix.tsc[20],
                          matrix.tsc[21], matrix.tsc[22], matrix.tsc[23],
                          matrix.tsc[24], matrix.tsc[25], matrix.tsc[26],
                          matrix.tsc[27]);
}

static PyMethodDef khmer_ReadAligner_methods[] = {
    {"align", (PyCFunction)readaligner_align, METH_VARARGS, ""},
    {"align_forward", (PyCFunction)readaligner_align_forward, METH_VARARGS, ""},
    {
        "get_scoring_matrix", (PyCFunction)khmer_ReadAligner_get_scoring_matrix,
        METH_VARARGS,
        "Get the scoring matrix in use.\n\n\
Returns a tuple of floats: (trusted_match, trusted_mismatch, untrusted_match, \
untrusted_mismatch)"
    },
    {
        "get_transition_probabilities",
        (PyCFunction)khmer_ReadAligner_get_transition_probabilities,
        METH_VARARGS,
        "Get the transition probabilties in use.\n\n\
HMM state notation abbreviations:\n\
    M_t - trusted match; M_u - untrusted match\n\
    Ir_t - trusted read insert; Ir_u - untrusted read insert\n\
    Ig_t - trusted graph insert; Ig_u - untrusted graph insert\n\
\
Returns a sparse matrix as a tuple of six tuples.\n\
The inner tuples contain 6, 4, 4, 6, 4, and 4 floats respectively.\n\
Transition are notated as 'StartState-NextState':\n\
(\n\
  ( M_t-M_t,  M_t-Ir_t,  M_t-Ig_t,  M_t-M_u,  M_t-Ir_u,  M_t-Ig_u),\n\
  (Ir_t-M_t, Ir_t-Ir_t,            Ir_t-M_u, Ir_t-Ir_u           ),\n\
  (Ig_t-M_t,          , Ig_t-Ig_t, Ig_t-M_u,            Ig_t-Ig_u),\n\
  ( M_u-M_t,  M_u-Ir_t,  M_u-Ig_t,  M_u-M_u,  M_u-Ir_u,  M_u-Ig_u),\n\
  (Ir_u-M_t, Ir_u-Ir_t,            Ir_u-M_u, Ir_u-Ir_u           ),\n\
  (Ig_u-M_t,          , Ig_u-Ig_t, Ig_u-M_u,            Ig_u-Ig_u)\n\
)"
    },
    {NULL} /* Sentinel */
};

//
// khmer_readaligner_dealloc -- clean up readaligner object
// GRAPHALIGN addition
//
static void khmer_readaligner_dealloc(khmer_ReadAligner_Object* obj)
{
    delete obj->aligner;
    obj->aligner = NULL;
    Py_TYPE(obj)->tp_free((PyObject*)obj);
}

//
// new_readaligner
//
static PyObject* khmer_ReadAligner_new(PyTypeObject *type, PyObject * args,
                                       PyObject *kwds)
{
    khmer_ReadAligner_Object * self;

    self = (khmer_ReadAligner_Object *)type->tp_alloc(type, 0);

    if (self != NULL) {
        khmer_KCountgraph_Object * ch = NULL;
        unsigned short int trusted_cov_cutoff = 2;
        double bits_theta = 1;
        double scoring_matrix[] = { 0, 0, 0, 0 };
        double * transitions = new double[28];

        if(!PyArg_ParseTuple(
                    args,
                    "O!Hd|(dddd)((dddddd)(dddd)(dddd)(dddddd)(dddd)(dddd))",
                    &khmer_KCountgraph_Type, &ch, &trusted_cov_cutoff,
                    &bits_theta, &scoring_matrix[0], &scoring_matrix[1],
                    &scoring_matrix[2], &scoring_matrix[3], &transitions[0],
                    &transitions[1], &transitions[2], &transitions[3],
                    &transitions[4], &transitions[5], &transitions[6],
                    &transitions[7], &transitions[8], &transitions[9],
                    &transitions[10], &transitions[11], &transitions[12],
                    &transitions[13], &transitions[14], &transitions[15],
                    &transitions[16], &transitions[17], &transitions[18],
                    &transitions[19], &transitions[20], &transitions[21],
                    &transitions[22], &transitions[23], &transitions[24],
                    &transitions[25], &transitions[26], &transitions[27])) {
            Py_DECREF(self);
            return NULL;
        }

        self->aligner = new ReadAligner(ch->countgraph, trusted_cov_cutoff,
                                        bits_theta, scoring_matrix,
                                        transitions);
    }

    return (PyObject *) self;
}

static PyTypeObject khmer_ReadAlignerType = {
    PyVarObject_HEAD_INIT(NULL, 0) /* init & ob_size */
    "_khmer.ReadAligner",		    /*tp_name*/
    sizeof(khmer_ReadAligner_Object),	    /*tp_basicsize*/
    0,					    /*tp_itemsize*/
    (destructor)khmer_readaligner_dealloc,  /*tp_dealloc*/
    0,                          /*tp_print*/
    0,                          /*tp_getattr*/
    0,                          /*tp_setattr*/
    0,                          /*tp_compare*/
    0,                          /*tp_repr*/
    0,                          /*tp_as_number*/
    0,                          /*tp_as_sequence*/
    0,                          /*tp_as_mapping*/
    0,                          /*tp_hash */
    0,                          /*tp_call*/
    0,                          /*tp_str*/
    0,                          /*tp_getattro*/
    0,                          /*tp_setattro*/
    0,                          /*tp_as_buffer*/
    Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE,         /*tp_flags*/
    "ReadAligner object",           /* tp_doc */
    0,                         /* tp_traverse */
    0,                         /* tp_clear */
    0,                         /* tp_richcompare */
    0,                         /* tp_weaklistoffset */
    0,                         /* tp_iter */
    0,                         /* tp_iternext */
    khmer_ReadAligner_methods, /* tp_methods */
    0,                         /* tp_members */
    0,                         /* tp_getset */
    0,                         /* tp_base */
    0,                         /* tp_dict */
    0,                         /* tp_descr_get */
    0,                         /* tp_descr_set */
    0,                         /* tp_dictoffset */
    0,			               /* tp_init */
    0,                         /* tp_alloc */
    khmer_ReadAligner_new,     /* tp_new */
};

//
// khmer_countgraph_dealloc -- clean up a countgraph hash object.
//

static void khmer_countgraph_dealloc(khmer_KCountgraph_Object * obj)
{
    delete obj->countgraph;
    obj->countgraph = NULL;
    Py_TYPE(obj)->tp_free((PyObject*)obj);
}

//
// khmer_nodegraph_dealloc -- clean up a nodegraph object.
//
static void khmer_nodegraph_dealloc(khmer_KNodegraph_Object * obj)
{
    delete obj->nodegraph;
    obj->nodegraph = NULL;

    Py_TYPE(obj)->tp_free((PyObject*)obj);
}


//
// khmer_subset_dealloc -- clean up a subset object.
//

static void khmer_subset_dealloc(khmer_KSubsetPartition_Object * obj)
{
    delete obj->subset;
    obj->subset = NULL;
    Py_TYPE(obj)->tp_free((PyObject*)obj);
}


/***********************************************************************/

//
// KHLLCounter object
//

typedef struct {
    PyObject_HEAD
    HLLCounter * hllcounter;
} khmer_KHLLCounter_Object;

static PyObject* khmer_hllcounter_new(PyTypeObject * type, PyObject * args,
                                      PyObject * kwds)
{
    khmer_KHLLCounter_Object * self;
    self = (khmer_KHLLCounter_Object *)type->tp_alloc(type, 0);

    if (self != NULL) {
        double error_rate = 0.01;
        WordLength ksize = 20;

        if (!PyArg_ParseTuple(args, "|db", &error_rate, &ksize)) {
            Py_DECREF(self);
            return NULL;
        }

        try {
            self->hllcounter = new HLLCounter(error_rate, ksize);
        } catch (InvalidValue &e) {
            Py_DECREF(self);
            PyErr_SetString(PyExc_ValueError, e.what());
            return NULL;
        }
    }

    return (PyObject *) self;
}

//
// khmer_hllcounter_dealloc -- clean up a hllcounter object.
//

static void khmer_hllcounter_dealloc(khmer_KHLLCounter_Object * obj)
{
    delete obj->hllcounter;
    obj->hllcounter = NULL;

    Py_TYPE(obj)->tp_free((PyObject*)obj);
}

static
PyObject *
hllcounter_add(khmer_KHLLCounter_Object * me, PyObject * args)
{
    const char * kmer_str;

    if (!PyArg_ParseTuple(args, "s", &kmer_str)) {
        return NULL;
    }

    try {
        me->hllcounter->add(kmer_str);
    } catch (khmer_exception &e) {
        PyErr_SetString(PyExc_ValueError, e.what());
        return NULL;
    }

    Py_RETURN_NONE;
}

static
PyObject *
hllcounter_estimate_cardinality(khmer_KHLLCounter_Object * me, PyObject * args)
{
    if (!PyArg_ParseTuple( args, "" )) {
        return NULL;
    }

    return PyLong_FromLong(me->hllcounter->estimate_cardinality());
}

static
PyObject *
hllcounter_consume_string(khmer_KHLLCounter_Object * me, PyObject * args)
{
    const char * kmer_str;
    unsigned long long n_consumed;

    if (!PyArg_ParseTuple(args, "s", &kmer_str)) {
        return NULL;
    }

    try {
        n_consumed = me->hllcounter->consume_string(kmer_str);
    } catch (khmer_exception &e) {
        PyErr_SetString(PyExc_ValueError, e.what());
        return NULL;
    }

    return PyLong_FromLong(n_consumed);
}

static PyObject * hllcounter_consume_seqfile(khmer_KHLLCounter_Object * me,
        PyObject * args, PyObject * kwds)
{
    const char * filename;
    PyObject * stream_records_o = NULL;

    static const char* const_kwlist[] = {"filename", "stream_records", NULL};
    static char** kwlist = const_cast<char**>(const_kwlist);

    bool stream_records = false;

    if (!PyArg_ParseTupleAndKeywords(args, kwds, "s|O", kwlist,
                                     &filename, &stream_records_o)) {
        return NULL;
    }

    if (stream_records_o != NULL && PyObject_IsTrue(stream_records_o)) {
        stream_records = true;
    }

    // call the C++ function, and trap signals => Python
    unsigned long long  n_consumed    = 0;
    unsigned int        total_reads   = 0;
    try {
        me->hllcounter->consume_seqfile<FastxReader>(filename, stream_records,
                total_reads,
                n_consumed);
    } catch (khmer_file_exception &exc) {
        PyErr_SetString(PyExc_OSError, exc.what());
        return NULL;
    } catch (khmer_value_exception &exc) {
        PyErr_SetString(PyExc_ValueError, exc.what());
        return NULL;
    }

    return Py_BuildValue("IK", total_reads, n_consumed);
}

static PyObject * hllcounter_merge(khmer_KHLLCounter_Object * me,
                                   PyObject * args);

static
PyObject *
hllcounter_get_erate(khmer_KHLLCounter_Object * me)
{
    return PyFloat_FromDouble(me->hllcounter->get_erate());
}

static
PyObject *
hllcounter_get_ksize(khmer_KHLLCounter_Object * me)
{
    return PyLong_FromLong(me->hllcounter->get_ksize());
}

static
int
hllcounter_set_ksize(khmer_KHLLCounter_Object * me, PyObject *value,
                     void *closure)
{
    if (value == NULL) {
        PyErr_SetString(PyExc_TypeError, "Cannot delete attribute");
        return -1;
    }

    long ksize = 0;
    if (PyLong_Check(value)) {
        ksize = PyLong_AsLong(value);
    } else if (PyInt_Check(value)) {
        ksize = PyInt_AsLong(value);
    } else {
        PyErr_SetString(PyExc_TypeError,
                        "Please use an integer value for k-mer size");
        return -1;
    }

    if (ksize <= 0) {
        PyErr_SetString(PyExc_ValueError, "Please set k-mer size to a value "
                        "greater than zero");
        return -1;
    }

    try {
        me->hllcounter->set_ksize(ksize);
    } catch (ReadOnlyAttribute &e) {
        PyErr_SetString(PyExc_AttributeError, e.what());
        return -1;
    }

    return 0;
}

static
int
hllcounter_set_erate(khmer_KHLLCounter_Object * me, PyObject *value,
                     void *closure)
{
    if (value == NULL) {
        PyErr_SetString(PyExc_TypeError, "Cannot delete attribute");
        return -1;
    }

    if (!PyFloat_Check(value)) {
        PyErr_SetString(PyExc_TypeError,
                        "Please use a float value for k-mer size");
        return -1;
    }

    double erate = PyFloat_AsDouble(value);
    try {
        me->hllcounter->set_erate(erate);
    } catch (InvalidValue &e) {
        PyErr_SetString(PyExc_ValueError, e.what());
        return -1;
    } catch (ReadOnlyAttribute &e) {
        PyErr_SetString(PyExc_AttributeError, e.what());
        return -1;
    }

    return 0;
}

static
PyObject *
hllcounter_getalpha(khmer_KHLLCounter_Object * me)
{
    return PyFloat_FromDouble(me->hllcounter->get_alpha());
}

static
PyObject *
hllcounter_getcounters(khmer_KHLLCounter_Object * me)
{
    std::vector<int> counters = me->hllcounter->get_M();

    PyObject * x = PyList_New(counters.size());
    for (size_t i = 0; i < counters.size(); i++) {
        PyList_SET_ITEM(x, i, PyLong_FromLong(counters[i]));
    }

    return x;
}

static PyMethodDef khmer_hllcounter_methods[] = {
    {
        "add", (PyCFunction)hllcounter_add,
        METH_VARARGS,
        "Add a k-mer to the counter."
    },
    {
        "estimate_cardinality", (PyCFunction)hllcounter_estimate_cardinality,
        METH_VARARGS,
        "Return the current estimation."
    },
    {
        "consume_string", (PyCFunction)hllcounter_consume_string,
        METH_VARARGS,
        "Break a sequence into k-mers and add each k-mer to the counter."
    },
    {
        "consume_seqfile", (PyCFunction)hllcounter_consume_seqfile,
        METH_VARARGS | METH_KEYWORDS,
        "Read sequences from file, break into k-mers, "
        "and add each k-mer to the counter. If optional keyword 'stream_out' "
        "is True, also prints each sequence to stdout."
    },
    {
        "merge", (PyCFunction)hllcounter_merge,
        METH_VARARGS,
        "Merge other counter into this one."
    },
    {NULL} /* Sentinel */
};

static PyGetSetDef khmer_hllcounter_getseters[] = {
    {
        (char *)"alpha",
        (getter)hllcounter_getalpha, NULL,
        (char *)"alpha constant for this HLL counter.",
        NULL
    },
    {
        (char *)"error_rate",
        (getter)hllcounter_get_erate, (setter)hllcounter_set_erate,
        (char *)"Error rate for this HLL counter. "
        "Can be changed prior to first counting, but becomes read-only after "
        "that (raising AttributeError)",
        NULL
    },
    {
        (char *)"ksize",
        (getter)hllcounter_get_ksize, (setter)hllcounter_set_ksize,
        (char *)"k-mer size for this HLL counter."
        "Can be changed prior to first counting, but becomes read-only after "
        "that (raising AttributeError)",
        NULL
    },
    {
        (char *)"counters",
        (getter)hllcounter_getcounters, NULL,
        (char *)"Read-only internal counters.",
        NULL
    },
    {NULL} /* Sentinel */
};

static PyTypeObject khmer_KHLLCounter_Type = {
    PyVarObject_HEAD_INIT(NULL, 0)
    "_khmer.KHLLCounter",                       /* tp_name */
    sizeof(khmer_KHLLCounter_Object),          /* tp_basicsize */
    0,                                         /* tp_itemsize */
    (destructor)khmer_hllcounter_dealloc,      /* tp_dealloc */
    0,                                         /* tp_print */
    0,                                         /* tp_getattr */
    0,                                         /* tp_setattr */
    0,                                         /* tp_compare */
    0,                                         /* tp_repr */
    0,                                         /* tp_as_number */
    0,                                         /* tp_as_sequence */
    0,                                         /* tp_as_mapping */
    0,                                         /* tp_hash */
    0,                                         /* tp_call */
    0,                                         /* tp_str */
    0,                                         /* tp_getattro */
    0,                                         /* tp_setattro */
    0,                                         /* tp_as_buffer */
    Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE,  /* tp_flags */
    "HyperLogLog counter",                     /* tp_doc */
    0,                                         /* tp_traverse */
    0,                                         /* tp_clear */
    0,                                         /* tp_richcompare */
    0,                                         /* tp_weaklistoffset */
    0,                                         /* tp_iter */
    0,                                         /* tp_iternext */
    khmer_hllcounter_methods,                  /* tp_methods */
    0,                                         /* tp_members */
    khmer_hllcounter_getseters,                /* tp_getset */
    0,                                         /* tp_base */
    0,                                         /* tp_dict */
    0,                                         /* tp_descr_get */
    0,                                         /* tp_descr_set */
    0,                                         /* tp_dictoffset */
    0,                                         /* tp_init */
    0,                                         /* tp_alloc */
    khmer_hllcounter_new,                      /* tp_new */
};

static PyObject * hllcounter_merge(khmer_KHLLCounter_Object * me,
                                   PyObject * args)
{
    khmer_KHLLCounter_Object * other;

    if (!PyArg_ParseTuple(args, "O!", &khmer_KHLLCounter_Type, &other)) {
        return NULL;
    }

    try {
        me->hllcounter->merge(*(other->hllcounter));
    } catch (khmer_exception &e) {
        PyErr_SetString(PyExc_ValueError, e.what());
        return NULL;
    }

    Py_RETURN_NONE;
}

/********************************
 * Assembler classes
 ********************************/


typedef struct {
    PyObject_HEAD
    LinearAssembler * assembler;
} khmer_KLinearAssembler_Object;

static void khmer_linearassembler_dealloc(khmer_KLinearAssembler_Object * obj)
{
    delete obj->assembler;
    obj->assembler = NULL;

    Py_TYPE(obj)->tp_free((PyObject*)obj);
}

static PyObject * khmer_linearassembler_new(PyTypeObject *type, PyObject *args,
        PyObject *kwds)
{
    khmer_KLinearAssembler_Object *self;
    self = (khmer_KLinearAssembler_Object*)type->tp_alloc(type, 0);

    if (self != NULL) {
        PyObject * hashgraph_o;
        Hashgraph * hashgraph = NULL;

        if (!PyArg_ParseTuple(args, "O", &hashgraph_o)) {
            Py_DECREF(self);
            return NULL;
        }

        if (PyObject_TypeCheck(hashgraph_o, &khmer_KNodegraph_Type)) {
            khmer_KNodegraph_Object * kho = (khmer_KNodegraph_Object *) hashgraph_o;
            hashgraph = kho->nodegraph;
        } else if (PyObject_TypeCheck(hashgraph_o, &khmer_KCountgraph_Type)) {
            khmer_KCountgraph_Object * cho = (khmer_KCountgraph_Object *) hashgraph_o;
            hashgraph = cho->countgraph;
        } else {
            PyErr_SetString(PyExc_ValueError,
                            "graph object must be a NodeGraph or CountGraph");
            Py_DECREF(self);
            return NULL;
        }

        try {
            self->assembler = new LinearAssembler(hashgraph);
        } catch (std::bad_alloc &e) {
            Py_DECREF(self);
            return PyErr_NoMemory();
        }

    }

    return (PyObject *) self;
}


static
PyObject *
linearassembler_assemble(khmer_KLinearAssembler_Object * me,
                         PyObject * args, PyObject *kwargs)
{
    LinearAssembler * assembler= me->assembler;

    PyObject * val_o;
    khmer_KNodegraph_Object * nodegraph_o = NULL;
    Nodegraph * stop_bf = NULL;
    const char * dir_str = NULL;
    char dir = NULL;

    const char *kwnames[] = {"seed_kmer", "stop_filter", "direction", NULL};

    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "O|O!s",
                                     const_cast<char **>(kwnames),
                                     &val_o, &khmer_KNodegraph_Type,
                                     &nodegraph_o, &dir_str)) {
        return NULL;
    }
    if (dir_str != NULL) {
        dir = dir_str[0];
    } else {
        dir = 'B';
    }

    Kmer start_kmer;
    if (!ht_convert_PyObject_to_Kmer(val_o, start_kmer, assembler->graph)) {
        return NULL;
    }

    if (nodegraph_o) {
        stop_bf = nodegraph_o->nodegraph;
    }

    std::string contig;
    if (dir == 'B') {
        contig = assembler->assemble(start_kmer, stop_bf);
    } else if (dir == 'L') {
        contig = assembler->assemble_left(start_kmer, stop_bf);
    } else if (dir == 'R') {
        contig = assembler->assemble_right(start_kmer, stop_bf);
    } else {
        PyErr_SetString(PyExc_ValueError, "Direction must be B (both), L (left),"
                        " or R (right).");
        return NULL;
    }

    PyObject * ret = Py_BuildValue("s", contig.c_str());
    return ret;
}


static PyMethodDef khmer_linearassembler_methods[] = {
    {
        "assemble",
        (PyCFunction)linearassembler_assemble, METH_VARARGS | METH_KEYWORDS,
        "Assemble a path linearly until a branch is reached."
    },
    {NULL, NULL, 0, NULL}           /* sentinel */
};

static PyTypeObject khmer_KLinearAssembler_Type = {
    PyVarObject_HEAD_INIT(NULL, 0)  /* init & ob_size */
    "_khmer.LinearAssembler",            /* tp_name */
    sizeof(khmer_KLinearAssembler_Object), /* tp_basicsize */
    0,                       /* tp_itemsize */
    (destructor)khmer_linearassembler_dealloc, /* tp_dealloc */
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
    Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE,   /* tp_flags */
    0,                       /* tp_doc */
    0,                       /* tp_traverse */
    0,                       /* tp_clear */
    0,                       /* tp_richcompare */
    0,                       /* tp_weaklistoffset */
    0,                       /* tp_iter */
    0,                       /* tp_iternext */
    khmer_linearassembler_methods, /* tp_methods */
    0,                       /* tp_members */
    0,                       /* tp_getset */
    0,                       /* tp_base */
    0,                       /* tp_dict */
    0,                       /* tp_descr_get */
    0,                       /* tp_descr_set */
    0,                       /* tp_dictoffset */
    0,                       /* tp_init */
    0,                       /* tp_alloc */
    khmer_linearassembler_new,      /* tp_new */
};



typedef struct {
    PyObject_HEAD
    SimpleLabeledAssembler * assembler;
} khmer_KSimpleLabeledAssembler_Object;


static void khmer_simplelabeledassembler_dealloc(khmer_KLinearAssembler_Object *
        obj)
{
    delete obj->assembler;
    obj->assembler = NULL;

    Py_TYPE(obj)->tp_free((PyObject*)obj);
}

static PyObject * khmer_simplelabeledassembler_new(PyTypeObject *type,
        PyObject *args,
        PyObject *kwds)
{
    khmer_KSimpleLabeledAssembler_Object *self;
    self = (khmer_KSimpleLabeledAssembler_Object*)type->tp_alloc(type, 0);

    if (self != NULL) {
        PyObject * labelhash_o;
        LabelHash * labelhash = NULL;

        if (!PyArg_ParseTuple(args, "O", &labelhash_o)) {
            Py_DECREF(self);
            return NULL;
        }

        if (PyObject_TypeCheck(labelhash_o, &khmer_KGraphLabels_Type)) {
            khmer_KGraphLabels_Object * klo = (khmer_KGraphLabels_Object *) labelhash_o;
            labelhash = klo->labelhash;
        } else {
            PyErr_SetString(PyExc_ValueError,
                            "SimpleLabeledAssembler needs a GraphLabels object.");
            Py_DECREF(self);
            return NULL;
        }

        try {
            self->assembler = new SimpleLabeledAssembler(labelhash);
        } catch (std::bad_alloc &e) {
            Py_DECREF(self);
            return PyErr_NoMemory();
        }

    }

    return (PyObject *) self;
}


static
PyObject *
simplelabeledassembler_assemble(khmer_KSimpleLabeledAssembler_Object * me,
                                PyObject * args, PyObject *kwargs)
{
    SimpleLabeledAssembler * assembler = me->assembler;

    PyObject * val_o;
    khmer_KNodegraph_Object * nodegraph_o = NULL;
    Nodegraph * stop_bf = NULL;

    const char *kwnames[] = {"seed_kmer", "stop_filter", NULL};

    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "O|O!",
                                     const_cast<char **>(kwnames),
                                     &val_o, &khmer_KNodegraph_Type,
                                     &nodegraph_o)) {
        return NULL;
    }


    Kmer start_kmer;
    if (!ht_convert_PyObject_to_Kmer(val_o, start_kmer, assembler->graph)) {
        return NULL;
    }

    if (nodegraph_o) {
        stop_bf = nodegraph_o->nodegraph;
    }

    std::vector<std::string> contigs = assembler->assemble(start_kmer, stop_bf);

    PyObject * ret = PyList_New(contigs.size());
    for (unsigned int i = 0; i < contigs.size(); i++) {
        PyList_SET_ITEM(ret, i, PyUnicode_FromString(contigs[i].c_str()));
    }

    return ret;
}


static PyMethodDef khmer_simplelabeledassembler_methods[] = {
    {
        "assemble",
        (PyCFunction)simplelabeledassembler_assemble, METH_VARARGS | METH_KEYWORDS,
        "Assemble paths, using labels to jump branches."
    },
    {NULL, NULL, 0, NULL}           /* sentinel */
};

static PyTypeObject khmer_KSimpleLabeledAssembler_Type = {
    PyVarObject_HEAD_INIT(NULL, 0)  /* init & ob_size */
    "_khmer.SimpleLabeledAssembler",            /* tp_name */
    sizeof(khmer_KSimpleLabeledAssembler_Object), /* tp_basicsize */
    0,                       /* tp_itemsize */
    (destructor)khmer_simplelabeledassembler_dealloc, /* tp_dealloc */
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
    Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE,   /* tp_flags */
    0,                       /* tp_doc */
    0,                       /* tp_traverse */
    0,                       /* tp_clear */
    0,                       /* tp_richcompare */
    0,                       /* tp_weaklistoffset */
    0,                       /* tp_iter */
    0,                       /* tp_iternext */
    khmer_simplelabeledassembler_methods, /* tp_methods */
    0,                       /* tp_members */
    0,                       /* tp_getset */
    0,                       /* tp_base */
    0,                       /* tp_dict */
    0,                       /* tp_descr_get */
    0,                       /* tp_descr_set */
    0,                       /* tp_dictoffset */
    0,                       /* tp_init */
    0,                       /* tp_alloc */
    khmer_simplelabeledassembler_new,      /* tp_new */
};



/********************************
 * JunctionCountAssembler
 ********************************/


typedef struct {
    PyObject_HEAD
    JunctionCountAssembler * assembler;
} khmer_KJunctionCountAssembler_Object;

static void khmer_junctioncountassembler_dealloc(
    khmer_KJunctionCountAssembler_Object * obj)
{
    delete obj->assembler;
    obj->assembler = NULL;

    Py_TYPE(obj)->tp_free((PyObject*)obj);
}

static PyObject * khmer_junctioncountassembler_new(PyTypeObject *type,
        PyObject *args,
        PyObject *kwds)
{
    khmer_KJunctionCountAssembler_Object *self;
    self = (khmer_KJunctionCountAssembler_Object*)type->tp_alloc(type, 0);

    if (self != NULL) {
        PyObject * hashgraph_o;
        Hashgraph * hashgraph = NULL;

        if (!PyArg_ParseTuple(args, "O", &hashgraph_o)) {
            Py_DECREF(self);
            return NULL;
        }

        if (PyObject_TypeCheck(hashgraph_o, &khmer_KNodegraph_Type)) {
            khmer_KNodegraph_Object * kho = (khmer_KNodegraph_Object *) hashgraph_o;
            hashgraph = kho->nodegraph;
        } else if (PyObject_TypeCheck(hashgraph_o, &khmer_KCountgraph_Type)) {
            khmer_KCountgraph_Object * cho = (khmer_KCountgraph_Object *) hashgraph_o;
            hashgraph = cho->countgraph;
        } else {
            PyErr_SetString(PyExc_ValueError,
                            "graph object must be a NodeGraph or CountGraph");
            Py_DECREF(self);
            return NULL;
        }

        try {
            self->assembler = new JunctionCountAssembler(hashgraph);
        } catch (std::bad_alloc &e) {
            Py_DECREF(self);
            return PyErr_NoMemory();
        }

    }

    return (PyObject *) self;
}


static
PyObject *
junctioncountassembler_assemble(khmer_KJunctionCountAssembler_Object * me,
                                PyObject * args, PyObject *kwargs)
{
    JunctionCountAssembler * assembler = me->assembler;

    PyObject * val_o;
    khmer_KNodegraph_Object * nodegraph_o = NULL;
    Nodegraph * stop_bf = NULL;

    const char *kwnames[] = {"seed_kmer", "stop_filter", NULL};

    if (!PyArg_ParseTupleAndKeywords(args, kwargs, "O|O!",
                                     const_cast<char **>(kwnames),
                                     &val_o, &khmer_KNodegraph_Type,
                                     &nodegraph_o)) {
        return NULL;
    }

    Kmer start_kmer;
    if (!ht_convert_PyObject_to_Kmer(val_o, start_kmer, assembler->graph)) {
        return NULL;
    }

    if (nodegraph_o) {
        stop_bf = nodegraph_o->nodegraph;
    }

    std::vector<std::string> contigs = assembler->assemble(start_kmer, stop_bf);

    PyObject * ret = PyList_New(contigs.size());
    for (unsigned int i = 0; i < contigs.size(); i++) {
        PyList_SET_ITEM(ret, i, PyUnicode_FromString(contigs[i].c_str()));
    }

    return ret;
}


static
PyObject *
junctioncountassembler_consume(khmer_KJunctionCountAssembler_Object * me,
                               PyObject * args)
{
    JunctionCountAssembler * assembler = me->assembler;
    const char * long_str;

    if (!PyArg_ParseTuple(args, "s", &long_str)) {
        return NULL;
    }

    if (strlen(long_str) < assembler->_ksize) {
        PyErr_SetString(PyExc_ValueError,
                        "string length must >= the hashgraph k-mer size");
        return NULL;
    }

    uint16_t n_junctions = assembler->consume(long_str);

    return PyLong_FromLong((HashIntoType) n_junctions);
}


static PyMethodDef khmer_junctioncountassembler_methods[] = {
    {
        "assemble",
        (PyCFunction)junctioncountassembler_assemble, METH_VARARGS | METH_KEYWORDS,
        "Assemble paths, using recorded junctions to jump branches."
    },
    {
        "consume",
        (PyCFunction)junctioncountassembler_consume, METH_VARARGS,
        "Consume a string and count its branch junctions."
    },
    {NULL, NULL, 0, NULL}           /* sentinel */
};

static PyTypeObject khmer_KJunctionCountAssembler_Type = {
    PyVarObject_HEAD_INIT(NULL, 0)  /* init & ob_size */
    "_khmer.JunctionCountAssembler",            /* tp_name */
    sizeof(khmer_KJunctionCountAssembler_Object), /* tp_basicsize */
    0,                       /* tp_itemsize */
    (destructor)khmer_junctioncountassembler_dealloc, /* tp_dealloc */
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
    Py_TPFLAGS_DEFAULT | Py_TPFLAGS_BASETYPE,   /* tp_flags */
    0,                       /* tp_doc */
    0,                       /* tp_traverse */
    0,                       /* tp_clear */
    0,                       /* tp_richcompare */
    0,                       /* tp_weaklistoffset */
    0,                       /* tp_iter */
    0,                       /* tp_iternext */
    khmer_junctioncountassembler_methods, /* tp_methods */
    0,                       /* tp_members */
    0,                       /* tp_getset */
    0,                       /* tp_base */
    0,                       /* tp_dict */
    0,                       /* tp_descr_get */
    0,                       /* tp_descr_set */
    0,                       /* tp_dictoffset */
    0,                       /* tp_init */
    0,                       /* tp_alloc */
    khmer_junctioncountassembler_new,      /* tp_new */
};

////

static
PyObject *
hashtable_abundance_distribution_with_reads_parser(khmer_KHashtable_Object * me,
        PyObject * args)
{
    Hashtable * hashtable = me->hashtable;

    khmer :: python :: khmer_ReadParser_Object * rparser_obj = NULL;
    khmer_KHashtable_Object * tracking_obj = NULL;

    if (!PyArg_ParseTuple(args, "O!O!", &python::khmer_ReadParser_Type,
                          &rparser_obj, &khmer_KHashtable_Type, &tracking_obj)) {
        return NULL;
    }

    FastxParserPtr& rparser = rparser_obj->parser;
    Hashtable          *tracking        = tracking_obj->hashtable;
    uint64_t           *dist            = NULL;
    const char         *value_exception = NULL;
    const char         *file_exception  = NULL;
    std::string exc_string;

    Py_BEGIN_ALLOW_THREADS
    try {
        dist = hashtable->abundance_distribution<FastxReader>(rparser, tracking);
    } catch (khmer_file_exception &exc) {
        exc_string = exc.what();
        file_exception = exc_string.c_str();
    } catch (khmer_value_exception &exc) {
        exc_string = exc.what();
        value_exception = exc_string.c_str();
    }
    Py_END_ALLOW_THREADS

    if (file_exception != NULL) {
        PyErr_SetString(PyExc_OSError, file_exception);
        return NULL;
    }
    if (value_exception != NULL) {
        PyErr_SetString(PyExc_ValueError, value_exception);
        return NULL;
    }

    PyObject * x = PyList_New(MAX_BIGCOUNT + 1);
    if (x == NULL) {
        delete[] dist;
        return NULL;
    }
    for (int i = 0; i < MAX_BIGCOUNT + 1; i++) {
        PyList_SET_ITEM(x, i, PyLong_FromUnsignedLongLong(dist[i]));
    }

    delete[] dist;
    return x;
}

static
PyObject *
hashtable_abundance_distribution(khmer_KHashtable_Object * me, PyObject * args)
{
    Hashtable * hashtable = me->hashtable;

    const char * filename = NULL;
    khmer_KHashtable_Object * tracking_obj = NULL;
    if (!PyArg_ParseTuple(args, "sO!", &filename, &khmer_KHashtable_Type,
                          &tracking_obj)) {
        return NULL;
    }

    Hashtable          *tracking       = tracking_obj->hashtable;
    uint64_t           *dist            = NULL;
    const char         *value_exception = NULL;
    const char         *file_exception  = NULL;
    std::string exc_string;

    Py_BEGIN_ALLOW_THREADS
    try {
        dist = hashtable->abundance_distribution<FastxReader>(filename, tracking);
    } catch (khmer_file_exception &exc) {
        exc_string = exc.what();
        file_exception = exc_string.c_str();
    } catch (khmer_value_exception &exc) {
        exc_string = exc.what();
        value_exception = exc_string.c_str();
    }
    Py_END_ALLOW_THREADS

    if (file_exception != NULL) {
        PyErr_SetString(PyExc_OSError, file_exception);
        if (dist != NULL) {
            delete []dist;
        }
        return NULL;
    }
    if (value_exception != NULL) {
        PyErr_SetString(PyExc_ValueError, value_exception);
        if (dist != NULL) {
            delete []dist;
        }
        return NULL;
    }

    PyObject * x = PyList_New(MAX_BIGCOUNT + 1);
    if (x == NULL) {
        if (dist != NULL) {
            delete []dist;
        }
        return NULL;
    }
    for (int i = 0; i < MAX_BIGCOUNT + 1; i++) {
        PyList_SET_ITEM(x, i, PyLong_FromUnsignedLongLong(dist[i]));
    }

    if (dist != NULL) {
        delete []dist;
    }

    return x;
}

//////////////////////////////
// standalone functions

static PyObject * forward_hash(PyObject * self, PyObject * args)
{
    const char * kmer;
    WordLength ksize;

    if (!PyArg_ParseTuple(args, "sb", &kmer, &ksize)) {
        return NULL;
    }

    if (ksize > KSIZE_MAX) {
        PyErr_Format(PyExc_ValueError, "k-mer size must be <= %u", KSIZE_MAX);
        return NULL;
    }

    if (strlen(kmer) != ksize) {
        PyErr_Format(PyExc_ValueError, "k-mer size different from ksize");
        return NULL;
    }

    try {
        PyObject * hash = nullptr;
        const HashIntoType h(_hash(kmer, ksize));
        convert_HashIntoType_to_PyObject(h, &hash);
        return hash;
    } catch (khmer_exception &e) {
        PyErr_SetString(PyExc_ValueError, e.what());
        return NULL;
    }

}

static PyObject * forward_hash_no_rc(PyObject * self, PyObject * args)
{
    const char * kmer;
    WordLength ksize;

    if (!PyArg_ParseTuple(args, "sb", &kmer, &ksize)) {
        return NULL;
    }

    if (ksize > KSIZE_MAX) {
        PyErr_Format(PyExc_ValueError, "k-mer size must be <= %u", KSIZE_MAX);
        return NULL;
    }

    if (strlen(kmer) != ksize) {
        PyErr_SetString(PyExc_ValueError,
                        "k-mer length must equal the k-size");
        return NULL;
    }

    PyObject * hash = nullptr;
    const HashIntoType h(_hash_forward(kmer, ksize));
    convert_HashIntoType_to_PyObject(h, &hash);
    return hash;
}

static PyObject * reverse_hash(PyObject * self, PyObject * args)
{
    PyObject * val;
    HashIntoType hash;
    WordLength ksize;

    if (!PyArg_ParseTuple(args, "Ob", &val, &ksize)) {
        return NULL;
    }

    if (PyLong_Check(val) || PyInt_Check(val)) {
        if (!convert_PyLong_to_HashIntoType(val, hash)) {
            return NULL;
        }
    } else {
        PyErr_SetString(PyExc_TypeError,
                        "Hash value must be an integer.");
        return NULL;
    }

    if (ksize > KSIZE_MAX) {
        PyErr_Format(PyExc_ValueError, "k-mer size must be <= %u", KSIZE_MAX);
        return NULL;
    }

    return PyUnicode_FromString(_revhash(hash, ksize).c_str());
}

static PyObject * murmur3_forward_hash(PyObject * self, PyObject * args)
{
    const char * kmer;

    if (!PyArg_ParseTuple(args, "s", &kmer)) {
        return NULL;
    }

    PyObject * hash = nullptr;
    const HashIntoType h(_hash_murmur(kmer, strlen(kmer)));
    convert_HashIntoType_to_PyObject(h, &hash);
    return hash;
}

static PyObject * murmur3_forward_hash_no_rc(PyObject * self, PyObject * args)
{
    const char * kmer;

    if (!PyArg_ParseTuple(args, "s", &kmer)) {
        return NULL;
    }

    PyObject * hash = nullptr;
    const HashIntoType h(_hash_murmur_forward(kmer, strlen(kmer)));
    convert_HashIntoType_to_PyObject(h, &hash);
    return hash;
}

static PyObject * reverse_complement(PyObject * self, PyObject * args)
{
    const char * sequence;
    if (!PyArg_ParseTuple(args, "s", &sequence)) {
        return NULL;
    }

    std::string s(sequence);
    try {
        s = _revcomp(s);
    } catch (khmer_exception &e) {
        PyErr_SetString(PyExc_RuntimeError, e.what());
        return NULL;
    }
    return PyUnicode_FromString(s.c_str());
}

//
// technique for resolving literal below found here:
// https://gcc.gnu.org/onlinedocs/gcc-4.9.1/cpp/Stringification.html
//

static
PyObject *
get_version_cpp( PyObject * self, PyObject * args )
{
#define xstr(s) str(s)
#define str(s) #s
    std::string dVersion = xstr(VERSION);
    return PyUnicode_FromString(dVersion.c_str());
}


//
// Module machinery.
//

static PyMethodDef KhmerMethods[] = {
    {
        "forward_hash",     forward_hash,
        METH_VARARGS,       "",
    },
    {
        "forward_hash_no_rc",   forward_hash_no_rc,
        METH_VARARGS,       "",
    },
    {
        "reverse_hash",     reverse_hash,
        METH_VARARGS,       "",
    },
    {
        "hash_murmur3",
        murmur3_forward_hash,
        METH_VARARGS,
        "Calculate the hash value of a k-mer using MurmurHash3 "
        "(with reverse complement)",
    },
    {
        "hash_no_rc_murmur3",
        murmur3_forward_hash_no_rc,
        METH_VARARGS,
        "Calculate the hash value of a k-mer using MurmurHash3 "
        "(no reverse complement)",
    },
    {
        "reverse_complement",
        reverse_complement,
        METH_VARARGS,
        "Calculate the reverse-complement of the DNA sequence "
        "with alphabet ACGT",
    },
    {
        "get_version_cpp", get_version_cpp,
        METH_VARARGS, "return the VERSION c++ compiler option"
    },
    { NULL, NULL, 0, NULL } // sentinel
};

MOD_INIT(_khmer)
{
    using namespace python;

    if (PyType_Ready(&khmer_KHashtable_Type) < 0) {
        return MOD_ERROR_VAL;
    }

    khmer_KCounttable_Type.tp_base = &khmer_KHashtable_Type;
    if (PyType_Ready(&khmer_KCounttable_Type) < 0) {
        return MOD_ERROR_VAL;
    }

    khmer_KSmallCounttable_Type.tp_base = &khmer_KHashtable_Type;
    if (PyType_Ready(&khmer_KSmallCounttable_Type) < 0) {
        return MOD_ERROR_VAL;
    }

    khmer_KNodetable_Type.tp_base = &khmer_KHashtable_Type;
    if (PyType_Ready(&khmer_KNodetable_Type) < 0) {
        return MOD_ERROR_VAL;
    }

    khmer_KHashgraph_Type.tp_base = &khmer_KHashtable_Type;
    khmer_KHashgraph_Type.tp_methods = khmer_hashgraph_methods;
    if (PyType_Ready(&khmer_KHashgraph_Type) < 0) {
        return MOD_ERROR_VAL;
    }

    khmer_KCountgraph_Type.tp_base = &khmer_KHashgraph_Type;
    if (PyType_Ready(&khmer_KCountgraph_Type) < 0) {
        return MOD_ERROR_VAL;
    }

    khmer_KSmallCountgraph_Type.tp_base = &khmer_KHashgraph_Type;
    if (PyType_Ready(&khmer_KSmallCountgraph_Type) < 0) {
        return MOD_ERROR_VAL;
    }

    if (PyType_Ready(&khmer_PrePartitionInfo_Type) < 0) {
        return MOD_ERROR_VAL;
    }

    khmer_KSubsetPartition_Type.tp_methods = khmer_subset_methods;
    if (PyType_Ready(&khmer_KSubsetPartition_Type) < 0) {
        return MOD_ERROR_VAL;
    }

    khmer_KNodegraph_Type.tp_base = &khmer_KHashgraph_Type;
    khmer_KNodegraph_Type.tp_methods = khmer_nodegraph_methods;
    if (PyType_Ready(&khmer_KNodegraph_Type) < 0) {
        return MOD_ERROR_VAL;
    }

    if (PyType_Ready(&khmer_KLinearAssembler_Type) < 0) {
        return MOD_ERROR_VAL;
    }
    if (PyType_Ready(&khmer_KSimpleLabeledAssembler_Type) < 0) {
        return MOD_ERROR_VAL;
    }
    if (PyType_Ready(&khmer_KJunctionCountAssembler_Type) < 0) {
        return MOD_ERROR_VAL;
    }

    khmer_KGraphLabels_Type.tp_base = &khmer_KNodegraph_Type;
    khmer_KGraphLabels_Type.tp_methods = khmer_graphlabels_methods;
    khmer_KGraphLabels_Type.tp_new = khmer_graphlabels_new;
    if (PyType_Ready(&khmer_KGraphLabels_Type) < 0) {
        return MOD_ERROR_VAL;
    }

    if (PyType_Ready(&khmer_KHLLCounter_Type) < 0) {
        return MOD_ERROR_VAL;
    }
    if (PyType_Ready(&khmer_ReadAlignerType) < 0) {
        return MOD_ERROR_VAL;
    }

    _init_ReadParser_Type_constants();
    if (PyType_Ready( &khmer_ReadParser_Type ) < 0) {
        return MOD_ERROR_VAL;
    }

    if (PyType_Ready(&khmer_Read_Type ) < 0) {
        return MOD_ERROR_VAL;
    }

    if (PyType_Ready(&khmer_ReadPairIterator_Type ) < 0) {
        return MOD_ERROR_VAL;
    }

    PyObject * m;

    MOD_DEF(m, "_khmer", "interface for the khmer module low-level extensions",
            KhmerMethods);

    if (m == NULL) {
        return MOD_ERROR_VAL;
    }

    PyObject * filetype_dict = Py_BuildValue("{s,i,s,i,s,i,s,i,s,i,s,i,s,i}",
                               "COUNTING_HT", SAVED_COUNTING_HT,
                               "HASHBITS", SAVED_HASHBITS,
                               "TAGS", SAVED_TAGS,
                               "STOPTAGS", SAVED_STOPTAGS,
                               "SUBSET", SAVED_SUBSET,
                               "LABELSET", SAVED_LABELSET,
                               "SMALLCOUNT", SAVED_SMALLCOUNT);
    if (PyModule_AddObject( m, "FILETYPES", filetype_dict ) < 0) {
        return MOD_ERROR_VAL;
    }

    Py_INCREF(&khmer_Read_Type);
    if (PyModule_AddObject( m, "Read",
                            (PyObject *)&khmer_Read_Type ) < 0) {
        return MOD_ERROR_VAL;
    }

    Py_INCREF(&khmer_ReadParser_Type);
    if (PyModule_AddObject( m, "ReadParser",
                            (PyObject *)&khmer_ReadParser_Type ) < 0) {
        return MOD_ERROR_VAL;
    }

    Py_INCREF(&khmer_KCounttable_Type);
    if (PyModule_AddObject( m, "Counttable",
                            (PyObject *)&khmer_KCounttable_Type ) < 0) {
        return MOD_ERROR_VAL;
    }

    Py_INCREF(&khmer_KSmallCounttable_Type);
    if (PyModule_AddObject( m, "SmallCounttable",
                            (PyObject *)&khmer_KSmallCounttable_Type ) < 0) {
        return MOD_ERROR_VAL;
    }

    Py_INCREF(&khmer_KNodetable_Type);
    if (PyModule_AddObject( m, "Nodetable",
                            (PyObject *)&khmer_KNodetable_Type ) < 0) {
        return MOD_ERROR_VAL;
    }

    Py_INCREF(&khmer_KCountgraph_Type);
    if (PyModule_AddObject( m, "Countgraph",
                            (PyObject *)&khmer_KCountgraph_Type ) < 0) {
        return MOD_ERROR_VAL;
    }

    Py_INCREF(&khmer_KSmallCountgraph_Type);
    if (PyModule_AddObject( m, "SmallCountgraph",
                            (PyObject *)&khmer_KSmallCountgraph_Type ) < 0) {
        return MOD_ERROR_VAL;
    }

    Py_INCREF(&khmer_KNodegraph_Type);
    if (PyModule_AddObject(m, "Nodegraph",
                           (PyObject *)&khmer_KNodegraph_Type) < 0) {
        return MOD_ERROR_VAL;
    }

    Py_INCREF(&khmer_KGraphLabels_Type);
    if (PyModule_AddObject(m, "GraphLabels",
                           (PyObject *)&khmer_KGraphLabels_Type) < 0) {
        return MOD_ERROR_VAL;
    }

    Py_INCREF(&khmer_KLinearAssembler_Type);
    if (PyModule_AddObject(m, "LinearAssembler",
                           (PyObject *)&khmer_KLinearAssembler_Type) < 0) {
        return MOD_ERROR_VAL;
    }

    Py_INCREF(&khmer_KSimpleLabeledAssembler_Type);
    if (PyModule_AddObject(m, "SimpleLabeledAssembler",
                           (PyObject *)&khmer_KSimpleLabeledAssembler_Type) < 0) {
        return MOD_ERROR_VAL;
    }

    Py_INCREF(&khmer_KJunctionCountAssembler_Type);
    if (PyModule_AddObject(m, "JunctionCountAssembler",
                           (PyObject *)&khmer_KJunctionCountAssembler_Type) < 0) {
        return MOD_ERROR_VAL;
    }

    if (PyType_Ready(&_HashSet_iter_Type) < 0) {
        return MOD_ERROR_VAL;
    }

    khmer_HashSet_Type.tp_new = khmer_HashSet_new;
    if (PyType_Ready(&khmer_HashSet_Type) < 0) {
        return MOD_ERROR_VAL;
    }

    Py_INCREF(&khmer_KHLLCounter_Type);
    if (PyModule_AddObject(m, "HLLCounter",
                           (PyObject *)&khmer_KHLLCounter_Type) < 0) {
        return MOD_ERROR_VAL;
    }

    Py_INCREF(&khmer_ReadAlignerType);
    if (PyModule_AddObject(m, "ReadAligner",
                           (PyObject *)&khmer_ReadAlignerType) < 0) {
        return MOD_ERROR_VAL;
    }

    Py_INCREF(&khmer_HashSet_Type);
    if (PyModule_AddObject(m, "HashSet",
                           (PyObject *)&khmer_HashSet_Type) < 0) {
        return MOD_ERROR_VAL;
    }

    return MOD_SUCCESS_VAL(m);
}

// vim: set ft=cpp sts=4 sw=4 tw=79:
