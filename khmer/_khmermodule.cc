//
// This file is part of khmer, http://github.com/ged-lab/khmer/, and is
// Copyright (C) Michigan State University, 2009-2015. It is licensed under
// the three-clause BSD license; see doc/LICENSE.txt.
// Contact: khmer-project@idyll.org
//

//
// A module for Python that exports khmer C++ library functions.
//

// Must be first.
#include <Python.h>

#include <iostream>

#include "khmer.hh"
#include "kmer_hash.hh"
#include "hashtable.hh"
#include "hashbits.hh"
#include "counting.hh"
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
#endif

//
// Python 2/3 compatibility: PyBytes and PyString
// https://docs.python.org/2/howto/cporting.html#str-unicode-unification
//

#include "bytesobject.h"

using namespace khmer;

//
// Function necessary for Python loading:
//

extern "C" {
    void init_khmer();
}

// Configure module logging.
//#define WITH_INTERNAL_TRACING
namespace khmer
{

namespace python
{

#ifdef WITH_INTERNAL_TRACING
#warning "Internal tracing of Python extension module is enabled."
static uint8_t const    _MODULE_TRACE_LEVEL = TraceLogger:: TLVL_DEBUG9;
static void     _trace_logger(
    uint8_t level, char const * format, ...
)
{
    static FILE *   _stream_handle  = NULL;

    if (NULL == _stream_handle) {
        _stream_handle = fopen( "pymod.log", "w" );
    }

    va_list varargs;

    if (_MODULE_TRACE_LEVEL <= level) {
        va_start( varargs, format );
        vfprintf( _stream_handle, format, varargs );
        va_end( varargs );
        fflush( _stream_handle );
    }

}
#endif


} // namespace python

} // namespace khmer


class _khmer_exception
{
private:
    std::string _message;
public:
    _khmer_exception(std::string message) : _message(message) { };
    inline const std::string get_message() const
    {
        return _message;
    };
};

class _khmer_signal : public _khmer_exception
{
public:
    _khmer_signal(std::string message) : _khmer_exception(message) { };
};

typedef pre_partition_info _pre_partition_info;

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
void
khmer_Read_dealloc(khmer_Read_Object * obj)
{
    delete obj->read;
    obj->read = NULL;
    Py_TYPE(obj)->tp_free((PyObject*)obj);
}


static
PyObject *
Read_get_name(khmer_Read_Object * obj, void * closure )
{
    return PyBytes_FromString(obj->read->name.c_str()) ;
}


static
PyObject *
Read_get_sequence(khmer_Read_Object * obj, void * closure)
{
    return PyBytes_FromString(obj->read->sequence.c_str()) ;
}


static
PyObject *
Read_get_quality(khmer_Read_Object * obj, void * closure)
{
    return PyBytes_FromString(obj->read->quality.c_str()) ;
}


static
PyObject *
Read_get_annotations(khmer_Read_Object * obj, void * closure)
{
    return PyBytes_FromString(obj->read->annotations.c_str()) ;
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
        (char *)"annotations",
        (getter)Read_get_annotations, (setter)NULL,
        (char *)"Annotations.", NULL
    },

    { NULL, NULL, NULL, NULL, NULL } // sentinel
};


static PyTypeObject khmer_Read_Type = {
    PyVarObject_HEAD_INIT(NULL, 0)        /* init & ob_size */
    "_khmer.Read",                         /* tp_name */
    sizeof(khmer_Read_Object),            /* tp_basicsize */
    0,                                    /* tp_itemsize */
    (destructor)khmer_Read_dealloc,       /* tp_dealloc */
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
};

/***********************************************************************/

//
// ReadParser object -- parse reads directly from streams
// ReadPairIterator -- return pairs of Read objects
//


typedef struct {
    PyObject_HEAD
    //! Pointer to the low-level parser object.
    read_parsers:: IParser *  parser;
} khmer_ReadParser_Object;


typedef struct {
    PyObject_HEAD
    //! Pointer to Python parser object for reference counting purposes.
    PyObject *  parent;
    //! Persistent value of pair mode across invocations.
    int pair_mode;
} khmer_ReadPairIterator_Object;


static
void
_ReadParser_dealloc(khmer_ReadParser_Object * obj)
{
    Py_DECREF(obj->parser);
    obj->parser = NULL;
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


static
PyObject *
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
        myself->parser =
            IParser:: get_parser( ifile_name );
    } catch (InvalidStreamHandle &exc) {
        PyErr_SetString( PyExc_ValueError, exc.what() );
        return NULL;
    }
    return self;
}


static
PyObject *
_ReadParser_iternext( PyObject * self )
{
    khmer_ReadParser_Object * myself  = (khmer_ReadParser_Object *)self;
    IParser *       parser  = myself->parser;

    bool    stop_iteration = false;
    char    const * exc = NULL;
    Read *  the_read_PTR;
    try {
        the_read_PTR = new Read( );
    } catch (std::bad_alloc &e) {
        return PyErr_NoMemory();
    }

    Py_BEGIN_ALLOW_THREADS
    stop_iteration = parser->is_complete( );
    if (!stop_iteration) {
        try {
            parser->imprint_next_read( *the_read_PTR );
        } catch (NoMoreReadsAvailable &e) {
            stop_iteration = true;
        } catch (StreamReadError &e) {
            exc = e.what();
        }
    }
    Py_END_ALLOW_THREADS

    // Note: Can simply return NULL instead of setting the StopIteration
    //       exception.
    if (stop_iteration) {
        delete the_read_PTR;
        return NULL;
    }

    if (exc != NULL) {
        delete the_read_PTR;
        PyErr_SetString(PyExc_IOError, exc);
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
    IParser *           parser    = parent->parser;
    uint8_t         pair_mode = myself->pair_mode;

    ReadPair    the_read_pair;
    bool    stop_iteration      = false;
    bool    unknown_pair_reading_mode   = false;
    bool    invalid_read_pair       = false;
    bool    stream_read_error = false;
    Py_BEGIN_ALLOW_THREADS
    stop_iteration = parser->is_complete( );
    if (!stop_iteration)
        try {
            parser->imprint_next_read_pair( the_read_pair, pair_mode );
        } catch (UnknownPairReadingMode &exc) {
            unknown_pair_reading_mode = true;
        } catch (InvalidReadPair &exc) {
            invalid_read_pair = true;
        } catch (StreamReadError &exc) {
            stream_read_error = true;
        } catch (NoMoreReadsAvailable &exc) {
            stop_iteration = true;
        }
    Py_END_ALLOW_THREADS

    // Note: Can return NULL instead of setting the StopIteration exception.
    if (stop_iteration) {
        return NULL;
    }

    if (unknown_pair_reading_mode) {
        PyErr_SetString(
            PyExc_ValueError, "Unknown pair reading mode supplied."
        );
        return NULL;
    }
    if (invalid_read_pair) {
        PyErr_SetString( PyExc_IOError, "Invalid read pair detected." );
        return NULL;
    }

    if (stream_read_error) {
        PyErr_SetString( PyExc_IOError, "Input file error.");
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
ReadParser_iter_read_pairs(PyObject * self, PyObject * args )
{
    int  pair_mode  = IParser:: PAIR_MODE_ERROR_ON_UNPAIRED;

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


static PyMethodDef _ReadParser_methods [ ] = {
    {
        "iter_reads",       (PyCFunction)ReadParser_iter_reads,
        METH_NOARGS,        "Iterates over reads."
    },
    {
        "iter_read_pairs",  (PyCFunction)ReadParser_iter_read_pairs,
        METH_VARARGS,       "Iterates over paired reads as pairs."
    },

    { NULL, NULL, 0, NULL } // sentinel
};


static PyTypeObject khmer_ReadParser_Type = {
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
    0,                                         /* tp_getset */
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

    PyObject * value = PyLong_FromLong( IParser:: PAIR_MODE_ALLOW_UNPAIRED );
    result = PyDict_SetItemString(cls_attrs_DICT,
                                  "PAIR_MODE_ALLOW_UNPAIRED", value);
    Py_XDECREF(value);
    if (!result) {
        Py_DECREF(cls_attrs_DICT);
        return;
    }

    value = PyLong_FromLong( IParser:: PAIR_MODE_IGNORE_UNPAIRED );
    result = PyDict_SetItemString(cls_attrs_DICT,
                                  "PAIR_MODE_IGNORE_UNPAIRED", value );
    Py_XDECREF(value);
    if (!result) {
        Py_DECREF(cls_attrs_DICT);
        return;
    }

    value = PyLong_FromLong( IParser:: PAIR_MODE_ERROR_ON_UNPAIRED );
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


static
read_parsers:: IParser *
_PyObject_to_khmer_ReadParser( PyObject * py_object )
{
    // TODO: Add type-checking.

    return ((python:: khmer_ReadParser_Object *)py_object)->parser;
}


/***********************************************************************/

//
// KCountingHash object
//

void free_pre_partition_info(void * p)
{
    _pre_partition_info * ppi = (_pre_partition_info *) p;
    delete ppi;
}

void free_subset_partition_info(void * p)
{
    SubsetPartition * subset_p = (SubsetPartition *) p;
    delete subset_p;
}

typedef struct {
    PyObject_HEAD
    CountingHash * counting;
} khmer_KCountingHash_Object;

typedef struct {
    PyObject_HEAD
    SubsetPartition * subset;
} khmer_KSubsetPartition_Object;

typedef struct {
    PyObject_HEAD
    Hashbits * hashbits;
} khmer_KHashbits_Object;

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
    PyObject_HEAD
    ReadAligner * aligner;
} khmer_ReadAligner_Object;

static void khmer_counting_dealloc(khmer_KCountingHash_Object * obj);

static
PyObject *
hash_abundance_distribution(khmer_KCountingHash_Object * me, PyObject * args);

static
PyObject *
hash_abundance_distribution_with_reads_parser(khmer_KCountingHash_Object * me,
        PyObject * args);

static
PyObject *
hash_get_raw_tables(khmer_KCountingHash_Object * self, PyObject * args)
{
    CountingHash * counting = self->counting;

    Byte ** table_ptrs = counting->get_raw_tables();
    std::vector<HashIntoType> sizes = counting->get_tablesizes();

    PyObject * raw_tables = PyList_New(sizes.size());
    for (unsigned int i=0; i<sizes.size(); ++i) {
        /*
        Py_buffer * buf = (Py_buffer * ) table_ptrs[i];
        buf->obj = NULL;
        buf->len = sizes[i];
        buf->readonly = 1;
        buf->ndim = 1;
        buf->format = NULL;
        buf->shape = NULL;
        buf->strides = NULL;
        buf->suboffsets = NULL;
        buf->internal = NULL;
        bufs.push_back(buf);
        */
        PyObject * buf = PyBuffer_FromMemory(table_ptrs[i], sizes[i]);
        if(!PyBuffer_Check(buf))
            return NULL;
        PyList_SET_ITEM(raw_tables, i, buf);
        //Py_XDECREF(buf);
    }

    return raw_tables;
}

static
PyObject *
hash_set_use_bigcount(khmer_KCountingHash_Object * me, PyObject * args)
{
    CountingHash * counting = me->counting;

    PyObject * x;
    if (!PyArg_ParseTuple(args, "O", &x)) {
        return NULL;
    }
    int setme = PyObject_IsTrue(x);
    if (setme < 0) {
        return NULL;
    }
    counting->set_use_bigcount((bool)setme);

    Py_RETURN_NONE;
}

static
PyObject *
hash_get_use_bigcount(khmer_KCountingHash_Object * me, PyObject * args)
{
    CountingHash * counting = me->counting;

    if (!PyArg_ParseTuple(args, "")) {
        return NULL;
    }

    bool val = counting->get_use_bigcount();

    return PyBool_FromLong((int)val);
}

static
PyObject *
hash_n_occupied(khmer_KCountingHash_Object * me, PyObject * args)
{
    CountingHash * counting = me->counting;

    HashIntoType start = 0, stop = 0;

    if (!PyArg_ParseTuple(args, "|KK", &start, &stop)) {
        return NULL;
    }

    HashIntoType n = counting->n_occupied(start, stop);

    return PyLong_FromUnsignedLongLong(n);
}

static
PyObject *
hash_n_unique_kmers(khmer_KCountingHash_Object * me, PyObject * args)
{
    CountingHash * counting = me->counting;

    HashIntoType n = counting->n_unique_kmers();

    return PyLong_FromUnsignedLongLong(n);
}

static
PyObject *
hash_n_entries(khmer_KCountingHash_Object * me, PyObject * args)
{
    CountingHash * counting = me->counting;

    if (!PyArg_ParseTuple(args, "")) {
        return NULL;
    }

    return PyLong_FromUnsignedLongLong(counting->n_entries());
}

static
PyObject *
hash_count(khmer_KCountingHash_Object * me, PyObject * args)
{
    CountingHash * counting = me->counting;

    const char * kmer;

    if (!PyArg_ParseTuple(args, "s", &kmer)) {
        return NULL;
    }

    if (strlen(kmer) != counting->ksize()) {
        PyErr_SetString(PyExc_ValueError,
                        "k-mer length must be the same as the hashtable k-size");
        return NULL;
    }

    counting->count(kmer);

    return PyLong_FromLong(1);
}

static
PyObject *
hash_output_fasta_kmer_pos_freq(khmer_KCountingHash_Object * me,
                                PyObject * args)
{
    CountingHash * counting = me->counting;

    const char * infile;
    const char * outfile;

    if (!PyArg_ParseTuple(args, "ss", &infile, &outfile)) {
        return NULL;
    }

    counting->output_fasta_kmer_pos_freq(infile, outfile);

    return PyLong_FromLong(0);
}

static
PyObject *
hash_consume_fasta(khmer_KCountingHash_Object * me, PyObject * args)
{
    CountingHash * counting  = me->counting;

    const char * filename;

    if (!PyArg_ParseTuple(args, "s", &filename)) {
        return NULL;
    }

    // call the C++ function, and trap signals => Python
    unsigned long long  n_consumed    = 0;
    unsigned int          total_reads   = 0;
    try {
        counting->consume_fasta(filename, total_reads, n_consumed);
    } catch (_khmer_signal &e) {
        PyErr_SetString(PyExc_IOError, e.get_message().c_str());
        return NULL;
    } catch (khmer_file_exception &e) {
        PyErr_SetString(PyExc_IOError, e.what());
        return NULL;
    }

    return Py_BuildValue("IK", total_reads, n_consumed);
}

static
PyObject *
hash_consume_fasta_with_reads_parser(khmer_KCountingHash_Object * me,
                                     PyObject * args)
{
    CountingHash * counting  = me->counting;

    PyObject * rparser_obj = NULL;

    if (!PyArg_ParseTuple(args, "O", &rparser_obj)) {
        return NULL;
    }

    read_parsers:: IParser * rparser =
        _PyObject_to_khmer_ReadParser( rparser_obj );

    char const * exc = "";
    // call the C++ function, and trap signals => Python
    unsigned long long  n_consumed  = 0;
    unsigned int    total_reads = 0;
    bool        exc_raised  = false;
    Py_BEGIN_ALLOW_THREADS
    try {
        counting->consume_fasta(rparser, total_reads, n_consumed);
    } catch (_khmer_signal &e) {
        exc = e.get_message().c_str();
        exc_raised = true;
    }
    Py_END_ALLOW_THREADS
    if (exc_raised) {
        PyErr_SetString(PyExc_IOError, exc);
        return NULL;
    }

    return Py_BuildValue("IK", total_reads, n_consumed);
}

static
PyObject *
hash_consume(khmer_KCountingHash_Object * me, PyObject * args)
{
    CountingHash * counting = me->counting;

    const char * long_str;

    if (!PyArg_ParseTuple(args, "s", &long_str)) {
        return NULL;
    }

    if (strlen(long_str) < counting->ksize()) {
        PyErr_SetString(PyExc_ValueError,
                        "string length must >= the hashtable k-mer size");
        return NULL;
    }

    unsigned int n_consumed;
    n_consumed = counting->consume_string(long_str);

    return PyLong_FromLong(n_consumed);
}

static
PyObject *
hash_get_min_count(khmer_KCountingHash_Object * me, PyObject * args)
{
    CountingHash * counting = me->counting;

    const char * long_str;

    if (!PyArg_ParseTuple(args, "s", &long_str)) {
        return NULL;
    }

    if (strlen(long_str) < counting->ksize()) {
        PyErr_SetString(PyExc_ValueError,
                        "string length must >= the hashtable k-mer size");
        return NULL;
    }

    BoundedCounterType c = counting->get_min_count(long_str);
    unsigned int N = c;

    return PyLong_FromLong(N);
}

static
PyObject *
hash_get_max_count(khmer_KCountingHash_Object * me, PyObject * args)
{
    CountingHash * counting = me->counting;

    const char * long_str;

    if (!PyArg_ParseTuple(args, "s", &long_str)) {
        return NULL;
    }

    if (strlen(long_str) < counting->ksize()) {
        PyErr_SetString(PyExc_ValueError,
                        "string length must >= the hashtable k-mer size");
        return NULL;
    }

    BoundedCounterType c = counting->get_max_count(long_str);
    unsigned int N = c;

    return PyLong_FromLong(N);
}

static
PyObject *
hash_get_median_count(khmer_KCountingHash_Object * me, PyObject * args)
{
    CountingHash * counting = me->counting;

    const char * long_str;

    if (!PyArg_ParseTuple(args, "s", &long_str)) {
        return NULL;
    }

    if (strlen(long_str) < counting->ksize()) {
        PyErr_SetString(PyExc_ValueError,
                        "string length must >= the hashtable k-mer size");
        return NULL;
    }

    BoundedCounterType med = 0;
    float average = 0, stddev = 0;

    counting->get_median_count(long_str, med, average, stddev);

    return Py_BuildValue("iff", med, average, stddev);
}

static
PyObject *
hash_get_kadian_count(khmer_KCountingHash_Object * me, PyObject * args)
{
    CountingHash * counting = me->counting;

    const char * long_str;
    unsigned int nk = 1;

    if (!PyArg_ParseTuple(args, "s|I", &long_str, &nk)) {
        return NULL;
    }

    if (strlen(long_str) < counting->ksize()) {
        PyErr_SetString(PyExc_ValueError,
                        "string length must >= the hashtable k-mer size");
        return NULL;
    }

    BoundedCounterType kad = 0;

    counting->get_kadian_count(long_str, kad, nk);

    return Py_BuildValue("i", kad);
}

static
PyObject *
hash_get(khmer_KCountingHash_Object * me, PyObject * args)
{
    CountingHash * counting = me->counting;

    PyObject * arg;

    if (!PyArg_ParseTuple(args, "O", &arg)) {
        return NULL;
    }

    unsigned long count = 0;

    if (PyInt_Check(arg)) {
        long pos = PyInt_AsLong(arg);
        count = counting->get_count((unsigned int) pos);
    } else if (PyBytes_Check(arg)) {
        std::string s = PyBytes_AsString(arg);

        if (strlen(s.c_str()) != counting->ksize()) {
            PyErr_SetString(PyExc_ValueError,
                            "k-mer size must equal the counting table k-mer size");
            return NULL;
        }

        count = counting->get_count(s.c_str());
    }

    return PyLong_FromLong(count);
}

static
PyObject *
count_trim_on_abundance(khmer_KCountingHash_Object * me, PyObject * args)
{
    CountingHash * counting = me->counting;

    const char * seq = NULL;
    unsigned int min_count_i = 0;

    if (!PyArg_ParseTuple(args, "sI", &seq, &min_count_i)) {
        return NULL;
    }

    unsigned long trim_at;
    Py_BEGIN_ALLOW_THREADS

    BoundedCounterType min_count = min_count_i;

    trim_at = counting->trim_on_abundance(seq, min_count);

    Py_END_ALLOW_THREADS;

    PyObject * trim_seq = PyBytes_FromStringAndSize(seq, trim_at);
    if (trim_seq == NULL) {
        return NULL;
    }
    PyObject * ret = Py_BuildValue("Ok", trim_seq, trim_at);
    Py_DECREF(trim_seq);

    return ret;
}

static
PyObject *
count_trim_below_abundance(khmer_KCountingHash_Object * me, PyObject * args)
{
    CountingHash * counting = me->counting;

    const char * seq = NULL;
    BoundedCounterType max_count_i = 0;

    if (!PyArg_ParseTuple(args, "sH", &seq, &max_count_i)) {
        return NULL;
    }

    unsigned long trim_at;
    Py_BEGIN_ALLOW_THREADS

    BoundedCounterType max_count = max_count_i;

    trim_at = counting->trim_below_abundance(seq, max_count);

    Py_END_ALLOW_THREADS;

    PyObject * trim_seq = PyBytes_FromStringAndSize(seq, trim_at);
    if (trim_seq == NULL) {
        return NULL;
    }
    PyObject * ret = Py_BuildValue("Ok", trim_seq, trim_at);
    Py_DECREF(trim_seq);

    return ret;
}

static
PyObject *
count_find_spectral_error_positions(khmer_KCountingHash_Object * me,
                                    PyObject * args)
{
    khmer::CountingHash * counting = me->counting;

    char * seq = NULL;
    khmer::BoundedCounterType max_count = 0; // unsigned short int

    if (!PyArg_ParseTuple(args, "sH", &seq, &max_count)) {
        return NULL;
    }

    std::vector<unsigned int> posns;

    try {
        posns = counting->find_spectral_error_positions(seq, max_count);
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
hash_fasta_count_kmers_by_position(khmer_KCountingHash_Object * me,
                                   PyObject * args)
{
    CountingHash * counting = me->counting;

    const char * inputfile;
    unsigned int max_read_len = 0;
    long max_read_len_long;
    int limit_by_count_int;

    if (!PyArg_ParseTuple(args, "sli", &inputfile, &max_read_len_long,
                          &limit_by_count_int)) {
        return NULL;
    }
    if (max_read_len_long < 0 || max_read_len_long >= pow(2, 32)) {
        PyErr_SetString(
            PyExc_ValueError,
            "The 2nd argument must be positive and less than 2^32");
        return NULL;
    }
    if (limit_by_count_int < 0 || limit_by_count_int >= pow(2, 16)) {
        PyErr_SetString(
            PyExc_ValueError,
            "The 3rd argument must be positive and less than 2^16");
        return NULL;
    }
    max_read_len = (unsigned int) max_read_len_long;

    unsigned long long * counts;
    counts = counting->fasta_count_kmers_by_position(inputfile, max_read_len,
             (unsigned short) limit_by_count_int);

    PyObject * x = PyList_New(max_read_len);
    if (x == NULL) {
        delete[] counts;
        return NULL;
    }

    for (unsigned int i = 0; i < max_read_len; i++) {
        int ret = PyList_SetItem(x, i, PyLong_FromUnsignedLongLong(counts[i]));
        if (ret < 0) {
            delete[] counts;
            return NULL;
        }
    }

    delete[] counts;

    return x;
}

static
PyObject *
hash_fasta_dump_kmers_by_abundance(khmer_KCountingHash_Object * me,
                                   PyObject * args)
{
    CountingHash * counting = me->counting;

    const char * inputfile;
    int limit_by = 0;

    if (!PyArg_ParseTuple(args, "si", &inputfile, &limit_by)) {
        return NULL;
    }

    counting->fasta_dump_kmers_by_abundance(inputfile,
                                            limit_by);

    Py_RETURN_NONE;
}

static
PyObject *
hash_load(khmer_KCountingHash_Object * me, PyObject * args)
{
    CountingHash * counting = me->counting;

    const char * filename = NULL;

    if (!PyArg_ParseTuple(args, "s", &filename)) {
        return NULL;
    }

    try {
        counting->load(filename);
    } catch (khmer_file_exception &e) {
        PyErr_SetString(PyExc_IOError, e.what());
        return NULL;
    }

    Py_RETURN_NONE;
}

static
PyObject *
hash_save(khmer_KCountingHash_Object * me, PyObject * args)
{
    CountingHash * counting = me->counting;

    const char * filename = NULL;

    if (!PyArg_ParseTuple(args, "s", &filename)) {
        return NULL;
    }

    counting->save(filename);

    Py_RETURN_NONE;
}

static
PyObject *
hash_get_ksize(khmer_KCountingHash_Object * me, PyObject * args)
{
    CountingHash * counting = me->counting;

    if (!PyArg_ParseTuple(args, "")) {
        return NULL;
    }

    unsigned int k = counting->ksize();

    return PyLong_FromLong(k);
}

static
PyObject *
hash_get_hashsizes(khmer_KCountingHash_Object * me, PyObject * args)
{
    CountingHash * counting = me->counting;


    if (!PyArg_ParseTuple(args, "")) {
        return NULL;
    }

    std::vector<HashIntoType> ts = counting->get_tablesizes();

    PyObject * x = PyList_New(ts.size());
    for (size_t i = 0; i < ts.size(); i++) {
        PyList_SET_ITEM(x, i, PyLong_FromUnsignedLongLong(ts[i]));
    }

    return x;
}

static
PyObject *
hash_collect_high_abundance_kmers(khmer_KCountingHash_Object * me,
                                  PyObject * args);

static
PyObject *
hash_consume_and_tag(khmer_KCountingHash_Object * me, PyObject * args)
{
    CountingHash * counting = me->counting;

    const char * seq;

    if (!PyArg_ParseTuple(args, "s", &seq)) {
        return NULL;
    }

    // call the C++ function, and trap signals => Python

    unsigned long long n_consumed = 0;
    try {
        // @CTB needs to normalize
        counting->consume_sequence_and_tag(seq, n_consumed);
    } catch (_khmer_signal &e) {
        PyErr_SetString(PyExc_ValueError, e.get_message().c_str());
        return NULL;
    }

    return Py_BuildValue("K", n_consumed);
}

static
PyObject *
hash_get_tags_and_positions(khmer_KCountingHash_Object * me, PyObject * args)
{
    CountingHash * counting = me->counting;

    const char * seq;

    if (!PyArg_ParseTuple(args, "s", &seq)) {
        return NULL;
    }

    // call the C++ function, and trap signals => Python

    std::vector<unsigned int> posns;
    std::vector<HashIntoType> tags;

    unsigned int pos = 1;
    KMerIterator kmers(seq, counting->ksize());

    while (!kmers.done()) {
        HashIntoType kmer = kmers.next();
        if (set_contains(counting->all_tags, kmer)) {
            posns.push_back(pos);
            tags.push_back(kmer);
        }
        pos++;
    }

    PyObject * posns_list = PyList_New(posns.size());
    for (size_t i = 0; i < posns.size(); i++) {
        PyObject * tup = Py_BuildValue("IK", posns[i], tags[i]);
        PyList_SET_ITEM(posns_list, i, tup);
    }

    return posns_list;
}

static
PyObject *
hash_find_all_tags_list(khmer_KCountingHash_Object * me, PyObject * args)
{
    CountingHash * counting = me->counting;

    const char * kmer_s = NULL;

    if (!PyArg_ParseTuple(args, "s", &kmer_s)) {
        return NULL;
    }

    if (strlen(kmer_s) != counting->ksize()) {
        PyErr_SetString(PyExc_ValueError,
                        "k-mer length must equal the counting table k-mer size");
        return NULL;
    }

    SeenSet tags;

    Py_BEGIN_ALLOW_THREADS

    HashIntoType kmer_f, kmer_r;
    _hash(kmer_s, counting->ksize(), kmer_f, kmer_r);

    counting->partition->find_all_tags(kmer_f, kmer_r, tags,
                                       counting->all_tags);

    Py_END_ALLOW_THREADS

    PyObject * x =  PyList_New(tags.size());
    if (x == NULL) {
        return NULL;
    }
    SeenSet::iterator si;
    unsigned long long i = 0;
    for (si = tags.begin(); si != tags.end(); ++si) {
        // type K for python unsigned long long
        PyList_SET_ITEM(x, i, Py_BuildValue("K", *si));
        i++;
    }

    return x;
}

static
PyObject *
hash_consume_fasta_and_tag(khmer_KCountingHash_Object * me, PyObject * args)
{
    CountingHash * counting = me->counting;

    const char * filename;

    if (!PyArg_ParseTuple(args, "s", &filename)) {
        return NULL;
    }

    // call the C++ function, and trap signals => Python

    unsigned long long n_consumed;
    unsigned int total_reads;

    try {
        counting->consume_fasta_and_tag(filename, total_reads, n_consumed);
    } catch (_khmer_signal &e) {
        PyErr_SetString(PyExc_IOError, e.get_message().c_str());
        return NULL;
    }

    return Py_BuildValue("IK", total_reads, n_consumed);
}

static
PyObject *
hash_find_all_tags_truncate_on_abundance(khmer_KCountingHash_Object * me,
        PyObject * args)
{
    CountingHash * counting = me->counting;

    const char * kmer_s = NULL;
    BoundedCounterType min_count, max_count;

    if (!PyArg_ParseTuple(args, "sHH", &kmer_s, &min_count, &max_count)) {
        return NULL;
    }

    if (strlen(kmer_s) != counting->ksize()) {
        PyErr_SetString(PyExc_ValueError,
                        "k-mer size must equal the k-mer size of the counting table");
        return NULL;
    }

    _pre_partition_info * ppi = NULL;

    Py_BEGIN_ALLOW_THREADS

    HashIntoType kmer, kmer_f, kmer_r;
    kmer = _hash(kmer_s, counting->ksize(), kmer_f, kmer_r);

    try {
        ppi = new _pre_partition_info(kmer);
    } catch (std::bad_alloc &e) {
        return PyErr_NoMemory();
    }
    counting->partition->find_all_tags_truncate_on_abundance(kmer_f, kmer_r,
            ppi->tagged_kmers,
            counting->all_tags,
            min_count,
            max_count);
    counting->add_kmer_to_tags(kmer);

    Py_END_ALLOW_THREADS

    return PyCObject_FromVoidPtr(ppi, free_pre_partition_info);
}

static
PyObject *
hash_do_subset_partition_with_abundance(khmer_KCountingHash_Object * me,
                                        PyObject * args)
{
    CountingHash * counting = me->counting;

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
        subset_p = new SubsetPartition(counting);
        subset_p->do_partition_with_abundance(start_kmer, end_kmer,
                                              min_count, max_count,
                                              break_on_stop_tags,
                                              stop_big_traversals);
        Py_END_ALLOW_THREADS
    } catch (_khmer_signal &e) {
        return NULL;
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

static PyMethodDef khmer_counting_methods[] = {
    {
        "ksize",
        (PyCFunction)hash_get_ksize,
        METH_VARARGS,
        ""
    },
    { "hashsizes", (PyCFunction)hash_get_hashsizes, METH_VARARGS, "" },
    { "set_use_bigcount", (PyCFunction)hash_set_use_bigcount, METH_VARARGS, "" },
    { "get_use_bigcount", (PyCFunction)hash_get_use_bigcount, METH_VARARGS, "" },
    { "n_unique_kmers", (PyCFunction)hash_n_unique_kmers, METH_VARARGS, "Count the number of unique kmers" },
    { "n_occupied", (PyCFunction)hash_n_occupied, METH_VARARGS, "Count the number of occupied bins" },
    { "n_entries", (PyCFunction)hash_n_entries, METH_VARARGS, "" },
    { "count", (PyCFunction)hash_count, METH_VARARGS, "Count the given kmer" },
    { "consume", (PyCFunction)hash_consume, METH_VARARGS, "Count all k-mers in the given string" },
    { "consume_fasta", (PyCFunction)hash_consume_fasta, METH_VARARGS, "Count all k-mers in a given file" },
    {
        "consume_fasta_with_reads_parser", (PyCFunction)hash_consume_fasta_with_reads_parser,
        METH_VARARGS, "Count all k-mers using a given reads parser"
    },
    { "output_fasta_kmer_pos_freq", (PyCFunction)hash_output_fasta_kmer_pos_freq, METH_VARARGS, "" },
    { "get", (PyCFunction)hash_get, METH_VARARGS, "Get the count for the given k-mer" },
    { "get_raw_tables", (PyCFunction)hash_get_raw_tables,
       METH_VARARGS, "Get a list of the raw tables as memoryview objects"
    },
    { "get_min_count", (PyCFunction)hash_get_min_count, METH_VARARGS, "Get the smallest count of all the k-mers in the string" },
    { "get_max_count", (PyCFunction)hash_get_max_count, METH_VARARGS, "Get the largest count of all the k-mers in the string" },
    { "get_median_count", (PyCFunction)hash_get_median_count, METH_VARARGS, "Get the median, average, and stddev of the k-mer counts in the string" },
    { "get_kadian_count", (PyCFunction)hash_get_kadian_count, METH_VARARGS, "Get the kadian (abundance of k-th rank-ordered k-mer) of the k-mer counts in the string" },
    { "trim_on_abundance", (PyCFunction)count_trim_on_abundance, METH_VARARGS, "Trim on >= abundance" },
    { "trim_below_abundance", (PyCFunction)count_trim_below_abundance, METH_VARARGS, "Trim on >= abundance" },
    { "find_spectral_error_positions", (PyCFunction)count_find_spectral_error_positions, METH_VARARGS, "Identify positions of low-abundance k-mers" },
    { "abundance_distribution", (PyCFunction)hash_abundance_distribution, METH_VARARGS, "" },
    { "abundance_distribution_with_reads_parser", (PyCFunction)hash_abundance_distribution_with_reads_parser, METH_VARARGS, "" },
    { "fasta_count_kmers_by_position", (PyCFunction)hash_fasta_count_kmers_by_position, METH_VARARGS, "" },
    { "fasta_dump_kmers_by_abundance", (PyCFunction)hash_fasta_dump_kmers_by_abundance, METH_VARARGS, "" },
    { "load", (PyCFunction)hash_load, METH_VARARGS, "" },
    { "save", (PyCFunction)hash_save, METH_VARARGS, "" },
    {
        "collect_high_abundance_kmers", (PyCFunction)hash_collect_high_abundance_kmers,
        METH_VARARGS, ""
    },
    { "consume_and_tag", (PyCFunction)hash_consume_and_tag, METH_VARARGS, "Consume a sequence and tag it" },
    { "get_tags_and_positions", (PyCFunction)hash_get_tags_and_positions, METH_VARARGS, "Retrieve tags and their positions in a sequence." },
    { "find_all_tags_list", (PyCFunction)hash_find_all_tags_list, METH_VARARGS, "Find all tags within range of the given k-mer, return as list" },
    { "consume_fasta_and_tag", (PyCFunction)hash_consume_fasta_and_tag, METH_VARARGS, "Count all k-mers in a given file" },
    { "do_subset_partition_with_abundance", (PyCFunction)hash_do_subset_partition_with_abundance, METH_VARARGS, "" },
    { "find_all_tags_truncate_on_abundance", (PyCFunction)hash_find_all_tags_truncate_on_abundance, METH_VARARGS, "" },

    {NULL, NULL, 0, NULL}           /* sentinel */
};

static PyObject* _new_counting_hash(PyTypeObject * type, PyObject * args,
                                    PyObject * kwds);

static PyTypeObject khmer_KCountingHash_Type
CPYCHECKER_TYPE_OBJECT_FOR_TYPEDEF("khmer_KCountingHash_Object")
= {
    PyVarObject_HEAD_INIT(NULL, 0)       /* init & ob_size */
    "_khmer.KCountingHash",              /*tp_name*/
    sizeof(khmer_KCountingHash_Object),  /*tp_basicsize*/
    0,                                   /*tp_itemsize*/
    (destructor)khmer_counting_dealloc,  /*tp_dealloc*/
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
    "counting hash object",              /* tp_doc */
    0,                                   /* tp_traverse */
    0,                                   /* tp_clear */
    0,                                   /* tp_richcompare */
    0,                                   /* tp_weaklistoffset */
    0,                                   /* tp_iter */
    0,                                   /* tp_iternext */
    khmer_counting_methods,              /* tp_methods */
    0,                                   /* tp_members */
    0,                                   /* tp_getset */
    0,                                   /* tp_base */
    0,                                   /* tp_dict */
    0,                                   /* tp_descr_get */
    0,                                   /* tp_descr_set */
    0,                                   /* tp_dictoffset */
    0,                                   /* tp_init */
    0,                                   /* tp_alloc */
    _new_counting_hash,                  /* tp_new */
};

#define is_counting_obj(v)  (Py_TYPE(v) == &khmer_KCountingHash_Type)

//
// new_hashtable
//

static PyObject* new_hashtable(PyObject * self, PyObject * args)
{
    unsigned int k = 0;
    unsigned long long size = 0;

    if (!PyArg_ParseTuple(args, "IK", &k, &size)) {
        return NULL;
    }

    khmer_KCountingHash_Object * kcounting_obj = (khmer_KCountingHash_Object *) \
            PyObject_New(khmer_KCountingHash_Object, &khmer_KCountingHash_Type);

    if (kcounting_obj == NULL) {
        return NULL;
    }

    try {
        kcounting_obj->counting = new CountingHash(k, size);
    } catch (std::bad_alloc &e) {
        return PyErr_NoMemory();
    }

    return (PyObject *) kcounting_obj;
}

//
// new_counting_hash
//

static PyObject* _new_counting_hash(PyTypeObject * type, PyObject * args,
                                    PyObject * kwds)
{
    khmer_KCountingHash_Object * self;

    self = (khmer_KCountingHash_Object *)type->tp_alloc(type, 0);

    if (self != NULL) {
        WordLength k = 0;
        PyListObject * sizes_list_o = NULL;

        if (!PyArg_ParseTuple(args, "bO!", &k, &PyList_Type, &sizes_list_o)) {
            Py_DECREF(self);
            return NULL;
        }

        std::vector<HashIntoType> sizes;
        Py_ssize_t sizes_list_o_length = PyList_GET_SIZE(sizes_list_o);
        if (sizes_list_o_length == -1) {
            Py_DECREF(self);
            PyErr_SetString(PyExc_ValueError, "error with hashtable primes!");
            return NULL;
        }
        for (Py_ssize_t i = 0; i < sizes_list_o_length; i++) {
            PyObject * size_o = PyList_GET_ITEM(sizes_list_o, i);
            if (PyLong_Check(size_o)) {
                sizes.push_back((HashIntoType) PyLong_AsUnsignedLongLong(size_o));
            } else if (PyInt_Check(size_o)) {
                sizes.push_back((HashIntoType) PyInt_AsLong(size_o));
            } else if (PyFloat_Check(size_o)) {
                sizes.push_back((HashIntoType) PyFloat_AS_DOUBLE(size_o));
            } else {
                Py_DECREF(self);
                PyErr_SetString(PyExc_TypeError,
                                "2nd argument must be a list of ints, longs, or floats");
                return NULL;
            }
        }

        try {
            self->counting = new CountingHash(k, sizes);
        } catch (std::bad_alloc &e) {
            return PyErr_NoMemory();
        }
    }

    return (PyObject *) self;
}

//
// hashbits stuff
//

static void khmer_hashbits_dealloc(khmer_KHashbits_Object * obj);
static PyObject* khmer_hashbits_new(PyTypeObject * type, PyObject * args,
                                    PyObject * kwds);
static int khmer_hashbits_init(khmer_KHashbits_Object * self, PyObject * args,
                               PyObject * kwds);

static PyTypeObject khmer_KHashbits_Type
CPYCHECKER_TYPE_OBJECT_FOR_TYPEDEF("khmer_KHashbits_Object")
= {
    PyVarObject_HEAD_INIT(NULL, 0) /* init & ob_size */
    "_khmer.Hashbits",             /* tp_name */
    sizeof(khmer_KHashbits_Object), /* tp_basicsize */
    0,                             /* tp_itemsize */
    (destructor)khmer_hashbits_dealloc, /*tp_dealloc*/
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
    "hashbits object",           /* tp_doc */
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
    (initproc)khmer_hashbits_init,   /* tp_init */
    0,                       /* tp_alloc */
    khmer_hashbits_new,                  /* tp_new */
};

static
PyObject *
hash_abundance_distribution_with_reads_parser(khmer_KCountingHash_Object * me,
        PyObject * args)
{
    CountingHash * counting = me->counting;

    khmer :: python :: khmer_ReadParser_Object * rparser_obj = NULL;
    khmer_KHashbits_Object *tracking_obj = NULL;

    if (!PyArg_ParseTuple(args, "O!O!", &python::khmer_ReadParser_Type,
                          &rparser_obj, &khmer_KHashbits_Type, &tracking_obj)) {
        return NULL;
    }

    read_parsers:: IParser * rparser = rparser_obj->parser;
    Hashbits * hashbits = tracking_obj->hashbits;

    HashIntoType * dist = NULL;

    const char * exception = NULL;
    Py_BEGIN_ALLOW_THREADS
    try {
        dist = counting->abundance_distribution(rparser, hashbits);
    } catch (khmer::read_parsers::NoMoreReadsAvailable &exc ) {
        exception = exc.what();
    }
    Py_END_ALLOW_THREADS
    if (exception != NULL) {
        delete[] dist;
        PyErr_SetString(PyExc_IOError, exception);
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
hash_abundance_distribution(khmer_KCountingHash_Object * me, PyObject * args)
{
    CountingHash * counting = me->counting;

    const char * filename = NULL;
    khmer_KHashbits_Object * tracking_obj = NULL;
    if (!PyArg_ParseTuple(args, "sO!", &filename, &khmer_KHashbits_Type,
                          &tracking_obj)) {
        return NULL;
    }

    Hashbits * hashbits = tracking_obj->hashbits;
    HashIntoType * dist;

    char const * result = "";
    bool exception = false;
    Py_BEGIN_ALLOW_THREADS
    try {
        dist = counting->abundance_distribution(filename, hashbits);
    } catch (khmer_file_exception &e) {
        exception = true;
        result = e.what();
    }
    Py_END_ALLOW_THREADS

    if (exception) {
        PyErr_SetString(PyExc_IOError, result);
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
hashbits_n_unique_kmers(khmer_KHashbits_Object * me, PyObject * args)
{
    Hashbits * hashbits = me->hashbits;

    HashIntoType n = hashbits->n_unique_kmers();

    return PyLong_FromUnsignedLongLong(n);
}


static
PyObject *
hashbits_count_overlap(khmer_KHashbits_Object * me, PyObject * args)
{
    Hashbits * hashbits = me->hashbits;
    khmer_KHashbits_Object * ht2_argu;
    const char * filename;
    Hashbits * ht2;

    if (!PyArg_ParseTuple(args, "sO!", &filename, &khmer_KHashbits_Type,
                          &ht2_argu)) {
        return NULL;
    }

    ht2 = ht2_argu->hashbits;

// call the C++ function, and trap signals => Python

    unsigned long long n_consumed;
    unsigned int total_reads;
    HashIntoType curve[2][100];

    try {
        hashbits->consume_fasta_overlap(filename, curve, *ht2, total_reads, n_consumed);
    } catch (_khmer_signal &e) {
        PyErr_SetString(PyExc_IOError, e.get_message().c_str());
        return NULL;
    }

    HashIntoType n = hashbits->n_unique_kmers();
    HashIntoType n_overlap = hashbits->n_overlap_kmers();

    PyObject * x = PyList_New(200);

    for (unsigned int i = 0; i < 100; i++) {
        PyList_SetItem(x, i, Py_BuildValue("K", curve[0][i]));
    }
    for (unsigned int i = 0; i < 100; i++) {
        PyList_SetItem(x, i + 100, Py_BuildValue("K", curve[1][i]));
    }
    return Py_BuildValue("KKO", n, n_overlap, x);
}

static
PyObject *
hashbits_n_occupied(khmer_KHashbits_Object * me, PyObject * args)
{
    Hashbits * hashbits = me->hashbits;

    HashIntoType start = 0, stop = 0;

    if (!PyArg_ParseTuple(args, "|KK", &start, &stop)) {
        return NULL;
    }

    HashIntoType n = hashbits->n_occupied(start, stop);

    return PyLong_FromUnsignedLongLong(n);
}

static
PyObject *
hashbits_n_tags(khmer_KHashbits_Object * me, PyObject * args)
{
    Hashbits * hashbits = me->hashbits;

    if (!PyArg_ParseTuple(args, "")) {
        return NULL;
    }

    return PyLong_FromSize_t(hashbits->n_tags());
}

static
PyObject *
hashbits_count(khmer_KHashbits_Object * me, PyObject * args)
{
    Hashbits * hashbits = me->hashbits;

    const char * kmer;

    if (!PyArg_ParseTuple(args, "s", &kmer)) {
        return NULL;
    }

    if (strlen(kmer) != hashbits->ksize()) {
        PyErr_SetString(PyExc_ValueError,
                        "k-mer length must equal the presence table k-mer size");
        return NULL;
    }

    hashbits->count(kmer);

    return PyLong_FromLong(1);
}

static
PyObject *
hashbits_consume(khmer_KHashbits_Object * me, PyObject * args)
{
    Hashbits * hashbits = me->hashbits;

    const char * long_str;

    if (!PyArg_ParseTuple(args, "s", &long_str)) {
        return NULL;
    }

    if (strlen(long_str) < hashbits->ksize()) {
        PyErr_SetString(PyExc_ValueError,
                        "string length must >= the hashbits k-mer size");
        return NULL;
    }

    unsigned int n_consumed;
    n_consumed = hashbits->consume_string(long_str);

    return PyLong_FromLong(n_consumed);
}

static
PyObject *
hashbits_print_stop_tags(khmer_KHashbits_Object * me, PyObject * args)
{
    Hashbits * hashbits = me->hashbits;

    const char * filename = NULL;

    if (!PyArg_ParseTuple(args, "s", &filename)) {
        return NULL;
    }

    hashbits->print_stop_tags(filename);

    Py_RETURN_NONE;
}

static
PyObject *
hashbits_print_tagset(khmer_KHashbits_Object * me, PyObject * args)
{
    Hashbits * hashbits = me->hashbits;

    const char * filename = NULL;

    if (!PyArg_ParseTuple(args, "s", &filename)) {
        return NULL;
    }

    hashbits->print_tagset(filename);

    Py_RETURN_NONE;
}

static
PyObject *
hashbits_load_stop_tags(khmer_KHashbits_Object * me, PyObject * args)
{
    Hashbits * hashbits = me->hashbits;

    const char * filename = NULL;
    PyObject * clear_tags_o = NULL;

    if (!PyArg_ParseTuple(args, "s|O", &filename, &clear_tags_o)) {
        return NULL;
    }

    bool clear_tags = true;
    if (clear_tags_o && !PyObject_IsTrue(clear_tags_o)) {
        clear_tags = false;
    }


    try {
        hashbits->load_stop_tags(filename, clear_tags);
    } catch (khmer_file_exception &e) {
        PyErr_SetString(PyExc_IOError, e.what());
        return NULL;
    }

    Py_RETURN_NONE;
}


static
PyObject *
hashbits_save_stop_tags(khmer_KHashbits_Object * me, PyObject * args)
{
    Hashbits * hashbits = me->hashbits;

    const char * filename = NULL;

    if (!PyArg_ParseTuple(args, "s", &filename)) {
        return NULL;
    }

    hashbits->save_stop_tags(filename);

    Py_RETURN_NONE;
}

static
PyObject *
hashbits_traverse_from_tags(khmer_KHashbits_Object * me, PyObject * args)
{
    Hashbits * hashbits = me->hashbits;

    khmer_KCountingHash_Object * counting_o = NULL;
    unsigned int distance, threshold, frequency;

    if (!PyArg_ParseTuple(args, "O!III", &khmer_KCountingHash_Type, &counting_o,
                          &distance, &threshold, &frequency)) {
        return NULL;
    }

    hashbits->traverse_from_tags(distance, threshold, frequency,
                                 * counting_o->counting);

    Py_RETURN_NONE;
}

static
PyObject *
hashbits_repartition_largest_partition(khmer_KHashbits_Object * me,
                                       PyObject * args)
{
    Hashbits * hashbits = me->hashbits;

    khmer_KCountingHash_Object * counting_o = NULL;
    PyObject * subset_o = NULL;
    unsigned int distance, threshold, frequency;

    if (!PyArg_ParseTuple(args, "OO!III", &subset_o, &khmer_KCountingHash_Type,
                          &counting_o, &distance, &threshold, &frequency)) {
        return NULL;
    }

    SubsetPartition * subset_p;
    if (subset_o != Py_None) {
        subset_p = (SubsetPartition *) PyCObject_AsVoidPtr(subset_o);
    } else {
        subset_p = hashbits->partition;
    }

    CountingHash * counting = counting_o->counting;

    unsigned long next_largest = subset_p->repartition_largest_partition(distance,
                                 threshold, frequency, *counting);

    return PyLong_FromLong(next_largest);
}

static
PyObject *
hashbits_get(khmer_KHashbits_Object * me, PyObject * args)
{
    Hashbits * hashbits = me->hashbits;

    PyObject * arg;

    if (!PyArg_ParseTuple(args, "O", &arg)) {
        return NULL;
    }

    unsigned long count = 0;

    if (PyInt_Check(arg)) {
        long pos = PyInt_AsLong(arg);
        count = hashbits->get_count((unsigned int) pos);
    } else if (PyBytes_Check(arg)) {
        std::string s = PyBytes_AsString(arg);

        if (strlen(s.c_str()) < hashbits->ksize()) {
            PyErr_SetString(PyExc_ValueError,
                            "string length must equal the presence table k-mer size");
            return NULL;
        }

        count = hashbits->get_count(s.c_str());
    } else {
        PyErr_SetString(PyExc_ValueError, "must pass in an int or string");
        return NULL;
    }

    return PyLong_FromLong(count);
}

static
PyObject *
hashbits_calc_connected_graph_size(khmer_KHashbits_Object * me, PyObject * args)
{
    Hashbits * hashbits = me->hashbits;

    const char * _kmer;
    unsigned int max_size = 0;
    PyObject * break_on_circum_o = NULL;
    if (!PyArg_ParseTuple(args, "s|IO", &_kmer, &max_size, &break_on_circum_o)) {
        return NULL;
    }

    bool break_on_circum = false;
    if (break_on_circum_o && PyObject_IsTrue(break_on_circum_o)) {
        break_on_circum = true;
    }

    unsigned long long size = 0;

    Py_BEGIN_ALLOW_THREADS
    SeenSet keeper;
    hashbits->calc_connected_graph_size(_kmer, size, keeper, max_size,
                                        break_on_circum);
    Py_END_ALLOW_THREADS

    return PyLong_FromUnsignedLongLong(size);
}

static
PyObject *
hashbits_kmer_degree(khmer_KHashbits_Object * me, PyObject * args)
{
    Hashbits * hashbits = me->hashbits;

    const char * kmer_s = NULL;

    if (!PyArg_ParseTuple(args, "s", &kmer_s)) {
        return NULL;
    }

    return PyLong_FromLong(hashbits->kmer_degree(kmer_s));
}

static
PyObject *
hashbits_trim_on_stoptags(khmer_KHashbits_Object * me, PyObject * args)
{
    Hashbits * hashbits = me->hashbits;

    const char * seq = NULL;

    if (!PyArg_ParseTuple(args, "s", &seq)) {
        return NULL;
    }

    size_t trim_at;
    Py_BEGIN_ALLOW_THREADS

    trim_at = hashbits->trim_on_stoptags(seq);

    Py_END_ALLOW_THREADS;

    PyObject * trim_seq = PyBytes_FromStringAndSize(seq, trim_at);
    if (trim_seq == NULL) {
        return NULL;
    }
    PyObject * ret = Py_BuildValue("Ok", trim_seq, (unsigned long) trim_at);
    Py_DECREF(trim_seq);

    return ret;
}

static
PyObject *
hashbits_identify_stoptags_by_position(khmer_KHashbits_Object * me,
                                       PyObject * args)
{
    Hashbits * hashbits = me->hashbits;

    const char * seq = NULL;

    if (!PyArg_ParseTuple(args, "s", &seq)) {
        return NULL;
    }

    std::vector<unsigned int> posns;
    Py_BEGIN_ALLOW_THREADS

    hashbits->identify_stop_tags_by_position(seq, posns);

    Py_END_ALLOW_THREADS;

    PyObject * x = PyList_New(posns.size());

    for (unsigned int i = 0; i < posns.size(); i++) {
        PyList_SET_ITEM(x, i, Py_BuildValue("I", posns[i]));
    }

    return x;
}

static
PyObject *
hashbits_do_subset_partition(khmer_KHashbits_Object * me, PyObject * args)
{
    Hashbits * hashbits = me->hashbits;

    HashIntoType start_kmer = 0, end_kmer = 0;
    PyObject * break_on_stop_tags_o = NULL;
    PyObject * stop_big_traversals_o = NULL;

    if (!PyArg_ParseTuple(args, "|KKOO", &start_kmer, &end_kmer,
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
        subset_p = new SubsetPartition(hashbits);
        subset_p->do_partition(start_kmer, end_kmer, break_on_stop_tags,
                               stop_big_traversals);
        Py_END_ALLOW_THREADS
    } catch (_khmer_signal &e) {
        return NULL;
    } catch (std::bad_alloc &e) {
        return PyErr_NoMemory();
    }

    return PyCObject_FromVoidPtr(subset_p, free_subset_partition_info);
}

static
PyObject *
hashbits_join_partitions_by_path(khmer_KHashbits_Object * me, PyObject * args)
{
    Hashbits * hashbits = me->hashbits;

    const char * sequence = NULL;
    if (!PyArg_ParseTuple(args, "s", &sequence)) {
        return NULL;
    }

    hashbits->partition->join_partitions_by_path(sequence);

    Py_RETURN_NONE;
}

static
PyObject *
hashbits_merge_subset(khmer_KHashbits_Object * me, PyObject * args)
{
    Hashbits * hashbits = me->hashbits;

    PyObject * subset_obj;
    if (!PyArg_ParseTuple(args, "O", &subset_obj)) {
        return NULL;
    }

    if (!PyCObject_Check(subset_obj)) {
        PyErr_SetString( PyExc_ValueError, "invalid subset");
        return NULL;
    }

    SubsetPartition * subset_p;
    subset_p = (SubsetPartition *) PyCObject_AsVoidPtr(subset_obj);

    hashbits->partition->merge(subset_p);

    Py_RETURN_NONE;
}

static
PyObject *
hashbits_merge_from_disk(khmer_KHashbits_Object * me, PyObject * args)
{
    Hashbits * hashbits = me->hashbits;

    const char * filename = NULL;
    if (!PyArg_ParseTuple(args, "s", &filename)) {
        return NULL;
    }

    try {
        hashbits->partition->merge_from_disk(filename);
    } catch (khmer_file_exception &e) {
        PyErr_SetString(PyExc_IOError, e.what());
        return NULL;
    }

    Py_RETURN_NONE;
}

static
PyObject *
hashbits_consume_fasta(khmer_KHashbits_Object * me, PyObject * args)
{
    Hashbits * hashbits = me->hashbits;

    const char * filename;

    if (!PyArg_ParseTuple(args, "s", &filename)) {
        return NULL;
    }

// call the C++ function, and trap signals => Python

    unsigned long long n_consumed = 0;
    unsigned int total_reads = 0;

    try {
        hashbits->consume_fasta(filename, total_reads, n_consumed);
    } catch (_khmer_signal &e) {
        PyErr_SetString(PyExc_IOError, e.get_message().c_str());
        return NULL;
    } catch (khmer_file_exception &e) {
        PyErr_SetString(PyExc_IOError, e.what());
        return NULL;
    }

    return Py_BuildValue("IK", total_reads, n_consumed);
}

static
PyObject *
hashbits_consume_fasta_with_reads_parser(khmer_KHashbits_Object * me,
        PyObject * args)
{
    Hashbits * hashbits = me->hashbits;

    PyObject * rparser_obj = NULL;

    if (!PyArg_ParseTuple(
                args, "O", &rparser_obj)) {
        return NULL;
    }

    read_parsers:: IParser * rparser =
        _PyObject_to_khmer_ReadParser( rparser_obj );

// call the C++ function, and trap signals => Python
    unsigned long long  n_consumed  = 0;
    unsigned int          total_reads = 0;
    char const * exc = NULL;
    Py_BEGIN_ALLOW_THREADS
    try {
        hashbits->consume_fasta(rparser, total_reads, n_consumed);
    } catch (_khmer_signal &e) {
        exc = e.get_message().c_str();
    }

    Py_END_ALLOW_THREADS
    if (exc != NULL) {
        PyErr_SetString(PyExc_IOError, exc);
        return NULL;
    }

    return Py_BuildValue("IK", total_reads, n_consumed);
}

static
PyObject *
hashbits_consume_fasta_and_traverse(khmer_KHashbits_Object * me,
                                    PyObject * args)
{
    Hashbits * hashbits = me->hashbits;

    const char * filename;
    unsigned int radius, big_threshold, transfer_threshold;
    khmer_KCountingHash_Object * counting_o = NULL;

    if (!PyArg_ParseTuple(args, "sIIIO!", &filename,
                          &radius, &big_threshold, &transfer_threshold,
                          &khmer_KCountingHash_Type, &counting_o)) {
        return NULL;
    }

    CountingHash * counting = counting_o->counting;

    hashbits->consume_fasta_and_traverse(filename, radius, big_threshold,
                                         transfer_threshold, *counting);


    Py_RETURN_NONE;
}

void sig(unsigned int total_reads, unsigned int n_consumed)
{
    std::cout << total_reads << " " << n_consumed << std::endl;
}

static
PyObject *
hashbits_consume_fasta_and_tag(khmer_KHashbits_Object * me, PyObject * args)
{
    Hashbits * hashbits = me->hashbits;

    const char * filename;

    if (!PyArg_ParseTuple(args, "s", &filename)) {
        return NULL;
    }

    // call the C++ function, and trap signals => Python

    unsigned long long n_consumed;
    unsigned int total_reads;

    try {
        hashbits->consume_fasta_and_tag(filename, total_reads, n_consumed);
    } catch (_khmer_signal &e) {
        PyErr_SetString(PyExc_IOError, e.get_message().c_str());
        return NULL;
    } catch (khmer_file_exception &e) {
        PyErr_SetString(PyExc_IOError, e.what());
        return NULL;
    }

    return Py_BuildValue("IK", total_reads, n_consumed);
}

static
PyObject *
hashbits_consume_fasta_and_tag_with_reads_parser(khmer_KHashbits_Object * me,
        PyObject * args)
{
    Hashbits * hashbits = me->hashbits;

    python::khmer_ReadParser_Object * rparser_obj = NULL;

    if (!PyArg_ParseTuple( args, "O!", &python::khmer_ReadParser_Type,
                           &rparser_obj)) {
        return NULL;
    }

    read_parsers:: IParser * rparser = rparser_obj-> parser;

    // call the C++ function, and trap signals => Python
    unsigned long long  n_consumed  = 0;
    unsigned int          total_reads = 0;
    char const * exc = NULL;
    Py_BEGIN_ALLOW_THREADS
    try {
        hashbits->consume_fasta_and_tag(
            rparser, total_reads, n_consumed
        );
    } catch (_khmer_signal &e) {
        exc = e.get_message().c_str();
    } catch (khmer::read_parsers::NoMoreReadsAvailable &e) {
        exc = e.what();
    }
    Py_END_ALLOW_THREADS
    if (exc != NULL) {
        PyErr_SetString(PyExc_IOError, exc);
        return NULL;
    }

    return Py_BuildValue("IK", total_reads, n_consumed);
}

static
PyObject *
hashbits_consume_fasta_and_tag_with_stoptags(khmer_KHashbits_Object * me,
        PyObject * args)
{
    Hashbits * hashbits = me->hashbits;

    const char * filename;

    if (!PyArg_ParseTuple(args, "s", &filename)) {
        return NULL;
    }

    // call the C++ function, and trap signals => Python

    unsigned long long n_consumed;
    unsigned int total_reads;

    try {
        hashbits->consume_fasta_and_tag_with_stoptags(filename,
                total_reads, n_consumed);
    } catch (_khmer_signal &e) {
        PyErr_SetString(PyExc_IOError, e.get_message().c_str());
        return NULL;
    } catch (khmer_file_exception &e) {
        PyErr_SetString(PyExc_IOError, e.what());
        return NULL;
    }

    return Py_BuildValue("IK", total_reads, n_consumed);
}

static
PyObject *
hashbits_consume_partitioned_fasta(khmer_KHashbits_Object * me, PyObject * args)
{
    Hashbits * hashbits = me->hashbits;

    const char * filename;

    if (!PyArg_ParseTuple(args, "s", &filename)) {
        return NULL;
    }

    // call the C++ function, and trap signals => Python

    unsigned long long n_consumed;
    unsigned int total_reads;

    try {
        hashbits->consume_partitioned_fasta(filename, total_reads, n_consumed);
    } catch (_khmer_signal &e) {
        PyErr_SetString(PyExc_IOError, e.get_message().c_str());
        return NULL;
    } catch (khmer_file_exception &e) {
        PyErr_SetString(PyExc_IOError, e.what());
        return NULL;
    }

    return Py_BuildValue("IK", total_reads, n_consumed);
}

static
PyObject *
hashbits_find_all_tags(khmer_KHashbits_Object * me, PyObject * args)
{
    Hashbits * hashbits = me->hashbits;

    const char * kmer_s = NULL;

    if (!PyArg_ParseTuple(args, "s", &kmer_s)) {
        return NULL;
    }

    if (strlen(kmer_s) != hashbits->ksize()) {
        PyErr_SetString( PyExc_ValueError,
                         "k-mer size must equal the k-mer size of the presence table");
        return NULL;
    }

    _pre_partition_info * ppi = NULL;

    Py_BEGIN_ALLOW_THREADS

    HashIntoType kmer, kmer_f, kmer_r;
    kmer = _hash(kmer_s, hashbits->ksize(), kmer_f, kmer_r);

    try {
        ppi = new _pre_partition_info(kmer);
    } catch (std::bad_alloc &e) {
        return PyErr_NoMemory();
    }
    hashbits->partition->find_all_tags(kmer_f, kmer_r, ppi->tagged_kmers,
                                       hashbits->all_tags);
    hashbits->add_kmer_to_tags(kmer);

    Py_END_ALLOW_THREADS

    return PyCObject_FromVoidPtr(ppi, free_pre_partition_info);
}

static
PyObject *
hashbits_assign_partition_id(khmer_KHashbits_Object * me, PyObject * args)
{
    Hashbits * hashbits = me->hashbits;

    PyObject * ppi_obj;
    if (!PyArg_ParseTuple(args, "O", &ppi_obj)) {
        return NULL;
    }

    if (!PyCObject_Check(ppi_obj)) {
        PyErr_SetString( PyExc_ValueError, "invalid pre_partition_info");
        return NULL;
    }

    _pre_partition_info * ppi;
    ppi = (_pre_partition_info *) PyCObject_AsVoidPtr(ppi_obj);

    PartitionID p;
    p = hashbits->partition->assign_partition_id(ppi->kmer,
            ppi->tagged_kmers);

    return PyLong_FromLong(p);
}

static
PyObject *
hashbits_add_tag(khmer_KHashbits_Object * me, PyObject * args)
{
    Hashbits * hashbits = me->hashbits;

    const char * kmer_s = NULL;
    if (!PyArg_ParseTuple(args, "s", &kmer_s)) {
        return NULL;
    }

    HashIntoType kmer = _hash(kmer_s, hashbits->ksize());
    hashbits->add_tag(kmer);

    Py_RETURN_NONE;
}

static
PyObject *
hashbits_add_stop_tag(khmer_KHashbits_Object * me, PyObject * args)
{
    Hashbits * hashbits = me->hashbits;

    const char * kmer_s = NULL;
    if (!PyArg_ParseTuple(args, "s", &kmer_s)) {
        return NULL;
    }

    HashIntoType kmer = _hash(kmer_s, hashbits->ksize());
    hashbits->add_stop_tag(kmer);

    Py_RETURN_NONE;
}

static
PyObject *
hashbits_get_stop_tags(khmer_KHashbits_Object * me, PyObject * args)
{
    Hashbits * hashbits = me->hashbits;

    if (!PyArg_ParseTuple(args, "")) {
        return NULL;
    }

    WordLength k = hashbits->ksize();
    SeenSet::const_iterator si;

    PyObject * x = PyList_New(hashbits->stop_tags.size());
    unsigned long long i = 0;
    for (si = hashbits->stop_tags.begin(); si != hashbits->stop_tags.end(); si++) {
        std::string s = _revhash(*si, k);
        PyList_SET_ITEM(x, i, Py_BuildValue("s", s.c_str()));
        i++;
    }

    return x;
}

static
PyObject *
hashbits_get_tagset(khmer_KHashbits_Object * me, PyObject * args)
{
    Hashbits * hashbits = me->hashbits;

    if (!PyArg_ParseTuple(args, "")) {
        return NULL;
    }

    WordLength k = hashbits->ksize();
    SeenSet::const_iterator si;

    PyObject * x = PyList_New(hashbits->all_tags.size());
    unsigned long long i = 0;
    for (si = hashbits->all_tags.begin(); si != hashbits->all_tags.end(); si++) {
        std::string s = _revhash(*si, k);
        PyList_SET_ITEM(x, i, Py_BuildValue("s", s.c_str()));
        i++;
    }

    return x;
}

static
PyObject *
hashbits_output_partitions(khmer_KHashbits_Object * me, PyObject * args)
{
    Hashbits * hashbits = me->hashbits;

    const char * filename = NULL;
    const char * output = NULL;
    PyObject * output_unassigned_o = NULL;

    if (!PyArg_ParseTuple(args, "ss|O", &filename, &output,
                          &output_unassigned_o)) {
        return NULL;
    }

    bool output_unassigned = false;
    if (output_unassigned_o != NULL && PyObject_IsTrue(output_unassigned_o)) {
        output_unassigned = true;
    }

    size_t n_partitions = 0;

    try {
        SubsetPartition * subset_p = hashbits->partition;
        n_partitions = subset_p->output_partitioned_file(filename,
                       output,
                       output_unassigned);
    } catch (_khmer_signal &e) {
        PyErr_SetString(PyExc_IOError, e.get_message().c_str());
        return NULL;
    } catch (khmer_file_exception &e) {
        PyErr_SetString(PyExc_IOError, e.what());
        return NULL;
    }

    return PyLong_FromLong(n_partitions);
}

static
PyObject *
hashbits_find_unpart(khmer_KHashbits_Object * me, PyObject * args)
{
    Hashbits * hashbits = me->hashbits;

    const char * filename = NULL;
    PyObject * traverse_o = NULL;
    PyObject * stop_big_traversals_o = NULL;

    if (!PyArg_ParseTuple(args, "sOO", &filename, &traverse_o,
                          &stop_big_traversals_o)) {
        return NULL;
    }

    bool traverse = PyObject_IsTrue(traverse_o);
    bool stop_big_traversals = PyObject_IsTrue(stop_big_traversals_o);
    unsigned int n_singletons = 0;

    try {
        SubsetPartition * subset_p = hashbits->partition;
        n_singletons = subset_p->find_unpart(filename, traverse,
                                             stop_big_traversals);
    } catch (_khmer_signal &e) {
        return NULL;
    }

    return PyLong_FromLong(n_singletons);

    // Py_INCREF(Py_None);
    // return Py_None;
}

static
PyObject *
hashbits_filter_if_present(khmer_KHashbits_Object * me, PyObject * args)
{
    Hashbits * hashbits = me->hashbits;

    const char * filename = NULL;
    const char * output = NULL;

    if (!PyArg_ParseTuple(args, "ss", &filename, &output)) {
        return NULL;
    }

    try {
        hashbits->filter_if_present(filename, output);
    } catch (_khmer_signal &e) {
        return NULL;
    }

    Py_RETURN_NONE;
}

static
PyObject *
hashbits_save_partitionmap(khmer_KHashbits_Object * me, PyObject * args)
{
    Hashbits * hashbits = me->hashbits;

    const char * filename = NULL;

    if (!PyArg_ParseTuple(args, "s", &filename)) {
        return NULL;
    }

    hashbits->partition->save_partitionmap(filename);

    Py_RETURN_NONE;
}

static
PyObject *
hashbits_load_partitionmap(khmer_KHashbits_Object * me, PyObject * args)
{
    Hashbits * hashbits = me->hashbits;

    const char * filename = NULL;

    if (!PyArg_ParseTuple(args, "s", &filename)) {
        return NULL;
    }

    hashbits->partition->load_partitionmap(filename);

    Py_RETURN_NONE;
}

static
PyObject *
hashbits__validate_partitionmap(khmer_KHashbits_Object * me, PyObject * args)
{
    Hashbits * hashbits = me->hashbits;

    if (!PyArg_ParseTuple(args, "")) {
        return NULL;
    }

    hashbits->partition->_validate_pmap();

    Py_RETURN_NONE;
}

static
PyObject *
hashbits_count_partitions(khmer_KHashbits_Object * me, PyObject * args)
{
    Hashbits * hashbits = me->hashbits;

    if (!PyArg_ParseTuple(args, "")) {
        return NULL;
    }

    size_t n_partitions = 0, n_unassigned = 0;
    hashbits->partition->count_partitions(n_partitions, n_unassigned);

    return Py_BuildValue("nn", (Py_ssize_t) n_partitions,
                         (Py_ssize_t) n_unassigned);
}

static
PyObject *
hashbits_subset_count_partitions(khmer_KHashbits_Object * me, PyObject * args)
{
    PyObject * subset_obj = NULL;

    if (!PyArg_ParseTuple(args, "O", &subset_obj)) {
        return NULL;
    }

    SubsetPartition * subset_p;
    subset_p = (SubsetPartition *) PyCObject_AsVoidPtr(subset_obj);

    size_t n_partitions = 0, n_unassigned = 0;
    subset_p->count_partitions(n_partitions, n_unassigned);

    return Py_BuildValue("nn", (Py_ssize_t) n_partitions,
                         (Py_ssize_t) n_unassigned);
}

static
PyObject *
hashbits_subset_partition_size_distribution(khmer_KHashbits_Object * me,
        PyObject * args)
{
    PyObject * subset_obj = NULL;
    if (!PyArg_ParseTuple(args, "O", &subset_obj)) {
        return NULL;
    }

    SubsetPartition * subset_p;
    subset_p = (SubsetPartition *) PyCObject_AsVoidPtr(subset_obj);

    PartitionCountDistribution d;

    unsigned int n_unassigned = 0;
    subset_p->partition_size_distribution(d, n_unassigned);

    PyObject * x = PyList_New(d.size());
    if (x == NULL) {
        return NULL;
    }
    PartitionCountDistribution::iterator di;

    unsigned int i;
    for (i = 0, di = d.begin(); di != d.end(); di++, i++) {
        PyObject * value =  Py_BuildValue("KK", di->first, di->second);
        if (value == NULL) {
            Py_DECREF(x);
            return NULL;
        }
        PyList_SET_ITEM(x, i, value);
    }
    if (!(i == d.size())) {
        throw khmer_exception();
    }

    PyObject * returnValue = Py_BuildValue("NI", x, n_unassigned);
    if (returnValue == NULL) {
        Py_DECREF(x);
        return NULL;
    }
    return returnValue;
}

static
PyObject *
hashbits_load(khmer_KHashbits_Object * me, PyObject * args)
{
    Hashbits * hashbits = me->hashbits;

    const char * filename = NULL;

    if (!PyArg_ParseTuple(args, "s", &filename)) {
        return NULL;
    }

    try {
        hashbits->load(filename);
    } catch (khmer_file_exception &e) {
        PyErr_SetString(PyExc_IOError, e.what());
        return NULL;
    }

    Py_RETURN_NONE;
}

static
PyObject *
hashbits_save(khmer_KHashbits_Object * me, PyObject * args)
{
    Hashbits * hashbits = me->hashbits;

    const char * filename = NULL;

    if (!PyArg_ParseTuple(args, "s", &filename)) {
        return NULL;
    }

    hashbits->save(filename);

    Py_RETURN_NONE;
}

static
PyObject *
hashbits_load_tagset(khmer_KHashbits_Object * me, PyObject * args)
{
    Hashbits * hashbits = me->hashbits;

    const char * filename = NULL;
    PyObject * clear_tags_o = NULL;

    if (!PyArg_ParseTuple(args, "s|O", &filename, &clear_tags_o)) {
        return NULL;
    }

    bool clear_tags = true;
    if (clear_tags_o && !PyObject_IsTrue(clear_tags_o)) {
        clear_tags = false;
    }

    try {
        hashbits->load_tagset(filename, clear_tags);
    } catch (khmer_file_exception &e) {
        PyErr_SetString(PyExc_IOError, e.what());
        return NULL;
    }

    Py_RETURN_NONE;
}

static
PyObject *
hashbits_save_tagset(khmer_KHashbits_Object * me, PyObject * args)
{
    Hashbits * hashbits = me->hashbits;

    const char * filename = NULL;

    if (!PyArg_ParseTuple(args, "s", &filename)) {
        return NULL;
    }

    hashbits->save_tagset(filename);

    Py_RETURN_NONE;
}

static
PyObject *
hashbits_save_subset_partitionmap(khmer_KHashbits_Object * me, PyObject * args)
{
    const char * filename = NULL;
    PyObject * subset_obj = NULL;

    if (!PyArg_ParseTuple(args, "Os", &subset_obj, &filename)) {
        return NULL;
    }

    SubsetPartition * subset_p;
    subset_p = (SubsetPartition *) PyCObject_AsVoidPtr(subset_obj);

    Py_BEGIN_ALLOW_THREADS

    subset_p->save_partitionmap(filename);

    Py_END_ALLOW_THREADS

    Py_RETURN_NONE;
}

static
PyObject *
hashbits_load_subset_partitionmap(khmer_KHashbits_Object * me, PyObject * args)
{
    Hashbits * hashbits = me->hashbits;

    const char * filename = NULL;

    if (!PyArg_ParseTuple(args, "s", &filename)) {
        return NULL;
    }

    SubsetPartition * subset_p;
    try {
        subset_p = new SubsetPartition(hashbits);
    } catch (std::bad_alloc &e) {
        return PyErr_NoMemory();
    }

    bool fail = false;
    std::string err;

    Py_BEGIN_ALLOW_THREADS

    try {
        subset_p->load_partitionmap(filename);
    } catch (khmer_file_exception &e) {
        fail = true;
        err = e.what();
    }

    Py_END_ALLOW_THREADS

    if (fail) {
        PyErr_SetString(PyExc_IOError, err.c_str());
        delete subset_p;
        return NULL;
    } else {
        return PyCObject_FromVoidPtr(subset_p, free_subset_partition_info);
    }
}

static
PyObject *
hashbits__set_tag_density(khmer_KHashbits_Object * me, PyObject * args)
{
    Hashbits * hashbits = me->hashbits;

    unsigned int d;
    if (!PyArg_ParseTuple(args, "I", &d)) {
        return NULL;
    }

    hashbits->_set_tag_density(d);

    Py_RETURN_NONE;
}

static
PyObject *
hashbits__get_tag_density(khmer_KHashbits_Object * me, PyObject * args)
{
    Hashbits * hashbits = me->hashbits;

    if (!PyArg_ParseTuple(args, "")) {
        return NULL;
    }

    unsigned int d = hashbits->_get_tag_density();

    return PyLong_FromLong(d);
}

static
PyObject *
hashbits__validate_subset_partitionmap(khmer_KHashbits_Object * me,
                                       PyObject * args)
{
    PyObject * subset_obj = NULL;

    if (!PyArg_ParseTuple(args, "O", &subset_obj)) {
        return NULL;
    }

    SubsetPartition * subset_p;
    subset_p = (SubsetPartition *) PyCObject_AsVoidPtr(subset_obj);
    subset_p->_validate_pmap();

    Py_RETURN_NONE;
}

static
PyObject *
hashbits_set_partition_id(khmer_KHashbits_Object * me, PyObject * args)
{
    Hashbits * hashbits = me->hashbits;

    const char * kmer = NULL;
    PartitionID p = 0;

    if (!PyArg_ParseTuple(args, "sI", &kmer, &p)) {
        return NULL;
    }

    hashbits->partition->set_partition_id(kmer, p);

    Py_RETURN_NONE;
}

static
PyObject *
hashbits_join_partitions(khmer_KHashbits_Object * me, PyObject * args)
{
    Hashbits * hashbits = me->hashbits;

    PartitionID p1 = 0, p2 = 0;

    if (!PyArg_ParseTuple(args, "II", &p1, &p2)) {
        return NULL;
    }

    p1 = hashbits->partition->join_partitions(p1, p2);

    return PyLong_FromLong(p1);
}

static
PyObject *
hashbits_get_partition_id(khmer_KHashbits_Object * me, PyObject * args)
{
    Hashbits * hashbits = me->hashbits;

    const char * kmer = NULL;

    if (!PyArg_ParseTuple(args, "s", &kmer)) {
        return NULL;
    }

    PartitionID partition_id;
    partition_id = hashbits->partition->get_partition_id(kmer);

    return PyLong_FromLong(partition_id);
}

static
PyObject *
hashbits_is_single_partition(khmer_KHashbits_Object * me, PyObject * args)
{
    Hashbits * hashbits = me->hashbits;

    const char * seq = NULL;

    if (!PyArg_ParseTuple(args, "s", &seq)) {
        return NULL;
    }

    bool v = hashbits->partition->is_single_partition(seq);

    PyObject * val;
    if (v) {
        val = Py_True;
    } else {
        val = Py_False;
    }
    Py_INCREF(val);

    return val;
}

static
PyObject *
hashbits_divide_tags_into_subsets(khmer_KHashbits_Object * me, PyObject * args)
{
    Hashbits * hashbits = me->hashbits;

    unsigned int subset_size = 0;

    if (!PyArg_ParseTuple(args, "I", &subset_size)) {
        return NULL;
    }

    SeenSet divvy;
    hashbits->divide_tags_into_subsets(subset_size, divvy);

    PyObject * x = PyList_New(divvy.size());
    unsigned int i = 0;
    for (SeenSet::const_iterator si = divvy.begin(); si != divvy.end();
            si++, i++) {
        PyList_SET_ITEM(x, i, PyLong_FromUnsignedLongLong(*si));
    }

    return x;
}

static
PyObject *
hashbits_count_kmers_within_radius(khmer_KHashbits_Object * me, PyObject * args)
{
    Hashbits * hashbits = me->hashbits;

    const char * kmer = NULL;
    unsigned int radius = 0;
    unsigned int max_count = 0;

    if (!PyArg_ParseTuple(args, "sI|I", &kmer, &radius, &max_count)) {
        return NULL;
    }

    unsigned int n;

    Py_BEGIN_ALLOW_THREADS

    HashIntoType kmer_f, kmer_r;
    _hash(kmer, hashbits->ksize(), kmer_f, kmer_r);
    n = hashbits->count_kmers_within_radius(kmer_f, kmer_r, radius,
                                            max_count);

    Py_END_ALLOW_THREADS

    return PyLong_FromUnsignedLong(n);
}

static
PyObject *
hashbits_get_ksize(khmer_KHashbits_Object * me, PyObject * args)
{
    Hashbits * hashbits = me->hashbits;

    if (!PyArg_ParseTuple(args, "")) {
        return NULL;
    }

    unsigned int k = hashbits->ksize();

    return PyLong_FromLong(k);
}


static
PyObject *
hashbits_get_hashsizes(khmer_KHashbits_Object * me, PyObject * args)
{
    Hashbits * hashbits = me->hashbits;

    if (!PyArg_ParseTuple(args, "")) {
        return NULL;
    }

    std::vector<HashIntoType> ts = hashbits->get_tablesizes();

    PyObject * x = PyList_New(ts.size());
    for (size_t i = 0; i < ts.size(); i++) {
        PyList_SET_ITEM(x, i, PyLong_FromUnsignedLongLong(ts[i]));
    }

    return x;
}

static
PyObject *
hashbits_extract_unique_paths(khmer_KHashbits_Object * me, PyObject * args)
{
    Hashbits * hashbits = me->hashbits;

    const char * sequence = NULL;
    unsigned int min_length = 0;
    float min_unique_f = 0;
    if (!PyArg_ParseTuple(args, "sIf", &sequence, &min_length, &min_unique_f)) {
        return NULL;
    }

    std::vector<std::string> results;
    hashbits->extract_unique_paths(sequence, min_length, min_unique_f, results);

    PyObject * x = PyList_New(results.size());
    if (x == NULL) {
        return NULL;
    }

    for (unsigned int i = 0; i < results.size(); i++) {
        PyList_SET_ITEM(x, i, PyBytes_FromString(results[i].c_str()));
    }

    return x;
}

static
PyObject *
hashbits_get_median_count(khmer_KHashbits_Object * me, PyObject * args)
{
    Hashbits * hashbits = me->hashbits;

    const char * long_str;

    if (!PyArg_ParseTuple(args, "s", &long_str)) {
        return NULL;
    }

    if (strlen(long_str) < hashbits->ksize()) {
        PyErr_SetString(PyExc_ValueError,
                        "string length must >= the hashtable k-mer size");
        return NULL;
    }

    BoundedCounterType med = 0;
    float average = 0, stddev = 0;

    hashbits->get_median_count(long_str, med, average, stddev);

    return Py_BuildValue("iff", med, average, stddev);
}

static PyMethodDef khmer_hashbits_methods[] = {
    { "extract_unique_paths", (PyCFunction)hashbits_extract_unique_paths, METH_VARARGS, "" },
    { "ksize", (PyCFunction)hashbits_get_ksize, METH_VARARGS, "" },
    { "hashsizes", (PyCFunction)hashbits_get_hashsizes, METH_VARARGS, "" },
    { "n_occupied", (PyCFunction)hashbits_n_occupied, METH_VARARGS, "Count the number of occupied bins" },
    { "n_unique_kmers", (PyCFunction)hashbits_n_unique_kmers,  METH_VARARGS, "Count the number of unique kmers" },
    { "count", (PyCFunction)hashbits_count, METH_VARARGS, "Count the given kmer" },
    { "count_overlap", (PyCFunction)hashbits_count_overlap, METH_VARARGS, "Count overlap kmers in two datasets" },
    { "consume", (PyCFunction)hashbits_consume, METH_VARARGS, "Count all k-mers in the given string" },
    { "load_stop_tags", (PyCFunction)hashbits_load_stop_tags, METH_VARARGS, "" },
    { "save_stop_tags", (PyCFunction)hashbits_save_stop_tags, METH_VARARGS, "" },
    { "print_stop_tags", (PyCFunction)hashbits_print_stop_tags, METH_VARARGS, "" },
    { "print_tagset", (PyCFunction)hashbits_print_tagset, METH_VARARGS, "" },
    { "get", (PyCFunction)hashbits_get, METH_VARARGS, "Get the count for the given k-mer" },
    { "calc_connected_graph_size", (PyCFunction)hashbits_calc_connected_graph_size, METH_VARARGS, "" },
    { "kmer_degree", (PyCFunction)hashbits_kmer_degree, METH_VARARGS, "" },
    { "trim_on_stoptags", (PyCFunction)hashbits_trim_on_stoptags, METH_VARARGS, "" },
    { "identify_stoptags_by_position", (PyCFunction)hashbits_identify_stoptags_by_position, METH_VARARGS, "" },
    { "do_subset_partition", (PyCFunction)hashbits_do_subset_partition, METH_VARARGS, "" },
    { "find_all_tags", (PyCFunction)hashbits_find_all_tags, METH_VARARGS, "" },
    { "assign_partition_id", (PyCFunction)hashbits_assign_partition_id, METH_VARARGS, "" },
    { "output_partitions", (PyCFunction)hashbits_output_partitions, METH_VARARGS, "" },
    { "find_unpart", (PyCFunction)hashbits_find_unpart, METH_VARARGS, "" },
    { "filter_if_present", (PyCFunction)hashbits_filter_if_present, METH_VARARGS, "" },
    { "add_tag", (PyCFunction)hashbits_add_tag, METH_VARARGS, "" },
    { "add_stop_tag", (PyCFunction)hashbits_add_stop_tag, METH_VARARGS, "" },
    { "get_stop_tags", (PyCFunction)hashbits_get_stop_tags, METH_VARARGS, "" },
    { "get_tagset", (PyCFunction)hashbits_get_tagset, METH_VARARGS, "" },
    { "load", (PyCFunction)hashbits_load, METH_VARARGS, "" },
    { "save", (PyCFunction)hashbits_save, METH_VARARGS, "" },
    { "load_tagset", (PyCFunction)hashbits_load_tagset, METH_VARARGS, "" },
    { "save_tagset", (PyCFunction)hashbits_save_tagset, METH_VARARGS, "" },
    { "n_tags", (PyCFunction)hashbits_n_tags, METH_VARARGS, "" },
    { "divide_tags_into_subsets", (PyCFunction)hashbits_divide_tags_into_subsets, METH_VARARGS, "" },
    { "load_partitionmap", (PyCFunction)hashbits_load_partitionmap, METH_VARARGS, "" },
    { "save_partitionmap", (PyCFunction)hashbits_save_partitionmap, METH_VARARGS, "" },
    { "_validate_partitionmap", (PyCFunction)hashbits__validate_partitionmap, METH_VARARGS, "" },
    { "_get_tag_density", (PyCFunction)hashbits__get_tag_density, METH_VARARGS, "" },
    { "_set_tag_density", (PyCFunction)hashbits__set_tag_density, METH_VARARGS, "" },
    { "consume_fasta", (PyCFunction)hashbits_consume_fasta, METH_VARARGS, "Count all k-mers in a given file" },
    { "consume_fasta_with_reads_parser", (PyCFunction)hashbits_consume_fasta_with_reads_parser, METH_VARARGS, "Count all k-mers in a given file" },
    { "consume_fasta_and_tag", (PyCFunction)hashbits_consume_fasta_and_tag, METH_VARARGS, "Count all k-mers in a given file" },
    {
        "consume_fasta_and_tag_with_reads_parser", (PyCFunction)hashbits_consume_fasta_and_tag_with_reads_parser,
        METH_VARARGS, "Count all k-mers using a given reads parser"
    },
    { "consume_fasta_and_traverse", (PyCFunction)hashbits_consume_fasta_and_traverse, METH_VARARGS, "" },
    { "consume_fasta_and_tag_with_stoptags", (PyCFunction)hashbits_consume_fasta_and_tag_with_stoptags, METH_VARARGS, "Count all k-mers in a given file" },
    { "consume_partitioned_fasta", (PyCFunction)hashbits_consume_partitioned_fasta, METH_VARARGS, "Count all k-mers in a given file" },
    { "join_partitions_by_path", (PyCFunction)hashbits_join_partitions_by_path, METH_VARARGS, "" },
    { "merge_subset", (PyCFunction)hashbits_merge_subset, METH_VARARGS, "" },
    { "merge_subset_from_disk", (PyCFunction)hashbits_merge_from_disk, METH_VARARGS, "" },
    { "count_partitions", (PyCFunction)hashbits_count_partitions, METH_VARARGS, "" },
    { "subset_count_partitions", (PyCFunction)hashbits_subset_count_partitions, METH_VARARGS, "" },
    { "subset_partition_size_distribution", (PyCFunction)hashbits_subset_partition_size_distribution, METH_VARARGS, "" },
    { "save_subset_partitionmap", (PyCFunction)hashbits_save_subset_partitionmap, METH_VARARGS },
    { "load_subset_partitionmap", (PyCFunction)hashbits_load_subset_partitionmap, METH_VARARGS },
    { "_validate_subset_partitionmap", (PyCFunction)hashbits__validate_subset_partitionmap, METH_VARARGS, "" },
    { "set_partition_id", (PyCFunction)hashbits_set_partition_id, METH_VARARGS, "" },
    { "join_partitions", (PyCFunction)hashbits_join_partitions, METH_VARARGS, "" },
    { "get_partition_id", (PyCFunction)hashbits_get_partition_id, METH_VARARGS, "" },
    { "is_single_partition", (PyCFunction)hashbits_is_single_partition, METH_VARARGS, "" },
    { "count_kmers_within_radius", (PyCFunction)hashbits_count_kmers_within_radius, METH_VARARGS, "" },
    { "traverse_from_tags", (PyCFunction)hashbits_traverse_from_tags, METH_VARARGS, "" },
    { "repartition_largest_partition", (PyCFunction)hashbits_repartition_largest_partition, METH_VARARGS, "" },
    { "get_median_count", (PyCFunction)hashbits_get_median_count, METH_VARARGS, "Get the median, average, and stddev of the k-mer counts in the string" },
    {NULL, NULL, 0, NULL}           /* sentinel */
};

// __new__ for hashbits; necessary for proper subclassing
// This will essentially do what the old factory function did. Unlike many __new__
// methods, we take our arguments here, because there's no "uninitialized" hashbits
// object; we have to have k and the table sizes before creating the new objects
static PyObject* khmer_hashbits_new(PyTypeObject * type, PyObject * args,
                                    PyObject * kwds)
{
    khmer_KHashbits_Object * self;
    self = (khmer_KHashbits_Object *)type->tp_alloc(type, 0);

    if (self != NULL) {
        WordLength k = 0;
        PyListObject* sizes_list_o = NULL;

        if (!PyArg_ParseTuple(args, "bO!", &k, &PyList_Type, &sizes_list_o)) {
            Py_DECREF(self);
            return NULL;
        }

        std::vector<HashIntoType> sizes;
        Py_ssize_t sizes_list_o_length = PyList_GET_SIZE(sizes_list_o);
        for (Py_ssize_t i = 0; i < sizes_list_o_length; i++) {
            PyObject * size_o = PyList_GET_ITEM(sizes_list_o, i);
            if (PyLong_Check(size_o)) {
                sizes.push_back((HashIntoType) PyLong_AsUnsignedLongLong(size_o));
            } else if (PyInt_Check(size_o)) {
                sizes.push_back((HashIntoType) PyInt_AsLong(size_o));
            } else if (PyFloat_Check(size_o)) {
                sizes.push_back((HashIntoType) PyFloat_AS_DOUBLE(size_o));
            } else {
                Py_DECREF(self);
                PyErr_SetString(PyExc_TypeError,
                                "2nd argument must be a list of ints, longs, or floats");
                return NULL;
            }
        }

        try {
            self->hashbits = new Hashbits(k, sizes);
        } catch (std::bad_alloc &e) {
            return PyErr_NoMemory();
        }
    }
    return (PyObject *) self;
}

// there are no attributes that we need at this time, so we'll just return 0
static int khmer_hashbits_init(khmer_KHashbits_Object * self, PyObject * args,
                               PyObject * kwds)
{
    return 0;
}

#define is_hashbits_obj(v)  (Py_TYPE(v) == &khmer_KHashbits_Type)

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
subset_compare_partitions(khmer_KSubsetPartition_Object * me, PyObject * args)
{
    SubsetPartition * subset1_p = me->subset;

    PyObject * subset2_obj = NULL;
    PartitionID pid1, pid2; // @CTB ensure that these are unsigned?

    if (!PyArg_ParseTuple(args, "IOI",
                          &pid1, &subset2_obj, &pid2)) {
        return NULL;
    }

    khmer_KSubsetPartition_Object *other = (khmer_KSubsetPartition_Object *)
                                           subset2_obj;
    SubsetPartition * subset2_p = other->subset;

    unsigned int n_only1 = 0, n_only2 = 0, n_shared = 0;
    subset1_p->compare_to_partition(pid1, subset2_p, pid2,
                                    n_only1, n_only2, n_shared);

    return Py_BuildValue("III", n_only1, n_only2, n_shared);
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
    for (i = 0, di = d.begin(); di != d.end(); di++, i++) {
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
    for (mi = cm.begin(); mi != cm.end(); mi++) {
        if (mi->second >= min_size) {
            i++;
        }
    }

    PyObject * x = PyList_New(i);
    if (x == NULL) {
        return NULL;
    }

    // this should probably be a dict. @CTB
    for (i = 0, mi = cm.begin(); mi != cm.end(); mi++) {
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

    khmer_KCountingHash_Object * counting_o;

    if (!PyArg_ParseTuple(args, "O!", &khmer_KCountingHash_Type, &counting_o)) {
        return NULL;
    }

    PartitionCountMap cm;
    subset_p->partition_average_coverages(cm, counting_o -> counting);

    unsigned int i;
    PartitionCountMap::iterator mi;

    PyObject * x = PyList_New(cm.size());
    if (x == NULL) {
        return NULL;
    }

    // this should probably be a dict. @CTB
    for (i = 0, mi = cm.begin(); mi != cm.end(); mi++, i++) {
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
        "compare_partitions",
        (PyCFunction)subset_compare_partitions,
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

/////////////////
// LabelHash
/////////////////

// LabelHash addition
typedef struct {
    //PyObject_HEAD
    khmer_KHashbits_Object khashbits;
    LabelHash * labelhash;
} khmer_KLabelHash_Object;

static int khmer_labelhash_init(khmer_KLabelHash_Object * self, PyObject *args,
                                PyObject *kwds);
static PyObject * khmer_labelhash_new(PyTypeObject * type, PyObject *args,
                                      PyObject *kwds);

#define is_labelhash_obj(v)  (Py_TYPE(v) == &khmer_KLabelHash_Type)

//
// khmer_labelhash_dealloc -- clean up a labelhash object.
//

static void khmer_labelhash_dealloc(khmer_KLabelHash_Object * obj)
{
    delete obj->labelhash;
    obj->labelhash = NULL;

    Py_TYPE(obj)->tp_free((PyObject*)obj);
}

// a little weird; we don't actually want to call Hashbits' new method. Rather, we
// define our own new method, and redirect the base's hashbits object to point to our
// labelhash object
static PyObject * khmer_labelhash_new(PyTypeObject *type, PyObject *args,
                                      PyObject *kwds)
{
    khmer_KLabelHash_Object *self;
    self = (khmer_KLabelHash_Object*)type->tp_alloc(type, 0);

    if (self != NULL) {
        WordLength k = 0;
        PyListObject * sizes_list_o = NULL;

        if (!PyArg_ParseTuple(args, "bO!", &k, &PyList_Type, &sizes_list_o)) {
            Py_DECREF(self);
            return NULL;
        }

        std::vector<HashIntoType> sizes;
        Py_ssize_t sizes_list_o_length = PyList_GET_SIZE(sizes_list_o);
        for (Py_ssize_t i = 0; i < sizes_list_o_length; i++) {
            PyObject * size_o = PyList_GET_ITEM(sizes_list_o, i);
            if (PyLong_Check(size_o)) {
                sizes.push_back((HashIntoType) PyLong_AsUnsignedLongLong(size_o));
            } else if (PyInt_Check(size_o)) {
                sizes.push_back((HashIntoType) PyInt_AsLong(size_o));
            } else if (PyFloat_Check(size_o)) {
                sizes.push_back((HashIntoType) PyFloat_AS_DOUBLE(size_o));
            } else {
                Py_DECREF(self);
                PyErr_SetString(PyExc_TypeError,
                                "2nd argument must be a list of ints, longs, or floats");
                return NULL;
            }
        }


        // We want the hashbits pointer in the base class to point to our labelhash,
        // so that the KHashbits methods are called on the correct object (a LabelHash)
        try {
            self->labelhash = new LabelHash(k, sizes);
        } catch (std::bad_alloc &e) {
            Py_DECREF(self);
            return PyErr_NoMemory();
        }
        self->khashbits.hashbits = (Hashbits *)self->labelhash;
    }

    return (PyObject *) self;
}

static int khmer_labelhash_init(khmer_KLabelHash_Object * self, PyObject *args,
                                PyObject *kwds)
{
    if (khmer_KHashbits_Type.tp_init((PyObject *)self, args, kwds) < 0) {
        return -1;
    }
    //std::cout << "testing my pointer ref to hashbits: " << self->khashbits.hashbits->n_tags() << std::endl;
    //std::cout << "hashbits: " << self->khashbits.hashbits << std::endl;
    //std::cout << "labelhash: " << self->labelhash << std::endl;
    return 0;
}

static
PyObject *
labelhash_get_label_dict(khmer_KLabelHash_Object * me, PyObject * args)
{
    LabelHash * hb = me->labelhash;

    PyObject * d = PyDict_New();
    if (d == NULL) {
        return NULL;
    }
    LabelPtrMap::iterator it;

    for (it = hb->label_ptrs.begin(); it != hb->label_ptrs.end(); ++it) {
        PyObject * key = Py_BuildValue("K", it->first);
        PyObject * val = Py_BuildValue("K", it->second);
        if (key != NULL && val != NULL) {
            PyDict_SetItem(d, key, val);
        }
        Py_XDECREF(key);
        Py_XDECREF(val);
    }

    return d;
}

static
PyObject *
labelhash_consume_fasta_and_tag_with_labels(khmer_KLabelHash_Object * me,
        PyObject * args)
{
    LabelHash * hb = me->labelhash;

    std::ofstream outfile;

    const char * filename;

    if (!PyArg_ParseTuple(args, "s", &filename)) {
        return NULL;
    }

    unsigned long long n_consumed;
    unsigned int total_reads;
    char const * exc = NULL;
    //Py_BEGIN_ALLOW_THREADS
    try {
        hb->consume_fasta_and_tag_with_labels(filename, total_reads,
                                              n_consumed);
    } catch (_khmer_signal &e) {
        exc = e.get_message().c_str();
    } catch (khmer_file_exception &e) {
        exc = e.what();
    }
    //Py_END_ALLOW_THREADS
    if (exc != NULL) {
        PyErr_SetString(PyExc_IOError, exc);
        return NULL;
    }

    return Py_BuildValue("IK", total_reads, n_consumed);

}

static
PyObject *
labelhash_consume_partitioned_fasta_and_tag_with_labels(
    khmer_KLabelHash_Object * me, PyObject * args)
{
    LabelHash * labelhash = me->labelhash;

    const char * filename;

    if (!PyArg_ParseTuple(args, "s", &filename)) {
        return NULL;
    }

    // call the C++ function, and trap signals => Python

    unsigned long long n_consumed;
    unsigned int total_reads;

    try {
        labelhash->consume_partitioned_fasta_and_tag_with_labels(filename,
                total_reads, n_consumed);
    } catch (_khmer_signal &e) {
        PyErr_SetString(PyExc_IOError,
                        "error parsing in consume_partitioned_fasta_and_tag_with_labels");
        return NULL;
    } catch (khmer_file_exception &e) {
        PyErr_SetString(PyExc_IOError, e.what());
        return NULL;
    }
    return Py_BuildValue("IK", total_reads, n_consumed);
}

static
PyObject *
labelhash_consume_sequence_and_tag_with_labels(khmer_KLabelHash_Object * me,
        PyObject * args)
{
    LabelHash * hb = me->labelhash;
    const char * seq = NULL;
    unsigned long long c = 0;
    if (!PyArg_ParseTuple(args, "sK", &seq, &c)) {
        return NULL;
    }
    unsigned long long n_consumed = 0;
    Label * the_label = hb->check_and_allocate_label(c);

    try {
        hb->consume_sequence_and_tag_with_labels(seq, n_consumed, *the_label);
    } catch (_khmer_signal &e) {
        return NULL;
    }
    return Py_BuildValue("K", n_consumed);
}

static
PyObject *
labelhash_sweep_label_neighborhood(khmer_KLabelHash_Object * me,
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

    unsigned int range = (2 * hb->_get_tag_density()) + 1;
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

    if (strlen(seq) < hb->ksize()) {
        PyErr_SetString(PyExc_ValueError,
                        "string length must >= the hashtable k-mer size");
        return NULL;
    }

    //std::pair<TagLabelPtrPair::iterator, TagLabelPtrPair::iterator> ret;
    LabelPtrSet found_labels;

    bool exc_raised = false;
    //unsigned int num_traversed = 0;
    //Py_BEGIN_ALLOW_THREADS
    try {
        hb->sweep_label_neighborhood(seq, found_labels, range, break_on_stop_tags,
                                     stop_big_traversals);
    } catch (_khmer_signal &e) {
        exc_raised = true;
    }
    //Py_END_ALLOW_THREADS

    //printf("...%u kmers traversed\n", num_traversed);

    if (exc_raised) {
        return NULL;
    }

    PyObject * x =  PyList_New(found_labels.size());
    LabelPtrSet::const_iterator si;
    unsigned long long i = 0;
    for (si = found_labels.begin(); si != found_labels.end(); ++si) {
        PyList_SET_ITEM(x, i, Py_BuildValue("K", *(*si)));
        i++;
    }

    return x;
}

// Similar to find_all_tags, but returns tags in a way actually usable by python
// need a tags_in_sequence iterator or function in c++ land for reuse in all
// these functions

static
PyObject *
labelhash_sweep_tag_neighborhood(khmer_KLabelHash_Object * me, PyObject * args)
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

    unsigned int range = (2 * labelhash->_get_tag_density()) + 1;
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

    if (strlen(seq) < labelhash->ksize()) {
        PyErr_SetString(PyExc_ValueError,
                        "string length must >= the hashtable k-mer size");
        return NULL;
    }

    SeenSet tagged_kmers;

    //Py_BEGIN_ALLOW_THREADS

    labelhash->partition->sweep_for_tags(seq, tagged_kmers,
                                         labelhash->all_tags, range, break_on_stop_tags, stop_big_traversals);

    //Py_END_ALLOW_THREADS

    PyObject * x =  PyList_New(tagged_kmers.size());
    if (x == NULL) {
        return NULL;
    }
    SeenSet::iterator si;
    unsigned long long i = 0;
    for (si = tagged_kmers.begin(); si != tagged_kmers.end(); ++si) {
        //std::string kmer_s = _revhash(*si, labelhash->ksize());
        // type K for python unsigned long long
        PyList_SET_ITEM(x, i, Py_BuildValue("K", *si));
        i++;
    }

    return x;
}

static
PyObject *
labelhash_get_tag_labels(khmer_KLabelHash_Object * me, PyObject * args)
{
    LabelHash * labelhash = me->labelhash;

    HashIntoType tag;

    if (!PyArg_ParseTuple(args, "K", &tag)) {
        return NULL;
    }

    LabelPtrSet labels;

    labels = labelhash->get_tag_labels(tag);

    PyObject * x =  PyList_New(labels.size());
    LabelPtrSet::const_iterator si;
    unsigned long long i = 0;
    for (si = labels.begin(); si != labels.end(); ++si) {
        //std::string kmer_s = _revhash(*si, labelhash->ksize());
        PyList_SET_ITEM(x, i, Py_BuildValue("K", *(*si)));
        i++;
    }

    return x;
}

static
PyObject *
labelhash_n_labels(khmer_KLabelHash_Object * me, PyObject * args)
{
    LabelHash * labelhash = me->labelhash;

    if (!PyArg_ParseTuple(args, "")) {
        return NULL;
    }

    return PyLong_FromSize_t(labelhash->n_labels());
}

static PyMethodDef khmer_labelhash_methods[] = {
    { "consume_fasta_and_tag_with_labels", (PyCFunction)labelhash_consume_fasta_and_tag_with_labels, METH_VARARGS, "" },
    { "sweep_label_neighborhood", (PyCFunction)labelhash_sweep_label_neighborhood, METH_VARARGS, "" },
    {"consume_partitioned_fasta_and_tag_with_labels", (PyCFunction)labelhash_consume_partitioned_fasta_and_tag_with_labels, METH_VARARGS, "" },
    {"sweep_tag_neighborhood", (PyCFunction)labelhash_sweep_tag_neighborhood, METH_VARARGS, "" },
    {"get_tag_labels", (PyCFunction)labelhash_get_tag_labels, METH_VARARGS, ""},
    {"consume_sequence_and_tag_with_labels", (PyCFunction)labelhash_consume_sequence_and_tag_with_labels, METH_VARARGS, "" },
    {"n_labels", (PyCFunction)labelhash_n_labels, METH_VARARGS, ""},
    {"get_label_dict", (PyCFunction)labelhash_get_label_dict, METH_VARARGS, "" },
    {NULL, NULL, 0, NULL}           /* sentinel */
};

static PyTypeObject khmer_KLabelHash_Type = {
    PyVarObject_HEAD_INIT(NULL, 0)  /* init & ob_size */
    "_khmer.LabelHash",            /* tp_name */
    sizeof(khmer_KLabelHash_Object), /* tp_basicsize */
    0,                       /* tp_itemsize */
    (destructor)khmer_labelhash_dealloc, /* tp_dealloc */
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
    khmer_labelhash_methods, /* tp_methods */
    0,                       /* tp_members */
    0,                       /* tp_getset */
    0,                       /* tp_base */
    0,                       /* tp_dict */
    0,                       /* tp_descr_get */
    0,                       /* tp_descr_set */
    0,                       /* tp_dictoffset */
    (initproc)khmer_labelhash_init,   /* tp_init */
    0,                       /* tp_alloc */
    khmer_labelhash_new,      /* tp_new */
};

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

static PyMethodDef khmer_ReadAligner_methods[] = {
    {"align", (PyCFunction)readaligner_align, METH_VARARGS, ""},
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
        khmer_KCountingHash_Object * ch = NULL;
        unsigned short int trusted_cov_cutoff = 2;
        double bits_theta = 1;

        if(!PyArg_ParseTuple(args, "O!Hd", &khmer_KCountingHash_Type, &ch,
                             &trusted_cov_cutoff, &bits_theta)) {
            Py_DECREF(self);
            return NULL;
        }

        self->aligner = new ReadAligner(ch->counting, trusted_cov_cutoff,
                                        bits_theta);
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
    Py_TPFLAGS_DEFAULT,         /*tp_flags*/
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

static PyObject * hash_collect_high_abundance_kmers(khmer_KCountingHash_Object *
        me , PyObject * args)
{
    CountingHash * counting = me->counting;

    const char * filename = NULL;
    unsigned int lower_count, upper_count;

    if (!PyArg_ParseTuple(args, "sII", &filename, &lower_count, &upper_count)) {
        return NULL;
    }

    SeenSet found_kmers;
    counting->collect_high_abundance_kmers(filename, lower_count, upper_count,
                                           found_kmers);

    // create a new hashbits object...
    std::vector<HashIntoType> sizes;
    sizes.push_back(1);

    khmer_KHashbits_Object * khashbits_obj = (khmer_KHashbits_Object *) \
            PyObject_New(khmer_KHashbits_Object, &khmer_KHashbits_Type);
    if (khashbits_obj == NULL) {
        return NULL;
    }

    // ...and set the collected kmers as the stoptags.
    try {
        khashbits_obj->hashbits = new Hashbits(counting->ksize(), sizes);
    } catch (std::bad_alloc &e) {
        return PyErr_NoMemory();
    }
    khashbits_obj->hashbits->stop_tags.swap(found_kmers);

    return (PyObject *) khashbits_obj;
}

//
// khmer_counting_dealloc -- clean up a counting hash object.
//

static void khmer_counting_dealloc(khmer_KCountingHash_Object * obj)
{
    delete obj->counting;
    obj->counting = NULL;
    Py_TYPE(obj)->tp_free((PyObject*)obj);
}

//
// khmer_hashbits_dealloc -- clean up a hashbits object.
//
static void khmer_hashbits_dealloc(khmer_KHashbits_Object * obj)
{
    delete obj->hashbits;
    obj->hashbits = NULL;

    Py_TYPE(obj)->tp_free((PyObject*)obj);
}


//
// khmer_subset_dealloc -- clean up a hashbits object.
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
    khmer::HLLCounter * hllcounter;
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

static PyObject * hllcounter_consume_fasta(khmer_KHLLCounter_Object * me,
        PyObject * args)
{
    const char * filename;

    if (!PyArg_ParseTuple(args, "s", &filename)) {
        return NULL;
    }

    // call the C++ function, and trap signals => Python
    unsigned long long  n_consumed    = 0;
    unsigned int        total_reads   = 0;
    try {
        me->hllcounter->consume_fasta(filename, total_reads, n_consumed);
    } catch (_khmer_signal &e) {
        PyErr_SetString(PyExc_IOError, e.get_message().c_str());
        return NULL;
    } catch (khmer_file_exception &e) {
        PyErr_SetString(PyExc_IOError, e.what());
        return NULL;
    }

    return Py_BuildValue("IK", total_reads, n_consumed);
}

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
        "consume_fasta", (PyCFunction)hllcounter_consume_fasta,
        METH_VARARGS,
        "Read sequences from file, break into k-mers, "
        "and add each k-mer to the counter."
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

#define is_hllcounter_obj(v)  (Py_TYPE(v) == &khmer_KHLLCounter_Type)


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

    return PyLong_FromUnsignedLongLong(_hash(kmer, ksize));
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

    return PyLong_FromUnsignedLongLong(_hash_forward(kmer, ksize));
}

static PyObject * reverse_hash(PyObject * self, PyObject * args)
{
    HashIntoType val;
    WordLength ksize;

    if (!PyArg_ParseTuple(args, "Kb", &val, &ksize)) {
        return NULL;
    }

    if (ksize > KSIZE_MAX) {
        PyErr_Format(PyExc_ValueError, "k-mer size must be <= %u", KSIZE_MAX);
        return NULL;
    }

    return PyBytes_FromString(_revhash(val, ksize).c_str());
}

static PyObject * murmur3_forward_hash(PyObject * self, PyObject * args)
{
    const char * kmer;

    if (!PyArg_ParseTuple(args, "s", &kmer)) {
        return NULL;
    }

    return PyLong_FromUnsignedLongLong(_hash_murmur(kmer));
}

static PyObject * murmur3_forward_hash_no_rc(PyObject * self, PyObject * args)
{
    const char * kmer;

    if (!PyArg_ParseTuple(args, "s", &kmer)) {
        return NULL;
    }

    return PyLong_FromUnsignedLongLong(_hash_murmur_forward(kmer));
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
    return PyBytes_FromString(dVersion.c_str());
}


//
// Module machinery.
//

static PyMethodDef KhmerMethods[] = {
#if (0)
    {
        "new_config",       new_config,
        METH_VARARGS,       "Create a default internals config"
    },
#endif
#if (0)
    {
        "set_config",       set_active_config,
        METH_VARARGS,       "Set active khmer configuration object"
    },
#endif
    {
        "new_hashtable",        new_hashtable,
        METH_VARARGS,       "Create an empty single-table counting hash"
    },
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
        "get_version_cpp", get_version_cpp,
        METH_VARARGS, "return the VERSION c++ compiler option"
    },
    { NULL, NULL, 0, NULL } // sentinel
};

PyMODINIT_FUNC
init_khmer(void)
{
    using namespace python;

    if (PyType_Ready(&khmer_KCountingHash_Type) < 0) {
        return;
    }

    khmer_KSubsetPartition_Type.tp_methods = khmer_subset_methods;
    if (PyType_Ready(&khmer_KSubsetPartition_Type) < 0) {
        return;
    }

    khmer_KHashbits_Type.tp_methods = khmer_hashbits_methods;
    if (PyType_Ready(&khmer_KHashbits_Type) < 0) {
        return;
    }
    // add LabelHash

    khmer_KLabelHash_Type.tp_base = &khmer_KHashbits_Type;
    if (PyType_Ready(&khmer_KLabelHash_Type) < 0) {
        return;
    }

    if (PyType_Ready(&khmer_ReadAlignerType) < 0) {
        return;
    }

    if (PyType_Ready(&khmer_KHLLCounter_Type) < 0) {
        return;
    }
    if (PyType_Ready(&khmer_ReadAlignerType) < 0) {
        return;
    }

    _init_ReadParser_Type_constants();
    if (PyType_Ready( &khmer_ReadParser_Type ) < 0) {
        return;
    }

    if (PyType_Ready(&khmer_Read_Type ) < 0) {
        return;
    }

    if (PyType_Ready(&khmer_ReadPairIterator_Type ) < 0) {
        return;
    }

    PyObject * m;
    m = Py_InitModule3( "_khmer", KhmerMethods,
                        "interface for the khmer module low-level extensions" );
    if (m == NULL) {
        return;
    }

    Py_INCREF(&khmer_ReadParser_Type);
    if (PyModule_AddObject( m, "ReadParser",
                            (PyObject *)&khmer_ReadParser_Type ) < 0) {
        return;
    }

    Py_INCREF(&khmer_KCountingHash_Type);
    if (PyModule_AddObject( m, "CountingHash",
                            (PyObject *)&khmer_KCountingHash_Type ) < 0) {
        return;
    }

    Py_INCREF(&khmer_KHashbits_Type);
    if (PyModule_AddObject(m, "Hashbits", (PyObject *)&khmer_KHashbits_Type) < 0) {
        return;
    }

    Py_INCREF(&khmer_KLabelHash_Type);
    if (PyModule_AddObject(m, "LabelHash",
                           (PyObject *)&khmer_KLabelHash_Type) < 0) {
        return;
    }

    Py_INCREF(&khmer_KHLLCounter_Type);
    PyModule_AddObject(m, "HLLCounter", (PyObject *)&khmer_KHLLCounter_Type);
    Py_INCREF(&khmer_ReadAlignerType);
    PyModule_AddObject(m, "ReadAligner", (PyObject *)&khmer_ReadAlignerType);
}

// vim: set ft=cpp sts=4 sw=4 tw=79:
