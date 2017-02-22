#include "khmer/_cpy_readparsers.hh"
#include <string>

using namespace oxli;
using namespace oxli::read_parsers;

namespace khmer {

PyTypeObject khmer_Read_Type = {
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

PyGetSetDef khmer_Read_accessors [ ] = {
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


PySequenceMethods khmer_Read_sequence_methods = {
    (lenfunc)khmer_Read_len,                  /* sq_length */
};


PyTypeObject khmer_ReadPairIterator_Type = {
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

PyMethodDef _ReadParser_methods [ ] = {
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


PyGetSetDef khmer_ReadParser_accessors[] = {
    {
        (char *)"num_reads",
        (getter)ReadParser_get_num_reads, NULL,
        (char *)"count of reads processed thus far.",
        NULL
    },
    {NULL, NULL, NULL, NULL, NULL} /* Sentinel */
};


PyTypeObject khmer_ReadParser_Type
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



/***********************************************************************/

//
// Read object -- name, sequence, and FASTQ stuff
//

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

void
khmer_Read_dealloc(khmer_Read_Object * obj)
{
    delete obj->read;
    obj->read = NULL;
    Py_TYPE(obj)->tp_free((PyObject*)obj);
}


Py_ssize_t
khmer_Read_len(khmer_Read_Object* obj)
{
    return obj->read->sequence.size();
}


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


int
Read_set_cleaned_seq(khmer_Read_Object *obj, PyObject *value, void *closure)
{
    if (value == NULL) {
        obj->read->cleaned_seq = "";
        return 0;
    }

    if (! (PyUnicode_Check(value) | PyBytes_Check(value))) {
        PyErr_SetString(PyExc_TypeError,
                        "The 'cleaned_seq' attribute value must be a string");
        return -1;
    }

    if (PyUnicode_Check(value)) {
        PyObject* temp = PyUnicode_AsASCIIString(value);
        if (temp == NULL) {
            PyErr_SetString(PyExc_TypeError,
                            "Could not encode 'cleaned_seq' as ASCII.");
            return -1;
        }

        obj->read->cleaned_seq = std::string(PyBytes_AS_STRING(temp));
        Py_DECREF(temp);
    }
    // in python2 not everything is a unicode object
    else {
        obj->read->cleaned_seq = std::string(PyBytes_AS_STRING(value));
    }

    return 0;
}


// TODO? Implement setters.



/***********************************************************************/

//
// ReadParser object -- parse reads directly from streams
// ReadPairIterator -- return pairs of Read objects
//


void
_ReadParser_dealloc(khmer_ReadParser_Object * obj)
{
    Py_TYPE(obj)->tp_free((PyObject*)obj);
}


void
khmer_ReadPairIterator_dealloc(khmer_ReadPairIterator_Object * obj)
{
    Py_DECREF(obj->parent);
    obj->parent = NULL;
    Py_TYPE(obj)->tp_free((PyObject*)obj);
}


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
        myself->parser = get_parser<FastxReader>(ifile_name);
    } catch (oxli_file_exception &exc) {
        PyErr_SetString( PyExc_OSError, exc.what() );
        return NULL;
    }
    return self;
}

PyObject *
_ReadParser_iternext( PyObject * self )
{
    khmer_ReadParser_Object * myself  = (khmer_ReadParser_Object *)self;
    FastxParserPtr&       parser  = myself->parser;
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
        } catch (oxli_file_exception &exc) {
            exc_string = exc.what();
            file_exception = exc_string.c_str();
        } catch (oxli_value_exception &exc) {
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
        } catch (oxli_file_exception &exc) {
            exc_string = exc.what();
            file_exception = exc_string.c_str();
        } catch (oxli_value_exception &exc) {
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


PyObject *
ReadParser_iter_reads(PyObject * self, PyObject * args )
{
    return PyObject_SelfIter( self );
}

PyObject *
ReadParser_get_num_reads(khmer_ReadParser_Object * me)
{
    return PyLong_FromLong(me->parser->get_num_reads());
}

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


FastxParserPtr& _PyObject_to_khmer_ReadParser(PyObject * py_object)
{
    // TODO: Add type-checking.

    return ((khmer_ReadParser_Object *)py_object)->parser;
}


}
