
#ifndef _CPY_READPARSERS_HH
#define _CPY_READPARSERS_HH

#include <Python.h>
#include "_cpy_utils.hh"
#include "oxli/read_parsers.hh"

namespace khmer
{

typedef struct {
    PyObject_HEAD
    //! Pointer to the low-level genomic read object.
    oxli::read_parsers::Read *   read;
} khmer_Read_Object;


typedef struct {
    PyObject_HEAD
    oxli::read_parsers::FastxParserPtr parser;
} khmer_ReadParser_Object;


typedef struct {
    PyObject_HEAD
    //! Pointer to Python parser object for reference counting purposes.
    PyObject * parent;
    //! Persistent value of pair mode across invocations.
    int pair_mode;
} khmer_ReadPairIterator_Object;



PyObject* 
khmer_Read_new(PyTypeObject * type, PyObject * args, PyObject * kwds);


int
khmer_Read_init(khmer_Read_Object *self, PyObject *args, PyObject *kwds);


void
khmer_Read_dealloc(khmer_Read_Object * obj);


Py_ssize_t
khmer_Read_len(khmer_Read_Object* obj);


extern PySequenceMethods khmer_Read_sequence_methods;


PyObject *
Read_get_name(khmer_Read_Object * obj, void * closure );


PyObject *
Read_get_sequence(khmer_Read_Object * obj, void * closure);


PyObject *
Read_get_quality(khmer_Read_Object * obj, void * closure);


PyObject *
Read_get_description(khmer_Read_Object * obj, void * closure);


PyObject *
Read_get_cleaned_seq(khmer_Read_Object * obj, void * closure);


int
Read_set_cleaned_seq(khmer_Read_Object *obj, PyObject *value, void *closure);


extern PyGetSetDef khmer_Read_accessors [];


extern PyTypeObject khmer_Read_Type;


void
_ReadParser_dealloc(khmer_ReadParser_Object * obj);


void
khmer_ReadPairIterator_dealloc(khmer_ReadPairIterator_Object * obj);


PyObject *
_ReadParser_new( PyTypeObject * subtype, PyObject * args, PyObject * kwds );


PyObject *
_ReadParser_iternext( PyObject * self );


PyObject *
_ReadPairIterator_iternext(khmer_ReadPairIterator_Object * myself);


extern PyTypeObject khmer_ReadPairIterator_Type;


PyObject *
ReadParser_iter_reads(PyObject * self, PyObject * args );


PyObject *
ReadParser_get_num_reads(khmer_ReadParser_Object * me);


PyObject *
ReadParser_iter_read_pairs(PyObject * self, PyObject * args );


PyObject *
ReadParser_close(PyObject * self, PyObject * args);


extern PyMethodDef _ReadParser_methods [];


extern PyGetSetDef khmer_ReadParser_accessors[];


extern PyTypeObject khmer_ReadParser_Type
CPYCHECKER_TYPE_OBJECT_FOR_TYPEDEF("khmer_ReadParser_Object");


void _init_ReadParser_Type_constants();


oxli::read_parsers::FastxParserPtr& _PyObject_to_khmer_ReadParser(PyObject * py_object);

} // namespace khmer

#endif
