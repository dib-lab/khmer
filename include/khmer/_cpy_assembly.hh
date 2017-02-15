#ifndef _CPY_ASSEMBLY_HH
#define _CPY_ASSEMBLY_HH

#include <Python.h>
#include "_cpy_utils.hh"
#include "oxli/assembler.hh"


namespace khmer {

typedef struct {
    PyObject_HEAD
    oxli::SimpleLabeledAssembler * assembler;
} khmer_KSimpleLabeledAssembler_Object;

extern PyTypeObject khmer_KSimpleLabeledAssembler_Type;

extern PyMethodDef khmer_simplelabeledassembler_methods[];

void khmer_simplelabeledassembler_dealloc(khmer_KSimpleLabeledAssembler_Object *
        obj);

PyObject * khmer_simplelabeledassembler_new(PyTypeObject *type,
        PyObject *args,
        PyObject *kwds);

PyObject *
simplelabeledassembler_assemble(khmer_KSimpleLabeledAssembler_Object * me,
                                PyObject * args, PyObject *kwargs);


/********************************
 * JunctionCountAssembler
 ********************************/


typedef struct {
    PyObject_HEAD
    oxli::JunctionCountAssembler * assembler;
} khmer_KJunctionCountAssembler_Object;


extern PyTypeObject khmer_KJunctionCountAssembler_Type;

extern PyMethodDef khmer_junctioncountassembler_methods[];


void khmer_junctioncountassembler_dealloc(khmer_KJunctionCountAssembler_Object * obj);

PyObject * khmer_junctioncountassembler_new(PyTypeObject *type,
        PyObject *args,
        PyObject *kwds);


PyObject *
junctioncountassembler_assemble(khmer_KJunctionCountAssembler_Object * me,
                                PyObject * args, PyObject *kwargs);


PyObject *
junctioncountassembler_consume(khmer_KJunctionCountAssembler_Object * me,
                               PyObject * args);

}

#endif
