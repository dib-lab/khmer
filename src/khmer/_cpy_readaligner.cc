#include "khmer/_cpy_countgraph.hh"
#include "khmer/_cpy_readaligner.hh"
#include <string>

using namespace oxli;

namespace khmer {




PyObject * readaligner_align_forward(khmer_ReadAligner_Object * me,
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


