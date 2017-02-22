#include "khmer/_cpy_countgraph.hh"
#include "khmer/_cpy_readaligner.hh"
#include <string>

using namespace oxli;

namespace khmer {

PyTypeObject khmer_ReadAlignerType = {
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



PyMethodDef khmer_ReadAligner_methods[] = {
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

PyObject * readaligner_align(khmer_ReadAligner_Object * me,
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

PyObject* khmer_ReadAligner_get_scoring_matrix(
    khmer_ReadAligner_Object * me, PyObject * args)
{

    if (!PyArg_ParseTuple(args, "")) {
        return NULL;
    }
    ScoringMatrix matrix = me->aligner->getScoringMatrix();

    return Py_BuildValue( "dddd", matrix.trusted_match, matrix.trusted_mismatch,
                          matrix.untrusted_match, matrix.untrusted_mismatch);
}

PyObject* khmer_ReadAligner_get_transition_probabilities(
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



//
// khmer_readaligner_dealloc -- clean up readaligner object
// GRAPHALIGN addition
//
void khmer_readaligner_dealloc(khmer_ReadAligner_Object* obj)
{
    delete obj->aligner;
    obj->aligner = NULL;
    Py_TYPE(obj)->tp_free((PyObject*)obj);
}

//
// new_readaligner
//
PyObject* khmer_ReadAligner_new(PyTypeObject *type, PyObject * args,
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

}
