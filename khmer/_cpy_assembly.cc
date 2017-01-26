#include "_cpy_assembly.hh"
#include "_cpy_nodegraph.hh"
#include "_cpy_countgraph.hh"
#include "_cpy_graphlabels.hh"

namespace khmer {

PyTypeObject khmer_KSimpleLabeledAssembler_Type = {
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

PyMethodDef khmer_simplelabeledassembler_methods[] = {
    {
        "assemble",
        (PyCFunction)simplelabeledassembler_assemble, METH_VARARGS | METH_KEYWORDS,
        "Assemble paths, using labels to jump branches."
    },
    {NULL, NULL, 0, NULL}           /* sentinel */
};


void khmer_simplelabeledassembler_dealloc(khmer_KSimpleLabeledAssembler_Object *
        obj)
{
    delete obj->assembler;
    obj->assembler = NULL;

    Py_TYPE(obj)->tp_free((PyObject*)obj);
}


PyObject * khmer_simplelabeledassembler_new(PyTypeObject *type,
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



/********************************
 * JunctionCountAssembler
 ********************************/


PyTypeObject khmer_KJunctionCountAssembler_Type = {
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

PyMethodDef khmer_junctioncountassembler_methods[] = {
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


void khmer_junctioncountassembler_dealloc(
    khmer_KJunctionCountAssembler_Object * obj)
{
    delete obj->assembler;
    obj->assembler = NULL;

    Py_TYPE(obj)->tp_free((PyObject*)obj);
}

PyObject * khmer_junctioncountassembler_new(PyTypeObject *type,
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

}
