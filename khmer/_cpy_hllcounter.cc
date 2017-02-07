#include "_cpy_hllcounter.hh"

namespace khmer {


PyTypeObject khmer_KHLLCounter_Type = {
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


PyMethodDef khmer_hllcounter_methods[] = {
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

PyGetSetDef khmer_hllcounter_getseters[] = {
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


PyObject* khmer_hllcounter_new(PyTypeObject * type, PyObject * args,
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


void khmer_hllcounter_dealloc(khmer_KHLLCounter_Object * obj)
{
    delete obj->hllcounter;
    obj->hllcounter = NULL;

    Py_TYPE(obj)->tp_free((PyObject*)obj);
}


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


PyObject *
hllcounter_estimate_cardinality(khmer_KHLLCounter_Object * me, PyObject * args)
{
    if (!PyArg_ParseTuple( args, "" )) {
        return NULL;
    }

    return PyLong_FromLong(me->hllcounter->estimate_cardinality());
}


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

PyObject * hllcounter_consume_fasta(khmer_KHLLCounter_Object * me,
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
        me->hllcounter->consume_fasta(filename, stream_records, total_reads,
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

 PyObject * hllcounter_merge(khmer_KHLLCounter_Object * me,
                                   PyObject * args);


PyObject *
hllcounter_get_erate(khmer_KHLLCounter_Object * me)
{
    return PyFloat_FromDouble(me->hllcounter->get_erate());
}


PyObject *
hllcounter_get_ksize(khmer_KHLLCounter_Object * me)
{
    return PyLong_FromLong(me->hllcounter->get_ksize());
}


int hllcounter_set_ksize(khmer_KHLLCounter_Object * me, PyObject *value,
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


int hllcounter_set_erate(khmer_KHLLCounter_Object * me, PyObject *value,
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


PyObject *
hllcounter_getalpha(khmer_KHLLCounter_Object * me)
{
    return PyFloat_FromDouble(me->hllcounter->get_alpha());
}


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


PyObject * hllcounter_merge(khmer_KHLLCounter_Object * me,
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

}
