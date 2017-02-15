#include "_cpy_nodegraph.hh"
#include "hashgraph.hh"

namespace khmer {


PyTypeObject khmer_KNodegraph_Type
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

PyMethodDef khmer_nodegraph_methods[] = {
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

void khmer_nodegraph_dealloc(khmer_KNodegraph_Object * obj)
{
    delete obj->nodegraph;
    obj->nodegraph = NULL;

    Py_TYPE(obj)->tp_free((PyObject*)obj);
}


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


// __new__ for nodegraph; necessary for proper subclassing
// This will essentially do what the old factory function did. Unlike many __new__
// methods, we take our arguments here, because there's no "uninitialized" nodegraph
// object; we have to have k and the table sizes before creating the new objects
PyObject* khmer_nodegraph_new(PyTypeObject * type, PyObject * args,
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

}
