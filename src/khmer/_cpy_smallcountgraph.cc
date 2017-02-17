#include "khmer/_cpy_smallcountgraph.hh"
#include "oxli/hashgraph.hh"

using namespace oxli;


namespace khmer {

PyTypeObject khmer_KSmallCountgraph_Type
CPYCHECKER_TYPE_OBJECT_FOR_TYPEDEF("khmer_KSmallCountgraph_Object")
= {
    PyVarObject_HEAD_INIT(NULL, 0)       /* init & ob_size */
    "_khmer.SmallCountgraph",                 /*tp_name*/
    sizeof(khmer_KSmallCountgraph_Object),  /*tp_basicsize*/
    0,                                   /*tp_itemsize*/
    (destructor)khmer_smallcountgraph_dealloc,  /*tp_dealloc*/
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
    "SmallCountgraph hash object",              /* tp_doc */
    0,                                   /* tp_traverse */
    0,                                   /* tp_clear */
    0,                                   /* tp_richcompare */
    0,                                   /* tp_weaklistoffset */
    0,                                   /* tp_iter */
    0,                                   /* tp_iternext */
    khmer_smallcountgraph_methods,              /* tp_methods */
    0,                                   /* tp_members */
    0,                                   /* tp_getset */
    0,                                   /* tp_base */
    0,                                   /* tp_dict */
    0,                                   /* tp_descr_get */
    0,                                   /* tp_descr_set */
    0,                                   /* tp_dictoffset */
    0,                                   /* tp_init */
    0,                                   /* tp_alloc */
    khmer_smallcountgraph_new,                /* tp_new */
};

PyMethodDef khmer_smallcountgraph_methods[] = {
    {
        "get_raw_tables",
        (PyCFunction)smallcount_get_raw_tables, METH_VARARGS,
        "Get a list of the raw storage tables as memoryview objects."
    },
    {NULL, NULL, 0, NULL}           /* sentinel */
};

PyObject *
smallcount_get_raw_tables(khmer_KSmallCountgraph_Object * self, PyObject * args)
{
    SmallCountgraph * countgraph = self->countgraph;

    oxli::Byte ** table_ptrs = countgraph->get_raw_tables();
    std::vector<uint64_t> sizes = countgraph->get_tablesizes();

    PyObject * raw_tables = PyList_New(sizes.size());
    for (unsigned int i=0; i<sizes.size(); ++i) {
        Py_buffer buffer;
        int res = PyBuffer_FillInfo(&buffer, NULL, table_ptrs[i],
                                    sizes[i] / 2 +1, 0,
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


//
// khmer_smallcountgraph_dealloc -- clean up a countgraph hash object.
//

void khmer_smallcountgraph_dealloc(khmer_KSmallCountgraph_Object * obj)
{
    delete obj->countgraph;
    obj->countgraph = NULL;
    Py_TYPE(obj)->tp_free((PyObject*)obj);
}

//
// khmer_smallcountgraph_new
//

PyObject* khmer_smallcountgraph_new(PyTypeObject * type, PyObject * args,
        PyObject * kwds)
{
    khmer_KSmallCountgraph_Object * self;

    self = (khmer_KSmallCountgraph_Object *)type->tp_alloc(type, 0);

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
            self->countgraph = new SmallCountgraph(k, sizes);
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

}
