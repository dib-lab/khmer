typedef struct {
    khmer_KHashbits_Object khashbits;
    Hashbits * nodegraph;
} khmer_KNodegraph_Object;

static void khmer_nodegraph_dealloc(khmer_KNodegraph_Object * obj);
static PyObject* khmer_nodegraph_new(PyTypeObject * type, PyObject * args,
                                     PyObject * kwds);

static PyTypeObject khmer_KNodegraph_Type
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


static PyMethodDef khmer_nodegraph_methods[] = {
    {NULL, NULL, 0, NULL}           /* sentinel */
};

// __new__ for nodegraph; necessary for proper subclassing
// This will essentially do what the old factory function did. Unlike many __new__
// methods, we take our arguments here, because there's no "uninitialized" nodegraph
// object; we have to have k and the table sizes before creating the new objects
static PyObject* khmer_nodegraph_new(PyTypeObject * type, PyObject * args,
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
        Py_ssize_t sizes_list_o_length = PyList_GET_SIZE(sizes_list_o);
        for (Py_ssize_t i = 0; i < sizes_list_o_length; i++) {
            PyObject * size_o = PyList_GET_ITEM(sizes_list_o, i);
            if (PyLong_Check(size_o)) {
                sizes.push_back(PyLong_AsUnsignedLongLong(size_o));
            } else if (PyInt_Check(size_o)) {
                sizes.push_back(PyInt_AsLong(size_o));
            } else if (PyFloat_Check(size_o)) {
                sizes.push_back(PyFloat_AS_DOUBLE(size_o));
            } else {
                Py_DECREF(self);
                PyErr_SetString(PyExc_TypeError,
                                "2nd argument must be a list of ints, longs, or floats");
                return NULL;
            }
        }

        try {
            self->nodegraph = new Hashbits(k, sizes);
        } catch (std::bad_alloc &e) {
            Py_DECREF(self);
            return PyErr_NoMemory();
        }
        self->khashbits.khashtable.hashtable = self->nodegraph;
        self->khashbits.hashbits = self->nodegraph;
    }
    return (PyObject *) self;
}

//
// khmer_nodegraph_dealloc -- clean up a nodegraph object.
//

static void khmer_nodegraph_dealloc(khmer_KNodegraph_Object * obj)
{
    delete obj->nodegraph;
    obj->nodegraph = NULL;

    Py_TYPE(obj)->tp_free((PyObject*)obj);
}

#define is_nodegraph_obj(v)  (Py_TYPE(v) == &khmer_KNodegraph_Type)
