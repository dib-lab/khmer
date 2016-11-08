typedef struct {
    khmer_KHashtable_Object khashtable;
    Hashgraph * hashgraph;
} khmer_KHashgraph_Object;

static PyTypeObject khmer_KHashgraph_Type
CPYCHECKER_TYPE_OBJECT_FOR_TYPEDEF("khmer_KHashgraph_Object")
= {
    PyVarObject_HEAD_INIT(NULL, 0) /* init & ob_size */
    "_khmer.Hashgraph",             /* tp_name */
    sizeof(khmer_KHashgraph_Object), /* tp_basicsize */
    0,                             /* tp_itemsize */
    0,                             /*tp_dealloc*/
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
    "hashgraph object",           /* tp_doc */
    0,                       /* tp_traverse */
    0,                       /* tp_clear */
    0,                       /* tp_richcompare */
    0,                       /* tp_weaklistoffset */
    0,                       /* tp_iter */
    0,                       /* tp_iternext */
    0,                       /* tp_methods */
    0,                       /* tp_members */
    0,                       /* tp_getset */
    0,                       /* tp_base */
    0,                       /* tp_dict */
    0,                       /* tp_descr_get */
    0,                       /* tp_descr_set */
    0,                       /* tp_dictoffset */
    0,                       /* tp_init */
    0,                       /* tp_alloc */
    0,                       /* tp_new */
};


static PyMethodDef khmer_hashgraph_methods[] = {
    {NULL, NULL, 0, NULL}           /* sentinel */
};

#define is_hashgraph_obj(v)  (Py_TYPE(v) == &khmer_KHashgraph_Type)
