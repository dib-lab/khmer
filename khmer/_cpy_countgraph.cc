#include "_cpy_countgraph.hh"
#include "_cpy_subsetpartition.hh"
#include "khmer.hh"
#include "hashgraph.hh"

namespace khmer {

PyTypeObject khmer_KCountgraph_Type
CPYCHECKER_TYPE_OBJECT_FOR_TYPEDEF("khmer_KCountgraph_Object")
= {
    PyVarObject_HEAD_INIT(NULL, 0)       /* init & ob_size */
    "_khmer.Countgraph",                 /*tp_name*/
    sizeof(khmer_KCountgraph_Object),  /*tp_basicsize*/
    0,                                   /*tp_itemsize*/
    (destructor)khmer_countgraph_dealloc,  /*tp_dealloc*/
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
    "countgraph hash object",              /* tp_doc */
    0,                                   /* tp_traverse */
    0,                                   /* tp_clear */
    0,                                   /* tp_richcompare */
    0,                                   /* tp_weaklistoffset */
    0,                                   /* tp_iter */
    0,                                   /* tp_iternext */
    khmer_countgraph_methods,              /* tp_methods */
    0,                                   /* tp_members */
    0,                                   /* tp_getset */
    0,                                   /* tp_base */
    0,                                   /* tp_dict */
    0,                                   /* tp_descr_get */
    0,                                   /* tp_descr_set */
    0,                                   /* tp_dictoffset */
    0,                                   /* tp_init */
    0,                                   /* tp_alloc */
    khmer_countgraph_new,                /* tp_new */
};


PyMethodDef khmer_countgraph_methods[] = {
    {
        "get_raw_tables",
        (PyCFunction)count_get_raw_tables, METH_VARARGS,
        "Get a list of the raw storage tables as memoryview objects."
    },
    { "do_subset_partition_with_abundance", (PyCFunction)count_do_subset_partition_with_abundance, METH_VARARGS, "" },
    {NULL, NULL, 0, NULL}           /* sentinel */
};


void khmer_countgraph_dealloc(khmer_KCountgraph_Object * obj)
{
    delete obj->countgraph;
    obj->countgraph = NULL;
    Py_TYPE(obj)->tp_free((PyObject*)obj);
}


PyObject *
count_get_raw_tables(khmer_KCountgraph_Object * self, PyObject * args)
{
    Countgraph * countgraph = self->countgraph;

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


PyObject *
count_do_subset_partition_with_abundance(khmer_KCountgraph_Object * me,
        PyObject * args)
{
    Countgraph * countgraph = me->countgraph;

    HashIntoType start_kmer = 0, end_kmer = 0;
    PyObject * break_on_stop_tags_o = NULL;
    PyObject * stop_big_traversals_o = NULL;
    BoundedCounterType min_count, max_count;

    if (!PyArg_ParseTuple(args, "HH|KKOO",
                          &min_count, &max_count,
                          &start_kmer, &end_kmer,
                          &break_on_stop_tags_o,
                          &stop_big_traversals_o)) {
        return NULL;
    }

    bool break_on_stop_tags = false;
    if (break_on_stop_tags_o && PyObject_IsTrue(break_on_stop_tags_o)) {
        break_on_stop_tags = true;
    }
    bool stop_big_traversals = false;
    if (stop_big_traversals_o && PyObject_IsTrue(stop_big_traversals_o)) {
        stop_big_traversals = true;
    }

    SubsetPartition * subset_p = NULL;
    try {
        Py_BEGIN_ALLOW_THREADS
        subset_p = new SubsetPartition(countgraph);
        subset_p->do_partition_with_abundance(start_kmer, end_kmer,
                                              min_count, max_count,
                                              break_on_stop_tags,
                                              stop_big_traversals);
        Py_END_ALLOW_THREADS
    } catch (std::bad_alloc &e) {
        return PyErr_NoMemory();
    }

    khmer_KSubsetPartition_Object * subset_obj = (khmer_KSubsetPartition_Object *)\
            PyObject_New(khmer_KSubsetPartition_Object, &khmer_KSubsetPartition_Type);

    if (subset_obj == NULL) {
        delete subset_p;
        return NULL;
    }

    subset_obj->subset = subset_p;

    return (PyObject *) subset_obj;
}


//
// khmer_countgraph_new
//

 PyObject* khmer_countgraph_new(PyTypeObject * type, PyObject * args,
                                      PyObject * kwds)
{
    khmer_KCountgraph_Object * self;

    self = (khmer_KCountgraph_Object *)type->tp_alloc(type, 0);

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
            self->countgraph = new Countgraph(k, sizes);
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
