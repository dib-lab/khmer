#include "khmer/_cpy_subsetpartition.hh"
#include "khmer/_cpy_countgraph.hh"

using namespace oxli;

namespace khmer {


PyTypeObject khmer_PrePartitionInfo_Type = {
    PyVarObject_HEAD_INIT(NULL, 0)        /* init & ob_size */
    "_khmer.PrePartitionInfo",            /* tp_name */
    sizeof(khmer_PrePartitionInfo_Object),/* tp_basicsize */
    0,                                    /* tp_itemsize */
    (destructor)khmer_PrePartitionInfo_dealloc,       /* tp_dealloc */
    0,                                    /* tp_print */
    0,                                    /* tp_getattr */
    0,                                    /* tp_setattr */
    0,                                    /* tp_compare */
    0,                                    /* tp_repr */
    0,                                    /* tp_as_number */
    0,                                    /* tp_as_sequence */
    0,                                    /* tp_as_mapping */
    0,                                    /* tp_hash */
    0,                                    /* tp_call */
    0,                                    /* tp_str */
    0,                                    /* tp_getattro */
    0,                                    /* tp_setattro */
    0,                                    /* tp_as_buffer */
    Py_TPFLAGS_DEFAULT,                   /* tp_flags */
    "Stores a k-kmer and a set of tagged seen k-mers.", /* tp_doc */
};

void
khmer_PrePartitionInfo_dealloc(khmer_PrePartitionInfo_Object * obj)
{
    delete obj->PrePartitionInfo;
    obj->PrePartitionInfo = NULL;
    Py_TYPE(obj)->tp_free((PyObject*)obj);
}



PyTypeObject khmer_KSubsetPartition_Type = {
    PyVarObject_HEAD_INIT(NULL, 0)         /* init & ob_size */
    "_khmer.KSubsetPartition",              /* tp_name */
    sizeof(khmer_KSubsetPartition_Object), /* tp_basicsize */
    0,                                     /* tp_itemsize */
    (destructor)khmer_subset_dealloc,      /*tp_dealloc*/
    0,                                     /*tp_print*/
    0,                                     /*tp_getattr*/
    0,                                     /*tp_setattr*/
    0,                                     /*tp_compare*/
    0,                                     /*tp_repr*/
    0,                                     /*tp_as_number*/
    0,                                     /*tp_as_sequence*/
    0,                                     /*tp_as_mapping*/
    0,                                     /*tp_hash */
    0,                                     /*tp_call*/
    0,                                     /*tp_str*/
    0,                                     /*tp_getattro*/
    0,                                     /*tp_setattro*/
    0,                                     /*tp_as_buffer*/
    Py_TPFLAGS_DEFAULT,                    /*tp_flags*/
    "subset object",                       /* tp_doc */
};

PyMethodDef khmer_subset_methods[] = {
    {
        "count_partitions",
        (PyCFunction)subset_count_partitions,
        METH_VARARGS,
        ""
    },
    {
        "report_on_partitions",
        (PyCFunction)subset_report_on_partitions,
        METH_VARARGS,
        ""
    },
    {
        "partition_size_distribution",
        (PyCFunction)subset_partition_size_distribution,
        METH_VARARGS,
        ""
    },
    {
        "partition_sizes",
        (PyCFunction)subset_partition_sizes,
        METH_VARARGS,
        ""
    },
    {
        "partition_average_coverages",
        (PyCFunction)subset_partition_average_coverages,
        METH_VARARGS,
        ""
    },
    {NULL, NULL, 0, NULL}           /* sentinel */
};

void khmer_subset_dealloc(khmer_KSubsetPartition_Object * obj)
{
    delete obj->subset;
    obj->subset = NULL;
    Py_TYPE(obj)->tp_free((PyObject*)obj);
}



PyObject *
subset_count_partitions(khmer_KSubsetPartition_Object * me, PyObject * args)
{
    SubsetPartition * subset_p = me->subset;

    if (!PyArg_ParseTuple(args, "")) {
        return NULL;
    }

    size_t n_partitions = 0, n_unassigned = 0;
    subset_p->count_partitions(n_partitions, n_unassigned);

    return Py_BuildValue("nn", (Py_ssize_t) n_partitions,
                         (Py_ssize_t) n_unassigned);
}


PyObject *
subset_report_on_partitions(khmer_KSubsetPartition_Object * me, PyObject * args)
{
    SubsetPartition * subset_p = me->subset;

    if (!PyArg_ParseTuple(args, "")) {
        return NULL;
    }

    subset_p->report_on_partitions();

    Py_RETURN_NONE;
}


PyObject *
subset_partition_size_distribution(khmer_KSubsetPartition_Object * me,
                                   PyObject * args)
{
    SubsetPartition * subset_p = me->subset;

    if (!PyArg_ParseTuple(args, "")) {
        return NULL;
    }

    PartitionCountDistribution d;

    unsigned int n_unassigned = 0;
    subset_p->partition_size_distribution(d, n_unassigned);

    PyObject * x = PyList_New(d.size());
    if (x == NULL) {
        return NULL;
    }
    PartitionCountDistribution::iterator di;

    unsigned int i;
    for (i = 0, di = d.begin(); di != d.end(); ++di, i++) {
        PyObject * tup = Py_BuildValue("KK", di->first, di->second);
        if (tup != NULL) {
            PyList_SET_ITEM(x, i, tup);
        }
        Py_XDECREF(tup);
    }
    if (!(i == d.size())) {
        throw oxli_exception();
    }

    PyObject * ret = Py_BuildValue("OI", x, n_unassigned);
    Py_DECREF(x);
    return ret;
}


PyObject *
subset_partition_sizes(khmer_KSubsetPartition_Object * me, PyObject * args)
{
    SubsetPartition * subset_p = me->subset;

    unsigned int min_size = 0;

    if (!PyArg_ParseTuple(args, "|I", &min_size)) {
        return NULL;
    }

    PartitionCountMap cm;
    unsigned int n_unassigned = 0;
    subset_p->partition_sizes(cm, n_unassigned);

    unsigned int i = 0;
    PartitionCountMap::const_iterator mi;
    for (mi = cm.begin(); mi != cm.end(); ++mi) {
        if (mi->second >= min_size) {
            i++;
        }
    }

    PyObject * x = PyList_New(i);
    if (x == NULL) {
        return NULL;
    }

    // this should probably be a dict. @CTB
    for (i = 0, mi = cm.begin(); mi != cm.end(); ++mi) {
        if (mi->second >= min_size) {
            PyObject * tup = Py_BuildValue("II", mi->first, mi->second);
            if (tup != NULL) {
                PyList_SET_ITEM(x, i, tup);
            }
            i++;
        }
    }

    PyObject * ret = Py_BuildValue("OI", x, n_unassigned);
    Py_DECREF(x);

    return ret;
}


PyObject *
subset_partition_average_coverages(khmer_KSubsetPartition_Object * me,
                                   PyObject * args)
{
    SubsetPartition * subset_p = me->subset;

    khmer_KCountgraph_Object * countgraph_o;

    if (!PyArg_ParseTuple(args, "O!", &khmer_KCountgraph_Type, &countgraph_o)) {
        return NULL;
    }

    PartitionCountMap cm;
    subset_p->partition_average_coverages(cm, countgraph_o -> countgraph);

    unsigned int i;
    PartitionCountMap::iterator mi;

    PyObject * x = PyList_New(cm.size());
    if (x == NULL) {
        return NULL;
    }

    // this should probably be a dict. @CTB
    for (i = 0, mi = cm.begin(); mi != cm.end(); ++mi, i++) {
        PyObject * tup = Py_BuildValue("II", mi->first, mi->second);
        if (tup != NULL) {
            PyList_SET_ITEM(x, i, tup);
        }
    }

    return x;
}



}
