#include "khmer/_cpy_hashset.hh"
#include "khmer/_cpy_hashgraph.hh"
#include "khmer/_cpy_nodegraph.hh"
#include "khmer/_cpy_countgraph.hh"
#include "khmer/_cpy_readparsers.hh"

#include <vector>
#include "oxli/oxli.hh"
#include "oxli/kmer_hash.hh"
#include "oxli/read_parsers.hh"
#include "oxli/assembler.hh"
#include "oxli/traversal.hh"

using namespace oxli;
using namespace oxli::read_parsers;

namespace khmer {

PyTypeObject khmer_KHashgraph_Type
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
    "hashgraph object"            /* tp_doc */
};


PyMethodDef khmer_hashgraph_methods[] = {


    //
    // tagging / sparse graph functionality
    //

    {NULL, NULL, 0, NULL}           /* sentinel */
};







PyObject *
hashgraph_merge_from_disk(khmer_KHashgraph_Object * me, PyObject * args)
{
    Hashgraph * hashgraph = me->hashgraph;

    const char * filename = NULL;
    if (!PyArg_ParseTuple(args, "s", &filename)) {
        return NULL;
    }

    try {
        hashgraph->partition->merge_from_disk(filename);
    } catch (oxli_file_exception &e) {
        PyErr_SetString(PyExc_OSError, e.what());
        return NULL;
    }

    Py_RETURN_NONE;
}


PyObject *
hashgraph_find_all_tags(khmer_KHashgraph_Object * me, PyObject * args)
{
    Hashgraph * hashgraph = me->hashgraph;

    const char * kmer_s = NULL;

    if (!PyArg_ParseTuple(args, "s", &kmer_s)) {
        return NULL;
    }

    if (strlen(kmer_s) != hashgraph->ksize()) {
        PyErr_SetString( PyExc_ValueError,
                         "k-mer size must equal the k-mer size of the presence table");
        return NULL;
    }

    pre_partition_info * ppi = NULL;

    Kmer kmer = hashgraph->build_kmer(kmer_s);

    Py_BEGIN_ALLOW_THREADS

    try {
        ppi = new pre_partition_info(kmer);
    } catch (std::bad_alloc &e) {
        return PyErr_NoMemory();
    }
    hashgraph->partition->find_all_tags(kmer, ppi->tagged_kmers,
                                        hashgraph->all_tags);
    hashgraph->add_kmer_to_tags(kmer);

    Py_END_ALLOW_THREADS

    khmer_PrePartitionInfo_Object * ppi_obj = (khmer_PrePartitionInfo_Object *) \
            PyObject_New(khmer_PrePartitionInfo_Object, &khmer_PrePartitionInfo_Type);

    ppi_obj->PrePartitionInfo = ppi;

    return (PyObject*)ppi_obj;
}


PyObject *
hashgraph_assign_partition_id(khmer_KHashgraph_Object * me, PyObject * args)
{
    Hashgraph * hashgraph = me->hashgraph;

    khmer_PrePartitionInfo_Object * ppi_obj;
    if (!PyArg_ParseTuple(args, "O!", &khmer_PrePartitionInfo_Type, &ppi_obj)) {
        return NULL;
    }

    pre_partition_info * ppi;
    ppi = ppi_obj->PrePartitionInfo;

    PartitionID p;
    p = hashgraph->partition->assign_partition_id(ppi->kmer,
            ppi->tagged_kmers);

    return PyLong_FromLong(p);
}




PyObject *
hashgraph_output_partitions(khmer_KHashgraph_Object * me, PyObject * args)
{
    Hashgraph * hashgraph = me->hashgraph;

    const char * filename = NULL;
    const char * output = NULL;
    PyObject * output_unassigned_o = NULL;

    if (!PyArg_ParseTuple(args, "ss|O", &filename, &output,
                          &output_unassigned_o)) {
        return NULL;
    }

    bool output_unassigned = false;
    if (output_unassigned_o != NULL && PyObject_IsTrue(output_unassigned_o)) {
        output_unassigned = true;
    }

    size_t n_partitions = 0;

    try {
        SubsetPartition * subset_p = hashgraph->partition;
        n_partitions = subset_p->output_partitioned_file(filename,
                       output,
                       output_unassigned);
    } catch (oxli_file_exception &e) {
        PyErr_SetString(PyExc_OSError, e.what());
        return NULL;
    } catch (oxli_value_exception &exc) {
        PyErr_SetString(PyExc_ValueError, exc.what());
        return NULL;
    }

    return PyLong_FromLong(n_partitions);
}


PyObject *
hashgraph_save_partitionmap(khmer_KHashgraph_Object * me, PyObject * args)
{
    Hashgraph * hashgraph = me->hashgraph;

    const char * filename = NULL;

    if (!PyArg_ParseTuple(args, "s", &filename)) {
        return NULL;
    }

    try {
        hashgraph->partition->save_partitionmap(filename);
    } catch (oxli_file_exception &e) {
        PyErr_SetString(PyExc_OSError, e.what());
        return NULL;
    }

    Py_RETURN_NONE;
}


PyObject *
hashgraph_load_partitionmap(khmer_KHashgraph_Object * me, PyObject * args)
{
    Hashgraph * hashgraph = me->hashgraph;

    const char * filename = NULL;

    if (!PyArg_ParseTuple(args, "s", &filename)) {
        return NULL;
    }

    try {
        hashgraph->partition->load_partitionmap(filename);
    } catch (oxli_file_exception &e) {
        PyErr_SetString(PyExc_OSError, e.what());
        return NULL;
    }

    Py_RETURN_NONE;
}


PyObject *
hashgraph__validate_partitionmap(khmer_KHashgraph_Object * me, PyObject * args)
{
    Hashgraph * hashgraph = me->hashgraph;

    if (!PyArg_ParseTuple(args, "")) {
        return NULL;
    }

    hashgraph->partition->_validate_pmap();

    Py_RETURN_NONE;
}



PyObject *
hashgraph_subset_count_partitions(khmer_KHashgraph_Object * me, PyObject * args)
{
    khmer_KSubsetPartition_Object * subset_obj = NULL;

    if (!PyArg_ParseTuple(args, "O!", &khmer_KSubsetPartition_Type,
                          &subset_obj)) {
        return NULL;
    }


    size_t n_partitions = 0, n_unassigned = 0;
    subset_obj->subset->count_partitions(n_partitions, n_unassigned);

    return Py_BuildValue("nn", (Py_ssize_t) n_partitions,
                         (Py_ssize_t) n_unassigned);
}


PyObject *
hashgraph_subset_partition_size_distribution(khmer_KHashgraph_Object * me,
        PyObject * args)
{
    khmer_KSubsetPartition_Object * subset_obj = NULL;
    if (!PyArg_ParseTuple(args, "O!", &khmer_KSubsetPartition_Type,
                          &subset_obj)) {
        return NULL;
    }

    SubsetPartition * subset_p = subset_obj->subset;

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
        PyObject * value =  Py_BuildValue("KK", di->first, di->second);
        if (value == NULL) {
            Py_DECREF(x);
            return NULL;
        }
        PyList_SET_ITEM(x, i, value);
    }
    if (!(i == d.size())) {
        throw oxli_exception();
    }

    PyObject * returnValue = Py_BuildValue("NI", x, n_unassigned);
    if (returnValue == NULL) {
        Py_DECREF(x);
        return NULL;
    }
    return returnValue;
}


PyObject *
hashgraph_save_subset_partitionmap(khmer_KHashgraph_Object * me,
                                   PyObject * args)
{
    const char * filename = NULL;
    khmer_KSubsetPartition_Object * subset_obj = NULL;

    if (!PyArg_ParseTuple(args, "O!s", &khmer_KSubsetPartition_Type,
                          &subset_obj, &filename)) {
        return NULL;
    }

    SubsetPartition * subset_p = subset_obj->subset;

    Py_BEGIN_ALLOW_THREADS

    try {
        subset_p->save_partitionmap(filename);
    } catch (oxli_file_exception &e) {
        PyErr_SetString(PyExc_OSError, e.what());
        return NULL;
    }

    Py_END_ALLOW_THREADS

    Py_RETURN_NONE;
}


PyObject *
hashgraph_load_subset_partitionmap(khmer_KHashgraph_Object * me,
                                   PyObject * args)
{
    Hashgraph * hashgraph = me->hashgraph;

    const char * filename = NULL;

    if (!PyArg_ParseTuple(args, "s", &filename)) {
        return NULL;
    }

    SubsetPartition * subset_p;
    try {
        subset_p = new SubsetPartition(hashgraph);
    } catch (std::bad_alloc &e) {
        return PyErr_NoMemory();
    }

    const char         *file_exception  = NULL;

    std::string exc_string ;
    Py_BEGIN_ALLOW_THREADS
    try {
        subset_p->load_partitionmap(filename);
    } catch (oxli_file_exception &exc) {
        exc_string = exc.what();
        file_exception = exc_string.c_str();
    }
    Py_END_ALLOW_THREADS

    if (file_exception != NULL) {
        PyErr_SetString(PyExc_OSError, file_exception);
        delete subset_p;
        return NULL;
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


PyObject *
hashgraph__validate_subset_partitionmap(khmer_KHashgraph_Object * me,
                                        PyObject * args)
{
    khmer_KSubsetPartition_Object * subset_obj = NULL;

    if (!PyArg_ParseTuple(args, "O!", &khmer_KSubsetPartition_Type,
                          &subset_obj)) {
        return NULL;
    }

    SubsetPartition * subset_p = subset_obj->subset;

    subset_p->_validate_pmap();

    Py_RETURN_NONE;
}




PyObject *
hashgraph_join_partitions(khmer_KHashgraph_Object * me, PyObject * args)
{
    Hashgraph * hashgraph = me->hashgraph;

    PartitionID p1 = 0, p2 = 0;

    if (!PyArg_ParseTuple(args, "II", &p1, &p2)) {
        return NULL;
    }

    p1 = hashgraph->partition->join_partitions(p1, p2);

    return PyLong_FromLong(p1);
}




}
