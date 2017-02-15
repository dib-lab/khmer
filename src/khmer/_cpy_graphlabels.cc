#include "_cpy_graphlabels.hh"
#include "_cpy_nodegraph.hh"
#include "_cpy_countgraph.hh"
#include "_cpy_hashset.hh"
#include "khmer.hh"
#include "assembler.hh"

namespace khmer {

PyTypeObject khmer_KGraphLabels_Type = {
    PyVarObject_HEAD_INIT(NULL, 0)  /* init & ob_size */
    "_khmer.LabelHash",            /* tp_name */
    sizeof(khmer_KGraphLabels_Object), /* tp_basicsize */
    0,                       /* tp_itemsize */
    (destructor)khmer_graphlabels_dealloc, /* tp_dealloc */
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
    khmer_graphlabels_methods, /* tp_methods */
    0,                       /* tp_members */
    0,                       /* tp_getset */
    0,                       /* tp_base */
    0,                       /* tp_dict */
    0,                       /* tp_descr_get */
    0,                       /* tp_descr_set */
    0,                       /* tp_dictoffset */
    0,                       /* tp_init */
    0,                       /* tp_alloc */
    khmer_graphlabels_new,      /* tp_new */
};


PyMethodDef khmer_graphlabels_methods[] = {
    { "consume_fasta_and_tag_with_labels", (PyCFunction)labelhash_consume_fasta_and_tag_with_labels, METH_VARARGS, "" },
    { "sweep_label_neighborhood", (PyCFunction)labelhash_sweep_label_neighborhood, METH_VARARGS, "" },
    {"consume_partitioned_fasta_and_tag_with_labels", (PyCFunction)labelhash_consume_partitioned_fasta_and_tag_with_labels, METH_VARARGS, "" },
    {"sweep_tag_neighborhood", (PyCFunction)labelhash_sweep_tag_neighborhood, METH_VARARGS, "" },
    {"get_tag_labels", (PyCFunction)labelhash_get_tag_labels, METH_VARARGS, ""},
    {"consume_sequence_and_tag_with_labels", (PyCFunction)labelhash_consume_sequence_and_tag_with_labels, METH_VARARGS, "" },
    {"n_labels", (PyCFunction)labelhash_n_labels, METH_VARARGS, ""},
    {"get_all_labels", (PyCFunction)labelhash_get_all_labels, METH_VARARGS, "" },
    {
        "label_across_high_degree_nodes",
        (PyCFunction)labelhash_label_across_high_degree_nodes, METH_VARARGS,
        "Connect graph across high degree nodes using labels.",
    },
    {
        "assemble_labeled_path",
        (PyCFunction)labelhash_assemble_labeled_path, METH_VARARGS,
        "Assemble all paths, using labels to negotiate tricky bits."
    },
    { "save_labels_and_tags", (PyCFunction)labelhash_save_labels_and_tags, METH_VARARGS, "" },
    { "load_labels_and_tags", (PyCFunction)labelhash_load_labels_and_tags, METH_VARARGS, "" },    {NULL, NULL, 0, NULL}           /* sentinel */
};


void khmer_graphlabels_dealloc(khmer_KGraphLabels_Object * obj)
{
    delete obj->labelhash;
    obj->labelhash = NULL;

    Py_TYPE(obj)->tp_free((PyObject*)obj);
}

 PyObject * khmer_graphlabels_new(PyTypeObject *type, PyObject *args,
                                        PyObject *kwds)
{
    khmer_KGraphLabels_Object *self;
    self = (khmer_KGraphLabels_Object*)type->tp_alloc(type, 0);

    if (self != NULL) {
        PyObject * hashgraph_o;
        khmer::Hashgraph * hashgraph = NULL; // @CTB

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
            self->labelhash = new LabelHash(hashgraph);
        } catch (std::bad_alloc &e) {
            Py_DECREF(self);
            return PyErr_NoMemory();
        }
    }

    return (PyObject *) self;
}


PyObject *
labelhash_get_all_labels(khmer_KGraphLabels_Object * me, PyObject * args)
{
    LabelHash * hb = me->labelhash;

    PyObject * d = PyList_New(hb->all_labels.size());
    if (d == NULL) {
        return NULL;
    }
    LabelSet::iterator it;

    unsigned long long i = 0;
    for (it = hb->all_labels.begin(); it != hb->all_labels.end(); ++it) {
        PyObject * val = Py_BuildValue("K", *it);
        if (val != NULL) {
            PyList_SetItem(d, i, val);
        }
        i++;
    }

    return d;
}


PyObject *
labelhash_consume_fasta_and_tag_with_labels(khmer_KGraphLabels_Object * me,
        PyObject * args)
{
    LabelHash * hb = me->labelhash;

    const char * filename;

    if (!PyArg_ParseTuple(args, "s", &filename)) {
        return NULL;
    }

    const char         *value_exception = NULL;
    const char         *file_exception  = NULL;
    unsigned long long  n_consumed      = 0;
    unsigned int        total_reads     = 0;
    std::string exc_string;

    //Py_BEGIN_ALLOW_THREADS
    try {
        hb->consume_fasta_and_tag_with_labels(filename, total_reads,
                                              n_consumed);
    } catch (khmer_file_exception &exc) {
        exc_string = exc.what();
        file_exception = exc_string.c_str();
    } catch (khmer_value_exception &exc) {
        exc_string = exc.what();
        value_exception = exc_string.c_str();
    }
    //Py_END_ALLOW_THREADS

    if (file_exception != NULL) {
        PyErr_SetString(PyExc_OSError, file_exception);
        return NULL;
    }
    if (value_exception != NULL) {
        PyErr_SetString(PyExc_ValueError, value_exception);
        return NULL;
    }

    return Py_BuildValue("IK", total_reads, n_consumed);
}


PyObject *
labelhash_consume_partitioned_fasta_and_tag_with_labels(
    khmer_KGraphLabels_Object * me, PyObject * args)
{
    LabelHash * labelhash = me->labelhash;

    const char * filename;

    if (!PyArg_ParseTuple(args, "s", &filename)) {
        return NULL;
    }

    // call the C++ function, and trap signals => Python

    unsigned long long  n_consumed  = 0;
    unsigned int        total_reads = 0;

    try {
        labelhash->consume_partitioned_fasta_and_tag_with_labels(filename,
                total_reads, n_consumed);
    } catch (khmer_file_exception &exc) {
        PyErr_SetString(PyExc_OSError, exc.what());
        return NULL;
    } catch (khmer_value_exception &exc) {
        PyErr_SetString(PyExc_ValueError, exc.what());
        return NULL;
    }

    return Py_BuildValue("IK", total_reads, n_consumed);
}


PyObject *
labelhash_consume_sequence_and_tag_with_labels(khmer_KGraphLabels_Object * me,
        PyObject * args)
{
    LabelHash * hb = me->labelhash;
    const char * seq = NULL;
    unsigned long long c = 0;
    if (!PyArg_ParseTuple(args, "sK", &seq, &c)) {
        return NULL;
    }
    unsigned long long n_consumed = 0;

    hb->consume_sequence_and_tag_with_labels(seq, n_consumed, c);
    return Py_BuildValue("K", n_consumed);
}


PyObject *
labelhash_sweep_label_neighborhood(khmer_KGraphLabels_Object * me,
                                   PyObject * args)
{
    LabelHash * hb = me->labelhash;

    const char * seq = NULL;
    int r = 0;
    PyObject * break_on_stop_tags_o = NULL;
    PyObject * stop_big_traversals_o = NULL;

    if (!PyArg_ParseTuple(args, "s|iOO", &seq, &r,
                          &break_on_stop_tags_o,
                          &stop_big_traversals_o)) {
        return NULL;
    }

    unsigned int range = (2 * hb->graph->_get_tag_density()) + 1;
    if (r >= 0) {
        range = r;
    }

    bool break_on_stop_tags = false;
    if (break_on_stop_tags_o && PyObject_IsTrue(break_on_stop_tags_o)) {
        break_on_stop_tags = true;
    }
    bool stop_big_traversals = false;
    if (stop_big_traversals_o && PyObject_IsTrue(stop_big_traversals_o)) {
        stop_big_traversals = true;
    }

    if (strlen(seq) < hb->graph->ksize()) {
        PyErr_SetString(PyExc_ValueError,
                        "string length must >= the hashtable k-mer size");
        return NULL;
    }

    //std::pair<TagLabelPair::iterator, TagLabelPair::iterator> ret;
    LabelSet found_labels;

    //unsigned int num_traversed = 0;
    //Py_BEGIN_ALLOW_THREADS
    hb->sweep_label_neighborhood(seq, found_labels, range, break_on_stop_tags,
                                 stop_big_traversals);
    //Py_END_ALLOW_THREADS

    //printf("...%u kmers traversed\n", num_traversed);

    PyObject * x =  PyList_New(found_labels.size());
    LabelSet::const_iterator si;
    unsigned long long i = 0;
    for (si = found_labels.begin(); si != found_labels.end(); ++si) {
        PyList_SET_ITEM(x, i, Py_BuildValue("K", *si));
        i++;
    }

    return x;
}

// Similar to find_all_tags, but returns tags in a way actually usable by python
// need a tags_in_sequence iterator or function in c++ land for reuse in all
// these functions


PyObject *
labelhash_sweep_tag_neighborhood(khmer_KGraphLabels_Object * me,
                                 PyObject * args)
{
    LabelHash * labelhash = me->labelhash;

    const char * seq = NULL;
    int r = 0;
    PyObject * break_on_stop_tags_o = NULL;
    PyObject * stop_big_traversals_o = NULL;

    if (!PyArg_ParseTuple(args, "s|iOO", &seq, &r,
                          &break_on_stop_tags_o,
                          &stop_big_traversals_o)) {
        return NULL;
    }

    unsigned int range = (2 * labelhash->graph->_get_tag_density()) + 1;
    if (r >= 0) {
        range = r;
    }

    bool break_on_stop_tags = false;
    if (break_on_stop_tags_o && PyObject_IsTrue(break_on_stop_tags_o)) {
        break_on_stop_tags = true;
    }
    bool stop_big_traversals = false;
    if (stop_big_traversals_o && PyObject_IsTrue(stop_big_traversals_o)) {
        stop_big_traversals = true;
    }

    if (strlen(seq) < labelhash->graph->ksize()) {
        PyErr_SetString(PyExc_ValueError,
                        "string length must >= the hashtable k-mer size");
        return NULL;
    }

    SeenSet * tagged_kmers = new SeenSet;

    //Py_BEGIN_ALLOW_THREADS

    labelhash->graph->partition->sweep_for_tags(seq, *tagged_kmers,
            labelhash->graph->all_tags,
            range, break_on_stop_tags,
            stop_big_traversals);

    //Py_END_ALLOW_THREADS

    PyObject * x = (PyObject *) create_HashSet_Object(tagged_kmers,
                   labelhash->graph->ksize());
    return x;
}


PyObject *
labelhash_get_tag_labels(khmer_KGraphLabels_Object * me, PyObject * args)
{
    LabelHash * labelhash = me->labelhash;

    PyObject * tag_o;
    HashIntoType tag;

    if (!PyArg_ParseTuple(args, "O", &tag_o)) {
        return NULL;
    }
    if (!ht_convert_PyObject_to_HashIntoType(tag_o, tag,
            labelhash->graph)) {
        return NULL;
    }

    LabelSet labels;
    labelhash->get_tag_labels(tag, labels);

    PyObject * x =  PyList_New(labels.size());
    LabelSet::const_iterator si;
    unsigned long long i = 0;
    for (si = labels.begin(); si != labels.end(); ++si) {
        PyList_SET_ITEM(x, i, Py_BuildValue("K", *si));
        i++;
    }

    return x;
}


PyObject *
labelhash_n_labels(khmer_KGraphLabels_Object * me, PyObject * args)
{
    LabelHash * labelhash = me->labelhash;

    if (!PyArg_ParseTuple(args, "")) {
        return NULL;
    }

    return PyLong_FromSize_t(labelhash->n_labels());
}


PyObject *
labelhash_label_across_high_degree_nodes(khmer_KGraphLabels_Object * me,
        PyObject * args)
{
    LabelHash * labelhash = me->labelhash;

    const char * long_str;
    khmer_HashSet_Object * hdn_o = NULL;
    Label label;

    if (!PyArg_ParseTuple(args, "sO!K", &long_str,
                          &khmer_HashSet_Type, &hdn_o, &label)) {
        return NULL;
    }

    if (strlen(long_str) < labelhash->graph->ksize()) {
        Py_INCREF(Py_None);
        return Py_None;
    }

    labelhash->label_across_high_degree_nodes(long_str, *hdn_o->hashes, label);

    Py_INCREF(Py_None);
    return Py_None;
}


PyObject *
labelhash_assemble_labeled_path(khmer_KGraphLabels_Object * me,
                                PyObject * args)
{
    LabelHash* labelhash = me->labelhash;

    PyObject * val_o;
    khmer_KNodegraph_Object * nodegraph_o = NULL;
    Nodegraph * stop_bf = NULL;

    if (!PyArg_ParseTuple(args, "O|O!", &val_o,
                          &khmer_KNodegraph_Type, &nodegraph_o)) {
        return NULL;
    }

    Kmer start_kmer;
    if (!ht_convert_PyObject_to_Kmer(val_o, start_kmer, labelhash->graph)) {
        return NULL;
    }

    if (nodegraph_o) {
        stop_bf = nodegraph_o->nodegraph;
    }

    SimpleLabeledAssembler assembler(labelhash);
    std::vector<std::string> contigs = assembler.assemble(start_kmer, stop_bf);

    PyObject * ret = PyList_New(contigs.size());
    for (unsigned int i = 0; i < contigs.size(); i++) {
        PyList_SET_ITEM(ret, i, PyUnicode_FromString(contigs[i].c_str()));
    }

    return ret;
}


PyObject *
labelhash_save_labels_and_tags(khmer_KGraphLabels_Object * me, PyObject * args)
{
    const char * filename = NULL;
    LabelHash * labelhash = me->labelhash;

    if (!PyArg_ParseTuple(args, "s", &filename)) {
        return NULL;
    }

    try {
        labelhash->save_labels_and_tags(filename);
    } catch (khmer_file_exception &e) {
        PyErr_SetString(PyExc_OSError, e.what());
        return NULL;
    }

    Py_RETURN_NONE;
}


PyObject *
labelhash_load_labels_and_tags(khmer_KGraphLabels_Object * me, PyObject * args)
{
    const char * filename = NULL;
    LabelHash * labelhash = me->labelhash;

    if (!PyArg_ParseTuple(args, "s", &filename)) {
        return NULL;
    }

    try {
        labelhash->load_labels_and_tags(filename);
    } catch (khmer_file_exception &e) {
        PyErr_SetString(PyExc_OSError, e.what());
        return NULL;
    }

    Py_RETURN_NONE;
}

}
