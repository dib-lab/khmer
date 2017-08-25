#include "khmer/_cpy_graphlabels.hh"
#include "khmer/_cpy_nodegraph.hh"
#include "khmer/_cpy_countgraph.hh"
#include "khmer/_cpy_hashset.hh"
#include "oxli/oxli.hh"
#include "oxli/assembler.hh"

using namespace oxli;
using namespace oxli::read_parsers;

 
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
    } catch (oxli_file_exception &e) {
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
    } catch (oxli_file_exception &e) {
        PyErr_SetString(PyExc_OSError, e.what());
        return NULL;
    }

    Py_RETURN_NONE;
}

}
