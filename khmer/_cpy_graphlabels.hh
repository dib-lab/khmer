#ifndef _CPY_GRAPHLABELS_HH
#define _CPY_GRAPHLABELS_HH

#include <Python.h>
#include "_cpy_utils.hh"
#include "labelhash.hh"

namespace khmer {

typedef struct {
    PyObject_HEAD
    LabelHash * labelhash;
} khmer_KGraphLabels_Object;


extern PyTypeObject khmer_KGraphLabels_Type;

extern PyMethodDef khmer_graphlabels_methods[];

PyObject * khmer_graphlabels_new(PyTypeObject * type, PyObject *args,
                                        PyObject *kwds);

void khmer_graphlabels_dealloc(khmer_KGraphLabels_Object * obj);


PyObject *
labelhash_get_all_labels(khmer_KGraphLabels_Object * me, PyObject * args);


PyObject *
labelhash_consume_fasta_and_tag_with_labels(khmer_KGraphLabels_Object * me,
        PyObject * args);


PyObject *
labelhash_consume_partitioned_fasta_and_tag_with_labels(
    khmer_KGraphLabels_Object * me, PyObject * args);


PyObject *
labelhash_consume_sequence_and_tag_with_labels(khmer_KGraphLabels_Object * me,
        PyObject * args);


PyObject *
labelhash_sweep_label_neighborhood(khmer_KGraphLabels_Object * me,
                                   PyObject * args);

// Similar to find_all_tags, but returns tags in a way actually usable by python
// need a tags_in_sequence iterator or function in c++ land for reuse in all
// these functions


PyObject *
labelhash_sweep_tag_neighborhood(khmer_KGraphLabels_Object * me,
                                 PyObject * args);


PyObject *
labelhash_get_tag_labels(khmer_KGraphLabels_Object * me, PyObject * args);


PyObject *
labelhash_n_labels(khmer_KGraphLabels_Object * me, PyObject * args);


PyObject *
labelhash_label_across_high_degree_nodes(khmer_KGraphLabels_Object * me,
        PyObject * args);


PyObject *
labelhash_assemble_labeled_path(khmer_KGraphLabels_Object * me,
                                PyObject * args);


PyObject *
labelhash_save_labels_and_tags(khmer_KGraphLabels_Object * me, PyObject * args);


PyObject *
labelhash_load_labels_and_tags(khmer_KGraphLabels_Object * me, PyObject * args);

}

#endif
