#ifndef _CPY_SUBSET_HH
#define _CPY_SUBSET_HH

#include <Python.h>
#include "_cpy_utils.hh"
#include "subset.hh"

namespace khmer {

typedef struct {
    PyObject_HEAD
    pre_partition_info *   PrePartitionInfo;
} khmer_PrePartitionInfo_Object;

typedef struct {
    PyObject_HEAD
    SubsetPartition * subset;
} khmer_KSubsetPartition_Object;



extern PyTypeObject khmer_PrePartitionInfo_Type;

void khmer_PrePartitionInfo_dealloc(khmer_PrePartitionInfo_Object * obj);


extern PyTypeObject khmer_KSubsetPartition_Type;

extern PyMethodDef khmer_subset_methods[];


void khmer_subset_dealloc(khmer_KSubsetPartition_Object * obj);

PyObject *
subset_count_partitions(khmer_KSubsetPartition_Object * me, PyObject * args);


PyObject *
subset_report_on_partitions(khmer_KSubsetPartition_Object * me, PyObject * args);


PyObject *
subset_partition_size_distribution(khmer_KSubsetPartition_Object * me,
                                   PyObject * args);


PyObject *
subset_partition_sizes(khmer_KSubsetPartition_Object * me, PyObject * args);


PyObject *
subset_partition_average_coverages(khmer_KSubsetPartition_Object * me,
                                   PyObject * args);


}

#endif
