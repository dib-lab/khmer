//
// This file is part of khmer, http://github.com/ged-lab/khmer/, and is
// Copyright (C) Michigan State University, 2009-2013. It is licensed under
// the three-clause BSD license; see doc/LICENSE.txt.
// Contact: khmer-project@idyll.org
//

//
// A module for Python that exports khmer C++ library functions.
//

// Must be first.

#ifndef KHMERMODULE_HH
#define KHMERMODULE_HH

#include <Python.h>
#include "hashtable.hh"

#include <iostream>

#include "khmer.hh"
#include "kmer_hash.hh"
#include "hashbits.hh"
#include "counting.hh"
#include "read_aligner.hh"
#include "labelhash.hh"
#include "khmer_async.hh"
#include "khmer_exception.hh"

using namespace khmer;

// Configure module logging.
//#define WITH_INTERNAL_TRACING
namespace khmer
{

namespace python
{

template < typename OBJECT >
void
_common_init_Type(PyTypeObject &tobj, char const * name, char const * doc);


} // namespace python

} // namespace khmer


class _khmer_exception
{
private:
    std::string _message;
public:
    _khmer_exception(std::string message) : _message(message) { };
    inline const std::string get_message() const
    {
        return _message;
    };
};

class _khmer_signal : public _khmer_exception
{
public:
    _khmer_signal(std::string message) : _khmer_exception(message) { };
};

typedef pre_partition_info _pre_partition_info;

// default callback obj;
static PyObject *_callback_obj = NULL;

// callback function to pass into C++ functions

void _report_fn(const char * info, void * data, unsigned long long n_reads,
                unsigned long long other);
/***********************************************************************/

//
// Read object -- name, sequence, and FASTQ stuff
//

namespace khmer
{

namespace python
{

static PyTypeObject Read_Type = { PyObject_HEAD_INIT( NULL ) };

typedef struct {
    PyObject_HEAD
    //! Pointer to the low-level genomic read object.
    read_parsers:: Read *   read;
} Read_Object;

static
void
_Read_dealloc( PyObject * self );

#define KHMER_READ_STRING_GETTER( SELF, ATTR_NAME ) \
    PyString_FromString( \
    ((((Read_Object *)(SELF))->read)->ATTR_NAME).c_str( ) \
    )


static
PyObject *
Read_get_name( PyObject * self, void * closure );

static
PyObject *
Read_get_sequence( PyObject * self, void * closure );

static
PyObject *
Read_get_accuracy( PyObject * self, void * closure );

static
PyObject *
Read_get_annotations( PyObject * self, void * closure );

static PyGetSetDef _Read_accessors [ ] = {
    {
        (char *)"name",
        (getter)Read_get_name, (setter)NULL,
        (char *)"Read identifier.", NULL
    },
    {
        (char *)"sequence",
        (getter)Read_get_sequence, (setter)NULL,
        (char *)"Genomic sequence.", NULL
    },
    {
        (char *)"accuracy",
        (getter)Read_get_accuracy, (setter)NULL,
        (char *)"Quality scores.", NULL
    },
    {
        (char *)"annotations",
        (getter)Read_get_annotations, (setter)NULL,
        (char *)"Annotations.", NULL
    },

    { NULL, NULL, NULL, NULL, NULL } // sentinel
};

static void _init_Read_Type( );

/***********************************************************************/

//
// ReadParser object -- parse reads directly from streams
// ReadPairIterator -- return pairs of Read objects
//


static PyTypeObject ReadParser_Type
CPYCHECKER_TYPE_OBJECT_FOR_TYPEDEF("ReadParser_Object")
    = { PyObject_HEAD_INIT( NULL ) };
static PyTypeObject ReadPairIterator_Type = { PyObject_HEAD_INIT( NULL ) };

typedef struct {
    PyObject_HEAD
    //! Pointer to the low-level parser object.
    read_parsers:: IParser *  parser;
} ReadParser_Object;


typedef struct {
    PyObject_HEAD
    //! Pointer to Python parser object for reference counting purposes.
    PyObject *  parent;
    //! Persistent value of pair mode across invocations.
    int pair_mode;
} ReadPairIterator_Object;


static
void
_ReadParser_dealloc( PyObject * self );


static
void
_ReadPairIterator_dealloc( PyObject * self );

static
PyObject *
_ReadParser_new( PyTypeObject * subtype, PyObject * args, PyObject * kwds );

static
PyObject *
_ReadParser_iternext( PyObject * self );

static
PyObject *
_ReadPairIterator_iternext( PyObject * self );


static
PyObject *
ReadParser_iter_reads( PyObject * self, PyObject * args );

static
PyObject *
ReadParser_iter_read_pairs( PyObject * self, PyObject * args );

static PyMethodDef _ReadParser_methods [ ] = {
    {
        "iter_reads",       (PyCFunction)ReadParser_iter_reads,
        METH_NOARGS,        "Iterates over reads."
    },
    {
        "iter_read_pairs",  (PyCFunction)ReadParser_iter_read_pairs,
        METH_VARARGS,       "Iterates over paired reads as pairs."
    },

    { NULL, NULL, 0, NULL } // sentinel
};


static
void
_init_ReadParser_Type( );

static
void
_init_ReadPairIterator_Type( );

} // namespace python

} // namespace khmer


static
read_parsers:: IParser *
_PyObject_to_khmer_ReadParser( PyObject * py_object );

/***********************************************************************/

//
// KCountingHash object
//

void free_pre_partition_info(void * p);

void free_subset_partition_info(void * p);

typedef struct {
    PyObject_HEAD
    CountingHash * counting;
} khmer_KCountingHashObject;

typedef struct {
    PyObject_HEAD
    SubsetPartition * subset;
} khmer_KSubsetPartitionObject;

typedef struct {
    PyObject_HEAD
    Hashbits * hashbits;
} khmer_KHashbitsObject;

static void khmer_subset_dealloc(PyObject *);
static PyObject * khmer_subset_getattr(PyObject * obj, char * name);

static PyTypeObject khmer_KSubsetPartitionType = {
    PyObject_HEAD_INIT(NULL)
    0,
    "KSubset", sizeof(khmer_KSubsetPartitionObject),
    0,
    khmer_subset_dealloc,   /*tp_dealloc*/
    0,              /*tp_print*/
    khmer_subset_getattr,   /*tp_getattr*/
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
    Py_TPFLAGS_DEFAULT,     /*tp_flags*/
    "subset object",           /* tp_doc */
};

typedef struct {
    PyObject_HEAD
    ReadAligner * aligner;
} khmer_ReadAlignerObject;

static void khmer_counting_dealloc(PyObject *);
static PyObject * hash_abundance_distribution(PyObject * self,
        PyObject * args);
static PyObject * hash_abundance_distribution_with_reads_parser(
        PyObject * self, PyObject * args);
static PyObject * hash_set_use_bigcount(PyObject * self, PyObject * args);
static PyObject * hash_get_use_bigcount(PyObject * self, PyObject * args);
static PyObject * hash_n_occupied(PyObject * self, PyObject * args);
static PyObject * hash_n_entries(PyObject * self, PyObject * args);
static PyObject * hash_count(PyObject * self, PyObject * args);
static PyObject * hash_output_fasta_kmer_pos_freq(PyObject * self,
        PyObject *args);
static PyObject * hash_consume_fasta(PyObject * self, PyObject * args);
static PyObject * hash_consume_fasta_with_reads_parser(
        PyObject * self, PyObject * args);
static PyObject * hash_consume(PyObject * self, PyObject * args);
static PyObject * hash_get_min_count(PyObject * self, PyObject * args);
static PyObject * hash_get_max_count(PyObject * self, PyObject * args);
static PyObject * hash_get_median_count(PyObject * self, PyObject * args);
static PyObject * hash_get_kadian_count(PyObject * self, PyObject * args);
static PyObject * hash_get(PyObject * self, PyObject * args);
static PyObject * count_trim_on_abundance(PyObject * self, PyObject * args);
static PyObject * count_trim_below_abundance(PyObject * self, PyObject * args);
static PyObject * hash_fasta_count_kmers_by_position(PyObject * self,
        PyObject * args);
static PyObject * hash_fasta_dump_kmers_by_abundance(PyObject * self,
        PyObject * args);
static PyObject * hash_load(PyObject * self, PyObject * args);
static PyObject * hash_save(PyObject * self, PyObject * args);
static PyObject * hash_get_ksize(PyObject * self, PyObject * args);
static PyObject * hash_get_hashsizes(PyObject * self, PyObject * args);
static PyObject * hash_collect_high_abundance_kmers(PyObject * self,
        PyObject * args);
static PyObject * hash_consume_and_tag(PyObject * self, PyObject * args);
static PyObject * hash_consume_fasta_and_tag(PyObject * self, PyObject * args);
static PyObject * hash_find_all_tags_truncate_on_abundance(PyObject * self,
        PyObject *args);
static PyObject * hash_do_subset_partition_with_abundance(PyObject * self,
        PyObject * args);

static PyMethodDef khmer_counting_methods[] = {
    { "ksize", hash_get_ksize, METH_VARARGS, "" },
    { "hashsizes", hash_get_hashsizes, METH_VARARGS, "" },
    { "set_use_bigcount", hash_set_use_bigcount, METH_VARARGS, "" },
    { "get_use_bigcount", hash_get_use_bigcount, METH_VARARGS, "" },
    { "n_occupied", hash_n_occupied, METH_VARARGS, "Count the number of occupied bins" },
    { "n_entries", hash_n_entries, METH_VARARGS, "" },
    { "count", hash_count, METH_VARARGS, "Count the given kmer" },
    { "consume", hash_consume, METH_VARARGS, "Count all k-mers in the given string" },
    { "consume_fasta", hash_consume_fasta, METH_VARARGS, "Count all k-mers in a given file" },
    {
        "consume_fasta_with_reads_parser", hash_consume_fasta_with_reads_parser,
        METH_VARARGS, "Count all k-mers using a given reads parser"
    },
    { "output_fasta_kmer_pos_freq", hash_output_fasta_kmer_pos_freq, METH_VARARGS, "" },
    { "get", hash_get, METH_VARARGS, "Get the count for the given k-mer" },
    { "get_min_count", hash_get_min_count, METH_VARARGS, "Get the smallest count of all the k-mers in the string" },
    { "get_max_count", hash_get_max_count, METH_VARARGS, "Get the largest count of all the k-mers in the string" },
    { "get_median_count", hash_get_median_count, METH_VARARGS, "Get the median, average, and stddev of the k-mer counts in the string" },
    { "get_kadian_count", hash_get_kadian_count, METH_VARARGS, "Get the kadian (abundance of k-th rank-ordered k-mer) of the k-mer counts in the string" },
    { "trim_on_abundance", count_trim_on_abundance, METH_VARARGS, "Trim on >= abundance" },
    { "trim_below_abundance", count_trim_below_abundance, METH_VARARGS, "Trim on >= abundance" },
    { "abundance_distribution", hash_abundance_distribution, METH_VARARGS, "" },
    { "abundance_distribution_with_reads_parser", hash_abundance_distribution_with_reads_parser, METH_VARARGS, "" },
    { "fasta_count_kmers_by_position", hash_fasta_count_kmers_by_position, METH_VARARGS, "" },
    { "fasta_dump_kmers_by_abundance", hash_fasta_dump_kmers_by_abundance, METH_VARARGS, "" },
    { "load", hash_load, METH_VARARGS, "" },
    { "save", hash_save, METH_VARARGS, "" },
    {
        "collect_high_abundance_kmers", hash_collect_high_abundance_kmers,
        METH_VARARGS, ""
    },
    { "consume_and_tag", hash_consume_and_tag, METH_VARARGS, "Consume a sequence and tag it" },
    { "consume_fasta_and_tag", hash_consume_fasta_and_tag, METH_VARARGS, "Count all k-mers in a given file" },
    { "do_subset_partition_with_abundance", hash_do_subset_partition_with_abundance, METH_VARARGS, "" },
    { "find_all_tags_truncate_on_abundance", hash_find_all_tags_truncate_on_abundance, METH_VARARGS, "" },

    {NULL, NULL, 0, NULL}           /* sentinel */
};

static PyObject *
khmer_counting_getattr(PyObject * obj, char * name);
#define is_counting_obj(v)  ((v)->ob_type == &khmer_KCountingHashType)

static PyTypeObject khmer_KCountingHashType
CPYCHECKER_TYPE_OBJECT_FOR_TYPEDEF("khmer_KCountingHashObject")
= {
    PyObject_HEAD_INIT(NULL)
    0,
    "KCountingHash", sizeof(khmer_KCountingHashObject),
    0,
    khmer_counting_dealloc, /*tp_dealloc*/
    0,              /*tp_print*/
    khmer_counting_getattr, /*tp_getattr*/
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
    Py_TPFLAGS_DEFAULT,     /*tp_flags*/
    "counting hash object",           /* tp_doc */
};

//
// new_hashtable
//

static PyObject* new_hashtable(PyObject * self, PyObject * args);
//
// new_counting_hash
//

static PyObject* _new_counting_hash(PyObject * self, PyObject * args);
//
// hashbits stuff
//

static void khmer_hashbits_dealloc(PyObject * obj);
static PyObject* khmer_hashbits_new(PyTypeObject * type, PyObject * args,
                                    PyObject * kwds);
static int khmer_hashbits_init(khmer_KHashbitsObject * self, PyObject * args,
                               PyObject * kwds);
static PyObject * khmer_hashbits_getattr(PyObject * obj, char * name);

static PyTypeObject khmer_KHashbitsType
CPYCHECKER_TYPE_OBJECT_FOR_TYPEDEF("khmer_KHashbitsObject")
= {
    PyObject_HEAD_INIT(NULL)
    0,
    "Hashbits", sizeof(khmer_KHashbitsObject),
    0,
    (destructor)khmer_hashbits_dealloc, /*tp_dealloc*/
    0,              /*tp_print*/
    khmer_hashbits_getattr, /*tp_getattr*/
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
    "hashbits object",           /* tp_doc */
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
    (initproc)khmer_hashbits_init,   /* tp_init */
    0,                       /* tp_alloc */
};

static PyObject * hash_abundance_distribution_with_reads_parser(
    PyObject * self, PyObject * args);
static PyObject * hash_abundance_distribution(PyObject * self, PyObject * args);
static PyObject * hashbits_n_unique_kmers(PyObject * self, PyObject * args);
static PyObject * hashbits_count_overlap(PyObject * self, PyObject * args);
static PyObject * hashbits_n_occupied(PyObject * self, PyObject * args);
static PyObject * hashbits_n_tags(PyObject * self, PyObject * args);
static PyObject * hashbits_count(PyObject * self, PyObject * args);
static PyObject * hashbits_consume(PyObject * self, PyObject * args);
static PyObject * hashbits_print_stop_tags(PyObject * self, PyObject * args);
static PyObject * hashbits_print_tagset(PyObject * self, PyObject * args);
static PyObject * hashbits_load_stop_tags(PyObject * self, PyObject * args);
static PyObject * hashbits_save_stop_tags(PyObject * self, PyObject * args);
static PyObject * hashbits_traverse_from_tags(PyObject * self, PyObject * args);
static PyObject * hashbits_repartition_largest_partition(PyObject * self,
        PyObject * args);
static PyObject * hashbits_get(PyObject * self, PyObject * args);
static PyObject * hashbits_calc_connected_graph_size(PyObject * self,
        PyObject * args);
static PyObject * hashbits_kmer_degree(PyObject * self, PyObject * args);
static PyObject * hashbits_trim_on_stoptags(PyObject * self, PyObject * args);
static PyObject * hashbits_identify_stoptags_by_position(PyObject * self,
        PyObject * args);
static PyObject * hashbits_do_subset_partition(PyObject * self,
        PyObject * args);
static PyObject * hashbits_join_partitions_by_path(PyObject * self,
        PyObject *args);
static PyObject * hashbits_merge_subset(PyObject * self, PyObject *args);
static PyObject * hashbits_merge_from_disk(PyObject * self, PyObject *args);
static PyObject * hashbits_consume_fasta(PyObject * self, PyObject * args);
static PyObject * hashbits_consume_fasta_with_reads_parser(
        PyObject * self, PyObject * args);
static PyObject * hashbits_consume_fasta_and_traverse(PyObject * self,
        PyObject * args);
void sig(unsigned int total_reads, unsigned int n_consumed);
static PyObject * hashbits_consume_fasta_and_tag(PyObject * self,
        PyObject * args);
static PyObject * hashbits_consume_fasta_and_tag_with_reads_parser(
        PyObject * self, PyObject * args);
static PyObject * hashbits_consume_fasta_and_tag_with_stoptags(PyObject * self,
        PyObject * args);
static PyObject * hashbits_consume_partitioned_fasta(PyObject * self,
        PyObject * args);
static PyObject * hashbits_find_all_tags(PyObject * self, PyObject *args);
static PyObject * hashbits_assign_partition_id(PyObject * self, PyObject *args);
static PyObject * hashbits_add_tag(PyObject * self, PyObject *args);
static PyObject * hashbits_add_stop_tag(PyObject * self, PyObject *args);
static PyObject * hashbits_get_stop_tags(PyObject * self, PyObject * args);
static PyObject * hashbits_get_tagset(PyObject * self, PyObject * args);
static PyObject * hashbits_output_partitions(PyObject * self, PyObject * args);
static PyObject * hashbits_find_unpart(PyObject * self, PyObject * args);
static PyObject * hashbits_filter_if_present(PyObject * self, PyObject * args);
static PyObject * hashbits_save_partitionmap(PyObject * self, PyObject * args);
static PyObject * hashbits_load_partitionmap(PyObject * self, PyObject * args);
static PyObject * hashbits__validate_partitionmap(PyObject * self,
        PyObject * args);
static PyObject * hashbits_count_partitions(PyObject * self, PyObject * args);
static PyObject * hashbits_subset_count_partitions(PyObject * self,
        PyObject * args);
static PyObject * hashbits_subset_partition_size_distribution(PyObject * self,
        PyObject * args);
static PyObject * hashbits_load(PyObject * self, PyObject * args);
static PyObject * hashbits_save(PyObject * self, PyObject * args);
static PyObject * hashbits_load_tagset(PyObject * self, PyObject * args);
static PyObject * hashbits_save_tagset(PyObject * self, PyObject * args);
static PyObject * hashbits_save_subset_partitionmap(PyObject * self,
        PyObject * args);
static PyObject * hashbits_load_subset_partitionmap(PyObject * self,
        PyObject * args);
static PyObject * hashbits__set_tag_density(PyObject * self, PyObject * args);
static PyObject * hashbits__get_tag_density(PyObject * self, PyObject * args);
static PyObject * hashbits__validate_subset_partitionmap(PyObject * self,
        PyObject * args);
static PyObject * hashbits_set_partition_id(PyObject * self, PyObject * args);
static PyObject * hashbits_join_partitions(PyObject * self, PyObject * args);
static PyObject * hashbits_get_partition_id(PyObject * self, PyObject * args);
static PyObject * hashbits_is_single_partition(PyObject * self,
        PyObject * args);
static PyObject * hashbits_divide_tags_into_subsets(PyObject * self,
        PyObject * args);
static PyObject * hashbits_count_kmers_within_radius(PyObject * self,
        PyObject * args);
static PyObject * hashbits_get_ksize(PyObject * self, PyObject * args);
static PyObject * hashbits_get_hashsizes(PyObject * self, PyObject * args);
static PyObject * hashbits_extract_unique_paths(PyObject * self,
        PyObject * args);
static PyObject * hashbits_get_median_count(PyObject * self, PyObject * args);

static PyMethodDef khmer_hashbits_methods[] = {
    { "extract_unique_paths", hashbits_extract_unique_paths, METH_VARARGS, "" },
    { "ksize", hashbits_get_ksize, METH_VARARGS, "" },
    { "hashsizes", hashbits_get_hashsizes, METH_VARARGS, "" },
    { "n_occupied", hashbits_n_occupied, METH_VARARGS, "Count the number of occupied bins" },
    { "n_unique_kmers", hashbits_n_unique_kmers,  METH_VARARGS, "Count the number of unique kmers" },
    { "count", hashbits_count, METH_VARARGS, "Count the given kmer" },
    { "count_overlap", hashbits_count_overlap, METH_VARARGS, "Count overlap kmers in two datasets" },
    { "consume", hashbits_consume, METH_VARARGS, "Count all k-mers in the given string" },
    { "load_stop_tags", hashbits_load_stop_tags, METH_VARARGS, "" },
    { "save_stop_tags", hashbits_save_stop_tags, METH_VARARGS, "" },
    { "print_stop_tags", hashbits_print_stop_tags, METH_VARARGS, "" },
    { "print_tagset", hashbits_print_tagset, METH_VARARGS, "" },
    { "get", hashbits_get, METH_VARARGS, "Get the count for the given k-mer" },
    { "calc_connected_graph_size", hashbits_calc_connected_graph_size, METH_VARARGS, "" },
    { "kmer_degree", hashbits_kmer_degree, METH_VARARGS, "" },
    { "trim_on_stoptags", hashbits_trim_on_stoptags, METH_VARARGS, "" },
    { "identify_stoptags_by_position", hashbits_identify_stoptags_by_position, METH_VARARGS, "" },
    { "do_subset_partition", hashbits_do_subset_partition, METH_VARARGS, "" },
    { "find_all_tags", hashbits_find_all_tags, METH_VARARGS, "" },
    { "assign_partition_id", hashbits_assign_partition_id, METH_VARARGS, "" },
    { "output_partitions", hashbits_output_partitions, METH_VARARGS, "" },
    { "find_unpart", hashbits_find_unpart, METH_VARARGS, "" },
    { "filter_if_present", hashbits_filter_if_present, METH_VARARGS, "" },
    { "add_tag", hashbits_add_tag, METH_VARARGS, "" },
    { "add_stop_tag", hashbits_add_stop_tag, METH_VARARGS, "" },
    { "get_stop_tags", hashbits_get_stop_tags, METH_VARARGS, "" },
    { "get_tagset", hashbits_get_tagset, METH_VARARGS, "" },
    { "load", hashbits_load, METH_VARARGS, "" },
    { "save", hashbits_save, METH_VARARGS, "" },
    { "load_tagset", hashbits_load_tagset, METH_VARARGS, "" },
    { "save_tagset", hashbits_save_tagset, METH_VARARGS, "" },
    { "n_tags", hashbits_n_tags, METH_VARARGS, "" },
    { "divide_tags_into_subsets", hashbits_divide_tags_into_subsets, METH_VARARGS, "" },
    { "load_partitionmap", hashbits_load_partitionmap, METH_VARARGS, "" },
    { "save_partitionmap", hashbits_save_partitionmap, METH_VARARGS, "" },
    { "_validate_partitionmap", hashbits__validate_partitionmap, METH_VARARGS, "" },
    { "_get_tag_density", hashbits__get_tag_density, METH_VARARGS, "" },
    { "_set_tag_density", hashbits__set_tag_density, METH_VARARGS, "" },
    { "consume_fasta", hashbits_consume_fasta, METH_VARARGS, "Count all k-mers in a given file" },
    { "consume_fasta_with_reads_parser", hashbits_consume_fasta_with_reads_parser, METH_VARARGS, "Count all k-mers in a given file" },
    { "consume_fasta_and_tag", hashbits_consume_fasta_and_tag, METH_VARARGS, "Count all k-mers in a given file" },
    {
        "consume_fasta_and_tag_with_reads_parser", hashbits_consume_fasta_and_tag_with_reads_parser,
        METH_VARARGS, "Count all k-mers using a given reads parser"
    },
    { "consume_fasta_and_traverse", hashbits_consume_fasta_and_traverse, METH_VARARGS, "" },
    { "consume_fasta_and_tag_with_stoptags", hashbits_consume_fasta_and_tag_with_stoptags, METH_VARARGS, "Count all k-mers in a given file" },
    { "consume_partitioned_fasta", hashbits_consume_partitioned_fasta, METH_VARARGS, "Count all k-mers in a given file" },
    { "join_partitions_by_path", hashbits_join_partitions_by_path, METH_VARARGS, "" },
    { "merge_subset", hashbits_merge_subset, METH_VARARGS, "" },
    { "merge_subset_from_disk", hashbits_merge_from_disk, METH_VARARGS, "" },
    { "count_partitions", hashbits_count_partitions, METH_VARARGS, "" },
    { "subset_count_partitions", hashbits_subset_count_partitions, METH_VARARGS, "" },
    { "subset_partition_size_distribution", hashbits_subset_partition_size_distribution, METH_VARARGS, "" },
    { "save_subset_partitionmap", hashbits_save_subset_partitionmap, METH_VARARGS },
    { "load_subset_partitionmap", hashbits_load_subset_partitionmap, METH_VARARGS },
    { "_validate_subset_partitionmap", hashbits__validate_subset_partitionmap, METH_VARARGS, "" },
    { "set_partition_id", hashbits_set_partition_id, METH_VARARGS, "" },
    { "join_partitions", hashbits_join_partitions, METH_VARARGS, "" },
    { "get_partition_id", hashbits_get_partition_id, METH_VARARGS, "" },
    { "is_single_partition", hashbits_is_single_partition, METH_VARARGS, "" },
    { "count_kmers_within_radius", hashbits_count_kmers_within_radius, METH_VARARGS, "" },
    { "traverse_from_tags", hashbits_traverse_from_tags, METH_VARARGS, "" },
    { "repartition_largest_partition", hashbits_repartition_largest_partition, METH_VARARGS, "" },
    { "get_median_count", hashbits_get_median_count, METH_VARARGS, "Get the median, average, and stddev of the k-mer counts in the string" },
    {NULL, NULL, 0, NULL}           /* sentinel */
};

static PyObject * khmer_hashbits_getattr(PyObject * obj, char * name);
static PyObject* khmer_hashbits_new(PyTypeObject * type, 
        PyObject * args, PyObject * kwds);
static int khmer_hashbits_init(khmer_KHashbitsObject * self, PyObject * args,
        PyObject * kwds);
#define is_hashbits_obj(v)  ((v)->ob_type == &khmer_KHashbitsType)

////////////////////////////////////////////////////////////////////////////

static PyObject * subset_count_partitions(PyObject * self,
        PyObject * args);
static PyObject * subset_report_on_partitions(PyObject * self,
        PyObject * args);
static PyObject * subset_compare_partitions(PyObject * self,
        PyObject * args);
static PyObject * subset_partition_size_distribution(PyObject * self,
        PyObject * args);
static PyObject * subset_partition_sizes(PyObject * self,
        PyObject * args);
static PyObject * subset_partition_average_coverages(PyObject * self,
        PyObject * args);

static PyMethodDef khmer_subset_methods[] = {
    { "count_partitions", subset_count_partitions, METH_VARARGS, "" },
    { "report_on_partitions", subset_report_on_partitions, METH_VARARGS, "" },
    { "compare_partitions", subset_compare_partitions, METH_VARARGS, "" },
    { "partition_size_distribution", subset_partition_size_distribution, METH_VARARGS, "" },
    { "partition_sizes", subset_partition_sizes, METH_VARARGS, "" },
    { "partition_average_coverages", subset_partition_average_coverages, METH_VARARGS, "" },
    {NULL, NULL, 0, NULL}           /* sentinel */
};

static PyObject *
khmer_subset_getattr(PyObject * obj, char * name);

/////////////////
// LabelHash
/////////////////

// LabelHash addition
typedef struct {
    //PyObject_HEAD
    khmer_KHashbitsObject khashbits;
    LabelHash * labelhash;
} khmer_KLabelHashObject;

static void khmer_labelhash_dealloc(PyObject *);
static int khmer_labelhash_init(khmer_KLabelHashObject * self, PyObject *args,
                                PyObject *kwds);
static PyObject * khmer_labelhash_new(PyTypeObject * type, PyObject *args,
                                      PyObject *kwds);

#define is_labelhash_obj(v)  ((v)->ob_type == &khmer_KLabelHashType)

//
// khmer_labelhash_dealloc -- clean up a labelhash object.
//

static void khmer_labelhash_dealloc(PyObject* obj);
static PyObject * khmer_labelhash_new(PyTypeObject *type, PyObject *args,
                                      PyObject *kwds);
static PyObject * labelhash_get_label_dict(PyObject * self, PyObject * args);
static PyObject * labelhash_consume_fasta_and_tag_with_labels(
    PyObject * self, PyObject * args);
static PyObject * labelhash_consume_partitioned_fasta_and_tag_with_labels(
    PyObject * self, PyObject * args);
static PyObject * labelhash_consume_sequence_and_tag_with_labels(
    PyObject * self, PyObject * args);
static PyObject * labelhash_sweep_label_neighborhood(PyObject * self,
        PyObject * args);
static PyObject * labelhash_sweep_tag_neighborhood(PyObject * self,
        PyObject *args);
static PyObject * labelhash_get_tag_labels(PyObject * self, PyObject * args);
static PyObject * labelhash_n_labels(PyObject * self, PyObject * args);

static PyMethodDef khmer_labelhash_methods[] = {
    { "consume_fasta_and_tag_with_labels", labelhash_consume_fasta_and_tag_with_labels, METH_VARARGS, "" },
    { "sweep_label_neighborhood", labelhash_sweep_label_neighborhood, METH_VARARGS, "" },
    {"consume_partitioned_fasta_and_tag_with_labels", labelhash_consume_partitioned_fasta_and_tag_with_labels, METH_VARARGS, "" },
    {"sweep_tag_neighborhood", labelhash_sweep_tag_neighborhood, METH_VARARGS, "" },
    {"get_tag_labels", labelhash_get_tag_labels, METH_VARARGS, ""},
    {"consume_sequence_and_tag_with_labels", labelhash_consume_sequence_and_tag_with_labels, METH_VARARGS, "" },
    {"n_labels", labelhash_n_labels, METH_VARARGS, ""},
    {"get_label_dict", labelhash_get_label_dict, METH_VARARGS, "" },

    {NULL, NULL, 0, NULL}           /* sentinel */
};

static PyTypeObject khmer_KLabelHashType = {
    PyObject_HEAD_INIT(NULL)
    0,                       /* ob_size */
    "_LabelHash",            /* tp_name */
    sizeof(khmer_KLabelHashObject), /* tp_basicsize */
    0,                       /* tp_itemsize */
    (destructor)khmer_labelhash_dealloc, /* tp_dealloc */
    0,                       /* tp_print */
    0,  /* khmer_labelhash_getattr, tp_getattr */
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
    khmer_labelhash_methods, /* tp_methods */
    0,                       /* tp_members */
    0,                       /* tp_getset */
    0,                       /* tp_base */
    0,                       /* tp_dict */
    0,                       /* tp_descr_get */
    0,                       /* tp_descr_set */
    0,                       /* tp_dictoffset */
    (initproc)khmer_labelhash_init,   /* tp_init */
    0,                       /* tp_alloc */
};

///////////////
// ReadAligner
///////////////

static PyObject * readaligner_align(PyObject * self, PyObject * args);

static PyMethodDef khmer_ReadAligner_methods[] = {
    {"align", readaligner_align, METH_VARARGS, ""},
    {NULL, NULL, 0, NULL}
};

static PyObject *
khmer_readaligner_getattr(PyObject * obj, char * name);
static void khmer_readaligner_dealloc(PyObject* self);

static PyTypeObject khmer_ReadAlignerType = {
    PyObject_HEAD_INIT(NULL)
    0,
    "ReadAligner", sizeof(khmer_ReadAlignerObject),
    0,
    khmer_readaligner_dealloc,     /*tp_dealloc*/
    0,                          /*tp_print*/
    khmer_readaligner_getattr,     /*tp_getattr*/
    0,                          /*tp_setattr*/
    0,                          /*tp_compare*/
    0,                          /*tp_repr*/
    0,                          /*tp_as_number*/
    0,                          /*tp_as_sequence*/
    0,                          /*tp_as_mapping*/
    0,                          /*tp_hash */
    0,                          /*tp_call*/
    0,                          /*tp_str*/
    0,                          /*tp_getattro*/
    0,                          /*tp_setattro*/
    0,                          /*tp_as_buffer*/
    Py_TPFLAGS_DEFAULT,         /*tp_flags*/
    "ReadAligner object",           /* tp_doc */
};

static PyObject* new_readaligner(PyObject * self, PyObject * args);
static PyObject* _new_hashbits(PyObject * self, PyObject * args);
static PyObject * hash_collect_high_abundance_kmers(PyObject * self,
        PyObject * args);
static void khmer_hashbits_dealloc(PyObject* obj);
static void khmer_subset_dealloc(PyObject* self);

//////////////////////
// module methods
//////////////////////

static PyObject * forward_hash(PyObject * self, PyObject * args);
static PyObject * forward_hash_no_rc(PyObject * self, PyObject * args);
static PyObject * reverse_hash(PyObject * self, PyObject * args);
static PyObject * set_reporting_callback(PyObject * self, PyObject * args);

static
PyObject *
get_version_cpp( PyObject * self, PyObject * args );

static PyMethodDef KhmerMethods[] = {
#if (0)
    {
        "new_config",       new_config,
        METH_VARARGS,       "Create a default internals config"
    },
#endif
#if (0)
    {
        "set_config",       set_active_config,
        METH_VARARGS,       "Set active khmer configuration object"
    },
#endif
    {
        "new_hashtable",        new_hashtable,
        METH_VARARGS,       "Create an empty single-table counting hash"
    },
    {
        "_new_counting_hash",   _new_counting_hash,
        METH_VARARGS,       "Create an empty counting hash"
    },
    {
        "_new_hashbits",        _new_hashbits,
        METH_VARARGS,       "Create an empty hashbits table"
    },
    {
        "new_readaligner",        new_readaligner,
        METH_VARARGS,             "Create a read aligner object"
    },
    {
        "forward_hash",     forward_hash,
        METH_VARARGS,       "",
    },
    {
        "forward_hash_no_rc",   forward_hash_no_rc,
        METH_VARARGS,       "",
    },
    {
        "reverse_hash",     reverse_hash,
        METH_VARARGS,       "",
    },
    {
        "set_reporting_callback",   set_reporting_callback,
        METH_VARARGS,       "",
    },
    {
        "get_version_cpp", get_version_cpp,
        METH_VARARGS, "return the VERSION c++ compiler option"
    },
    { NULL, NULL, 0, NULL } // sentinel
};

#endif

// vim: set ft=cpp sts=4 sw=4 tw=79:
