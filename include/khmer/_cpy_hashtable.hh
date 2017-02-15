#ifndef _CPY_HASHTABLE_HH
#define _CPY_HASHTABLE_HH

#include <Python.h>
#include "_cpy_utils.hh"
#include "oxli/oxli.hh"
#include "oxli/hashtable.hh"


namespace khmer {

typedef struct {
    PyObject_HEAD
    oxli::Hashtable * hashtable;
} khmer_KHashtable_Object;

extern PyMethodDef khmer_hashtable_methods[];

extern PyTypeObject khmer_KHashtable_Type
CPYCHECKER_TYPE_OBJECT_FOR_TYPEDEF("khmer_KHashtable_Object");


PyObject *
hashtable_ksize(khmer_KHashtable_Object * me, PyObject * args);


PyObject *
hashtable_hash(khmer_KHashtable_Object * me, PyObject * args);


PyObject *
hashtable_reverse_hash(khmer_KHashtable_Object * me, PyObject * args);


PyObject *
hashtable_n_occupied(khmer_KHashtable_Object * me, PyObject * args);


PyObject *
hashtable_n_unique_kmers(khmer_KHashtable_Object * me, PyObject * args);


PyObject *
hashtable_count(khmer_KHashtable_Object * me, PyObject * args);


PyObject *
hashtable_consume_fasta(khmer_KHashtable_Object * me, PyObject * args);


PyObject *
hashtable_consume_fasta_with_reads_parser(khmer_KHashtable_Object * me,
        PyObject * args);


PyObject *
hashtable_consume(khmer_KHashtable_Object * me, PyObject * args);


PyObject *
hashtable_get(khmer_KHashtable_Object * me, PyObject * args);


PyObject *
hashtable_set_use_bigcount(khmer_KHashtable_Object * me, PyObject * args);


PyObject *
hashtable_get_use_bigcount(khmer_KHashtable_Object * me, PyObject * args);


PyObject *
hashtable_get_min_count(khmer_KHashtable_Object * me, PyObject * args);


PyObject *
hashtable_get_max_count(khmer_KHashtable_Object * me, PyObject * args);


PyObject *
hashtable_abundance_distribution_with_reads_parser(khmer_KHashtable_Object * me,
        PyObject * args);


PyObject *
hashtable_trim_on_abundance(khmer_KHashtable_Object * me, PyObject * args);


PyObject *
hashtable_trim_below_abundance(khmer_KHashtable_Object * me, PyObject * args);


PyObject *
hashtable_find_spectral_error_positions(khmer_KHashtable_Object * me,
                                        PyObject * args);


PyObject *
hashtable_abundance_distribution(khmer_KHashtable_Object * me, PyObject * args);


PyObject *
hashtable_load(khmer_KHashtable_Object * me, PyObject * args);


PyObject *
hashtable_save(khmer_KHashtable_Object * me, PyObject * args);


PyObject *
hashtable_get_hashsizes(khmer_KHashtable_Object * me, PyObject * args);


PyObject *
hashtable_get_median_count(khmer_KHashtable_Object * me, PyObject * args);


PyObject *
hashtable_median_at_least(khmer_KHashtable_Object * me, PyObject * args);


PyObject *
hashtable_get_kmers(khmer_KHashtable_Object * me, PyObject * args);


PyObject *
hashtable_get_kmer_counts(khmer_KHashtable_Object * me, PyObject * args);



PyObject *
hashtable_get_kmer_hashes(khmer_KHashtable_Object * me, PyObject * args);



PyObject *
hashtable_get_kmer_hashes_as_hashset(khmer_KHashtable_Object * me,
                                     PyObject * args);


}

#endif
