#include <Python.h>

#include <iostream>

#include "khmer.hh"
#include "kmer_hash.hh"
#include "hashtable.hh"
#include "hashbits.hh"
#include "counting.hh"
#include "assembler.hh"
#include "read_aligner.hh"
#include "labelhash.hh"
#include "khmer_exception.hh"
#include "hllcounter.hh"


namespace khmer {
namespace python {


typedef struct {
    PyObject_HEAD
    //! Pointer to the low-level genomic read object.
    read_parsers:: Read *   read;
} khmer_Read_Object;


typedef struct {
    PyObject_HEAD
    //! Pointer to the low-level parser object.
    read_parsers:: IParser *  parser;
} khmer_ReadParser_Object;


typedef struct {
    PyObject_HEAD
    //! Pointer to Python parser object for reference counting purposes.
    PyObject *  parent;
    //! Persistent value of pair mode across invocations.
    int pair_mode;
} khmer_ReadPairIterator_Object;

}; //python

typedef struct {
    PyObject_HEAD
    pre_partition_info *   PrePartitionInfo;
} khmer_PrePartitionInfo_Object;



typedef struct {
    PyObject_HEAD
    SeenSet * hashes;
    WordLength ksize;
} khmer_HashSet_Object;


typedef struct {
    PyObject_HEAD
    khmer_HashSet_Object * parent;
    SeenSet::iterator * it;
} _HashSet_iterobj;


typedef struct {
    PyObject_HEAD
    Hashtable * hashtable;
} khmer_KHashtable_Object;


typedef struct {
    PyObject_HEAD
    SubsetPartition * subset;
} khmer_KSubsetPartition_Object;


typedef struct {
    khmer_KHashtable_Object khashtable;
    Hashbits * hashbits;
} khmer_KHashbits_Object;


typedef struct {
    khmer_KHashtable_Object khashtable;
    CountingHash * counting;
} khmer_KCountingHash_Object;

typedef struct {
    PyObject_HEAD
    ReadAligner * aligner;
} khmer_ReadAligner_Object;


typedef struct {
    PyObject_HEAD
    LabelHash * labelhash;
} khmer_KGraphLabels_Object;


typedef struct {
    PyObject_HEAD
    HLLCounter * hllcounter;
} khmer_KHLLCounter_Object;


typedef struct {
    PyObject_HEAD
    LinearAssembler * assembler;
} khmer_KLinearAssembler_Object;


typedef struct {
    PyObject_HEAD
    SimpleLabeledAssembler * assembler;
} khmer_KSimpleLabeledAssembler_Object;


typedef struct {
    PyObject_HEAD
    JunctionCountAssembler * assembler;
} khmer_KJunctionCountAssembler_Object;



};
