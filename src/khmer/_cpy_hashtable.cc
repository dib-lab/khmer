#include "khmer/_cpy_hashtable.hh"
#include "khmer/_cpy_nodegraph.hh"
#include "khmer/_cpy_countgraph.hh"
#include "khmer/_cpy_hashset.hh"
#include "khmer/_cpy_readparsers.hh"

using namespace oxli;
using namespace oxli::read_parsers;

namespace khmer {

PyTypeObject khmer_KHashtable_Type
CPYCHECKER_TYPE_OBJECT_FOR_TYPEDEF("khmer_KHashtable_Object")
= {
    PyVarObject_HEAD_INIT(NULL, 0)       /* init & ob_size */
    "_khmer.KHashtable   ",              /*tp_name*/
    sizeof(khmer_KHashtable_Object),     /*tp_basicsize*/
    0,                                   /*tp_itemsize*/
    0,                                   /*tp_dealloc*/
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
    Py_TPFLAGS_DEFAULT,                  /*tp_flags*/
    "base hashtable object",             /* tp_doc */
    0,                                   /* tp_traverse */
    0,                                   /* tp_clear */
    0,                                   /* tp_richcompare */
    0,                                   /* tp_weaklistoffset */
    0,                                   /* tp_iter */
    0,                                   /* tp_iternext */
    khmer_hashtable_methods,             /* tp_methods */
    0,                                   /* tp_members */
    0,                                   /* tp_getset */
    0,                                   /* tp_base */
    0,                                   /* tp_dict */
    0,                                   /* tp_descr_get */
    0,                                   /* tp_descr_set */
    0,                                   /* tp_dictoffset */
    0,                                   /* tp_init */
    0,                                   /* tp_alloc */
    0,                                   /* tp_new */
};


PyMethodDef khmer_hashtable_methods[] = {
    //
    // Basic methods
    //

    {
        "ksize",
        (PyCFunction)hashtable_ksize, METH_VARARGS,
        "Returns the k-mer size of this graph."
    },
    {
        "hash",
        (PyCFunction)hashtable_hash, METH_VARARGS,
        "Returns the hash of this k-mer. For Nodetables and Counttables, this "
        "function will fail if the supplied k-mer contains non-ACGT "
        "characters."
    },
    {
        "reverse_hash",
        (PyCFunction)hashtable_reverse_hash, METH_VARARGS,
        "Turns a k-mer hash back into a DNA k-mer, if possible."
    },
    {
        "hashsizes",
        (PyCFunction)hashtable_get_hashsizes, METH_VARARGS,
        "" },
    {
        "n_unique_kmers",
        (PyCFunction)hashtable_n_unique_kmers, METH_VARARGS,
        "Count the number of unique kmers in this graph."
    },
    {
        "n_occupied", (PyCFunction)hashtable_n_occupied, METH_VARARGS,
        "Count the number of occupied bins."
    },
    {
        "count",
        (PyCFunction)hashtable_count, METH_VARARGS,
        "Increment the count of this k-mer."
    },
    {
        "add",
        (PyCFunction)hashtable_count, METH_VARARGS,
        "Increment the count of this k-mer. (Synonym for 'count'.)"
    },
    {
        "consume",
        (PyCFunction)hashtable_consume, METH_VARARGS,
        "Increment the counts of all of the k-mers in the string."
    },
    {
        "consume_seqfile",
        (PyCFunction)hashtable_consume_seqfile, METH_VARARGS,
        "Increment the counts of all the k-mers in the sequences in the "
        "given file"
    },
    {
        "consume_seqfile_banding",
        (PyCFunction)hashtable_consume_seqfile_banding, METH_VARARGS,
        "Consume sequences in k-mer banding mode"
    },
    {
        "consume_seqfile_with_mask",
        (PyCFunction)hashtable_consume_seqfile_with_mask, METH_VARARGS,
        "Consume any k-mers not present in the provided mask"
    },
    {
        "consume_seqfile_banding_with_mask",
        (PyCFunction)hashtable_consume_seqfile_banding_with_mask, METH_VARARGS,
        "Consume sequences in k-mer banding mode, with a mask"
    },
    {
        "consume_seqfile_with_reads_parser",
        (PyCFunction)hashtable_consume_seqfile_with_reads_parser, METH_VARARGS,
        "Count all k-mers retrieved with this reads parser object."
    },
    {
        "get",
        (PyCFunction)hashtable_get, METH_VARARGS,
        "Retrieve the count for the given k-mer. For Nodetables and "
        "Counttables, this function will fail if the supplied k-mer contains "
        "non-ACGT characters."
    },
    {
        "load",
        (PyCFunction)hashtable_load, METH_VARARGS,
        "Load the graph from the specified file."
    },
    {
        "save",
        (PyCFunction)hashtable_save, METH_VARARGS,
        "Save the graph to the specified file."
    },
    {
        "get_kmers",
        (PyCFunction)hashtable_get_kmers, METH_VARARGS,
        "Generate an ordered list of all substrings of length k in the string."
    },
    {
        "get_kmer_hashes",
        (PyCFunction)hashtable_get_kmer_hashes, METH_VARARGS,
        "Retrieve an ordered list of all hashes of all k-mers in the string."
    },
    {
        "get_kmer_hashes_as_hashset",
        (PyCFunction)hashtable_get_kmer_hashes_as_hashset, METH_VARARGS,
        "Retrieve a HashSet containing all the k-mers in the string."
    },
    {
        "get_kmer_counts",
        (PyCFunction)hashtable_get_kmer_counts, METH_VARARGS,
        "Retrieve an ordered list of the counts of all k-mers in the string."
    },

    {
        "set_use_bigcount",
        (PyCFunction)hashtable_set_use_bigcount, METH_VARARGS,
        "Count past maximum binsize of hashtable (set to T/F)"
    },
    {
        "get_use_bigcount",
        (PyCFunction)hashtable_get_use_bigcount, METH_VARARGS,
        "Get value of bigcount flag (T/F)"
    },
    {
        "get_min_count",
        (PyCFunction)hashtable_get_min_count, METH_VARARGS,
        "Get the smallest count of all the k-mers in the string"
    },
    {
        "get_max_count",
        (PyCFunction)hashtable_get_max_count, METH_VARARGS,
        "Get the largest count of all the k-mers in the string"
    },
    {
        "trim_on_abundance",
        (PyCFunction)hashtable_trim_on_abundance, METH_VARARGS,
        "Trim string at first k-mer below the given abundance"
    },
    {
        "trim_below_abundance",
        (PyCFunction)hashtable_trim_below_abundance, METH_VARARGS,
        "Trim string at first k-mer above the given abundance"
    },
    {
        "find_spectral_error_positions",
        (PyCFunction)hashtable_find_spectral_error_positions, METH_VARARGS,
        "Identify positions of low-abundance k-mers"
    },
    {
        "abundance_distribution",
        (PyCFunction)hashtable_abundance_distribution, METH_VARARGS,
        "Calculate the k-mer abundance distribution of the given file"
    },
    {
        "abundance_distribution_with_reads_parser",
        (PyCFunction)hashtable_abundance_distribution_with_reads_parser,
        METH_VARARGS,
        "Calculate the k-mer abundance distribution for a reads parser handle"
    },
    {
        "get_median_count",
        (PyCFunction)hashtable_get_median_count, METH_VARARGS,
        "Get the median, average, and stddev of the k-mer counts in the string"
    },
    {
        "median_at_least",
        (PyCFunction)hashtable_median_at_least, METH_VARARGS,
        "Return true if the median is at least the given cutoff"
    },
    {NULL, NULL, 0, NULL}           /* sentinel */
};



PyObject *
hashtable_ksize(khmer_KHashtable_Object * me, PyObject * args)
{
    Hashtable * hashtable = me->hashtable;

    if (!PyArg_ParseTuple(args, "")) {
        return NULL;
    }

    unsigned int k = hashtable->ksize();

    return PyLong_FromLong(k);
}


PyObject *
hashtable_hash(khmer_KHashtable_Object * me, PyObject * args)
{
    Hashtable * hashtable = me->hashtable;

    char * kmer;
    if (!PyArg_ParseTuple(args, "s", &kmer)) {
        return NULL;
    }

    if (strlen(kmer) != hashtable->ksize()) {
        PyErr_SetString(PyExc_ValueError,
                        "provided k-mer is wrong length");
        return NULL;
    }

    try {
        PyObject * hash = nullptr;
        const HashIntoType h(hashtable->hash_dna(kmer));
        convert_HashIntoType_to_PyObject(h, &hash);
        return hash;
    } catch (oxli_exception &e) {
        PyErr_SetString(PyExc_ValueError, e.what());
        return NULL;
    }
}


PyObject *
hashtable_reverse_hash(khmer_KHashtable_Object * me, PyObject * args)
{
    Hashtable * hashtable = me->hashtable;

    PyObject * val_o;
    HashIntoType val;
    if (!PyArg_ParseTuple(args, "O", &val_o)) {
        return NULL;
    }

    if (!ht_convert_PyObject_to_HashIntoType(val_o, val, hashtable)) {
        return NULL;
    }

    try {
        return PyUnicode_FromString(hashtable->unhash_dna(val).c_str());
    } catch (oxli_exception &e) {
        PyErr_SetString(PyExc_ValueError, e.what());
        return NULL;
    }
}


PyObject *
hashtable_n_occupied(khmer_KHashtable_Object * me, PyObject * args)
{
    Hashtable * hashtable = me->hashtable;

    if (!PyArg_ParseTuple(args, "")) {
        return NULL;
    }

    uint64_t n = hashtable->n_occupied();

    return PyLong_FromUnsignedLongLong(n);
}


PyObject *
hashtable_n_unique_kmers(khmer_KHashtable_Object * me, PyObject * args)
{
    Hashtable * hashtable = me->hashtable;

    uint64_t n = hashtable->n_unique_kmers();

    return PyLong_FromUnsignedLongLong(n);
}


PyObject *
hashtable_count(khmer_KHashtable_Object * me, PyObject * args)
{
    Hashtable * hashtable = me->hashtable;

    PyObject * v;
    if (!PyArg_ParseTuple(args, "O", &v)) {
        return NULL;
    }

    HashIntoType hashval;

    if (!ht_convert_PyObject_to_HashIntoType(v, hashval, hashtable)) {
        return NULL;
    }

    hashtable->count(hashval);

    return PyLong_FromLong(1);
}

PyObject *
hashtable_consume_seqfile(khmer_KHashtable_Object * me, PyObject * args)
{
    Hashtable * hashtable  = me->hashtable;

    const char * filename;

    if (!PyArg_ParseTuple(args, "s", &filename)) {
        return NULL;
    }

    // call the C++ function, and trap signals => Python
    unsigned long long  n_consumed    = 0;
    unsigned int          total_reads   = 0;
    try {
        hashtable->consume_seqfile<FastxReader>(filename, total_reads, n_consumed);
    } catch (oxli_file_exception &exc) {
        PyErr_SetString(PyExc_OSError, exc.what());
        return NULL;
    } catch (oxli_value_exception &exc) {
        PyErr_SetString(PyExc_ValueError, exc.what());
        return NULL;
    }

    return Py_BuildValue("IK", total_reads, n_consumed);
}

PyObject *
hashtable_consume_seqfile_banding(khmer_KHashtable_Object * me, PyObject * args)
{
    Hashtable * hashtable  = me->hashtable;

    const char * filename;
    unsigned int num_bands;
    unsigned int band;

    if (!PyArg_ParseTuple(args, "sII", &filename, &num_bands, &band)) {
        return NULL;
    }

    // call the C++ function, and trap signals => Python
    unsigned long long  n_consumed    = 0;
    unsigned int          total_reads   = 0;
    try {
        hashtable->consume_seqfile_banding<FastxReader>(filename, num_bands, band, total_reads, n_consumed);
    } catch (oxli_file_exception &exc) {
        PyErr_SetString(PyExc_OSError, exc.what());
        return NULL;
    } catch (oxli_value_exception &exc) {
        PyErr_SetString(PyExc_ValueError, exc.what());
        return NULL;
    }

    return Py_BuildValue("IK", total_reads, n_consumed);
}

PyObject *
hashtable_consume_seqfile_with_mask(khmer_KHashtable_Object * me, PyObject * args)
{
    Hashtable * hashtable  = me->hashtable;

    const char * filename;
    khmer_KHashtable_Object *mask = NULL;
    unsigned int threshold = 0;

    if (!PyArg_ParseTuple(args, "sO|I", &filename, &mask, &threshold)) {
        return NULL;
    }

    // call the C++ function, and trap signals => Python
    unsigned long long n_consumed = 0;
    unsigned int total_reads = 0;
    try {
        hashtable->consume_seqfile_with_mask<FastxReader>(
            filename, mask->hashtable, threshold, total_reads, n_consumed
        );
    } catch (oxli_file_exception &exc) {
        PyErr_SetString(PyExc_OSError, exc.what());
        return NULL;
    } catch (oxli_value_exception &exc) {
        PyErr_SetString(PyExc_ValueError, exc.what());
        return NULL;
    }

    return Py_BuildValue("IK", total_reads, n_consumed);
}

PyObject *
hashtable_consume_seqfile_banding_with_mask(khmer_KHashtable_Object * me, PyObject * args)
{
    Hashtable * hashtable  = me->hashtable;

    const char * filename;
    unsigned int num_bands;
    unsigned int band;
    khmer_KHashtable_Object *mask = NULL;
    unsigned int threshold = 0;

    if (!PyArg_ParseTuple(args, "sIIO|I", &filename, &num_bands, &band, &mask, &threshold)) {
        return NULL;
    }

    // call the C++ function, and trap signals => Python
    unsigned long long n_consumed = 0;
    unsigned int total_reads = 0;
    try {
        hashtable->consume_seqfile_banding_with_mask<FastxReader>(
            filename, num_bands, band, mask->hashtable, threshold, total_reads,
            n_consumed
        );
    } catch (oxli_file_exception &exc) {
        PyErr_SetString(PyExc_OSError, exc.what());
        return NULL;
    } catch (oxli_value_exception &exc) {
        PyErr_SetString(PyExc_ValueError, exc.what());
        return NULL;
    }

    return Py_BuildValue("IK", total_reads, n_consumed);
}

PyObject *
hashtable_consume_seqfile_with_reads_parser(khmer_KHashtable_Object * me,
        PyObject * args)
{
    Hashtable * hashtable = me->hashtable;

    PyObject * rparser_obj = NULL;

    if (!PyArg_ParseTuple(args, "O", &rparser_obj)) {
        return NULL;
    }

    FastxParserPtr& rparser = _PyObject_to_khmer_ReadParser( rparser_obj );

    // call the C++ function, and trap signals => Python
    unsigned long long  n_consumed      = 0;
    unsigned int        total_reads     = 0;
    const char         *value_exception = NULL;
    const char         *file_exception  = NULL;
    std::string exc_string;

    Py_BEGIN_ALLOW_THREADS
    try {
        hashtable->consume_seqfile<FastxReader>(rparser, total_reads, n_consumed);
    } catch (oxli_file_exception &exc) {
        exc_string = exc.what();
        file_exception = exc_string.c_str();
    } catch (oxli_value_exception &exc) {
        exc_string = exc.what();
        value_exception = exc_string.c_str();
    }
    Py_END_ALLOW_THREADS

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
hashtable_consume(khmer_KHashtable_Object * me, PyObject * args)
{
    Hashtable * hashtable = me->hashtable;

    const char * long_str;

    if (!PyArg_ParseTuple(args, "s", &long_str)) {
        return NULL;
    }

    if (strlen(long_str) < hashtable->ksize()) {
        PyErr_SetString(PyExc_ValueError,
                        "string length must >= the hashtable k-mer size");
        return NULL;
    }

    unsigned int n_consumed;
    n_consumed = hashtable->consume_string(long_str);

    return PyLong_FromLong(n_consumed);
}


PyObject *
hashtable_get(khmer_KHashtable_Object * me, PyObject * args)
{
    Hashtable * hashtable = me->hashtable;

    PyObject * arg;

    if (!PyArg_ParseTuple(args, "O", &arg)) {
        return NULL;
    }

    HashIntoType hashval;

    if (!ht_convert_PyObject_to_HashIntoType(arg, hashval, hashtable)) {
        return NULL;
    }

    unsigned int count = hashtable->get_count(hashval);
    return PyLong_FromLong(count);
}


PyObject *
hashtable_set_use_bigcount(khmer_KHashtable_Object * me, PyObject * args)
{
    Hashtable * hashtable = me->hashtable;

    PyObject * x;
    if (!PyArg_ParseTuple(args, "O", &x)) {
        return NULL;
    }
    int setme = PyObject_IsTrue(x);
    if (setme < 0) {
        return NULL;
    }
    try {
        hashtable->set_use_bigcount((bool)setme);
    } catch (oxli_exception &e) {
        PyErr_SetString(PyExc_ValueError, e.what());
        return NULL;
    }

    Py_RETURN_NONE;
}


PyObject *
hashtable_get_use_bigcount(khmer_KHashtable_Object * me, PyObject * args)
{
    Hashtable * hashtable = me->hashtable;

    if (!PyArg_ParseTuple(args, "")) {
        return NULL;
    }

    bool val = hashtable->get_use_bigcount();

    return PyBool_FromLong((int)val);
}


PyObject *
hashtable_get_min_count(khmer_KHashtable_Object * me, PyObject * args)
{
    Hashtable * hashtable = me->hashtable;

    const char * long_str;

    if (!PyArg_ParseTuple(args, "s", &long_str)) {
        return NULL;
    }

    if (strlen(long_str) < hashtable->ksize()) {
        PyErr_SetString(PyExc_ValueError,
                        "string length must >= the hashtable k-mer size");
        return NULL;
    }

    BoundedCounterType c = hashtable->get_min_count(long_str);
    unsigned int N = c;

    return PyLong_FromLong(N);
}


PyObject *
hashtable_get_max_count(khmer_KHashtable_Object * me, PyObject * args)
{
    Hashtable * hashtable = me->hashtable;

    const char * long_str;

    if (!PyArg_ParseTuple(args, "s", &long_str)) {
        return NULL;
    }

    if (strlen(long_str) < hashtable->ksize()) {
        PyErr_SetString(PyExc_ValueError,
                        "string length must >= the hashtable k-mer size");
        return NULL;
    }

    BoundedCounterType c = hashtable->get_max_count(long_str);
    unsigned int N = c;

    return PyLong_FromLong(N);
}


PyObject *
hashtable_abundance_distribution_with_reads_parser(khmer_KHashtable_Object * me,
        PyObject * args)
{
    Hashtable * hashtable = me->hashtable;

    khmer :: khmer_ReadParser_Object * rparser_obj = NULL;
    khmer_KHashtable_Object * tracking_obj = NULL;

    if (!PyArg_ParseTuple(args, "O!O!", &khmer_ReadParser_Type,
                          &rparser_obj, &khmer_KHashtable_Type, &tracking_obj)) {
        return NULL;
    }

    FastxParserPtr& rparser = rparser_obj->parser;
    Hashtable          *tracking        = tracking_obj->hashtable;
    uint64_t           *dist            = NULL;
    const char         *value_exception = NULL;
    const char         *file_exception  = NULL;
    std::string exc_string;

    Py_BEGIN_ALLOW_THREADS
    try {
        dist = hashtable->abundance_distribution<FastxReader>(rparser, tracking);
    } catch (oxli_file_exception &exc) {
        exc_string = exc.what();
        file_exception = exc_string.c_str();
    } catch (oxli_value_exception &exc) {
        exc_string = exc.what();
        value_exception = exc_string.c_str();
    }
    Py_END_ALLOW_THREADS

    if (file_exception != NULL) {
        PyErr_SetString(PyExc_OSError, file_exception);
        return NULL;
    }
    if (value_exception != NULL) {
        PyErr_SetString(PyExc_ValueError, value_exception);
        return NULL;
    }

    PyObject * x = PyList_New(MAX_BIGCOUNT + 1);
    if (x == NULL) {
        delete[] dist;
        return NULL;
    }
    for (int i = 0; i < MAX_BIGCOUNT + 1; i++) {
        PyList_SET_ITEM(x, i, PyLong_FromUnsignedLongLong(dist[i]));
    }

    delete[] dist;
    return x;
}


PyObject *
hashtable_trim_on_abundance(khmer_KHashtable_Object * me, PyObject * args)
{
    Hashtable * hashtable = me->hashtable;

    const char * seq = NULL;
    unsigned int min_count_i = 0;

    if (!PyArg_ParseTuple(args, "sI", &seq, &min_count_i)) {
        return NULL;
    }

    unsigned long trim_at;
    Py_BEGIN_ALLOW_THREADS

    BoundedCounterType min_count = min_count_i;

    trim_at = hashtable->trim_on_abundance(seq, min_count);

    Py_END_ALLOW_THREADS;

    PyObject * trim_seq = PyUnicode_FromStringAndSize(seq, trim_at);
    if (trim_seq == NULL) {
        return NULL;
    }
    PyObject * ret = Py_BuildValue("Ok", trim_seq, trim_at);
    Py_DECREF(trim_seq);

    return ret;
}


PyObject *
hashtable_trim_below_abundance(khmer_KHashtable_Object * me, PyObject * args)
{
    Hashtable * hashtable = me->hashtable;

    const char * seq = NULL;
    BoundedCounterType max_count_i = 0;

    if (!PyArg_ParseTuple(args, "sH", &seq, &max_count_i)) {
        return NULL;
    }

    unsigned long trim_at;
    Py_BEGIN_ALLOW_THREADS

    BoundedCounterType max_count = max_count_i;

    trim_at = hashtable->trim_below_abundance(seq, max_count);

    Py_END_ALLOW_THREADS;

    PyObject * trim_seq = PyUnicode_FromStringAndSize(seq, trim_at);
    if (trim_seq == NULL) {
        return NULL;
    }
    PyObject * ret = Py_BuildValue("Ok", trim_seq, trim_at);
    Py_DECREF(trim_seq);

    return ret;
}


PyObject *
hashtable_find_spectral_error_positions(khmer_KHashtable_Object * me,
                                        PyObject * args)
{
    Hashtable * hashtable = me->hashtable;

    const char * seq = NULL;
    BoundedCounterType max_count = 0; // unsigned short int

    if (!PyArg_ParseTuple(args, "sH", &seq, &max_count)) {
        return NULL;
    }

    std::vector<unsigned int> posns;

    try {
        posns = hashtable->find_spectral_error_positions(seq, max_count);
    } catch (oxli_exception &e) {
        PyErr_SetString(PyExc_ValueError, e.what());
        return NULL;
    }

    Py_ssize_t posns_size = posns.size();

    PyObject * x = PyList_New(posns_size);
    if (x == NULL) {
        return NULL;
    }
    for (Py_ssize_t i = 0; i < posns_size; i++) {
        PyList_SET_ITEM(x, i, PyLong_FromLong(posns[i]));
    }

    return x;
}


PyObject *
hashtable_abundance_distribution(khmer_KHashtable_Object * me, PyObject * args)
{
    Hashtable * hashtable = me->hashtable;

    const char * filename = NULL;
    khmer_KHashtable_Object * tracking_obj = NULL;
    if (!PyArg_ParseTuple(args, "sO!", &filename, &khmer_KHashtable_Type,
                          &tracking_obj)) {
        return NULL;
    }

    Hashtable          *tracking       = tracking_obj->hashtable;
    uint64_t           *dist            = NULL;
    const char         *value_exception = NULL;
    const char         *file_exception  = NULL;
    std::string exc_string;

    Py_BEGIN_ALLOW_THREADS
    try {
        dist = hashtable->abundance_distribution<FastxReader>(filename, tracking);
    } catch (oxli_file_exception &exc) {
        exc_string = exc.what();
        file_exception = exc_string.c_str();
    } catch (oxli_value_exception &exc) {
        exc_string = exc.what();
        value_exception = exc_string.c_str();
    }
    Py_END_ALLOW_THREADS

    if (file_exception != NULL) {
        PyErr_SetString(PyExc_OSError, file_exception);
        if (dist != NULL) {
            delete []dist;
        }
        return NULL;
    }
    if (value_exception != NULL) {
        PyErr_SetString(PyExc_ValueError, value_exception);
        if (dist != NULL) {
            delete []dist;
        }
        return NULL;
    }

    PyObject * x = PyList_New(MAX_BIGCOUNT + 1);
    if (x == NULL) {
        if (dist != NULL) {
            delete []dist;
        }
        return NULL;
    }
    for (int i = 0; i < MAX_BIGCOUNT + 1; i++) {
        PyList_SET_ITEM(x, i, PyLong_FromUnsignedLongLong(dist[i]));
    }

    if (dist != NULL) {
        delete []dist;
    }

    return x;
}


PyObject *
hashtable_load(khmer_KHashtable_Object * me, PyObject * args)
{
    Hashtable * hashtable = me->hashtable;

    const char * filename = NULL;

    if (!PyArg_ParseTuple(args, "s", &filename)) {
        return NULL;
    }

    try {
        hashtable->load(filename);
    } catch (oxli_file_exception &e) {
        PyErr_SetString(PyExc_OSError, e.what());
        return NULL;
    }

    Py_RETURN_NONE;
}


PyObject *
hashtable_save(khmer_KHashtable_Object * me, PyObject * args)
{
    Hashtable * hashtable = me->hashtable;

    const char * filename = NULL;

    if (!PyArg_ParseTuple(args, "s", &filename)) {
        return NULL;
    }

    try {
        hashtable->save(filename);
    } catch (oxli_file_exception &e) {
        PyErr_SetString(PyExc_OSError, e.what());
        return NULL;
    }

    Py_RETURN_NONE;
}


PyObject *
hashtable_get_hashsizes(khmer_KHashtable_Object * me, PyObject * args)
{
    Hashtable * hashtable = me->hashtable;


    if (!PyArg_ParseTuple(args, "")) {
        return NULL;
    }

    std::vector<uint64_t> ts = hashtable->get_tablesizes();

    PyObject * x = PyList_New(ts.size());
    for (size_t i = 0; i < ts.size(); i++) {
        PyList_SET_ITEM(x, i, PyLong_FromUnsignedLongLong(ts[i]));
    }

    return x;
}


PyObject *
hashtable_get_median_count(khmer_KHashtable_Object * me, PyObject * args)
{
    Hashtable * hashtable = me->hashtable;

    const char * long_str;

    if (!PyArg_ParseTuple(args, "s", &long_str)) {
        return NULL;
    }

    if (strlen(long_str) < hashtable->ksize()) {
        PyErr_SetString(PyExc_ValueError,
                        "string length must >= the hashtable k-mer size");
        return NULL;
    }

    BoundedCounterType med = 0;
    float average = 0, stddev = 0;

    hashtable->get_median_count(long_str, med, average, stddev);

    return Py_BuildValue("iff", med, average, stddev);
}


PyObject *
hashtable_median_at_least(khmer_KHashtable_Object * me, PyObject * args)
{
    Hashtable * hashtable = me->hashtable;

    const char * long_str;
    unsigned int cutoff;

    if (!PyArg_ParseTuple(args, "sI", &long_str, &cutoff)) {
        return NULL;
    }

    if (strlen(long_str) < hashtable->ksize()) {
        PyErr_SetString(PyExc_ValueError,
                        "string length must >= the hashtable k-mer size");
        return NULL;
    }

    if (hashtable->median_at_least(long_str, cutoff)) {
        Py_RETURN_TRUE;
    }
    Py_RETURN_FALSE;

}


PyObject *
hashtable_get_kmers(khmer_KHashtable_Object * me, PyObject * args)
{
    Hashtable * hashtable = me->hashtable;
    const char * sequence;

    if (!PyArg_ParseTuple(args, "s", &sequence)) {
        return NULL;
    }

    std::vector<std::string> kmers;

    hashtable->get_kmers(sequence, kmers);

    PyObject * x = PyList_New(kmers.size());
    for (unsigned int i = 0; i < kmers.size(); i++) {
        PyObject * obj = PyUnicode_FromString(kmers[i].c_str());
        PyList_SET_ITEM(x, i, obj);
    }

    return x;
}


PyObject *
hashtable_get_kmer_counts(khmer_KHashtable_Object * me, PyObject * args)
{
    Hashtable * hashtable = me->hashtable;
    const char * sequence;

    if (!PyArg_ParseTuple(args, "s", &sequence)) {
        return NULL;
    }

    std::vector<BoundedCounterType> counts;
    try {
        hashtable->get_kmer_counts(sequence, counts);
    } catch (oxli_exception &e) {
        PyErr_SetString(PyExc_ValueError, e.what());
        return NULL;
    }

    PyObject * x = PyList_New(counts.size());
    for (unsigned int i = 0; i <counts.size(); i++) {
        PyObject * obj = PyInt_FromLong(counts[i]);
        PyList_SET_ITEM(x, i, obj);
    }

    return x;
}



PyObject *
hashtable_get_kmer_hashes(khmer_KHashtable_Object * me, PyObject * args)
{
    Hashtable * hashtable = me->hashtable;
    const char * sequence;

    if (!PyArg_ParseTuple(args, "s", &sequence)) {
        return NULL;
    }

    std::vector<HashIntoType> hashes;
    try {
        hashtable->get_kmer_hashes(sequence, hashes);
    } catch (oxli_exception &e) {
        PyErr_SetString(PyExc_ValueError, e.what());
        return NULL;
    }

    PyObject * x = PyList_New(hashes.size());
    for (unsigned int i = 0; i < hashes.size(); i++) {
        PyObject * obj = nullptr;
        convert_HashIntoType_to_PyObject(hashes[i], &obj);
        PyList_SET_ITEM(x, i, obj);
    }

    return x;
}



PyObject *
hashtable_get_kmer_hashes_as_hashset(khmer_KHashtable_Object * me,
                                     PyObject * args)
{
    Hashtable * hashtable = me->hashtable;
    const char * sequence;

    if (!PyArg_ParseTuple(args, "s", &sequence)) {
        return NULL;
    }

    SeenSet * hashes = new SeenSet;
    try {
        hashtable->get_kmer_hashes_as_hashset(sequence, *hashes);
    } catch (oxli_exception &e) {
        PyErr_SetString(PyExc_ValueError, e.what());
        return NULL;
    }

    PyObject * x = (PyObject *) create_HashSet_Object(hashes,
                   hashtable->ksize());

    return x;
}

}
