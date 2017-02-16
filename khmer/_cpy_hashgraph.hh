/*
This file is part of khmer, https://github.com/dib-lab/khmer/, and is
Copyright (C) 2010-2015, Michigan State University.
Copyright (C) 2015-2016, The Regents of the University of California.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are
met:

    * Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.

    * Redistributions in binary form must reproduce the above
      copyright notice, this list of conditions and the following
      disclaimer in the documentation and/or other materials provided
      with the distribution.

    * Neither the name of the Michigan State University nor the names
      of its contributors may be used to endorse or promote products
      derived from this software without specific prior written
      permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
"AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
LICENSE (END)

Contact: khmer-project@idyll.org
*/

static PyTypeObject khmer_KHashgraph_Type
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

//
// Method definitions
//

static
PyObject *
hashgraph_find_high_degree_nodes(khmer_KHashgraph_Object * me, PyObject * args)
{
    Hashgraph * hashgraph = me->hashgraph;

    const char * long_str;

    if (!PyArg_ParseTuple(args, "s", &long_str)) {
        return NULL;
    }

    if (strlen(long_str) < hashgraph->ksize()) {
        PyErr_SetString(PyExc_ValueError,
                        "string length must >= the hashgraph k-mer size");
        return NULL;
    }

    SeenSet * hashes = new SeenSet;
    hashgraph->find_high_degree_nodes(long_str, *hashes);

    khmer_HashSet_Object * o;
    o = create_HashSet_Object(hashes, hashgraph->ksize());

    return (PyObject *) o;
}

static
PyObject *
hashgraph_neighbors(khmer_KHashgraph_Object * me, PyObject * args)
{
    Hashgraph * hashgraph = me->hashgraph;
    PyObject * val_obj;

    if (!PyArg_ParseTuple(args, "O", &val_obj)) {
        return NULL;
    }

    Kmer start_kmer;
    if (!ht_convert_PyObject_to_Kmer(val_obj, start_kmer, hashgraph)) {
        return NULL;
    }

    KmerQueue node_q;
    Traverser traverser(hashgraph);

    traverser.traverse(start_kmer, node_q);

    PyObject * x =  PyList_New(node_q.size());
    if (x == NULL) {
        return NULL;
    }

    unsigned int i;
    PyObject * value = nullptr;
    for (i = 0; node_q.size() > 0; i++) {
        const HashIntoType h = node_q.front();
        node_q.pop();
        convert_HashIntoType_to_PyObject(h, &value);
        PyList_SET_ITEM(x, i, value);
    }

    return x;
}

static
PyObject *
hashgraph_traverse_linear_path(khmer_KHashgraph_Object * me, PyObject * args)
{
    Hashgraph * hashgraph = me->hashgraph;

    PyObject * val_o;
    khmer_KNodegraph_Object * nodegraph_o = NULL;
    khmer_HashSet_Object * hdn_o = NULL;

    if (!PyArg_ParseTuple(args, "OO!O!", &val_o,
                          &khmer_HashSet_Type, &hdn_o,
                          &khmer_KNodegraph_Type, &nodegraph_o)) {
        return NULL;
    }
    Kmer start_kmer;
    if (!ht_convert_PyObject_to_Kmer(val_o, start_kmer, hashgraph)) {
        return NULL;
    }

    SeenSet * adj = new SeenSet;
    SeenSet * visited = new SeenSet;
    unsigned int size = hashgraph->traverse_linear_path(start_kmer,
                        *adj, *visited,
                        *nodegraph_o->nodegraph,
                        *hdn_o->hashes);

    khmer_HashSet_Object * adj_o = create_HashSet_Object(adj,
                                   hashgraph->ksize());
    khmer_HashSet_Object * visited_o = create_HashSet_Object(visited,
                                       hashgraph->ksize());

    PyObject * ret = Py_BuildValue("kOO", (unsigned long) size,
                                   (PyObject *) adj_o, (PyObject *) visited_o);
    Py_DECREF(adj_o);
    Py_DECREF(visited_o);

    return ret;
}

static
PyObject *
hashgraph_assemble_linear_path(khmer_KHashgraph_Object * me, PyObject * args)
{
    Hashgraph * hashgraph = me->hashgraph;

    PyObject * val_o;
    khmer_KNodegraph_Object * nodegraph_o = NULL;
    Nodegraph * stop_bf = NULL;

    if (!PyArg_ParseTuple(args, "O|O!", &val_o,
                          &khmer_KNodegraph_Type, &nodegraph_o)) {
        return NULL;
    }

    Kmer start_kmer;
    if (!ht_convert_PyObject_to_Kmer(val_o, start_kmer, hashgraph)) {
        return NULL;
    }

    if (nodegraph_o) {
        stop_bf = nodegraph_o->nodegraph;
    }
    LinearAssembler assembler(hashgraph);

    std::string contig = assembler.assemble(start_kmer, stop_bf);

    PyObject * ret = Py_BuildValue("s", contig.c_str());

    return ret;
}

static
PyObject *
hashgraph_n_tags(khmer_KHashgraph_Object * me, PyObject * args)
{
    Hashgraph * hashgraph = me->hashgraph;

    if (!PyArg_ParseTuple(args, "")) {
        return NULL;
    }

    return PyLong_FromSize_t(hashgraph->n_tags());
}

static
PyObject *
hashgraph_print_stop_tags(khmer_KHashgraph_Object * me, PyObject * args)
{
    Hashgraph * hashgraph = me->hashgraph;

    const char * filename = NULL;

    if (!PyArg_ParseTuple(args, "s", &filename)) {
        return NULL;
    }

    hashgraph->print_stop_tags(filename);

    Py_RETURN_NONE;
}

static
PyObject *
hashgraph_print_tagset(khmer_KHashgraph_Object * me, PyObject * args)
{
    Hashgraph * hashgraph = me->hashgraph;

    const char * filename = NULL;

    if (!PyArg_ParseTuple(args, "s", &filename)) {
        return NULL;
    }

    hashgraph->print_tagset(filename);

    Py_RETURN_NONE;
}

static
PyObject *
hashgraph_load_stop_tags(khmer_KHashgraph_Object * me, PyObject * args)
{
    Hashgraph * hashgraph = me->hashgraph;

    const char * filename = NULL;
    PyObject * clear_tags_o = NULL;

    if (!PyArg_ParseTuple(args, "s|O", &filename, &clear_tags_o)) {
        return NULL;
    }

    bool clear_tags = true;
    if (clear_tags_o && !PyObject_IsTrue(clear_tags_o)) {
        clear_tags = false;
    }


    try {
        hashgraph->load_stop_tags(filename, clear_tags);
    } catch (khmer_file_exception &e) {
        PyErr_SetString(PyExc_OSError, e.what());
        return NULL;
    }

    Py_RETURN_NONE;
}


static
PyObject *
hashgraph_save_stop_tags(khmer_KHashgraph_Object * me, PyObject * args)
{
    Hashgraph * hashgraph = me->hashgraph;

    const char * filename = NULL;

    if (!PyArg_ParseTuple(args, "s", &filename)) {
        return NULL;
    }

    try {
        hashgraph->save_stop_tags(filename);
    } catch (khmer_file_exception &e) {
        PyErr_SetString(PyExc_OSError, e.what());
        return NULL;
    }

    Py_RETURN_NONE;
}

static PyObject * hashgraph_repartition_largest_partition(
    khmer_KHashgraph_Object * me,
    PyObject * args);

static
PyObject *
hashgraph_calc_connected_graph_size(khmer_KHashgraph_Object * me,
                                    PyObject * args)
{
    Hashgraph * hashgraph = me->hashgraph;

    const char * _kmer;
    unsigned int max_size = 0;
    PyObject * break_on_circum_o = NULL;
    if (!PyArg_ParseTuple(args, "s|IO", &_kmer, &max_size, &break_on_circum_o)) {
        return NULL;
    }

    bool break_on_circum = false;
    if (break_on_circum_o && PyObject_IsTrue(break_on_circum_o)) {
        break_on_circum = true;
    }

    unsigned long long size = 0;
    Kmer start_kmer = hashgraph->build_kmer(_kmer);

    Py_BEGIN_ALLOW_THREADS
    KmerSet keeper;
    hashgraph->calc_connected_graph_size(start_kmer, size, keeper, max_size,
                                         break_on_circum);
    Py_END_ALLOW_THREADS

    return PyLong_FromUnsignedLongLong(size);
}

static
PyObject *
hashgraph_kmer_degree(khmer_KHashgraph_Object * me, PyObject * args)
{
    Hashgraph * hashgraph = me->hashgraph;

    const char * kmer_s = NULL;

    if (!PyArg_ParseTuple(args, "s", &kmer_s)) {
        return NULL;
    }

    return PyLong_FromLong(hashgraph->kmer_degree(kmer_s));
}

static
PyObject *
hashgraph_trim_on_stoptags(khmer_KHashgraph_Object * me, PyObject * args)
{
    Hashgraph * hashgraph = me->hashgraph;

    const char * seq = NULL;

    if (!PyArg_ParseTuple(args, "s", &seq)) {
        return NULL;
    }

    size_t trim_at;
    Py_BEGIN_ALLOW_THREADS

    trim_at = hashgraph->trim_on_stoptags(seq);

    Py_END_ALLOW_THREADS;

    PyObject * trim_seq = PyUnicode_FromStringAndSize(seq, trim_at);
    if (trim_seq == NULL) {
        return NULL;
    }
    PyObject * ret = Py_BuildValue("Ok", trim_seq, (unsigned long) trim_at);
    Py_DECREF(trim_seq);

    return ret;
}

static
PyObject *
hashgraph_do_subset_partition(khmer_KHashgraph_Object * me, PyObject * args)
{
    Hashgraph * hashgraph = me->hashgraph;

    PyObject * start_kmer_obj;
    PyObject * end_kmer_obj;
    HashIntoType start_kmer, end_kmer;
    PyObject * break_on_stop_tags_o = NULL;
    PyObject * stop_big_traversals_o = NULL;

    if (!PyArg_ParseTuple(args, "|OOOO", &start_kmer_obj, &end_kmer_obj,
                          &break_on_stop_tags_o,
                          &stop_big_traversals_o)) {
        return NULL;
    }
    if (!ht_convert_PyObject_to_HashIntoType(start_kmer_obj, start_kmer,
            hashgraph)) {
        return NULL;
    }
    if (!ht_convert_PyObject_to_HashIntoType(end_kmer_obj, end_kmer,
            hashgraph)) {
        return NULL;
    }

    bool break_on_stop_tags = false;
    if (break_on_stop_tags_o && PyObject_IsTrue(break_on_stop_tags_o)) {
        break_on_stop_tags = true;
    }
    bool stop_big_traversals = false;
    if (stop_big_traversals_o && PyObject_IsTrue(stop_big_traversals_o)) {
        stop_big_traversals = true;
    }

    SubsetPartition * subset_p = NULL;
    try {
        Py_BEGIN_ALLOW_THREADS
        subset_p = new SubsetPartition(hashgraph);
        subset_p->do_partition(start_kmer, end_kmer, break_on_stop_tags,
                               stop_big_traversals);
        Py_END_ALLOW_THREADS
    } catch (std::bad_alloc &e) {
        return PyErr_NoMemory();
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


static
PyObject *
hashgraph_merge_subset(khmer_KHashgraph_Object * me, PyObject * args)
{
    Hashgraph * hashgraph = me->hashgraph;

    khmer_KSubsetPartition_Object * subset_obj = NULL;
    if (!PyArg_ParseTuple(args, "O!", &khmer_KSubsetPartition_Type,
                          &subset_obj)) {
        return NULL;
    }

    SubsetPartition * subset_p = subset_obj->subset;

    hashgraph->partition->merge(subset_p);

    Py_RETURN_NONE;
}

static
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
    } catch (khmer_file_exception &e) {
        PyErr_SetString(PyExc_OSError, e.what());
        return NULL;
    }

    Py_RETURN_NONE;
}

static
PyObject *
hashgraph_consume_seqfile_and_tag_with_reads_parser(khmer_KHashgraph_Object * me,
        PyObject * args)
{
    Hashgraph * hashgraph = me->hashgraph;

    python::khmer_ReadParser_Object * rparser_obj = NULL;

    if (!PyArg_ParseTuple( args, "O!", &python::khmer_ReadParser_Type,
                           &rparser_obj)) {
        return NULL;
    }

    FastxParserPtr& rparser = rparser_obj->parser;

    // call the C++ function, and trap signals => Python
    const char         *value_exception = NULL;
    const char         *file_exception  = NULL;
    unsigned long long  n_consumed      = 0;
    unsigned int        total_reads     = 0;
    std::string exc_string;

    Py_BEGIN_ALLOW_THREADS
    try {
        hashgraph->consume_seqfile_and_tag<FastxReader>(rparser, total_reads, n_consumed);
    } catch (khmer_file_exception &exc) {
        exc_string = exc.what();
        file_exception = exc_string.c_str();
    } catch (khmer_value_exception &exc) {
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
    }

    return Py_BuildValue("IK", total_reads, n_consumed);
}

static
PyObject *
hashgraph_consume_partitioned_fasta(khmer_KHashgraph_Object * me,
                                    PyObject * args)
{
    Hashgraph * hashgraph = me->hashgraph;

    const char * filename;

    if (!PyArg_ParseTuple(args, "s", &filename)) {
        return NULL;
    }

    // call the C++ function, and trap signals => Python

    unsigned long long n_consumed;
    unsigned int total_reads;

    try {
        hashgraph->consume_partitioned_fasta<FastxReader>(filename, total_reads, n_consumed);
    } catch (khmer_file_exception &exc) {
        PyErr_SetString(PyExc_OSError, exc.what());
        return NULL;
    } catch (khmer_value_exception &exc) {
        PyErr_SetString(PyExc_ValueError, exc.what());
        return NULL;
    }

    return Py_BuildValue("IK", total_reads, n_consumed);
}

static
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

static
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

static
PyObject *
hashgraph_add_tag(khmer_KHashgraph_Object * me, PyObject * args)
{
    Hashgraph * hashgraph = me->hashgraph;

    const char * kmer_s = NULL;
    if (!PyArg_ParseTuple(args, "s", &kmer_s)) {
        return NULL;
    }

    HashIntoType kmer = hashgraph->hash_dna(kmer_s);
    hashgraph->add_tag(kmer);

    Py_RETURN_NONE;
}

static
PyObject *
hashgraph_add_stop_tag(khmer_KHashgraph_Object * me, PyObject * args)
{
    Hashgraph * hashgraph = me->hashgraph;

    const char * kmer_s = NULL;
    if (!PyArg_ParseTuple(args, "s", &kmer_s)) {
        return NULL;
    }

    HashIntoType kmer = hashgraph->hash_dna(kmer_s);
    hashgraph->add_stop_tag(kmer);

    Py_RETURN_NONE;
}

static
PyObject *
hashgraph_get_stop_tags(khmer_KHashgraph_Object * me, PyObject * args)
{
    Hashgraph * hashgraph = me->hashgraph;

    if (!PyArg_ParseTuple(args, "")) {
        return NULL;
    }

    SeenSet::const_iterator si;

    PyObject * x = PyList_New(hashgraph->stop_tags.size());
    unsigned long long i = 0;
    for (si = hashgraph->stop_tags.begin(); si != hashgraph->stop_tags.end();
            ++si) {
        std::string s = hashgraph->unhash_dna(*si);
        PyList_SET_ITEM(x, i, Py_BuildValue("s", s.c_str()));
        i++;
    }

    return x;
}

static
PyObject *
hashgraph_get_tagset(khmer_KHashgraph_Object * me, PyObject * args)
{
    Hashgraph * hashgraph = me->hashgraph;

    if (!PyArg_ParseTuple(args, "")) {
        return NULL;
    }

    SeenSet::const_iterator si;

    PyObject * x = PyList_New(hashgraph->all_tags.size());
    unsigned long long i = 0;
    for (si = hashgraph->all_tags.begin(); si != hashgraph->all_tags.end();
            ++si) {
        std::string s = hashgraph->unhash_dna(*si);
        PyList_SET_ITEM(x, i, Py_BuildValue("s", s.c_str()));
        i++;
    }

    return x;
}

static
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
    } catch (khmer_file_exception &e) {
        PyErr_SetString(PyExc_OSError, e.what());
        return NULL;
    } catch (khmer_value_exception &exc) {
        PyErr_SetString(PyExc_ValueError, exc.what());
        return NULL;
    }

    return PyLong_FromLong(n_partitions);
}

static
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
    } catch (khmer_file_exception &e) {
        PyErr_SetString(PyExc_OSError, e.what());
        return NULL;
    }

    Py_RETURN_NONE;
}

static
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
    } catch (khmer_file_exception &e) {
        PyErr_SetString(PyExc_OSError, e.what());
        return NULL;
    }

    Py_RETURN_NONE;
}

static
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

static
PyObject *
hashgraph_count_partitions(khmer_KHashgraph_Object * me, PyObject * args)
{
    Hashgraph * hashgraph = me->hashgraph;

    if (!PyArg_ParseTuple(args, "")) {
        return NULL;
    }

    size_t n_partitions = 0, n_unassigned = 0;
    hashgraph->partition->count_partitions(n_partitions, n_unassigned);

    return Py_BuildValue("nn", (Py_ssize_t) n_partitions,
                         (Py_ssize_t) n_unassigned);
}

static
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

static
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
        throw khmer_exception();
    }

    PyObject * returnValue = Py_BuildValue("NI", x, n_unassigned);
    if (returnValue == NULL) {
        Py_DECREF(x);
        return NULL;
    }
    return returnValue;
}

static
PyObject *
hashgraph_load_tagset(khmer_KHashgraph_Object * me, PyObject * args)
{
    Hashgraph * hashgraph = me->hashgraph;

    const char * filename = NULL;
    PyObject * clear_tags_o = NULL;

    if (!PyArg_ParseTuple(args, "s|O", &filename, &clear_tags_o)) {
        return NULL;
    }

    bool clear_tags = true;
    if (clear_tags_o && !PyObject_IsTrue(clear_tags_o)) {
        clear_tags = false;
    }

    try {
        hashgraph->load_tagset(filename, clear_tags);
    } catch (khmer_file_exception &e) {
        PyErr_SetString(PyExc_OSError, e.what());
        return NULL;
    }

    Py_RETURN_NONE;
}

static
PyObject *
hashgraph_save_tagset(khmer_KHashgraph_Object * me, PyObject * args)
{
    Hashgraph * hashgraph = me->hashgraph;

    const char * filename = NULL;

    if (!PyArg_ParseTuple(args, "s", &filename)) {
        return NULL;
    }

    try {
        hashgraph->save_tagset(filename);
    } catch (khmer_file_exception &e) {
        PyErr_SetString(PyExc_OSError, e.what());
        return NULL;
    }

    Py_RETURN_NONE;
}

static
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
    } catch (khmer_file_exception &e) {
        PyErr_SetString(PyExc_OSError, e.what());
        return NULL;
    }

    Py_END_ALLOW_THREADS

    Py_RETURN_NONE;
}

static
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
    } catch (khmer_file_exception &exc) {
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

static
PyObject *
hashgraph__set_tag_density(khmer_KHashgraph_Object * me, PyObject * args)
{
    Hashgraph * hashgraph = me->hashgraph;

    unsigned int d;
    if (!PyArg_ParseTuple(args, "I", &d)) {
        return NULL;
    }

    hashgraph->_set_tag_density(d);

    Py_RETURN_NONE;
}

static
PyObject *
hashgraph__get_tag_density(khmer_KHashgraph_Object * me, PyObject * args)
{
    Hashgraph * hashgraph = me->hashgraph;

    if (!PyArg_ParseTuple(args, "")) {
        return NULL;
    }

    unsigned int d = hashgraph->_get_tag_density();

    return PyLong_FromLong(d);
}

static
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

static
PyObject *
hashgraph_set_partition_id(khmer_KHashgraph_Object * me, PyObject * args)
{
    Hashgraph * hashgraph = me->hashgraph;

    const char * kmer = NULL;
    PartitionID p = 0;

    if (!PyArg_ParseTuple(args, "sI", &kmer, &p)) {
        return NULL;
    }

    hashgraph->partition->set_partition_id(kmer, p);

    Py_RETURN_NONE;
}

static
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

static
PyObject *
hashgraph_get_partition_id(khmer_KHashgraph_Object * me, PyObject * args)
{
    Hashgraph * hashgraph = me->hashgraph;

    const char * kmer = NULL;

    if (!PyArg_ParseTuple(args, "s", &kmer)) {
        return NULL;
    }

    PartitionID partition_id;
    partition_id = hashgraph->partition->get_partition_id(kmer);

    return PyLong_FromLong(partition_id);
}

static
PyObject *
hashgraph_divide_tags_into_subsets(khmer_KHashgraph_Object * me,
                                   PyObject * args)
{
    Hashgraph * hashgraph = me->hashgraph;

    unsigned int subset_size = 0;

    if (!PyArg_ParseTuple(args, "I", &subset_size)) {
        return NULL;
    }

    SeenSet * divvy = new SeenSet;
    hashgraph->divide_tags_into_subsets(subset_size, *divvy);

    PyObject * x = (PyObject *) create_HashSet_Object(divvy,
                   hashgraph->ksize());
    return x;
}

static
PyObject *
hashgraph_count_kmers_within_radius(khmer_KHashgraph_Object * me,
                                    PyObject * args)
{
    Hashgraph * hashgraph = me->hashgraph;

    const char * kmer = NULL;
    unsigned int radius = 0;
    unsigned int max_count = 0;

    if (!PyArg_ParseTuple(args, "sI|I", &kmer, &radius, &max_count)) {
        return NULL;
    }

    unsigned int n;

    Py_BEGIN_ALLOW_THREADS
    Kmer start_kmer = hashgraph->build_kmer(kmer);
    KmerSet seen;
    n = hashgraph->traverse_from_kmer(start_kmer, radius,
                                      seen, max_count);

    Py_END_ALLOW_THREADS

    return PyLong_FromUnsignedLong(n);
}

static
PyObject *
hashgraph_extract_unique_paths(khmer_KHashgraph_Object * me, PyObject * args)
{
    Hashgraph * hashgraph = me->hashgraph;

    const char * sequence = NULL;
    unsigned int min_length = 0;
    float min_unique_f = 0;
    if (!PyArg_ParseTuple(args, "sIf", &sequence, &min_length, &min_unique_f)) {
        return NULL;
    }

    std::vector<std::string> results;
    hashgraph->extract_unique_paths(sequence, min_length, min_unique_f, results);

    PyObject * x = PyList_New(results.size());
    if (x == NULL) {
        return NULL;
    }

    for (unsigned int i = 0; i < results.size(); i++) {
        PyList_SET_ITEM(x, i, PyUnicode_FromString(results[i].c_str()));
    }

    return x;
}

static
PyObject *
hashgraph_consume_and_tag(khmer_KHashgraph_Object * me, PyObject * args)
{
    Hashgraph * hashgraph = me->hashgraph;

    const char * seq;

    if (!PyArg_ParseTuple(args, "s", &seq)) {
        return NULL;
    }

    // call the C++ function, and trap signals => Python

    unsigned long long n_consumed = 0;

    // @CTB needs to normalize
    hashgraph->consume_sequence_and_tag(seq, n_consumed);

    return Py_BuildValue("K", n_consumed);
}

static
PyObject *
hashgraph_get_tags_and_positions(khmer_KHashgraph_Object * me, PyObject * args)
{
    Hashgraph * hashgraph = me->hashgraph;

    const char * seq;

    if (!PyArg_ParseTuple(args, "s", &seq)) {
        return NULL;
    }

    // call the C++ function, and trap signals => Python

    std::vector<unsigned int> posns;
    std::vector<HashIntoType> tags;

    unsigned int pos = 1;
    KmerIterator kmers(seq, hashgraph->ksize());

    while (!kmers.done()) {
        HashIntoType kmer = kmers.next();
        if (set_contains(hashgraph->all_tags, kmer)) {
            posns.push_back(pos);
            tags.push_back(kmer);
        }
        pos++;
    }

    PyObject * tag = nullptr;
    PyObject * posns_list = PyList_New(posns.size());
    for (size_t i = 0; i < posns.size(); i++) {
        convert_HashIntoType_to_PyObject(tags[i], &tag);
        PyObject * tup = Py_BuildValue("IO", posns[i], tag);
        PyList_SET_ITEM(posns_list, i, tup);
    }

    return posns_list;
}

static
PyObject *
hashgraph_find_all_tags_list(khmer_KHashgraph_Object * me, PyObject * args)
{
    Hashgraph * hashgraph = me->hashgraph;

    const char * kmer_s = NULL;

    if (!PyArg_ParseTuple(args, "s", &kmer_s)) {
        return NULL;
    }

    if (strlen(kmer_s) != hashgraph->ksize()) {
        PyErr_SetString(PyExc_ValueError,
                        "k-mer length must equal the counting table k-mer size");
        return NULL;
    }

    SeenSet * tags = new SeenSet;

    Kmer start_kmer = hashgraph->build_kmer(kmer_s);

    Py_BEGIN_ALLOW_THREADS

    hashgraph->partition->find_all_tags(start_kmer, *tags,
                                        hashgraph->all_tags);

    Py_END_ALLOW_THREADS

    PyObject * x = (PyObject *) create_HashSet_Object(tags,
                   hashgraph->ksize());
    return x;
}

static
PyObject *
hashgraph_consume_seqfile_and_tag(khmer_KHashgraph_Object * me, PyObject * args)
{
    Hashgraph * hashgraph = me->hashgraph;

    const char * filename;

    if (!PyArg_ParseTuple(args, "s", &filename)) {
        return NULL;
    }

    // call the C++ function, and trap signals => Python

    unsigned long long n_consumed;
    unsigned int total_reads;

    try {
        hashgraph->consume_seqfile_and_tag<FastxReader>(filename, total_reads, n_consumed);
    } catch (khmer_file_exception &exc) {
        PyErr_SetString(PyExc_OSError, exc.what());
        return NULL;
    } catch (khmer_value_exception &exc) {
        PyErr_SetString(PyExc_ValueError, exc.what());
        return NULL;
    }

    return Py_BuildValue("IK", total_reads, n_consumed);
}

static PyMethodDef khmer_hashgraph_methods[] = {
    //
    // graph/traversal functionality
    //

    {
        "neighbors",
        (PyCFunction)hashgraph_neighbors, METH_VARARGS,
        "Get a list of neighbor nodes for this k-mer.",
    },
    {
        "calc_connected_graph_size",
        (PyCFunction)hashgraph_calc_connected_graph_size, METH_VARARGS, ""
    },
    {
        "kmer_degree",
        (PyCFunction)hashgraph_kmer_degree, METH_VARARGS,
        "Calculate the number of immediate neighbors this k-mer has in "
        "the graph."
    },
    {
        "count_kmers_within_radius",
        (PyCFunction)hashgraph_count_kmers_within_radius, METH_VARARGS,
        "Calculate the number of neighbors with given radius in the graph."
    },
    {
        "find_high_degree_nodes",
        (PyCFunction)hashgraph_find_high_degree_nodes, METH_VARARGS,
        "Examine the given sequence for degree > 2 nodes and add to  "
        "list; used in graph contraction.",
    },
    {
        "traverse_linear_path",
        (PyCFunction)hashgraph_traverse_linear_path, METH_VARARGS,
        "Traverse the path through the graph starting with the given "
        "k-mer and avoiding high-degree nodes, finding (and returning) "
        "traversed k-mers and any encountered high-degree nodes.",
    },
    {
        "assemble_linear_path",
        (PyCFunction)hashgraph_assemble_linear_path, METH_VARARGS,
        "Assemble a purely linear path starting with the given "
        "k-mer, returning traversed k-mers and any encountered high-degree "
        "nodes.",
    },

    //
    // tagging / sparse graph functionality
    //

    {
        "consume_and_tag",
        (PyCFunction)hashgraph_consume_and_tag, METH_VARARGS,
        "Consume a sequence and tag it."
    },
    {
        "get_tags_and_positions",
        (PyCFunction)hashgraph_get_tags_and_positions, METH_VARARGS,
        "Retrieve tags and their positions in a sequence."
    },
    {
        "find_all_tags_list",
        (PyCFunction)hashgraph_find_all_tags_list, METH_VARARGS,
        "Find all tags within range of the given k-mer, return as list"
    },
    {
        "consume_seqfile_and_tag",
        (PyCFunction)hashgraph_consume_seqfile_and_tag, METH_VARARGS,
        "Consume all sequences in a FASTA/FASTQ file and tag the resulting "
        "graph."
    },
    {
        "extract_unique_paths",
        (PyCFunction)hashgraph_extract_unique_paths, METH_VARARGS,
        "@CTB remove."
    },
    {
        "print_tagset",
        (PyCFunction)hashgraph_print_tagset, METH_VARARGS,
        "Print out all of the tags."
    },
    {
        "add_tag",
        (PyCFunction)hashgraph_add_tag, METH_VARARGS,
        "Add a k-mer to the tagset."
    },
    {
        "get_tagset",
        (PyCFunction)hashgraph_get_tagset, METH_VARARGS,
        "Get all tagged k-mers as DNA strings."
    },
    {
        "load_tagset",
        (PyCFunction)hashgraph_load_tagset, METH_VARARGS,
        "Load tags from a file."
    },
    {
        "save_tagset",
        (PyCFunction)hashgraph_save_tagset, METH_VARARGS,
        "Save tags to a file."
    },
    {
        "n_tags",
        (PyCFunction)hashgraph_n_tags, METH_VARARGS,
        "Return the count of all tags."
    },
    {
        "divide_tags_into_subsets",
        (PyCFunction)hashgraph_divide_tags_into_subsets, METH_VARARGS,
        "Divide tags equally up into subsets of given size."
    },
    {
        "_get_tag_density",
        (PyCFunction)hashgraph__get_tag_density, METH_VARARGS,
        "Get the tagging density."
    },
    {
        "_set_tag_density",
        (PyCFunction)hashgraph__set_tag_density, METH_VARARGS,
        "Set the tagging density."
    },

    //
    // partitioning
    //
    {
        "do_subset_partition",
        (PyCFunction)hashgraph_do_subset_partition, METH_VARARGS,
        "Partition the graph starting from a given subset of tags."
    },
    {
        "find_all_tags",
        (PyCFunction)hashgraph_find_all_tags, METH_VARARGS,
        "Starting from the given k-mer, find all closely connected tags."
    },
    {
        "assign_partition_id",
        (PyCFunction)hashgraph_assign_partition_id, METH_VARARGS,
        "Assign a partition ID to a given tag."
    },
    {
        "output_partitions",
        (PyCFunction)hashgraph_output_partitions, METH_VARARGS,
        "Write out sequences in given filename to another file, annotating "
        "with partition IDs."
    },
    {
        "load_partitionmap",
        (PyCFunction)hashgraph_load_partitionmap, METH_VARARGS,
        "Load a partitionmap for a given subset."
    },
    {
        "save_partitionmap",
        (PyCFunction)hashgraph_save_partitionmap, METH_VARARGS,
        "Save a partitionmap for the given subset."
    },
    {
        "_validate_partitionmap",
        (PyCFunction)hashgraph__validate_partitionmap, METH_VARARGS,
        "Run internal validation checks."
    },
    {
        "consume_seqfile_and_tag_with_reads_parser",
        (PyCFunction)hashgraph_consume_seqfile_and_tag_with_reads_parser,
        METH_VARARGS,
        "Count all k-mers using the given reads parser"
    },
    {
        "consume_partitioned_fasta",
        (PyCFunction)hashgraph_consume_partitioned_fasta, METH_VARARGS,
        "Count all k-mers in a given file"
    },
    {
        "merge_subset",
        (PyCFunction)hashgraph_merge_subset, METH_VARARGS,
        "Merge the given subset into this one."
    },
    {
        "merge_subset_from_disk",
        (PyCFunction)hashgraph_merge_from_disk, METH_VARARGS,
        "Merge the given subset (filename) into this one."
    },
    {
        "count_partitions",
        (PyCFunction)hashgraph_count_partitions, METH_VARARGS,
        "Count the number of partitions in the master partitionmap."
    },
    {
        "subset_count_partitions",
        (PyCFunction)hashgraph_subset_count_partitions, METH_VARARGS,
        "Count the number of partitions in this subset partitionmap."
    },
    {
        "subset_partition_size_distribution",
        (PyCFunction)hashgraph_subset_partition_size_distribution,
        METH_VARARGS,
        "Get the size distribution of partitions in this subset."
    },
    {
        "save_subset_partitionmap",
        (PyCFunction)hashgraph_save_subset_partitionmap, METH_VARARGS,
        "Save the partition map for this subset."
    },
    {
        "load_subset_partitionmap",
        (PyCFunction)hashgraph_load_subset_partitionmap, METH_VARARGS,
        "Save the partition map for this subset."
    },
    {
        "_validate_subset_partitionmap",
        (PyCFunction)hashgraph__validate_subset_partitionmap, METH_VARARGS,
        "Run internal validation checks on this subset."
    },
    {
        "set_partition_id",
        (PyCFunction)hashgraph_set_partition_id, METH_VARARGS,
        "Set the partition ID for this tag."
    },
    {
        "join_partitions",
        (PyCFunction)hashgraph_join_partitions, METH_VARARGS,
        "Join the partitions of these two tags."
    },
    {
        "get_partition_id",
        (PyCFunction)hashgraph_get_partition_id, METH_VARARGS,
        "Get the partition ID of this tag."
    },
    {
        "repartition_largest_partition",
        (PyCFunction)hashgraph_repartition_largest_partition, METH_VARARGS,
        "Repartition the largest partition (in the face of stop tags)."
    },

    // stop tags
    {
        "load_stop_tags",
        (PyCFunction)hashgraph_load_stop_tags, METH_VARARGS,
        "Load the set of stop tags."
    },
    {
        "save_stop_tags",
        (PyCFunction)hashgraph_save_stop_tags, METH_VARARGS,
        "Save the set of stop tags."
    },
    {
        "print_stop_tags",
        (PyCFunction)hashgraph_print_stop_tags, METH_VARARGS,
        "Print out the set of stop tags."
    },
    {
        "trim_on_stoptags",
        (PyCFunction)hashgraph_trim_on_stoptags, METH_VARARGS,
        "Trim the reads on the given stop tags."
    },
    {
        "add_stop_tag",
        (PyCFunction)hashgraph_add_stop_tag, METH_VARARGS,
        "Add this k-mer as a stop tag."
    },
    {
        "get_stop_tags",
        (PyCFunction)hashgraph_get_stop_tags, METH_VARARGS,
        "Return a DNA list of all of the stop tags."
    },
    {NULL, NULL, 0, NULL}           /* sentinel */
};
