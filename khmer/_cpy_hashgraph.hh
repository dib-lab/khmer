typedef struct {
    khmer_KHashtable_Object khashtable;
    Hashgraph * hashgraph;
} khmer_KHashgraph_Object;

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
    "hashgraph object",           /* tp_doc */
    0,                       /* tp_traverse */
    0,                       /* tp_clear */
    0,                       /* tp_richcompare */
    0,                       /* tp_weaklistoffset */
    0,                       /* tp_iter */
    0,                       /* tp_iternext */
    0,                       /* tp_methods */
    0,                       /* tp_members */
    0,                       /* tp_getset */
    0,                       /* tp_base */
    0,                       /* tp_dict */
    0,                       /* tp_descr_get */
    0,                       /* tp_descr_set */
    0,                       /* tp_dictoffset */
    0,                       /* tp_init */
    0,                       /* tp_alloc */
    0,                       /* tp_new */
};


static PyMethodDef khmer_hashgraph_methods[] = {
    //
    // graph/traversal functionality
    //

    {
        "neighbors",
        (PyCFunction)hashtable_neighbors, METH_VARARGS,
        "Get a list of neighbor nodes for this k-mer.",
    },
    {
        "calc_connected_graph_size",
        (PyCFunction)hashtable_calc_connected_graph_size, METH_VARARGS, ""
    },
    {
        "kmer_degree",
        (PyCFunction)hashtable_kmer_degree, METH_VARARGS,
        "Calculate the number of immediate neighbors this k-mer has in "
        "the graph."
    },
    {
        "count_kmers_within_radius",
        (PyCFunction)hashtable_count_kmers_within_radius, METH_VARARGS,
        "Calculate the number of neighbors with given radius in the graph."
    },

    {
        "find_high_degree_nodes",
        (PyCFunction)hashtable_find_high_degree_nodes, METH_VARARGS,
        "Examine the given sequence for degree > 2 nodes and add to  "
        "list; used in graph contraction.",
    },
    {
        "traverse_linear_path",
        (PyCFunction)hashtable_traverse_linear_path, METH_VARARGS,
        "Traverse the path through the graph starting with the given "
        "k-mer and avoiding high-degree nodes, finding (and returning) "
        "traversed k-mers and any encountered high-degree nodes.",
    },
    {
        "assemble_linear_path",
        (PyCFunction)hashtable_assemble_linear_path, METH_VARARGS,
        "Assemble a purely linear path starting with the given "
        "k-mer, returning traversed k-mers and any encountered high-degree "
        "nodes.",
    },

    //
    // tagging / sparse graph functionality
    //

    { "consume_and_tag", (PyCFunction)hashtable_consume_and_tag, METH_VARARGS, "Consume a sequence and tag it" },
    { "get_tags_and_positions", (PyCFunction)hashtable_get_tags_and_positions, METH_VARARGS, "Retrieve tags and their positions in a sequence." },
    { "find_all_tags_list", (PyCFunction)hashtable_find_all_tags_list, METH_VARARGS, "Find all tags within range of the given k-mer, return as list" },
    { "consume_fasta_and_tag", (PyCFunction)hashtable_consume_fasta_and_tag, METH_VARARGS, "Count all k-mers in a given file" },
    { "get_median_count", (PyCFunction)hashtable_get_median_count, METH_VARARGS, "Get the median, average, and stddev of the k-mer counts in the string" },
    { "median_at_least", (PyCFunction)hashtable_median_at_least, METH_VARARGS, "Return true if the median is at least the given cutoff" },
    { "extract_unique_paths", (PyCFunction)hashtable_extract_unique_paths, METH_VARARGS, "" },
    { "print_tagset", (PyCFunction)hashtable_print_tagset, METH_VARARGS, "" },
    { "add_tag", (PyCFunction)hashtable_add_tag, METH_VARARGS, "" },
    { "get_tagset", (PyCFunction)hashtable_get_tagset, METH_VARARGS, "" },
    { "load_tagset", (PyCFunction)hashtable_load_tagset, METH_VARARGS, "" },
    { "save_tagset", (PyCFunction)hashtable_save_tagset, METH_VARARGS, "" },
    { "n_tags", (PyCFunction)hashtable_n_tags, METH_VARARGS, "" },
    { "divide_tags_into_subsets", (PyCFunction)hashtable_divide_tags_into_subsets, METH_VARARGS, "" },
    { "_get_tag_density", (PyCFunction)hashtable__get_tag_density, METH_VARARGS, "" },
    { "_set_tag_density", (PyCFunction)hashtable__set_tag_density, METH_VARARGS, "" },

    // partitioning
    { "do_subset_partition", (PyCFunction)hashtable_do_subset_partition, METH_VARARGS, "" },
    { "find_all_tags", (PyCFunction)hashtable_find_all_tags, METH_VARARGS, "" },
    { "assign_partition_id", (PyCFunction)hashtable_assign_partition_id, METH_VARARGS, "" },
    { "output_partitions", (PyCFunction)hashtable_output_partitions, METH_VARARGS, "" },
    { "load_partitionmap", (PyCFunction)hashtable_load_partitionmap, METH_VARARGS, "" },
    { "save_partitionmap", (PyCFunction)hashtable_save_partitionmap, METH_VARARGS, "" },
    { "_validate_partitionmap", (PyCFunction)hashtable__validate_partitionmap, METH_VARARGS, "" },
    {
        "consume_fasta_and_tag_with_reads_parser", (PyCFunction)hashtable_consume_fasta_and_tag_with_reads_parser,
        METH_VARARGS, "Count all k-mers using a given reads parser"
    },
    { "consume_partitioned_fasta", (PyCFunction)hashtable_consume_partitioned_fasta, METH_VARARGS, "Count all k-mers in a given file" },
    { "merge_subset", (PyCFunction)hashtable_merge_subset, METH_VARARGS, "" },
    { "merge_subset_from_disk", (PyCFunction)hashtable_merge_from_disk, METH_VARARGS, "" },
    { "count_partitions", (PyCFunction)hashtable_count_partitions, METH_VARARGS, "" },
    { "subset_count_partitions", (PyCFunction)hashtable_subset_count_partitions, METH_VARARGS, "" },
    { "subset_partition_size_distribution", (PyCFunction)hashtable_subset_partition_size_distribution, METH_VARARGS, "" },
    { "save_subset_partitionmap", (PyCFunction)hashtable_save_subset_partitionmap, METH_VARARGS },
    { "load_subset_partitionmap", (PyCFunction)hashtable_load_subset_partitionmap, METH_VARARGS },
    { "_validate_subset_partitionmap", (PyCFunction)hashtable__validate_subset_partitionmap, METH_VARARGS, "" },
    { "set_partition_id", (PyCFunction)hashtable_set_partition_id, METH_VARARGS, "" },
    { "join_partitions", (PyCFunction)hashtable_join_partitions, METH_VARARGS, "" },
    { "get_partition_id", (PyCFunction)hashtable_get_partition_id, METH_VARARGS, "" },
    { "repartition_largest_partition", (PyCFunction)hashtable_repartition_largest_partition, METH_VARARGS, "" },

    // stop tags
    { "load_stop_tags", (PyCFunction)hashtable_load_stop_tags, METH_VARARGS, "" },
    { "save_stop_tags", (PyCFunction)hashtable_save_stop_tags, METH_VARARGS, "" },
    { "print_stop_tags", (PyCFunction)hashtable_print_stop_tags, METH_VARARGS, "" },
    { "trim_on_stoptags", (PyCFunction)hashtable_trim_on_stoptags, METH_VARARGS, "" },
    { "add_stop_tag", (PyCFunction)hashtable_add_stop_tag, METH_VARARGS, "" },
    { "get_stop_tags", (PyCFunction)hashtable_get_stop_tags, METH_VARARGS, "" },
    {NULL, NULL, 0, NULL}           /* sentinel */
};

#define is_hashgraph_obj(v)  (Py_TYPE(v) == &khmer_KHashgraph_Type)
