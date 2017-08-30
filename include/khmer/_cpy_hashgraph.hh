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
#ifndef _CPY_HASHGRAPH_HH
#define _CPY_HASHGRAPH_HH

#include <Python.h>
#include "_cpy_utils.hh"
#include "_cpy_hashtable.hh"
#include "oxli/hashgraph.hh"


namespace khmer {


typedef struct {
    khmer_KHashtable_Object khashtable;
    oxli::Hashgraph * hashgraph;
} khmer_KHashgraph_Object;


extern PyTypeObject khmer_KHashgraph_Type
CPYCHECKER_TYPE_OBJECT_FOR_TYPEDEF("khmer_KHashgraph_Object");


extern PyMethodDef khmer_hashgraph_methods[];


//
// Method definitions
//


PyObject *
hashgraph_find_high_degree_nodes(khmer_KHashgraph_Object * me, PyObject * args);


PyObject *
hashgraph_neighbors(khmer_KHashgraph_Object * me, PyObject * args);


PyObject *
hashgraph_traverse_linear_path(khmer_KHashgraph_Object * me, PyObject * args);


PyObject *
hashgraph_assemble_linear_path(khmer_KHashgraph_Object * me, PyObject * args);


PyObject *
hashgraph_n_tags(khmer_KHashgraph_Object * me, PyObject * args);


PyObject *
hashgraph_print_stop_tags(khmer_KHashgraph_Object * me, PyObject * args);


PyObject *
hashgraph_print_tagset(khmer_KHashgraph_Object * me, PyObject * args);


PyObject *
hashgraph_load_stop_tags(khmer_KHashgraph_Object * me, PyObject * args);



PyObject *
hashgraph_save_stop_tags(khmer_KHashgraph_Object * me, PyObject * args);



PyObject *
hashgraph_repartition_largest_partition(khmer_KHashgraph_Object * me,
                                        PyObject * args);


PyObject *
hashgraph_calc_connected_graph_size(khmer_KHashgraph_Object * me,
                                    PyObject * args);


PyObject *
hashgraph_kmer_degree(khmer_KHashgraph_Object * me, PyObject * args);


PyObject *
hashgraph_trim_on_stoptags(khmer_KHashgraph_Object * me, PyObject * args);


PyObject *
hashgraph_do_subset_partition(khmer_KHashgraph_Object * me, PyObject * args);



PyObject *
hashgraph_merge_subset(khmer_KHashgraph_Object * me, PyObject * args);


PyObject *
hashgraph_merge_from_disk(khmer_KHashgraph_Object * me, PyObject * args);


PyObject *
hashgraph_consume_seqfile_and_tag_with_reads_parser(khmer_KHashgraph_Object * me,
        PyObject * args);


PyObject *
hashgraph_consume_partitioned_fasta(khmer_KHashgraph_Object * me,
                                    PyObject * args);


PyObject *
hashgraph_find_all_tags(khmer_KHashgraph_Object * me, PyObject * args);


PyObject *
hashgraph_assign_partition_id(khmer_KHashgraph_Object * me, PyObject * args);


PyObject *
hashgraph_add_tag(khmer_KHashgraph_Object * me, PyObject * args);


PyObject *
hashgraph_add_stop_tag(khmer_KHashgraph_Object * me, PyObject * args);


PyObject *
hashgraph_get_stop_tags(khmer_KHashgraph_Object * me, PyObject * args);


PyObject *
hashgraph_get_tagset(khmer_KHashgraph_Object * me, PyObject * args);


PyObject *
hashgraph_output_partitions(khmer_KHashgraph_Object * me, PyObject * args);


PyObject *
hashgraph_save_partitionmap(khmer_KHashgraph_Object * me, PyObject * args);


PyObject *
hashgraph_load_partitionmap(khmer_KHashgraph_Object * me, PyObject * args);


PyObject *
hashgraph__validate_partitionmap(khmer_KHashgraph_Object * me, PyObject * args);


PyObject *
hashgraph_count_partitions(khmer_KHashgraph_Object * me, PyObject * args);


PyObject *
hashgraph_subset_count_partitions(khmer_KHashgraph_Object * me, PyObject * args);


PyObject *
hashgraph_subset_partition_size_distribution(khmer_KHashgraph_Object * me,
        PyObject * args);


PyObject *
hashgraph_load_tagset(khmer_KHashgraph_Object * me, PyObject * args);


PyObject *
hashgraph_save_tagset(khmer_KHashgraph_Object * me, PyObject * args);


PyObject *
hashgraph_save_subset_partitionmap(khmer_KHashgraph_Object * me,
                                   PyObject * args);


PyObject *
hashgraph_load_subset_partitionmap(khmer_KHashgraph_Object * me,
                                   PyObject * args);


PyObject *
hashgraph__set_tag_density(khmer_KHashgraph_Object * me, PyObject * args);


PyObject *
hashgraph__get_tag_density(khmer_KHashgraph_Object * me, PyObject * args);


PyObject *
hashgraph__validate_subset_partitionmap(khmer_KHashgraph_Object * me,
                                        PyObject * args);


PyObject *
hashgraph_set_partition_id(khmer_KHashgraph_Object * me, PyObject * args);


PyObject *
hashgraph_join_partitions(khmer_KHashgraph_Object * me, PyObject * args);


PyObject *
hashgraph_get_partition_id(khmer_KHashgraph_Object * me, PyObject * args);


PyObject *
hashgraph_divide_tags_into_subsets(khmer_KHashgraph_Object * me,
                                   PyObject * args);


PyObject *
hashgraph_count_kmers_within_radius(khmer_KHashgraph_Object * me,
                                    PyObject * args);


PyObject *
hashgraph_extract_unique_paths(khmer_KHashgraph_Object * me, PyObject * args);


PyObject *
hashgraph_consume_and_tag(khmer_KHashgraph_Object * me, PyObject * args);


PyObject *
hashgraph_get_tags_for_sequence(khmer_KHashgraph_Object * me, PyObject * args);


PyObject *
hashgraph_get_tags_and_positions(khmer_KHashgraph_Object * me, PyObject * args);


PyObject *
hashgraph_find_all_tags_list(khmer_KHashgraph_Object * me, PyObject * args);


PyObject *
hashgraph_consume_seqfile_and_tag(khmer_KHashgraph_Object * me, PyObject * args);

}

#endif
