Types --

Table types
-----------

C++ class name:

Python methods:

* ksize
* hash
* reverse_hash (@CTB)
* hashsizes
* n_unique_kmers
* n_occupied
* count
* add (see count)
* consume (should be consume string?)
* consume_fasta
* consume_fasta_with_reads_parser
* get
* load
* save
* get_kmers
* get_kmer_hashes
* get_kmer_hashes_as_hashset
* get_kmer_counts
* set_use_bigcount
* get_use_bigcount
* get_min_count
* get_median_count
* get_max_count
* trim_on_abundance
* trim_below_abundance
* find_spectral_error_positions
* abundance_distribution
* abundance_distribution_with_reads_paresr

Counting types
--------------

Countgraph:

* get_raw_tables
* do_subset_partition_with_abundance

Nodegraph:

* update
* get_raw_tables

Graph types
-----------

All the methods of table types, and in addition:

* neighbors
* calc_connected_graph_size
* kmer_degree
* count_kmers_within_radius
* find_high_degree_nodes
* traverse_linear_path
* assemble_linear_path
* consume_and_tag
* get_tags_and_positions
* find_all_tags_list
* consume_fasta_and_tag
* extract_unique_paths
* print_tagset
* add_tag
* get_tagset
* load_tagset
* save_tagset
* n_tags
* divide_tags_into_subsets
* _get_tag_density
* _set_tag_density
* do_subset_partition
* find_all_tags
* assign_partition_id
* output_partitions
* load_partitionmap
* save_partitionmap
* _validate_partitionmap
* consume_fasta_and_tag_with_reads_parser
* consume_partitioned_fasta
* merge_subset
* merge_subset_from_disk
* count_partitions
* subset_count_partitions
* subset_partition_size_distribution
* save_subset_partitionmap
* load_subset_partitionmap
* _validate_subset_partitionmap
* set_partition_id
* join_partitions
* get_partition_id
* repartition_latest_partition
* load_stop_tags
* save_stop_tags
* print_stop_tags
* trim_on_stoptags
* add_stop_tags
* get_stop_tags

Smallcountgraph:

* get_raw_tables
