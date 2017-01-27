Types --

Type names consist of two parts. The first part indicates how far the type
can count and the second part whether it is a table or a graph.

Possible choices for the first part:
* Node, uses 1bit counter
* SmallCount, uses a 4bit counter
* Count, uses a 8bit counter

Possible choices for the second part:
* Table, keep track of kmers
* Graph, navigate, tag, etc the de Bruijn graph formed by the khmers


Table types
-----------

C++ class name:

Python methods:

* k = ksize() - return the k-mer size (a positive integer).

* hashval = hash(dna_kmer) - return the result of hashing ``dna_kmer``, which will be a non-negative integer. ``len(dna_kmer)`` must be exactly the k-mer size.  Which hash function is used is dependent on the table type (@document).

* dna_kmer = reverse_hash(hashval) - return a DNA string that will hash to ``hashval``.  If there are multiple such strings, return only one.  May be unimplemented for particular table types in which case a ValueError will be returned.

* sizelist = hashsizes() - return the list of table sizes used in construction.

* n_unique_kmers - foo.
* n_occupied - foo.

* add(dna_kmer_or_hashval) - increment the count associated with either a DNA k-mer or a hashval.  Depending on max count for the tabletype and bigcount settings, the count may top out at 1, 15, 255, or 65535. (@CTB add method for retrieving max_count)
* count (synonym for add)
* get(dna_kmer or hashval) - retrieve the count associated with a DNA k-mer or a hashval.

* save
* get_kmers
* get_kmer_hashes
* get_kmer_hashes_as_hashset
* get_kmer_counts
* get_min_count
* get_median_count
* get_max_count
* consume (should be consume string?)
* consume_fasta
* consume_fasta_with_reads_parser
* load
* set_use_bigcount
* get_use_bigcount
* abundance_distribution
* abundance_distribution_with_reads_paresr
* trim_on_abundance
* trim_below_abundance
* find_spectral_error_positions

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
