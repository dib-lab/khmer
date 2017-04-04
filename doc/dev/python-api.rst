Some notes on the Python API
============================

Some programming guidelines
---------------------------

...

Valid and invalid DNA
~~~~~~~~~~~~~~~~~~~~~

Generally speaking, the Python API for khmer and oxli assume that
they are receiving valid DNA (ACGT).  Low-level hash functions like
``hash(kmer)`` and mid-level hash functions like ``hash_kmer_hashes(str)``
will neither check for correct DNA nor is their output with respect to
incorrect DNA characters specified.

However, bulk loading functions provide a ``cleaned_seq`` attribute that
will ... document me here.

Table types
-----------

Type names consist of two parts. The first part indicates how far the type
can count and the second part whether it is a table or a graph.

Possible choices for the first part:
* Node, uses 1bit counter
* SmallCount, uses a 4bit counter
* Count, uses a 8bit counter

Possible choices for the second part:
* Table, keep track of kmers
* Graph, navigate, tag, etc the de Bruijn graph formed by the k-mers

C++ class name:

Python methods:

* k = ksize() - return the k-mer size (a positive integer).

* hashval = hash(dna_kmer) - return the result of hashing ``dna_kmer``, which will be a non-negative integer. ``len(dna_kmer)`` must be exactly the k-mer size.  Which hash function is used is dependent on the table type (@document).

* dna_kmer = reverse_hash(hashval) - return a DNA string that will hash to ``hashval``.  If there are multiple such strings, return only one.  May be unimplemented for particular table types in which case a ValueError will be returned.

* sizelist = hashsizes() - return the list of table sizes used in construction.

* n = n_unique_kmers() - retrieve an estimate for the number of unique k-mers inserted into the table. Note, this may be order dependent.
  
* n = n_occupied() - retrieve the fraction of bins occupied in the table.

* add(dna_kmer_or_hashval) - increment the count associated with either a DNA k-mer or a hashval.  Depending on max count for the tabletype and bigcount settings, the count may top out at 1, 15, 255, or 65535. (@CTB add method for retrieving max_count)
  
* count (synonym for add)
  
* get(dna_kmer or hashval) - retrieve the count associated with a DNA k-mer or a hashval.

* list_of_strings = get_kmers(seq) - return the list of k-mer strings in the given sequence.
* list_of_hashes = get_kmer_hashes(seq) - return the list of the hashed k-mers in the given sequence.
* hashset = get_kmer_hashes_as_hashset(seq) - return the hashset of hashed k-mers in the given sequence.

* list_of_counts = get_kmer_counts(seq) - return the list of the counts of the k-mers in the given sequence.

* save(filename) - save the data to a file on disk.
* load(filename) - load the data from a file on disk.

* min_count = get_min_count(seq) - return the minimum count for k-mers in seq
* med_count, avg_count, stddev_count = get_median_count(seq) - return the median, average, and stddev of the counts for k-mers in the sequence.
  
* max_count = get_max_count(seq) - return the maximum count for k-mers in the sequence.
  
* num_kmers = consume(seq) - count all the k-mers in the given DNA string. @CTB should be consume string
* consume_fasta(filename) - count all the k-mers in a (DNA) FASTA/FASTQ file.
* consume_fasta_with_reads_parser(khmer.ReadParser object) - count all the k-mers in (DNA) sequences returned from a ReadParser object.  This can be used to ensure various forms of pairing are present, etc.

* (trim_seq, trim_pos) = trim_on_abundance(seq, abund) - trim the sequence at the first k-mer with an abundance strictly below (<) the provided abundance.
* (trim_seq, trim_pos) = trim_below_abundance(seq, abund) - trim the sequence at the first k-mer with an abundance strictly above (>) the provided abundance.
* list_of_posns = find_spectral_error_positions(seq, abund) - find potential locations of errors in input sequence, where k-mers are below given abundance.

* set_use_bigcount(bool) - some table types (Counttable and Countgraph) support counting past their max value, using a (memory intensive) C++ stl::map. That turns on this "bigcount" behavior.  This will raise a ValueError when called on a table type that does not support it.
* get_use_bigcount - return the bigcount value (False by default).

* dist = abundance_distribution(filename, tracking_obj) - generate an abundance distribution for the k-mers in the given file, using the tracking_obj to avoid double-counting identical k-mers.
* dist = abundance_distribution_with_reads_parser(readparser_obj, tracking_obj) - generate an abundance distribution for the k-mers loaded from the given readparser, using the tracking_obj to avoid double-counting identical k-mers.

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

Countgraph:

* get_raw_tables
* do_subset_partition_with_abundance

Nodegraph:

* update
* get_raw_tables
