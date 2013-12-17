//
// This file is part of khmer, http://github.com/ged-lab/khmer/, and is
// Copyright (C) Michigan State University, 2009-2013. It is licensed under
// the three-clause BSD license; see doc/LICENSE.txt. Contact: ctb@msu.edu
//

#include "labelhash.hh"

using namespace khmer;

/*
 * @camillescott
 * Might be time for a refactor: could do a general consume_fasta
 * function which accepts a consume_sequence function pointer as a parameter
 */

void
LabelHash::consume_fasta_and_tag_with_labels(
  std:: string const  &filename,
  unsigned int	      &total_reads, unsigned long long	&n_consumed,
  CallbackFn	      callback,	    void *		callback_data
)
{
  khmer:: Config    &the_config	  = khmer:: get_active_config( );

  // Note: Always assume only 1 thread if invoked this way.
  IParser *	  parser = 
  IParser::get_parser(
    filename, 1, the_config.get_reads_input_buffer_size( ),
    the_config.get_reads_parser_trace_level( )
  );


  consume_fasta_and_tag_with_labels(
    parser,
    total_reads, n_consumed,
    callback, callback_data
  );

  delete parser;
}

void
LabelHash::consume_fasta_and_tag_with_labels(
    read_parsers:: IParser *  parser,
    unsigned int		    &total_reads,   unsigned long long	&n_consumed,
    CallbackFn		    callback,	    void *		callback_data
  )
  {
    Hasher		  &hasher		= 
    _get_hasher( parser->uuid( ) );
    unsigned int		  total_reads_LOCAL	= 0;
  #if (0) // Note: Used with callback - currently disabled.
    unsigned long long int  n_consumed_LOCAL	= 0;
  #endif
    Read			  read;

    // TODO? Delete the following assignments.
    total_reads = 0;
    n_consumed = 0;
    
    hasher.trace_logger(
      TraceLogger:: TLVL_DEBUG2,
      "Starting trace of 'consume_fasta_and_tag_with_labels'....\n"
    );
    
    Label _tag_label = 0;

    Label * the_label;
    // Iterate through the reads and consume their k-mers.
    while (!parser->is_complete( ))
    {
      unsigned long long this_n_consumed   = 0;

      read = parser->get_next_read( );

      if (check_and_normalize_read( read.sequence ))
      {
        // TODO: make threadsafe!
        the_label = check_and_allocate_label(_tag_label);
        consume_sequence_and_tag_with_labels( read.sequence,
					      this_n_consumed,
					      *the_label );
	    _tag_label++;

  #ifdef WITH_INTERNAL_METRICS
        hasher.pmetrics.start_timers( );
  #endif
  #if (0) // Note: Used with callback - currently disabled.
        n_consumed_LOCAL  = __sync_add_and_fetch( &n_consumed, this_n_consumed );
  #else
        __sync_add_and_fetch( &n_consumed, this_n_consumed );
  #endif
        total_reads_LOCAL = __sync_add_and_fetch( &total_reads, 1 );
  #ifdef WITH_INTERNAL_METRICS
        hasher.pmetrics.stop_timers( );
        hasher.pmetrics.accumulate_timer_deltas(
	  (uint32_t)HashTablePerformanceMetrics:: MKEY_TIME_UPDATE_TALLIES
        );
  #endif
      }

      if (0 == (total_reads_LOCAL % 10000))
        hasher.trace_logger(
	  TraceLogger:: TLVL_DEBUG3,
	  "Total number of reads processed: %llu\n",
	  (unsigned long long int)total_reads_LOCAL
        );

      // TODO: Figure out alternative to callback into Python VM
      //       Cannot use in multi-threaded operation.
  #if (0)
        // run callback, if specified
        if (total_reads_TL % CALLBACK_PERIOD == 0 && callback) {
	  std::cout << "n tags: " << all_tags.size() << "\n";
	  try {
	    callback("consume_fasta_and_tag_with_labels", callback_data, total_reads_TL,
		     n_consumed);
	  } catch (...) {
	    delete parser;
	    throw;
	  }
        }
  #endif // 0

    } // while reads left for parser

  }

void LabelHash::consume_partitioned_fasta_and_tag_with_labels(const std::string &filename,
					  unsigned int &total_reads,
					  unsigned long long &n_consumed,
					  CallbackFn callback,
					  void * callback_data)
{
  total_reads = 0;
  n_consumed = 0;

  IParser* parser = IParser::get_parser(filename.c_str());
  Read read;

  string seq = "";

  // reset the master subset partition
  delete partition;
  partition = new SubsetPartition(this);

  //
  // iterate through the FASTA file & consume the reads.
  //
  Label * c;
  PartitionID p;
  while(!parser->is_complete())  {
    read = parser->get_next_read();
    seq = read.sequence;

    if (check_and_normalize_read(seq)) {
      // First, figure out what the partition is (if non-zero), and save that.
      p = _parse_partition_id(read.name);
      c = check_and_allocate_label(p);

      consume_sequence_and_tag_with_labels( seq,
					      n_consumed,
					      *c );
    }
	       
    // reset the sequence info, increment read number
    total_reads++;

    // run callback, if specified
    if (total_reads % CALLBACK_PERIOD == 0 && callback) {
      try {
        callback("consume_partitioned_fasta_and_tag_with_labels", callback_data, 
        total_reads, n_consumed);
      } catch (...) {
	delete parser;
        throw;
      }
    }
  }

  // @cswelcher TODO: check that deallocate LabelPtrMap is correct
  delete parser;
}

// @cswelcher: double-check -- is it valid to pull the address from a reference?
void LabelHash::link_tag_and_label(HashIntoType& kmer, Label& kmer_label) {
  tag_labels.insert(TagLabelPtrPair(kmer, &kmer_label));
  label_tag_ptrs.insert(LabelTagPtrPair(kmer_label, &kmer));
}

void LabelHash::consume_sequence_and_tag_with_labels(const std::string& seq,
					unsigned long long& n_consumed,
					Label& current_label,
					SeenSet * found_tags)
  {
    bool is_new_kmer;
    bool kmer_tagged;

    KMerIterator kmers(seq.c_str(), _ksize);
    HashIntoType kmer;

    unsigned int since = _tag_density / 2 + 1;

    while(!kmers.done()) {
      kmer = kmers.next();

      if ((is_new_kmer = test_and_set_bits( kmer )))
        ++n_consumed;

  #if (1)
      if (is_new_kmer) {
        ++since;
      } else {
        ACQUIRE_ALL_TAGS_SPIN_LOCK
        kmer_tagged = set_contains(all_tags, kmer);
        RELEASE_ALL_TAGS_SPIN_LOCK
        if (kmer_tagged) {
	      since = 1;
	      
	      // Labeling code
	      // TODO: MAKE THREADSAFE!
	      
	      if (!_cmap_contains_label(tag_labels, kmer, current_label)) {
	        ACQUIRE_TAG_COLORS_SPIN_LOCK
	        link_tag_and_label(kmer, current_label);
	        RELEASE_TAG_COLORS_SPIN_LOCK
	      }
	      if (found_tags) {
	        found_tags->insert(kmer);
	      }
        }  else ++since;
      }
  #else
      if (!is_new_kmer && set_contains(all_tags, kmer)) {
        since = 1;
        if (found_tags) { found_tags->insert(kmer); }
      } else {
        since++;
      }
  #endif
      //
      if (since >= _tag_density) {
        ACQUIRE_ALL_TAGS_SPIN_LOCK
        all_tags.insert(kmer);
        RELEASE_ALL_TAGS_SPIN_LOCK
        
        // Labeling code
        // TODO: MAKE THREADSAFE!
        ACQUIRE_TAG_COLORS_SPIN_LOCK
        link_tag_and_label(kmer, current_label);
        RELEASE_TAG_COLORS_SPIN_LOCK
        
        if (found_tags) { found_tags->insert(kmer); }
        since = 1;
      }

    } // iteration over kmers

    if (since >= _tag_density/2 - 1) {
      ACQUIRE_ALL_TAGS_SPIN_LOCK
      all_tags.insert(kmer);	// insert the last k-mer, too.
      RELEASE_ALL_TAGS_SPIN_LOCK
      
      // Label code: TODO: MAKE THREADSAFE!
      link_tag_and_label(kmer, current_label);
      
      if (found_tags) { found_tags->insert(kmer); }
    }
  }
/*
 * Find all labels associated with the sequence
 * For now, check /every/ k-mer with find_all_tags
 * THIS SUCKS AND IT'S YOUR FAULT @CTB
 */
unsigned int LabelHash::sweep_sequence_for_labels(const std::string& seq,
					LabelPtrSet& found_labels,
					bool break_on_stoptags,
					bool stop_big_traversals) {
					
    SeenSet tagged_kmers;
    //LabelPtrSet found_labels;
    
    HashIntoType kmer_f, kmer_r, kmer;
    
    KMerIterator kmers(seq.c_str(), _ksize);
    std::string kmer_s;
    // keep a list of kmers which have already been traversed
    SeenSet traversed_kmers;
    while (!kmers.done()) {
      kmer = kmers.next();
      kmer_s = _revhash(kmer, _ksize);
      _hash(kmer_s.c_str(), _ksize, kmer_f, kmer_r);
      
      // don't even try traversing from k-mers not in the hashtable
      //traversed_kmers.clear();
      if (get_count(uniqify_rc(kmer_f,kmer_r))) {
        partition->find_all_tags(kmer_f, kmer_r, tagged_kmers,
                   all_tags, break_on_stoptags, stop_big_traversals);
        traverse_labels_and_resolve(tagged_kmers, found_labels);
      }
    }
    return traversed_kmers.size();
}

unsigned int LabelHash::sweep_label_neighborhood(const std::string& seq,
                                                  LabelPtrSet& found_labels,
                                                  unsigned int range,
                                                  bool break_on_stoptags,
                                                  bool stop_big_traversals) {

    SeenSet tagged_kmers;
    unsigned int num_traversed;
    num_traversed = partition->sweep_for_tags(seq, tagged_kmers, all_tags, 
                              range, break_on_stoptags, stop_big_traversals);
    traverse_labels_and_resolve(tagged_kmers, found_labels);
    //printf("range=%u ", range);
    if (range == 0) {
      assert(num_traversed == seq.length()-ksize()+1);
    }
    tagged_kmers.clear();
    return num_traversed;
}

LabelPtrSet LabelHash::get_tag_labels(const HashIntoType& tag) {
  LabelPtrSet labels;
  //unsigned int num_labels;
  _get_tag_labels(tag, tag_labels, labels);
  return labels;
}

TagPtrSet LabelHash::get_label_tags(const Label& label) {
  TagPtrSet tags;
  //unsigned int num_tags;
  _get_tags_from_label(label, label_tag_ptrs, tags);
  return tags;
}

void LabelHash::traverse_labels_and_resolve(const SeenSet& tagged_kmers,
                                              LabelPtrSet& found_labels) {
  
  SeenSet::const_iterator si;
  unsigned int num_labels = 0;
  for (si=tagged_kmers.begin(); si!=tagged_kmers.end(); ++si) {
    HashIntoType tag = *si;
    // get the labels associated with this tag
    num_labels = _get_tag_labels(tag, tag_labels, found_labels);
    if (num_labels > 1) {
      // reconcile labels
      // for now do nothing ha
    }
  }
}


