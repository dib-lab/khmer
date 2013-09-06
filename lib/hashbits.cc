//
// This file is part of khmer, http://github.com/ged-lab/khmer/, and is
// Copyright (C) Michigan State University, 2009-2013. It is licensed under
// the three-clause BSD license; see doc/LICENSE.txt. Contact: ctb@msu.edu
//

#include <iostream>
#include "hashtable.hh"
#include "hashbits.hh"
#include "read_parsers.hh"

using namespace std;
using namespace khmer;
using namespace khmer:: read_parsers;

void Hashbits::save(std::string outfilename)
{
  assert(_counts[0]);

  unsigned int save_ksize = _ksize;
  unsigned char save_n_tables = _n_tables;
  unsigned long long save_tablesize;

  ofstream outfile(outfilename.c_str(), ios::binary);

  unsigned char version = SAVED_FORMAT_VERSION;
  outfile.write((const char *) &version, 1);

  unsigned char ht_type = SAVED_HASHBITS;
  outfile.write((const char *) &ht_type, 1);

  outfile.write((const char *) &save_ksize, sizeof(save_ksize));
  outfile.write((const char *) &save_n_tables, sizeof(save_n_tables));

  for (unsigned int i = 0; i < _n_tables; i++) {
    save_tablesize = _tablesizes[i];
    unsigned long long tablebytes = save_tablesize / 8 + 1;

    outfile.write((const char *) &save_tablesize, sizeof(save_tablesize));

    outfile.write((const char *) _counts[i], tablebytes);
  }
  outfile.close();
}

void Hashbits::load(std::string infilename)
{
  if (_counts) {
    for (unsigned int i = 0; i < _n_tables; i++) {
      delete _counts[i]; _counts[i] = NULL;
    }
    delete _counts; _counts = NULL;
  }
  _tablesizes.clear();
  
  unsigned int save_ksize = 0;
  unsigned char save_n_tables = 0;
  unsigned long long save_tablesize = 0;
  unsigned char version, ht_type;

  ifstream infile(infilename.c_str(), ios::binary);
  assert(infile.is_open());

  infile.read((char *) &version, 1);
  infile.read((char *) &ht_type, 1);
  assert(version == SAVED_FORMAT_VERSION);
  assert(ht_type == SAVED_HASHBITS);

  infile.read((char *) &save_ksize, sizeof(save_ksize));
  infile.read((char *) &save_n_tables, sizeof(save_n_tables));

  _ksize = (WordLength) save_ksize;
  _n_tables = (unsigned int) save_n_tables;
  _init_bitstuff();

  _counts = new Byte*[_n_tables];
  for (unsigned int i = 0; i < _n_tables; i++) {
    HashIntoType tablesize;
    unsigned long long tablebytes;

    infile.read((char *) &save_tablesize, sizeof(save_tablesize));

    tablesize = (HashIntoType) save_tablesize;
    _tablesizes.push_back(tablesize);

    tablebytes = tablesize / 8 + 1;
    _counts[i] = new Byte[tablebytes];

    unsigned long long loaded = 0;
    while (loaded != tablebytes) {
      infile.read((char *) _counts[i], tablebytes - loaded);
      loaded += infile.gcount();	// do I need to do this loop?
    }
  }
  infile.close();
}

// for counting overlap k-mers specifically!!

//
// check_and_process_read: checks for non-ACGT characters before consuming
//

unsigned int Hashbits::check_and_process_read_overlap(std::string &read,
					    bool &is_valid,
                                            Hashbits &ht2)
{
   is_valid = check_and_normalize_read(read);

   if (!is_valid) { return 0; }

   return consume_string_overlap(read, ht2);
}

/*
 * Pretty much copy-pasta from the above functions
 * Might be time for a refactor: could do a general consume_fasta
 * function which accepts a consume_sequence function pointer as a parameter
 */

void
Hashbits::
consume_fasta_and_tag_with_colors(
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


  consume_fasta_and_tag_with_colors(
    parser,
    total_reads, n_consumed,
    callback, callback_data
  );

  delete parser;
}

void
Hashbits::
consume_fasta_and_tag_with_colors(
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
    "Starting trace of 'consume_fasta_and_tag'....\n"
  );

  // Iterate through the reads and consume their k-mers.
  while (!parser->is_complete( ))
  {
    unsigned long long this_n_consumed   = 0;

    read = parser->get_next_read( );

    if (check_and_normalize_read( read.sequence ))
    {
      // TODO: make threadsafe!
      consume_sequence_and_tag_with_colors( read.sequence,
					    this_n_consumed,
					    _tag_color );
      ++_tag_color;

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
	  callback("consume_fasta_and_tag", callback_data, total_reads_TL,
		   n_consumed);
	} catch (...) {
	  delete parser;
	  throw;
	}
      }
#endif // 0

  } // while reads left for parser

}


/* This is essentially the same code as above, only it assigns colors to the
 * tags through multimap TagColorMap defined in hashtable.hh, declared in
 * hashbits.hh
 */
void Hashbits::consume_sequence_and_tag_with_colors(const std::string& seq,
					unsigned long long& n_consumed,
					Color& current_color,
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
	    
	    // Coloring code
	    // TODO: MAKE THREADSAFE!
	    
	    if (!_map_contains(color_map, kmer, current_color)) {
	      color_map.insert(TagColorPair(kmer, current_color))
	    }
	    if (found_tags) {
	      found_tags->insert(kmer);
	    }
      }  else ++since;
    }
    // Should I bother adding new code down here?
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
      
      // Coloring code
      // TODO: MAKE THREADSAFE!
      color_map.insert(TagColorPair(kmer, current_color))
      
      if (found_tags) { found_tags->insert(kmer); }
      since = 1;
    }

  } // iteration over kmers

  if (since >= _tag_density/2 - 1) {
    ACQUIRE_ALL_TAGS_SPIN_LOCK
    all_tags.insert(kmer);	// insert the last k-mer, too.
    RELEASE_ALL_TAGS_SPIN_LOCK
    
    // Color code: TODO: MAKE THREADSAFE!
    color_map.insert(TagColorPair(kmer, current_color))
    
    if (found_tags) { found_tags->insert(kmer); }
  }
}


//
// consume_fasta: consume a FASTA file of reads
//

void Hashbits::consume_fasta_overlap(const std::string &filename,
                              HashIntoType curve[2][100],Hashbits &ht2,
			      unsigned int &total_reads,
			      unsigned long long &n_consumed,
			      CallbackFn callback,
			      void * callback_data)
{
  total_reads = 0;
  n_consumed = 0;
  Read read;

//get total number of reads in dataset

  IParser* parser = IParser::get_parser(filename.c_str());
  while(!parser->is_complete())  {
    read = parser->get_next_read();
    total_reads++;
  }
//block size for curve
  int block_size = total_reads/100;
  
  total_reads = 0;
  khmer::HashIntoType start = 0, stop = 0;
  
  delete parser;
  parser = IParser::get_parser(filename.c_str());



  string currName = "";
  string currSeq = "";

  //
  // iterate through the FASTA file & consume the reads.
  //

  while(!parser->is_complete())  {
    read = parser->get_next_read();
    currSeq = read.sequence;
    currName = read.name; 

      unsigned int this_n_consumed;
      bool is_valid;

      this_n_consumed = check_and_process_read_overlap(currSeq,
					       is_valid, ht2);

        n_consumed += this_n_consumed;
	       
    // reset the sequence info, increment read number

    total_reads++;

    if (total_reads%block_size == 0) {
        curve[0][total_reads/block_size-1] = n_overlap_kmers(start,stop);
        curve[1][total_reads/block_size-1] = n_kmers(start,stop);
    }
    // run callback, if specified
    if (total_reads % CALLBACK_PERIOD == 0 && callback) {
      try {
        callback("consume_fasta", callback_data, total_reads, n_consumed);
      } catch (...) {
        throw;
      }
    }

  } // while
  
  delete parser;
}

//
// consume_string: run through every k-mer in the given string, & hash it.
//

unsigned int Hashbits::consume_string_overlap(const std::string &s,
					      Hashbits &ht2)
{
  const char * sp = s.c_str();
  unsigned int n_consumed = 0;

  KMerIterator kmers(sp, _ksize);
  HashIntoType kmer;

  while(!kmers.done()) {
    kmer = kmers.next();
  
    count_overlap(kmer,ht2);
    n_consumed++;
  }

  return n_consumed;
}

// vim: set sts=2 sw=2:
