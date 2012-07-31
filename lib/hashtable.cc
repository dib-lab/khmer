#include "khmer.hh"
#include "hashtable.hh"
#include "parsers.hh"

#include <omp.h>

// Note: This simple inlined code should be quicker than a call to a 
//	 'toupper' function in the C library.
//	 This should boil down to one integer compare, one branch, and 
//	 maybe one integer subtraction.
//	 This will be potentially called on trillions of bytes, 
//	 so efficiency probably matters some.
#define quick_toupper( c )  (0x61 <= (c) ? (c) - 0x20 : (c))

using namespace khmer;
using namespace std;


Hashtable:: Hasher::
Hasher( uint32_t const thread_id, uint8_t const	trace_level )
: thread_id( thread_id ),
  trace_logger(
    TraceLogger(
      trace_level, "hastable-%lu.log", (unsigned long int)thread_id
    )
  )
{ }


Hashtable:: Hasher::
~Hasher( )
{ }


//
// check_and_process_read: checks for non-ACGT characters before consuming
//

unsigned int Hashtable::check_and_process_read(std::string &read,
					    bool &is_valid,
                                            HashIntoType lower_bound,
                                            HashIntoType upper_bound)
{
   is_valid = check_and_normalize_read(read);

   if (!is_valid) { return 0; }

   return consume_string(read, lower_bound, upper_bound);
}

//
// check_and_normalize_read: checks for non-ACGT characters
//			     converts lowercase characters to uppercase one
// Note: Usually it is desirable to keep checks and mutations separate.
//	 However, in the interests of efficiency (we are potentially working 
//	 with TB of data), a check and mutation have been placed inside the 
//	 same loop. Potentially trillions fewer fetches from memory would 
//	 seem to be a worthwhile goal.
//

bool Hashtable::check_and_normalize_read(std::string &read) const
{
  if (read.length() < _ksize) {
    return false;
  }

  for (unsigned int i = 0; i < read.length(); i++)  {
    read[ i ] = quick_toupper( read[ i ] ); // normalize to uppercase letters
    if (!is_valid_dna(read[i])) {
      return false;
    }
  }

  return true;
}

//
// consume_fasta: consume a FASTA file of reads
//

#define USE_NEW_PARSER	  1

void Hashtable::consume_fasta(const std::string &filename,
			      unsigned int &total_reads,
			      unsigned long long &n_consumed,
			      HashIntoType lower_bound,
			      HashIntoType upper_bound,
			      CallbackFn callback,
			      void * callback_data)
{
#ifndef USE_NEW_PARSER
  using namespace khmer:: parsers;
#else
  using namespace khmer:: read_parsers;
#endif

  total_reads = 0;
  n_consumed = 0;

  // TODO? Lookup parser by filename.
  //	   Use an empty string for STDIN.
  //	   If no parser exists for filename, then create one.
  //	   Else, re-use existing one.
  //	   Alternatively, we might be able to get clever in the Python wrapper.
#ifndef USE_NEW_PARSER
  IParser* parser = IParser::get_parser(filename);
#else
  IParser *	  parser = 
  IParser::get_parser(
    filename, omp_get_max_threads( ), 2*1024*1024*1024U,
    TraceLogger:: TLVL_NONE
  );
#endif

  //
  // iterate through the FASTA file & consume the reads.
  //
#ifdef USE_NEW_PARSER
# pragma omp parallel default( shared ) 
  {
#endif
    Read read;
    string currName = "";
    string currSeq = "";
  
  // TODO: Move TraceLogger instances into Hasher.
#if (0)
  // DEBUG
  TraceLogger	    trace_logger(
    TraceLogger:: TLVL_DEBUG5,
    "consume_fasta-%lu.log",
    (unsigned long int)omp_get_thread_num( )
  );

  // DEBUG
  trace_logger( TraceLogger:: TLVL_DEBUG2, "Starting trace...\n" );
#endif

  while(!parser->is_complete())  {
    
#if (0)
    // DEBUG
    if (0 == (total_reads % 1000))
      trace_logger(
	TraceLogger:: TLVL_DEBUG3,
	"Total number of reads processed: %llu\n",
	(unsigned long long int)total_reads
      );
#endif

    read = parser->get_next_read();
    currSeq = read.seq;
    currName = read.name; 

    unsigned int this_n_consumed;
    bool is_valid;

//#pragma omp critical (process_read)
    this_n_consumed = check_and_process_read(currSeq,
					     is_valid,
					     lower_bound,
					     upper_bound);

    // TODO: Get back return values with atomic operations
    //	     to use in reporting below.
    __sync_add_and_fetch( &n_consumed, this_n_consumed );
    __sync_add_and_fetch( &total_reads, 1 );
    // run callback, if specified
    if (total_reads % CALLBACK_PERIOD == 0 && callback) {
      try {
	callback("consume_fasta", callback_data, total_reads, n_consumed);
      } catch (...) {
	throw;
      }
    }

  } // while reads left for parser

#ifdef USE_NEW_PARSER
  } // parallel region
#endif
}

//
// consume_string: run through every k-mer in the given string, & hash it.
//

unsigned int Hashtable::consume_string(const std::string &s,
				       HashIntoType lower_bound,
				       HashIntoType upper_bound)
{
  const char * sp = s.c_str();
  unsigned int n_consumed = 0;

  bool bounded = true;

  KMerIterator kmers(sp, _ksize);
  HashIntoType kmer;

#if (0)
  // DEBUG
  TraceLogger	trace_logger(
    TraceLogger:: TLVL_DEBUG5,
    "consume_string-%lu.log",
    (unsigned long int)omp_get_thread_num( )
  );
#endif

  if (lower_bound == upper_bound && upper_bound == 0) {
    bounded = false;
  }

#if (0)
  // DEBUG
  trace_logger( TraceLogger:: TLVL_DEBUG3, "Starting trace...\n" );
#endif

  while(!kmers.done()) {
    kmer = kmers.next();
  
    if (!bounded || (kmer >= lower_bound && kmer < upper_bound)) {

#if (0)
      // DEBUG
      trace_logger(
	TraceLogger:: TLVL_DEBUG4, "Processing kmer: %llu\n",
	(unsigned long long int)kmer
      );
#endif

// #pragma omp critical (count_kmer)
      count(kmer);
      n_consumed++;

#if (0)
      // DEBUG
      trace_logger(
	TraceLogger:: TLVL_DEBUG4, "Processed kmer.\n"
      );
#endif

    }
  }

  return n_consumed;
}

// vim: set sts=2 sw=2:
