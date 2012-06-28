#include "khmer.hh"
#include "hashtable.hh"
#include "parsers.hh"
#include "read_parsers.hh"

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

void Hashtable::consume_fasta(const std::string &filename,
			      unsigned int &total_reads,
			      unsigned long long &n_consumed,
			      HashIntoType lower_bound,
			      HashIntoType upper_bound,
			      ReadMaskTable ** orig_readmask,
			      bool update_readmask,
			      CallbackFn callback,
			      void * callback_data)
{
#if (0)
  using namespace khmer:: parsers;
#else
  using namespace khmer:: read_parsers;
#endif

  total_reads = 0;
  n_consumed = 0;

  //IParser* parser = IParser::get_parser(filename);
  IParser *	  parser = 
  IParser::get_parser(
    filename, omp_get_max_threads( ), 2*1024*1024*1024U,
    TraceLogger:: TLVL_DEBUG5
  );

  //
  // readmask stuff: were we given one? do we want to update it?
  // 

  ReadMaskTable * readmask = NULL;
  std::list<unsigned int> masklist;

  if (orig_readmask && *orig_readmask) {
    readmask = *orig_readmask;
  }

  //
  // iterate through the FASTA file & consume the reads.
  //
#pragma omp parallel default( shared ) 
  {
    Read read;
    string currName = "";
    string currSeq = "";
  
#if (1)
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
    //if (0 == (total_reads % 100000))
      trace_logger(
	TraceLogger:: TLVL_DEBUG3,
	"Total number of reads processed: %llu\n",
	(unsigned long long int)total_reads
      );
#endif

    read = parser->get_next_read();
    currSeq = read.seq;
    currName = read.name; 

#if (1)
    // do we want to process it?
    if (!readmask || readmask->get(total_reads)) {

      // yep! process.
      unsigned int this_n_consumed;
      bool is_valid;

//#pragma omp critical (process_read)
      this_n_consumed = check_and_process_read(currSeq,
					       is_valid,
					       lower_bound,
					       upper_bound);

      // was this an invalid sequence -> mark as bad?
      if (!is_valid && update_readmask) {
	if (readmask) {
#pragma omp critical (set_read_mask)
	  readmask->set(total_reads, false);
	} else {
#pragma omp critical (append_to_masklist)
	  masklist.push_back(total_reads);
	}
      } else {		// nope -- count it!
	__sync_add_and_fetch( &n_consumed, this_n_consumed );
      }
    } // check masked read
	       
    // reset the sequence info, increment read number
    __sync_add_and_fetch( &total_reads, 1 );

    // run callback, if specified
    if (total_reads % CALLBACK_PERIOD == 0 && callback) {
      try {
	callback("consume_fasta", callback_data, total_reads, n_consumed);
      } catch (...) {
	throw;
      }
    }
#endif
  } // while reads left for parser

  } // parallel region

  //
  // We've either updated the readmask in place, OR we need to create a
  // new one.
  //

  if (orig_readmask && update_readmask && readmask == NULL) {
    // allocate, fill in from masklist
    readmask = new ReadMaskTable(total_reads);

    list<unsigned int>::const_iterator it;
    for(it = masklist.begin(); it != masklist.end(); ++it) {
      readmask->set(*it, false);
    }
    *orig_readmask = readmask;
  }
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

  // DEBUG
  TraceLogger	trace_logger(
    TraceLogger:: TLVL_DEBUG5,
    "consume_string-%lu.log",
    (unsigned long int)omp_get_thread_num( )
  );

  if (lower_bound == upper_bound && upper_bound == 0) {
    bounded = false;
  }

  // DEBUG
  trace_logger( TraceLogger:: TLVL_DEBUG3, "Starting trace...\n" );

  while(!kmers.done()) {
    kmer = kmers.next();
  
    if (!bounded || (kmer >= lower_bound && kmer < upper_bound)) {

      // DEBUG
      trace_logger(
	TraceLogger:: TLVL_DEBUG4, "Processing kmer: %llu\n",
	(unsigned long long int)kmer
      );

// #pragma omp critical (count_kmer)
      count(kmer);
      n_consumed++;

      // DEBUG
      trace_logger(
	TraceLogger:: TLVL_DEBUG4, "Processed kmer.\n"
      );

    }
  }

  return n_consumed;
}

// vim: set sts=2 sw=2:
