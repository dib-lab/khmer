#include "khmer.hh"
#include "hashtable.hh"
#include "parsers.hh"

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

// TODO? Inline in header.
void
Hashtable::
consume_fasta(
  std:: string const  &filename,
  unsigned int	      &total_reads, unsigned long long	&n_consumed,
  HashIntoType	      lower_bound,  HashIntoType	upper_bound,
  CallbackFn	      callback,	    void *		callback_data
)
{
  using namespace khmer:: read_parsers;

  // TODO: Get defaults from config.
  // Note: Always assume only 1 thread if invoked this way.
  IParser *	  parser = 
  IParser::get_parser(
    filename, 1, 2*1024*1024*1024U, TraceLogger:: TLVL_NONE
  );

  consume_fasta(
    parser, 
    total_reads, n_consumed, 
    lower_bound, upper_bound, 
    callback, callback_data
  );
}

void
Hashtable::
consume_fasta(
  read_parsers:: IParser *	      parser,
  unsigned int	      &total_reads, unsigned long long  &n_consumed,
  HashIntoType	      lower_bound,  HashIntoType	upper_bound,
  CallbackFn	      callback,	    void *		callback_data
)
{
  using namespace khmer:: read_parsers;

  Hasher		  &hasher		= _get_hasher( );
  unsigned int		  total_reads_LOCAL	= 0;
  unsigned long long int  n_consumed_LOCAL	= 0;
  Read			  read;
  
  hasher.trace_logger(
    TraceLogger:: TLVL_DEBUG2, "Starting trace of 'consume_fasta'....\n"
  );

  // Iterate through the reads and consume their kmers.
  while (!parser->is_complete( ))
  {
    unsigned int  this_n_consumed;
    bool	  is_valid;

    read      = parser->get_next_read();

    this_n_consumed = 
    check_and_process_read(read.sequence, is_valid, lower_bound, upper_bound);

    n_consumed_LOCAL  = __sync_add_and_fetch( &n_consumed, this_n_consumed );
    total_reads_LOCAL = __sync_add_and_fetch( &total_reads, 1 );

    if (0 == (total_reads_LOCAL % 10000))
      hasher.trace_logger(
	TraceLogger:: TLVL_DEBUG3,
	"Total number of reads processed: %llu\n",
	(unsigned long long int)total_reads_LOCAL
      );

    // run callback, if specified
    if (callback && (0 == (total_reads_LOCAL % CALLBACK_PERIOD)))
    {
      try 
      {
	callback(
	  "consume_fasta", callback_data,
	  total_reads_LOCAL, n_consumed_LOCAL
	);
      }
      catch (...) { throw; }
    }

  } // while reads left for parser

} // consume_fasta

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

  if (lower_bound == upper_bound && upper_bound == 0) {
    bounded = false;
  }

  while(!kmers.done()) {
    kmer = kmers.next();
  
    if (!bounded || (kmer >= lower_bound && kmer < upper_bound)) {

      count(kmer);
      n_consumed++;

    }
  }

  return n_consumed;
}

// vim: set sts=2 sw=2:
