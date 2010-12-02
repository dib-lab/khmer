#include "khmer.hh"
#include "hashtable.hh"
#include "parsers.hh"

#define CALLBACK_PERIOD 10000

using namespace khmer;
using namespace std;

//
// check_and_process_read: checks for non-ACGT characters before consuming
//

unsigned int Hashtable::check_and_process_read(const std::string &read,
					    bool &is_valid,
                                            HashIntoType lower_bound,
                                            HashIntoType upper_bound)
{
   is_valid = check_read(read);

   if (!is_valid) { return 0; }

   return consume_string(read, lower_bound, upper_bound);
}

//
// check_read: checks for non-ACGT characters
//

bool Hashtable::check_read(const std::string &read) const
{
  if (read.length() < _ksize) {
    return false;
  }

  for (unsigned int i = 0; i < read.length(); i++)  {
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
  total_reads = 0;
  n_consumed = 0;

  IParser* parser = IParser::get_parser(filename.c_str());
  Read read;

  string currName = "";
  string currSeq = "";

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

  while(!parser->is_complete())  {
    read = parser->get_next_read();
    currSeq = read.seq;
    currName = read.name; 

    // do we want to process it?
    if (!readmask || readmask->get(total_reads)) {

      // yep! process.
      unsigned int this_n_consumed;
      bool is_valid;

      this_n_consumed = check_and_process_read(currSeq,
					       is_valid,
					       lower_bound,
					       upper_bound);

      // was this an invalid sequence -> mark as bad?
      if (!is_valid && update_readmask) {
        if (readmask) {
	  readmask->set(total_reads, false);
	} else {
	  masklist.push_back(total_reads);
	}
      } else {		// nope -- count it!
        n_consumed += this_n_consumed;
      }
    }
	       
    // reset the sequence info, increment read number
    total_reads++;

    // run callback, if specified
    if (total_reads % CALLBACK_PERIOD == 0 && callback) {
      try {
        callback("consume_fasta", callback_data, total_reads, n_consumed);
      } catch (...) {
        throw;
      }
    }
  }


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
  const unsigned int length = s.length();
  unsigned int n_consumed = 0;

  HashIntoType h = 0, r = 0;
  bool bounded = true;

  if (lower_bound == upper_bound && upper_bound == 0) {
    bounded = false;
  }
  
  HashIntoType bin = _hash(sp, _ksize, h, r);

  try {
    if (!bounded || (bin >= lower_bound && bin < upper_bound)) {
      count(bin);
      n_consumed++;
    }

    for (unsigned int i = _ksize; i < length; i++) {
      bin = _next_hash(sp[i], h, r);

      if (!bounded || (bin >= lower_bound && bin < upper_bound)) {
	count(bin);
	n_consumed++;
      }
    }
  } catch (...) {
    throw;
  }

  return n_consumed;
}

