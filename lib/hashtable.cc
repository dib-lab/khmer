#include <iostream>
#include <list>
#include <queue>

#include "khmer.hh"
#include "hashtable.hh"
#include "parsers.hh"

#define IO_BUF_SIZE 50*1000*1000
#define CALLBACK_PERIOD 10000
#define PARTITION_ALL_TAG_DEPTH 500

using namespace khmer;
using namespace std;

MinMaxTable * Hashtable::fasta_file_to_minmax(const std::string &inputfile,
					      unsigned int total_reads,
					      ReadMaskTable * readmask,
					      CallbackFn callback,
					      void * callback_data)
{
   IParser* parser = IParser::get_parser(inputfile.c_str());
   Read read;
   string seq = "";
   unsigned int read_num = 0;

   MinMaxTable * mmt = new MinMaxTable(total_reads);

   while(!parser->is_complete()) {
     read = parser->get_next_read();
     seq = read.seq;

     bool valid_read = true;
     if (!readmask || readmask->get(read_num)) {
       valid_read = check_read(seq);

       if (valid_read) {
         BoundedCounterType minval = get_min_count(seq);
	 BoundedCounterType maxval = get_max_count(seq);

	 mmt->add_min(read_num, minval);
	 mmt->add_max(read_num, maxval);
       }
     }

     seq.clear();
     read_num += 1;

     // run callback, if specified
     if (read_num % CALLBACK_PERIOD == 0 && callback) {
       try {
         callback("fasta_file_to_minmax", callback_data, read_num, 0);
       } catch (...) {
         delete mmt;
	 throw;
       }
     }
   }

  return mmt;
}

//
// filter_fasta_file_any: filters a FASTA file based on whether any (at
// least one) k-mer in a sequence has 'threshold' counts in the hashtable.
//

ReadMaskTable * Hashtable::filter_fasta_file_any(MinMaxTable &minmax,
						 BoundedCounterType threshold,
						 ReadMaskTable * old_readmask,
						 CallbackFn callback,
						 void * callback_data)

{
   unsigned int read_num;
   const unsigned int tablesize = minmax.get_tablesize();
   ReadMaskTable * readmask = new ReadMaskTable(tablesize);

   if (old_readmask) {
     readmask->merge(*old_readmask);
   }

   for (read_num = 0; read_num < tablesize; read_num++) {
     if (readmask->get(read_num)) {
       BoundedCounterType maxval = minmax.get_max(read_num);

       if (maxval < threshold) {
	 readmask->set(read_num, false);
       }

       // run callback, if specified
       if (read_num % CALLBACK_PERIOD == 0 && callback) {
	 try {
	   callback("filter_fasta_file_any", callback_data, read_num, 0);
	 } catch (...) {
	   delete readmask;
	   throw;
	 }
       }
     }
   }
  
   return readmask;
}

//
// filter_fasta_file_limit_n: filters a FASTA file based on whether
// a read has at least n kmers that meet a 'threshold' count
//

ReadMaskTable * Hashtable::filter_fasta_file_limit_n(const std::string &readsfile,
                                                     MinMaxTable &minmax,
                                                     BoundedCounterType threshold,
                                                     BoundedCounterType n,
                                                     ReadMaskTable * old_readmask,
                                                     CallbackFn callback,
                                                     void * callback_data)
{
   IParser* parser = IParser::get_parser(readsfile.c_str());
   string seq;
   Read read;
   unsigned int read_num = 0;
   const unsigned int tablesize = minmax.get_tablesize();

   ReadMaskTable * readmask = new ReadMaskTable(tablesize);

   if (old_readmask) {
     readmask->merge(*old_readmask);
   }

   while(!parser->is_complete()) {
      read = parser->get_next_read();
      seq = read.seq;
     
      if (readmask->get(read_num)) {
         int numPos = seq.length() - _ksize + 1;
         unsigned int n_met = 0;
 
         for (int i = 0; i < numPos; i++)  {
            string kmer = seq.substr(i, _ksize);
            if ((int)this->get_count(kmer.c_str()) >= threshold)  {
               n_met++;
            }
         }
 
         if (n_met < n)  {
            readmask->set(read_num, false);
         }
 
         read_num++;
 
         // run callback, if specified
         if (read_num % CALLBACK_PERIOD == 0 && callback) {
            try {
               callback("filter_fasta_file_limit_n", callback_data, read_num, 0);
            } catch (...) {
               delete readmask;
               throw;
            }
         }
      }
   }
         
   return readmask;
}

//
// filter_fasta_file_all: filters a FASTA file based on whether all
// k-mers in a sequence have 'threshold' counts in the hashtable.
//

ReadMaskTable * Hashtable::filter_fasta_file_all(MinMaxTable &minmax,
						 BoundedCounterType threshold,
						 ReadMaskTable * old_readmask,
						 CallbackFn callback,
						 void * callback_data)
{
   unsigned int read_num;
   const unsigned int tablesize = minmax.get_tablesize();

   ReadMaskTable * readmask = new ReadMaskTable(tablesize);

   if (old_readmask) {
     readmask->merge(*old_readmask);
   }

   for (read_num = 0; read_num < tablesize; read_num++) {
     if (readmask->get(read_num)) {
       BoundedCounterType minval = minmax.get_min(read_num);

       if (minval < threshold) {
	 readmask->set(read_num, false);
       }

       // run callback, if specified
       if (read_num % CALLBACK_PERIOD == 0 && callback) {
	 try {
	   callback("filter_fasta_file_all", callback_data, read_num, 0);
	 } catch (...) {
	   delete readmask;
	   throw;
	 }
       }
     }
   }
  
   return readmask;
}

//
// filter_fasta_file_run: filters a FASTA file based on whether a run of
// k-mers in a sequence al have 'threshold' counts in the hashtable.
//

ReadMaskTable * Hashtable::filter_fasta_file_run(const std::string &inputfile,
						 unsigned int total_reads,
						 BoundedCounterType threshold,
						 unsigned int runlength,
						 ReadMaskTable * old_readmask,
						 CallbackFn callback,
						 void * callback_data)

{
   IParser* parser = IParser::get_parser(inputfile.c_str());
   string seq;
   Read read;
   unsigned int read_num = 0;
   unsigned int n_kept = 0;
   ReadMaskTable * readmask = new ReadMaskTable(total_reads);

   if (old_readmask) {
     readmask->merge(*old_readmask);
   }

   while(parser->is_complete()) {
      read = parser->get_next_read();
      seq = read.seq;

      if (readmask->get(read_num)) {
         bool keep = false;
	   
         const unsigned int length = seq.length();
         const char * s = seq.c_str();
         unsigned int this_run = 0;

         for (unsigned int i = 0; i < length - _ksize + 1; i++) {
            HashIntoType count = this->get_count(s);
	    this_run++;
	    if (count < threshold) {
	       this_run = 0;
	    } else if (this_run >= runlength) {
	       keep = true;
	       break;
	    }
	    s++;
         }

         if (!keep) {
            readmask->set(read_num, false);
         } else {
            n_kept++;
         }
      }

      seq.clear();

      read_num += 1;

      // run callback, if specified
      if (read_num % CALLBACK_PERIOD == 0 && callback) {
         try {
            callback("filter_fasta_file_run", callback_data, read_num,n_kept);
         } catch (...) {
            throw;
         }
      }
   }

   return readmask;
}

///
/// output_fasta_kmer_pos_freq: outputs the kmer frequencies for each read
///

void Hashtable::output_fasta_kmer_pos_freq(const std::string &inputfile,
                                           const std::string &outputfile)
{
   IParser* parser = IParser::get_parser(inputfile.c_str());
   ofstream outfile;
   outfile.open(outputfile.c_str());
   string seq;
   Read read;

   while(!parser->is_complete()) {
      read = parser->get_next_read();
      seq = read.seq;

      int numPos = seq.length() - _ksize + 1;

      for (int i = 0; i < numPos; i++)  {
         string kmer = seq.substr(i, _ksize);
         outfile << (int)get_count(kmer.c_str()) << " ";
      }
      outfile << endl;
   }

   outfile.close();
}


unsigned int khmer::output_filtered_fasta_file(const std::string &inputfile,
					       const std::string &outputfile,
					       ReadMaskTable * readmask,
					       CallbackFn callback,
					       void * callback_data)
{
   IParser* parser = IParser::get_parser(inputfile.c_str());
   ofstream outfile;
   outfile.open(outputfile.c_str());
   Read read;
   string name;
   string seq;
   unsigned int n_kept = 0;
   unsigned int read_num = 0;


   while(!parser->is_complete()) {
      read = parser->get_next_read();

      seq = read.seq;
      name = read.name;

      if (readmask->get(read_num)) {
         outfile << ">" << name << endl;
	 outfile << seq << endl;

	 n_kept++;
      }

      name.clear();
      seq.clear();

      read_num++;

      // run callback, if specified
      if (read_num % CALLBACK_PERIOD == 0 && callback) {
         try {
            callback("output_filtered_fasta_file", callback_data,
		      read_num, n_kept);
	 } catch (...) {
	     outfile.close();
	     throw;
	 }
      }
   }
  
   outfile.close();
   return n_kept;
}

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

bool Hashtable::check_read(const std::string &read)
{
   unsigned int i;
   bool is_valid = true;

   if (read.length() < _ksize) {
     is_valid = false;
     return 0;
   }

   for (i = 0; i < read.length(); i++)  {
     if (!is_valid_dna(read[i])) {
         is_valid = false;
	 return 0;
      }
   }

   return is_valid;
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
      bin = bin % _tablesize;
      if (_counts[bin] != MAX_COUNT) {
	_counts[bin]++;
      }
      n_consumed++;
    }

    for (unsigned int i = _ksize; i < length; i++) {
      // left-shift the previous hash over
      h = h << 2;

      // 'or' in the current nt
      h |= twobit_repr(sp[i]);

      // mask off the 2 bits we shifted over.
      h &= bitmask;

      // now handle reverse complement
      r = r >> 2;
      r |= (twobit_comp(sp[i]) << (_ksize*2 - 2));

      bin = uniqify_rc(h, r);

      if (!bounded || (bin >= lower_bound && bin < upper_bound)) {
	bin = bin % _tablesize;
	if (_counts[bin] != MAX_COUNT) {
	  _counts[bin]++;
	}
	n_consumed++;
      }
    }
  } catch (...) {
    throw;
  }

  return n_consumed;
}


BoundedCounterType Hashtable::get_min_count(const std::string &s,
					    HashIntoType lower_bound,
					    HashIntoType upper_bound)
{
  const unsigned int length = s.length();
  const char * sp = s.c_str();
  BoundedCounterType min_count = MAX_COUNT, count;

  HashIntoType h = 0, r = 0;
  bool bounded = true;

  if (lower_bound == upper_bound && upper_bound == 0) {
    bounded = false;
  }

  HashIntoType bin;
  
  bin = _hash(sp, _ksize, h, r);
  if (!bounded || (bin >= lower_bound && bin < upper_bound)) {
    min_count = this->get_count(bin);
  }

  for (unsigned int i = _ksize; i < length; i++) {
    // left-shift the previous hash over
    h = h << 2;

    // 'or' in the current nt
    h |= twobit_repr(sp[i]);

    // mask off the 2 bits we shifted over.
    h &= bitmask;

    // now handle reverse complement
    r = r >> 2;
    r |= (twobit_comp(sp[i]) << (_ksize*2 - 2));

    bin = uniqify_rc(h, r);

    if (!bounded || (bin >= lower_bound && bin < upper_bound)) {
      count = this->get_count(bin);
    
      if (count < min_count) {
	min_count = count;
      }
    }
  }
  return min_count;
}

BoundedCounterType Hashtable::get_max_count(const std::string &s,
					    HashIntoType lower_bound,
					    HashIntoType upper_bound)
{
  const unsigned int length = s.length();
  const char * sp = s.c_str();
  BoundedCounterType max_count = 0, count;

  HashIntoType h = 0, r = 0;
  bool bounded = true;

  if (lower_bound == upper_bound && upper_bound == 0) {
    bounded = false;
  }

  HashIntoType bin = _hash(sp, _ksize, h, r);
  if (!bounded || (bin >= lower_bound && bin < upper_bound)) {
    max_count = this->get_count(bin);
  }

  for (unsigned int i = _ksize; i < length; i++) {
    // left-shift the previous hash over
    h = h << 2;

    // 'or' in the current nt
    h |= twobit_repr(sp[i]);

    // mask off the 2 bits we shifted over.
    h &= bitmask;

    // now handle reverse complement
    r = r >> 2;
    r |= (twobit_comp(sp[i]) << (_ksize*2-2));

    bin = uniqify_rc(h, r);
    if (!bounded || (bin >= lower_bound && bin < upper_bound)) {
      count = this->get_count(bin);

      if (count > max_count) {
	max_count = count;
      }
    }
  }
  return max_count;
}

HashIntoType * Hashtable::abundance_distribution() const
{
  HashIntoType * dist = new HashIntoType[256];
  HashIntoType i;
  
  for (i = 0; i < 256; i++) {
    dist[i] = 0;
  }

  for (i = 0; i < _tablesize; i++) {
    dist[_counts[i]]++;
  }

  return dist;
}

HashIntoType * Hashtable::fasta_count_kmers_by_position(const std::string &inputfile,
					     const unsigned int max_read_len,
					     ReadMaskTable * readmask,
					     BoundedCounterType limit_by_count,
					     CallbackFn callback,
					     void * callback_data)
{
   unsigned long long *counts = new unsigned long long[max_read_len];

   for (unsigned int i = 0; i < max_read_len; i++) {
     counts[i] = 0;
   }

   Read read;
   IParser* parser = IParser::get_parser(inputfile.c_str());
   string name;
   string seq;
   unsigned int read_num = 0;

   while(!parser->is_complete()) {
      read = parser->get_next_read();

      seq = read.seq;
      bool valid_read = true;
	 
      if (!readmask || readmask->get(read_num)) {
	valid_read = check_read(seq);

	if (valid_read) {
	  for (unsigned int i = 0; i < seq.length() - _ksize + 1; i++) {
	    string kmer = seq.substr(i, i + _ksize);
	    BoundedCounterType n = get_count(kmer.c_str());
	    
	    if (limit_by_count == 0 || n == limit_by_count) {
	      if (i < max_read_len) {
		counts[i]++;
	      }
	    }
	  }
	}
 
	name.clear();
	seq.clear();
      }

      read_num += 1;

      // run callback, if specified
      if (read_num % CALLBACK_PERIOD == 0 && callback) {
         try {
	    callback("fasta_file_count_kmers_by_position", callback_data, read_num, 0);
         } catch (...) {
	    throw;
         }
      }
   }

   return counts;
}

void Hashtable::fasta_dump_kmers_by_abundance(const std::string &inputfile,
					      ReadMaskTable * readmask,
					      BoundedCounterType limit_by_count,
					      CallbackFn callback,
					      void * callback_data)
{
  Read read;
  IParser* parser = IParser::get_parser(inputfile.c_str());
  string name;
  string seq;
  unsigned int read_num = 0;

  while(!parser->is_complete()) {
    read = parser->get_next_read();
    bool valid_read = true;
    seq = read.seq;

    if (!readmask || readmask->get(read_num)) {
      valid_read = check_read(seq);

      if (valid_read) {
        for (unsigned int i = 0; i < seq.length() - _ksize + 1; i++) {
	  string kmer = seq.substr(i, i + _ksize);
	  BoundedCounterType n = get_count(kmer.c_str());
	  char ss[_ksize + 1];
	  strncpy(ss, kmer.c_str(), _ksize);
	  ss[_ksize] = 0;

	  if (n == limit_by_count) {
            cout << ss << endl;
	  }
	}
      }

      name.clear();
      seq.clear();
    }

    read_num += 1;

    // run callback, if specified
    if (read_num % CALLBACK_PERIOD == 0 && callback) {
      try {
        callback("fasta_file_dump_kmers_by_abundance", callback_data, read_num, 0);
      } catch (...) {
        throw;
      }
    }
  }
}

//////////////////////////////////////////////////////////////////////
// graph stuff

ReadMaskTable * Hashtable::filter_file_connected(const std::string &est,
                                                 const std::string &readsfile,
                                                 unsigned int total_reads)
{
   unsigned int read_num = 0;
   unsigned int n_kept = 0;
   unsigned long long int cluster_size;
   ReadMaskTable * readmask = new ReadMaskTable(total_reads);
   IParser* parser = IParser::get_parser(readsfile.c_str());


   std::string first_kmer = est.substr(0, _ksize);
   SeenSet keeper;
   calc_connected_graph_size(first_kmer.c_str(),
                             cluster_size,
                             keeper);

   while(!parser->is_complete())
   {
      std::string seq = parser->get_next_read().seq;

      if (readmask->get(read_num))
      {
         bool keep = false;

         HashIntoType h = 0, r = 0, kmer;
         kmer = _hash(seq.substr(0, _ksize).c_str(), _ksize, h, r);
         kmer = uniqify_rc(h, r);

         SeenSet::iterator i = keeper.find(kmer);
         if (i != keeper.end()) {
            keep = true;
         }

         if (!keep) {
            readmask->set(read_num, false);
         } else {
            n_kept++;
         }
      }

      read_num++;
   }

   return readmask;
}

void Hashtable::calc_connected_graph_size(const HashIntoType kmer_f,
					  const HashIntoType kmer_r,
					  unsigned long long& count,
					  SeenSet& keeper,
					  const unsigned long long threshold)
const
{
  HashIntoType kmer = uniqify_rc(kmer_f, kmer_r);
  const BoundedCounterType val = _counts[kmer % _tablesize];

  if (val == 0) {
    return;
  }

  // have we already seen me? don't count; exit.
  SeenSet::iterator i = keeper.find(kmer);
  if (i != keeper.end()) {
    return;
  }

  // keep track of both seen kmers, and counts.
  keeper.insert(kmer);
  count += 1;

  // are we past the threshold? truncate search.
  if (threshold && count >= threshold) {
    return;
  }

  // otherwise, explore in all directions.

  // NEXT.

  HashIntoType f, r;
  const unsigned int rc_left_shift = _ksize*2 - 2;

  f = ((kmer_f << 2) & bitmask) | twobit_repr('A');
  r = kmer_r >> 2 | (twobit_comp('A') << rc_left_shift);
  calc_connected_graph_size(f, r, count, keeper, threshold);

  f = ((kmer_f << 2) & bitmask) | twobit_repr('C');
  r = kmer_r >> 2 | (twobit_comp('C') << rc_left_shift);
  calc_connected_graph_size(f, r, count, keeper, threshold);

  f = ((kmer_f << 2) & bitmask) | twobit_repr('G');
  r = kmer_r >> 2 | (twobit_comp('G') << rc_left_shift);
  calc_connected_graph_size(f, r, count, keeper, threshold);

  f = ((kmer_f << 2) & bitmask) | twobit_repr('T');
  r = kmer_r >> 2 | (twobit_comp('T') << rc_left_shift);
  calc_connected_graph_size(f, r, count, keeper, threshold);

  // PREVIOUS.

  r = ((kmer_r << 2) & bitmask) | twobit_comp('A');
  f = kmer_f >> 2 | (twobit_repr('A') << rc_left_shift);
  calc_connected_graph_size(f, r, count, keeper, threshold);

  r = ((kmer_r << 2) & bitmask) | twobit_comp('C');
  f = kmer_f >> 2 | (twobit_repr('C') << rc_left_shift);
  calc_connected_graph_size(f, r, count, keeper, threshold);

  r = ((kmer_r << 2) & bitmask) | twobit_comp('G');
  f = kmer_f >> 2 | (twobit_repr('G') << rc_left_shift);
  calc_connected_graph_size(f, r, count, keeper, threshold);

  r = ((kmer_r << 2) & bitmask) | twobit_comp('T');
  f = kmer_f >> 2 | (twobit_repr('T') << rc_left_shift);
  calc_connected_graph_size(f, r, count, keeper, threshold);
}

void Hashtable::trim_graphs(const std::string infilename,
			    const std::string outfilename,
			    unsigned int min_size,
			    CallbackFn callback,
			    void * callback_data)
{
  IParser* parser = IParser::get_parser(infilename.c_str());
  unsigned int total_reads = 0;
  unsigned int reads_kept = 0;
  Read read;
  string seq;

  string line;
  ofstream outfile(outfilename.c_str());

  //
  // iterate through the FASTA file & consume the reads.
  //

  while(!parser->is_complete())  {
    read = parser->get_next_read();
    seq = read.seq;

    bool is_valid = check_read(seq);

    if (is_valid) {
      std::string first_kmer = seq.substr(0, _ksize);
      unsigned long long clustersize = 0;
      SeenSet keeper;
      calc_connected_graph_size(first_kmer.c_str(), clustersize, keeper,
				min_size);

      if (clustersize >= min_size) {
	outfile << ">" << read.name << endl;
	outfile << seq << endl;
	reads_kept++;
      }
    }
	       
    total_reads++;

    // run callback, if specified
    if (total_reads % CALLBACK_PERIOD == 0 && callback) {
      try {
	callback("trim_graphs", callback_data, total_reads, reads_kept);
      } catch (...) {
	delete parser;
	throw;
      }
    }
  }

  delete parser;
}

HashIntoType * Hashtable::graphsize_distribution(const unsigned int &max_size)
{
  HashIntoType * p = new HashIntoType[max_size];
  const unsigned char seen = 1 << 7;
  unsigned long long size;

  for (unsigned int i = 0; i < max_size; i++) {
    p[i] = 0;
  }

  for (HashIntoType i = 0; i < _tablesize; i++) {
    BoundedCounterType count = _counts[i];
    if (count && !(count & seen)) {
      std::string kmer = _revhash(i, _ksize);
      size = 0;

      SeenSet keeper;
      calc_connected_graph_size(kmer.c_str(), size, keeper, max_size);
      if (size) {
	if (size < max_size) {
	  p[size] += 1;
	}
      }
    }
  }

  return p;
}

void Hashtable::save(std::string outfilename)
{
  assert(_counts);

  unsigned int save_ksize = _ksize;
  unsigned long long save_tablesize = _tablesize;

  ofstream outfile(outfilename.c_str(), ios::binary);

  outfile.write((const char *) &save_ksize, sizeof(save_ksize));
  outfile.write((const char *) &save_tablesize, sizeof(save_tablesize));

  outfile.write((const char *) _counts,
		sizeof(BoundedCounterType) * _tablesize);
  outfile.close();
}

void Hashtable::load(std::string infilename)
{
  if (_counts) { delete _counts; _counts = NULL; }
  
  unsigned int save_ksize = 0;
  unsigned long long save_tablesize = 0;

  ifstream infile(infilename.c_str(), ios::binary);
  infile.read((char *) &save_ksize, sizeof(save_ksize));
  infile.read((char *) &save_tablesize, sizeof(save_tablesize));

  _ksize = (WordLength) save_ksize;
  _tablesize = (HashIntoType) save_tablesize;
  _counts = new BoundedCounterType[_tablesize];

  unsigned long long loaded = 0;
  while (loaded != _tablesize) {
    infile.read((char *) _counts, _tablesize - loaded);
    loaded += infile.gcount();	// do I need to do this loop?
  }
  infile.close();
}

void Hashtable::save_tagset(std::string outfilename)
{
  ofstream outfile(outfilename.c_str(), ios::binary);
  const unsigned int tagset_size = all_tags.size();
  HashIntoType * buf = new HashIntoType[tagset_size];

  outfile.write((const char *) &tagset_size, sizeof(tagset_size));

  unsigned int i = 0;
  for (SeenSet::iterator pi = all_tags.begin(); pi != all_tags.end();
	 pi++, i++) {
    buf[i] = *pi;
  }

  outfile.write((const char *) buf, sizeof(HashIntoType) * tagset_size);
  outfile.close();

  delete buf;
}

void Hashtable::load_tagset(std::string infilename)
{
  ifstream infile(infilename.c_str(), ios::binary);
  all_tags.clear();

  unsigned int tagset_size = 0;
  infile.read((char *) &tagset_size, sizeof(tagset_size));
  HashIntoType * buf = new HashIntoType[tagset_size];

  infile.read((char *) buf, sizeof(HashIntoType) * tagset_size);

  for (unsigned int i = 0; i < tagset_size; i++) {
    all_tags.insert(buf[i]);
  }
  
  delete buf;
}


// do_truncated_partition: progressive partitioning.
//   1) load all sequences, tagging first kmer of each
//   2) do a truncated BFS search for all connected tagged kmers & assign
//      partition ID.
//
// CTB note: for unlimited PARTITION_ALL_TAG_DEPTH, yields perfect clustering.

void Hashtable::do_truncated_partition(const std::string infilename,
				       CallbackFn callback,
				       void * callback_data)
{
  unsigned int total_reads = 0;

  IParser* parser = IParser::get_parser(infilename);
  Read read;
  string seq;
  bool is_valid;

  std::string first_kmer;
  HashIntoType kmer_f, kmer_r;
  SeenSet tagged_kmers;
  bool surrender;

  if (!partition) {
    partition = new SubsetPartition(this);
  }

  while(!parser->is_complete()) {
    // increment read number
    read = parser->get_next_read();
    total_reads++;

    seq = read.seq;

    check_and_process_read(seq, is_valid);

    if (is_valid) {
      first_kmer = seq.substr(0, _ksize);
      _hash(first_kmer.c_str(), _ksize, kmer_f, kmer_r);

      // find all tagged kmers within range.
      tagged_kmers.clear();
      surrender = false;
      partition->find_all_tags(kmer_f, kmer_r, tagged_kmers, surrender, true);

      // assign the partition ID
      partition->assign_partition_id(kmer_f, tagged_kmers, surrender);
      all_tags.insert(kmer_f);

      // run callback, if specified
      if (total_reads % CALLBACK_PERIOD == 0 && callback) {
	try {
	  callback("do_truncated_partition/read", callback_data, total_reads,
		   partition->next_partition_id);
	} catch (...) {
	  delete parser;
	  throw;
	}
      }
    }
  }

  delete parser;
}

void SubsetPartition::count_partitions(unsigned int& n_partitions,
				       unsigned int& n_unassigned,
				       unsigned int& n_surrendered)
{
  n_partitions = 0;
  n_unassigned = 0;
  n_surrendered = 0;

  PartitionSet partitions;

  //
  // now, go through all the reads, and take those with assigned partitions
  // and output them.
  //

  for (PartitionMap::const_iterator pi = partition_map.begin();
       pi != partition_map.end(); pi++) {
    PartitionID * partition_p = pi->second;
    if (partition_p) {
      partitions.insert(*partition_p);

      if (*partition_p == SURRENDER_PARTITION) {
	n_surrendered++;
      }
    }
    else {
      n_unassigned++;
    }
  }
  n_partitions = partitions.size();
}


unsigned int SubsetPartition::output_partitioned_file(const std::string infilename,
						      const std::string outputfile,
						      bool output_unassigned,
						      CallbackFn callback,
						      void * callback_data)
{
  IParser* parser = IParser::get_parser(infilename);
  ofstream outfile(outputfile.c_str());

  unsigned int total_reads = 0;
  unsigned int reads_kept = 0;
  unsigned int n_singletons = 0;

  PartitionSet partitions;

  Read read;
  string seq;

  std::string first_kmer;
  HashIntoType kmer_f, kmer_r;

  const unsigned char ksize = _ht->ksize();

  //
  // go through all the reads, and take those with assigned partitions
  // and output them.
  //

  while(!parser->is_complete()) {
    read = parser->get_next_read();
    seq = read.seq;

    if (_ht->check_read(seq)) {
      first_kmer = seq.substr(0, ksize);
      _hash(first_kmer.c_str(), ksize, kmer_f, kmer_r);

      PartitionID * partition_p = partition_map[kmer_f];
      PartitionID partition_id;
      if (partition_p == NULL ){
	partition_id = 0;
	n_singletons++;
      } else {
	partition_id = *partition_p;
	partitions.insert(partition_id);
      }

      // Is this a partition that has not been entirely explored? If so,
      // mark it.
      char surrender_flag = ' ';
      if (partition_id == SURRENDER_PARTITION) {
	surrender_flag = '*';
      }

      if (partition_id > 0 || output_unassigned) {
	outfile << ">" << read.name << "\t" << partition_id
		<< surrender_flag << "\n" 
		<< seq << "\n";
      }
	       
      total_reads++;

      // run callback, if specified
      if (total_reads % CALLBACK_PERIOD == 0 && callback) {
	try {
	  callback("output_partitions", callback_data,
		   total_reads, reads_kept);
	} catch (...) {
	  delete parser; parser = NULL;
	  outfile.close();
	  throw;
	}
      }
    }
  }

  delete parser; parser = NULL;

  return partitions.size() + n_singletons;
}

///

bool SubsetPartition::_do_continue(const HashIntoType kmer,
				   const SeenSet& keeper)
{
  // have we already seen me? don't count; exit.
  SeenSet::iterator i = keeper.find(kmer);

  return (i == keeper.end());
}

bool SubsetPartition::_is_tagged_kmer(const HashIntoType kmer_f,
				      const HashIntoType kmer_r,
				      HashIntoType& tagged_kmer)
{
  SeenSet * tags = &_ht->all_tags;
  SeenSet::const_iterator fi = tags->find(kmer_f);
  if (fi != tags->end()) {
    tagged_kmer = kmer_f;
    return true;
  }

  fi = tags->find(kmer_r);
  if (fi != tags->end()) {
    tagged_kmer = kmer_r;
    return true;
  }

  return false;
}
		     

// used by do_truncated_partition

void SubsetPartition::find_all_tags(HashIntoType kmer_f,
				    HashIntoType kmer_r,
				    SeenSet& tagged_kmers,
				    bool& surrender,
				    bool do_initial_check)
{
  const HashIntoType bitmask = _ht->bitmask;

  HashIntoType tagged_kmer;
  if (do_initial_check && _is_tagged_kmer(kmer_f, kmer_r, tagged_kmer)) {
    if (partition_map[kmer_f] != NULL) {
      tagged_kmers.insert(tagged_kmer); // this might connect kmer_r and kmer_f
      return;
    }
  }

  HashIntoType f, r;
  bool first = true;
  NodeQueue node_q;
  const unsigned int rc_left_shift = _ht->ksize()*2 - 2;
  unsigned int total = 0;

  SeenSet keeper;		// keep track of traversed kmers

  // start breadth-first search.

  node_q.push(kmer_f);
  node_q.push(kmer_r);

  while(!node_q.empty()) {
    if (total > PARTITION_ALL_TAG_DEPTH && \
	node_q.size() > PARTITION_ALL_TAG_DEPTH) {
      surrender = true;
      break;
    }

    kmer_f = node_q.front();
    node_q.pop();
    kmer_r = node_q.front();
    node_q.pop();

    HashIntoType kmer = uniqify_rc(kmer_f, kmer_r);
    if (!_do_continue(kmer, keeper)) {
      continue;
    }

    // keep track of seen kmers
    keeper.insert(kmer);
    //    cout << "INSERT: " << _revhash(kmer, _ht->ksize()) << "=" << (int) (_ht->get_count(kmer)) << " xx " << kmer % _ht->n_entries() << " =\n";
    total++;

    // Is this a kmer-to-tag, and have we put this tag in a partition already?
    // Search no further in this direction.
    if (!first && _is_tagged_kmer(kmer_f, kmer_r, tagged_kmer)) {
      tagged_kmers.insert(tagged_kmer);
      first = false;
      continue;
    }

    //
    // Enqueue next set of nodes.
    //

    // NEXT.
    f = ((kmer_f << 2) & bitmask) | twobit_repr('A');
    r = kmer_r >> 2 | (twobit_comp('A') << rc_left_shift);
    if (_ht->get_count(uniqify_rc(f,r))) {
      node_q.push(f); node_q.push(r);
    }

    f = ((kmer_f << 2) & bitmask) | twobit_repr('C');
    r = kmer_r >> 2 | (twobit_comp('C') << rc_left_shift);
    if (_ht->get_count(uniqify_rc(f,r))) {
      node_q.push(f); node_q.push(r);
    }

    f = ((kmer_f << 2) & bitmask) | twobit_repr('G');
    r = kmer_r >> 2 | (twobit_comp('G') << rc_left_shift);
    if (_ht->get_count(uniqify_rc(f,r))) {
      node_q.push(f); node_q.push(r);
    }

    f = ((kmer_f << 2) & bitmask) | twobit_repr('T');
    r = kmer_r >> 2 | (twobit_comp('T') << rc_left_shift);
    if (_ht->get_count(uniqify_rc(f,r))) {
      node_q.push(f); node_q.push(r);
    }

    // PREVIOUS.
    r = ((kmer_r << 2) & bitmask) | twobit_comp('A');
    f = kmer_f >> 2 | (twobit_repr('A') << rc_left_shift);
    if (_ht->get_count(uniqify_rc(f,r))) {
      node_q.push(f); node_q.push(r);
    }

    r = ((kmer_r << 2) & bitmask) | twobit_comp('C');
    f = kmer_f >> 2 | (twobit_repr('C') << rc_left_shift);
    if (_ht->get_count(uniqify_rc(f,r))) {
      node_q.push(f); node_q.push(r);
    }
    
    r = ((kmer_r << 2) & bitmask) | twobit_comp('G');
    f = kmer_f >> 2 | (twobit_repr('G') << rc_left_shift);
    if (_ht->get_count(uniqify_rc(f,r))) {
      node_q.push(f); node_q.push(r);
    }

    r = ((kmer_r << 2) & bitmask) | twobit_comp('T');
    f = kmer_f >> 2 | (twobit_repr('T') << rc_left_shift);
    if (_ht->get_count(uniqify_rc(f,r))) {
      node_q.push(f); node_q.push(r);
    }

    first = false;
  }
}

///////////////////////////////////////////////////////////////////////

void SubsetPartition::do_partition(HashIntoType first_kmer,
				   HashIntoType last_kmer,
				   CallbackFn callback,
				   void * callback_data)
{
  unsigned int total_reads = 0;

  std::string kmer_s;
  HashIntoType kmer_f, kmer_r;
  SeenSet tagged_kmers;
  bool surrender;
  const unsigned char ksize = _ht->ksize();

  SeenSet::const_iterator si, end;

  si = _ht->all_tags.find(first_kmer);
  if (last_kmer) {
    end = _ht->all_tags.find(last_kmer);
  } else {
    end = _ht->all_tags.end();
  }

  for (; si != end; si++) {
    total_reads++;

    kmer_s = _revhash(*si, ksize);
    _hash(kmer_s.c_str(), ksize, kmer_f, kmer_r);

    // find all tagged kmers within range.
    tagged_kmers.clear();
    surrender = false;
    find_all_tags(kmer_f, kmer_r, tagged_kmers, surrender, false);

    // assign the partition ID
    assign_partition_id(kmer_f, tagged_kmers, surrender);

    // run callback, if specified
    if (total_reads % CALLBACK_PERIOD == 0 && callback) {
#if 0 // @CTB
	try {
	  callback("do_subset_partition/read", callback_data, total_reads,
		   next_partition_id);
	} catch (...) {
	  delete parser;
	  throw;
	}
#endif // 0
      }
  }
}

//

PartitionID SubsetPartition::assign_partition_id(HashIntoType kmer_f,
						 SeenSet& tagged_kmers,
						 bool surrender)

{
  PartitionID return_val = 0; 
  PartitionID * pp = NULL;

  for (SeenSet::iterator si = tagged_kmers.begin(); si != tagged_kmers.end();
       si++) {
    PartitionID * pp = partition_map[*si];
    if (pp && *pp == SURRENDER_PARTITION) {
      surrender = true;
    }
  }

  // did we find a tagged kmer?
  if (surrender) {
    PartitionPtrSet * s = reverse_pmap[SURRENDER_PARTITION];
    PartitionPtrSet::iterator si = s->begin();

    PartitionID * surrender_pp = *si;
    

    SeenSet::iterator ii = tagged_kmers.begin();
    for (; ii != tagged_kmers.end(); ii++) {
      PartitionID * pp = partition_map[*ii];
      if (pp) {
	if (*pp != SURRENDER_PARTITION) {
	  _add_partition_ptr(surrender_pp, pp);
	}
      } else {
	partition_map[*ii] = surrender_pp;
      }
    }
    partition_map[kmer_f] = surrender_pp;
  }
  else if (tagged_kmers.size() >= 1) {
    pp = _reassign_partition_ids(tagged_kmers, kmer_f);
    return_val = *pp;
  } else {
    partition_map[kmer_f] = NULL;
    return_val = 0;
  }

  return return_val;
}

PartitionID * SubsetPartition::_reassign_partition_ids(SeenSet& tagged_kmers,
						     const HashIntoType kmer_f)
{
  SeenSet::iterator it = tagged_kmers.begin();
  unsigned int * this_partition_p = NULL;

  // find first assigned partition ID in tagged set
  while (it != tagged_kmers.end()) {
    this_partition_p = partition_map[*it];
    if (this_partition_p != NULL) {
      break;
    }
    it++;
  }

  // none? allocate!
  if (this_partition_p == NULL) {
    this_partition_p = new PartitionID(next_partition_id);
    next_partition_id++;

    PartitionPtrSet * s = new PartitionPtrSet();
    s->insert(this_partition_p);
    reverse_pmap[*this_partition_p] = s;
  }

  it = tagged_kmers.begin();
  for (; it != tagged_kmers.end(); ++it) {
    PartitionID * pp_id = partition_map[*it];
    if (pp_id == NULL) {
      partition_map[*it] = this_partition_p;
    } else {
      if (*pp_id != *this_partition_p) { // join partitions
	_add_partition_ptr(this_partition_p, pp_id);
      }
    }
  }

  assert(this_partition_p != NULL);
  partition_map[kmer_f] = this_partition_p;

  return this_partition_p;
}

///

void SubsetPartition::_add_partition_ptr(PartitionID *orig_pp, PartitionID *new_pp)
{
  PartitionPtrSet * s = reverse_pmap[*orig_pp];
  PartitionPtrSet * t = reverse_pmap[*new_pp];
  reverse_pmap.erase(*new_pp);
    
  for (PartitionPtrSet::iterator pi = t->begin(); pi != t->end(); pi++) {
    PartitionID * iter_pp;
    iter_pp = *pi;

    *iter_pp = *orig_pp;
    s->insert(iter_pp);
  }
  delete t;
}

static void make_partitions_to_tags(PartitionMap& pmap,
				    PartitionsToTagsMap& pttmap)
{
  SeenSet * sp;
  HashIntoType tag;
  PartitionID p;

  for (PartitionMap::const_iterator pi = pmap.begin(); pi != pmap.end();
       pi++) {
    if (pi->second) {
      tag = pi->first;
      p = *(pi->second);

      if (p != SURRENDER_PARTITION) {
	sp = pttmap[p];
	if (sp == NULL) {
	  sp = new SeenSet();
	  pttmap[p] = sp;
	}
	sp->insert(tag);
      }
    }
  }
}

static void del_partitions_to_tags(PartitionsToTagsMap& pttmap)
{
  for (PartitionsToTagsMap::iterator pt = pttmap.begin();
       pt != pttmap.end(); pt++) {
    SeenSet * sp = pt->second;
    if (sp != NULL) {
      delete sp;
      pt->second = NULL;
    }
  }
}

#if 0

static void print_partition_set(PartitionSet& p)
{
  cout << "\tpartition set: ";
  for (PartitionSet::iterator pi = p.begin(); pi != p.end(); pi++) {
    cout << *pi << ", ";
  }
  cout << "\n";
}

static void print_tag_set(SeenSet& p)
{
  cout << "\ttag set: ";
  for (SeenSet::iterator si = p.begin(); si != p.end(); si++) {
    cout << *si << ", ";
  }
  cout << "\n";
}

#endif //0

static void get_new_tags_from_partitions(SeenSet& old_tags,
					 SeenSet& new_tags,
					 PartitionSet& new_partitions,
					 PartitionsToTagsMap& pttm)
{
  for (PartitionSet::const_iterator psi = new_partitions.begin();
       psi != new_partitions.end(); psi++) {
    PartitionID p = *psi;
    SeenSet * s = pttm[p];

    // should only happen for the incompletely traversed sets
    if (s == NULL) {
      assert(p == SURRENDER_PARTITION);
      continue;
    }
    
    for (SeenSet::const_iterator si = s->begin(); si != s->end(); si++) {
      SeenSet::const_iterator test = old_tags.find(*si);
      if (test == old_tags.end()) {
	new_tags.insert(*si);
      }
    }
  }
}

static void get_new_partitions_from_tags(PartitionSet& old_parts,
					 PartitionSet& new_parts,
					 SeenSet& new_tags,
					 PartitionMap& pmap)
{
  for (SeenSet::const_iterator si = new_tags.begin(); si != new_tags.end();
       si++) {
    PartitionID * pp = pmap[*si];
    if (pp != NULL) {
      PartitionID p = *pp;

      PartitionSet::const_iterator test = old_parts.find(p);
      if (test == old_parts.end()) {
	new_parts.insert(p);
      }
    }
  }
}

static void transfer_tags(SeenSet& from, SeenSet& to)
{
  for (SeenSet::const_iterator si = from.begin(); si != from.end(); si++) {
    to.insert(*si);
  }
  from.clear();
}

static void transfer_partitions(PartitionSet& from, PartitionSet& to)
{
  for (PartitionSet::const_iterator pi = from.begin(); pi != from.end();
       pi++) {
    to.insert(*pi);
  }
  from.clear();
}

void SubsetPartition::merge(SubsetPartition * other)
{
  assert (this != other);

  PartitionsToTagsMap subset_pttm, master_pttm;
  PartitionID * pp;
  PartitionID p;

  //
  // Convert the partition maps.
  //

  make_partitions_to_tags(other->partition_map, subset_pttm);
  make_partitions_to_tags(this->partition_map, master_pttm);

  for (PartitionsToTagsMap::iterator ptti = subset_pttm.begin();
       ptti != subset_pttm.end(); ptti++) {
    p = ptti->first;
    SeenSet * these_tags = ptti->second;

    if (these_tags == NULL) {
      continue;
    }

    // this partition has not yet been transferred.  DEAL.

    //
    // Here, we want to get all of the partitions connected to this
    // one in the master map and the subset map -- and do
    // transitively.  Loop until no more are found.
    //

    SeenSet old_tags;
    SeenSet new_tags;

    PartitionSet old_subset_partitions, new_subset_partitions;
    PartitionSet old_master_partitions, new_master_partitions;

    old_subset_partitions.insert(p);
    transfer_tags(*these_tags, new_tags);

    while(new_tags.size()) {
      // first, get partitions (and then tags) for partitions in the *master*
      // that overlap with any of the tags in this subset partition.
      get_new_partitions_from_tags(old_master_partitions,
				   new_master_partitions,
				   new_tags,
				   this->partition_map);
      transfer_tags(new_tags, old_tags);
      get_new_tags_from_partitions(old_tags, new_tags,
				   new_master_partitions, master_pttm);
      transfer_partitions(new_master_partitions, old_master_partitions);

      // ok, now get partitions (and then tags) for partitions in *subset*
      // that overlap with any of the tags from the master.

      get_new_partitions_from_tags(old_subset_partitions,
				   new_subset_partitions,
				   new_tags,
				   other->partition_map);
      transfer_tags(new_tags, old_tags);
      get_new_tags_from_partitions(old_tags, new_tags,
				   new_subset_partitions, subset_pttm);
      transfer_partitions(new_subset_partitions, old_subset_partitions);

      // aaaaaand.... iterate until no more tags show up!
    }

    // Deal with merging incompletely traversed sets...

    bool surrender = false;
    if (old_subset_partitions.find(SURRENDER_PARTITION) != 
	old_subset_partitions.end() ||
	old_master_partitions.find(SURRENDER_PARTITION) != 
	old_master_partitions.end()) {
      surrender = true;
    }

    // All right!  We've now got all the tags that we want to be part of
    // the master map partition; create or merge.

    PartitionPtrSet * pp_set = NULL;

    if (old_master_partitions.size() == 0 && !surrender) {
      pp = this->get_new_partition();
	
      pp_set = new PartitionPtrSet();
      pp_set->insert(pp);
      this->reverse_pmap[*pp] = pp_set;
    } else {
      if (surrender) {
	pp_set = this->reverse_pmap[SURRENDER_PARTITION];
      } else {
	PartitionSet::iterator psi = old_master_partitions.begin();
	pp_set = this->reverse_pmap[*psi];
      }
      assert(pp_set != NULL);
      pp = *(pp_set->begin());
    }

    if (surrender) {
      assert(*pp == SURRENDER_PARTITION);
    }

    // Remove all of the SeenSets in the subset_pttm map for these
    // now-connected partitions.
    for (PartitionSet::iterator psi2 = old_subset_partitions.begin();
	 psi2 != old_subset_partitions.end(); psi2++) {
      SeenSet * sp = subset_pttm[*psi2];
      if (sp != NULL) {
	subset_pttm[*psi2] = NULL;
	delete sp;
      }
    }

    // Remap all of the tags in the master map. This has the side
    // benefit of condensing pretty much everything into a single
    // pointer...

    PartitionPtrSet remove_me;
    for (SeenSet::iterator si = old_tags.begin(); si != old_tags.end(); si++) {
      PartitionID * old_pp = this->partition_map[*si];

      if (old_pp && *old_pp == SURRENDER_PARTITION) {
	assert(surrender);
	continue;
      }

      if (old_pp == NULL) {
	this->partition_map[*si] = pp;
      } else if (old_pp != pp) {
	remove_me.insert(old_pp);
	this->partition_map[*si] = pp;
      }
    }

    // reset the reverse_pmap for this entry, too.
    if (this->reverse_pmap[*pp]->size() > 1) {
      pp_set = new PartitionPtrSet();
      pp_set->insert(pp);
      this->reverse_pmap[*pp] = pp_set;
    }

    // finally: remove remove_me & associated reverse_pmap entries.
    for (PartitionPtrSet::iterator si = remove_me.begin();
	 si != remove_me.end(); si++) {
      PartitionID p = *(*si);
      assert (p != SURRENDER_PARTITION);
      pp_set = this->reverse_pmap[p];
      if (pp_set) {
	delete pp_set;
	this->reverse_pmap.erase(p);
      }
      delete *si;
    }
  }

  del_partitions_to_tags(subset_pttm);
  del_partitions_to_tags(master_pttm);

  // @CTB deal with surrender partitions
  PartitionID * surr_pp = *(this->reverse_pmap[SURRENDER_PARTITION]->begin());

  for (PartitionMap::iterator pi = other->partition_map.begin();
       pi != other->partition_map.end(); pi++) {
    if (pi->second && *(pi->second) == SURRENDER_PARTITION) {
      PartitionID * this_pp = this->partition_map[pi->first];
      if (!this_pp) {
	this->partition_map[pi->first] = surr_pp;
      } else if (this_pp) {
	if (*this_pp != SURRENDER_PARTITION) {
	  this->_add_partition_ptr(surr_pp, this_pp);
	}
      }
    }
  }
}

//
// consume_fasta: consume a FASTA file of reads
//

void Hashtable::consume_fasta_and_tag(const std::string &filename,
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

  //
  // iterate through the FASTA file & consume the reads.
  //

  while(!parser->is_complete())  {
    read = parser->get_next_read();
    seq = read.seq;

    // yep! process.
    unsigned int this_n_consumed = 0;
    bool is_valid;

    this_n_consumed = check_and_process_read(seq, is_valid);
    n_consumed += this_n_consumed;
    if (is_valid) {
      string first_kmer = seq.substr(0, _ksize);
      HashIntoType kmer_f, kmer_r;
      _hash(first_kmer.c_str(), _ksize, kmer_f, kmer_r);

      all_tags.insert(kmer_f);
    }
	       
    // reset the sequence info, increment read number
    total_reads++;

    // run callback, if specified
    if (total_reads % CALLBACK_PERIOD == 0 && callback) {
      try {
        callback("consume_fasta_and_tag", callback_data, total_reads,
		 n_consumed);
      } catch (...) {
	delete parser;
        throw;
      }
    }
  }
}

// load partition maps from & save to disk 

void SubsetPartition::save_partitionmap(string pmap_filename)
{
  ofstream outfile(pmap_filename.c_str(), ios::binary);
  char * buf = NULL;
  buf = new char[IO_BUF_SIZE];
  unsigned int n_bytes = 0;

  HashIntoType * kmer_p = NULL;
  PartitionID * pp;

  PartitionMap::const_iterator pi = partition_map.begin();
  for (; pi != partition_map.end(); pi++) {
    PartitionID p_id;

    HashIntoType kmer = pi->first;
    if (pi->second != NULL) {
      p_id = *(pi->second);

      kmer_p = (HashIntoType *) (buf + n_bytes);
      *kmer_p = kmer;
      n_bytes += sizeof(HashIntoType);
      pp = (PartitionID *) (buf + n_bytes);
      *pp = p_id;
      n_bytes += sizeof(PartitionID);

      if (n_bytes >= IO_BUF_SIZE - sizeof(HashIntoType) - sizeof(PartitionID)) {
	outfile.write(buf, n_bytes);
	n_bytes = 0;
      }
    }
  }
  if (n_bytes) {
    outfile.write(buf, n_bytes);
  }
  outfile.close();

  delete buf;
}
					 
void SubsetPartition::load_partitionmap(string infilename)
{
  ifstream infile(infilename.c_str(), ios::binary);
  char * buf = NULL;
  buf = new char[IO_BUF_SIZE];

  unsigned int n_bytes = 0;
  unsigned int loaded = 0;
  unsigned int remainder;

  assert(infile.is_open());

  PartitionSet partitions;

  HashIntoType * kmer_p = NULL;
  PartitionID * pp = NULL;

  remainder = 0;
  while (!infile.eof()) {
    unsigned int i;

    infile.read(buf + remainder, IO_BUF_SIZE - remainder);
    n_bytes = infile.gcount() + remainder;
    remainder = n_bytes % (sizeof(PartitionID) + sizeof(HashIntoType));
    n_bytes -= remainder;

    for (i = 0; i < n_bytes;) {
      // ignore kmer for this loop.
      i += sizeof(HashIntoType);
      pp = (PartitionID *) (buf + i);
      i += sizeof(PartitionID);

      partitions.insert(*pp);

      loaded++;
    }
    assert(i == n_bytes);
    memcpy(buf, buf + n_bytes, remainder);
  }

  PartitionID max_p_id = 1;
  PartitionPtrMap ppmap;
  for (PartitionSet::const_iterator si = partitions.begin();
      si != partitions.end(); si++) {
    if (*si == 0) {
      continue;
    }

    PartitionID * p = new PartitionID;
    *p = *(si);
    ppmap[*p] = p;

    PartitionPtrSet * s = new PartitionPtrSet();
    s->insert(p);
    reverse_pmap[*p] = s;

    if (max_p_id < *p) {
      max_p_id = *p;
    }
  }
  next_partition_id = max_p_id + 1;

  infile.clear();
  infile.seekg(0, ios::beg);

  remainder = 0;
  while (!infile.eof()) {
    unsigned int i;

    infile.read(buf + remainder, IO_BUF_SIZE - remainder);
    n_bytes = infile.gcount() + remainder;
    remainder = n_bytes % (sizeof(PartitionID) + sizeof(HashIntoType));
    n_bytes -= remainder;

    for (i = 0; i < n_bytes;) {
      kmer_p = (HashIntoType *) (buf + i);
      i += sizeof(HashIntoType);
      pp = (PartitionID *) (buf + i);
      i += sizeof(PartitionID);

      if (*pp == 0) {
	partition_map[*kmer_p] = NULL;
      } else {
	partition_map[*kmer_p] = ppmap[*pp];
      } 
    }
    assert(i == n_bytes);
    memcpy(buf, buf + n_bytes, remainder);
  }

  infile.close();
  delete buf; buf = NULL;
}


void SubsetPartition::_validate_pmap()
{
  // cout << "validating partition_map\n";

  for (PartitionMap::const_iterator pi = partition_map.begin();
       pi != partition_map.end(); pi++) {
    //HashIntoType kmer = (*pi).first;
    PartitionID * pp_id = (*pi).second;

    if (pp_id != NULL) {
      assert(*pp_id >= 1);
      assert(*pp_id < next_partition_id);
    }
  }

  // cout << "validating reverse_pmap -- st 1\n";
  for (ReversePartitionMap::const_iterator ri = reverse_pmap.begin();
       ri != reverse_pmap.end(); ri++) {
    PartitionID p = (*ri).first;
    PartitionPtrSet *s = (*ri).second;

    assert(s != NULL);

    for (PartitionPtrSet::const_iterator si = s->begin(); si != s->end();
	 si++) {
      PartitionID * pp;
      pp = *si;

      assert (p == *pp);
    }
  }
}

void SubsetPartition::_clear_partitions()
{
  for (ReversePartitionMap::iterator ri = reverse_pmap.begin();
       ri != reverse_pmap.end(); ri++) {
    PartitionPtrSet * s = (*ri).second;

    for (PartitionPtrSet::iterator pi = s->begin(); pi != s->end(); pi++) {
      PartitionID * pp = (*pi);
      delete pp;
    }
    delete s;
  }
  partition_map.clear();
  next_partition_id = 1;
}


void Hashtable::divide_tags_into_subsets(unsigned int subset_size,
					 SeenSet& divvy)
{
  unsigned int i = 0;

  for (SeenSet::const_iterator si = all_tags.begin(); si != all_tags.end();
       si++) {
    if (i % subset_size == 0) {
      divvy.insert(*si);
      i = 0;
    }
    i++;
  }
}

void Hashtable::tags_to_map(TagCountMap& tag_map)
{
  for (SeenSet::const_iterator si = all_tags.begin(); si != all_tags.end();
       si++) {
    tag_map[*si] = 0;
  }
  cout << "TM size: " << tag_map.size() << "\n";
}

void SubsetPartition::maxify_partition_size(TagCountMap& tag_map)
{
  PartitionCountMap partition_count;

  for (PartitionMap::const_iterator pi = partition_map.begin();
       pi != partition_map.end(); pi++) {
    PartitionID * pp = pi->second;
    if (pp) {
      partition_count[*pp]++;
    }
  }

  for (PartitionMap::const_iterator pi = partition_map.begin();
       pi != partition_map.end(); pi++) {
    PartitionID * pp = pi->second;
    if (pp) {    
      if (tag_map[pi->first] < partition_count[*pp]) {
	tag_map[pi->first] = partition_count[*pp];
      }
    }
  }
}

void Hashtable::discard_tags(TagCountMap& tag_map,
				   unsigned int threshold)
{
  SeenSet delete_me;

  for (TagCountMap::const_iterator ti = tag_map.begin(); ti != tag_map.end();
       ti++) {
    if (ti->second < threshold) {
      delete_me.insert(ti->first);
    }
  }

  for (SeenSet::const_iterator si = delete_me.begin(); si != delete_me.end();
       si++) {
    tag_map.erase(*si);
  }
}

void SubsetPartition::filter_against_tags(TagCountMap& tag_map)
{
  PartitionMap new_pmap;

  for (TagCountMap::const_iterator ti = tag_map.begin(); ti != tag_map.end();
       ti++) {
    new_pmap[ti->first] = partition_map[ti->first];
  }

  cout << "OLD partition map size: " << partition_map.size() << "\n";
  cout << "NEW partition map size: " << new_pmap.size() << "\n";

  partition_map.swap(new_pmap);
  new_pmap.clear();
}
