#include <iostream>
#include <list>

#include "khmer.hh"
#include "hashtable.hh"
#include "parsers.hh"

#define CALLBACK_PERIOD 10000
#define PARTITION_FIRST_TAG_DEPTH 50
#define PARTITION_ALL_TAG_DEPTH 500
#define PARTITION_MAX_TAG_EXAMINED 1e6

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
       for (unsigned int i = 0; i < seq.length(); i++)  {
         if (!is_valid_dna(seq[i])) {
           valid_read = false;
           break;
         }
       }

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
   unsigned int i;

   is_valid = true;

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

   return consume_string(read, lower_bound, upper_bound);
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

  writelock_acquire();

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
    writelock_release();
    throw;
  }

  writelock_release();

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
         for (unsigned int i = 0; i < seq.length(); i++)  {
	    if (!is_valid_dna(seq[i])) {
	       valid_read = false;
	       break;
	    }
	 }

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
      for (unsigned int i = 0; i < seq.length(); i++)  {
        if (!is_valid_dna(seq[i])) {
	  valid_read = false;
	  break;
	}
      }

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

#if 0
  // @@@
  if (val != MAX_COUNT) {
    _counts[kmer % _tablesize] = val - 1;
  }
#endif // 0

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

    // yep! process.
    bool is_valid;

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

// used by do_partition to explore an entire graph and assign the tagged
// kmers to a particular partition.

void Hashtable::partition_set_id(const HashIntoType kmer_f,
				 const HashIntoType kmer_r,
				 SeenSet& keeper,
				 const unsigned int partition_id,
				 PartitionMap& partition_map)

{
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

    // keep track of seen kmers
    keeper.insert(kmer);

    // Is this a kmer-to-tag, and have we tagged it already?
    PartitionMap::iterator fi = partition_map.find(kmer_f);
    if (fi != partition_map.end()) {
      unsigned int existing = partition_map[kmer_f];
      if (existing != 0) {	// can eliminate once it works :) @CTB
	assert(existing == partition_id);
      } else {
	partition_map[kmer_f] = partition_id;
      }
    }

    fi = partition_map.find(kmer_r);
    if (fi != partition_map.end()) {
      unsigned int existing = partition_map[kmer_r];
      if (existing != 0) {
	assert(existing == partition_id);
      } else {
	partition_map[kmer_r] = partition_id;
      }
    }
  }

  // NEXT.

  HashIntoType f, r;
  const unsigned int rc_left_shift = _ksize*2 - 2;

  f = ((kmer_f << 2) & bitmask) | twobit_repr('A');
  r = kmer_r >> 2 | (twobit_comp('A') << rc_left_shift);
  partition_set_id(f, r, keeper, partition_id, partition_map);

  f = ((kmer_f << 2) & bitmask) | twobit_repr('C');
  r = kmer_r >> 2 | (twobit_comp('C') << rc_left_shift);
  partition_set_id(f, r, keeper, partition_id, partition_map);

  f = ((kmer_f << 2) & bitmask) | twobit_repr('G');
  r = kmer_r >> 2 | (twobit_comp('G') << rc_left_shift);
  partition_set_id(f, r, keeper, partition_id, partition_map);

  f = ((kmer_f << 2) & bitmask) | twobit_repr('T');
  r = kmer_r >> 2 | (twobit_comp('T') << rc_left_shift);
  partition_set_id(f, r, keeper, partition_id, partition_map);

  // PREVIOUS.

  r = ((kmer_r << 2) & bitmask) | twobit_comp('A');
  f = kmer_f >> 2 | (twobit_repr('A') << rc_left_shift);
  partition_set_id(f, r, keeper, partition_id, partition_map);

  r = ((kmer_r << 2) & bitmask) | twobit_comp('C');
  f = kmer_f >> 2 | (twobit_repr('C') << rc_left_shift);
  partition_set_id(f, r, keeper, partition_id, partition_map);

  r = ((kmer_r << 2) & bitmask) | twobit_comp('G');
  f = kmer_f >> 2 | (twobit_repr('G') << rc_left_shift);
  partition_set_id(f, r, keeper, partition_id, partition_map);

  r = ((kmer_r << 2) & bitmask) | twobit_comp('T');
  f = kmer_f >> 2 | (twobit_repr('T') << rc_left_shift);
  partition_set_id(f, r, keeper, partition_id, partition_map);
}

// do_exact_partition: simple partitioning, done once per cluster.
//   1) load in all the sequences, tagging the first kmer of each sequence.
//   2) then, for each tag, explore entre cluster & set tag
//         using partition_set_id.
//
// slow for big clusters, because it has to find all tagged k-mers in each
// cluster.  no provision for giving up and retagging.

unsigned int Hashtable::do_exact_partition(const std::string infilename,
					   CallbackFn callback,
					   void * callback_data)
{
  PartitionMap partition_map;

  unsigned int total_reads = 0;
  unsigned int reads_kept = 0;

  IParser* parser = IParser::get_parser(infilename);
  Read read;
  string seq;

  while(!parser->is_complete())  {
    bool is_valid;
    read = parser->get_next_read();
    seq = read.seq;

    check_and_process_read(seq, is_valid);

    if (is_valid) {
      std::string first_kmer = seq.substr(0, _ksize);
      HashIntoType kmer_f = _hash_forward(first_kmer.c_str(), _ksize);

      partition_map[kmer_f] = 0;
    }
	       
    total_reads++;

    // run callback, if specified
    if (total_reads % CALLBACK_PERIOD == 0 && callback) {
      try {
	callback("do_exact_partition", callback_data, total_reads, reads_kept);
      } catch (...) {
	delete parser;
	throw;
      }
    }
   }

  delete parser;

  // now, build the partition maps
   
  unsigned int next_partition_id = 1;
  for(PartitionMap::iterator i=partition_map.begin();
      i != partition_map.end();
      ++i) {

    HashIntoType kmer_f = (*i).first;
    unsigned int partition_id = (*i).second;

    if (partition_id == 0) {
      HashIntoType kmer_r;
      SeenSet keeper;

      partition_id = next_partition_id;
      next_partition_id++;

      std::string kmer_s = _revhash(kmer_f, _ksize);

      _hash(kmer_s.c_str(), _ksize, kmer_f, kmer_r);
      partition_set_id(kmer_f, kmer_r, keeper, partition_id, partition_map);
    }
  }

  return next_partition_id - 1;
}

// used by do_truncated_partition

void Hashtable::partition_find_first_tag(const HashIntoType kmer_f,
					 const HashIntoType kmer_r,
					 SeenSet& keeper,
					 HashIntoType& tagged_kmer,
					 const PartitionMap& partition_map,
					 bool& done,
					 unsigned int depth)
{
  if (done || depth == 0) return;
  depth -= 1;

  HashIntoType kmer = uniqify_rc(kmer_f, kmer_r);
  if (!_do_continue(kmer, keeper)) {
    return;
  }

  // keep track of seen kmers
  keeper.insert(kmer);

  if (_is_tagged_kmer(kmer_f, kmer_r, partition_map, tagged_kmer)) {
    done = true;
    return;
  }

  // NEXT.

  HashIntoType f, r;
  const unsigned int rc_left_shift = _ksize*2 - 2;

  f = ((kmer_f << 2) & bitmask) | twobit_repr('A');
  r = kmer_r >> 2 | (twobit_comp('A') << rc_left_shift);
  partition_find_first_tag(f, r, keeper, tagged_kmer, partition_map, done, depth);

  f = ((kmer_f << 2) & bitmask) | twobit_repr('C');
  r = kmer_r >> 2 | (twobit_comp('C') << rc_left_shift);
  partition_find_first_tag(f, r, keeper, tagged_kmer, partition_map, done, depth);

  f = ((kmer_f << 2) & bitmask) | twobit_repr('G');
  r = kmer_r >> 2 | (twobit_comp('G') << rc_left_shift);
  partition_find_first_tag(f, r, keeper, tagged_kmer, partition_map, done, depth);

  f = ((kmer_f << 2) & bitmask) | twobit_repr('T');
  r = kmer_r >> 2 | (twobit_comp('T') << rc_left_shift);
  partition_find_first_tag(f, r, keeper, tagged_kmer, partition_map, done, depth);

  // PREVIOUS.

  r = ((kmer_r << 2) & bitmask) | twobit_comp('A');
  f = kmer_f >> 2 | (twobit_repr('A') << rc_left_shift);
  partition_find_first_tag(f, r, keeper, tagged_kmer, partition_map, done, depth);

  r = ((kmer_r << 2) & bitmask) | twobit_comp('C');
  f = kmer_f >> 2 | (twobit_repr('C') << rc_left_shift);
  partition_find_first_tag(f, r, keeper, tagged_kmer, partition_map, done, depth);

  r = ((kmer_r << 2) & bitmask) | twobit_comp('G');
  f = kmer_f >> 2 | (twobit_repr('G') << rc_left_shift);
  partition_find_first_tag(f, r, keeper, tagged_kmer, partition_map, done, depth);

  r = ((kmer_r << 2) & bitmask) | twobit_comp('T');
  f = kmer_f >> 2 | (twobit_repr('T') << rc_left_shift);
  partition_find_first_tag(f, r, keeper, tagged_kmer, partition_map, done, depth);
}

// do_truncated_partition: less truncated progressive partitioning.
//   1) load all sequences, tagging first kmer of each
//   2) do a truncated DFS search for the *first* tagged kmer; assign cluster
//         (partition_find_first_tag, to PARTITION_FIRST_TAG_DEPTH)
//   3) after loading all sequences, do a truncated DFS search for *all* tagged
//         kmers (partition_find_i4, to PARTITION_ALL_TAG_DEPTH),
//         reassigning now-connected clusters.
//
// #2 is technically not required, but for high-abundance clusters
// with closely spaced tags, it should reduce the amount of
// reassignment in #3.  Benchmarking on well-connected data sets is
// needed.  @CTB
//
// CTB note: for unlimited PARTITION_ALL_TAG_DEPTH, yields perfect clustering.

unsigned int Hashtable::do_truncated_partition(const std::string infilename,
					       const std::string outputfile,
					       const unsigned int threshold,
					       CallbackFn callback,
					       void * callback_data)
{
  PartitionMap partition_map;
  ReversePartitionMap rev_pmap;

  unsigned int total_reads = 0;
  unsigned int reads_kept = 0;
  unsigned int next_partition_id = 1;
  SeenSet keeper;
  HashIntoType tagged_kmer;

  IParser* parser = IParser::get_parser(infilename);
  Read read;
  string seq;

  while(!parser->is_complete()) {
    read = parser->get_next_read();
    seq = read.seq;

    bool is_valid;
    check_and_process_read(seq, is_valid);

    if (is_valid) {
      std::string first_kmer = seq.substr(0, _ksize);

      HashIntoType kmer_f, kmer_r;
      _hash(first_kmer.c_str(), _ksize, kmer_f, kmer_r);
      bool done = false;

      keeper.empty();
      tagged_kmer = 0;
      partition_find_first_tag(kmer_f, kmer_r, keeper, tagged_kmer,
			       partition_map, done,
			       PARTITION_FIRST_TAG_DEPTH);

      unsigned int partition_id;
      if (!done) {		// no tagged_kmer found.
	partition_id = next_partition_id;
	next_partition_id++;

	partition_map[kmer_f] = partition_id;

	SeenSet * x = new SeenSet();
	x->insert(kmer_f);
	rev_pmap[partition_id] = x;
      } else {
	// get graph ID of first tagged kmer
	partition_id = partition_map[tagged_kmer];
	partition_map[kmer_f] = partition_id;
	rev_pmap[partition_id]->insert(kmer_f);
      }
	       
      // reset the sequence info, increment read number
      total_reads++;

      // run callback, if specified
      if (total_reads % CALLBACK_PERIOD == 0 && callback) {
	try {
	  callback("do_truncated_partition/read", callback_data, total_reads,
		   reads_kept);
	} catch (...) {
	  delete parser;
	}
      }
    }
  }

  delete parser;

  // join partitions that should overlap, based on discovery of *all* 
  // tagged kmers within given range (partition_find_all_tags).

  PartitionMap::iterator pi;
  unsigned int n = 0;
  std::set<unsigned int> surrender_set;

  for (pi = partition_map.begin(); pi != partition_map.end(); ++pi) {
    n++;
    if (n % 10000 == 0 && callback) {
      callback("do_truncated_partition/dfs", callback_data, n,
	       rev_pmap.size());
    }
    std::string kmer = _revhash((*pi).first, _ksize);
    unsigned int this_pid = (*pi).second;
    
    std::set<unsigned int>::const_iterator ii;
    ii = surrender_set.find(this_pid);
	   
    if (ii != surrender_set.end()) { // we've already surrendered!
      continue;
    }


    HashIntoType kmer_f, kmer_r;
    _hash(kmer.c_str(), _ksize, kmer_f, kmer_r);

    SeenSet keeper;
    SeenSet tagged_kmers;
    bool surrender = false;

    // find all tagged kmers within range.
    unsigned int total = 0;
    partition_find_all_tags(kmer_f, kmer_r, keeper, tagged_kmers,
			    partition_map, true,
			    PARTITION_ALL_TAG_DEPTH, surrender, total);

    if (surrender) {
      std::cout << "SURRENDER on partition: " << this_pid << "\n";
      surrender_set.insert(this_pid);
    }

    // did we find more than one tagged kmer?
    if (tagged_kmers.size() >= 1) {
      // collate the partitions from the tagged kmers
      set<unsigned int> other_partition_ids;
      SeenSet::iterator it = tagged_kmers.begin();
      for (; it != tagged_kmers.end(); ++it) {
	unsigned int id = partition_map[*it];

	if (id != this_pid) {
	  other_partition_ids.insert(id);
	}
      }

      // did we find a different partition id linked to this partition?
      if (other_partition_ids.size()) {

	// yep -- reassign.
	for (set<unsigned int>::iterator si = other_partition_ids.begin();
	     si != other_partition_ids.end(); si++) {
	  SeenSet * x = rev_pmap[*si];
	  for (SeenSet::iterator pi = x->begin(); pi != x->end(); ++pi){
	    partition_map[*pi] = this_pid;
	    rev_pmap[this_pid]->insert(*pi);
	  }

	  rev_pmap.erase(*si);
	  delete x;
	}
      }
    }
  }

  // restart!
  parser = IParser::get_parser(infilename);

  ofstream outfile(outputfile.c_str());

  total_reads = 0;
  std::set<unsigned int> partitions;

  while(!parser->is_complete()) {
    read = parser->get_next_read();
    seq = read.seq;

    bool is_valid;
    check_and_process_read(seq, is_valid);

    if (is_valid) {
      std::string first_kmer = seq.substr(0, _ksize);

      HashIntoType kmer_f, kmer_r;
      _hash(first_kmer.c_str(), _ksize, kmer_f, kmer_r);

      unsigned int partition_id = partition_map[kmer_f];
      unsigned int cluster_size = rev_pmap[partition_id]->size();

      std::set<unsigned int>::const_iterator ii;
      ii = surrender_set.find(partition_id);
	   
      char surrender_flag = ' ';
      if (ii != surrender_set.end()) {
	surrender_flag = '*';
      }

      if (cluster_size >= threshold) {
	outfile << ">" << read.name << "\t" << partition_id
		<< surrender_flag << "\t" << cluster_size << "\n" 
		<< seq << "\n";
	partitions.insert(partition_id);
      }
	       
      // reset the sequence info, increment read number
      total_reads++;

      // run callback, if specified
      if (total_reads % CALLBACK_PERIOD == 0 && callback) {
	try {
	  callback("do_truncated_partition/output", callback_data,
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

  return partitions.size();
}

bool Hashtable::_do_continue(const HashIntoType kmer,
			     const SeenSet& keeper)
{
  const BoundedCounterType val = _counts[kmer % _tablesize];

  if (val == 0) {
    return false;
  }

  // have we already seen me? don't count; exit.
  SeenSet::iterator i = keeper.find(kmer);
  if (i != keeper.end()) {
    return false;
  }

  return true;
}

bool Hashtable::_is_tagged_kmer(const HashIntoType kmer_f,
				const HashIntoType kmer_r,
				const PartitionMap& partition_map,
				HashIntoType& tagged_kmer)
{
  bool found = false;
  PartitionMap::const_iterator fi = partition_map.find(kmer_f);
  if (fi != partition_map.end()) {
    tagged_kmer = kmer_f;
    found = true;
  }

  if (!found) {
    fi = partition_map.find(kmer_r);
    if (fi != partition_map.end()) {
      tagged_kmer = kmer_r;
      found = true;
    }
  }

  return found;
}
		     

// used by do_truncated_partition

void Hashtable::partition_find_all_tags(const HashIntoType kmer_f,
					const HashIntoType kmer_r,
					SeenSet& keeper,
					SeenSet& tagged_kmers,
					const PartitionMap& partition_map,
					bool first,
					unsigned int depth,
					bool& surrender,
					unsigned int& total)
{
  if (depth == 0 || total > PARTITION_MAX_TAG_EXAMINED) {
    surrender = true;
    return;
  }
  depth -= 1;
  total++;

  HashIntoType kmer = uniqify_rc(kmer_f, kmer_r);
  if (!_do_continue(kmer, keeper)) {
    return;
  }

  // keep track of seen kmers
  keeper.insert(kmer);

  // Is this a kmer-to-tag, and have we tagged it already? then exit.
  HashIntoType tagged_kmer;
  if (!first && _is_tagged_kmer(kmer_f, kmer_r, partition_map, tagged_kmer)) {
    tagged_kmers.insert(tagged_kmer);
    return;
  }

  // NEXT.

  HashIntoType f, r;
  const unsigned int rc_left_shift = _ksize*2 - 2;

  f = ((kmer_f << 2) & bitmask) | twobit_repr('A');
  r = kmer_r >> 2 | (twobit_comp('A') << rc_left_shift);
  partition_find_all_tags(f, r, keeper, tagged_kmers, partition_map, false, depth, surrender, total);

  f = ((kmer_f << 2) & bitmask) | twobit_repr('C');
  r = kmer_r >> 2 | (twobit_comp('C') << rc_left_shift);
  partition_find_all_tags(f, r, keeper, tagged_kmers, partition_map, false, depth, surrender, total);

  f = ((kmer_f << 2) & bitmask) | twobit_repr('G');
  r = kmer_r >> 2 | (twobit_comp('G') << rc_left_shift);
  partition_find_all_tags(f, r, keeper, tagged_kmers, partition_map, false, depth, surrender, total);

  f = ((kmer_f << 2) & bitmask) | twobit_repr('T');
  r = kmer_r >> 2 | (twobit_comp('T') << rc_left_shift);
  partition_find_all_tags(f, r, keeper, tagged_kmers, partition_map, false, depth, surrender, total);

  // PREVIOUS.

  r = ((kmer_r << 2) & bitmask) | twobit_comp('A');
  f = kmer_f >> 2 | (twobit_repr('A') << rc_left_shift);
  partition_find_all_tags(f, r, keeper, tagged_kmers, partition_map, false, depth, surrender, total);

  r = ((kmer_r << 2) & bitmask) | twobit_comp('C');
  f = kmer_f >> 2 | (twobit_repr('C') << rc_left_shift);
  partition_find_all_tags(f, r, keeper, tagged_kmers, partition_map, false, depth, surrender, total);

  r = ((kmer_r << 2) & bitmask) | twobit_comp('G');
  f = kmer_f >> 2 | (twobit_repr('G') << rc_left_shift);
  partition_find_all_tags(f, r, keeper, tagged_kmers, partition_map, false, depth, surrender, total);

  r = ((kmer_r << 2) & bitmask) | twobit_comp('T');
  f = kmer_f >> 2 | (twobit_repr('T') << rc_left_shift);
  partition_find_all_tags(f, r, keeper, tagged_kmers, partition_map, false, depth, surrender, total);
}
