#include <iostream>
#include <list>

#include "khmer.hh"
#include "hashtable.hh"

#define CALLBACK_PERIOD 10000
#define MAX_CLUSTER_EXPLORE 250

using namespace khmer;
using namespace std;

MinMaxTable * Hashtable::fasta_file_to_minmax(const std::string &inputfile,
					      unsigned int total_reads,
					      ReadMaskTable * readmask,
					      CallbackFn callback,
					      void * callback_data)
{
   string line;
   ifstream infile(inputfile.c_str());
   string name = "";
   string seq = "";
   unsigned int read_num = 0;

   MinMaxTable * mmt = new MinMaxTable(total_reads);

   if (infile.is_open()) {
     while(!infile.eof()) {
       getline(infile, line);

       if (line[0] == '>' || infile.eof()) {
	 if (seq != "") {

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
	   name.clear();
	   seq.clear();

	   read_num += 1;

	   // run callback, if specified
	   if (read_num % CALLBACK_PERIOD == 0 && callback) {
	     try {
	       callback("fasta_file_to_minmax", callback_data, read_num, 0);
	     } catch (...) {
	       infile.close();
	       delete mmt;
	       throw;
	     }
	   }
	 }
	 if (line[0] == '>') {
	   name = line.substr(1, line.length() - 1);
	 }
       }
       else {
	 seq += line;
       }
     }
   }
  
   infile.close();

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
   string line;
   ifstream infile(inputfile.c_str());
   int isRead = 0;
   string name;
   string seq;
   unsigned int read_num = 0;
   unsigned int n_kept = 0;
   ReadMaskTable * readmask = new ReadMaskTable(total_reads);

   if (old_readmask) {
     readmask->merge(*old_readmask);
   }

   if (infile.is_open()) {
     while(!infile.eof()) {
       getline(infile, line);
       if (line.length() == 0) {
	 break;
       }

       if (isRead) {
	 seq = line;
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
	   name.clear();
	   seq.clear();
	 }

	 read_num += 1;

	 // run callback, if specified
	 if (read_num % CALLBACK_PERIOD == 0 && callback) {
	   try {
	     callback("filter_fasta_file_run", callback_data, read_num,n_kept);
	   } catch (...) {
	     infile.close();
	     throw;
	   }
	 }
       }
       else {
	 name = line.substr(1, line.length()-1);
       }

       isRead = isRead ? 0 : 1;
     }
   }
  
   infile.close();

   return readmask;
}

///
/// output_fasta_kmer_pos_freq: outputs the kmer frequencies for each read
///

void Hashtable::output_fasta_kmer_pos_freq(const std::string &inputfile,
                                           const std::string &outputfile)
{
  string line;
  ifstream infile(inputfile.c_str());
  ofstream outfile;
  outfile.open(outputfile.c_str());
  int isRead = 0;
  string name;
  string seq;

  if (infile.is_open()) {
    while(!infile.eof()) {
      getline(infile, line);
      if (line.length() == 0) {
        break;
      }

      if (isRead) {
        seq = line;

        int numPos = seq.length() - Hashtable::_ksize + 1;

        for (int i = 0; i < numPos; i++)  {
          string kmer = seq.substr(i, Hashtable::_ksize);
          outfile << (int)Hashtable::get_count(kmer.c_str()) << " ";
        }
        outfile << endl;
      }

      isRead = isRead? 0 : 1;
    }
  }

  infile.close();
  outfile.close();
}


unsigned int khmer::output_filtered_fasta_file(const std::string &inputfile,
					       const std::string &outputfile,
					       ReadMaskTable * readmask,
					       CallbackFn callback,
					       void * callback_data)
{
   string line;
   ifstream infile(inputfile.c_str());
   ofstream outfile;
   outfile.open(outputfile.c_str());
   int isRead = 0;
   string name;
   string seq;
   unsigned int n_kept = 0;
   unsigned int read_num = 0;

   if (infile.is_open()) {
     while(!infile.eof()) {
       getline(infile, line);
       if (line.length() == 0) {
	 break;
       }

       if (isRead) {
	 seq = line;

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
	     infile.close();
	     outfile.close();
	     throw;
	   }
	 }

       } else {
	 name = line.substr(1, line.length()-1);
       }

       isRead = isRead? 0 : 1;
     }
   }
  
   infile.close();
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

   for (i = 0; i < read.length(); i++)  {
     if (!is_valid_dna(read[i])) {
         is_valid = false;
         break;
      }
   }

   if (is_valid) {
      return consume_string(read, lower_bound, upper_bound);
   }
   else  {
      return 0;
   }
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

   string line;
   ifstream infile(filename.c_str());

   if (!infile.is_open())  {
     return;
   }

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

   while(1)  {
     getline(infile, line);

     if (line[0] == '>' || infile.eof())  {

       // do we have a sequence to process?
       if (currSeq != "")  {

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
	 currSeq = "";
	 total_reads++;

	 // run callback, if specified
	 if (total_reads % CALLBACK_PERIOD == 0 && callback) {
	   try {
	     callback("consume_fasta", callback_data, total_reads, n_consumed);
	   } catch (...) {
	     infile.close();
	     throw;
	   }
	 }
       }

       // new sequence => new sequence name
       if (line[0] == '>') {
	 currName = line.substr(1, line.length()-1);
       }
     }
     else  {			// additional line for sequence
       currSeq += line;
     }
     
     // @ end of file? break out.
     if (infile.eof()) {
       break;
     }
   }

   infile.close();

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

  HashIntoType mask = 0;
  for (unsigned int i = 0; i < _ksize; i++) {
    mask = mask << 2;
    mask |= 3;
  }

  HashIntoType h = 0, r = 0;
  bool bounded = true;

  if (lower_bound == upper_bound && upper_bound == 0) {
    bounded = false;
  }
  
  HashIntoType bin = _hash(sp, _ksize, &h, &r);

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
    h &= mask;

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

  return n_consumed;
}


BoundedCounterType Hashtable::get_min_count(const std::string &s,
					    HashIntoType lower_bound,
					    HashIntoType upper_bound)
{
  const unsigned int length = s.length();
  const char * sp = s.c_str();
  BoundedCounterType min_count = MAX_COUNT, count;

  HashIntoType mask = 0;
  for (unsigned int i = 0; i < (unsigned int) _ksize; i++) {
    mask = mask << 2;
    mask |= 3;
  }

  HashIntoType h = 0, r = 0;
  bool bounded = true;

  if (lower_bound == upper_bound && upper_bound == 0) {
    bounded = false;
  }

  HashIntoType bin;
  
  bin = _hash(sp, _ksize, &h, &r);
  if (!bounded || (bin >= lower_bound && bin < upper_bound)) {
    min_count = this->get_count(bin);
  }

  for (unsigned int i = _ksize; i < length; i++) {
    // left-shift the previous hash over
    h = h << 2;

    // 'or' in the current nt
    h |= twobit_repr(sp[i]);

    // mask off the 2 bits we shifted over.
    h &= mask;

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

  HashIntoType mask = 0;
  for (unsigned int i = 0; i < (unsigned int) _ksize; i++) {
    mask = mask << 2;
    mask |= 3;
  }

  HashIntoType h = 0, r = 0;
  bool bounded = true;

  if (lower_bound == upper_bound && upper_bound == 0) {
    bounded = false;
  }

  HashIntoType bin = _hash(sp, _ksize, &h, &r);
  if (!bounded || (bin >= lower_bound && bin < upper_bound)) {
    max_count = this->get_count(bin);
  }

  for (unsigned int i = _ksize; i < length; i++) {
    // left-shift the previous hash over
    h = h << 2;

    // 'or' in the current nt
    h |= twobit_repr(sp[i]);

    // mask off the 2 bits we shifted over.
    h &= mask;

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

  string line;
  ifstream infile(inputfile.c_str());
  int isRead = 0;
  string name;
  string seq;
  unsigned int read_num = 0;

   if (infile.is_open()) {
     while(!infile.eof()) {
       getline(infile, line);
       if (line.length() == 0) {
	 break;
       }

       if (isRead) {
	 bool valid_read = true;
	 seq = line;
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
	     infile.close();
	     throw;
	   }
	 }
       }
       else {
	 name = line.substr(1, line.length()-1);
       }

       isRead = isRead ? 0 : 1;
     }
   }
  
   infile.close();

  return counts;
}

void Hashtable::fasta_dump_kmers_by_abundance(const std::string &inputfile,
					      ReadMaskTable * readmask,
					      BoundedCounterType limit_by_count,
					      CallbackFn callback,
					      void * callback_data)
{
  string line;
  ifstream infile(inputfile.c_str());
  int isRead = 0;
  string name;
  string seq;
  unsigned int read_num = 0;

   if (infile.is_open()) {
     while(!infile.eof()) {
       getline(infile, line);
       if (line.length() == 0) {
	 break;
       }

       if (isRead) {
	 bool valid_read = true;
	 seq = line;
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
	     infile.close();
	     throw;
	   }
	 }
       }
       else {
	 name = line.substr(1, line.length()-1);
       }

       isRead = isRead ? 0 : 1;
     }
   }
  
   infile.close();
}

//////////////////////////////////////////////////////////////////////
// graph stuff

void Hashtable::mark_connected_graph(const char * kmer)
{
  const unsigned char seen = 1 << 7;

  HashIntoType bin = _hash(kmer, _ksize); // % _tablesize;
  const BoundedCounterType val = _counts[bin];

  if (empty(val) || marked(val)) {
    return;
  }
  _counts[bin] |= seen;

  // std::cout << kmer << std::endl;

  char new_kmer[_ksize + 1];
  new_kmer[_ksize] = 0;		// NULL terminate
  strncpy(new_kmer, kmer + 1, _ksize - 1);

  new_kmer[_ksize - 1] = 'A';
  mark_connected_graph(new_kmer);
  new_kmer[_ksize - 1] = 'C';
  mark_connected_graph(new_kmer);
  new_kmer[_ksize - 1] = 'G';
  mark_connected_graph(new_kmer);
  new_kmer[_ksize - 1] = 'T';
  mark_connected_graph(new_kmer);

  strncpy(new_kmer + 1, kmer, _ksize - 1);

  new_kmer[0] = 'A';
  mark_connected_graph(new_kmer);
  new_kmer[0] = 'C';
  mark_connected_graph(new_kmer);
  new_kmer[0] = 'G';
  mark_connected_graph(new_kmer);
  new_kmer[0] = 'T';
  mark_connected_graph(new_kmer);
}

unsigned long long Hashtable::zero_connected_graph(const char * kmer)
{
  HashIntoType bin = _hash(kmer, _ksize); // % _tablesize;
  const BoundedCounterType val = _counts[bin];

  if (empty(val)) {
    return 0;
  }
  _counts[bin] = 0;
  unsigned long long removed = 1;

  // std::cout << kmer << std::endl;

  char new_kmer[_ksize + 1];
  new_kmer[_ksize] = 0;		// NULL terminate
  strncpy(new_kmer, kmer + 1, _ksize - 1);

  new_kmer[_ksize - 1] = 'A';
  removed += zero_connected_graph(new_kmer);
  new_kmer[_ksize - 1] = 'C';
  removed += zero_connected_graph(new_kmer);
  new_kmer[_ksize - 1] = 'G';
  removed += zero_connected_graph(new_kmer);
  new_kmer[_ksize - 1] = 'T';
  removed += zero_connected_graph(new_kmer);

  strncpy(new_kmer + 1, kmer, _ksize - 1);

  new_kmer[0] = 'A';
  removed += zero_connected_graph(new_kmer);
  new_kmer[0] = 'C';
  removed += zero_connected_graph(new_kmer);
  new_kmer[0] = 'G';
  removed += zero_connected_graph(new_kmer);
  new_kmer[0] = 'T';
  removed += zero_connected_graph(new_kmer);

  return removed;
}

void Hashtable::clear_marks_for_connected_graph(const char * kmer)
{
  const unsigned char seen = 1 << 7;

  HashIntoType bin = _hash(kmer, _ksize); // % _tablesize;
  const BoundedCounterType val = _counts[bin];

  if (empty(val) || !(_counts[bin] & seen)) {
    return;
  }
  _counts[bin] &= 127;

  // std::cout << kmer << std::endl;

  char new_kmer[_ksize + 1];
  new_kmer[_ksize] = 0;		// NULL terminate
  strncpy(new_kmer, kmer + 1, _ksize - 1);

  new_kmer[_ksize - 1] = 'A';
  clear_marks_for_connected_graph(new_kmer);
  new_kmer[_ksize - 1] = 'C';
  clear_marks_for_connected_graph(new_kmer);
  new_kmer[_ksize - 1] = 'G';
  clear_marks_for_connected_graph(new_kmer);
  new_kmer[_ksize - 1] = 'T';
  clear_marks_for_connected_graph(new_kmer);

  strncpy(new_kmer + 1, kmer, _ksize - 1);

  new_kmer[0] = 'A';
  clear_marks_for_connected_graph(new_kmer);
  new_kmer[0] = 'C';
  clear_marks_for_connected_graph(new_kmer);
  new_kmer[0] = 'G';
  clear_marks_for_connected_graph(new_kmer);
  new_kmer[0] = 'T';
  clear_marks_for_connected_graph(new_kmer);
}


void Hashtable::calc_connected_graph_size(const char * kmer,
					  unsigned long long& count,
					  unsigned long long threshold)
{
  const unsigned char seen = 1 << 7;

  HashIntoType bin = _hash(kmer, _ksize); // % _tablesize;
  const BoundedCounterType val = _counts[bin];

  if (empty(val) || marked(val)) {
    return;
  }
  _counts[bin] |= seen;
  count += 1;

  if (threshold && count >= threshold) {
    return;
  }

  char new_kmer[_ksize + 1];
  new_kmer[_ksize] = 0;		// NULL terminate
  strncpy(new_kmer, kmer + 1, _ksize - 1);

  new_kmer[_ksize - 1] = 'A';
  calc_connected_graph_size(new_kmer, count);
  new_kmer[_ksize - 1] = 'C';
  calc_connected_graph_size(new_kmer, count);
  new_kmer[_ksize - 1] = 'G';
  calc_connected_graph_size(new_kmer, count);
  new_kmer[_ksize - 1] = 'T';
  calc_connected_graph_size(new_kmer, count);

  strncpy(new_kmer + 1, kmer, _ksize - 1);

  new_kmer[0] = 'A';
  calc_connected_graph_size(new_kmer, count);
  new_kmer[0] = 'C';
  calc_connected_graph_size(new_kmer, count);
  new_kmer[0] = 'G';
  calc_connected_graph_size(new_kmer, count);
  new_kmer[0] = 'T';
  calc_connected_graph_size(new_kmer, count);
}

void Hashtable::calc_connected_graph_size2(const char * kmer,
					   unsigned long long& count,
					   unsigned long long threshold,
					   const HashIntoType watermark)
{
  const unsigned char seen = 1 << 7;

  HashIntoType bin = _hash(kmer, _ksize); // % _tablesize;
  const BoundedCounterType val = _counts[bin];

  // nothing here, go home.
  if (empty(val)) {
    return;
  }

  // below the watermark, already been here & kept it, go away.
  if (bin < watermark) {
    // std::cout << "watermark\n";
    count = threshold + 1;
    return;
  }

  // have we already seen me? don't count; exit.
  if (val & seen) {
    // std::cout << "seen.\n";
    return;
  }

  // mark as seen, to prevent cycles... and increment count.
  _counts[bin] |= seen;
  count += 1;

  // now that we've counted you, are we above the threshold? if so, exit.
  if (count < threshold) {
    char new_kmer[_ksize + 1];
    new_kmer[_ksize] = 0;		// NULL terminate

    strncpy(new_kmer + 1, kmer, _ksize - 1);

    new_kmer[0] = 'A';
    calc_connected_graph_size2(new_kmer, count, threshold, watermark);
    new_kmer[0] = 'C';
    calc_connected_graph_size2(new_kmer, count, threshold, watermark);
    new_kmer[0] = 'G';
    calc_connected_graph_size2(new_kmer, count, threshold, watermark);
    new_kmer[0] = 'T';
    calc_connected_graph_size2(new_kmer, count, threshold, watermark);

    strncpy(new_kmer, kmer + 1, _ksize - 1);

    new_kmer[_ksize - 1] = 'T';
    calc_connected_graph_size2(new_kmer, count, threshold, watermark);
    new_kmer[_ksize - 1] = 'A';
    calc_connected_graph_size2(new_kmer, count, threshold, watermark);
    new_kmer[_ksize - 1] = 'C';
    calc_connected_graph_size2(new_kmer, count, threshold, watermark);
    new_kmer[_ksize - 1] = 'G';
    calc_connected_graph_size2(new_kmer, count, threshold, watermark);
  }
}

bool Hashtable::is_graph_size_larger(const char * kmer,
				     const unsigned long long threshold)
{
  return false;
}

void Hashtable::empty_bins(bool empty_marked)
{
  for (HashIntoType i = 0; i < _tablesize; i++) {
    if (!_counts[i]) {		// no counts, not marked
      continue;
    }

    if (marked(_counts[i])) {
      if (empty_marked) {
	_counts[i] = 0;
      }
    }
    else {
      if (!empty_marked) {
	_counts[i] = 0;
      }
    }
  }
}

void Hashtable::trim_graphs(unsigned int min_size)
{
  const unsigned char seen = 1 << 7;
  unsigned int max_depth = MAX_CLUSTER_EXPLORE;
  unsigned long long queried = 0;
  unsigned long long removed = 0;

  if (min_size > max_depth) {
    max_depth = min_size + 2;
  }

  for (HashIntoType i = 0; i < _tablesize; i++) {
    if (_counts[i] && !(_counts[i] & seen)) {
      std::cout << "at: " << i / 1000000 << "m of " << _tablesize / 1000000 << "m queried; " << removed/1e9 << "b removed\n";
      std::string kmer = _revhash(i, _ksize);

      // ASSUME: we are at the lowest-valued entry in the hashtable.
      HashIntoType watermark = i;
      unsigned long long clustersize = 0;
      calc_connected_graph_size2(kmer.c_str(), clustersize, max_depth,
				 watermark);

      queried += 1;

      if (clustersize >= min_size) {
	if (clustersize >= max_depth) {
	  clear_marks_for_connected_graph(kmer.c_str());
	}
      } else {
	removed += zero_connected_graph(kmer.c_str());
      }
    }
  }

  clear_marks();
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
      calc_connected_graph_size(kmer.c_str(), size);
      if (size) {
	if (size > 5000) { std::cout << "GRAPH SIZE: " << size << "\n"; }
	if (size >= max_size) {
	  size = max_size;
	}
	p[size] += 1;
      }
    }
  }

  clear_marks();

  return p;
}
