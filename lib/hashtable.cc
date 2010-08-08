#include <iostream>
#include <list>

#include "khmer.hh"
#include "hashtable.hh"

#define CALLBACK_PERIOD 10000
#define ID3_DEPTH 500
#define ID4_DEPTH 100

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
   string line;
   ifstream infile(readsfile.c_str());
   int isRead = 0;
   string name;
   string seq;

   unsigned int read_num = 0;
   const unsigned int tablesize = minmax.get_tablesize();

   ReadMaskTable * readmask = new ReadMaskTable(tablesize);

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
 
       isRead = isRead? 0 : 1;
     }
   }

   infile.close();

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

        int numPos = seq.length() - _ksize + 1;

        for (int i = 0; i < numPos; i++)  {
          string kmer = seq.substr(i, _ksize);
          outfile << (int)get_count(kmer.c_str()) << " ";
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
    }

    // new sequence => new sequence name
    if (line[0] == '>') {
      currName = line.substr(1, line.length()-1);
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
  unsigned int total_reads = 0;
  unsigned int reads_kept = 0;

  string line;
  ifstream infile(infilename.c_str());
  ofstream outfile(outfilename.c_str());

  if (!infile.is_open())  {
    return;
  }

   string currName = "";
   string currSeq = "";

   //
   // iterate through the FASTA file & consume the reads.
   //

   while(1)  {
     getline(infile, line);

     if (line[0] == '>' || infile.eof())  {

       // do we have a sequence to process?
       if (currSeq != "")  {

	 // yep! process.

	 bool is_valid;

	 check_and_process_read(currSeq, is_valid);

	 if (is_valid) {
	   std::string first_kmer = currSeq.substr(0, _ksize);
	   unsigned long long clustersize = 0;
	   SeenSet keeper;
	   calc_connected_graph_size(first_kmer.c_str(), clustersize, keeper,
				     min_size);

	   if (clustersize >= min_size) {
	     outfile << ">" << currName << endl;
	     outfile << currSeq << endl;
	     reads_kept++;
	   }
	 }
	       
	 // reset the sequence info, increment read number
	 currSeq = "";
	 total_reads++;

	 // run callback, if specified
	 if (total_reads % CALLBACK_PERIOD == 0 && callback) {
	   try {
	     callback("trim_graphs", callback_data, total_reads, reads_kept);
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
	// if (size > 5000) { std::cout << "GRAPH SIZE: " << size << "\n"; }
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

// do_partition: simple partitioning, done once per cluster.
//   1) load in all the sequences, tagging the first kmer of each sequence.
//   2) then, for each tag, explore entre cluster & set tag
//         using partition_set_id.
//
// slow for big clusters, because it has to find all tagged k-mers in each
// cluster.  no provision for giving up and retagging.

unsigned int Hashtable::do_partition(const std::string infilename,
				      CallbackFn callback,
				      void * callback_data)
{
  PartitionMap partition_map;

  unsigned int total_reads = 0;
  unsigned int reads_kept = 0;

  string line;
  ifstream infile(infilename.c_str());

  if (!infile.is_open())  {
    return 0;
  }

  string currName = "";
  string currSeq = "";

  while(1)  {
     getline(infile, line);

     if (line[0] == '>' || infile.eof())  {

       // do we have a sequence to process?
       if (currSeq != "")  {

	 // yep! process.

	 bool is_valid;

	 check_and_process_read(currSeq, is_valid);

	 if (is_valid) {
	   std::string first_kmer = currSeq.substr(0, _ksize);
	   HashIntoType kmer_f = _hash_forward(first_kmer.c_str(), _ksize);

	   partition_map[kmer_f] = 0;
	 }
	       
	 // reset the sequence info, increment read number
	 currSeq = "";
	 total_reads++;

	 // run callback, if specified
	 if (total_reads % CALLBACK_PERIOD == 0 && callback) {
	   try {
	     callback("AAAAA", callback_data, total_reads, reads_kept);
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
   
   unsigned int next_partition_id = 1;
   for(PartitionMap::iterator i=partition_map.begin();
       i != partition_map.end(); ++i) {
     HashIntoType kmer_f = (*i).first;
     unsigned int partition_id = (*i).second;

     if (partition_id == 0) {
       partition_id = next_partition_id;
       next_partition_id++;

       std::string kmer_s = _revhash(kmer_f, _ksize);
       HashIntoType kmer_r;

       _hash(kmer_s.c_str(), _ksize, kmer_f, kmer_r);

#if PARTITION2_DEBUG
       std::cout << "traversing partition " << partition_id << "\n";
#endif // PARTITION2_DEBUG

       SeenSet keeper;
       partition_set_id(kmer_f, kmer_r, keeper, partition_id, partition_map);
#if PARTITION2_DEBUG
       std::cout << "graph size: " << keeper.size() << "\n";
#endif // PARTITION2_DEBUG
     }
   }

   return next_partition_id - 1;
}

// used by do_partition3 and do_partition4.

void Hashtable::partition_find_id3(const HashIntoType kmer_f,
				   const HashIntoType kmer_r,
				   SeenSet& keeper,
				   SeenSet& tagged_kmers,
				   PartitionMap& partition_map,
				   ReversePartitionMap& rev_pmap,
				   bool& done,
				   bool first,
				   unsigned int depth)
{
  if (done || depth == 0) return;

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

    if (!first) {
      // Is this a kmer-to-tag, and have we tagged it already?
      bool found = false;
      PartitionMap::iterator fi = partition_map.find(kmer_f);
      if (fi != partition_map.end()) {
	tagged_kmers.insert(kmer_f);
	found = true;
      }

      fi = partition_map.find(kmer_r);
      if (!found && fi != partition_map.end()) {
	tagged_kmers.insert(kmer_r);
	found = true;
      }

      if (found) {
	done = true;
	return;
      }
    }
  }

  // NEXT.

  HashIntoType f, r;
  const unsigned int rc_left_shift = _ksize*2 - 2;

  f = ((kmer_f << 2) & bitmask) | twobit_repr('A');
  r = kmer_r >> 2 | (twobit_comp('A') << rc_left_shift);
  partition_find_id3(f, r, keeper, tagged_kmers, partition_map, rev_pmap, done, false, depth - 1);

  f = ((kmer_f << 2) & bitmask) | twobit_repr('C');
  r = kmer_r >> 2 | (twobit_comp('C') << rc_left_shift);
  partition_find_id3(f, r, keeper, tagged_kmers, partition_map, rev_pmap, done, false, depth - 1);

  f = ((kmer_f << 2) & bitmask) | twobit_repr('G');
  r = kmer_r >> 2 | (twobit_comp('G') << rc_left_shift);
  partition_find_id3(f, r, keeper, tagged_kmers, partition_map, rev_pmap, done, false, depth - 1);

  f = ((kmer_f << 2) & bitmask) | twobit_repr('T');
  r = kmer_r >> 2 | (twobit_comp('T') << rc_left_shift);
  partition_find_id3(f, r, keeper, tagged_kmers, partition_map, rev_pmap, done, false, depth - 1);

  // PREVIOUS.

  r = ((kmer_r << 2) & bitmask) | twobit_comp('A');
  f = kmer_f >> 2 | (twobit_repr('A') << rc_left_shift);
  partition_find_id3(f, r, keeper, tagged_kmers, partition_map, rev_pmap, done, false, depth - 1);

  r = ((kmer_r << 2) & bitmask) | twobit_comp('C');
  f = kmer_f >> 2 | (twobit_repr('C') << rc_left_shift);
  partition_find_id3(f, r, keeper, tagged_kmers, partition_map, rev_pmap, done, false, depth - 1);

  r = ((kmer_r << 2) & bitmask) | twobit_comp('G');
  f = kmer_f >> 2 | (twobit_repr('G') << rc_left_shift);
  partition_find_id3(f, r, keeper, tagged_kmers, partition_map, rev_pmap, done, false, depth - 1);

  r = ((kmer_r << 2) & bitmask) | twobit_comp('T');
  f = kmer_f >> 2 | (twobit_repr('T') << rc_left_shift);
  partition_find_id3(f, r, keeper, tagged_kmers, partition_map, rev_pmap, done, false, depth - 1);
}

// do_partition4: less truncated progressive partitioning.
//   1) load all sequences, tagging first kmer of each
//   2) do a truncated DFS search for the *first* tagged kmer; assign cluster
//         (partition_find_id3, to ID3_DEPTH)
//   3) after loading all sequences, do a truncated DFS search for *all* tagged
//         kmers (partition_find_i4, to ID4_DEPTH), reassigning now-connected
//         clusters.
//
// CTB note: for unlimited ID4_DEPTH, yields perfect clustering.

unsigned int Hashtable::do_partition4(const std::string infilename,
				      CallbackFn callback,
				      void * callback_data)
{
  PartitionMap partition_map;
  ReversePartitionMap rev_pmap;

  unsigned int total_reads = 0;
  unsigned int reads_kept = 0;
  unsigned int next_partition_id = 1;
  SeenSet keeper;
  SeenSet tagged_kmers;
	   

  string line;
  ifstream infile(infilename.c_str());

  if (!infile.is_open())  {
    return 0;
  }

  string currName = "";
  string currSeq = "";

  while(1)  {
     getline(infile, line);

     if (line[0] == '>' || infile.eof())  {

       // do we have a sequence to process?
       if (currSeq != "")  {

	 // yep! process.

	 bool is_valid;

	 check_and_process_read(currSeq, is_valid);

	 if (is_valid) {
	   std::string first_kmer = currSeq.substr(0, _ksize);

	   HashIntoType kmer_f, kmer_r;
	   _hash(first_kmer.c_str(), _ksize, kmer_f, kmer_r);
	   bool done = false;

	   keeper.empty();
	   tagged_kmers.empty();
	   partition_find_id3(kmer_f, kmer_r, keeper, tagged_kmers,
			      partition_map, rev_pmap, done);

	   unsigned int partition_id;

	   if (tagged_kmers.size() == 0) {
	     assert(!done);
	     partition_id = next_partition_id;
	     next_partition_id++;

	     partition_map[kmer_f] = partition_id;

	     SeenSet * x = new SeenSet();
	     x->insert(kmer_f);
	     rev_pmap[partition_id] = x;
	   } else {
	     assert(tagged_kmers.size() == 1);
	     SeenSet::iterator it = tagged_kmers.begin();
	     partition_id = partition_map[*it]; // get graph ID of first tagged kmer

	     partition_map[kmer_f] = partition_id;
	     rev_pmap[partition_id]->insert(kmer_f);
	   }
	 }
	       
	 // reset the sequence info, increment read number
	 currSeq = "";
	 total_reads++;

	 // run callback, if specified
	 if (total_reads % CALLBACK_PERIOD == 0 && callback) {
	   try {
	     callback("AAAAA", callback_data, total_reads, reads_kept);
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

   ///

   PartitionMap::iterator pi;
   unsigned int n_done = 0;
   unsigned int n = 0;
   for (pi = partition_map.begin(); pi != partition_map.end(); ++pi) {
     n++;
     if (n % 10000 == 0) {
       std::cout << "...fixing " << n << " of " << partition_map.size() << "\n";
     }
     std::string kmer = _revhash((*pi).first, _ksize);

     HashIntoType kmer_f, kmer_r;
     _hash(kmer.c_str(), _ksize, kmer_f, kmer_r);
     SeenSet keeper;
     SeenSet tagged_kmers;
     bool done = false;

     partition_find_id4(kmer_f, kmer_r, keeper, tagged_kmers, partition_map,
			rev_pmap, done, true, ID4_DEPTH);

     if (tagged_kmers.size() >= 1) {
       unsigned int this_pid = partition_map[kmer_f];

       set<unsigned int> other_partition_ids;
       SeenSet::iterator it = tagged_kmers.begin();
       for (; it != tagged_kmers.end(); ++it) {
	 unsigned int id = partition_map[*it];

	 if (id != this_pid) {
#ifdef DEBUG_PRINT_4
	   std::cout << "inserting " << id << "\n";
#endif // DEBUG_PRINT_4
	   other_partition_ids.insert(id);
	 }
       }

       if (other_partition_ids.size()) {
	 for (set<unsigned int>::iterator si = other_partition_ids.begin();
	      si != other_partition_ids.end(); si++) {
	   SeenSet * x = rev_pmap[*si];
#ifdef DEBUG_PRINT_4
	   std::cout << "looking at opid: " << *si << " - " << x->size() << "\n";	 
#endif // DEBUG_PRINT_4
	   for (SeenSet::iterator pi = x->begin(); pi != x->end(); ++pi){
#ifdef DEBUG_PRINT_4
	     std::cout << "reassigning " << *pi << " to " << this_pid << "--" << n_done << "\n";
#endif // DEBUG_PRINT_4

	     partition_map[*pi] = this_pid;
	     rev_pmap[this_pid]->insert(*pi);
	   }

	   rev_pmap.erase(*si);
	   delete x;
	 }
       }
     }
     n_done++;
   }

   ///

#if 0
   ReversePartitionMap::iterator ri;
   for (ri = rev_pmap.begin(); ri != rev_pmap.end(); ri++) {
     SeenSet * x = (*ri).second;
     if (x->size() >= 2) {
#ifdef DEBUG_PRINT_4
       std::cout << "partition: " << (*ri).first << "\n";
#endif // DEBUG_PRINT_4
       SeenSet::iterator j;
       for (j = x->begin(); j != x->end(); j++) {
	 HashIntoType kmer_f = *j;
#ifdef DEBUG_PRINT_4
	 std::cout << x->size() << " <- " << _revhash(kmer_f, _ksize) << "\n";
#endif // DEBUG_PRINT_4
       }
     }
   }
#endif

   SeenSet unique_partitions;
   unsigned int i = 0;
   for (PartitionMap::iterator pi = partition_map.begin();
	pi != partition_map.end(); ++pi, i++) {
     unique_partitions.insert((*pi).second);
     // std::cout << i << " - " << (*pi).second << "\n";
   }

   return unique_partitions.size();
}

// used by do_partition4.

void Hashtable::partition_find_id4(const HashIntoType kmer_f,
				   const HashIntoType kmer_r,
				   SeenSet& keeper,
				   SeenSet& tagged_kmers,
				   PartitionMap& partition_map,
				   ReversePartitionMap& rev_pmap,
				   bool& done,
				   bool first,
				   unsigned int depth)
{
  // if (done) return;
  if (depth == 0) return;

  depth -= 1;

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

    if (!first) {
      // Is this a kmer-to-tag, and have we tagged it already?
      bool found = false;
      PartitionMap::iterator fi = partition_map.find(kmer_f);
      if (fi != partition_map.end()) {
	tagged_kmers.insert(kmer_f);
	found = true;
      }

      fi = partition_map.find(kmer_r);
      if (fi != partition_map.end()) {
	tagged_kmers.insert(kmer_r);
	found = true;
      }

      if (found) {
	done = true;
	return;
      }
    }
  }

  // NEXT.

  HashIntoType f, r;
  const unsigned int rc_left_shift = _ksize*2 - 2;

  f = ((kmer_f << 2) & bitmask) | twobit_repr('A');
  r = kmer_r >> 2 | (twobit_comp('A') << rc_left_shift);
  partition_find_id4(f, r, keeper, tagged_kmers, partition_map, rev_pmap, done, false, depth);

  f = ((kmer_f << 2) & bitmask) | twobit_repr('C');
  r = kmer_r >> 2 | (twobit_comp('C') << rc_left_shift);
  partition_find_id4(f, r, keeper, tagged_kmers, partition_map, rev_pmap, done, false, depth);

  f = ((kmer_f << 2) & bitmask) | twobit_repr('G');
  r = kmer_r >> 2 | (twobit_comp('G') << rc_left_shift);
  partition_find_id4(f, r, keeper, tagged_kmers, partition_map, rev_pmap, done, false, depth);

  f = ((kmer_f << 2) & bitmask) | twobit_repr('T');
  r = kmer_r >> 2 | (twobit_comp('T') << rc_left_shift);
  partition_find_id4(f, r, keeper, tagged_kmers, partition_map, rev_pmap, done, false, depth);

  // PREVIOUS.

  r = ((kmer_r << 2) & bitmask) | twobit_comp('A');
  f = kmer_f >> 2 | (twobit_repr('A') << rc_left_shift);
  partition_find_id4(f, r, keeper, tagged_kmers, partition_map, rev_pmap, done, false, depth);

  r = ((kmer_r << 2) & bitmask) | twobit_comp('C');
  f = kmer_f >> 2 | (twobit_repr('C') << rc_left_shift);
  partition_find_id4(f, r, keeper, tagged_kmers, partition_map, rev_pmap, done, false, depth);

  r = ((kmer_r << 2) & bitmask) | twobit_comp('G');
  f = kmer_f >> 2 | (twobit_repr('G') << rc_left_shift);
  partition_find_id4(f, r, keeper, tagged_kmers, partition_map, rev_pmap, done, false, depth);

  r = ((kmer_r << 2) & bitmask) | twobit_comp('T');
  f = kmer_f >> 2 | (twobit_repr('T') << rc_left_shift);
  partition_find_id4(f, r, keeper, tagged_kmers, partition_map, rev_pmap, done, false, depth);
}
