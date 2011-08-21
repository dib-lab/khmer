#include "hashtable.hh"
#include "counting.hh"
#include "hashbits.hh"
#include "parsers.hh"

#include <math.h>
#include <algorithm>

using namespace std;
using namespace khmer;

MinMaxTable * CountingHash::fasta_file_to_minmax(const std::string &inputfile,
					      unsigned long long total_reads,
					      ReadMaskTable * readmask,
					      CallbackFn callback,
					      void * callback_data)
{
   IParser* parser = IParser::get_parser(inputfile.c_str());
   Read read;
   string seq = "";
   unsigned long long read_num = 0;

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

ReadMaskTable * CountingHash::filter_fasta_file_any(MinMaxTable &minmax,
						 BoundedCounterType threshold,
						 ReadMaskTable * old_readmask,
						 CallbackFn callback,
						 void * callback_data)

{
   unsigned long long read_num;
   const unsigned long long tablesize = minmax.get_tablesize();
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

ReadMaskTable * CountingHash::filter_fasta_file_limit_n(const std::string &readsfile,
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
   unsigned long long read_num = 0;
   const unsigned long long tablesize = minmax.get_tablesize();

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

ReadMaskTable * CountingHash::filter_fasta_file_all(MinMaxTable &minmax,
						 BoundedCounterType threshold,
						 ReadMaskTable * old_readmask,
						 CallbackFn callback,
						 void * callback_data)
{
   unsigned long long read_num;
   const unsigned long long tablesize = minmax.get_tablesize();

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

ReadMaskTable * CountingHash::filter_fasta_file_run(const std::string &inputfile,
						 unsigned long long total_reads,
						 BoundedCounterType threshold,
						 unsigned int runlength,
						 ReadMaskTable * old_readmask,
						 CallbackFn callback,
						 void * callback_data)

{
   IParser* parser = IParser::get_parser(inputfile.c_str());
   string seq;
   Read read;
   unsigned long long read_num = 0;
   unsigned long long n_kept = 0;
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

void CountingHash::output_fasta_kmer_pos_freq(const std::string &inputfile,
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


unsigned long long khmer::output_filtered_fasta_file(const std::string &inputfile,
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
   unsigned long long n_kept = 0;
   unsigned long long read_num = 0;


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

BoundedCounterType CountingHash::get_min_count(const std::string &s,
					    HashIntoType lower_bound,
					    HashIntoType upper_bound)
{
  KMerIterator kmers(s.c_str(), _ksize);
  HashIntoType kmer;

  BoundedCounterType min_count = MAX_COUNT, count;

  bool bounded = true;
  if (lower_bound == upper_bound && upper_bound == 0) {
    bounded = false;
  }

  while(!kmers.done()) {
    kmer = kmers.next();

    if (!bounded || (kmer >= lower_bound && kmer < upper_bound)) {
      count = this->get_count(kmer);
    
      if (count < min_count) {
	min_count = count;
      }
    }
  }
  return min_count;
}

BoundedCounterType CountingHash::get_max_count(const std::string &s,
					    HashIntoType lower_bound,
					    HashIntoType upper_bound)
{
  KMerIterator kmers(s.c_str(), _ksize);

  BoundedCounterType max_count = 0, count;

  bool bounded = true;
  if (lower_bound == upper_bound && upper_bound == 0) {
    bounded = false;
  }

  HashIntoType kmer;
  while(!kmers.done()) {
    kmer = kmers.next();

    if (!bounded || (kmer >= lower_bound && kmer < upper_bound)) {
      count = this->get_count(kmer);

      if (count > max_count) {
	max_count = count;
      }
    }
  }
  return max_count;
}

HashIntoType * CountingHash::abundance_distribution(std::string filename,
						    Hashbits * tracking,
			    CallbackFn callback,
			    void * callback_data) const
{
  HashIntoType * dist = new HashIntoType[MAX_BIGCOUNT + 1];
  HashIntoType i;
  
  for (i = 0; i <= MAX_BIGCOUNT; i++) {
    dist[i] = 0;
  }

  Read read;
  IParser* parser = IParser::get_parser(filename.c_str());
  string name;
  string seq;
  unsigned long long read_num = 0;

  // if not, could lead to overflow.
  assert(sizeof(BoundedCounterType) == 2);

  while(!parser->is_complete()) {
    read = parser->get_next_read();
    seq = read.seq;

    if (check_read(seq)) {
      HashIntoType kmer;
      KMerIterator kmers(seq.c_str(), _ksize);

      while(!kmers.done()) {
	kmer = kmers.next();

	if (!tracking->get_count(kmer)) {
	  tracking->count(kmer);

	  BoundedCounterType n = get_count(kmer);
	  dist[n]++;
	}
      }

      name.clear();
      seq.clear();
    }

    read_num += 1;

    // run callback, if specified
    if (read_num % CALLBACK_PERIOD == 0 && callback) {
      try {
        callback("abundance_distribution", callback_data, read_num, 0);
      } catch (...) {
        throw;
      }
    }
  }

  return dist;
}

HashIntoType * CountingHash::fasta_count_kmers_by_position(const std::string &inputfile,
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
   unsigned long long read_num = 0;

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

void CountingHash::fasta_dump_kmers_by_abundance(const std::string &inputfile,
					      ReadMaskTable * readmask,
					      BoundedCounterType limit_by_count,
					      CallbackFn callback,
					      void * callback_data)
{
  Read read;
  IParser* parser = IParser::get_parser(inputfile.c_str());
  string name;
  string seq;
  unsigned long long read_num = 0;

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

void CountingHash::save(std::string outfilename)
{
  assert(_counts[0]);

  unsigned int save_ksize = _ksize;
  unsigned char save_n_tables = _n_tables;
  unsigned long long save_tablesize;

  ofstream outfile(outfilename.c_str(), ios::binary);

  unsigned char version = SAVED_FORMAT_VERSION;
  outfile.write((const char *) &version, 1);

  unsigned char ht_type = SAVED_COUNTING_HT;
  outfile.write((const char *) &ht_type, 1);

  unsigned char use_bigcount = 0;
  if (_use_bigcount) {
    use_bigcount = 1;
  }
  outfile.write((const char *) &use_bigcount, 1);

  outfile.write((const char *) &save_ksize, sizeof(save_ksize));
  outfile.write((const char *) &save_n_tables, sizeof(save_n_tables));

  for (unsigned int i = 0; i < _n_tables; i++) {
    save_tablesize = _tablesizes[i];

    outfile.write((const char *) &save_tablesize, sizeof(save_tablesize));
    outfile.write((const char *) _counts[i], save_tablesize);
  }

  HashIntoType n_counts = _bigcounts.size();
  outfile.write((const char *) &n_counts, sizeof(n_counts));

  if (n_counts) {
    KmerCountMap::const_iterator it = _bigcounts.begin();
    
    for (; it != _bigcounts.end(); it++) {
      outfile.write((const char *) &it->first, sizeof(it->first));
      outfile.write((const char *) &it->second, sizeof(it->second));
    }
  }

  outfile.close();
}

void CountingHash::load(std::string infilename)
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
  unsigned char version, ht_type, use_bigcount;

  ifstream infile(infilename.c_str(), ios::binary);

  infile.read((char *) &version, 1);
  infile.read((char *) &ht_type, 1);
  assert(version == SAVED_FORMAT_VERSION);
  assert(ht_type == SAVED_COUNTING_HT);

  infile.read((char *) &use_bigcount, 1);
  infile.read((char *) &save_ksize, sizeof(save_ksize));
  infile.read((char *) &save_n_tables, sizeof(save_n_tables));

  _ksize = (WordLength) save_ksize;
  _n_tables = (unsigned int) save_n_tables;
  _init_bitstuff();

  _use_bigcount = use_bigcount;

  _counts = new Byte*[_n_tables];
  for (unsigned int i = 0; i < _n_tables; i++) {
    HashIntoType tablesize;

    infile.read((char *) &save_tablesize, sizeof(save_tablesize));

    tablesize = (HashIntoType) save_tablesize;
    _tablesizes.push_back(tablesize);

    _counts[i] = new Byte[tablesize];

    unsigned long long loaded = 0;
    while (loaded != tablesize) {
      infile.read((char *) _counts[i], tablesize - loaded);
      loaded += infile.gcount();	// do I need to do this loop?
    }
  }

  HashIntoType n_counts = 0;
  infile.read((char *) &n_counts, sizeof(n_counts));

  if (n_counts) {
    _bigcounts.clear();

    HashIntoType kmer;
    BoundedCounterType count;

    for (HashIntoType n = 0; n < n_counts; n++) {
      infile.read((char *) &kmer, sizeof(kmer));
      infile.read((char *) &count, sizeof(count));
      _bigcounts[kmer] = count;
    }
  }

  infile.close();
}

// technically, get medioid count... our "median" is always a member of the
// population.

void CountingHash::get_median_count(const std::string &s,
				    BoundedCounterType &median,
				    float &average,
				    float &stddev)
{
  BoundedCounterType count;
  std::vector<BoundedCounterType> counts;
  KMerIterator kmers(s.c_str(), _ksize);

  while(!kmers.done()) {
    HashIntoType kmer = kmers.next();
    count = this->get_count(kmer);
    counts.push_back(count);
  }

  assert(counts.size());

  if (!counts.size()) {
    median = 0;
    average = 0;
    stddev = 0;

    return;
  }

  average = 0;
  for (std::vector<BoundedCounterType>::const_iterator i = counts.begin();
       i != counts.end(); i++) {
    average += *i;
  }
  average /= float(counts.size());

  stddev = 0;
  for (std::vector<BoundedCounterType>::const_iterator i = counts.begin();
       i != counts.end(); i++) {
    stddev += (float(*i) - average) * (float(*i) - average);
  }
  stddev /= float(counts.size());
  stddev = sqrt(stddev);

  sort(counts.begin(), counts.end());
  median = counts[counts.size() / 2]; // rounds down
}


void CountingHash::get_kmer_abund_mean(const std::string &filename,
				       unsigned long long &total,
				       unsigned long long &count,
				       float &mean) const
{
  total = 0;
  count = 0;
  mean = 0.0;

  Read read;
  IParser* parser = IParser::get_parser(filename.c_str());
  string name;
  string seq;
  unsigned long long read_num = 0;

  while(!parser->is_complete()) {
    read = parser->get_next_read();
    seq = read.seq;

    if (check_read(seq)) {
      for (unsigned int i = 0; i < seq.length() - _ksize + 1; i++) {
	string kmer = seq.substr(i, i + _ksize);
	BoundedCounterType n = get_count(kmer.c_str());

	total += n;
	count ++;
      }

      name.clear();
      seq.clear();
    }

    read_num += 1;

#if 0
    // run callback, if specified
    if (read_num % CALLBACK_PERIOD == 0 && callback) {
      try {
        callback("abundance_distribution", callback_data, read_num, 0);
      } catch (...) {
        throw;
      }
    }
#endif // 0
  }

  mean = float(total) / float(count);
}

void CountingHash::get_kmer_abund_abs_deviation(const std::string &filename,
						float mean,
						float &abs_deviation) const
{
  float total = 0.0;
  unsigned long long count = 0;

  Read read;
  IParser* parser = IParser::get_parser(filename.c_str());
  string name;
  string seq;
  unsigned long long read_num = 0;

  while(!parser->is_complete()) {
    read = parser->get_next_read();
    seq = read.seq;

    if (check_read(seq)) {
      for (unsigned int i = 0; i < seq.length() - _ksize + 1; i++) {
	string kmer = seq.substr(i, i + _ksize);
	BoundedCounterType n = get_count(kmer.c_str());

	float diff = mean - (unsigned int)n;
	if (diff < 0) { diff = -diff; }
	total += diff;
	count ++;
      }

      name.clear();
      seq.clear();
    }

    read_num += 1;

#if 0
    // run callback, if specified
    if (read_num % CALLBACK_PERIOD == 0 && callback) {
      try {
        callback("abundance_distribution", callback_data, read_num, 0);
      } catch (...) {
        throw;
      }
    }
#endif // 0
  }

  abs_deviation = total / float(count);
}

unsigned int CountingHash::max_hamming1_count(const std::string kmer_s)
{
  std::string ksub;

  unsigned int max_count = 0;
  for (unsigned int i = 0; i < _ksize; i++) {
    unsigned int the_count;

    ksub = kmer_s;
    ksub[i] = 'A';
    the_count = get_count(ksub.c_str());
    if (the_count > max_count) { max_count = the_count; }

    ksub[i] = 'C';
    the_count = get_count(ksub.c_str());
    if (the_count > max_count) { max_count = the_count; }

    ksub[i] = 'G';
    the_count = get_count(ksub.c_str());
    if (the_count > max_count) { max_count = the_count; }

    ksub[i] = 'T';
    the_count = get_count(ksub.c_str());
    if (the_count > max_count) { max_count = the_count; }
  }

  return max_count;
}

unsigned int CountingHash::trim_on_abundance(std::string seq,
					     BoundedCounterType min_abund)
  const
{
  if (!check_read(seq)) {
    return 0;
  }

  KMerIterator kmers(seq.c_str(), _ksize);

  SeenSet path;

  HashIntoType kmer;

  if (kmers.done()) { return 0; }
  kmer = kmers.next();

  if (kmers.done() || get_count(kmer) < min_abund) {
    return 0;
  }

  unsigned int i = _ksize;
  while (!kmers.done()) {
    kmer = kmers.next();

    if (get_count(kmer) < min_abund) {
      return i;
    }
    i++;
  }

  return seq.length();
}

