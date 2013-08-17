#include "hashtable.hh"
#include "counting.hh"
#include "hashbits.hh"
#include "read_parsers.hh"

#include "zlib/zlib.h"
#include <math.h>
#include <algorithm>

using namespace std;
using namespace khmer;
using namespace khmer:: read_parsers;

MinMaxTable * CountingHash::fasta_file_to_minmax(
    const std::string &inputfile,
    unsigned long long total_reads,
    CallbackFn callback,
    void * callback_data
)
{
    IParser* parser = IParser::get_parser(inputfile.c_str());
    Read read;
    string seq = "";
    unsigned long long read_num = 0;

    MinMaxTable * mmt = new MinMaxTable(total_reads);

    while(!parser->is_complete()) {
        read = parser->get_next_read();
        seq = read.sequence;

        bool valid_read = true;
        valid_read = check_and_normalize_read(seq);
        if (valid_read) {
            BoundedCounterType minval = get_min_count(seq);
            BoundedCounterType maxval = get_max_count(seq);

            mmt->add_min(read_num, minval);
            mmt->add_max(read_num, maxval);
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

    delete parser;

    return mmt;
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
      seq = read.sequence;

      int numPos = seq.length() - _ksize + 1;

      for (int i = 0; i < numPos; i++)  {
         string kmer = seq.substr(i, _ksize);
         outfile << (int)get_count(kmer.c_str()) << " ";
      }
      outfile << endl;
   }

   delete parser;

   outfile.close();
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

HashIntoType * 
CountingHash::abundance_distribution(read_parsers::IParser * parser,
				     Hashbits * tracking)
{
  HashIntoType * dist = new HashIntoType[MAX_BIGCOUNT + 1];
  HashIntoType i;
  
  for (i = 0; i <= MAX_BIGCOUNT; i++) {
    dist[i] = 0;
  }

  Read read;
  
  string name;
  string seq;
  unsigned long long read_num = 0;

  // if not, could lead to overflow.
  assert(sizeof(BoundedCounterType) == 2);

  while(!parser->is_complete()) {
    read = parser->get_next_read();
    seq = read.sequence;

    if (check_and_normalize_read(seq)) {
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
  }
  return dist;
}


HashIntoType * CountingHash::abundance_distribution(std::string filename,
						    Hashbits * tracking)
{
  IParser* parser = IParser::get_parser(filename.c_str());

  return abundance_distribution(parser, tracking);
}

HashIntoType * CountingHash::fasta_count_kmers_by_position(const std::string &inputfile,
					     const unsigned int max_read_len,
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

        seq = read.sequence;
        bool valid_read = true;
         
        valid_read = check_and_normalize_read(seq);

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

          read_num += 1;

          // run callback, if specified
          if (read_num % CALLBACK_PERIOD == 0 && callback) {
             try {
                callback("fasta_file_count_kmers_by_position", callback_data, read_num, 0);
             } catch (...) {
                throw;
             }
          }
   } // while reads

   delete parser;

   return counts;
}

void CountingHash::fasta_dump_kmers_by_abundance(const std::string &inputfile,
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
        seq = read.sequence;

        valid_read = check_and_normalize_read(seq);

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

        read_num += 1;

        // run callback, if specified
        if (read_num % CALLBACK_PERIOD == 0 && callback) {
            try {
                callback("fasta_file_dump_kmers_by_abundance", callback_data, read_num, 0);
            } catch (...) {
                throw;
            }
        }
    } // while reads

    delete parser;
}

void CountingHash::save(std::string outfilename)
{
  CountingHashFile::save(outfilename, *this);
}

void CountingHash::load(std::string infilename)
{
  CountingHashFile::load(infilename, *this);
}

void CountingHash::get_kadian_count(const std::string &s,
				    BoundedCounterType &kadian,
				    unsigned int nk)
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
  unsigned int kpos = nk*_ksize;
  
  if (counts.size() < kpos) {
    kadian = 0;

    return;
  }

  sort(counts.begin(), counts.end());
  kadian = counts[kpos - 1];

#if 0
  std::cout << "k " << kpos << ": ";
  for (unsigned int i = 0; i < counts.size(); i++) {
    std::cout << i << "-" << counts[i] << " ";
  }
  std::cout << "\n";
#endif // 0
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
    seq = read.sequence;

    if (check_and_normalize_read(seq)) {
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

  delete parser;

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
    seq = read.sequence;

    if (check_and_normalize_read(seq)) {
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

  delete parser;

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
  if (!check_and_normalize_read(seq)) {
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


unsigned int CountingHash::trim_below_abundance(std::string seq,
						BoundedCounterType max_abund)
  const
{
  if (!check_and_normalize_read(seq)) {
    return 0;
  }

  KMerIterator kmers(seq.c_str(), _ksize);

  SeenSet path;

  HashIntoType kmer;

  if (kmers.done()) { return 0; }
  kmer = kmers.next();

  if (kmers.done() || get_count(kmer) > max_abund) {
    return 0;
  }

  unsigned int i = _ksize;
  while (!kmers.done()) {
    kmer = kmers.next();

    if (get_count(kmer) > max_abund) {
      return i;
    }
    i++;
  }

  return seq.length();
}


void CountingHashFile::load(const std::string &infilename, CountingHash &ht)
{
   std::string filename(infilename);
   int found = filename.find_last_of(".");
   std::string type = filename.substr(found+1);

   if (type == "gz") { CountingHashGzFileReader(filename, ht); }
   else { CountingHashFileReader(filename, ht); }
}


void CountingHashFile::save(const std::string &outfilename, const CountingHash &ht)
{
   std::string filename(outfilename);
   int found = filename.find_last_of(".");
   std::string type = filename.substr(found+1);

   if (type == "gz") { CountingHashGzFileWriter(filename, ht); }
   else { CountingHashFileWriter(filename, ht); }
}


CountingHashFileReader::CountingHashFileReader(const std::string &infilename, CountingHash &ht)
{
  if (ht._counts) {
    for (unsigned int i = 0; i < ht._n_tables; i++) {
      delete ht._counts[i]; ht._counts[i] = NULL;
    }
    delete ht._counts; ht._counts = NULL;
  }
  ht._tablesizes.clear();
  
  unsigned int save_ksize = 0;
  unsigned char save_n_tables = 0;
  unsigned long long save_tablesize = 0;
  unsigned char version, ht_type, use_bigcount;

  ifstream infile(infilename.c_str(), ios::binary);
  assert(infile.is_open());

  infile.read((char *) &version, 1);
  infile.read((char *) &ht_type, 1);
  assert(version == SAVED_FORMAT_VERSION);
  assert(ht_type == SAVED_COUNTING_HT);

  infile.read((char *) &use_bigcount, 1);
  infile.read((char *) &save_ksize, sizeof(save_ksize));
  infile.read((char *) &save_n_tables, sizeof(save_n_tables));

  ht._ksize = (WordLength) save_ksize;
  ht._n_tables = (unsigned int) save_n_tables;
  ht._init_bitstuff();

  ht._use_bigcount = use_bigcount;

  ht._counts = new Byte*[ht._n_tables];
  for (unsigned int i = 0; i < ht._n_tables; i++) {
    HashIntoType tablesize;

    infile.read((char *) &save_tablesize, sizeof(save_tablesize));

    tablesize = (HashIntoType) save_tablesize;
    ht._tablesizes.push_back(tablesize);

    ht._counts[i] = new Byte[tablesize];

    unsigned long long loaded = 0;
    while (loaded != tablesize) {
      infile.read((char *) ht._counts[i], tablesize - loaded);
      loaded += infile.gcount();	// do I need to do this loop?
    }
  }

  HashIntoType n_counts = 0;
  infile.read((char *) &n_counts, sizeof(n_counts));

  if (n_counts) {
    ht._bigcounts.clear();

    HashIntoType kmer;
    BoundedCounterType count;

    for (HashIntoType n = 0; n < n_counts; n++) {
      infile.read((char *) &kmer, sizeof(kmer));
      infile.read((char *) &count, sizeof(count));
      ht._bigcounts[kmer] = count;
    }
  }

  infile.close();
}

CountingHashGzFileReader::CountingHashGzFileReader(const std::string &infilename, CountingHash &ht)
{
  if (ht._counts) {
    for (unsigned int i = 0; i < ht._n_tables; i++) {
      delete ht._counts[i]; ht._counts[i] = NULL;
    }
    delete ht._counts; ht._counts = NULL;
  }
  ht._tablesizes.clear();
  
  unsigned int save_ksize = 0;
  unsigned char save_n_tables = 0;
  unsigned long long save_tablesize = 0;
  unsigned char version, ht_type, use_bigcount;

  gzFile infile = gzopen(infilename.c_str(), "rb");

  gzread(infile, (char *) &version, 1);
  gzread(infile, (char *) &ht_type, 1);
  assert(version == SAVED_FORMAT_VERSION);
  assert(ht_type == SAVED_COUNTING_HT);

  gzread(infile, (char *) &use_bigcount, 1);
  gzread(infile, (char *) &save_ksize, sizeof(save_ksize));
  gzread(infile, (char *) &save_n_tables, sizeof(save_n_tables));

  ht._ksize = (WordLength) save_ksize;
  ht._n_tables = (unsigned int) save_n_tables;
  ht._init_bitstuff();

  ht._use_bigcount = use_bigcount;

  ht._counts = new Byte*[ht._n_tables];
  for (unsigned int i = 0; i < ht._n_tables; i++) {
    HashIntoType tablesize;

    gzread(infile, (char *) &save_tablesize, sizeof(save_tablesize));

    tablesize = (HashIntoType) save_tablesize;
    ht._tablesizes.push_back(tablesize);

    ht._counts[i] = new Byte[tablesize];

    unsigned long long loaded = 0;
    while (loaded != tablesize) {
      loaded += gzread(infile, (char *) ht._counts[i], tablesize - loaded);
    }
  }

  HashIntoType n_counts = 0;
  gzread(infile, (char *) &n_counts, sizeof(n_counts));

  if (n_counts) {
    ht._bigcounts.clear();

    HashIntoType kmer;
    BoundedCounterType count;

    for (HashIntoType n = 0; n < n_counts; n++) {
      gzread(infile, (char *) &kmer, sizeof(kmer));
      gzread(infile, (char *) &count, sizeof(count));
      ht._bigcounts[kmer] = count;
    }
  }

  gzclose(infile);
}

CountingHashFileWriter::CountingHashFileWriter(const std::string &outfilename, const CountingHash &ht)
{
  assert(ht._counts[0]);

  unsigned int save_ksize = ht._ksize;
  unsigned char save_n_tables = ht._n_tables;
  unsigned long long save_tablesize;

  ofstream outfile(outfilename.c_str(), ios::binary);

  unsigned char version = SAVED_FORMAT_VERSION;
  outfile.write((const char *) &version, 1);

  unsigned char ht_type = SAVED_COUNTING_HT;
  outfile.write((const char *) &ht_type, 1);

  unsigned char use_bigcount = 0;
  if (ht._use_bigcount) {
    use_bigcount = 1;
  }
  outfile.write((const char *) &use_bigcount, 1);

  outfile.write((const char *) &save_ksize, sizeof(save_ksize));
  outfile.write((const char *) &save_n_tables, sizeof(save_n_tables));

  for (unsigned int i = 0; i < save_n_tables; i++) {
    save_tablesize = ht._tablesizes[i];

    outfile.write((const char *) &save_tablesize, sizeof(save_tablesize));
    outfile.write((const char *) ht._counts[i], save_tablesize);
  }

  HashIntoType n_counts = ht._bigcounts.size();
  outfile.write((const char *) &n_counts, sizeof(n_counts));

  if (n_counts) {
    KmerCountMap::const_iterator it = ht._bigcounts.begin();
    
    for (; it != ht._bigcounts.end(); it++) {
      outfile.write((const char *) &it->first, sizeof(it->first));
      outfile.write((const char *) &it->second, sizeof(it->second));
    }
  }

  outfile.close();
}

CountingHashGzFileWriter::CountingHashGzFileWriter(const std::string &outfilename, const CountingHash &ht)
{
  assert(ht._counts[0]);

  unsigned int save_ksize = ht._ksize;
  unsigned char save_n_tables = ht._n_tables;
  unsigned long long save_tablesize;

  gzFile outfile = gzopen(outfilename.c_str(), "wb");

  unsigned char version = SAVED_FORMAT_VERSION;
  gzwrite(outfile, (const char *) &version, 1);

  unsigned char ht_type = SAVED_COUNTING_HT;
  gzwrite(outfile, (const char *) &ht_type, 1);

  unsigned char use_bigcount = 0;
  if (ht._use_bigcount) {
    use_bigcount = 1;
  }
  gzwrite(outfile, (const char *) &use_bigcount, 1);

  gzwrite(outfile, (const char *) &save_ksize, sizeof(save_ksize));
  gzwrite(outfile, (const char *) &save_n_tables, sizeof(save_n_tables));

  for (unsigned int i = 0; i < save_n_tables; i++) {
    save_tablesize = ht._tablesizes[i];

    gzwrite(outfile, (const char *) &save_tablesize, sizeof(save_tablesize));
    gzwrite(outfile, (const char *) ht._counts[i], save_tablesize);
  }

  HashIntoType n_counts = ht._bigcounts.size();
  gzwrite(outfile, (const char *) &n_counts, sizeof(n_counts));

  if (n_counts) {
    KmerCountMap::const_iterator it = ht._bigcounts.begin();
    
    for (; it != ht._bigcounts.end(); it++) {
      gzwrite(outfile, (const char *) &it->first, sizeof(it->first));
      gzwrite(outfile, (const char *) &it->second, sizeof(it->second));
    }
  }

  gzclose(outfile);
}

void CountingHash::collect_high_abundance_kmers(const std::string &filename,
						unsigned int lower_count,
						unsigned int upper_count,
						SeenSet& found_kmers)
{
  unsigned long long total_reads = 0;

  IParser* parser = IParser::get_parser(filename.c_str());
  Read read;

  string currName = "";
  string currSeq = "";

  //
  // iterate through the FASTA file & consume the reads, until we hit
  // upper_count.
  //

  bool done = false;
  while(!parser->is_complete() && !done)  {
    read = parser->get_next_read();
    currSeq = read.sequence;
    currName = read.name; 

    // do we want to process it?
    if (check_and_normalize_read(currSeq)) {
      const char * sp = currSeq.c_str();

      KMerIterator kmers(sp, _ksize);
      HashIntoType kmer;

      while(!kmers.done()) {
	kmer = kmers.next();

	count(kmer);
	if (get_count(kmer) >= upper_count) {
	  done = true;
	}
      }
    }
	       
    // increment read number
    total_reads++;

    if (total_reads % 100000 == 0) {
      std::cout << "..." << total_reads << "\n";
    }
  }

  delete parser; parser = NULL;

  unsigned long long stop_at_read = total_reads;

  //
  // go back through the file again, and store all k-mers >= lower_count
  //

  parser = IParser::get_parser(filename.c_str());

  total_reads = 0;
  while(!parser->is_complete() && total_reads != stop_at_read)  {
    read = parser->get_next_read();
    currSeq = read.sequence;
    currName = read.name; 

    // do we want to process it?
    if (check_and_normalize_read(currSeq)) {
      const char * sp = currSeq.c_str();

      KMerIterator kmers(sp, _ksize);
      HashIntoType kmer;

      while(!kmers.done()) {
	kmer = kmers.next();

	if (get_count(kmer) >= lower_count) {
	  found_kmers.insert(kmer);
	}
      }
    }
	       
    // increment read number
    total_reads++;

    if (total_reads % 100000 == 0) {
      std::cout << "... x 2 " << total_reads << "\n";
    }
  }
  delete parser; parser = NULL;
}

/* vim: set ft=cpp ts=8 sts=4 sw=4 et tw=79 */
