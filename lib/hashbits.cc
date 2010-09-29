#include "hashtable.hh"
#include "hashbits.hh"
#include "parsers.hh"

using namespace std;
using namespace khmer;

void Hashbits::save(std::string outfilename)
{
  assert(_counts[0]);

  unsigned int save_ksize = _ksize;
  unsigned long long save_tablesize;

  ofstream outfile(outfilename.c_str(), ios::binary);

  outfile.write((const char *) &save_ksize, sizeof(save_ksize));

  for (unsigned int i = 0; i < n_tables; i++) {
    save_tablesize = _tablesizes[i];
    unsigned long long tablebytes = save_tablesize / 8 + 1;

    outfile.write((const char *) &save_tablesize, sizeof(save_tablesize));

    outfile.write((const char *) _counts[i], tablebytes);
  }
  outfile.close();
}

void Hashbits::load(std::string infilename)
{
  if (_counts) {
    for (unsigned int i = 0; i < n_tables; i++) {
      delete _counts[i]; _counts[i] = NULL;
    }
    delete _counts; _counts = NULL;
  }
  _tablesizes.clear();
  
  unsigned int save_ksize = 0;
  unsigned long long save_tablesize = 0;

  ifstream infile(infilename.c_str(), ios::binary);
  infile.read((char *) &save_ksize, sizeof(save_ksize));
  _ksize = (WordLength) save_ksize;

  _counts = new BoundedCounterType*[n_tables];
  for (unsigned int i = 0; i < n_tables; i++) {
    HashIntoType tablesize;
    unsigned long long tablebytes;

    infile.read((char *) &save_tablesize, sizeof(save_tablesize));

    tablesize = (HashIntoType) save_tablesize;
    _tablesizes.push_back(tablesize);

    tablebytes = tablesize / 8 + 1;
    _counts[i] = new BoundedCounterType[tablebytes];

    unsigned long long loaded = 0;
    while (loaded != tablebytes) {
      infile.read((char *) _counts[i], tablebytes - loaded);
      loaded += infile.gcount();	// do I need to do this loop?
    }
  }
  infile.close();
}

//////////////////////////////////////////////////////////////////////
// graph stuff

ReadMaskTable * Hashbits::filter_file_connected(const std::string &est,
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

void Hashbits::calc_connected_graph_size(const HashIntoType kmer_f,
					  const HashIntoType kmer_r,
					  unsigned long long& count,
					  SeenSet& keeper,
					  const unsigned long long threshold)
const
{
  HashIntoType kmer = uniqify_rc(kmer_f, kmer_r);
  const BoundedCounterType val = get_count(kmer);

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

void Hashbits::trim_graphs(const std::string infilename,
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

HashIntoType * Hashbits::graphsize_distribution(const unsigned int &max_size)
{
  HashIntoType * p = new HashIntoType[max_size];
  const unsigned char seen = 1 << 7;
  unsigned long long size;

  for (unsigned int i = 0; i < max_size; i++) {
    p[i] = 0;
  }

  for (HashIntoType i = 0; i < _tablesizes[0]; i++) {
    BoundedCounterType count = get_count(i);
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

void Hashbits::save_tagset(std::string outfilename)
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

void Hashbits::load_tagset(std::string infilename)
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

void Hashbits::do_truncated_partition(const std::string infilename,
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
      tagged_kmers.clear();
      for (unsigned int i = 0; i < seq.length() - _ksize + 1;
	   i += TAG_DENSITY) {
	first_kmer = seq.substr(i, _ksize);
	_hash(first_kmer.c_str(), _ksize, kmer_f, kmer_r);
	HashIntoType kmer = uniqify_rc(kmer_f, kmer_r);
	tagged_kmers.insert(kmer);
      }

      for (unsigned int i = 0; i < seq.length() - _ksize + 1;
	   i += TAG_DENSITY) {
	first_kmer = seq.substr(i, _ksize);
	_hash(first_kmer.c_str(), _ksize, kmer_f, kmer_r);

	// find all tagged kmers within range.
	surrender = false;
	partition->find_all_tags(kmer_f, kmer_r, tagged_kmers, surrender, false);
	for (SeenSet::iterator si = tagged_kmers.begin(); si != tagged_kmers.end(); si++) {
	  PartitionID * pmap = partition->partition_map[*si];
	  PartitionID pid = 0;
	  if (pmap) pid = *pmap;
	}

	// assign the partition ID
	HashIntoType kmer = uniqify_rc(kmer_f, kmer_r);
	partition->assign_partition_id(kmer, tagged_kmers, surrender);
	all_tags.insert(kmer);
      }

      // run callback, if specified
      if (total_reads % CALLBACK_PERIOD == 0 && callback) {
	try {
	  callback("do_truncated_partition", callback_data, total_reads,
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

// do_threaded_partition:

void Hashbits::do_threaded_partition(const std::string infilename,
				      CallbackFn callback,
				      void * callback_data)
{
  unsigned int total_reads = 0;

  IParser* parser = IParser::get_parser(infilename);
  Read read;
  string seq;
  bool is_valid;

  HashIntoType kmer_f, kmer_r;

  if (!partition) {
    partition = new SubsetPartition(this);
  }

  while(!parser->is_complete()) {
    // increment read number
    read = parser->get_next_read();
    total_reads++;

    seq = read.seq;

    is_valid = check_read(seq);
    if (is_valid) {
      const char * kmer_s = seq.c_str();
      HashIntoType kmer;
      HashIntoType last_overlap;

      bool found = false;

      for (unsigned int i = 0; i < seq.length() - _ksize + 1; i++) {
	_hash(kmer_s + i, _ksize, kmer_f, kmer_r);
	kmer = uniqify_rc(kmer_f, kmer_r);

	if (get_count(kmer)) {
	  if (!found) {
	    if (all_tags.find(kmer) == all_tags.end()) { // tag first intersect
	      all_tags.insert(kmer);
	      last_overlap = kmer;
	    }
	    found = true;
	  } else {
	    last_overlap = kmer;
	  }
	} else {		// no overlap
	  count(kmer);
	}
      }

      if (found) {		// tag last intersect
	all_tags.insert(last_overlap);
      } else {			// no intersect? insert first kmer.
	_hash(kmer_s, _ksize, kmer_f, kmer_r);
	kmer = uniqify_rc(kmer_f, kmer_r);
	all_tags.insert(kmer);
      }

      // run callback, if specified
      if (total_reads % CALLBACK_PERIOD == 0 && callback) {
	try {
	  callback("do_threaded_partition", callback_data, total_reads,
		   all_tags.size());
	} catch (...) {
	  delete parser;
	  throw;
	}
      }
    }
  }

  delete parser;
}

void Hashbits::connectivity_distribution(const std::string infilename,
					 HashIntoType dist[9],
					 CallbackFn callback,
					 void * callback_data)
{
  const unsigned int rc_left_shift = _ksize*2 - 2;
  unsigned int total_reads = 0;
  for (unsigned int i = 0; i < 9; i++) {
    dist[i] = 0;
  }

  IParser* parser = IParser::get_parser(infilename);
  Read read;
  string seq;
  bool is_valid;

  HashIntoType kmer_f, kmer_r;

  while(!parser->is_complete()) {
    // increment read number
    read = parser->get_next_read();
    total_reads++;

    seq = read.seq;

    is_valid = check_read(seq);
    if (is_valid) {
      const char * kmer_s = seq.c_str();

      for (unsigned int i = 0; i < seq.length() - _ksize + 1; i++) {
	unsigned int neighbors = 0;
	_hash(kmer_s + i, _ksize, kmer_f, kmer_r);

	HashIntoType f, r;

	// NEXT.
	f = ((kmer_f << 2) & bitmask) | twobit_repr('A');
	r = kmer_r >> 2 | (twobit_comp('A') << rc_left_shift);
	if (get_count(uniqify_rc(f, r))) { neighbors++; }
	  
	f = ((kmer_f << 2) & bitmask) | twobit_repr('C');
	r = kmer_r >> 2 | (twobit_comp('C') << rc_left_shift);
	if (get_count(uniqify_rc(f, r))) { neighbors++; }

	f = ((kmer_f << 2) & bitmask) | twobit_repr('G');
	r = kmer_r >> 2 | (twobit_comp('G') << rc_left_shift);
	if (get_count(uniqify_rc(f, r))) { neighbors++; }

	f = ((kmer_f << 2) & bitmask) | twobit_repr('T');
	r = kmer_r >> 2 | (twobit_comp('T') << rc_left_shift);
	if (get_count(uniqify_rc(f, r))) { neighbors++; }

	// PREVIOUS.
	r = ((kmer_r << 2) & bitmask) | twobit_comp('A');
	f = kmer_f >> 2 | (twobit_repr('A') << rc_left_shift);
	if (get_count(uniqify_rc(f, r))) { neighbors++; }

	r = ((kmer_r << 2) & bitmask) | twobit_comp('C');
	f = kmer_f >> 2 | (twobit_repr('C') << rc_left_shift);
	if (get_count(uniqify_rc(f, r))) { neighbors++; }
    
	r = ((kmer_r << 2) & bitmask) | twobit_comp('G');
	f = kmer_f >> 2 | (twobit_repr('G') << rc_left_shift);
	if (get_count(uniqify_rc(f, r))) { neighbors++; }

	r = ((kmer_r << 2) & bitmask) | twobit_comp('T');
	f = kmer_f >> 2 | (twobit_repr('T') << rc_left_shift);
	if (get_count(uniqify_rc(f, r))) { neighbors++; }

	dist[neighbors]++;
      }

      // run callback, if specified
      if (total_reads % CALLBACK_PERIOD == 0 && callback) {
	try {
	  callback("connectivity_dist", callback_data, total_reads, 0);
	} catch (...) {
	  delete parser;
	  throw;
	}
      }
    }
  }

  delete parser;
}

//
// consume_fasta: consume a FASTA file of reads
//

void Hashbits::consume_fasta_and_tag(const std::string &filename,
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
      const char * first_kmer = seq.c_str();
      for (unsigned int i = 0; i < seq.length() - _ksize + 1;
	   i += TAG_DENSITY) {
	HashIntoType kmer = _hash(first_kmer + i, _ksize);
	all_tags.insert(kmer);
      }
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
  delete parser;
}

//
// divide_tags_into_subsets - take all of the tags in 'all_tags', and
//   divide them into subsets (based on starting tag) of <= given size.
//

void Hashbits::divide_tags_into_subsets(unsigned int subset_size,
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

//
// tags_to_map - convert the 'all_tags' set into a TagCountMap, connecting
//    each tag to a number (defaulting to zero).
//

void Hashbits::tags_to_map(TagCountMap& tag_map)
{
  for (SeenSet::const_iterator si = all_tags.begin(); si != all_tags.end();
       si++) {
    tag_map[*si] = 0;
  }
  cout << "TM size: " << tag_map.size() << "\n";
}

//
// discard_tags - remove tags from a TagCountMap if they have fewer than
//   threshold count tags in their partition.  Used to eliminate tags belonging
//   to small partitions.
//

void Hashbits::discard_tags(TagCountMap& tag_map, unsigned int threshold)
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

//
// consume_partitioned_fasta: consume a FASTA file of reads
//

void Hashbits::consume_partitioned_fasta(const std::string &filename,
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

  // reset the master subset partition
  delete partition;
  partition = new SubsetPartition(this);

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
      // First, compute the tag (first k-mer)
      string first_kmer = seq.substr(0, _ksize);
      HashIntoType kmer_f, kmer_r;
      _hash(first_kmer.c_str(), _ksize, kmer_f, kmer_r);

      all_tags.insert(kmer_f);

      // Next, figure out the partition is (if non-zero), and save that.
      const char * s = read.name.c_str() + read.name.length() - 1;
      assert(*(s + 1) == (unsigned int) NULL);

      while(*s != '\t' && s >= read.name.c_str()) {
	s--;
      }

      if (*s == '\t') {
	PartitionID p = (PartitionID) atoi(s + 1);
	if (p > 0) {
	  partition->set_partition_id(kmer_f, p);
	}
      }
    }
	       
    // reset the sequence info, increment read number
    total_reads++;

    // run callback, if specified
    if (total_reads % CALLBACK_PERIOD == 0 && callback) {
      try {
        callback("consume_partitioned_fasta", callback_data, total_reads,
		 n_consumed);
      } catch (...) {
	delete parser;
        throw;
      }
    }
  }

  delete parser;
}
