#include "hashtable.hh"
#include "hashbits.hh"
#include "parsers.hh"

#define MAX_KEEPER_SIZE int(1e6)

using namespace std;
using namespace khmer;

void Hashbits::save(std::string outfilename)
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

  outfile.write((const char *) &save_ksize, sizeof(save_ksize));
  outfile.write((const char *) &save_n_tables, sizeof(save_n_tables));

  for (unsigned int i = 0; i < _n_tables; i++) {
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
    for (unsigned int i = 0; i < _n_tables; i++) {
      delete _counts[i]; _counts[i] = NULL;
    }
    delete _counts; _counts = NULL;
  }
  _tablesizes.clear();
  
  unsigned int save_ksize = 0;
  unsigned char save_n_tables = 0;
  unsigned long long save_tablesize = 0;
  unsigned char version, ht_type;

  ifstream infile(infilename.c_str(), ios::binary);
  assert(infile.is_open());

  infile.read((char *) &version, 1);
  infile.read((char *) &ht_type, 1);
  assert(version == SAVED_FORMAT_VERSION);
  assert(ht_type == SAVED_COUNTING_HT);

  infile.read((char *) &save_ksize, sizeof(save_ksize));
  infile.read((char *) &save_n_tables, sizeof(save_n_tables));

  _ksize = (WordLength) save_ksize;
  _n_tables = (unsigned int) save_n_tables;
  _init_bitstuff();

  _counts = new Byte*[_n_tables];
  for (unsigned int i = 0; i < _n_tables; i++) {
    HashIntoType tablesize;
    unsigned long long tablebytes;

    infile.read((char *) &save_tablesize, sizeof(save_tablesize));

    tablesize = (HashIntoType) save_tablesize;
    _tablesizes.push_back(tablesize);

    tablebytes = tablesize / 8 + 1;
    _counts[i] = new Byte[tablebytes];

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

void Hashbits::calc_connected_graph_size(const HashIntoType kmer_f,
					 const HashIntoType kmer_r,
					 unsigned long long& count,
					 SeenSet& keeper,
					 const unsigned long long threshold,
					 bool break_on_circum)
const
{
  HashIntoType kmer = uniqify_rc(kmer_f, kmer_r);
  const BoundedCounterType val = get_count(kmer);

  if (val == 0) {
    return;
  }

  // have we already seen me? don't count; exit.
  if (set_contains(keeper, kmer)) {
    return;
  }

  // is this in stop_tags?
  if (set_contains(stop_tags, kmer)) {
    return;
  }

  // keep track of both seen kmers, and counts.
  keeper.insert(kmer);

  // is this a high-circumference k-mer? if so, don't count it; get outta here!
  if (break_on_circum && \
      kmer_degree(kmer_f, kmer_r) > 4) {
    return;
  }

  count += 1;

  // are we past the threshold? truncate search.
  if (threshold && count >= threshold) {
    return;
  }

  // otherwise, explore in all directions.

  // NEXT.

  HashIntoType f, r;
  const unsigned int rc_left_shift = _ksize*2 - 2;

  f = next_f(kmer_f, 'A');
  r = next_r(kmer_r, 'A');
  calc_connected_graph_size(f, r, count, keeper, threshold, break_on_circum);

  f = next_f(kmer_f, 'C');
  r = next_r(kmer_r, 'C');
  calc_connected_graph_size(f, r, count, keeper, threshold, break_on_circum);

  f = next_f(kmer_f, 'G');
  r = next_r(kmer_r, 'G');
  calc_connected_graph_size(f, r, count, keeper, threshold, break_on_circum);

  f = next_f(kmer_f, 'T');
  r = next_r(kmer_r, 'T');
  calc_connected_graph_size(f, r, count, keeper, threshold, break_on_circum);

  // PREVIOUS.

  
  r = prev_r(kmer_r, 'A');
  f = prev_f(kmer_f, 'A');
  calc_connected_graph_size(f, r, count, keeper, threshold, break_on_circum);

  r = prev_r(kmer_r, 'C');
  f = prev_f(kmer_f, 'C');
  calc_connected_graph_size(f, r, count, keeper, threshold, break_on_circum);

  r = prev_r(kmer_r, 'G');
  f = prev_f(kmer_f, 'G');
  calc_connected_graph_size(f, r, count, keeper, threshold, break_on_circum);

  r = prev_r(kmer_r, 'T');
  f = prev_f(kmer_f, 'T');
  calc_connected_graph_size(f, r, count, keeper, threshold, break_on_circum);
}

void Hashbits::save_tagset(std::string outfilename)
{
  ofstream outfile(outfilename.c_str(), ios::binary);
  const unsigned int tagset_size = all_tags.size();

  HashIntoType * buf = new HashIntoType[tagset_size];

  outfile.write((const char *) &tagset_size, sizeof(tagset_size));
  outfile.write((const char *) &_tag_density, sizeof(_tag_density));

  unsigned int i = 0;
  for (SeenSet::iterator pi = all_tags.begin(); pi != all_tags.end();
	 pi++, i++) {
    buf[i] = *pi;
  }

  outfile.write((const char *) buf, sizeof(HashIntoType) * tagset_size);
  outfile.close();

  delete buf;
}

void Hashbits::load_tagset(std::string infilename, bool clear_tags)
{
  ifstream infile(infilename.c_str(), ios::binary);
  assert(infile.is_open());

  if (clear_tags) {
    all_tags.clear();
  }

  unsigned int tagset_size = 0;
  infile.read((char *) &tagset_size, sizeof(tagset_size));
  infile.read((char *) &_tag_density, sizeof(_tag_density));

  HashIntoType * buf = new HashIntoType[tagset_size];

  infile.read((char *) buf, sizeof(HashIntoType) * tagset_size);

  for (unsigned int i = 0; i < tagset_size; i++) {
    all_tags.insert(buf[i]);
  }

  delete buf;
}

unsigned int Hashbits::kmer_degree(HashIntoType kmer_f, HashIntoType kmer_r)
const
{
  unsigned int neighbors = 0;

  const unsigned int rc_left_shift = _ksize*2 - 2;

  HashIntoType f, r;

  // NEXT.
  f = next_f(kmer_f, 'A');
  r = next_r(kmer_r, 'A');
  if (get_count(uniqify_rc(f, r))) { neighbors++; }
	  
  f = next_f(kmer_f, 'C');
  r = next_r(kmer_r, 'C');
  if (get_count(uniqify_rc(f, r))) { neighbors++; }

  f = next_f(kmer_f, 'G');
  r = next_r(kmer_r, 'G');
  if (get_count(uniqify_rc(f, r))) { neighbors++; }

  f = next_f(kmer_f, 'T');
  r = next_r(kmer_r, 'T');
  if (get_count(uniqify_rc(f, r))) { neighbors++; }

  // PREVIOUS.
  r = prev_r(kmer_r, 'A');
  f = prev_f(kmer_f, 'A');
  if (get_count(uniqify_rc(f, r))) { neighbors++; }

  r = prev_r(kmer_r, 'C');
  f = prev_f(kmer_f, 'C');
  if (get_count(uniqify_rc(f, r))) { neighbors++; }
    
  r = prev_r(kmer_r, 'G');
  f = prev_f(kmer_f, 'G');
  if (get_count(uniqify_rc(f, r))) { neighbors++; }

  r = prev_r(kmer_r, 'T');
  f = prev_f(kmer_f, 'T');
  if (get_count(uniqify_rc(f, r))) { neighbors++; }

  return neighbors;
}


//
// consume_fasta_and_tag: consume a FASTA file of reads, tagging reads every
//     so often.
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

    // n_consumed += this_n_consumed;

    if (check_read(seq)) {	// process?
      bool is_new_kmer;

      KMerIterator kmers(seq.c_str(), _ksize);
      HashIntoType kmer;

      unsigned int since = _tag_density / 2 + 1;

      while(!kmers.done()) {
	kmer = kmers.next();

	is_new_kmer = (bool) !get_count(kmer);
	if (is_new_kmer) {
	  count(kmer);
	  n_consumed++;
	}

	if (!is_new_kmer && set_contains(all_tags, kmer)) {
	  since = 1;
	} else {
	  since++;
	}

	if (since >= _tag_density) {
	  all_tags.insert(kmer);
	  since = 1;
	}
      }

      if (since >= _tag_density/2 - 1) {
	all_tags.insert(kmer);	// insert the last k-mer, too.
      }
    }
	       
    // reset the sequence info, increment read number
    total_reads++;

    // run callback, if specified
    if (total_reads % CALLBACK_PERIOD == 0 && callback) {
      std::cout << "n tags: " << all_tags.size() << "\n";
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
// consume_fasta_and_tag_with_stoptags: consume a FASTA file of reads,
//     tagging reads every so often.  Do not insert matches to stoptags,
//     and join the tags across those gaps.
//

void Hashbits::consume_fasta_and_tag_with_stoptags(const std::string &filename,
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

  SeenSet read_tags;

  //
  // iterate through the FASTA file & consume the reads.
  //

  while(!parser->is_complete())  {
    read = parser->get_next_read();
    seq = read.seq;

    read_tags.clear();

    // n_consumed += this_n_consumed;

    if (check_read(seq)) {	// process?
      bool is_new_kmer;
      KMerIterator kmers(seq.c_str(), _ksize);

      HashIntoType kmer, last_kmer;
      bool is_first_kmer = true;

      unsigned int since = _tag_density / 2 + 1;
      while (!kmers.done()) {
	kmer = kmers.next();

	if (!set_contains(stop_tags, kmer)) { // NOT a stop tag... ok.
	  is_new_kmer = (bool) !get_count(kmer);
	  if (is_new_kmer) {
	    count(kmer);
	    n_consumed++;
	  }

	  if (!is_new_kmer && set_contains(all_tags, kmer)) {
	    read_tags.insert(kmer);
	    since = 1;
	  } else {
	    since++;
	  }

	  if (since >= _tag_density) {
	    all_tags.insert(kmer);
	    read_tags.insert(kmer);
	    since = 1;
	  }
	} else {		// stop tag!  do not insert, but connect.
	  // before first tag insertion; insert last kmer.
	  if (!is_first_kmer && read_tags.size() == 0) {
	    read_tags.insert(last_kmer);
	    all_tags.insert(last_kmer);
	  }
	  
	  since = _tag_density - 1; // insert next kmer, too.
	}

	last_kmer = kmer;
	is_first_kmer = false;
      }

      if (!set_contains(stop_tags, kmer)) { // NOT a stop tag... ok.
	is_new_kmer = (bool) !get_count(kmer);
	if (is_new_kmer) {
	  count(kmer);
	  n_consumed++;
	}

	if (since >= _tag_density/2 - 1) {
	  all_tags.insert(kmer);	// insert the last k-mer, too.
	  read_tags.insert(kmer);
	}
      }
    }

    if (read_tags.size() > 1) {
      partition->assign_partition_id(*(read_tags.begin()), read_tags);
    }
	       
    // reset the sequence info, increment read number
    total_reads++;

    // run callback, if specified
    if (total_reads % CALLBACK_PERIOD == 0 && callback) {
      std::cout << "n tags: " << all_tags.size() << "\n";
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

static PartitionID _parse_partition_id(string name)
{
  PartitionID p = 0;
  const char * s = name.c_str() + name.length() - 1;
  assert(*(s + 1) == (unsigned int) NULL);

  while(*s != '\t' && s >= name.c_str()) {
    s--;
  }

  if (*s == '\t') {
    p = (PartitionID) atoi(s + 1);
  } else {
    cerr << "consume_partitioned_fasta barfed on read "  << name << "\n";
    assert(0);
  }

  return p;
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

    if (check_read(seq)) {
      // First, figure out what the partition is (if non-zero), and save that.
      PartitionID p = _parse_partition_id(read.name);

      // Then consume the sequence
      n_consumed += consume_string(seq); // @CTB why are we doing this?

      // Next, compute the tag & set the partition, if nonzero
      HashIntoType kmer = _hash(seq.c_str(), _ksize);
      all_tags.insert(kmer);
      if (p > 0) {
	partition->set_partition_id(kmer, p);
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

void Hashbits::filter_if_present(const std::string infilename,
				 const std::string outputfile,
				 CallbackFn callback,
				 void * callback_data)
{
  IParser* parser = IParser::get_parser(infilename);
  ofstream outfile(outputfile.c_str());

  unsigned int total_reads = 0;
  unsigned int reads_kept = 0;

  Read read;
  string seq;

  std::string first_kmer;
  HashIntoType kmer;

  while(!parser->is_complete()) {
    read = parser->get_next_read();
    seq = read.seq;

    if (check_read(seq)) {
      KMerIterator kmers(seq.c_str(), _ksize);
      bool keep = true;

      while (!kmers.done()) {
	kmer = kmers.next();
	if (get_count(kmer)) {
	  keep = false;
	  break;
	}
      }

      if (keep) {
	outfile << ">" << read.name;
	outfile << "\n" << seq << "\n";
	reads_kept++;
      }
	       
      total_reads++;

      // run callback, if specified
      if (total_reads % CALLBACK_PERIOD == 0 && callback) {
	try {
	  callback("filter_if_present", callback_data,total_reads, reads_kept);
	} catch (...) {
	  delete parser; parser = NULL;
	  outfile.close();
	  throw;
	}
      }
    }
  }

  delete parser; parser = NULL;

  return;
}


unsigned int Hashbits::count_kmers_within_radius(HashIntoType kmer_f,
						 HashIntoType kmer_r,
						 unsigned int radius,
						 unsigned int max_count,
						 const SeenSet * seen)
const
{
  HashIntoType f, r;
  NodeQueue node_q;
  std::queue<unsigned int> breadth_q;
  unsigned int cur_breadth = 0;
  unsigned int breadth = 0;

  const unsigned int rc_left_shift = _ksize*2 - 2;
  unsigned int total = 0;

  SeenSet keeper;		// keep track of traversed kmers
  if (seen) { keeper = *seen; }

  // start breadth-first search.

  node_q.push(kmer_f);
  node_q.push(kmer_r);
  breadth_q.push(0);

  while(!node_q.empty()) {
    kmer_f = node_q.front();
    node_q.pop();
    kmer_r = node_q.front();
    node_q.pop();
    breadth = breadth_q.front();
    breadth_q.pop();

    if (breadth > radius) {
      break;
    }

    HashIntoType kmer = uniqify_rc(kmer_f, kmer_r);
    if (set_contains(keeper, kmer)) {
      continue;
    }

    // keep track of seen kmers
    keeper.insert(kmer);
    total++;

    if (max_count && total > max_count) {
      break;
    }

    assert(breadth >= cur_breadth); // keep track of watermark, for debugging.
    if (breadth > cur_breadth) { cur_breadth = breadth; }

    //
    // Enqueue next set of nodes.
    //

    // NEXT.
    f = next_f(kmer_f, 'A');
    r = next_r(kmer_r, 'A');
    if (get_count(uniqify_rc(f,r)) && !set_contains(keeper, uniqify_rc(f, r))){
      node_q.push(f); node_q.push(r);
      breadth_q.push(breadth + 1);
    }

    f = next_f(kmer_f, 'C');
    r = next_r(kmer_r, 'C');
    if (get_count(uniqify_rc(f,r)) && !set_contains(keeper, uniqify_rc(f, r))){
      node_q.push(f); node_q.push(r);
      breadth_q.push(breadth + 1);
    }

    f = next_f(kmer_f, 'G');
    r = next_r(kmer_r, 'G');
    if (get_count(uniqify_rc(f,r)) && !set_contains(keeper, uniqify_rc(f, r))){
      node_q.push(f); node_q.push(r);
      breadth_q.push(breadth + 1);
    }

    f = next_f(kmer_f, 'T');
    r = next_r(kmer_r, 'T');
    if (get_count(uniqify_rc(f,r)) && !set_contains(keeper, uniqify_rc(f, r))){
      node_q.push(f); node_q.push(r);
      breadth_q.push(breadth + 1);
    }

    // PREVIOUS.
    r = prev_r(kmer_r, 'A');
    f = prev_f(kmer_f, 'A');
    if (get_count(uniqify_rc(f,r)) && !set_contains(keeper, uniqify_rc(f, r))){
      node_q.push(f); node_q.push(r);
      breadth_q.push(breadth + 1);
    }

    r = prev_r(kmer_r, 'C');
    f = prev_f(kmer_f, 'C');
    if (get_count(uniqify_rc(f,r)) && !set_contains(keeper, uniqify_rc(f, r))){
      node_q.push(f); node_q.push(r);
      breadth_q.push(breadth + 1);
    }
    
    r = prev_r(kmer_r, 'G');
    f = prev_f(kmer_f, 'G');
    if (get_count(uniqify_rc(f,r)) && !set_contains(keeper, uniqify_rc(f, r))){
      node_q.push(f); node_q.push(r);
      breadth_q.push(breadth + 1);
    }

    r = prev_r(kmer_r, 'T');
    f = prev_f(kmer_f, 'T');
    if (get_count(uniqify_rc(f,r)) && !set_contains(keeper, uniqify_rc(f, r))){
      node_q.push(f); node_q.push(r);
      breadth_q.push(breadth + 1);
    }
  }

  return total;
}

unsigned int Hashbits::count_kmers_within_depth(HashIntoType kmer_f,
						HashIntoType kmer_r,
						unsigned int depth,
						unsigned int max_count,
						SeenSet * seen)
const
{
  HashIntoType f, r;
  unsigned int count = 1;

  if (depth == 0) { return 0; }

  const unsigned int rc_left_shift = _ksize*2 - 2;

  seen->insert(uniqify_rc(kmer_f, kmer_r));

  // NEXT.
  f = next_f(kmer_f, 'A');
  r = next_r(kmer_r, 'A');
  if (get_count(uniqify_rc(f,r)) && !set_contains(*seen, uniqify_rc(f, r))) {
    count += count_kmers_within_depth(f, r, depth - 1, max_count - count,
				      seen);
    if (count >= max_count) { return count; }
  }

  f = next_f(kmer_f, 'C');
  r = next_r(kmer_r, 'C');
  if (get_count(uniqify_rc(f,r)) && !set_contains(*seen, uniqify_rc(f, r))) {
    count += count_kmers_within_depth(f, r, depth -1, max_count - count,
				      seen);
    if (count >= max_count) { return count; }
    ;
  }

  f = next_f(kmer_f, 'G');
  r = next_r(kmer_r, 'G');
  if (get_count(uniqify_rc(f,r)) && !set_contains(*seen, uniqify_rc(f, r))) {
    count += count_kmers_within_depth(f, r, depth -1, max_count - count,
				      seen);
    if (count >= max_count) { return count; }
    ;
  }

  f = next_f(kmer_f, 'T');
  r = next_r(kmer_r, 'T');
  if (get_count(uniqify_rc(f,r)) && !set_contains(*seen, uniqify_rc(f, r))) {
    count += count_kmers_within_depth(f, r, depth -1, max_count - count,
				      seen);
    if (count >= max_count) { return count; }
    ;
  }

  // PREVIOUS.
  r = prev_r(kmer_r, 'A');
  f = prev_f(kmer_f, 'A');
  if (get_count(uniqify_rc(f,r)) && !set_contains(*seen, uniqify_rc(f, r))) {
    count += count_kmers_within_depth(f, r, depth -1, max_count - count,
				      seen);
    if (count >= max_count) { return count; }
    ;
  }

  r = prev_r(kmer_r, 'C');
  f = prev_f(kmer_f, 'C');
  if (get_count(uniqify_rc(f,r)) && !set_contains(*seen, uniqify_rc(f, r))) {
    count += count_kmers_within_depth(f, r, depth -1, max_count - count,
				      seen);
    if (count >= max_count) { return count; }
    ;
  }
    
  r = prev_r(kmer_r, 'G');
  f = prev_f(kmer_f, 'G');
  if (get_count(uniqify_rc(f,r)) && !set_contains(*seen, uniqify_rc(f, r))) {
    count += count_kmers_within_depth(f, r, depth -1, max_count - count,
				      seen);
    if (count >= max_count) { return count; }
    ;
  }

  r = prev_r(kmer_r, 'T');
  f = prev_f(kmer_f, 'T');
  if (get_count(uniqify_rc(f,r)) && !set_contains(*seen, uniqify_rc(f, r))) {
    count += count_kmers_within_depth(f, r, depth -1, max_count - count,
				      seen);
    if (count >= max_count) { return count; }
    ;
  }

  return count;
}

unsigned int Hashbits::find_radius_for_volume(HashIntoType kmer_f,
					      HashIntoType kmer_r,
					      unsigned int max_count,
					      unsigned int max_radius)
const
{
  HashIntoType f, r;
  NodeQueue node_q;
  std::queue<unsigned int> breadth_q;
  unsigned int breadth = 0;

  const unsigned int rc_left_shift = _ksize*2 - 2;
  unsigned int total = 0;

  SeenSet keeper;		// keep track of traversed kmers

  // start breadth-first search.

  node_q.push(kmer_f);
  node_q.push(kmer_r);
  breadth_q.push(0);

  while(!node_q.empty()) {
    kmer_f = node_q.front();
    node_q.pop();
    kmer_r = node_q.front();
    node_q.pop();
    breadth = breadth_q.front();
    breadth_q.pop();

    HashIntoType kmer = uniqify_rc(kmer_f, kmer_r);
    if (set_contains(keeper, kmer)) {
      continue;
    }

    // keep track of seen kmers
    keeper.insert(kmer);
    total++;

    if (total >= max_count || breadth >= max_radius) {
      break;
    }

    //
    // Enqueue next set of nodes.
    //

    // NEXT.
    f = next_f(kmer_f, 'A');
    r = next_r(kmer_r, 'A');
    if (get_count(uniqify_rc(f,r)) && !set_contains(keeper, uniqify_rc(f,r))) {
      node_q.push(f); node_q.push(r);
      breadth_q.push(breadth + 1);
    }

    f = next_f(kmer_f, 'C');
    r = next_r(kmer_r, 'C');
    if (get_count(uniqify_rc(f,r)) && !set_contains(keeper, uniqify_rc(f,r))) {
      node_q.push(f); node_q.push(r);
      breadth_q.push(breadth + 1);
    }

    f = next_f(kmer_f, 'G');
    r = next_r(kmer_r, 'G');
    if (get_count(uniqify_rc(f,r)) && !set_contains(keeper, uniqify_rc(f,r))) {
      node_q.push(f); node_q.push(r);
      breadth_q.push(breadth + 1);
    }

    f = next_f(kmer_f, 'T');
    r = next_r(kmer_r, 'T');
    if (get_count(uniqify_rc(f,r)) && !set_contains(keeper, uniqify_rc(f,r))) {
      node_q.push(f); node_q.push(r);
      breadth_q.push(breadth + 1);
    }

    // PREVIOUS.
    r = prev_r(kmer_r, 'A');
    f = prev_f(kmer_f, 'A');
    if (get_count(uniqify_rc(f,r)) && !set_contains(keeper, uniqify_rc(f,r))) {
      node_q.push(f); node_q.push(r);
      breadth_q.push(breadth + 1);
    }

    r = prev_r(kmer_r, 'C');
    f = prev_f(kmer_f, 'C');
    if (get_count(uniqify_rc(f,r)) && !set_contains(keeper, uniqify_rc(f,r))) {
      node_q.push(f); node_q.push(r);
      breadth_q.push(breadth + 1);
    }
    
    r = prev_r(kmer_r, 'G');
    f = prev_f(kmer_f, 'G');
    if (get_count(uniqify_rc(f,r)) && !set_contains(keeper, uniqify_rc(f,r))) {
      node_q.push(f); node_q.push(r);
      breadth_q.push(breadth + 1);
    }

    r = prev_r(kmer_r, 'T');
    f = prev_f(kmer_f, 'T');
    if (get_count(uniqify_rc(f,r)) && !set_contains(keeper, uniqify_rc(f,r))) {
      node_q.push(f); node_q.push(r);
      breadth_q.push(breadth + 1);
    }

    if (node_q.empty()) {
      breadth = max_radius;
      break;
    }
  }

  return breadth;
}

unsigned int Hashbits::count_kmers_on_radius(HashIntoType kmer_f,
					     HashIntoType kmer_r,
					     unsigned int radius,
					     unsigned int max_volume)
const
{
  HashIntoType f, r;
  NodeQueue node_q;
  std::queue<unsigned int> breadth_q;
  unsigned int cur_breadth = 0;
  unsigned int breadth = 0;
  unsigned int count = 0;

  const unsigned int rc_left_shift = _ksize*2 - 2;
  unsigned int total = 0;

  SeenSet keeper;		// keep track of traversed kmers

  // start breadth-first search.

  node_q.push(kmer_f);
  node_q.push(kmer_r);
  breadth_q.push(0);

  while(!node_q.empty()) {
    kmer_f = node_q.front();
    node_q.pop();
    kmer_r = node_q.front();
    node_q.pop();
    breadth = breadth_q.front();
    breadth_q.pop();

    if (breadth > radius) {
      break;
    }

    HashIntoType kmer = uniqify_rc(kmer_f, kmer_r);
    if (set_contains(keeper, kmer)) {
      continue;
    }

    if (breadth == radius) {
      count++;
    }

    // keep track of seen kmers
    keeper.insert(kmer);
    total++;

    if (max_volume && total > max_volume) {
      break;
    }

    assert(breadth >= cur_breadth); // keep track of watermark, for debugging.
    if (breadth > cur_breadth) { cur_breadth = breadth; }

    //
    // Enqueue next set of nodes.
    //

    // NEXT.
    f = next_f(kmer_f, 'A');
    r = next_r(kmer_r, 'A');
    if (get_count(uniqify_rc(f,r)) && !set_contains(keeper, uniqify_rc(f,r))) {
      node_q.push(f); node_q.push(r);
      breadth_q.push(breadth + 1);
    }

    f = next_f(kmer_f, 'C');
    r = next_r(kmer_r, 'C');
    if (get_count(uniqify_rc(f,r)) && !set_contains(keeper, uniqify_rc(f,r))) {
      node_q.push(f); node_q.push(r);
      breadth_q.push(breadth + 1);
    }

    f = next_f(kmer_f, 'G');
    r = next_r(kmer_r, 'G');
    if (get_count(uniqify_rc(f,r)) && !set_contains(keeper, uniqify_rc(f,r))) {
      node_q.push(f); node_q.push(r);
      breadth_q.push(breadth + 1);
    }

    f = next_f(kmer_f, 'T');
    r = next_r(kmer_r, 'T');
    if (get_count(uniqify_rc(f,r)) && !set_contains(keeper, uniqify_rc(f,r))) {
      node_q.push(f); node_q.push(r);
      breadth_q.push(breadth + 1);
    }

    // PREVIOUS.
    r = prev_r(kmer_r, 'A');
    f = prev_f(kmer_f, 'A');
    if (get_count(uniqify_rc(f,r)) && !set_contains(keeper, uniqify_rc(f,r))) {
      node_q.push(f); node_q.push(r);
      breadth_q.push(breadth + 1);
    }

    r = prev_r(kmer_r, 'C');
    f = prev_f(kmer_f, 'C');
    if (get_count(uniqify_rc(f,r)) && !set_contains(keeper, uniqify_rc(f,r))) {
      node_q.push(f); node_q.push(r);
      breadth_q.push(breadth + 1);
    }
    
    r = prev_r(kmer_r, 'G');
    f = prev_f(kmer_f, 'G');
    if (get_count(uniqify_rc(f,r)) && !set_contains(keeper, uniqify_rc(f,r))) {
      node_q.push(f); node_q.push(r);
      breadth_q.push(breadth + 1);
    }

    r = prev_r(kmer_r, 'T');
    f = prev_f(kmer_f, 'T');
    if (get_count(uniqify_rc(f,r)) && !set_contains(keeper, uniqify_rc(f,r))) {
      node_q.push(f); node_q.push(r);
      breadth_q.push(breadth + 1);
    }
  }

  return count;
}

unsigned int Hashbits::trim_on_degree(std::string seq, unsigned int max_degree)
const
{
  if (!check_read(seq)) {
    return 0;

  }

  HashIntoType kmer_f = 0, kmer_r = 0;
  KMerIterator kmers(seq.c_str(), _ksize);

  unsigned int i = _ksize;
  while(!kmers.done()) {
    kmers.next(kmer_f, kmer_r);
  
    if (kmer_degree(kmer_f, kmer_r) > max_degree) {
      return i;
    }
    i++;
  }

  return seq.length();
}

unsigned int Hashbits::trim_on_sodd(std::string seq, unsigned int max_degree)
const
{
  if (!check_read(seq)) {
    return 0;
  }

  const unsigned int RADIUS = 2;
  const unsigned int INCR = 2*RADIUS;
  const char * first_kmer = seq.c_str();

  HashIntoType kmer_f, kmer_r;
  _hash(first_kmer, _ksize, kmer_f, kmer_r);
  if (count_kmers_on_radius(kmer_f, kmer_r, RADIUS, 20) > max_degree) {
    return _ksize - 1;
  }

  for (unsigned int i = INCR; i < seq.length() - _ksize + 1; i += INCR) {
    _hash(first_kmer + i, _ksize, kmer_f, kmer_r);
    if (count_kmers_on_radius(kmer_f, kmer_r, RADIUS, 20) > max_degree) {

      i -= INCR;
      unsigned int pos = 1;

      for (; pos < INCR; pos++) {
	_hash(first_kmer + i + pos, _ksize, kmer_f, kmer_r);
	if (count_kmers_on_radius(kmer_f, kmer_r, RADIUS, 20) > max_degree) {
	  break;
	}
      }

      if (pos == INCR) pos--;
      return i + pos + _ksize - 1;
    }
  }

  return seq.length();
}

unsigned int Hashbits::trim_on_density_explosion(std::string seq,
						 unsigned int radius,
						 unsigned int max_volume)
  const
{
  if (!check_read(seq)) {
    return 0;
  }
  unsigned int count;
  SeenSet path;

  HashIntoType kmer_f = 0, kmer_r = 0;
  SeenSet seen;

  KMerIterator kmers(seq.c_str(), _ksize);

  unsigned int i = _ksize - 2;
  while(!kmers.done()) {
    kmers.next(kmer_f, kmer_r);
    count = count_kmers_within_depth(kmer_f, kmer_r, radius,
				     max_volume, &seen);
    if (count >= max_volume) {
      return i;
    }
    
    i++;
  }

  
  return seq.length();
}

unsigned int Hashbits::trim_on_stoptags(std::string seq) const
{
  if (!check_read(seq)) {
    return 0;
  }

  SeenSet path;
  HashIntoType kmer;

  KMerIterator kmers(seq.c_str(), _ksize);
  
  unsigned int i = _ksize - 2;
  while (!kmers.done()) {
    kmer = kmers.next();
    if (set_contains(stop_tags, kmer)) {
      return i;
    }
    i++;
  }

  return seq.length();
}

void Hashbits::traverse_from_tags(unsigned int distance,
				  unsigned int threshold,
				  unsigned int frequency,
				  CountingHash &counting)
{
  unsigned int i = 0;
  unsigned int n = 0;
  unsigned int count;
  unsigned int n_big = 0;
  SeenSet keeper;

#if VERBOSE_REPARTITION
  std::cout << all_tags.size() << " tags...\n";
#endif // 0
  SeenSet::const_iterator si = all_tags.begin();

  for (; si != all_tags.end(); si++, i++) {
    n++;
    count = traverse_from_kmer(*si, distance, keeper);

    if (count >= threshold) {
      n_big++;
	
      SeenSet::const_iterator ti;
      for (ti = keeper.begin(); ti != keeper.end(); ti++) {
	if (counting.get_count(*ti) > frequency) {
	  stop_tags.insert(*ti);
	} else {
	  counting.count(*ti);
	}
      }
#if VERBOSE_REPARTITION
      std::cout << "traversed from " << n << " tags total; "
		<< n_big << " big; " << keeper.size() << "\n";
#endif // 0
    }
    keeper.clear();

    if (n % 100 == 0) {
#if VERBOSE_REPARTITION
      std::cout << "traversed " << n << " " << n_big << " " <<
	all_tags.size() << " " << stop_tags.size() << "\n";
#endif // 0
    }
  }
}

unsigned int Hashbits::traverse_from_kmer(HashIntoType start,
					  unsigned int radius,
					  SeenSet &keeper)
const
{
  std::string kmer_s = _revhash(start, _ksize);
  HashIntoType kmer, kmer_f, kmer_r;
  kmer = _hash(kmer_s.c_str(), _ksize, kmer_f, kmer_r);

  HashIntoType f, r;
  NodeQueue node_q;
  std::queue<unsigned int> breadth_q;
  unsigned int cur_breadth = 0;
  unsigned int breadth = 0;
  bool is_first_kmer = true;

  const unsigned int rc_left_shift = _ksize*2 - 2;
  unsigned int total = 0;

  // start breadth-first search.

  node_q.push(kmer_f);
  node_q.push(kmer_r);
  breadth_q.push(0);

  while(!node_q.empty()) {
    kmer_f = node_q.front();
    node_q.pop();
    kmer_r = node_q.front();
    node_q.pop();
    breadth = breadth_q.front();
    breadth_q.pop();

    if (breadth > radius) {
      break;
    }

    if (total > MAX_KEEPER_SIZE) {
      break;
    }

    HashIntoType kmer = uniqify_rc(kmer_f, kmer_r);
    if (set_contains(keeper, kmer)) {
      continue;
    }

    if (set_contains(stop_tags, kmer)) {
      continue;
    }

    // keep track of seen kmers
    keeper.insert(kmer);
    total++;

    if (false && !is_first_kmer && set_contains(all_tags, kmer)) {
      continue;
    }

    assert(breadth >= cur_breadth); // keep track of watermark, for debugging.
    if (breadth > cur_breadth) { cur_breadth = breadth; }

    //
    // Enqueue next set of nodes.
    //

    // NEXT.
    f = next_f(kmer_f, 'A');
    r = next_r(kmer_r, 'A');
      
    // f = ((kmer_f << 2) & bitmask) | twobit_repr('A');
    r = next_r(kmer_r, 'A');
    if (get_count(uniqify_rc(f,r)) && !set_contains(keeper, uniqify_rc(f,r))) {
      node_q.push(f); node_q.push(r);
      breadth_q.push(breadth + 1);
    }

    f = next_f(kmer_f, 'C');
    r = next_r(kmer_r, 'C');
    if (get_count(uniqify_rc(f,r)) && !set_contains(keeper, uniqify_rc(f,r))) {
      node_q.push(f); node_q.push(r);
      breadth_q.push(breadth + 1);
    }

    f = next_f(kmer_f, 'G');
    r = next_r(kmer_r, 'G');
    if (get_count(uniqify_rc(f,r)) && !set_contains(keeper, uniqify_rc(f,r))) {
      node_q.push(f); node_q.push(r);
      breadth_q.push(breadth + 1);
    }

    f = next_f(kmer_f, 'T');
    r = next_r(kmer_r, 'T');
    if (get_count(uniqify_rc(f,r)) && !set_contains(keeper, uniqify_rc(f,r))) {
      node_q.push(f); node_q.push(r);
      breadth_q.push(breadth + 1);
    }

    // PREVIOUS.
    r = prev_r(kmer_r, 'A');
    f = prev_f(kmer_f, 'A');
    if (get_count(uniqify_rc(f,r)) && !set_contains(keeper, uniqify_rc(f,r))) {
      node_q.push(f); node_q.push(r);
      breadth_q.push(breadth + 1);
    }

    r = prev_r(kmer_r, 'C');
    f = prev_f(kmer_f, 'C');
    if (get_count(uniqify_rc(f,r)) && !set_contains(keeper, uniqify_rc(f,r))) {
      node_q.push(f); node_q.push(r);
      breadth_q.push(breadth + 1);
    }
    
    r = prev_r(kmer_r, 'G');
    f = prev_f(kmer_f, 'G');
    if (get_count(uniqify_rc(f,r)) && !set_contains(keeper, uniqify_rc(f,r))) {
      node_q.push(f); node_q.push(r);
      breadth_q.push(breadth + 1);
    }

    r = prev_r(kmer_r, 'T');
    f = prev_f(kmer_f, 'T');
    if (get_count(uniqify_rc(f,r)) && !set_contains(keeper, uniqify_rc(f,r))) {
      node_q.push(f); node_q.push(r);
      breadth_q.push(breadth + 1);
    }

    is_first_kmer = false;
  }

  return total;
}

void Hashbits::hitraverse_to_stoptags(std::string filename,
				      CountingHash &counting,
				      unsigned int cutoff)
{
  Read read;
  IParser* parser = IParser::get_parser(filename);
  string name;
  string seq;
  unsigned int read_num = 0;

  while(!parser->is_complete()) {
    read = parser->get_next_read();
    seq = read.seq;

    if (check_read(seq)) {
      for (unsigned int i = 0; i < seq.length() - _ksize + 1; i++) {
	string kmer = seq.substr(i, i + _ksize); // @CTB this wrong!
	HashIntoType kmer_n = _hash(kmer.c_str(), _ksize);
	BoundedCounterType n = counting.get_count(kmer_n);

	if (n >= cutoff) {
	  stop_tags.insert(kmer_n);
	}
      }

      name.clear();
      seq.clear();
    }

    read_num += 1;
  }

#if VERBOSE_REPARTITION
  std::cout << "Inserted " << stop_tags.size() << " stop tags\n";
#endif // 0
}

void Hashbits::load_stop_tags(std::string infilename, bool clear_tags)
{
  ifstream infile(infilename.c_str(), ios::binary);
  assert(infile.is_open());

  if (clear_tags) {
    stop_tags.clear();
  }

  unsigned int tagset_size = 0;
  infile.read((char *) &tagset_size, sizeof(tagset_size));

  HashIntoType * buf = new HashIntoType[tagset_size];

  infile.read((char *) buf, sizeof(HashIntoType) * tagset_size);

  for (unsigned int i = 0; i < tagset_size; i++) {
    stop_tags.insert(buf[i]);
  }

  delete buf;
}

void Hashbits::save_stop_tags(std::string outfilename)
{
  ofstream outfile(outfilename.c_str(), ios::binary);
  const unsigned int tagset_size = stop_tags.size();

  HashIntoType * buf = new HashIntoType[tagset_size];

  outfile.write((const char *) &tagset_size, sizeof(tagset_size));

  unsigned int i = 0;
  for (SeenSet::iterator pi = stop_tags.begin(); pi != stop_tags.end();
	 pi++, i++) {
    buf[i] = *pi;
  }

  outfile.write((const char *) buf, sizeof(HashIntoType) * tagset_size);
  outfile.close();

  delete buf;
}

void Hashbits::print_stop_tags(std::string infilename)
{
  ofstream printfile(infilename.c_str());

  unsigned int i = 0;
  for (SeenSet::iterator pi = stop_tags.begin(); pi != stop_tags.end();
	 pi++, i++) {
    std::string kmer = _revhash(*pi, _ksize);
    printfile << kmer << "\n";
  }
  
  printfile.close();
}

void Hashbits::print_tagset(std::string infilename)
{
  ofstream printfile(infilename.c_str());

  unsigned int i = 0;
  for (SeenSet::iterator pi = all_tags.begin(); pi != all_tags.end();
	 pi++, i++) {
    std::string kmer = _revhash(*pi, _ksize);
    printfile << kmer << "\n";
  }
  
  printfile.close();
}

unsigned int Hashbits::count_and_transfer_to_stoptags(SeenSet &keeper,
						      unsigned int threshold,
						      CountingHash &counting)
{
  unsigned int n_inserted = 0;

  SeenSet::const_iterator ti;
  for (ti = keeper.begin(); ti != keeper.end(); ti++) {
    if (counting.get_count(*ti) >= threshold) {
      stop_tags.insert(*ti);
      n_inserted++;
    } else {
      counting.count(*ti);
    }
  }

  return n_inserted;
}

void Hashbits::traverse_from_reads(std::string filename,
				   unsigned int radius,
				   unsigned int big_threshold,
				   unsigned int transfer_threshold,
				   CountingHash &counting)
{
  unsigned long long total_reads = 0;
  unsigned long long total_stop = 0;

  IParser* parser = IParser::get_parser(filename.c_str());
  Read read;
  SeenSet keeper;

  string seq = "";

  //
  // iterate through the FASTA file & consume the reads.
  //

  while(!parser->is_complete())  {
    read = parser->get_next_read();
    seq = read.seq;

    if (check_read(seq)) {	// process?
      const char * last_kmer = seq.c_str() + seq.length() - _ksize;
      HashIntoType kmer = _hash(last_kmer, _ksize);

      unsigned int n = traverse_from_kmer(kmer, radius, keeper);

      if (n >= big_threshold) {
#if VERBOSE_REPARTITION
	std::cout << "lump: " << n << "; added: " << total_stop << "\n";
#endif
	total_stop += count_and_transfer_to_stoptags(keeper,
						     transfer_threshold,
						     counting);
      }

      keeper.clear();
    }
	       
    // reset the sequence info, increment read number
    total_reads++;

    // run callback, if specified
    if (total_reads % CALLBACK_PERIOD == 0) {
      std::cout << "n reads: " << total_reads << "; n tags: " << stop_tags.size() << "\n";
    }
  }
  delete parser;
}

//
// consume_fasta: consume a FASTA file of reads
//

void Hashbits::consume_fasta_and_traverse(const std::string &filename,
					  unsigned int radius,
					  unsigned int big_threshold,
					  unsigned int transfer_threshold,
					  CountingHash &counting)
{
  unsigned long long total_reads = 0;

  IParser* parser = IParser::get_parser(filename.c_str());
  Read read;

  string seq = "";

  //
  // iterate through the FASTA file & consume the reads.
  //

  while(!parser->is_complete())  {
    read = parser->get_next_read();
    seq = read.seq;

    if (check_read(seq)) {	// process?
      KMerIterator kmers(seq.c_str(), _ksize);

      HashIntoType kmer = 0;
      bool is_first_kmer = true;
      while (!kmers.done()) {
	kmer = kmers.next();

	if (set_contains(stop_tags, kmer)) {
	  break;
	}
	count(kmer);
	is_first_kmer = false;
      }

      if (!is_first_kmer) {	// traverse
	SeenSet keeper;

	unsigned int n = traverse_from_kmer(kmer, radius, keeper);
	if (n >= big_threshold) {
#if VERBOSE_REPARTITION
	  std::cout << "lmp: " << n << "; added: " << stop_tags.size() << "\n";
#endif // VERBOSE_REPARTITION
	  count_and_transfer_to_stoptags(keeper, transfer_threshold, counting);
	}
      }
    }

    // reset the sequence info, increment read number
    total_reads++;

    // run callback, if specified
    if (total_reads % CALLBACK_PERIOD == 0) {
      std::cout << "n reads: " << total_reads << "\n";
    }
  }
  delete parser;
}

void Hashbits::identify_stop_tags_by_position(std::string seq,
					      std::vector<unsigned int> &posns)
const
{
  if (!check_read(seq)) {
    return;
  }

  SeenSet path;
  HashIntoType kmer;

  KMerIterator kmers(seq.c_str(), _ksize);
  
  unsigned int i = 0;
  while(!kmers.done()) {
    kmer = kmers.next();

    if (set_contains(stop_tags, kmer)) {
      posns.push_back(i);
    }
    i++;
  }

  return;
}

void Hashbits::extract_unique_paths(std::string seq,
				    unsigned int min_length,
				    float min_unique_f,
				    std::vector<std::string> &results)
{
  if (seq.size() < min_length) {
    return;
  }

  float max_seen = 1.0 - min_unique_f;

  min_length = min_length - _ksize + 1; // adjust for k-mer size.

  KMerIterator kmers(seq.c_str(), _ksize);
  HashIntoType kmer;

  std::deque<bool> seen_queue;
  unsigned int n_already_seen = 0;
  unsigned int n_kmers = 0;

  // first, put together an array for presence/absence of the k-mer
  // at each given position.
  while (!kmers.done()) {
    kmer = kmers.next();

    if (get_count(kmer)) {
      seen_queue.push_back(true);
      n_already_seen++;
    } else {
      seen_queue.push_back(false);
    }
    n_kmers++;
  }

  // next, run through this array with 'i'.

  unsigned int i = 0;
  while (i < n_kmers - min_length) {
    unsigned int seen_counter, j;

    // For each starting 'i', count the number of 'seen' k-mers in the
    // given window.

    // yes, inefficient n^2 algorithm.  sue me.
    for (seen_counter = 0, j = 0; j < min_length; j++) {
      if (seen_queue[i + j]) {
	seen_counter++;
      }
    }

    // If the fraction seen is small enough to be interesting, suggesting
    // that this, in fact, a "new" window -- extend until it isn't, and
    // then extract.

    assert(j == min_length);
    if ( ((float)seen_counter / (float) j) <= max_seen) {
      unsigned int start = i;

      // extend the window until the end of the sequence...
      while ((start + min_length) < n_kmers) {
	if (seen_queue[start]) {
	  seen_counter--;
	}
	if (seen_queue[start + min_length]) {
	  seen_counter++;
	}
	start++;

	// ...or until we've seen too many of the k-mers.
	if (((float)seen_counter / (float) min_length) > max_seen) {
	  break;
	}
      }

      // adjust for ending point.
      if (start + min_length == n_kmers) {	// potentially decrement twice at end
	if (((float)seen_counter / (float) min_length) > max_seen) {
	  start--;
	}
	start--;
      }
      else {
	start -= 2;
      }

      // ...and now extract the relevant portion of the sequence, and adjust
      // starting pos'n.
      results.push_back(seq.substr(i, start + min_length + _ksize - i));

      i = start + min_length + 1;
    } else {
      i++;
    }
  }
}
