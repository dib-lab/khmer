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
  unsigned long long save_tablesize = 0;

  ifstream infile(infilename.c_str(), ios::binary);
  infile.read((char *) &save_ksize, sizeof(save_ksize));
  _ksize = (WordLength) save_ksize;
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
  SeenSet::iterator i = keeper.find(kmer);
  if (i != keeper.end()) {
    return;
  }

  // is this in stop_tags?
  i = stop_tags.find(kmer);
  if (i != stop_tags.end()) {
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

void Hashbits::load_tagset(std::string infilename)
{
  ifstream infile(infilename.c_str(), ios::binary);
  all_tags.clear();

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

  return neighbors;
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

    // n_consumed += this_n_consumed;

    if (check_read(seq)) {	// process?
      bool is_new_kmer;
      const char * first_kmer = seq.c_str();
      HashIntoType kmer_f = 0, kmer_r = 0;
      HashIntoType kmer = _hash(first_kmer, _ksize, kmer_f, kmer_r);

      unsigned int since = _tag_density / 2 + 1;
      for (unsigned int i = _ksize; i < seq.length(); i++) {

	is_new_kmer = (bool) !get_count(kmer);
	if (is_new_kmer) {
	  count(kmer);
	  n_consumed++;
	}

	if (!is_new_kmer && all_tags.find(kmer) != all_tags.end()) {
	  since = 1;
	} else {
	  since++;
	}

	if (since >= _tag_density) {
	  all_tags.insert(kmer);
	  since = 1;
	}

	kmer = _next_hash(seq[i], kmer_f, kmer_r);
      }

      is_new_kmer = (bool) !get_count(kmer);
      if (is_new_kmer) {
	count(kmer);
	n_consumed++;
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
      const char * kmer_s = seq.c_str();
      bool keep = true;
      
      for (unsigned int i = 0; i < seq.length() - _ksize + 1; i++) {
	kmer = _hash(kmer_s + i, _ksize);

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
    if (keeper.find(kmer) != keeper.end()) {
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
    f = ((kmer_f << 2) & bitmask) | twobit_repr('A');
    r = kmer_r >> 2 | (twobit_comp('A') << rc_left_shift);
    if (get_count(uniqify_rc(f,r)) && 
	keeper.find(uniqify_rc(f,r)) == keeper.end()) {
      node_q.push(f); node_q.push(r);
      breadth_q.push(breadth + 1);
    }

    f = ((kmer_f << 2) & bitmask) | twobit_repr('C');
    r = kmer_r >> 2 | (twobit_comp('C') << rc_left_shift);
    if (get_count(uniqify_rc(f,r)) && 
	keeper.find(uniqify_rc(f,r)) == keeper.end()) {
      node_q.push(f); node_q.push(r);
      breadth_q.push(breadth + 1);
    }

    f = ((kmer_f << 2) & bitmask) | twobit_repr('G');
    r = kmer_r >> 2 | (twobit_comp('G') << rc_left_shift);
    if (get_count(uniqify_rc(f,r)) && 
	keeper.find(uniqify_rc(f,r)) == keeper.end()) {
      node_q.push(f); node_q.push(r);
      breadth_q.push(breadth + 1);
    }

    f = ((kmer_f << 2) & bitmask) | twobit_repr('T');
    r = kmer_r >> 2 | (twobit_comp('T') << rc_left_shift);
    if (get_count(uniqify_rc(f,r)) && 
	keeper.find(uniqify_rc(f,r)) == keeper.end()) {
      node_q.push(f); node_q.push(r);
      breadth_q.push(breadth + 1);
    }

    // PREVIOUS.
    r = ((kmer_r << 2) & bitmask) | twobit_comp('A');
    f = kmer_f >> 2 | (twobit_repr('A') << rc_left_shift);
    if (get_count(uniqify_rc(f,r)) && 
	keeper.find(uniqify_rc(f,r)) == keeper.end()) {
      node_q.push(f); node_q.push(r);
      breadth_q.push(breadth + 1);
    }

    r = ((kmer_r << 2) & bitmask) | twobit_comp('C');
    f = kmer_f >> 2 | (twobit_repr('C') << rc_left_shift);
    if (get_count(uniqify_rc(f,r)) && 
	keeper.find(uniqify_rc(f,r)) == keeper.end()) {
      node_q.push(f); node_q.push(r);
      breadth_q.push(breadth + 1);
    }
    
    r = ((kmer_r << 2) & bitmask) | twobit_comp('G');
    f = kmer_f >> 2 | (twobit_repr('G') << rc_left_shift);
    if (get_count(uniqify_rc(f,r)) && 
	keeper.find(uniqify_rc(f,r)) == keeper.end()) {
      node_q.push(f); node_q.push(r);
      breadth_q.push(breadth + 1);
    }

    r = ((kmer_r << 2) & bitmask) | twobit_comp('T');
    f = kmer_f >> 2 | (twobit_repr('T') << rc_left_shift);
    if (get_count(uniqify_rc(f,r)) && 
	keeper.find(uniqify_rc(f,r)) == keeper.end()) {
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
  f = ((kmer_f << 2) & bitmask) | twobit_repr('A');
  r = kmer_r >> 2 | (twobit_comp('A') << rc_left_shift);
  if (get_count(uniqify_rc(f,r)) && 
      seen->find(uniqify_rc(f,r)) == seen->end()) {
    count += count_kmers_within_depth(f, r, depth - 1, max_count - count,
				      seen);
    if (count >= max_count) { return count; }
  }

  f = ((kmer_f << 2) & bitmask) | twobit_repr('C');
  r = kmer_r >> 2 | (twobit_comp('C') << rc_left_shift);
  if (get_count(uniqify_rc(f,r)) && 
      seen->find(uniqify_rc(f,r)) == seen->end()) {
    count += count_kmers_within_depth(f, r, depth -1, max_count - count,
				      seen);
    if (count >= max_count) { return count; }
    ;
  }

  f = ((kmer_f << 2) & bitmask) | twobit_repr('G');
  r = kmer_r >> 2 | (twobit_comp('G') << rc_left_shift);
  if (get_count(uniqify_rc(f,r)) && 
      seen->find(uniqify_rc(f,r)) == seen->end()) {
    count += count_kmers_within_depth(f, r, depth -1, max_count - count,
				      seen);
    if (count >= max_count) { return count; }
    ;
  }

  f = ((kmer_f << 2) & bitmask) | twobit_repr('T');
  r = kmer_r >> 2 | (twobit_comp('T') << rc_left_shift);
  if (get_count(uniqify_rc(f,r)) && 
      seen->find(uniqify_rc(f,r)) == seen->end()) {
    count += count_kmers_within_depth(f, r, depth -1, max_count - count,
				      seen);
    if (count >= max_count) { return count; }
    ;
  }

  // PREVIOUS.
  r = ((kmer_r << 2) & bitmask) | twobit_comp('A');
  f = kmer_f >> 2 | (twobit_repr('A') << rc_left_shift);
  if (get_count(uniqify_rc(f,r)) && 
      seen->find(uniqify_rc(f,r)) == seen->end()) {
    count += count_kmers_within_depth(f, r, depth -1, max_count - count,
				      seen);
    if (count >= max_count) { return count; }
    ;
  }

  r = ((kmer_r << 2) & bitmask) | twobit_comp('C');
  f = kmer_f >> 2 | (twobit_repr('C') << rc_left_shift);
  if (get_count(uniqify_rc(f,r)) && 
      seen->find(uniqify_rc(f,r)) == seen->end()) {
    count += count_kmers_within_depth(f, r, depth -1, max_count - count,
				      seen);
    if (count >= max_count) { return count; }
    ;
  }
    
  r = ((kmer_r << 2) & bitmask) | twobit_comp('G');
  f = kmer_f >> 2 | (twobit_repr('G') << rc_left_shift);
  if (get_count(uniqify_rc(f,r)) && 
      seen->find(uniqify_rc(f,r)) == seen->end()) {
    count += count_kmers_within_depth(f, r, depth -1, max_count - count,
				      seen);
    if (count >= max_count) { return count; }
    ;
  }

  r = ((kmer_r << 2) & bitmask) | twobit_comp('T');
  f = kmer_f >> 2 | (twobit_repr('T') << rc_left_shift);
  if (get_count(uniqify_rc(f,r)) && 
      seen->find(uniqify_rc(f,r)) == seen->end()) {
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
    if (keeper.find(kmer) != keeper.end()) {
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
    f = ((kmer_f << 2) & bitmask) | twobit_repr('A');
    r = kmer_r >> 2 | (twobit_comp('A') << rc_left_shift);
    if (get_count(uniqify_rc(f,r)) && 
	keeper.find(uniqify_rc(f,r)) == keeper.end()) {
      node_q.push(f); node_q.push(r);
      breadth_q.push(breadth + 1);
    }

    f = ((kmer_f << 2) & bitmask) | twobit_repr('C');
    r = kmer_r >> 2 | (twobit_comp('C') << rc_left_shift);
    if (get_count(uniqify_rc(f,r)) && 
	keeper.find(uniqify_rc(f,r)) == keeper.end()) {
      node_q.push(f); node_q.push(r);
      breadth_q.push(breadth + 1);
    }

    f = ((kmer_f << 2) & bitmask) | twobit_repr('G');
    r = kmer_r >> 2 | (twobit_comp('G') << rc_left_shift);
    if (get_count(uniqify_rc(f,r)) && 
	keeper.find(uniqify_rc(f,r)) == keeper.end()) {
      node_q.push(f); node_q.push(r);
      breadth_q.push(breadth + 1);
    }

    f = ((kmer_f << 2) & bitmask) | twobit_repr('T');
    r = kmer_r >> 2 | (twobit_comp('T') << rc_left_shift);
    if (get_count(uniqify_rc(f,r)) && 
	keeper.find(uniqify_rc(f,r)) == keeper.end()) {
      node_q.push(f); node_q.push(r);
      breadth_q.push(breadth + 1);
    }

    // PREVIOUS.
    r = ((kmer_r << 2) & bitmask) | twobit_comp('A');
    f = kmer_f >> 2 | (twobit_repr('A') << rc_left_shift);
    if (get_count(uniqify_rc(f,r)) && 
	keeper.find(uniqify_rc(f,r)) == keeper.end()) {
      node_q.push(f); node_q.push(r);
      breadth_q.push(breadth + 1);
    }

    r = ((kmer_r << 2) & bitmask) | twobit_comp('C');
    f = kmer_f >> 2 | (twobit_repr('C') << rc_left_shift);
    if (get_count(uniqify_rc(f,r)) && 
	keeper.find(uniqify_rc(f,r)) == keeper.end()) {
      node_q.push(f); node_q.push(r);
      breadth_q.push(breadth + 1);
    }
    
    r = ((kmer_r << 2) & bitmask) | twobit_comp('G');
    f = kmer_f >> 2 | (twobit_repr('G') << rc_left_shift);
    if (get_count(uniqify_rc(f,r)) && 
	keeper.find(uniqify_rc(f,r)) == keeper.end()) {
      node_q.push(f); node_q.push(r);
      breadth_q.push(breadth + 1);
    }

    r = ((kmer_r << 2) & bitmask) | twobit_comp('T');
    f = kmer_f >> 2 | (twobit_repr('T') << rc_left_shift);
    if (get_count(uniqify_rc(f,r)) && 
	keeper.find(uniqify_rc(f,r)) == keeper.end()) {
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
    if (keeper.find(kmer) != keeper.end()) {
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
    f = ((kmer_f << 2) & bitmask) | twobit_repr('A');
    r = kmer_r >> 2 | (twobit_comp('A') << rc_left_shift);
    if (get_count(uniqify_rc(f,r)) && 
	keeper.find(uniqify_rc(f,r)) == keeper.end()) {
      node_q.push(f); node_q.push(r);
      breadth_q.push(breadth + 1);
    }

    f = ((kmer_f << 2) & bitmask) | twobit_repr('C');
    r = kmer_r >> 2 | (twobit_comp('C') << rc_left_shift);
    if (get_count(uniqify_rc(f,r)) && 
	keeper.find(uniqify_rc(f,r)) == keeper.end()) {
      node_q.push(f); node_q.push(r);
      breadth_q.push(breadth + 1);
    }

    f = ((kmer_f << 2) & bitmask) | twobit_repr('G');
    r = kmer_r >> 2 | (twobit_comp('G') << rc_left_shift);
    if (get_count(uniqify_rc(f,r)) && 
	keeper.find(uniqify_rc(f,r)) == keeper.end()) {
      node_q.push(f); node_q.push(r);
      breadth_q.push(breadth + 1);
    }

    f = ((kmer_f << 2) & bitmask) | twobit_repr('T');
    r = kmer_r >> 2 | (twobit_comp('T') << rc_left_shift);
    if (get_count(uniqify_rc(f,r)) && 
	keeper.find(uniqify_rc(f,r)) == keeper.end()) {
      node_q.push(f); node_q.push(r);
      breadth_q.push(breadth + 1);
    }

    // PREVIOUS.
    r = ((kmer_r << 2) & bitmask) | twobit_comp('A');
    f = kmer_f >> 2 | (twobit_repr('A') << rc_left_shift);
    if (get_count(uniqify_rc(f,r)) && 
	keeper.find(uniqify_rc(f,r)) == keeper.end()) {
      node_q.push(f); node_q.push(r);
      breadth_q.push(breadth + 1);
    }

    r = ((kmer_r << 2) & bitmask) | twobit_comp('C');
    f = kmer_f >> 2 | (twobit_repr('C') << rc_left_shift);
    if (get_count(uniqify_rc(f,r)) && 
	keeper.find(uniqify_rc(f,r)) == keeper.end()) {
      node_q.push(f); node_q.push(r);
      breadth_q.push(breadth + 1);
    }
    
    r = ((kmer_r << 2) & bitmask) | twobit_comp('G');
    f = kmer_f >> 2 | (twobit_repr('G') << rc_left_shift);
    if (get_count(uniqify_rc(f,r)) && 
	keeper.find(uniqify_rc(f,r)) == keeper.end()) {
      node_q.push(f); node_q.push(r);
      breadth_q.push(breadth + 1);
    }

    r = ((kmer_r << 2) & bitmask) | twobit_comp('T');
    f = kmer_f >> 2 | (twobit_repr('T') << rc_left_shift);
    if (get_count(uniqify_rc(f,r)) && 
	keeper.find(uniqify_rc(f,r)) == keeper.end()) {
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
  const char * first_kmer = seq.c_str();
  HashIntoType kmer_f = 0, kmer_r = 0;
  _hash(first_kmer, _ksize, kmer_f, kmer_r);

  if (kmer_degree(kmer_f, kmer_r) > max_degree) {
    return _ksize;
  }

  for (unsigned int i = _ksize; i < seq.length(); i++) {
    _next_hash(seq[i], kmer_f, kmer_r);

    if (kmer_degree(kmer_f, kmer_r) > max_degree) {
      return i;
    }
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

  const char * first_kmer = seq.c_str();
  SeenSet path;

  HashIntoType kmer_f = 0, kmer_r = 0;
  HashIntoType kmer;
  SeenSet seen;

#if 0
  kmer = _hash(first_kmer, _ksize, kmer_f, kmer_r);
  path.insert(kmer);

  for (unsigned int i = _ksize; i < seq.length(); i++) {
    kmer = _next_hash(seq[i], kmer_f, kmer_r);
    path.insert(kmer);
  }
#endif // 0

  kmer = _hash(first_kmer, _ksize, kmer_f, kmer_r);
  count = count_kmers_within_depth(kmer_f, kmer_r, radius, max_volume, &seen);
  if (count >= max_volume) {
    return 0;
  }
  
  for (unsigned int i = _ksize; i < seq.length(); i++) {
    SeenSet seen;
    kmer = _next_hash(seq[i], kmer_f, kmer_r);
    count = count_kmers_within_depth(kmer_f, kmer_r, radius, max_volume,&seen);
    if (count >= max_volume) {
      return i - 1;
    }
  }

  return seq.length();
}

unsigned int Hashbits::trim_on_stoptags(std::string seq) const
{
  if (!check_read(seq)) {
    return 0;
  }

  const char * first_kmer = seq.c_str();
  SeenSet path;

  HashIntoType kmer_f = 0, kmer_r = 0;
  HashIntoType kmer;

  kmer = _hash(first_kmer, _ksize, kmer_f, kmer_r);

  if (stop_tags.find(kmer) != stop_tags.end()) {
    return 0;
  }

  for (unsigned int i = _ksize; i < seq.length(); i++) {
    kmer = _next_hash(seq[i], kmer_f, kmer_r);
    if (stop_tags.find(kmer) != stop_tags.end()) {
      return i - 1;
    }
  }
  return seq.length();
}

#if 0

void Hashbits::load_stop_tags(std::string filename)
{
  std::ifstream infile(filename.c_str());
  assert(infile.is_open());

  std::string line;
  HashIntoType kmer;

  getline(infile, line);
  while (!infile.eof()) {
    if (line[_ksize] != ' ') {
      std::cerr << "incorrect format\n";
      std::cerr << line << "\n";
      assert(0);
    }
    kmer = _hash(line.c_str(), _ksize);
    stop_tags.insert(kmer);

    getline(infile, line);
  }
  std::cout << "stop tags: " << stop_tags.size() << "\n";
}

#endif // 0

void Hashbits::traverse_from_tags(unsigned int distance,
				  unsigned int threshold,
				  unsigned int num_high_todo,
				  CountingHash &counting) const
{
  unsigned int i = 0;
  unsigned int n = 0;
  unsigned int count;
  unsigned int n_big = 0;

  std::cout << all_tags.size() << " tags...\n";
  SeenSet::const_iterator si = all_tags.end();
  si--;

  for (; si != all_tags.begin(); si--, i++) {
      n++;
      count = _traverse_from_tag(*si, distance, counting);

      if (count >= threshold) {
	n_big++;
      }
      if (n_big >= num_high_todo) {
	break;
      }
      std::cout << "traverse-counting from " << *si << " " << i
		<< " " << n << " " << n_big << "\n";
  }
  std::cout << "traversed from " << n << " tags total.\n";
}

unsigned int Hashbits::_traverse_from_tag(HashIntoType start,
					  unsigned int radius,
					  CountingHash &counting)
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
    if (keeper.find(kmer) != keeper.end()) {
      continue;
    }

    // keep track of seen kmers
    keeper.insert(kmer);
    total++;

    counting.count(kmer);

    assert(breadth >= cur_breadth); // keep track of watermark, for debugging.
    if (breadth > cur_breadth) { cur_breadth = breadth; }

    //
    // Enqueue next set of nodes.
    //

    // NEXT.
    f = next_f(kmer_f, 'A');
    r = next_r(kmer_r, 'A');
      
    // f = ((kmer_f << 2) & bitmask) | twobit_repr('A');
    // r = kmer_r >> 2 | (twobit_comp('A') << rc_left_shift);
    if (get_count(uniqify_rc(f,r)) && 
	keeper.find(uniqify_rc(f,r)) == keeper.end()) {
      node_q.push(f); node_q.push(r);
      breadth_q.push(breadth + 1);
    }

    f = ((kmer_f << 2) & bitmask) | twobit_repr('C');
    r = kmer_r >> 2 | (twobit_comp('C') << rc_left_shift);
    if (get_count(uniqify_rc(f,r)) && 
	keeper.find(uniqify_rc(f,r)) == keeper.end()) {
      node_q.push(f); node_q.push(r);
      breadth_q.push(breadth + 1);
    }

    f = ((kmer_f << 2) & bitmask) | twobit_repr('G');
    r = kmer_r >> 2 | (twobit_comp('G') << rc_left_shift);
    if (get_count(uniqify_rc(f,r)) && 
	keeper.find(uniqify_rc(f,r)) == keeper.end()) {
      node_q.push(f); node_q.push(r);
      breadth_q.push(breadth + 1);
    }

    f = ((kmer_f << 2) & bitmask) | twobit_repr('T');
    r = kmer_r >> 2 | (twobit_comp('T') << rc_left_shift);
    if (get_count(uniqify_rc(f,r)) && 
	keeper.find(uniqify_rc(f,r)) == keeper.end()) {
      node_q.push(f); node_q.push(r);
      breadth_q.push(breadth + 1);
    }

    // PREVIOUS.
    r = ((kmer_r << 2) & bitmask) | twobit_comp('A');
    f = kmer_f >> 2 | (twobit_repr('A') << rc_left_shift);
    if (get_count(uniqify_rc(f,r)) && 
	keeper.find(uniqify_rc(f,r)) == keeper.end()) {
      node_q.push(f); node_q.push(r);
      breadth_q.push(breadth + 1);
    }

    r = ((kmer_r << 2) & bitmask) | twobit_comp('C');
    f = kmer_f >> 2 | (twobit_repr('C') << rc_left_shift);
    if (get_count(uniqify_rc(f,r)) && 
	keeper.find(uniqify_rc(f,r)) == keeper.end()) {
      node_q.push(f); node_q.push(r);
      breadth_q.push(breadth + 1);
    }
    
    r = ((kmer_r << 2) & bitmask) | twobit_comp('G');
    f = kmer_f >> 2 | (twobit_repr('G') << rc_left_shift);
    if (get_count(uniqify_rc(f,r)) && 
	keeper.find(uniqify_rc(f,r)) == keeper.end()) {
      node_q.push(f); node_q.push(r);
      breadth_q.push(breadth + 1);
    }

    r = ((kmer_r << 2) & bitmask) | twobit_comp('T');
    f = kmer_f >> 2 | (twobit_repr('T') << rc_left_shift);
    if (get_count(uniqify_rc(f,r)) && 
	keeper.find(uniqify_rc(f,r)) == keeper.end()) {
      node_q.push(f); node_q.push(r);
      breadth_q.push(breadth + 1);
    }
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
	string kmer = seq.substr(i, i + _ksize);
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

  std::cout << "Inserted " << stop_tags.size() << " stop tags\n";
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

void Hashbits::load_stop_tags(std::string infilename)
{
  ifstream infile(infilename.c_str(), ios::binary);
  stop_tags.clear();

  unsigned int tagset_size = 0;
  infile.read((char *) &tagset_size, sizeof(tagset_size));

  HashIntoType * buf = new HashIntoType[tagset_size];

  infile.read((char *) buf, sizeof(HashIntoType) * tagset_size);

  for (unsigned int i = 0; i < tagset_size; i++) {
    stop_tags.insert(buf[i]);
  }
  
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
