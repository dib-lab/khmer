#include "hashbits.hh"
#include "subset.hh"
#include "parsers.hh"

#define IO_BUF_SIZE 1000*1000*1000

// #define VALIDATE_PARTITIONS

using namespace khmer;
using namespace std;

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

void SubsetPartition::count_partitions(unsigned int& n_partitions,
				       unsigned int& n_unassigned)
{
  n_partitions = 0;
  n_unassigned = 0;

  PartitionSet partitions;

  //
  // go through all the tagged kmers and count partitions/orphan.
  //

  for (PartitionMap::const_iterator pi = partition_map.begin();
       pi != partition_map.end(); pi++) {
    PartitionID * partition_p = pi->second;
    if (partition_p) {
      partitions.insert(*partition_p);
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
  HashIntoType kmer = 0;

  const unsigned int ksize = _ht->ksize();

  //
  // go through all the reads, and take those with assigned partitions
  // and output them.
  //

  while(!parser->is_complete()) {
    read = parser->get_next_read();
    seq = read.seq;

    if (_ht->check_read(seq)) {
      const char * kmer_s = seq.c_str();

      bool found_tag = false;
      for (unsigned int i = 0; i < seq.length() - ksize + 1; i++) {
	kmer = _hash(kmer_s + i, ksize);

	// is this a known tag?
	if (set_contains(partition_map, kmer)) {
	  found_tag = true;
	  break;
	}
      }

      // all sequences should have at least one tag in them.
      // assert(found_tag);  @CTB currently breaks tests.  give fn flag to
      // disable.

      PartitionID partition_id = 0;
      if (found_tag) {
	PartitionID * partition_p = partition_map[kmer];
	if (partition_p == NULL ){
	  partition_id = 0;
	  n_singletons++;
	} else {
	  partition_id = *partition_p;
	  partitions.insert(partition_id);
	}
      }

      if (partition_id > 0 || output_unassigned) {
	outfile << ">" << read.name << "\t" << partition_id;
	outfile << "\n" << seq << "\n";
      }
#ifdef VALIDATE_PARTITIONS
      std::cout << "checking: " << read.name << "\n";
      assert(is_single_partition(seq));
#endif // VALIDATE_PARTITIONS
	       
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

// find_all_tags: the core of the partitioning code.  finds all tagged k-mers
//    connected to kmer_f/kmer_r in the graph.

void SubsetPartition::find_all_tags(HashIntoType kmer_f,
				    HashIntoType kmer_r,
				    SeenSet& tagged_kmers,
				    const SeenSet& all_tags,
				    bool break_on_stop_tags)
{
  const HashIntoType bitmask = _ht->bitmask;

  HashIntoType f, r;
  bool first = true;
  NodeQueue node_q;
  std::queue<unsigned int> breadth_q;
  unsigned int cur_breadth = 0;
  unsigned int breadth = 0;
  const unsigned int max_breadth = (2 * _ht->_tag_density) + 1;

  const unsigned int rc_left_shift = _ht->ksize()*2 - 2;
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

    if (break_on_stop_tags && set_contains(_ht->stop_tags, kmer)) {
      continue;
    }

    // keep track of seen kmers
    keeper.insert(kmer);
    //    cout << "INSERT: " << _revhash(kmer, _ht->ksize()) << "=" << (int) (_ht->get_count(kmer)) << " xx " << kmer % _ht->n_entries() << " =\n";
    total++;

    // if this is a high-circumference k-mer, do nothing more with it;
    // definitely do not connect partitions across it!
#if 0
    if (_ht->count_kmers_on_radius(kmer_f, kmer_r, CIRCUM_RADIUS,
				   CIRCUM_MAX_VOL) > MAX_CIRCUM) {
      continue;
    }
#endif // 0

    // Is this a kmer-to-tag, and have we put this tag in a partition already?
    // Search no further in this direction.  (This is where we connect
    // partitions.)
    if (!first && set_contains(all_tags, kmer)) {
      tagged_kmers.insert(kmer);
      continue;
    }

    assert(breadth >= cur_breadth); // keep track of watermark, for debugging.
    if (breadth > cur_breadth) { cur_breadth = breadth; }

    if (breadth >= max_breadth) { continue; } // truncate search @CTB exit?

    //
    // Enqueue next set of nodes.
    //

    // NEXT
    f = next_f(kmer_f, 'A');
    r = next_r(kmer_r, 'A');
    if (_ht->get_count(uniqify_rc(f,r)) &&
	!set_contains(keeper, uniqify_rc(f,r))) {
      node_q.push(f); node_q.push(r);
      breadth_q.push(breadth + 1);
    }

    f = next_f(kmer_f, 'C');
    r = next_r(kmer_r, 'C');
    if (_ht->get_count(uniqify_rc(f,r)) && 
        !set_contains(keeper, uniqify_rc(f,r))) {
      node_q.push(f); node_q.push(r);
      breadth_q.push(breadth + 1);
    }

    f = next_f(kmer_f, 'G');
    r = next_r(kmer_r, 'G');
    if (_ht->get_count(uniqify_rc(f,r)) && 
        !set_contains(keeper, uniqify_rc(f,r))) {
      node_q.push(f); node_q.push(r);
      breadth_q.push(breadth + 1);
    }

    f = next_f(kmer_f, 'T');
    r = next_r(kmer_r, 'T');
    if (_ht->get_count(uniqify_rc(f,r)) && 
        !set_contains(keeper, uniqify_rc(f,r))) {
      node_q.push(f); node_q.push(r);
      breadth_q.push(breadth + 1);
    }

    // PREVIOUS.
    r = prev_r(kmer_r, 'A');
    f = prev_f(kmer_f, 'A');
    if (_ht->get_count(uniqify_rc(f,r)) && 
        !set_contains(keeper, uniqify_rc(f,r))) {
      node_q.push(f); node_q.push(r);
      breadth_q.push(breadth + 1);
    }

    r = prev_r(kmer_r, 'C');
    f = prev_f(kmer_f, 'C');
    if (_ht->get_count(uniqify_rc(f,r)) && 
        !set_contains(keeper, uniqify_rc(f,r))) {
      node_q.push(f); node_q.push(r);
      breadth_q.push(breadth + 1);
    }
    
    r = prev_r(kmer_r, 'G');
    f = prev_f(kmer_f, 'G');
    if (_ht->get_count(uniqify_rc(f,r)) && 
        !set_contains(keeper, uniqify_rc(f,r))) {
      node_q.push(f); node_q.push(r);
      breadth_q.push(breadth + 1);
    }

    r = prev_r(kmer_r, 'T');
    f = prev_f(kmer_f, 'T');
    if (_ht->get_count(uniqify_rc(f,r)) && 
        !set_contains(keeper, uniqify_rc(f,r))) {
      node_q.push(f); node_q.push(r);
      breadth_q.push(breadth + 1);
    }

    first = false;
  }
}

///////////////////////////////////////////////////////////////////////

void SubsetPartition::do_partition(HashIntoType first_kmer,
				   HashIntoType last_kmer,
				   bool break_on_stop_tags,
				   CallbackFn callback,
				   void * callback_data)
{
  unsigned int total_reads = 0;

  std::string kmer_s;
  HashIntoType kmer_f, kmer_r, kmer;
  SeenSet tagged_kmers;
  const unsigned char ksize = _ht->ksize();

  SeenSet::const_iterator si, end;

  if (first_kmer) {
    si = _ht->all_tags.find(first_kmer);
  } else {
    si = _ht->all_tags.begin();
  }
  if (last_kmer) {
    end = _ht->all_tags.find(last_kmer);
  } else {
    end = _ht->all_tags.end();
  }

  for (; si != end; si++) {
    total_reads++;

    kmer_s = _revhash(*si, ksize); // @CTB hackity hack hack!
    kmer = _hash(kmer_s.c_str(), ksize, kmer_f, kmer_r);

    // find all tagged kmers within range.
    tagged_kmers.clear();
    find_all_tags(kmer_f, kmer_r, tagged_kmers, _ht->all_tags,
		  break_on_stop_tags);

    // assign the partition ID
    assign_partition_id(kmer, tagged_kmers);

    // run callback, if specified
    if (total_reads % CALLBACK_PERIOD == 0 && callback) {
      cout << "...subset-part " << first_kmer << "-" << last_kmer << ": " << total_reads << " <- " << next_partition_id << "\n";
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

void SubsetPartition::set_partition_id(std::string kmer_s, PartitionID p)
{
  HashIntoType kmer;
  assert(kmer_s.length() >= _ht->ksize());
  kmer = _hash(kmer_s.c_str(), _ht->ksize());

  set_partition_id(kmer, p);
}

void SubsetPartition::set_partition_id(HashIntoType kmer, PartitionID p)
{
  PartitionPtrSet * s = reverse_pmap[p];
  PartitionID * pp = NULL;
  if (s == NULL) {
    s = new PartitionPtrSet();
    pp = new unsigned int(p);
    s->insert(pp);
    reverse_pmap[p] = s;
  } else {
    pp = *(s->begin());
  }
  partition_map[kmer] = pp;

  if (next_partition_id <= p) {
    next_partition_id = p + 1;
  }
}

PartitionID SubsetPartition::assign_partition_id(HashIntoType kmer,
						 SeenSet& tagged_kmers)

{
  PartitionID return_val = 0; 
  PartitionID * pp = NULL;

  // did we find a tagged kmer?
  if (tagged_kmers.size() >= 1) {
    pp = _join_partitions_by_tags(tagged_kmers, kmer);
    return_val = *pp;
  } else {
    PartitionMap::iterator pi = partition_map.find(kmer);
    if (pi != partition_map.end()) {
      assert(pi->second == NULL); // if it's not... reverse_pmap removal TBD.
    }

    partition_map.erase(kmer);
    return_val = 0;
  }

  return return_val;
}

// _join_partitions_by_tags combines the tags in 'tagged_kmers' into a single
// partition, creating or reassigning partitions as necessary.  Low level
// function!

PartitionID * SubsetPartition::_join_partitions_by_tags(
                   const SeenSet& tagged_kmers,
		   const HashIntoType kmer)
{
  SeenSet::const_iterator it = tagged_kmers.begin();
  unsigned int * this_partition_p = NULL;

  // find first assigned partition ID in tagged set
  while (it != tagged_kmers.end()) {
    this_partition_p = partition_map[*it];
    if (this_partition_p != NULL) {
      break;
    }
    it++;
  }

  // no partition ID? allocate new!
  if (this_partition_p == NULL) {
    this_partition_p = new PartitionID(next_partition_id);
    next_partition_id++;

    PartitionPtrSet * s = new PartitionPtrSet();
    s->insert(this_partition_p);
    reverse_pmap[*this_partition_p] = s;
  }
  
  // reassign all partitions individually.
  it = tagged_kmers.begin();
  for (; it != tagged_kmers.end(); ++it) {
    PartitionMap::iterator pi = partition_map.find(*it);

    if (pi == partition_map.end()) { // no entry? insert.
      partition_map[*it] = this_partition_p;
    } else {
      PartitionID * pp_id = pi->second;

      if (pp_id == NULL) {	// entry is null? set;
	pi->second = this_partition_p;
      } else if (*pp_id != *this_partition_p) { // != entry? join partitions.
	_merge_two_partitions(this_partition_p, pp_id);
      }
    }
  }

  assert(this_partition_p != NULL);
  partition_map[kmer] = this_partition_p;

  return this_partition_p;
}

// _merge_two_partitions merges the 'merge_pp' partition into the
// 'the_pp' partition.  It does this by joining the reverse pointer
// map structures for two partitions and resetting each partition
// pointer individually.

PartitionID * SubsetPartition::_merge_two_partitions(PartitionID *the_pp,
						     PartitionID *merge_pp)
{
  PartitionPtrSet * s = reverse_pmap[*the_pp];
  PartitionPtrSet * t = reverse_pmap[*merge_pp];

  // Choose the smaller of two sets to loop over.
  if (s->size() < t->size()) {
    PartitionPtrSet * tmp = s;  s = t; t = tmp;
    PartitionID * tmp2 = the_pp; the_pp = merge_pp; merge_pp = tmp2;
  }

  // Get rid of the reverse pointer for the old partition.
  reverse_pmap.erase(*merge_pp);

  // Merge all of the elements in the to-be-replaced PartitionPtrSet
  // into the merged partition.
  for (PartitionPtrSet::iterator pi = t->begin(); pi != t->end(); pi++) {
    PartitionID * iter_pp;
    iter_pp = *pi;

    *iter_pp = *the_pp;	// reset the partition ID to the new one.
    s->insert(iter_pp);
  }
  delete t;

  return the_pp;
}

PartitionID SubsetPartition::join_partitions(PartitionID orig, PartitionID join)
{
  if (orig == join) { return orig; }
  if (orig == 0 || join == 0) { return 0; }

  if (reverse_pmap.find(orig) == reverse_pmap.end() ||
      reverse_pmap.find(join) == reverse_pmap.end() ||
      reverse_pmap[orig] == NULL ||
      reverse_pmap[join] == NULL) {
    return 0;
  }

  PartitionID * orig_pp = *(reverse_pmap[orig]->begin());
  PartitionID * join_pp = *(reverse_pmap[join]->begin());

  _merge_two_partitions(orig_pp, join_pp);

  return orig;
}

PartitionID SubsetPartition::get_partition_id(std::string kmer_s)
{
  HashIntoType kmer;
  assert(kmer_s.length() >= _ht->ksize());
  kmer = _hash(kmer_s.c_str(), _ht->ksize());

  return get_partition_id(kmer);
}

PartitionID SubsetPartition::get_partition_id(HashIntoType kmer)
{
  if (partition_map.find(kmer) != partition_map.end()) {
    PartitionID * pp = partition_map[kmer];
    if (pp == NULL) {
      return 0;
    }
    return *pp;
  }
  return 0;
}

void SubsetPartition::merge(SubsetPartition * other)
{
  if (this == other) { return; }

  PartitionPtrMap other_to_this;

  PartitionMap::const_iterator pi = other->partition_map.begin();
  for (; pi != other->partition_map.end(); pi++) {
    if (pi->second) {
      _merge_other(pi->first, *(pi->second), other_to_this);
    }
  }
}


// @CTB dead code?
void SubsetPartition::_merge_from_disk_consolidate(PartitionPtrMap& diskp_to_pp)
{
  for (PartitionPtrMap::iterator pp = diskp_to_pp.begin();
       pp != diskp_to_pp.end(); pp++) {
    PartitionPtrSet * s = reverse_pmap[*(pp->second)];
    if (s->size() > 1) {
      pp->second = *(s->begin());
    }
  }

  for (PartitionMap::iterator pi = partition_map.begin();
       pi != partition_map.end(); pi++) {
    if (pi->second) {
      PartitionPtrSet * s = reverse_pmap[*(pi->second)];
      if (s->size() > 1) {
	pi->second = *(s->begin());
      }
    }
  }

  return;			// @@CTB??

  for (ReversePartitionMap::iterator ri = reverse_pmap.begin();
       ri != reverse_pmap.end(); ri++) {
    PartitionPtrSet * s = ri->second;
    if (s->size() > 1) {
      PartitionPtrSet::iterator si = s->begin();
      si++;

      while(si != s->end()) { free(*si); si++; }
      
      si = s->begin();
      si++;
      s->erase(si, s->end());
    }
  }
}

void SubsetPartition::_merge_other(HashIntoType tag,
				   PartitionID other_partition,
				   PartitionPtrMap& diskp_to_pp)
{
  if (set_contains(_ht->stop_tags, tag)) { // don't merge if it's a stop_tag
    return;
  }

  // OK.  Does our current partitionmap have this?
  PartitionID * pp_0;
  pp_0 = partition_map[tag];

  if (pp_0 == NULL) {	// No!  OK, map to new 'un.
    PartitionID * existing_pp_0 = diskp_to_pp[other_partition];

    if (existing_pp_0) {	// already seen this other_partition
      partition_map[tag] = existing_pp_0;
    }
    else {			// new other_partition! create a new partition.
      pp_0 = get_new_partition();

      PartitionPtrSet * pp_set = new PartitionPtrSet();
      pp_set->insert(pp_0);
      reverse_pmap[*pp_0] = pp_set;
      partition_map[tag] = pp_0;

      diskp_to_pp[other_partition] = pp_0;
    }
  }
  else {			// yes, we've seen this tag before...
    PartitionID * existing_pp_0 = diskp_to_pp[other_partition];

    if (existing_pp_0) {	// mapping exists.  copacetic?
      if (*pp_0 == *existing_pp_0) {
	;			// yep! nothing to do, yay!
      } else {
	// remapping must be done... we need to merge!
	// the two partitions to merge are *pp_0 and *existing_pp_0.
	// we also need to reset existing_pp_0 in diskp_to_pp to pp_0.

	pp_0 = _merge_two_partitions(pp_0, existing_pp_0);
	diskp_to_pp[other_partition] = pp_0;
      }
    }
    else {
      // no, does not exist in our mapping yet.  but that's ok,
      // we can fix that.
      diskp_to_pp[other_partition] = pp_0;
    }
  }
}

void SubsetPartition::merge_from_disk(string other_filename)
{
  ifstream infile(other_filename.c_str(), ios::binary);
  char * buf = NULL;
  buf = new char[IO_BUF_SIZE];

  unsigned int n_bytes = 0;
  unsigned int loaded = 0;
  unsigned int remainder;

  assert(infile.is_open());

  PartitionPtrMap diskp_to_pp;

  HashIntoType * kmer_p = NULL;
  PartitionID * diskp = NULL;

  //
  // Run through the entire partitionmap file, figuring out what partition IDs
  // are present.
  //

  remainder = 0;
  unsigned int iteration = 0;
  while (!infile.eof()) {
    unsigned int i;

    infile.read(buf + remainder, IO_BUF_SIZE - remainder);
    n_bytes = infile.gcount() + remainder;
    remainder = n_bytes % (sizeof(PartitionID) + sizeof(HashIntoType));
    n_bytes -= remainder;

    iteration++;

    for (i = 0; i < n_bytes;) {
      kmer_p = (HashIntoType *) (buf + i);
      i += sizeof(HashIntoType);
      diskp = (PartitionID *) (buf + i);
      i += sizeof(PartitionID);

      assert(*diskp != 0);		// sanity check.

      _merge_other(*kmer_p, *diskp, diskp_to_pp);

      loaded++;
    }
    assert(i == n_bytes);
    memcpy(buf, buf + n_bytes, remainder);

    // _merge_from_disk_consolidate(diskp_to_pp);
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

  // For each tag in the partition map, save the tag and the associated
  // partition ID.

  PartitionMap::const_iterator pi = partition_map.begin();
  for (; pi != partition_map.end(); pi++) {
    PartitionID p_id;

    HashIntoType kmer = pi->first;
    if (pi->second != NULL) {	// if a partition ID has been assigned... save.
      p_id = *(pi->second);

      // each record consists of one tag followed by one PartitionID.
      kmer_p = (HashIntoType *) (buf + n_bytes);
      *kmer_p = kmer;
      n_bytes += sizeof(HashIntoType);

      pp = (PartitionID *) (buf + n_bytes);
      *pp = p_id;
      n_bytes += sizeof(PartitionID);

      // flush to disk
      if (n_bytes >= IO_BUF_SIZE - sizeof(HashIntoType) - sizeof(PartitionID)) {
	outfile.write(buf, n_bytes);
	n_bytes = 0;
      }
    }
  }
  // save remainder.
  if (n_bytes) {
    outfile.write(buf, n_bytes);
  }
  outfile.close();

  delete buf;
}
					 
void SubsetPartition::load_partitionmap(string infilename)
{
  // @CTB make sure this is an empty partition...
  merge_from_disk(infilename);
}


void SubsetPartition::_validate_pmap()
{
  for (PartitionMap::const_iterator pi = partition_map.begin();
       pi != partition_map.end(); pi++) {
    //HashIntoType kmer = (*pi).first;
    PartitionID * pp_id = (*pi).second;

    if (pp_id != NULL) {
      assert(*pp_id >= 1);
      assert(*pp_id < next_partition_id);
    }
  }

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

void SubsetPartition::_clear_all_partitions()
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


bool SubsetPartition::is_single_partition(std::string seq)
{
  if (!_ht->check_read(seq)) {
    return 0;
  }

  const char * first_kmer = seq.c_str();

  HashIntoType kmer_f = 0, kmer_r = 0;
  HashIntoType kmer;

  kmer = _hash(first_kmer, _ht->ksize(), kmer_f, kmer_r);

  PartitionSet partitions;
  PartitionID *pp;

  if (partition_map.find(kmer) != partition_map.end()) {
    pp = partition_map[kmer];
    if (pp) {
      partitions.insert(*pp);
    }
  }

  for (unsigned int i = _ht->ksize(); i < seq.length(); i++) {
    kmer = _ht->_next_hash(seq[i], kmer_f, kmer_r);
    if (partition_map.find(kmer) != partition_map.end()) {
      pp = partition_map[kmer];
      if (pp) {
	partitions.insert(*pp);
      }
    }
  }

  if (partitions.size() > 1) { return false; }

  return true;
} 

void SubsetPartition::join_partitions_by_path(std::string seq)
{
  SeenSet tagged_kmers;
  HashIntoType kmer_f, kmer_r, kmer;
  const unsigned int ksize = _ht->ksize();

  kmer = _hash(seq.c_str(), ksize, kmer_f, kmer_r);
  if (_ht->all_tags.find(kmer) != _ht->all_tags.end()) {
    tagged_kmers.insert(kmer);
  }

  for (unsigned int i = ksize; i < seq.length(); i++) {
    kmer = _ht->_next_hash(seq[i], kmer_f, kmer_r);
    if (_ht->all_tags.find(kmer) != _ht->all_tags.end()) {
      tagged_kmers.insert(kmer);
    }
  }

  // assert(tagged_kmers.size());
  assign_partition_id(*(tagged_kmers.begin()), tagged_kmers);
}

void SubsetPartition::partition_size_distribution(PartitionCountDistribution &d,
						  unsigned int& n_unassigned)
const
{
  PartitionCountMap cm;
  n_unassigned = 0;

  for (PartitionMap::const_iterator pi = partition_map.begin();
       pi != partition_map.end(); pi++) {
    if (pi->second) {
      cm[*(pi->second)]++;
    } else {
      n_unassigned++;
    }
  }

  for (PartitionCountMap::const_iterator cmi = cm.begin(); cmi != cm.end();
       cmi++) {
    d[cmi->second]++;
  }
}

unsigned int SubsetPartition::repartition_largest_partition(unsigned int distance,
						    unsigned int threshold,
						    unsigned int frequency,
						    CountingHash &counting)
{
  PartitionCountMap cm;
  unsigned int n_unassigned = 0;
  PartitionID biggest_p = 0;
  unsigned int next_largest = 0;

  for (PartitionMap::const_iterator pi = partition_map.begin();
       pi != partition_map.end(); pi++) {
    if (pi->second) {
      cm[*(pi->second)]++;
    } else {
      n_unassigned++;
    }
  }

  PartitionCountDistribution d;

  for (PartitionCountMap::const_iterator cmi = cm.begin(); cmi != cm.end();
       cmi++) {
    d[cmi->second]++;
  }

  PartitionCountDistribution::const_iterator di = d.end();
  di--;

  assert(d.size());

  for (PartitionCountMap::const_iterator cmi = cm.begin(); cmi != cm.end();
       cmi++) {
    if (cmi->second == di->first) {
      biggest_p = cmi->first;	// find PID of largest partition
    }
  }
  assert(biggest_p != 0);

  std::cout << "biggest partition: " << di->first << "\n";
  di--;
  std::cout << "biggest partition ID: " << biggest_p << "\n";

  next_largest = di->first;
  std::cout << "next biggest partition: " << di->first << "\n";

  ///

  SeenSet bigtags;
  _clear_partition(biggest_p, bigtags);
  std::cout << "gathered/cleared " << bigtags.size() << " tags.\n";

  /// 

  unsigned int i = 0;
  unsigned int n = 0;
  unsigned int count;
  unsigned int n_big = 0;
  SeenSet keeper;

  SeenSet::const_iterator si = bigtags.begin();

  for (; si != bigtags.end(); si++, i++) {
    n++;
    count = _ht->traverse_from_kmer(*si, distance, keeper);

    if (count >= threshold) {
      n_big++;
	
      SeenSet::const_iterator ti;
      for (ti = keeper.begin(); ti != keeper.end(); ti++) {
	if (counting.get_count(*ti) > frequency) {
	  _ht->stop_tags.insert(*ti);
	} else {
	  counting.count(*ti);
	}
      }
      std::cout << "traversed from " << n << " tags total, of " 
	        << bigtags.size() << "; "
		<< n_big << " big; size is " << keeper.size() << "\n";
    }
    keeper.clear();

    if (n % 1000 == 0) {
      std::cout << "found big 'un!  traversed " << n << " tags, " << n_big << " big; " <<
	bigtags.size() << " total tags; " << _ht->stop_tags.size() 
		<< " stop tags\n";
    }
  }

  // return next_largest;
  std::cout << "repartitioning...\n";
  repartition_a_partition(bigtags);

  // 

  cm.clear();
  for (PartitionMap::const_iterator pi = partition_map.begin();
       pi != partition_map.end(); pi++) {
    if (pi->second) {
      cm[*(pi->second)]++;
    } else {
      n_unassigned++;
    }
  }

  d.clear();
  for (PartitionCountMap::const_iterator cmi = cm.begin(); cmi != cm.end();
       cmi++) {
    d[cmi->second]++;
  }

  di = d.end();
  di--;

  std::cout << "post-repart new biggest partition size: " << di->first << "\n";

  return next_largest;
}

void SubsetPartition::repartition_a_partition(const SeenSet& partition_tags)
{
  SeenSet tagged_kmers;
  std::string kmer_s;
  HashIntoType kmer_f, kmer_r, kmer;
  unsigned int ksize = _ht->ksize();

  SeenSet::const_iterator si;

  unsigned n = 0;
  for (si = partition_tags.begin(); si != partition_tags.end(); si++, n++) {
    if (n % 1000 == 0) {
      std::cout << "repartitioning... on " << n << " of " << partition_tags.size() << "\n";
    }

    kmer_s = _revhash(*si, ksize); // @CTB hackity hack hack!
    kmer = _hash(kmer_s.c_str(), ksize, kmer_f, kmer_r);

    tagged_kmers.clear();
    find_all_tags(kmer_f, kmer_r, tagged_kmers, _ht->all_tags, true);

    // only join things already in bigtags.
    for (SeenSet::iterator ssi = tagged_kmers.begin();
	 ssi != tagged_kmers.end(); ssi++) {
      if (!set_contains(partition_tags, *ssi)) {
	tagged_kmers.erase(ssi);
      }
    }
    // std::cout << "joining: " << tagged_kmers.size() << "\n";
    assign_partition_id(kmer, tagged_kmers);
  }
}

// _clear_partition: given a partition ID, identifies all tags that belong
//    to that partition & (a) clears their PID, and (b) adds them to
//    the SeenSet partition_tags.  partition_tags is cleared first.

void SubsetPartition::_clear_partition(PartitionID the_partition,
				       SeenSet& partition_tags)
{
  partition_tags.clear();

  for (PartitionMap::iterator pi = partition_map.begin();
       pi != partition_map.end(); pi++) {
    if (pi->second && *(pi->second) == the_partition) {
      partition_tags.insert(pi->first);
    }
  }

  for (SeenSet::const_iterator si = partition_tags.begin();
       si != partition_tags.end(); si++) {
    partition_map.erase(*si);
  }

  // clear out the reverse partition mapping, too.
  PartitionPtrSet * ps = reverse_pmap[the_partition];
  for (PartitionPtrSet::iterator psi = ps->begin(); psi != ps->end(); psi++) {
    delete *psi;
  }
  delete ps;

  reverse_pmap.erase(the_partition);
}
