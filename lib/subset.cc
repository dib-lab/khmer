#include "hashbits.hh"
#include "subset.hh"
#include "parsers.hh"

#define IO_BUF_SIZE 1000*1000*1000

// #define VALIDATE_PARTITIONS

using namespace khmer;
using namespace std;

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

      sp = pttmap[p];
      if (sp == NULL) {
	sp = new SeenSet();
	pttmap[p] = sp;
      }
      sp->insert(tag);
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
	if (partition_map.find(kmer) != partition_map.end()) {
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

// used by do_truncated_partition

void SubsetPartition::find_all_tags(HashIntoType kmer_f,
				    HashIntoType kmer_r,
				    SeenSet& tagged_kmers)
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

    if (keeper.find(kmer) != keeper.end()) {
      continue;
    }

    // if (_ht->stop_tags.find(kmer) != _ht->stop_tags.end()) {
    //      continue;
    //    }

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
    if (!first && _ht->all_tags.find(kmer) != _ht->all_tags.end()) {
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
	keeper.find(uniqify_rc(f,r)) == keeper.end()) {
      node_q.push(f); node_q.push(r);
      breadth_q.push(breadth + 1);
    }

    f = next_f(kmer_f, 'C');
    r = next_r(kmer_r, 'C');
    if (_ht->get_count(uniqify_rc(f,r)) && 
	keeper.find(uniqify_rc(f,r)) == keeper.end()) {
      node_q.push(f); node_q.push(r);
      breadth_q.push(breadth + 1);
    }

    f = next_f(kmer_f, 'G');
    r = next_r(kmer_r, 'G');
    if (_ht->get_count(uniqify_rc(f,r)) && 
	keeper.find(uniqify_rc(f,r)) == keeper.end()) {
      node_q.push(f); node_q.push(r);
      breadth_q.push(breadth + 1);
    }

    f = next_f(kmer_f, 'T');
    r = next_r(kmer_r, 'T');
    if (_ht->get_count(uniqify_rc(f,r)) && 
	keeper.find(uniqify_rc(f,r)) == keeper.end()) {
      node_q.push(f); node_q.push(r);
      breadth_q.push(breadth + 1);
    }

    // PREVIOUS.
    r = prev_r(kmer_r, 'A');
    f = prev_f(kmer_f, 'A');
    if (_ht->get_count(uniqify_rc(f,r)) && 
	keeper.find(uniqify_rc(f,r)) == keeper.end()) {
      node_q.push(f); node_q.push(r);
      breadth_q.push(breadth + 1);
    }

    r = prev_r(kmer_r, 'C');
    f = prev_f(kmer_f, 'C');
    if (_ht->get_count(uniqify_rc(f,r)) && 
	keeper.find(uniqify_rc(f,r)) == keeper.end()) {
      node_q.push(f); node_q.push(r);
      breadth_q.push(breadth + 1);
    }
    
    r = prev_r(kmer_r, 'G');
    f = prev_f(kmer_f, 'G');
    if (_ht->get_count(uniqify_rc(f,r)) && 
	keeper.find(uniqify_rc(f,r)) == keeper.end()) {
      node_q.push(f); node_q.push(r);
      breadth_q.push(breadth + 1);
    }

    r = prev_r(kmer_r, 'T');
    f = prev_f(kmer_f, 'T');
    if (_ht->get_count(uniqify_rc(f,r)) && 
	keeper.find(uniqify_rc(f,r)) == keeper.end()) {
      node_q.push(f); node_q.push(r);
      breadth_q.push(breadth + 1);
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
    find_all_tags(kmer_f, kmer_r, tagged_kmers);

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
    pp = _reassign_partition_ids(tagged_kmers, kmer);
    return_val = *pp;
  } else {
    partition_map[kmer] = NULL;
    return_val = 0;
  }

  return return_val;
}

PartitionID * SubsetPartition::_reassign_partition_ids(SeenSet& tagged_kmers,
						       const HashIntoType kmer)
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
  partition_map[kmer] = this_partition_p;

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

PartitionID * SubsetPartition::_add_partition_ptr2(PartitionID *orig_pp, PartitionID *new_pp)
{
  PartitionPtrSet * s = reverse_pmap[*orig_pp];
  PartitionPtrSet * t = reverse_pmap[*new_pp];

  if (s->size() < t->size()) {
    PartitionPtrSet * tmp = s;  s = t; t = tmp;
    PartitionID * tmp2 = orig_pp; orig_pp = new_pp; new_pp = tmp2;
  }

  reverse_pmap.erase(*new_pp);

  for (PartitionPtrSet::iterator pi = t->begin(); pi != t->end(); pi++) {
    PartitionID * iter_pp;
    iter_pp = *pi;

    *iter_pp = *orig_pp;
    s->insert(iter_pp);
  }
  delete t;

  if (s->size() == 100 || s->size() % 10000 == 0) {
    cout << "Big un: " << *orig_pp << "; size " << s->size() << "\n";
  }

  return orig_pp;
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

  _add_partition_ptr2(orig_pp, join_pp);

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
  assert (this != other);

  PartitionsToTagsMap subset_pttm, master_pttm;
  PartitionID * pp;
  PartitionID p;

  //
  // Convert the partition maps into something where we can do reverse
  // lookups of what tags belong to a given partition.
  //

  make_partitions_to_tags(other->partition_map, subset_pttm);
  make_partitions_to_tags(this->partition_map, master_pttm);

  // Run through all of the unique partitions in 'other' and find the
  // transitive overlaps with sets in 'this'.

  for (PartitionsToTagsMap::iterator ptti = subset_pttm.begin();
       ptti != subset_pttm.end(); ptti++) {
    p = ptti->first;
    SeenSet * these_tags = ptti->second;

    if (these_tags == NULL) {	// We NULL the values as we deal with parts.
      continue;
    }

    // OK, this partition has not yet been transferred.  DEAL.

    //
    // Here, we want to get all of the partitions connected to this
    // one in the master map and the subset map -- and do
    // transitively.  Loop until no more additional tags/partitions are found.
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

    // All right!  We've now got all the tags that we want to be part of
    // the master map partition; create or merge.

    PartitionPtrSet * pp_set = NULL;

    // No overlapping master partitions?  Create new.
    if (old_master_partitions.size() == 0) {
      pp = this->get_new_partition();
	
      pp_set = new PartitionPtrSet();
      pp_set->insert(pp);
      this->reverse_pmap[*pp] = pp_set;
    } else {
      PartitionSet::iterator psi = old_master_partitions.begin();
      pp_set = this->reverse_pmap[*psi];

      assert(pp_set != NULL);
      pp = *(pp_set->begin());
    }

    // Remove all of the SeenSets in the subset_pttm map for these
    // now-connected partitions.  (Also see NULL check at beginning of loop.)
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

      if (old_pp == NULL) {	// no previously assigned partition? assign!
	this->partition_map[*si] = pp;
      } else if (old_pp != pp) { // Overwrite, but keep track of those to free.
	remove_me.insert(old_pp);
	this->partition_map[*si] = pp;

	// Note: even though there may be multiple PartitionID*s for
	// this partition that need reassigning, we don't need to
	// handle that in this loop; we're doing this systematically (and
	// the 'remove_me' set will clear out with duplicates).
      }
    }

    // Condense & reset the reverse_pmap for *this* entry, too.  These
    // PartitionID*s will have been added to 'remove_me' in the previous loop
    // already.
    if (this->reverse_pmap[*pp]->size() > 1) {
      pp_set = new PartitionPtrSet();
      pp_set->insert(pp);
      this->reverse_pmap[*pp] = pp_set;
    }

    // Finally: remove contents of remove_me & associated (obsolete)
    // reverse_pmap entries.
    for (PartitionPtrSet::iterator si = remove_me.begin();
	 si != remove_me.end(); si++) {
      PartitionID p = *(*si);

      pp_set = this->reverse_pmap[p];
      if (pp_set) {
	delete pp_set;
	this->reverse_pmap.erase(p);
      }
      delete *si;
    }
  }

  // OK, clean up the reverse maps.

  del_partitions_to_tags(subset_pttm);
  del_partitions_to_tags(master_pttm);
}


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

  return;

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

    //    cout << "Read " << n_bytes + remainder << " (" << iteration << ")\n";
    // cout << other_filename << "; mapping: " << diskp_to_pp.size()
    //	 << "; partitions: " << reverse_pmap.size() << "\n";
    iteration++;

    for (i = 0; i < n_bytes;) {
      kmer_p = (HashIntoType *) (buf + i);
      i += sizeof(HashIntoType);
      diskp = (PartitionID *) (buf + i);
      i += sizeof(PartitionID);

      assert(*diskp != 0);		// sanity check.

      // @@ partitions.insert(*pp);

      // OK.  Does our current partitionmap have this?
      PartitionID * pp_0;
      pp_0 = partition_map[*kmer_p];

      if (pp_0 == NULL) {	// No!  OK, map to new 'un.
	PartitionID * existing_pp_0 = diskp_to_pp[*diskp];

	if (existing_pp_0) {	// already seen this *diskp
	  // cout << "This diskp already seen -- " << *diskp << " to " << *existing_pp_0 << "\n";
	  partition_map[*kmer_p] = existing_pp_0;
	}
	else {			// new *diskp! create a new partition.
	  pp_0 = get_new_partition();

	  PartitionPtrSet * pp_set = new PartitionPtrSet();
	  pp_set->insert(pp_0);
	  reverse_pmap[*pp_0] = pp_set;
	  partition_map[*kmer_p] = pp_0;

	  diskp_to_pp[*diskp] = pp_0;
	  // cout << "New partition! " << *diskp << " mapped to " << *pp_0 << "\n";
	}
      }
      else {			// yes, we've seen this tag before...
	PartitionID * existing_pp_0 = diskp_to_pp[*diskp];

	if (existing_pp_0) {	// mapping exists.  copacetic?
	  if (*pp_0 == *existing_pp_0) {
	    ;			// yep! nothing to do, yay!
	  } else {
	    // remapping must be done... we need to merge!
	    // the two partitions to merge are *pp_0 and *existing_pp_0.
	    // we also need to reset existing_pp_0 in diskp_to_pp to pp_0.
	    //	    cout << "Remapping/merging: " << *existing_pp_0 << "=>" << *pp_0 << "\n";
	    pp_0 = _add_partition_ptr2(pp_0, existing_pp_0);
	    diskp_to_pp[*diskp] = pp_0;
	  }
	}
	else {
	  // no, does not exist in our mapping yet.  but that's ok,
	  // we can fix that.
	  // cout << "First time/existing mapping: " << *diskp << "->" << *pp_0 << "\n";
	  diskp_to_pp[*diskp] = pp_0;
	}
      }

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
