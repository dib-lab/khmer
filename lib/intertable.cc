//

#include "ktable.hh"
#include "intertable.hh"
#include "parsers.hh"

#define CALLBACK_PERIOD 10000

using namespace std;
using namespace khmer;

bool IntersectTable::check_read(const std::string &read)
{
   unsigned int i;
   bool is_valid = true;

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

   return is_valid;
}

void IntersectTable::do_partition(const std::string infile, const std::string outfile)
{
  unsigned int total_reads = 0;

  IParser* parser = IParser::get_parser(infile);
  Read read;
  string seq;
  bool is_valid;

  cout << "KSIZE: " << (int) _ksize << "\n";

  while(!parser->is_complete()) {
    // increment read number
    read = parser->get_next_read();
    total_reads++;

    seq = read.seq;

    is_valid = check_read(seq);
    if (is_valid) {
      const char * kmer_s = seq.c_str();
      HashIntoType kmer;

      IntersectionSet crossed;
      IntersectionID iid;

      for (unsigned int i = 0; i < seq.length() - _ksize + 1; i++) {
	kmer = _hash(kmer_s + i, _ksize);
	iid = get(kmer);
	if (iid) {
	  crossed.insert(partitions[iid]);
	}
      }

      if (crossed.size() == 0) {
	iid = get_next_iid();

	IntersectionSet * s = new IntersectionSet();
	s->insert(iid);
	revmap[iid] = s;

	assert(partitions[iid] == iid); // pid == iid
      } else if (crossed.size() >= 1) {
	iid = *(crossed.begin()); // no merging or anything; just take one.
      }

      for (unsigned int i = 0; i < seq.length() - _ksize + 1; i++) {
	kmer = _hash(kmer_s + i, _ksize);
	set(kmer, iid);
      }
    }

    // run callback, if specified
    if (total_reads % CALLBACK_PERIOD == 0) {
      try {
	cout << "..." << total_reads << "\n";
      } catch (...) {
	delete parser;
	throw;
      }
    }
  }

  delete parser;

  total_reads = 0;

  parser = IParser::get_parser(infile);

  while(!parser->is_complete()) {
    // increment read number
    read = parser->get_next_read();
    total_reads++;

    seq = read.seq;

    is_valid = check_read(seq);
    if (is_valid) {
      const char * kmer_s = seq.c_str();
      HashIntoType kmer;

      IntersectionSet crossed;
      IntersectionID iid;

      for (unsigned int i = 0; i < seq.length() - _ksize + 1; i++) {
	kmer = _hash(kmer_s + i, _ksize);
	iid = get(kmer);
	crossed.insert(partitions[iid]);
	assert(partitions[iid] <= iid); // pid <= iid
      }

      if (crossed.size() > 1) {
	IntersectionSet::const_iterator is = crossed.begin();

#if 0
	cout << "overlap: ";
	for (; is != crossed.end(); is++) {
	  cout << *is << " ";
	}
	cout << "\n";
#endif

	is = crossed.begin();
	IntersectionID first_pid = *(is); // will be LOWEST.
	// cout << "starting at: " << first_pid;

	IntersectionSet * s = revmap[first_pid];
	assert (s != NULL);

	IntersectionSet * t;
	is++;
	for (; is != crossed.end(); is++) {
	  t = revmap[*is];
	  for (IntersectionSet::const_iterator ti = t->begin(); ti != t->end();
	       ti++) {

	    s->insert(*ti);
	    partitions[*ti] = first_pid;
	    // cout << "ASSIGN " << *ti << " " << first_pid << "\n";
	  }
	  revmap.erase(*is);
	  delete t;
	  // cout << "merged: " << *is << " into " << first_pid << "\n";
	}
      }
    }

    // run callback, if specified
    if (total_reads % CALLBACK_PERIOD == 0) {
      try {
	cout << "..." << total_reads << " x2\n";
      } catch (...) {
	delete parser;
	throw;
      }
    }
  }

  delete parser;

  cout << "remapping x 2...\n";
  remap();
  cout << "done\n";

  total_reads = 0;

  parser = IParser::get_parser(infile);
  ofstream outfp(outfile.c_str());

  while(!parser->is_complete()) {
    // increment read number
    read = parser->get_next_read();
    total_reads++;

    seq = read.seq;

    is_valid = check_read(seq);
    if (is_valid) {
      const char * kmer_s = seq.c_str();
      HashIntoType kmer;

      IntersectionSet crossed;
      IntersectionID iid;

      for (unsigned int i = 0; i < seq.length() - _ksize + 1; i++) {
	kmer = _hash(kmer_s + i, _ksize);
	iid = get(kmer);
	crossed.insert(partitions[iid]);
      }

      assert(crossed.size() == 1);

      outfp << ">" << read.name << "\t" << *(crossed.begin()) << "\n"
	      << seq << "\n";
    }

    // run callback, if specified
    if (total_reads % CALLBACK_PERIOD == 0) {
      try {
	cout << "..." << total_reads << " x3\n";
      } catch (...) {
	delete parser;
	throw;
      }
    }
  }

  delete parser;
}

void IntersectTable::remap()
{
  IntersectionSet seen;

  for (HashIntoType i = 0; i < _tablesize; i++) {
    IntersectionID iid = _table[i];
    if (iid) {
      _table[i] = partitions[iid];
      
      seen.insert(partitions[_table[i]]);
    }
  }

  for (ReverseIntersectionMap::iterator ri = revmap.begin();
       ri != revmap.end(); ri++) {
    if (ri->second->size() > 1) {
      delete ri->second;
    
      IntersectionSet * s = new IntersectionSet();
      s->insert(ri->first);
      ri->second = s;
    }
  }

  cout << "total intersection IDs seen: " << seen.size() << "\n";

  cout << "revmap: " << revmap.size() << "\n";
  //assert(seen.size() == revmap.size());
}
