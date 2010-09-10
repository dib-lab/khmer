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
	  crossed.insert(iid);
	}
      }

      if (crossed.size() == 0) {
	cout << "no crossing.\n";
	iid = get_next_iid();
      } else {
	iid = find_min_iid(crossed);
	cout << "crossed! " << iid << "\n";

	for (IntersectionSet::const_iterator si = crossed.begin();
	     si != crossed.end(); si++ ){
	  if (*si != iid) {
	    intermap[*si] = iid;
	  }
	}
      }

      for (unsigned int i = 0; i < seq.length() - _ksize + 1; i++) {
	kmer = _hash(kmer_s + i, _ksize);
	set(kmer, iid);
	// cout << "SET: "<< kmer << " " << iid << "\n";
	std::string s = seq.substr(i, _ksize);
	// cout << "S is " << s << "\n";
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

  return;

  cout << "remapping...\n";
  remap();
  cout << "done\n";

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
	crossed.insert(iid);
      }


      if (crossed.size() > 1) {
	iid = find_min_iid(crossed);
	for (IntersectionSet::const_iterator si = crossed.begin();
	     si != crossed.end(); si++ ){
	  if (*si != iid) {
	    intermap[*si] = iid;
	  }
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
	crossed.insert(iid);
      }

      assert(crossed.size() == 1);

      outfp << ">" << read.name << "\t" << *(crossed.begin()) << "\n"
	      << seq << "\n";
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
}

IntersectionID IntersectTable::find_min_iid(IntersectionSet s)
{
  IntersectionSet::const_iterator si = s.begin();
  IntersectionInterMap::const_iterator pi;

  IntersectionID min_iid = *si;

  si++;
  for (; si != s.end(); si++) {
    IntersectionID iid = *si;
	
    while (1) {
      pi = intermap.find(iid);
      if (pi == intermap.end()) { break; }
      iid = intermap[iid];
    }
	
    if (iid < min_iid) { min_iid = iid; }
  }

  return min_iid;
}

void IntersectTable::remap()
{
  IntersectionSet seen;

  for (HashIntoType i = 0; i < _tablesize; i++) {
    IntersectionID iid = _table[i];
    if (iid) {
      if (intermap.find(iid) != intermap.end()) {
	_table[i] = intermap[iid];
      }
      seen.insert(_table[i]);
    }
  }
  intermap.clear();

  cout << "total intersection IDs seen: " << seen.size() << "\n";
}
