#ifndef INTERTABLE_HH
#define INTERTABLE_HH

#include <fstream>
#include <string>
#include <set>
#include <map>
#include <queue>

#include "khmer.hh"

#define MAX_IID 50*1000*1000	// times 4 bytes for mem consumption

namespace khmer {
  typedef unsigned int IntersectionID;
  typedef unsigned int PartitionID;
  typedef std::map <IntersectionID, IntersectionID> IntersectionInterMap;
  typedef std::set<IntersectionID> IntersectionSet;
  typedef std::map<IntersectionID, IntersectionSet*> ReverseIntersectionMap;
  
  class IntersectTable {
  protected:
    WordLength _ksize;
    HashIntoType _tablesize;

    IntersectionID * _table;
    PartitionID * partitions;

    IntersectionID next_intersection_id;
    ReverseIntersectionMap revmap;

    bool rolled;
  public:
    IntersectTable(WordLength ksize, HashIntoType tablesize) :
      _ksize(ksize), _tablesize(tablesize) {
      next_intersection_id = 1;
      _allocate_counters();

      rolled = false;
    }

    ~IntersectTable() {
      delete _table;
    }

    void _allocate_counters() {
      _table = new IntersectionID[_tablesize];
      memset(_table, 0, _tablesize * sizeof(IntersectionID));

      partitions = new PartitionID[MAX_IID];
      for (unsigned int i = 0; i < MAX_IID; i++) {
	partitions[i] = i;
      }
    }

    // accessor to get 'k'
    const WordLength ksize() const { return _ksize; }

    // accessors to get table info
    const HashIntoType tablesize() const { return _tablesize; }

    IntersectionID get_next_iid() {
      IntersectionID iid = next_intersection_id;

      next_intersection_id++;
      if (next_intersection_id == MAX_IID) {
	next_intersection_id = 1;
      }
      return iid;
    }

    void set(HashIntoType kmer, IntersectionID iid) {
      _table[kmer % _tablesize] = iid;
    }

    IntersectionID get(HashIntoType kmer) {
      return _table[kmer % _tablesize];
    }

    bool check_read(const std::string &read);
    void do_partition(const std::string infile, const std::string outfile);

    void remap();
  };
};

#endif
