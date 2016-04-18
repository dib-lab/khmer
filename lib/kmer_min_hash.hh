#ifndef KMER_MIN_HASH_HH
#define KMER_MIN_HASH_HH

#include <set>

#include "MurmurHash3.h"
#include "khmer.hh"

////

typedef std::set<khmer::HashIntoType> CMinHashType;

class KmerMinHash
{
public:
  unsigned int num;
  unsigned int ksize;
  long int prime;
  bool is_protein;
  CMinHashType mins;

  KmerMinHash(unsigned int n, unsigned int k, long int p, bool prot) :
    num(n), ksize(k), prime(p), is_protein(prot) { };

  void add_hash(long int h)
  {
    h = ((h % prime) + prime) % prime;
    mins.insert(h);

    if (mins.size() > num) {
      CMinHashType::iterator mi = mins.end();
      mi--;
      mins.erase(mi);
    }
  }
  void add_kmer(std::string kmer)
  {
    long int hash = _hash_murmur32(kmer);
    add_hash(hash);
  }
  int _hash_murmur32(const std::string& kmer)
  {
    int out[2];
    uint32_t seed = 0;
    MurmurHash3_x86_32((void *)kmer.c_str(), kmer.size(), seed, &out);
    return out[0];
  }
};

#endif // KMER_MIN_HASH_HH
