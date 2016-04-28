#ifndef KMER_MIN_HASH_HH
#define KMER_MIN_HASH_HH

#include <set>
#include <map>

#include "MurmurHash3.h"
#include "khmer.hh"

////

#define DEFAULT_MINHASH_PRIME 9999999967

namespace khmer {

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

  void _shrink()
  {
    while (mins.size() > num) {
      CMinHashType::iterator mi = mins.end();
      mi--;
      mins.erase(mi);
    }
  }
  void add_hash(long int h)
  {
    h = ((h % prime) + prime) % prime;
    mins.insert(h);
    _shrink();
  }
  void add_kmer(std::string kmer)
  {
    long int hash = _hash_murmur32(kmer);
    add_hash(hash);
  }
  void add_sequence(const char * sequence) {
    if (!is_protein) {
      std::string seq = sequence;
      for (unsigned int i = 0; i < seq.length() - ksize + 1; i++) {
        std::string kmer = seq.substr(i, ksize);
        add_kmer(kmer);
      }
      std::string rc = _revcomp(seq);
      for (unsigned int i = 0; i < rc.length() - ksize + 1; i++) {
        std::string kmer = rc.substr(i, ksize);
        add_kmer(kmer);
      }
    } else {                      // protein
      std::string seq = sequence;
      for (unsigned int i = 0; i < seq.length() - ksize + 1; i ++) {
        std::string kmer = seq.substr(i, ksize);
        std::string aa = _dna_to_aa(kmer);

        add_kmer(aa);
        
        std::string rc = _revcomp(kmer);
        aa = _dna_to_aa(rc);

        add_kmer(aa);
      }
    }
  }
  int _hash_murmur32(const std::string& kmer)
  {
    int out[2];
    uint32_t seed = 0;
    MurmurHash3_x86_32((void *)kmer.c_str(), kmer.size(), seed, &out);
    return out[0];
  }

  std::string _dna_to_aa(const std::string& dna)
  {
    std::string aa;
    unsigned int dna_size = (dna.size() / 3) * 3; // floor it
    for (unsigned int j = 0; j < dna_size; j += 3) {
      std::string codon = dna.substr(j, 3);
      aa += (_codon_table)[codon];
    }
    return aa;
  }

  std::string _revcomp(const std::string& kmer)
  {
    std::string out = kmer;
    size_t ksize = out.size();

    for (size_t i=0; i < ksize; ++i) {
        char complement;

        switch(kmer[i]) {
        case 'A':
            complement = 'T';
            break;
        case 'C':
            complement = 'G';
            break;
        case 'G':
            complement = 'C';
            break;
        case 'T':
            complement = 'A';
            break;
        default:
            throw std::exception();
            break;
        }
        out[ksize - i - 1] = complement;
    }
    return out;
  }
  
  void merge(const KmerMinHash& other)
  {
    CMinHashType::iterator mi;
    for (mi = other.mins.begin(); mi != other.mins.end(); ++mi) {
      mins.insert(*mi);
    }
    _shrink();
  }
  unsigned int count_common(const KmerMinHash& other)
  {
    CMinHashType combined;
    
    CMinHashType::iterator mi;
    for (mi = mins.begin(); mi != mins.end(); ++mi) {
      combined.insert(*mi);
    }
    for (mi = other.mins.begin(); mi != other.mins.end(); ++mi) {
      combined.insert(*mi);
    }
    return mins.size() + other.mins.size() - combined.size();
  }

private:
  std::map<std::string, std::string> _codon_table = {
       {"TTT", "F"}, {"TTC", "F"},
       {"TTA", "L"}, {"TTG", "L"},

       {"TCT", "S"}, {"TCC", "S"}, {"TCA", "S"}, {"TCG", "S"},

       {"TAT", "Y"}, {"TAC", "Y"},
       {"TAA", "*"}, {"TAG", "*"},

       {"TGT", "C"}, {"TGC", "C"},
       {"TGA", "*"},
       {"TGG", "W"},

       {"CTT", "L"}, {"CTC", "L"}, {"CTA", "L"}, {"CTG", "L"},

       {"CCT", "P"}, {"CCC", "P"}, {"CCA", "P"}, {"CCG", "P"},

       {"CAT", "H"}, {"CAC", "H"},
       {"CAA", "Q"}, {"CAG", "Q"},

       {"CGT", "R"}, {"CGC", "R"}, {"CGA", "R"}, {"CGG", "R"},

       {"ATT", "I"}, {"ATC", "I"}, {"ATA", "I"},
       {"ATG", "M"},

       {"ACT", "T"}, {"ACC", "T"}, {"ACA", "T"}, {"ACG", "T"},

       {"AAT", "N"}, {"AAC", "N"},
       {"AAA", "K"}, {"AAG", "K"},

       {"AGT", "S"}, {"AGC", "S"},
       {"AGA", "R"}, {"AGG", "R"},

       {"GTT", "V"}, {"GTC", "V"}, {"GTA", "V"}, {"GTG", "V"},

       {"GCT", "A"}, {"GCC", "A"}, {"GCA", "A"}, {"GCG", "A"},

       {"GAT", "D"}, {"GAC", "D"},
       {"GAA", "E"}, {"GAG", "E"},

       {"GGT", "G"}, {"GGC", "G"}, {"GGA", "G"}, {"GGG", "G"}
  };
};

typedef std::map<khmer::HashIntoType, khmer::TagSet> TagToTagSet;
typedef std::map<khmer::HashIntoType, khmer::KmerMinHash *> TagToMinHash;

class NeighborhoodMinHash {
public:
  TagToTagSet tag_connections;
  TagToMinHash tag_to_mh;

  ~NeighborhoodMinHash() {
    cleanup_neighborhood_hash();
  }

  void cleanup_neighborhood_hash() {
    for (TagToMinHash::iterator mhi = tag_to_mh.begin();
         mhi != tag_to_mh.end(); mhi++) {
      delete mhi->second;
      mhi->second = NULL;
    }
  }
};

class CombinedMinHash {
public:
  TagSet tags;
  KmerMinHash * mh;

  CombinedMinHash() : mh(NULL) { }

  ~CombinedMinHash() { cleanup(); }
  
  void cleanup() {
    delete mh;
    mh = NULL;
  }
};

}
#endif // KMER_MIN_HASH_HH
