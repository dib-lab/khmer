#ifndef KMER_HH
#define KMER_HH

#include "khmer.hh"
#include "ktable.hh"
#include "hashbits.hh"

#include <set>
#include <string>

using namespace khmer;

class Kmer;

class Kmer {
private:
   HashIntoType h;
   HashIntoType r;
   unsigned int k;
   unsigned char direction; // 1 if h is forward

public:
   Kmer(std::string kmer);
   Kmer(HashIntoType, HashIntoType, unsigned char, unsigned int);
   Kmer() { }

   HashIntoType getUniqueHash() const;

   unsigned int getK();
   HashIntoType getH();
   HashIntoType getR();
   unsigned char getDir();

   std::string toString();
   std::string toStringH();
   std::string toStringR();
   std::set<Kmer> getNeighbors();
   std::set<Kmer> getForwardNeighbors();
   std::set<Kmer> getBackwardNeighbors();
   std::set<Kmer> getForwardStates(CountingHash*);
   std::set<Kmer> getBackwardStates(CountingHash*);

   bool operator== (const Kmer& b) const;
   bool operator< (const Kmer& b) const;
};

#endif // KMER_HH
