#ifndef NODEENUMERATOR_HH
#define NODEENUMERATOR_HH

#include <string>
#include <queue>
#include "astarsearchnode.hh"
#include "scoringmatrix.hh"
#include "hashbits.hh"

class NodeEnumerator {
   ScoringMatrix * sm;
   
   int bestMatch;
   AStarSearchNode * next;
   std::string nextKmer;
   int index;
   char forward;
   std::string seq;
   char nextNucl;

private:
   char getNextNucl(int);
   std::string makeNextKmer(char, char, AStarSearchNode*);

public:
   NodeEnumerator(char _forward, std::string _seq, ScoringMatrix* _sm) {
      forward = _forward;
      seq = _seq;
      sm = _sm;
      bestMatch = 5;
   }

   std::queue<AStarSearchNode*> enumerateNodes(AStarSearchNode * curr,
                                        khmer::Hashbits *);
};

#endif //NODEENUMERATOR_HH
