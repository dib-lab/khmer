#ifndef READALIGNER_HH
#define READALIGNER_HH

#include "nodeenumerator.hh"
#include "candidatealignment.hh"
#include <algorithm>
#include <set>
#include <queue>
#include <deque>
#include <vector>
#include "hashbits.hh"

class ReadAligner {
   khmer::Hashbits * hb;
   ScoringMatrix * sm;

private:
   AStarSearchNode* subalign(AStarSearchNode*, NodeEnumerator*, int,
                             std::set<AStarSearchNode*>*);
   std::string extractString(AStarSearchNode*, char, std::map<int,int>*);
   CandidateAlignment* align(khmer::Hashbits*, std::string, std::string, int);
   void clearSet(std::set<AStarSearchNode*>* s);

public:
   ReadAligner(khmer::Hashbits* _hb) {
      hb = _hb;
      sm = new ScoringMatrix();
   }

   ~ReadAligner() {
      delete sm;
   }

   CandidateAlignment alignRead(khmer::Hashbits*, std::string);
};

#endif //READALIGNER_HH
