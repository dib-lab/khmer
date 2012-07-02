#ifndef READALIGNER_HH
#define READALIGNER_HH

#include "nodeenumerator.hh"
#include "candidatealignment.hh"
#include <algorithm>
#include <set>
#include <queue>
#include <deque>
#include <vector>
#include "counting.hh"

class ReadAligner {
   khmer::CountingHash * ch;
   ScoringMatrix * sm;
   int k;

private:
//   AStarSearchNode* subalign(AStarSearchNode*, NodeEnumerator*, int,
//                             std::set<AStarSearchNode*>*);
   AStarSearchNode* subalign(AStarSearchNode*, NodeEnumerator*, int,
                             std::set<AStarSearchNode*>*, std::string);
   std::string extractString(AStarSearchNode*, char, std::map<int,int>*);
   CandidateAlignment* align(khmer::CountingHash*, std::string, std::string, int);
   void clearSet(std::set<AStarSearchNode*>* s);

public:
   ReadAligner(khmer::CountingHash* _ch) {
      ch = _ch;
      sm = new ScoringMatrix();
      k = ch->ksize();
   }

   int ksize() {
      return k;
   }

   ~ReadAligner() {
      delete sm;
   }

   CandidateAlignment alignRead(std::string);
};

#endif //READALIGNER_HH
