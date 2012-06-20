#include "nodeenumerator.hh"

enum { A, C, G, T };

char NodeEnumerator::getNextNucl(int base) {
   if (base == A) {
      return 'A';
   } else if (base == C) {
      return 'C';
   } else if (base == G) {
      return 'G';
   } else if (base == T) {
      return 'T';
   }
   return '-';
} 

std::queue<AStarSearchNode*> NodeEnumerator::enumerateNodes(AStarSearchNode * curr,
                                                    khmer::Hashbits * hb) {
   std::queue<AStarSearchNode*> ret;

   int remaining;
   if (forward) {
      index = curr->stateNo + 1;
      remaining = seq.length() - index;
   } else {
      index = curr->stateNo - 1;
      remaining = index;
   }

   // loop for matches and insertions
   for (int i = A; i <= T; i++) {
      nextNucl = getNextNucl(i);
      nextKmer = makeNextKmer(forward, nextNucl, curr);

      if (!hb->get_count(nextKmer.c_str())) {
         continue;
      }

      next = new AStarSearchNode::AStarSearchNode(curr, nextNucl, index, 'm', nextKmer);
      next->score = curr->score + sm->score(nextNucl, seq[index]);
      next->hcost = remaining * bestMatch;
      next->fval = (int) next->score + next->hcost;
      next->deletes = curr->deletes;
      ret.push(next);

      next = new AStarSearchNode(curr, nextNucl, curr->stateNo, 'i', nextKmer);
      next->score = curr->score + sm->score(nextNucl, '-');
      next->hcost = remaining * bestMatch;
      next->fval = (int) next->score + next->hcost;
      next->deletes = curr->deletes;
      ret.push(next); 
   }

   // now handle deletion
   next = new AStarSearchNode(curr, '-', index, 'd', curr->kmer);
   next->score = curr->score + sm->score('-', seq[index]);
   next->hcost = remaining * bestMatch;
   next->fval = (int) next->score + next->hcost;
   next->deletes = curr->deletes + 1;
   ret.push(next);

   return ret;
}

std::string NodeEnumerator::makeNextKmer(char forward, char b, AStarSearchNode * curr) {
   std::string kmer;

   if (forward) {
      kmer = curr->kmer.substr(1, curr->kmer.length() - 1) + b;
   } else {
      kmer = b + curr->kmer.substr(0, curr->kmer.length() - 1);
   }

   return kmer;
}
