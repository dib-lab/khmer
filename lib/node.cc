#include "node.hh"

enum { A, C, G, T };

char getNextNucl(int base) {
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

bool Node::operator== (const Node &param) const {
   if (!(kmer == param.kmer)) {
      return 0;
   }

   if (state != param.state) {
      return 0;
   }

   if (stateNo != param.stateNo) {
      return 0;
   }

   return 1;
}

bool Node::operator< (const Node &param) const {
   if (fval < param.fval) {
      return 1;
   }
   return 0;
}

Kmer Node::makeNextKmer(unsigned char forward, char b) {
   HashIntoType ret_h;
   HashIntoType ret_r;

   const unsigned int rc_left_shift = kmer.getK()*2 - 2;

   if (forward) {
      ret_h = next_f(kmer.getH(), b);
      ret_r = next_r(kmer.getR(), b);
   } else {
      ret_h = prev_f(kmer.getH(), b);
      ret_r = prev_r(kmer.getR(), b);
   }

   return Kmer(ret_h, ret_r, kmer.getDir(), kmer.getK());
}

Node::Node(Node* _parent, 
           char _emission, 
           unsigned int _stateNo, 
           char _state,
           Kmer _kmer) {
   parent = _parent;
   kmer = _kmer;
   emission = _emission;
   stateNo = _stateNo;
   state = _state;

   fval = 0;
   gval = 0;
   hval = 0;

   diff = 0;

   bitmask = 0;
   for (unsigned int i = 0; i < kmer.getK(); i++) {
      bitmask = (bitmask << 2) | 3;
   }
}

std::queue<Node*> Node::enumerate(CountingHash* ch,
                         ScoringMatrix* sm,
                         unsigned char forward,
                         const std::string& seq) {
   std::queue<Node*> ret;

   int index;
   int remaining;

   if (forward) {
      index = stateNo + 1;
      remaining = kmer.getK();      
   } else {
      index = stateNo - 1;
      remaining = index;
   }

   // loop for matches and insertions
   for (int i = A; i <= T; i++) {
      char nextNucl = getNextNucl(i);
      Kmer nextKmer = makeNextKmer(forward, nextNucl);

      if (!ch->get_count(nextKmer.getUniqueHash())) {
         continue;
      }

      // match
      Node * nextMatch = new Node(this,
                                nextNucl,
                                index,
                                'm',
                                nextKmer);
      nextMatch->gval = gval + sm->score(nextNucl, seq[index]);
      nextMatch->hval = 0;
      nextMatch->fval = nextMatch->gval + nextMatch->hval;
      
      if (nextNucl == seq[index]) {
         nextMatch->diff = diff;
      } else {
         nextMatch->diff = diff + 1;
      }

      if (nextMatch->diff <= 3) {
         ret.push(nextMatch);
      } else {
         delete nextMatch;
      }

      // insertion
      if (state != 'd') {
         Node * nextIns = new Node(this,
                                 nextNucl,
                                 stateNo,
                                 'i',
                                 nextKmer);
         nextIns->gval = gval + sm->score(nextNucl, '-');
         nextIns->hval = 0;
         nextIns->fval = nextIns->gval + nextIns->hval;

         nextIns->diff = diff + 1;

         if (nextIns->diff <= 3) {
            ret.push(nextIns);
         } else {
            delete nextIns;
         }
      }
   }

   // deletion
   if (state != 'i') {
      Node * nextDel = new Node(this,
                          '-',
                          index,
                          'd',
                          kmer);
      nextDel->gval = gval + sm->score('-', seq[index]);
      nextDel->hval = 0;
      nextDel->fval = nextDel->gval + nextDel->hval;

      nextDel->diff = diff + 1;

      if (nextDel->diff <= 3) {
         ret.push(nextDel);
      } else {
         delete nextDel;
      }
   }

   return ret;
}
