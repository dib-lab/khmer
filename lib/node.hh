#ifndef NODE_HH
#define NODE_HH

#include "kmer.hh"
#include "counting.hh"
#include "scoringmatrix.hh"
#include <queue>

class Node;

class Node {
public:
   Node *parent;
   Kmer kmer;
   char emission; // change to enum
   unsigned int stateNo;
   char state; //change to enum
   int fval;
   int hval;
   int gval;

   unsigned int diff;

   HashIntoType bitmask;

   bool operator== (const Node &param) const;
   bool operator< (const Node &param) const;

   Node(Node* _parent, char _emission, 
        unsigned int _stateNo, char _state, Kmer _kmer);

   Kmer makeNextKmer(unsigned char forward, char b);

   std::queue<Node*> enumerate(CountingHash *, ScoringMatrix * sm, 
                            unsigned char forward, const std::string&);
};

class NodeCompare {
public:
   bool operator()(Node* o1, Node* o2) {
      if (o1->fval < o2-> fval) {
         return true;
      } else {
         return false;
      }
   }
};

#endif //FILE_HH
