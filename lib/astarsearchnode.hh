#ifndef ASTARSEARCHNODE_HH
#define ASTARSEARCHNODE_HH

#include <string>
#include <iostream>

class AStarSearchNode;

class AStarSearchNode {
public:
   AStarSearchNode * discoveredFrom;
   char emission;
   std::string kmer;
   double score;
   char state;
   int stateNo;
   double hcost;
   int fval;
   char hasNewKmer; // boolean
   int deletes;
   int snps;
   double thisNodeScore;
   AStarSearchNode(AStarSearchNode*, char, int, char, std::string);
   bool operator== (const AStarSearchNode &param);
   bool operator< (const AStarSearchNode &param);
};

class ASSNCompare {
public:
   bool operator()(AStarSearchNode* o1, AStarSearchNode* o2) {
      if (o1->fval < o2->fval) {
         return true;
      }  else  {
         return false;
      }
   }
};

#endif //ASTARSEARCHNODE_HH
