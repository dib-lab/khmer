#include "astarsearchnode.hh"
#include <iostream>

AStarSearchNode::AStarSearchNode(AStarSearchNode * _discoveredFrom,
                                 char _emission,
                                 int _stateNo,
                                 char _state,
                                 std::string _kmer) {
   discoveredFrom = _discoveredFrom;
   emission = _emission;
   kmer = _kmer;
   stateNo = _stateNo;
   state = _state;
   score = 0.0;
   deletes = 0;
   snps = 0;
}

bool AStarSearchNode::operator== (const AStarSearchNode &param) {
   if (kmer != param.kmer) {
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

bool AStarSearchNode::operator< (const AStarSearchNode &param) {
   if (fval < param.fval) {
      return 1;
   }
   return 0;
}
