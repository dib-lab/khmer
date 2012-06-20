#include "astarsearchnode.hh"

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
}
