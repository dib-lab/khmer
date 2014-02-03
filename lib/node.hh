//
// This file is part of khmer, http://github.com/ged-lab/khmer/, and is
// Copyright (C) Michigan State University, 2009-2013. It is licensed under
// the three-clause BSD license; see doc/LICENSE.txt. Contact: ctb@msu.edu
//

#ifndef NODE_HH
#define NODE_HH

#include "kmer.hh"
#include "counting.hh"
#include "scoringmatrix.hh"
#include <queue>

#include <cmath>

namespace khmer
{

extern unsigned int factorial(unsigned int n);
extern double pois(double l, unsigned int k);
extern double weight_nonerror(unsigned int kCov, double lambOne, double lambTwo);
extern bool isCorrectKmer(unsigned int kCov, double lambOne, double lambTwo);

class Node;

class Node
{
public:
    Node *parent;
    Kmer kmer;
    char emission; // change to enum
    unsigned int stateNo;
    char state; //change to enum
    double fval;
    double hval;
    double gval;

    unsigned int diff;

    HashIntoType bitmask;

    bool operator== (const Node &param) const;
    bool operator< (const Node &param) const;

    Node(Node* _parent, char _emission,
         unsigned int _stateNo, char _state, Kmer _kmer);

    Kmer makeNextKmer(unsigned char forward, char b);

    std::queue<Node*> enumerate(CountingHash *, ScoringMatrix * sm,
                                unsigned char forward, const std::string&,
                                double lambdaOne=0, double lambdaTwo=0);
};

class NodeCompare
{
public:
    bool operator()(Node* o1, Node* o2) {
        if (o1->fval > o2->fval) {
            return true;
        } else {
            return false;
        }
    }
};

};
#endif //FILE_HH
