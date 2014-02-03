//
// This file is part of khmer, http://github.com/ged-lab/khmer/, and is
// Copyright (C) Michigan State University, 2009-2013. It is licensed under
// the three-clause BSD license; see doc/LICENSE.txt. Contact: ctb@msu.edu
//

#ifndef ALIGNER_HH
#define ALIGNER_HH

#include "node.hh"
#include "counting.hh"

#include <algorithm>
#include <set>
#include <vector>

namespace khmer {

class CandidateAlignment {
public:
   std::map<int,int> readDeletions;
   std::string alignment;

   CandidateAlignment(std::map<int,int> _readDel, std::string _aln)
   {
      readDeletions = _readDel;
      alignment = _aln;
   }

   CandidateAlignment() {
      alignment = "";
   }

   bool operator<(const CandidateAlignment& param) const {
      return alignment < param.alignment;
   }

   std::string getReadAlignment(std::string seq) {
      int tmpDels = 0;
      int readIndex = 0;
      std::string readAlign;

      for (int i = 0; i < (int)alignment.length(); i++) {
         if (tmpDels > 0) {
            readAlign += '-';
         } else {
            readAlign += seq[readIndex];
         }
   
         if (tmpDels == 0 && readDeletions.count(readIndex+1) == 0) {
            readIndex++;
         } else if (tmpDels > 0) {
            tmpDels--;

            if (tmpDels == 0) {
               readIndex++;
            }
         } else {
            tmpDels = readDeletions[readIndex+1];
         }
      }

      return readAlign;
   }

   void outputReadDels() {      
      std::map<int, int>::iterator it;

      for (it = readDeletions.begin(); it != readDeletions.end(); it++) {
         std::cout << "Key: " << (*it).first << " Value: " <<
                      (*it).second << std::endl;
      }
   }    
};

class Aligner {
   khmer::CountingHash * ch;
   ScoringMatrix * sm;
   int k;
   double lambdaOne; // error distribution parameter
   double lambdaTwo; // non-error distribution parameter
   unsigned int maxErrorRegion;

public:
   Node * subalign(Node *, unsigned int, unsigned char, std::set<Node*>&, 
                  std::vector<Node*>&, const std::string&);
   std::string extractString(Node*, unsigned char, std::map<int,int>*);
   CandidateAlignment align(khmer::CountingHash*, const std::string&, 
                            const std::string&, int);

   Aligner(khmer::CountingHash* _ch, 
           double lOne=0.0, double lTwo=0.0,
	   unsigned int maxErrorReg=UINT_MAX) {
      ch = _ch;
      sm = new ScoringMatrix();
      k = ch->ksize();
      lambdaOne=lOne;
      lambdaTwo=lTwo;
      maxErrorRegion = maxErrorReg;
   }

   int ksize() {
      return k;
   }

   ~Aligner() {
      delete sm;
   }

   void printErrorFootprint(const std::string& read);
   CandidateAlignment alignRead(const std::string&);
};

};

#endif // ALIGNER_HH
