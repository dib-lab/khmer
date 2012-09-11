#ifndef ALIGNER_HH
#define ALIGNER_HH

#include "node.hh"
#include "counting.hh"

#include <algorithm>
#include <set>
#include <vector>

class CandidateAlignment {
public:
   std::map<int,int> readDeletions;
   std::string alignment;
   double score;

   CandidateAlignment(std::map<int,int> _readDel, std::string _aln, 
                      double _score)
   {
      readDeletions = _readDel;
      alignment = _aln;
      score = _score;
   }

   CandidateAlignment() {
      alignment = "";
      score = -1.0; // a negative score isn't possible
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

public:
   Node * subalign(Node *, unsigned int, unsigned char, std::set<Node*>&, 
                  std::vector<Node*>&, const std::string&);
   std::string extractString(Node*, unsigned char, std::map<int,int>*);
   CandidateAlignment align(khmer::CountingHash*, const std::string&, 
                            const std::string&, int);

   Aligner(khmer::CountingHash* _ch, double lOne=0, double lTwo=0) {
      ch = _ch;
      sm = new ScoringMatrix();
      k = ch->ksize();
      lambdaOne=lOne;
      lambdaTwo=lTwo;
   }

   int ksize() {
      return k;
   }

   ~Aligner() {
      delete sm;
   }

   CandidateAlignment alignRead(const std::string&);
};

#endif // ALIGNER_HH
