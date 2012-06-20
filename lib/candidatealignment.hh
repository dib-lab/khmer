#ifndef CANDIDATEALIGNMENT_HH
#define CANDIDATEALIGNMENT_HH

#include <string>
#include <map>

class CandidateAlignment {
public:
   std::map<int,int> readDeletions;
   std::string alignment;
   int score;

   CandidateAlignment(std::map<int,int> _readDel, std::string _aln, int _score)
   {
      readDeletions = _readDel;
      alignment = _aln;
      score = _score;
   }

   std::string getReadAlignment(std::string seq) {
      int tmpDels = 0;
      int readIndex = 0;
      std::string readAlign;

      for (int i = 0; i < (int)alignment.length(); i++) {
         if (tmpDels > 0) {
            readAlign += '-';
            tmpDels--;
         } else if (readDeletions.count(i) == 0) {
            readAlign += seq[readIndex];
            readIndex++;
         } else {
            readAlign += seq[readIndex];
            tmpDels = readDeletions[i];
            readIndex++;
         }
      }

      return readAlign;
   }
};

#endif // CANDIDATEALIGNMENT_HH
