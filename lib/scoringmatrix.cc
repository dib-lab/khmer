#include "scoringmatrix.hh"
#include <cmath>

using namespace std;

int ScoringMatrix::assign(char b) {
   if (b == 'A') {
      return 0;
   } else if (b == 'C') {
      return 1;
   } else if (b == 'G') {
      return 2;
   } else if (b == 'T') {
      return 3;
   }
   return 4;
}

double ScoringMatrix::score(char ref, char qry) {
   //int r = assign(ref);
   //int q = assign(qry);

   if (ref == qry) { // match
      //return 0-log2(probs[MAT]); 
      return scores[MAT];
   } else if (ref == '-') { // deletion
      //return 0-log2(probs[DEL]);
      return scores[DEL];
   } else if (qry == '-') { // insertion
      //return 0-(log2(probs[INS])/4.0);
      return scores[INS];
   } else { // snp
      //return 0-(log2(probs[SNP])/3.0);
      return scores[SNP];
   }   

   //return ScoringMatrix::sm[r][q];
}
