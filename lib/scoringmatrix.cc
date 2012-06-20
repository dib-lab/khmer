#include "scoringmatrix.hh"

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

int ScoringMatrix::score(char ref, char qry) {
   int r = assign(ref);
   int q = assign(qry);

   return ScoringMatrix::sm[r][q];
}
