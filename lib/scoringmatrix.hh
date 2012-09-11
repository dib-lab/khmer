#ifndef SCORINGMATRIX_HH
#define SCORINGMATRIX_HH

#define MAT 0
#define SNP 1
#define INS 2
#define DEL 3

class ScoringMatrix {
   //double sm[5][5];
   double probs[4]; 

private:
   int assign(char);

public:
   ScoringMatrix() {
      probs[MAT] = .98;
      probs[SNP] = .0066;
      probs[INS] = .0066;
      probs[DEL] = .0066;
      /*
      sm =
      {
      {0, 9, 9, 9, 7},
      {9, 0, 9, 9, 7},
      {9, 9, 0, 9, 7},
      {9, 9, 9, 0, 7},
      {7, 7, 7, 7, 6}
      };
      */
//      sm = 
//      {
//      {0.03, 6.64, 6.64, 6.64, 8.23},
//      {6.64, 0.03, 6.64, 6.64, 8.23},
//      {6.64, 0.03, 6.64, 6.64, 8.23},
//      {6.64, 6.64, 6.64, 0.03, 8.23},
//      {8.23, 8.23, 8.23, 8.23, 10.0}
//      };
   }
   
   double score(char, char);
};

#endif // SCORINGMATRIX_HH
