#ifndef SCORINGMATRIX_HH
#define SCORINGMATRIX_HH

class ScoringMatrix {
   double sm[5][5];

private:
   int assign(char);

public:
   ScoringMatrix() {
      sm =
      {
      { 0, -9, -9, -9, -7},
      {-9,  0, -9, -9, -7},
      {-9, -9,  0, -9, -7},
      {-9, -9, -9,  0, -7},
      {-7, -7, -7, -7, -6}
      };
   }
   
   int score(char, char);
};

#endif // SCORINGMATRIX_HH
