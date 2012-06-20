class ScoringMatrix {
   int sm[5][5];

private:
   int assign(char);

public:
   ScoringMatrix() {
      sm = 
      {
      { 5, -4, -4, -4, -2},
      {-4,  5, -4, -4, -2},
      {-4, -4,  5, -4, -2},
      {-4, -4, -4,  5, -2},
      {-2, -2, -2, -2, -1}
      };
   }
   
   int score(char, char);
};
