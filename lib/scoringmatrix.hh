//
// This file is part of khmer, http://github.com/ged-lab/khmer/, and is
// Copyright (C) Michigan State University, 2009-2013. It is licensed under
// the three-clause BSD license; see doc/LICENSE.txt.
// Contact: khmer-project@idyll.org
//

#ifndef SCORINGMATRIX_HH
#define SCORINGMATRIX_HH

#define MAT 0
#define SNP 1
#define INS 2
#define DEL 3

class ScoringMatrix
{
    //double sm[5][5];
    double scores[4];

private:

public:
    ScoringMatrix() {
        /*
        scores[MAT] = .96;
        scores[SNP] = .02;
        scores[INS] = .01;
        scores[DEL] = .01;
        */
        scores[MAT] = 0.0;
        scores[SNP] = 7.0;
        scores[INS] = 4.0;
        scores[DEL] = 4.0;
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
