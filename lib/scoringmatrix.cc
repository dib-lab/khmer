//
// This file is part of khmer, http://github.com/ged-lab/khmer/, and is
// Copyright (C) Michigan State University, 2009-2013. It is licensed under
// the three-clause BSD license; see doc/LICENSE.txt.
// Contact: khmer-project@idyll.org
//

#include "scoringmatrix.hh"
#include <cmath>

using namespace std;

double ScoringMatrix::score(char ref, char qry)
{
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
