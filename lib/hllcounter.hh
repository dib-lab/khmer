//
// This file is part of khmer, http://github.com/ged-lab/khmer/, and is
// Copyright (C) Michigan State University, 2009-2014. It is licensed under
// the three-clause BSD license; see doc/LICENSE.txt. Contact: ctb@msu.edu
//

#ifndef HLLCOUNTER_HH
#define HLLCOUNTER_HH

#include <vector>
#include <string>
#include "khmer_config.hh"

/*
#ifdef __cplusplus
extern "C" {
#endif
extern const int THRESHOLD_DATA[] = {
  10, 20, 40, 80, 220, 400, 900, 1800, 3100,
  6500, 11500, 20000, 50000, 120000, 350000 };

#ifdef __cplusplus
}
#endif
*/

namespace khmer
{
class HLLCounter
{
public:
    HLLCounter(double error_rate);

    void add(const std::string &);
    HashIntoType estimate_cardinality();
    virtual ~HLLCounter() {}
private:
    double _Ep();
    double alpha;
    int p;
    int m;
    std::vector<int> M;
};
};

#endif // HLLCOUNTER_HH

// vim: set sts=2 sw=2:
