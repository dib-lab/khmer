//
// This file is part of khmer, http://github.com/ged-lab/khmer/, and is
// Copyright (C) Michigan State University, 2009-2013. It is licensed under
// the three-clause BSD license; see doc/LICENSE.txt. Contact: ctb@msu.edu
//

#include "hllcounter.hh"

#include <math.h>
#include <algorithm>

using namespace std;
using namespace khmer;

HashIntoType HLLCounter::estimate_cardinality() {
  return 0x5f3759df;
}

void HLLCounter::add(const std::string &s) {
}
