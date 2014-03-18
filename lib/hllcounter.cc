//
// This file is part of khmer, http://github.com/ged-lab/khmer/, and is
// Copyright (C) Michigan State University, 2009-2013. It is licensed under
// the three-clause BSD license; see doc/LICENSE.txt. Contact: ctb@msu.edu
//

#include "hllcounter.hh"

#include <math.h>
#include <algorithm>
#include <numeric>

#include "khmer.hh"
#include "ktable.hh"

using namespace std;
using namespace khmer;

double get_alpha(const int p)
{
    if ((p < 4) or (p > 16)) {
        return 0;    // TODO: treat error
    }

    switch (p) {
    case 4:
        return 0.673;
    case 5:
        return 0.697;
    case 6:
        return 0.709;
    default:
        return 0.7213 / (1.0 + 1.079 / (1 << p));
    }
}

double get_threshold(int p)
{
    static const int THRESHOLD_DATA[] = {
        10, 20, 40, 80, 220, 400, 900, 1800, 3100,
        6500, 11500, 20000, 50000, 120000, 350000
    };

    return THRESHOLD_DATA[p - 4];
}

double ep_sum(double acc, int b)
{
    return acc += pow(2.0, -b);
}

int get_rho(HashIntoType w, int max_width)
{
    int rho = max_width - floor(log2(w)) + 1; /* TODO find last bit set */

    if (rho <= 0) {
        return -1;    // TODO: treat error, w overflow
    }

    return rho;
}

HLLCounter::HLLCounter(double error_rate)
{
    // TODO: check if 0 < error_rate < 1
    int p = ceil(log2(pow(1.04 / error_rate, 2)));
    int m = 1 << p;
    std::vector<int> M(m, 0.0);

    this->alpha = get_alpha(p);
    this->p = p;
    this->m = 1 << p;
    this->M = M;
}

double HLLCounter::_Ep()
{
    double sum = std::accumulate(this->M.begin(), this->M.end(), 0.0, ep_sum);
    double E = this->alpha * pow(this->m, 2) / sum;

    /*
      if (E <= (5 * this->m))
        return E - estimate_bias(E, this->p)
    */
    return E;
}

HashIntoType HLLCounter::estimate_cardinality()
{
    double H;
    int V = std::count(this->M.begin(), this->M.end(), 0);

    if (V > 0) {
        H = this->m * log(this->m / V);
        if (H <= get_threshold(this->p)) {
            return H;
        }
    }

    return this->_Ep();
}

void HLLCounter::add(const std::string &value) {
  HashIntoType j;
  HashIntoType x = khmer::_hash(value.c_str(), 32);
  j = x & (this->m - 1);

    this->M[j] = max(this->M[j], get_rho(x >> this->p, 64 - this->p));
}
