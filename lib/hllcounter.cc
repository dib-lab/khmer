//
// This file is part of khmer, http://github.com/ged-lab/khmer/, and is
// Copyright (C) Michigan State University, 2009-2013. It is licensed under
// the three-clause BSD license; see doc/LICENSE.txt. Contact: ctb@msu.edu
//

#include "hllcounter.hh"

#include <math.h>
#include <algorithm>
#include <numeric>
#include <inttypes.h>

#include "khmer.hh"
#include "kmer_hash.hh"

#define arr_len(a) (a + sizeof a / sizeof a[0])

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

vector<double> raw_estimate_data(int p) {
    /* TODO avoid all this repetition/copying around */
    vector<double> v;
    switch(p) {
      case 4:
        v.assign(RAW_ESTIMATE_DATA_4, arr_len(RAW_ESTIMATE_DATA_4));
        break;
      case 5:
        v.assign(RAW_ESTIMATE_DATA_5, arr_len(RAW_ESTIMATE_DATA_5));
        break;
      case 6:
        v.assign(RAW_ESTIMATE_DATA_6, arr_len(RAW_ESTIMATE_DATA_6));
        break;
      case 7:
        v.assign(RAW_ESTIMATE_DATA_7, arr_len(RAW_ESTIMATE_DATA_7));
        break;
      case 8:
        v.assign(RAW_ESTIMATE_DATA_8, arr_len(RAW_ESTIMATE_DATA_8));
        break;
      case 9:
        v.assign(RAW_ESTIMATE_DATA_9, arr_len(RAW_ESTIMATE_DATA_9));
        break;
      case 10:
        v.assign(RAW_ESTIMATE_DATA_10, arr_len(RAW_ESTIMATE_DATA_10));
        break;
      case 11:
        v.assign(RAW_ESTIMATE_DATA_11, arr_len(RAW_ESTIMATE_DATA_11));
        break;
      case 12:
        v.assign(RAW_ESTIMATE_DATA_12, arr_len(RAW_ESTIMATE_DATA_12));
        break;
      case 13:
        v.assign(RAW_ESTIMATE_DATA_13, arr_len(RAW_ESTIMATE_DATA_13));
        break;
      case 14:
        v.assign(RAW_ESTIMATE_DATA_14, arr_len(RAW_ESTIMATE_DATA_14));
        break;
      case 15:
        v.assign(RAW_ESTIMATE_DATA_15, arr_len(RAW_ESTIMATE_DATA_15));
        break;
      case 16:
        v.assign(RAW_ESTIMATE_DATA_16, arr_len(RAW_ESTIMATE_DATA_16));
        break;
      case 17:
        v.assign(RAW_ESTIMATE_DATA_17, arr_len(RAW_ESTIMATE_DATA_17));
        break;
      case 18:
        v.assign(RAW_ESTIMATE_DATA_18, arr_len(RAW_ESTIMATE_DATA_18));
        break;
      default:
        /* TODO raise exception */
        break;
    }
    return v;
}

vector<double> raw_bias_data(int p) {
    /* TODO avoid all this repetition/copying around */
    vector<double> v;
    switch(p) {
      case 4:
        v.assign(RAW_BIAS_DATA_4, arr_len(RAW_BIAS_DATA_4));
        break;
      case 5:
        v.assign(RAW_BIAS_DATA_5, arr_len(RAW_BIAS_DATA_5));
        break;
      case 6:
        v.assign(RAW_BIAS_DATA_6, arr_len(RAW_BIAS_DATA_6));
        break;
      case 7:
        v.assign(RAW_BIAS_DATA_7, arr_len(RAW_BIAS_DATA_7));
        break;
      case 8:
        v.assign(RAW_BIAS_DATA_8, arr_len(RAW_BIAS_DATA_8));
        break;
      case 9:
        v.assign(RAW_BIAS_DATA_9, arr_len(RAW_BIAS_DATA_9));
        break;
      case 10:
        v.assign(RAW_BIAS_DATA_10, arr_len(RAW_BIAS_DATA_10));
        break;
      case 11:
        v.assign(RAW_BIAS_DATA_11, arr_len(RAW_BIAS_DATA_11));
        break;
      case 12:
        v.assign(RAW_BIAS_DATA_12, arr_len(RAW_BIAS_DATA_12));
        break;
      case 13:
        v.assign(RAW_BIAS_DATA_13, arr_len(RAW_BIAS_DATA_13));
        break;
      case 14:
        v.assign(RAW_BIAS_DATA_14, arr_len(RAW_BIAS_DATA_14));
        break;
      case 15:
        v.assign(RAW_BIAS_DATA_15, arr_len(RAW_BIAS_DATA_15));
        break;
      case 16:
        v.assign(RAW_BIAS_DATA_16, arr_len(RAW_BIAS_DATA_16));
        break;
      case 17:
        v.assign(RAW_BIAS_DATA_17, arr_len(RAW_BIAS_DATA_17));
        break;
      case 18:
        v.assign(RAW_BIAS_DATA_18, arr_len(RAW_BIAS_DATA_18));
        break;
      default:
        /* TODO raise exception */
        break;
    }
    return v;
}

double get_threshold(int p)
{
    return THRESHOLD_DATA[p - 4];
}

vector<int> get_nearest_neighbors(double E, vector<double> estimate)
{
    vector< pair<double,int> > distance_map;
    vector<int> nearest;

    int i = 0;
    for (vector<double>::iterator it = estimate.begin();
         it != estimate.end();
         ++it) {
      std::pair<double, int> p(pow(E - *it, 2.0), i);
      distance_map.push_back(p);
      i++;
    }

    sort(distance_map.begin(), distance_map.end());

    for(int k=0; k < 6; k++) {
        nearest.push_back(distance_map[k].second);
    }

    return nearest;
}

double estimate_bias(double E, int p)
{
    vector<double> bias = raw_bias_data(p);
    vector<double> raw_estimate = raw_estimate_data(p);
    vector<int> nearest = get_nearest_neighbors(E, raw_estimate);
    double estimate = 0.0;

    for (vector<int>::iterator it = nearest.begin();
         it != nearest.end();
         ++it) {
      estimate += bias[*it];
    }
    return estimate / nearest.size();
}

double ep_sum(double acc, int b)
{
    return acc += pow(2.0, float(-b));
}

int get_rho(HashIntoType w, int max_width)
{
    int rho = max_width - floor(log2(w)); /* TODO find last bit set */

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
    vector<int> M(m, 0.0);

    this->alpha = get_alpha(p);
    this->p = p;
    this->m = 1 << p;
    this->M = M;
}

double HLLCounter::_Ep()
{
    double sum = accumulate(this->M.begin(), this->M.end(), 0.0, ep_sum);
    double E = this->alpha * pow(this->m, 2.0) / sum;

    if (E <= (5 * (double)this->m))
      return E - estimate_bias(E, this->p);

    return E;
}

HashIntoType HLLCounter::estimate_cardinality()
{
    double H;
    int V = count(this->M.begin(), this->M.end(), 0);

    if (V > 0) {
        H = this->m * log((double)this->m / V);
        if (H <= get_threshold(this->p)) {
            return H;
        }
    }
    return this->_Ep();
}

void HLLCounter::add(const string &value)
{
    //HashIntoType x = khmer::_hash(value.c_str(), value.size());
    //HashIntoType x = khmer::_hash_forward(value.c_str(), value.size());

    //HashIntoType x = khmer::_hash_murmur(value);
    //HashIntoType x = khmer::_hash_murmur_forward(value);

    HashIntoType x = khmer::_hash_sha1(value);
    //HashIntoType x = khmer::_hash_sha1_forward(value);

    HashIntoType j = x & (this->m - 1);
    this->M[j] = max(this->M[j], get_rho(x >> this->p, 64 - this->p));
}

unsigned int HLLCounter::consume_string(const std::string &s, unsigned int ksize) {
    unsigned int n_consumed = 0;
    std::string kmer = "";

    for(std::string::const_iterator it = s.begin(); it != s.end(); ++it) {
        kmer.push_back(*it);
        if (kmer.size() < ksize) {
            continue;
        }
        this->add(kmer);

        kmer.erase(0, 1);
        n_consumed++;
    }
    return n_consumed;
}
