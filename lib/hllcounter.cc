//
// This file is part of khmer, https://github.com/dib-lab/khmer/, and is
// Copyright (C) Michigan State University, 2014-2015. It is licensed under
// the three-clause BSD license; see LICENSE.
// Contact: khmer-project@idyll.org
//

#include "hllcounter.hh"

#include <math.h>
#include <algorithm>
#include <numeric>
#include <inttypes.h>
#include <sstream>

#include "khmer.hh"
#include "kmer_hash.hh"
#include "read_parsers.hh"
#include "khmer_exception.hh"

#ifdef _OPENMP
#include <omp.h>
#else
#define omp_get_thread_num() 0
#define omp_get_num_threads() 1
#endif

#define arr_len(a) (a + sizeof a / sizeof a[0])

using namespace khmer;

std::map<int, std::vector<double> > rawEstimateData;
std::map<int, std::vector<double> > biasData;


double calc_alpha(const int p)
{
    if (p < 4) {
        // ceil(log2((1.04 / x) ^ 2)) = 4, solve for x
        throw InvalidValue("Please set error rate to a value "
                           "smaller than 0.367696");
    } else if (p > 16) {
        // ceil(log2((1.04 / x) ^ 2)) = 16, solve for x
        throw InvalidValue("Please set error rate to a value "
                           "greater than 0.0040624");
    }

    /*
       For a description of following constants see
       HyperLogLog in Practice: Algorithmic Engineering of a State of The Art
          Cardinality Estimation Algorithm
       Stefan Heule, Marc Nunkesser and Alex Hall
       dx.doi.org/10.1145/2452376.2452456
    */
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

void init_raw_estimate_data()
{
    if (rawEstimateData.empty()) {
        for(int i=4; i <= 18; i++) {
            std::vector<double> v;
            switch(i) {
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
            }
            rawEstimateData[i] = v;
        }
    }
}

void init_bias_data()
{
    if (biasData.empty()) {
        for(int i=4; i <= 18; i++) {
            std::vector<double> v;
            switch(i) {
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
            }
            biasData[i] = v;
        }
    }
}

double get_threshold(int p)
{
    return THRESHOLD_DATA[p - 4];
}

std::vector<int> get_nearest_neighbors(double E, std::vector<double> estimate)
{
    std::vector< std::pair<double,int> > distance_map;
    std::vector<int> nearest;

    int i = 0;
    for (std::vector<double>::iterator it = estimate.begin();
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
    std::vector<double> bias = biasData[p];
    std::vector<double> raw_estimate = rawEstimateData[p];

    std::vector<int> nearest = get_nearest_neighbors(E, raw_estimate);
    double estimate = 0.0;

    for (std::vector<int>::iterator it = nearest.begin();
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
    return max_width - floor(log2(w));
}

HLLCounter::HLLCounter(double error_rate, WordLength ksize)
{
    if (error_rate < 0) {
        throw InvalidValue("Please set error rate to a value "
                           "greater than zero");
    }
    int p = ceil(log2(pow(1.04 / error_rate, 2)));
    this->init(p, ksize);
}

HLLCounter::HLLCounter(int p, WordLength ksize)
{
    this->init(p, ksize);
}

void HLLCounter::init(int p, WordLength ksize)
{
    this->alpha = calc_alpha(p);
    this->p = p;
    this->_ksize = ksize;
    this->m = 1 << p;
    std::vector<int> M(this->m, 0.0);
    this->M = M;

    init_raw_estimate_data();
    init_bias_data();
}

double HLLCounter::get_erate()
{
    return 1.04 / sqrt(this->m);
}

void HLLCounter::set_erate(double error_rate)
{
    if (count(this->M.begin(), this->M.end(), 0) != this->m) {
        throw ReadOnlyAttribute("You can only change error rate prior to "
                                "first counting");
    }

    if (error_rate < 0) {
        throw InvalidValue("Please set error rate to a value "
                           "greater than zero");
    }
    int p = ceil(log2(pow(1.04 / error_rate, 2)));
    this->init(p, this->_ksize);
}

void HLLCounter::set_ksize(WordLength new_ksize)
{
    if (count(this->M.begin(), this->M.end(), 0) != this->m) {
        throw ReadOnlyAttribute("You can only change k-mer size prior to "
                                "first counting");
    }

    this->init(this->p, new_ksize);
}

double HLLCounter::_Ep()
{
    double sum = accumulate(this->M.begin(), this->M.end(), 0.0, ep_sum);
    double E = this->alpha * pow(this->m, 2.0) / sum;

    if (E <= (5 * (double)this->m)) {
        return E - estimate_bias(E, this->p);
    }

    return E;
}

HashIntoType HLLCounter::estimate_cardinality()
{
    long V = count(this->M.begin(), this->M.end(), 0);

    if (V > 0) {
        double H = this->m * log((double)this->m / V);
        if (H <= get_threshold(this->p)) {
            return H;
        }
    }
    return this->_Ep();
}

void HLLCounter::add(const std::string &value)
{
    HashIntoType x = khmer::_hash_murmur(value);
    HashIntoType j = x & (this->m - 1);
    this->M[j] = std::max(this->M[j], get_rho(x >> this->p, 64 - this->p));
}

unsigned int HLLCounter::consume_string(const std::string &inp)
{
    unsigned int n_consumed = 0;
    std::string kmer = "";
    std::string s = inp;

    for (unsigned int i = 0; i < s.length(); i++)  {
        s[i] &= 0xdf; // toupper - knock out the "lowercase bit"
    }

    for(std::string::const_iterator it = s.begin(); it != s.end(); ++it) {
        kmer.push_back(*it);
        if (kmer.size() < _ksize) {
            continue;
        }
        this->add(kmer);

        kmer.erase(0, 1);
        n_consumed++;
    }
    return n_consumed;
}

void HLLCounter::consume_fasta(
    std::string const &filename,
    unsigned int &total_reads,
    unsigned long long &n_consumed)
{
    read_parsers::IParser * parser = read_parsers::IParser::get_parser(filename);

    consume_fasta(parser, total_reads, n_consumed);

    delete parser;
}

void HLLCounter::consume_fasta(
    read_parsers::IParser *parser,
    unsigned int &      total_reads,
    unsigned long long &    n_consumed)
{

    read_parsers::Read read;
    HLLCounter** counters;
    unsigned int *n_consumed_partial;
    unsigned int *total_reads_partial;

    n_consumed = 0;

    #pragma omp parallel
    {
        #pragma omp single
        {
            counters = (HLLCounter**)calloc(omp_get_num_threads(),
            sizeof(HLLCounter*));
            n_consumed_partial = (unsigned int*)calloc(omp_get_num_threads(),
            sizeof(unsigned int));
            total_reads_partial = (unsigned int*)calloc(omp_get_num_threads(),
            sizeof(unsigned int));

            for (int i=0; i < omp_get_num_threads(); i++)
            {
                HLLCounter *newc = new HLLCounter(this->p, this->_ksize);
                counters[i] = newc;
            }

            while (!parser->is_complete())
            {
                // Iterate through the reads and consume their k-mers.
                try {
                    read = parser->get_next_read();
                } catch (read_parsers::NoMoreReadsAvailable) {
                    break;
                }

                #pragma omp task default(none) firstprivate(read) \
                shared(counters, n_consumed_partial, total_reads_partial)
                {
                    bool is_valid;
                    int n, t = omp_get_thread_num();
                    n = counters[t]->check_and_process_read(read.sequence,
                    is_valid);
                    n_consumed_partial[t] += n;
                    if (is_valid) {
                        total_reads_partial[t] += 1;
                    }
                }

            } // while reads left for parser
        }
        #pragma omp taskwait

        #pragma omp single
        {
            for (int i=0; i < omp_get_num_threads(); ++i)
            {
                this->merge(*counters[i]);
                delete counters[i];
                n_consumed += n_consumed_partial[i];
                total_reads += total_reads_partial[i];;
            }
            free(counters);
            free(n_consumed_partial);
            free(total_reads_partial);
        }
    }
}

unsigned int HLLCounter::check_and_process_read(std::string &read,
        bool &is_valid)
{
    is_valid = check_and_normalize_read(read);

    if (!is_valid) {
        return 0;
    }

    return consume_string(read);
}

bool HLLCounter::check_and_normalize_read(std::string &read) const
{
    bool is_valid = true;

    if (read.length() < this->_ksize) {
        return false;
    }

    for (unsigned int i = 0; i < read.length(); i++) {
        read[ i ] &= 0xdf; // toupper - knock out the "lowercase bit"
        if (!is_valid_dna( read[ i ] )) {
            is_valid = false;
            break;
        }
    }

    return is_valid;
}

void HLLCounter::merge(HLLCounter &other)
{
    if (this->p != other.p || this->_ksize != other._ksize) {
        throw khmer_exception("HLLCounters to be merged must be created with same parameters");
    }
    for(unsigned int i=0; i < this->M.size(); ++i) {
        this->M[i] = std::max(other.M[i], this->M[i]);
    }
}
