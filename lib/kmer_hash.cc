//
// This file is part of khmer, https://github.com/dib-lab/khmer/, and is
// Copyright (C) Michigan State University, 2009-2015. It is licensed under
// the three-clause BSD license; see LICENSE.
// Contact: khmer-project@idyll.org
//

#include <math.h>
#include <string>
#include <iostream>
#include <algorithm>

#include "khmer.hh"
#include "kmer_hash.hh"
#include "MurmurHash3.h"

using namespace std;

//
// _hash: hash a k-length DNA sequence into a 64-bit number.
//

namespace khmer
{

HashIntoType _hash(const char * kmer, const WordLength k,
                   HashIntoType& _h, HashIntoType& _r)
{
    // sizeof(HashIntoType) * 8 bits / 2 bits/base
    if (!(k <= sizeof(HashIntoType)*4) || !(strlen(kmer) >= k)) {
        throw khmer_exception("Supplied kmer string doesn't match the underlying k-size.");
    }

    HashIntoType h = 0, r = 0;

    h |= twobit_repr(kmer[0]);
    r |= twobit_comp(kmer[k-1]);

    for (WordLength i = 1, j = k - 2; i < k; i++, j--) {
        h = h << 2;
        r = r << 2;

        h |= twobit_repr(kmer[i]);
        r |= twobit_comp(kmer[j]);
    }

    _h = h;
    _r = r;

    return uniqify_rc(h, r);
}

// _hash: return the maximum of the forward and reverse hash.

HashIntoType _hash(const char * kmer, const WordLength k)
{
    HashIntoType h = 0;
    HashIntoType r = 0;

    return khmer::_hash(kmer, k, h, r);
}

// _hash_forward: return the hash from the forward direction only.

HashIntoType _hash_forward(const char * kmer, WordLength k)
{
    HashIntoType h = 0;
    HashIntoType r = 0;


    khmer::_hash(kmer, k, h, r);
    return h;			// return forward only
}

//
// _revhash: given an unsigned int, return the associated k-mer.
//

std::string _revhash(HashIntoType hash, WordLength k)
{
    std::string s = "";

    unsigned int val = hash & 3;
    s += revtwobit_repr(val);

    for (WordLength i = 1; i < k; i++) {
        hash = hash >> 2;
        val = hash & 3;
        s += revtwobit_repr(val);
    }

    reverse(s.begin(), s.end());

    return s;
}

std::string _revcomp(const std::string& kmer)
{
    std::string out = kmer;
    size_t ksize = out.size();

    for (size_t i=0; i < ksize; ++i) {
        char complement;

        switch(kmer[i]) {
        case 'A':
            complement = 'T';
            break;
        case 'C':
            complement = 'G';
            break;
        case 'G':
            complement = 'C';
            break;
        case 'T':
            complement = 'A';
            break;
        default:
            throw khmer::khmer_exception("Invalid base in read");
            break;
        }
        out[ksize - i - 1] = complement;
    }
    return out;
}

HashIntoType _hash_murmur(const std::string& kmer)
{
    HashIntoType h = 0;
    HashIntoType r = 0;

    return khmer::_hash_murmur(kmer, h, r);
}

HashIntoType _hash_murmur(const std::string& kmer,
                          HashIntoType& h, HashIntoType& r)
{
    HashIntoType out[2];
    uint32_t seed = 0;
    MurmurHash3_x64_128((void *)kmer.c_str(), kmer.size(), seed, &out);
    h = out[0];

    std::string rev = khmer::_revcomp(kmer);
    MurmurHash3_x64_128((void *)rev.c_str(), rev.size(), seed, &out);
    r = out[0];

    return h ^ r;
}

HashIntoType _hash_murmur_forward(const std::string& kmer)
{
    HashIntoType h = 0;
    HashIntoType r = 0;

    khmer::_hash_murmur(kmer, h, r);
    return h;
}

};
