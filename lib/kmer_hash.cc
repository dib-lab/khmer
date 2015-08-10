//
// This file is part of khmer, https://github.com/dib-lab/khmer/, and is
// Copyright (C) Michigan State University, 2009-2015. It is licensed under
// the three-clause BSD license; see LICENSE.
// Contact: khmer-project@idyll.org
//

#include <stddef.h>
#include <stdint.h>
#include <string.h>
#include <algorithm>
#include <string>

#include "MurmurHash3.h"
#include "khmer.hh"
#include "khmer_exception.hh"
#include "kmer_hash.hh"

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

KmerIterator::KmerIterator(const char * seq,
                           unsigned char k) :
                           KmerFactory(k), _seq(seq)
{
    bitmask = 0;
    for (unsigned char i = 0; i < _ksize; i++) {
        bitmask = (bitmask << 2) | 3;
    }
    _nbits_sub_1 = (_ksize*2 - 2);

    index = _ksize - 1;
    length = strlen(seq);
    _kmer_f = 0;
    _kmer_r = 0;

    initialized = false;
}

Kmer KmerIterator::first(HashIntoType& f, HashIntoType& r)
{
    HashIntoType x;
    x = _hash(_seq, _ksize, _kmer_f, _kmer_r);

    f = _kmer_f;
    r = _kmer_r;

    index = _ksize;

    return Kmer(_kmer_f, _kmer_r, x);
}

Kmer KmerIterator::next(HashIntoType& f, HashIntoType& r)
{
    if (done()) {
        throw khmer_exception();
    }

    if (!initialized) {
        initialized = true;
        return first(f, r);
    }

    unsigned char ch = _seq[index];
    index++;
    if (!(index <= length)) {
        throw khmer_exception();
    }

    // left-shift the previous hash over
    _kmer_f = _kmer_f << 2;

    // 'or' in the current nt
    _kmer_f |= twobit_repr(ch);

    // mask off the 2 bits we shifted over.
    _kmer_f &= bitmask;

    // now handle reverse complement
    _kmer_r = _kmer_r >> 2;
    _kmer_r |= (twobit_comp(ch) << _nbits_sub_1);

    f = _kmer_f;
    r = _kmer_r;

    return build_kmer(_kmer_f, _kmer_r);
}

};
