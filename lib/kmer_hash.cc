//
// This file is part of khmer, http://github.com/ged-lab/khmer/, and is
// Copyright (C) Michigan State University, 2009-2015. It is licensed under
// the three-clause BSD license; see doc/LICENSE.txt.
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
  // constants FTW
  
  // sizeof(HashIntoType) * 8 bits / 2 bits/base
  if (!(k <= sizeof(HashIntoType)*4) || !(strlen(kmer) >= k)) {
    throw khmer_exception("Supplied kmer string doesn't match the underlying k-size.");
  }

  HashIntoType h = 0, r = 0;

  h = _cyclichash(kmer, k);
  r = _revcyclichash(kmer, k);

  _h = h;
  _r = r;

  return uniqify_rc(h, r);
}

// Bitwise Shifting 
HashIntoType _cyclichash(std::string kmer_string, const WordLength k)
{
  // Our list of random numbers for each base value
  // Hard coding in random vals from array
  // Random values generated from www.random.org
  HashIntoType vals[] = {303837760, 392233993, 127908739, 491989698}; 
  //std::vector<HashIntoType> random_base_vals (vals, vals + sizeof(vals) / sizeof(vals));

  HashIntoType hashvalue = 0;

  // Iterate through the kmer, hashing each 
  for (int i = 0; i < kmer_string.length(); i++) 
  {
    // A circular bitshift (http://en.wikipedia.org/wiki/Circular_shift)
    hashvalue = (hashvalue << 1 | hashvalue >> (k - 1));
    
    // Hash each letter of the current kmer
    switch (kmer_string[i]) 
    {
      case 'A':
      case 'a':
        hashvalue ^= vals[0];
      case 'C':
      case 'c':
        hashvalue ^= vals[1];
      case 'G':
      case 'g':
        hashvalue ^= vals[2];
      case 'T':
      case 't':
      case 'U':
      case 'u':
        hashvalue ^= vals[3];
      default:
        throw khmer_exception(); 
    }// of Switch
  }
  return hashvalue;
}

HashIntoType _revcyclichash(std::string kmer_string, const WordLength k)
{
  // Our list of random numbers for each base value
  // Hard coding in random vals from array
  // Random values generated from www.random.org
  HashIntoType compliment_vals[] = {491989698, 127908739, 392233993, 303837760}; 
  //std::vector<HashIntoType> compliment_random_base_vals (vals, vals + sizeof(vals) / sizeof(vals));

  HashIntoType hashvalue = 0;

  // Iterate through the kmer backwards, hashing each base 
  for (int i = kmer_string.length(); i > -1; i--) 
  {
    // A circular bitshift (http://en.wikipedia.org/wiki/Circular_shift)
    hashvalue = (hashvalue << 1 | hashvalue >> (k - 1));
    
    // Hash each letter of the current kmer
    switch (kmer_string[i]) 
    {
      case 'A':
      case 'a':
        hashvalue ^= compliment_vals[0];
      case 'C':
      case 'c':
        hashvalue ^= compliment_vals[1];
      case 'G':
      case 'g':
        hashvalue ^= compliment_vals[2];
      case 'T':
      case 't':
      case 'U':
      case 'u':
        hashvalue ^= compliment_vals[3];
      default:
        throw khmer_exception(); 
    }// of Switch
  }
  return hashvalue;
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
