//
// This file is part of khmer, http://github.com/ged-lab/khmer/, and is
// Copyright (C) Michigan State University, 2009-2014. It is licensed under
// the three-clause BSD license; see doc/LICENSE.txt. 
//Contact: khmer-project@idyll.org
//
#include <assert.h>
#include <math.h>
#include <string>
#include <iostream>
#include <algorithm>

#include "khmer.h"
#include "ktable.h"

using namespace std;
using namespace khmer;

//
// _hash: hash a k-length DNA sequence into a 64-bit number.
//

HashIntoType khmer::_hash(const char * kmer, const WordLength k, 
			  HashIntoType& _h, HashIntoType& _r)
{
  // sizeof(HashIntoType) * 8 bits / 2 bits/base  
  assert(k <= sizeof(HashIntoType)*4);

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

HashIntoType khmer::_hash(const char * kmer, const WordLength k)
{
  HashIntoType h = 0;
  HashIntoType r = 0;

  return _hash(kmer, k, h, r);
}

// _hash_forward: return the hash from the forward direction only.

HashIntoType khmer::_hash_forward(const char * kmer, WordLength k)
{
  HashIntoType h = 0;
  HashIntoType r = 0;

  
  _hash(kmer, k, h, r);
  return h;			// return forward only
}

//
// _revhash: given an unsigned int, return the associated k-mer.
//

std::string khmer::_revhash(HashIntoType hash, WordLength k)
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

//
// consume_string: run through every k-mer in the given string, & hash it.
//

void KTable::consume_string(const std::string &s)
{
  const char * sp = s.c_str();
  unsigned int length = s.length();

#if 0
  const unsigned int length = s.length() - _ksize + 1;
  for (unsigned int i = 0; i < length; i++) {
    count(&sp[i]);
  }
#else

  unsigned long long int mask = 0;
  for (unsigned int i = 0; i < _ksize; i++) {
    mask = mask << 2;
    mask |= 3;
  }

  HashIntoType h;
  HashIntoType r;

  _hash(sp, _ksize, h, r);
  
  _counts[uniqify_rc(h, r)]++;

  for (unsigned int i = _ksize; i < length; i++) {
    // left-shift the previous hash over
    h = h << 2;

    // 'or' in the current nt
    h |= twobit_repr(sp[i]);

    // mask off the 2 bits we shifted over.
    h &= mask;

    // now handle reverse complement
    r = r >> 2;
    r |= (twobit_comp(sp[i]) << (_ksize*2 - 2));

    _counts[uniqify_rc(h, r)]++;
  }

#endif // 0
}

void KTable::update(const KTable &other)
{
  assert(_ksize == other._ksize);

  for (unsigned int i = 0; i < n_entries(); i++) {
    _counts[i] += other._counts[i];
  }
}

KTable * KTable::intersect(const KTable &other) const
{
  assert(_ksize == other._ksize);

  KTable * intersection = new KTable(_ksize);

  for (unsigned int i = 0; i < n_entries(); i++) {
    if (_counts[i] > 0 && other._counts[i] > 0) {
      intersection->_counts[i] = _counts[i] + other._counts[i];
    }
  }
  return intersection;
}
