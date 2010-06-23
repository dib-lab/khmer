#include <assert.h>
#include "ktable.hh"
#include <math.h>

using namespace khmer;

//
// _hash: hash a k-length DNA sequence into an unsigned int.
//

HashIntoType khmer::_hash(const char * kmer, WordLength k, 
			  HashIntoType * h, HashIntoType * r)
{
  *h |= twobit_repr(kmer[0]);
  *r |= twobit_comp(kmer[k-1]);

  for (WordLength i = 1; i < k; i++) {
    *h = *h << 2;
    *r = *r << 2;

    *h |= twobit_repr(kmer[i]);
    *r |= twobit_comp(kmer[k-1-i]);
  }

  return *h < *r ? *h : *r;
}

HashIntoType khmer::_hash(const char * kmer, WordLength k)
{
  HashIntoType h = 0;
  HashIntoType r = 0;

  
  return _hash(kmer, k, &h, &r);
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

#if 1
  const unsigned int length = s.length() - _ksize + 1;
  for (unsigned int i = 0; i < length; i++) {
    count(&sp[i]);
  }
#else

  unsigned int mask = 0;
  for (unsigned int i = 0; i < _ksize; i++) {
    mask = mask << 2;
    mask |= 3;
  }

  unsigned long long int h;
  unsigned long long int r;

  _hash(sp, _ksize, &h, &r);
  
  if (h < r)
     _counts[h]++;
  else
    _counts[r]++;

  for (unsigned int i = _ksize; i < length; i++) {
    short int repr = twobit_repr(sp[i]);

    // left-shift the previous hash over
    h = h << 2;

    // 'or' in the current nt
    h |= twobit_repr(sp[i]);

    // mask off the 2 bits we shifted over.
    h &= mask;

    // now handle reverse complement
    r = r << 2;
    r &= mask;
    r |= twobit_repr(sp[i]);

    if (h < r)
      _counts[h]++;
    else
      _counts[r]++;
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
