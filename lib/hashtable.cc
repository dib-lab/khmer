#include "khmer.hh"
#include "hashtable.hh"

using namespace khmer;

//
// consume_string: run through every k-mer in the given string, & hash it.
//

void Hashtable::consume_string(const std::string &s)
{
  const unsigned int length = s.length();
  const char * sp = s.c_str();

  unsigned int mask = 0;
  for (unsigned int i = 0; i < _ksize; i++) {
    mask = mask << 2;
    mask |= 3;
  }

#if 0
  for (unsigned int i = 0; i < s.length() - _ksize + 1; i++) {
    count(&sp[i]);
  }
#else

  unsigned int h = _hash(sp, _ksize);
  _counts[h % _tablesize]++;

  for (unsigned int i = _ksize; i < length; i++) {
    short int repr = twobit_repr(sp[i]);

    // left-shift the previous hash over
    h = h << 2;

    // 'or' in the current nt
    h |= twobit_repr(sp[i]);

    // mask off the 2 bits we shifted over.
    h &= mask;

    _counts[h % _tablesize]++;
  }

#endif // 0
}
