#include <iostream>

using namespace std;

#define twobit_repr(ch) ((toupper(ch)) == 'A' ? 0 : \
                          (toupper(ch)) == 'C' ? 1 : \
                          (toupper(ch)) == 'G' ? 2 : 3)

#define twobit_comp(ch) ((toupper(ch)) == 'A' ? 3 : \
                          (toupper(ch)) == 'C' ? 2 : \
                          (toupper(ch)) == 'G' ? 1 : 0)

//
// _hash: hash a k-length DNA sequence into an unsigned int.
//

unsigned long long int _hash(const char * kmer, unsigned int k)
{
  unsigned long long int h = 0;
  unsigned long long int r = 0;

  h |= twobit_repr(kmer[0]);
  r |= twobit_comp(kmer[k-1]);

  for (unsigned int i = 1; i < k; i++) {
    h = h << 2;
    r = r << 2;

    h |= twobit_repr(kmer[i]);
    r |= twobit_comp(kmer[k-1-i]);
  }

  cout << h << endl;
  cout << r << endl;

  return h < r ? h : r;
}

int main()
{
   cout << _hash("ACGTTGCAACGTTGCAA", 17) << endl;
   cout << _hash("TTGCAACGTTGCAACGT", 17) << endl;

   return 0;
}
