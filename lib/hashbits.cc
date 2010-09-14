#include "hashtable.hh"
#include "hashbits.hh"

using namespace std;
using namespace khmer;

#if 1
HashIntoType khmer::primes[] = { 22906493,
				 22906519,
				 22906561,
				 22906567,
				 22906619,
				 22906649,
				 22906657,
				 22906661 };
#endif

#if 0
HashIntoType khmer::primes[] = { 32000000017,
				 32000000023,
				 32000000059,
				 32000000093,
				 32000000113,
				 32000000177,
				 32000000213,
				 32000000237 };
#endif

#if 0
HashIntoType khmer::primes[] = { 16000000039,
				 16000000067,
				 16000000091,
				 16000000097,
				 16000000103,
				 16000000127,
				 16000000157,
				 16000000201 };
#endif

//

void Hashbits::save(std::string outfilename)
{
  assert(_counts[0]);

  unsigned int save_ksize = _ksize;
  unsigned long long save_tablesize = _tablesize;

  ofstream outfile(outfilename.c_str(), ios::binary);

  outfile.write((const char *) &save_ksize, sizeof(save_ksize));

  for (unsigned int i = 0; i < N_TABLES; i++) {
    outfile.write((const char *) &save_tablesize, sizeof(save_tablesize));

    outfile.write((const char *) _counts[i],
		  sizeof(BoundedCounterType) * _tablebytes);
  }
  outfile.close();
}

void Hashbits::load(std::string infilename)
{
  if (_counts) {
    for (unsigned int i = 0; i < N_TABLES; i++) {
      delete _counts[i]; _counts[i] = NULL;
    }
  }
  
  unsigned int save_ksize = 0;
  unsigned long long save_tablesize = 0;

  ifstream infile(infilename.c_str(), ios::binary);
  infile.read((char *) &save_ksize, sizeof(save_ksize));
  _ksize = (WordLength) save_ksize;

  for (unsigned int i = 0; i < N_TABLES; i++) {
    infile.read((char *) &save_tablesize, sizeof(save_tablesize));

    _tablesize = (HashIntoType) save_tablesize;
    _tablebytes = _tablesize / 8 + 1;
    _counts[i] = new BoundedCounterType[_tablebytes];

    unsigned long long loaded = 0;
    while (loaded != _tablebytes) {
      infile.read((char *) _counts[i], _tablebytes - loaded);
      loaded += infile.gcount();	// do I need to do this loop?
    }
  }
  infile.close();
}

