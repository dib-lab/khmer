#include "hashtable.hh"
#include "hashbits.hh"

using namespace std;
using namespace khmer;

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

