#include "ktable.hh"
#include "seqfuncs.hh"
#include <iostream>
#include <string>

using namespace std;

int main()
{
  khmer::KTable a(1);
  assert(a.n_entries() == 4);
  assert(a.max_hash() == 3);

  khmer::KTable b(2);
  assert(b.n_entries() == 16);
  assert(b.max_hash() == 15);

  khmer::KTable c(5);
  assert(c.n_entries() == 1024);
  assert(c.max_hash() == 1023);

  khmer::KTable d(6);
  assert(d.n_entries() == 4096);
  assert(d.max_hash() == 4096 - 1);

#if 0
  assert(d.get_count("AAAAAA") == 0);
  d.count("AAAAAA");
  assert(d.get_count("AAAAAA") == 1);

  assert(d.get_count("ATCGAT") == 0);
  d.count("ATCGAT");
  assert(d.get_count("ATCGAT") == 1);
#endif // 0

  unsigned int k = 6;
  assert(khmer::_hash(khmer::_revhash(52, k).c_str(), k) == 52);

  d.clear();

  char * s = "ATCGATATGAGGATCCAGGATCAAGATAGACCAGATATGACCAGAG";
  d.consume_string(s);

  for (unsigned int i = 0; i < strlen(s) - k + 1; i++) {
    assert(d.get_count(&s[i]) >= 1);
  }

  string str1 = getReverseComplement("CCGTA");
  string str2 = string("TACGG");
  assert(str1.compare(str2) == 0); 
  
  printf("Tests SUCCESSFUL.\n");

  return 0;
}
