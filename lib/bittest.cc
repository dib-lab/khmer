#include <iostream>
#include "ktable.hh"

using namespace std;

int main()
{
  cout << khmer::_hash("ACGTTGCAACGTTGCAA", 17) << endl;
  cout << khmer::_hash("TTGCAACGTTGCAACGT", 17) << endl;

  cout << khmer::_hash("ACGTTGCAACGTTGCAA", 16) << endl;
  cout << khmer::_hash("ACGTTGCAACGTTGCAA", 17) << endl;

  return 0;
}
