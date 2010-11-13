#include <iostream>
#include "ktable.hh"

using namespace std;

int main()
{
  cout << khmer::_hash("ATGGACCAGATG", 12) << endl;
  cout << khmer::_hash("TGGACCAGATGA", 12) << endl;
  cout << khmer::_hash("GGACCAGATGAC", 12) << endl;
  cout << khmer::_hash("GACCAGATGACA", 12) << endl;
  cout << khmer::_hash("ACCAGATGACAC", 12) << endl;
  return 0;
}
