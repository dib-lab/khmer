//
// This file is part of khmer, http://github.com/ged-lab/khmer/, and is
// Copyright (C) Michigan State University, 2009-2013. It is licensed under
// the three-clause BSD license; see doc/LICENSE.txt. Contact: ctb@msu.edu
//

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
