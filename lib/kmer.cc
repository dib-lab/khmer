//
// This file is part of khmer, http://github.com/ged-lab/khmer/, and is
// Copyright (C) Michigan State University, 2009-2013. It is licensed under
// the three-clause BSD license; see doc/LICENSE.txt. Contact: ctb@msu.edu
//

#include "kmer.hh"

#include <iostream>
using namespace std;

Kmer::Kmer(std::string kmer) {
   khmer::_hash(kmer.c_str(), kmer.size(), h, r);   
   k = kmer.size();

   if (kmer == khmer::_revhash(h, k)) {
      direction = 1;
   } else {
      direction = 0;
   }
}

Kmer::Kmer(khmer::HashIntoType _h, khmer::HashIntoType _r, unsigned char _direction, 
           unsigned int _k) {
   h = _h;
   r = _r;
   direction = _direction;
   k = _k;
}

unsigned char Kmer::getDir() {
   return direction;
}

khmer::HashIntoType Kmer::getUniqueHash() const {
   return uniqify_rc(h, r);
}

unsigned int Kmer::getK() {
   return k;
}

khmer::HashIntoType Kmer::getH() {
   return h;
}

khmer::HashIntoType Kmer::getR() {
   return r;
}

std::string Kmer::toString() {
   return khmer::_revhash(uniqify_rc(h, r), k); 
}

std::string Kmer::toStringH() {
   return khmer::_revhash(h, k);
}

std::string Kmer::toStringR() {
   return khmer::_revhash(r, k);
}

bool Kmer::operator== (const Kmer &b) const {
   if (getUniqueHash() == b.getUniqueHash()) {
      return 1;
   } else {
      return 0;
   }
}

bool Kmer::operator< (const Kmer &b) const {
   if (getUniqueHash() < b.getUniqueHash()) {
      return 1;
   } else {
      return 0;
   }
}

/*
int main() {
   Kmer a = Kmer("AAAACCCC");
   Kmer b = Kmer("GGGGTTTT");   

   cout << a.toString() << endl;
   cout << a.toStringH() << endl;
   cout << a.toStringR() << endl;
   cout << (a == b) << endl;
   cout << (a < b) << endl;

   return 0;
}
*/
