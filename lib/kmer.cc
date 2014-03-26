//
// This file is part of khmer, http://github.com/ged-lab/khmer/, and is
// Copyright (C) Michigan State University, 2009-2013. It is licensed under
// the three-clause BSD license; see doc/LICENSE.txt. 
// Contact: khmer-project@idyll.org
//

#include "kmer.hh"

#include <iostream>
using namespace std;

namespace khmer
{

Kmer::Kmer(std::string kmer)
{
    _hash(kmer.c_str(), kmer.size(), h, r);
    k = kmer.size();

    if (kmer == _revhash(h, k)) {
        direction = 1;
    } else {
        direction = 0;
    }
}

Kmer::Kmer(HashIntoType _h, HashIntoType _r, unsigned char _direction,
           unsigned int _k)
{
    h = _h;
    r = _r;
    direction = _direction;
    k = _k;
}

unsigned char Kmer::getDir()
{
    return direction;
}

HashIntoType Kmer::getUniqueHash() const
{
    return uniqify_rc(h, r);
}

unsigned int Kmer::getK()
{
    return k;
}

HashIntoType Kmer::getH()
{
    return h;
}

HashIntoType Kmer::getR()
{
    return r;
}

std::string Kmer::toString()
{
    return _revhash(uniqify_rc(h, r), k);
}

std::string Kmer::toStringH()
{
    return _revhash(h, k);
}

std::string Kmer::toStringR()
{
    return _revhash(r, k);
}

bool Kmer::operator== (const Kmer &b) const
{
    if (getUniqueHash() == b.getUniqueHash()) {
        return 1;
    } else {
        return 0;
    }
}

bool Kmer::operator< (const Kmer &b) const
{
    if (getUniqueHash() < b.getUniqueHash()) {
        return 1;
    } else {
        return 0;
    }
}
};
