/*
This file is part of khmer, https://github.com/dib-lab/khmer/, and is
Copyright (C) 2010-2015, Michigan State University.
Copyright (C) 2015-2016, The Regents of the University of California.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are
met:

    * Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.

    * Redistributions in binary form must reproduce the above
      copyright notice, this list of conditions and the following
      disclaimer in the documentation and/or other materials provided
      with the distribution.

    * Neither the name of the Michigan State University nor the names
      of its contributors may be used to endorse or promote products
      derived from this software without specific prior written
      permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
"AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
LICENSE (END)

Contact: khmer-project@idyll.org
*/
#include <stddef.h>
#include <stdint.h>
#include <string.h>
#include <algorithm>
#include <string>

#include "khmer.hh"
#include "kmer_hash.hh"

using namespace std;

namespace khmer
{


KmerIterator::KmerIterator(const char * seq,
                           unsigned char k) :
    _hash(k), _seq(seq)
{
    bitmask = 0;
    for (unsigned char i = 0; i < _ksize; i++) {
        bitmask = (bitmask << 2) | 3;
    }
    _nbits_sub_1 = (_ksize*2 - 2);

    index = _ksize - 1;
    length = strlen(_seq);
    _kmer_f = 0;
    _kmer_r = 0;

    initialized = false;
}

Kmer KmerIterator::first()
{
    index = _ksize;
    Kmer start = _hash(_seq);
    _kmer_f = start.kmer_f;
    _kmer_r = start.kmer_r;
    return start;
}

Kmer KmerIterator::next()
{
    if (done()) {
        throw khmer_exception();
    }

    if (!initialized) {
        initialized = true;
        return first();
    }

    unsigned char ch = _seq[index];
    index++;
    if (!(index <= length)) {
        throw khmer_exception();
    }

    // left-shift the previous hash over
    _kmer_f = _kmer_f << 2;

    // 'or' in the current nt
    _kmer_f |= twobit_repr(ch);

    // mask off the 2 bits we shifted over.
    _kmer_f &= bitmask;

    // now handle reverse complement
    _kmer_r = _kmer_r >> 2;
    _kmer_r |= (twobit_comp(ch) << _nbits_sub_1);

    return Kmer(_kmer_f, kmer_r);
}

}
