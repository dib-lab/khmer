/*
This file is part of khmer, https://github.com/dib-lab/khmer/, and is
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
#ifndef HIST_HH
#define HIST_HH

#include <functional>
#include <vector>
#include <cstddef>

#include "oxli.hh"
#include "kmer_hash.hh"

namespace oxli {

inline size_t highest_bit(uint64_t num)
{
    if (!num)
        return 0;

    int pos = 1;

    while (num >>= 1) {
        pos += 1;
    }

    return pos;
}


template <size_t n_bins>
class Histogram {

    public:

        uint64_t bins[n_bins];

        Histogram() {
            clear();
        }

        void add(uint64_t val) {
            size_t bin = highest_bit(val) - 1;
            if (bin >= n_bins) {
                bins[n_bins-1] += 1;
            } else {
                bins[bin] += 1;
            }
        }

        void clear() {
            for (auto&& b : bins) {
                b = 0;
            }
        }
};

template class Histogram<8>;
template class Histogram<16>;
template class Histogram<32>;
template class Histogram<64>;


}

#endif
