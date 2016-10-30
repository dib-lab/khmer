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
#ifndef COUNTING_HH
#define COUNTING_HH

#include <stddef.h>
#include <stdint.h>
#include <string.h>
#include <map>
#include <string>
#include <utility>
#include <vector>

#include "hashtable.hh"
#include "khmer.hh"
#include "kmer_hash.hh"

namespace khmer
{
class Hashbits;

namespace read_parsers
{
struct IParser;
}  // namespace read_parsers
}  // namespace khmer

namespace khmer
{
class CountingHashFile;
class CountingHashFileReader;
class CountingHashFileWriter;
class CountingHashGzFileReader;
class CountingHashGzFileWriter;
class CountingHashIntersect;

class CountingHash : public khmer::Hashgraph
{
    friend class CountingHashIntersect;
    friend class CountingHashFile;
    friend class CountingHashFileReader;
    friend class CountingHashFileWriter;
    friend class CountingHashGzFileReader;
    friend class CountingHashGzFileWriter;

public:
    explicit CountingHash(WordLength ksize, std::vector<uint64_t> sizes)
        : Hashgraph(ksize, new ByteStorage(sizes)) { } ;

    void set_use_bigcount(bool b) { store->set_use_bigcount(b); }
    bool get_use_bigcount() { return store->get_use_bigcount(); }

    BoundedCounterType get_min_count(const std::string &s);

    BoundedCounterType get_max_count(const std::string &s);

    uint64_t * abundance_distribution(read_parsers::IParser * parser,
                                      Hashbits * tracking);
    uint64_t * abundance_distribution(std::string filename,
                                      Hashbits * tracking);

    unsigned long trim_on_abundance(std::string seq,
                                    BoundedCounterType min_abund) const;
    unsigned long trim_below_abundance(std::string seq,
                                       BoundedCounterType max_abund) const;
    std::vector<unsigned int> find_spectral_error_positions(std::string seq,
            BoundedCounterType min_abund) const;
};


class CountingHashFile
{
public:
    static void load(const std::string &infilename, CountingHash &ht);
    static void save(const std::string &outfilename, const CountingHash &ht);
};

class CountingHashFileReader : public CountingHashFile
{
public:
    CountingHashFileReader(const std::string &infilename, CountingHash &ht);
};

class CountingHashGzFileReader : public CountingHashFile
{
public:
    CountingHashGzFileReader(const std::string &infilename, CountingHash &ht);
};


class CountingHashFileWriter : public CountingHashFile
{
public:
    CountingHashFileWriter(const std::string &outfilename, const CountingHash &ht);
};

class CountingHashGzFileWriter : public CountingHashFile
{
public:
    CountingHashGzFileWriter(const std::string &outfilename,
                             const CountingHash &ht);
};
}

#endif // COUNTING_HH

// vim: set sts=2 sw=2:
