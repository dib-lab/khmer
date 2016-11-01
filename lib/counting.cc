/*
This file is part of khmer, https://github.com/dib-lab/khmer/, and is
Copyright (C) 2010-2015, Michigan State University.
Copyright (C) 2015, The Regents of the University of California.

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

#include <errno.h>
#include <algorithm>
#include <iostream>
#include <sstream> // IWYU pragma: keep

#include "counting.hh"
#include "hashbits.hh"
#include "hashtable.hh"
#include "khmer_exception.hh"
#include "kmer_hash.hh"
#include "read_parsers.hh"
#include "zlib.h"

using namespace std;
using namespace khmer;
using namespace khmer:: read_parsers;

BoundedCounterType CountingHash::get_min_count(const std::string &s)
{
    KmerIterator kmers(s.c_str(), _ksize);

    BoundedCounterType min_count = MAX_KCOUNT;

    while(!kmers.done()) {
        HashIntoType kmer = kmers.next();

        BoundedCounterType count = this->get_count(kmer);

        if (this->get_count(kmer) < min_count) {
            min_count = count;
        }
    }
    return min_count;
}

BoundedCounterType CountingHash::get_max_count(const std::string &s)
{
    KmerIterator kmers(s.c_str(), _ksize);

    BoundedCounterType max_count = 0;

    while(!kmers.done()) {
        HashIntoType kmer = kmers.next();

        BoundedCounterType count = this->get_count(kmer);

        if (count > max_count) {
            max_count = count;
        }
    }
    return max_count;
}

uint64_t *
CountingHash::abundance_distribution(
    read_parsers::IParser * parser,
    Hashbits *          tracking)
{
    uint64_t * dist = new uint64_t[MAX_BIGCOUNT + 1];
    uint64_t i;

    for (i = 0; i <= MAX_BIGCOUNT; i++) {
        dist[i] = 0;
    }

    Read read;

    string name;
    string seq;

    // if not, could lead to overflow.
    if (sizeof(BoundedCounterType) != 2) {
        delete[] dist;
        throw khmer_exception();
    }

    while(!parser->is_complete()) {
        try {
            read = parser->get_next_read();
        } catch (NoMoreReadsAvailable &exc) {
            break;
        }
        seq = read.sequence;

        if (check_and_normalize_read(seq)) {
            KmerIterator kmers(seq.c_str(), _ksize);

            while(!kmers.done()) {
                HashIntoType kmer = kmers.next();

                if (!tracking->get_count(kmer)) {
                    tracking->count(kmer);

                    BoundedCounterType n = get_count(kmer);
                    dist[n]++;
                }
            }

            name.clear();
            seq.clear();
        }
    }
    return dist;
}


uint64_t * CountingHash::abundance_distribution(
    std::string filename,
    Hashbits *  tracking)
{
    IParser* parser = IParser::get_parser(filename.c_str());

    uint64_t * distribution = abundance_distribution(parser, tracking);
    delete parser;
    return distribution;
}

unsigned long CountingHash::trim_on_abundance(
    std::string     seq,
    BoundedCounterType  min_abund)
const
{
    if (!check_and_normalize_read(seq)) {
        return 0;
    }

    KmerIterator kmers(seq.c_str(), _ksize);

    HashIntoType kmer;

    if (kmers.done()) {
        return 0;
    }
    kmer = kmers.next();

    if (kmers.done() || get_count(kmer) < min_abund) {
        return 0;
    }

    unsigned long i = _ksize;
    while (!kmers.done()) {
        kmer = kmers.next();

        if (get_count(kmer) < min_abund) {
            return i;
        }
        i++;
    }

    return seq.length();
}

unsigned long CountingHash::trim_below_abundance(
    std::string     seq,
    BoundedCounterType  max_abund)
const
{
    if (!check_and_normalize_read(seq)) {
        return 0;
    }

    KmerIterator kmers(seq.c_str(), _ksize);

    HashIntoType kmer;

    if (kmers.done()) {
        return 0;
    }
    kmer = kmers.next();

    if (kmers.done() || get_count(kmer) > max_abund) {
        return 0;
    }

    unsigned long i = _ksize;
    while (!kmers.done()) {
        kmer = kmers.next();

        if (get_count(kmer) > max_abund) {
            return i;
        }
        i++;
    }

    return seq.length();
}

std::vector<unsigned int> CountingHash::find_spectral_error_positions(
    std::string seq,
    BoundedCounterType max_abund)
const
{
    std::vector<unsigned int> posns;
    if (!check_and_normalize_read(seq)) {
        throw khmer_exception("invalid read");
    }

    KmerIterator kmers(seq.c_str(), _ksize);

    HashIntoType kmer = kmers.next();
    if (kmers.done()) {
        return posns;
    }

    // find the first trusted k-mer
    while (!kmers.done()) {
        if (get_count(kmer) > max_abund) {
            break;
        }
        kmer = kmers.next();
    }

    if (kmers.done()) {
        return posns;
    }

    // did we bypass some erroneous k-mers? call the last one.
    if (kmers.get_start_pos() > 0) {
        // if we are well past the first k, forget the whole thing (!? @CTB)
        if (kmers.get_start_pos() >= _ksize && 0) {
            return posns;
        }
        posns.push_back(kmers.get_start_pos() - 1);
    }

    while (!kmers.done()) {
        kmer = kmers.next();
        if (get_count(kmer) <= max_abund) { // error!
            posns.push_back(kmers.get_end_pos() - 1);

            // find next good
            while (!kmers.done()) {
                kmer = kmers.next();
                if (get_count(kmer) > max_abund) { // a good stretch again.
                    break;
                }
            }
        }
    }

    return posns;
}

/* vim: set ft=cpp ts=8 sts=4 sw=4 et tw=79 */
