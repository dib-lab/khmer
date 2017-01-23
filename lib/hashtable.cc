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
#include <errno.h>
#include <math.h>
#include <algorithm>
#include <deque>
#include <fstream>
#include <iostream>
#include <sstream> // IWYU pragma: keep
#include <queue>
#include <set>
#include <memory>

#include "hashtable.hh"
#include "khmer.hh"
#include "traversal.hh"
#include "read_parsers.hh"

using namespace std;
using namespace khmer;
using namespace khmer:: read_parsers;

//
// check_and_process_read: checks for non-ACGT characters before consuming
//

unsigned int Hashtable::check_and_process_read(std::string &read,
        bool &is_valid)
{
    is_valid = check_and_normalize_read(read);

    if (!is_valid) {
        return 0;
    }

    return consume_string(read);
}

//
// check_and_normalize_read: checks for non-ACGT characters
//			     converts lowercase characters to uppercase one
// Note: Usually it is desirable to keep checks and mutations separate.
//	 However, in the interests of efficiency (we are potentially working
//	 with TB of data), a check and mutation have been placed inside the
//	 same loop. Potentially trillions fewer fetches from memory would
//	 seem to be a worthwhile goal.
//

bool Hashtable::check_and_normalize_read(std::string &read) const
{
    bool rc = true;

    if (read.length() < _ksize) {
        return false;
    }

    for (unsigned int i = 0; i < read.length(); i++)  {
        read[ i ] &= 0xdf; // toupper - knock out the "lowercase bit"
        if (!is_valid_dna( read[ i ] )) {
            rc = false;
            break;
        }
    }

    return rc;
}

//
// consume_fasta: consume a FASTA file of reads
//

// TODO? Inline in header.
void
Hashtable::
consume_fasta(
    std:: string const  &filename,
    unsigned int	      &total_reads, unsigned long long	&n_consumed
)
{
    IParser *	  parser =
        IParser::get_parser( filename );

    consume_fasta(
        parser,
        total_reads, n_consumed
    );

    delete parser;
}

void
Hashtable::
consume_fasta(
    read_parsers:: IParser *  parser,
    unsigned int		    &total_reads, unsigned long long  &n_consumed
)
{
    Read			  read;

    // Iterate through the reads and consume their k-mers.
    while (!parser->is_complete( )) {
        bool is_valid;
        try {
            read = parser->get_next_read( );
        } catch (NoMoreReadsAvailable) {
            break;
        }

        unsigned int this_n_consumed =
            check_and_process_read(read.sequence, is_valid);

        __sync_add_and_fetch( &n_consumed, this_n_consumed );
        __sync_add_and_fetch( &total_reads, 1 );

    } // while reads left for parser

} // consume_fasta

//
// consume_string: run through every k-mer in the given string, & hash it.
//

unsigned int Hashtable::consume_string(const std::string &s)
{
    unsigned int n_consumed = 0;

    unique_ptr<KmerHashIterator> kmers = new_kmer_iterator(s);

    while(!kmers->done()) {
        HashIntoType kmer = kmers->next();

        count(kmer);
        n_consumed++;
    }

    return n_consumed;
}

// technically, get medioid count... our "median" is always a member of the
// population.

void Hashtable::get_median_count(const std::string &s,
                                 BoundedCounterType &median,
                                 float &average,
                                 float &stddev)
{
    std::vector<BoundedCounterType> counts;
    this->get_kmer_counts(s, counts);

    if (!counts.size()) {
        throw khmer_exception("no k-mer counts for this string; too short?");
    }

    average = 0;
    for (std::vector<BoundedCounterType>::const_iterator i = counts.begin();
            i != counts.end(); ++i) {
        average += *i;
    }
    average /= float(counts.size());

    stddev = 0;
    for (std::vector<BoundedCounterType>::const_iterator i = counts.begin();
            i != counts.end(); ++i) {
        stddev += (float(*i) - average) * (float(*i) - average);
    }
    stddev /= float(counts.size());
    stddev = sqrt(stddev);

    sort(counts.begin(), counts.end());
    median = counts[counts.size() / 2]; // rounds down
}

//
// Optimized filter function for normalize-by-median
//
bool Hashtable::median_at_least(const std::string &s,
                                unsigned int cutoff)
{
    unique_ptr<KmerHashIterator> kmers = new_kmer_iterator(s);
    unsigned int min_req = 0.5 + float(s.size() - _ksize + 1) / 2;
    unsigned int num_cutoff_kmers = 0;

    // first loop:
    // accumulate at least min_req worth of counts before checking to see
    // if we have enough high-abundance k-mers to indicate success.
    for (unsigned int i = 0; i < min_req; ++i) {
        HashIntoType kmer = kmers->next();
        if (this->get_count(kmer) >= cutoff) {
            ++num_cutoff_kmers;
        }
    }

    // second loop: now check to see if we pass the threshold for each k-mer.
    if (num_cutoff_kmers >= min_req) {
        return true;
    }
    while(!kmers->done()) {
        HashIntoType kmer = kmers->next();
        if (this->get_count(kmer) >= cutoff) {
            ++num_cutoff_kmers;
            if (num_cutoff_kmers >= min_req) {
                return true;
            }
        }
    }
    return false;
}

void Hashtable::get_kmers(const std::string &s,
                          std::vector<std::string> &kmers_vec) const
{
    if (s.length() < _ksize) {
        return;
    }
    for (unsigned int i = 0; i < s.length() - _ksize + 1; i++) {
        std::string sub = s.substr(i, _ksize);
        kmers_vec.push_back(sub);
    }
}


void Hashtable::get_kmer_hashes(const std::string &s,
                                std::vector<HashIntoType> &kmers_vec) const
{
    unique_ptr<KmerHashIterator> kmers = new_kmer_iterator(s);

    while(!kmers->done()) {
        HashIntoType kmer = kmers->next();
        kmers_vec.push_back(kmer);
    }
}


void Hashtable::get_kmer_hashes_as_hashset(const std::string &s,
        SeenSet& hashes) const
{
    unique_ptr<KmerHashIterator> kmers = new_kmer_iterator(s);

    while(!kmers->done()) {
        HashIntoType kmer = kmers->next();
        hashes.insert(kmer);
    }
}


void Hashtable::get_kmer_counts(const std::string &s,
                                std::vector<BoundedCounterType> &counts) const
{
    unique_ptr<KmerHashIterator> kmers = new_kmer_iterator(s);

    while(!kmers->done()) {
        HashIntoType kmer = kmers->next();
        BoundedCounterType c = this->get_count(kmer);
        counts.push_back(c);
    }
}

BoundedCounterType Hashtable::get_min_count(const std::string &s)
{
    unique_ptr<KmerHashIterator> kmers = new_kmer_iterator(s);

    BoundedCounterType min_count = MAX_KCOUNT;

    while(!kmers->done()) {
        HashIntoType kmer = kmers->next();

        BoundedCounterType count = this->get_count(kmer);

        if (this->get_count(kmer) < min_count) {
            min_count = count;
        }
    }
    return min_count;
}

BoundedCounterType Hashtable::get_max_count(const std::string &s)
{
    BoundedCounterType max_count = 0;

    unique_ptr<KmerHashIterator> kmers = new_kmer_iterator(s);

    while(!kmers->done()) {
        HashIntoType kmer = kmers->next();

        BoundedCounterType count = this->get_count(kmer);

        if (count > max_count) {
            max_count = count;
        }
    }
    return max_count;
}

uint64_t *
Hashtable::abundance_distribution(
    read_parsers::IParser * parser,
    Hashtable *          tracking)
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
            unique_ptr<KmerHashIterator> kmers = new_kmer_iterator(seq);

            while(!kmers->done()) {
                HashIntoType kmer = kmers->next();

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


uint64_t * Hashtable::abundance_distribution(
    std::string filename,
    Hashtable *  tracking)
{
    IParser* parser = IParser::get_parser(filename.c_str());

    uint64_t * distribution = abundance_distribution(parser, tracking);
    delete parser;
    return distribution;
}

unsigned long Hashtable::trim_on_abundance(
    std::string     seq,
    BoundedCounterType  min_abund)
const
{
    if (!check_and_normalize_read(seq)) {
        return 0;
    }

    unique_ptr<KmerHashIterator> kmers = new_kmer_iterator(seq);

    HashIntoType kmer;

    if (kmers->done()) {
        return 0;
    }
    kmer = kmers->next();

    if (kmers->done() || get_count(kmer) < min_abund) {
        return 0;
    }

    unsigned long i = _ksize;
    while (!kmers->done()) {
        kmer = kmers->next();

        if (get_count(kmer) < min_abund) {
            return i;
        }
        i++;
    }

    return seq.length();
}

unsigned long Hashtable::trim_below_abundance(
    std::string     seq,
    BoundedCounterType  max_abund)
const
{
    if (!check_and_normalize_read(seq)) {
        return 0;
    }

    unique_ptr<KmerHashIterator> kmers = new_kmer_iterator(seq);

    HashIntoType kmer;

    if (kmers->done()) {
        return 0;
    }
    kmer = kmers->next();

    if (kmers->done() || get_count(kmer) > max_abund) {
        return 0;
    }

    unsigned long i = _ksize;
    while (!kmers->done()) {
        kmer = kmers->next();

        if (get_count(kmer) > max_abund) {
            return i;
        }
        i++;
    }

    return seq.length();
}

std::vector<unsigned int> Hashtable::find_spectral_error_positions(
    std::string seq,
    BoundedCounterType max_abund)
const
{
    std::vector<unsigned int> posns;
    if (!check_and_normalize_read(seq)) {
        throw khmer_exception("invalid read");
    }

    unique_ptr<KmerHashIterator> kmers = new_kmer_iterator(seq);

    HashIntoType kmer = kmers->next();
    if (kmers->done()) {
        return posns;
    }

    // find the first trusted k-mer
    while (!kmers->done()) {
        if (get_count(kmer) > max_abund) {
            break;
        }
        kmer = kmers->next();
    }

    if (kmers->done()) {
        return posns;
    }

    // did we bypass some erroneous k-mers? call the last one.
    if (kmers->get_start_pos() > 0) {
        // if we are well past the first k, forget the whole thing (!? @CTB)
        if (kmers->get_start_pos() >= _ksize && 0) {
            return posns;
        }
        posns.push_back(kmers->get_start_pos() - 1);
    }

    while (!kmers->done()) {
        kmer = kmers->next();
        if (get_count(kmer) <= max_abund) { // error!
            posns.push_back(kmers->get_end_pos() - 1);

            // find next good
            while (!kmers->done()) {
                kmer = kmers->next();
                if (get_count(kmer) > max_abund) { // a good stretch again.
                    break;
                }
            }
        }
    }

    return posns;
}

class MurmurKmerHashIterator : public KmerHashIterator
{
    const char * _seq;
    const char _ksize;
    unsigned int index;
    unsigned int length;
public:
    MurmurKmerHashIterator(const char * seq, unsigned char k) :
        _seq(seq), _ksize(k) {
        index = 0;
        length = strlen(_seq);
    };

    HashIntoType first() { return next(); }

    HashIntoType next() {
        if (done()) {
            throw khmer_exception("past end of iterator");
        }

        uint64_t out[2];
        uint32_t seed = 0;
        HashIntoType h = 0, r = 0;

        std::string kmer;
        kmer.assign(_seq + index, _ksize);
        MurmurHash3_x64_128((void *)kmer.c_str(), _ksize, seed, &out);
        h = out[0];

        std::string kmer_rc = khmer::_revcomp(kmer);
        MurmurHash3_x64_128((void *)kmer_rc.c_str(), _ksize, seed, &out);
        r = out[0];

        index += 1;
        return h ^ r;
    }

    bool done() const {
        return (index + _ksize > length);
    }

    unsigned int get_start_pos() const { return index; }
    unsigned int get_end_pos() const { return index + _ksize; }
};

unique_ptr<KmerHashIterator> Counttable::new_kmer_iterator(const char * sp) const {
    KmerHashIterator * ki = new MurmurKmerHashIterator(sp, _ksize);
    return unique_ptr<KmerHashIterator>(ki);
}

// vim: set sts=2 sw=2:
