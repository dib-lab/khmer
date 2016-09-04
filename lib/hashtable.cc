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

#include "counting.hh"
#include "hashtable.hh"
#include "khmer.hh"
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
    const char * sp = s.c_str();
    unsigned int n_consumed = 0;

    StringToHashIterator kmers(sp, _ksize);

    while(!kmers.done()) {
        HashIntoType kmer = kmers.next();

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
    StringToHashIterator kmers(s.c_str(), _ksize);
    unsigned int min_req = 0.5 + float(s.size() - _ksize + 1) / 2;
    unsigned int num_cutoff_kmers = 0;

    // first loop:
    // accumulate at least min_req worth of counts before checking to see
    // if we have enough high-abundance k-mers to indicate success.
    for (unsigned int i = 0; i < min_req; ++i) {
        HashIntoType kmer = kmers.next();
        if (this->get_count(kmer) >= cutoff) {
            ++num_cutoff_kmers;
        }
    }

    // second loop: now check to see if we pass the threshold for each k-mer.
    if (num_cutoff_kmers >= min_req) {
        return true;
    }
    while(!kmers.done()) {
        HashIntoType kmer = kmers.next();
        if (this->get_count(kmer) >= cutoff) {
            ++num_cutoff_kmers;
            if (num_cutoff_kmers >= min_req) {
                return true;
            }
        }
    }
    return false;
}

//////////////////////////////////////////////////////////////////////
// graph stuff

unsigned int Hashtable::kmer_degree(const char * kmer_s)
{
    Traverser traverser(this);
    Kmer node = build_kmer(kmer_s);
    return traverser.degree(node);
}

unsigned int Hashtable::traverse_from_kmer(Kmer start,
        unsigned int radius,
        KmerSet &keeper,
        unsigned int max_count)
const
{

    Traverser traverser(this);
    KmerQueue node_q;
    std::queue<unsigned int> breadth_q;
    unsigned int cur_breadth = 0;
    unsigned int total = 0;
    unsigned int nfound = 0;

    auto filter = [&] (Kmer& n) {
        return !set_contains(keeper, n);
    };

    node_q.push(start);
    breadth_q.push(0);

    while(!node_q.empty()) {
        Kmer node = node_q.front();
        node_q.pop();

        unsigned int breadth = breadth_q.front();
        breadth_q.pop();

        if (breadth > radius) {
            break;
        }

        if (max_count && total > max_count) {
            break;
        }

        if (set_contains(keeper, node)) {
            continue;
        }

        // keep track of seen kmers
        keeper.insert(node);
        total++;

        if (!(breadth >= cur_breadth)) { // keep track of watermark, for debugging.
            throw khmer_exception();
        }
        if (breadth > cur_breadth) {
            cur_breadth = breadth;
        }

        nfound = traverser.traverse_right(node, node_q, filter);
        for (unsigned int i = 0; i<nfound; ++i) {
            breadth_q.push(breadth + 1);
        }

        nfound = traverser.traverse_left(node, node_q, filter);
        for (unsigned int i = 0; i<nfound; ++i) {
            breadth_q.push(breadth + 1);
        }
    }

    return total;
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
    StringToHashIterator kmers(s.c_str(), _ksize);

    while(!kmers.done()) {
        HashIntoType kmer = kmers.next();
        kmers_vec.push_back(kmer);
    }
}


void Hashtable::get_kmer_hashes_as_hashset(const std::string &s,
        SeenSet& hashes) const
{
    StringToHashIterator kmers(s.c_str(), _ksize);

    while(!kmers.done()) {
        HashIntoType kmer = kmers.next();
        hashes.insert(kmer);
    }
}


void Hashtable::get_kmer_counts(const std::string &s,
                                std::vector<BoundedCounterType> &counts) const
{
    StringToHashIterator kmers(s.c_str(), _ksize);

    while(!kmers.done()) {
        HashIntoType kmer = kmers.next();
        BoundedCounterType c = this->get_count(kmer);
        counts.push_back(c);
    }
}

void Hashtable::find_high_degree_nodes(const char * s,
                                       SeenSet& high_degree_nodes)
const
{
    Traverser traverser(this);
    KmerIterator kmers(s, _ksize);

    unsigned long n = 0;
    while(!kmers.done()) {
        n++;
        if (n % 10000 == 0) {
            std::cout << "... find_high_degree_nodes: " << n << "\n";
            std::cout << std::flush;
        }

        Kmer kmer = kmers.next();
        if ((traverser.degree(kmer)) > 2) {
            high_degree_nodes.insert(kmer);
        }
    }
}


unsigned int Hashtable::traverse_linear_path(const Kmer seed_kmer,
        SeenSet &adjacencies,
        SeenSet &visited, Hashtable &bf,
        SeenSet &high_degree_nodes)
const
{
    unsigned int size = 0;

    Traverser traverser(this);

    // if this k-mer is in the Bloom filter, truncate search.
    // This prevents paths from being traversed in two directions.
    if (bf.get_count(seed_kmer)) {
        return 0;
    }

    std::vector<Kmer> to_be_visited;
    to_be_visited.push_back(seed_kmer);

    while (to_be_visited.size()) {
        Kmer kmer = to_be_visited.back();
        to_be_visited.pop_back();

        visited.insert(kmer);
        size += 1;

        KmerQueue node_q;
        traverser.traverse(kmer, node_q);

        while (node_q.size()) {
            Kmer node = node_q.front();
            node_q.pop();

            if (set_contains(high_degree_nodes, node)) {
                // if there are any adjacent high degree nodes, record;
                adjacencies.insert(node);
                // also, add this to the stop Bloom filter.
                bf.count(kmer);
            } else if (set_contains(visited, node)) {
                // do nothing - already visited
                ;
            } else {
                to_be_visited.push_back(node);
            }
        }
    }
    return size;
}

// vim: set sts=2 sw=2:
