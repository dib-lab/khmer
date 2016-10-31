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

    KmerIterator kmers(sp, _ksize);

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
    KmerIterator kmers(s.c_str(), _ksize);
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

void Hashgraph::save_tagset(std::string outfilename)
{
    ofstream outfile(outfilename.c_str(), ios::binary);
    const size_t tagset_size = n_tags();
    unsigned int save_ksize = _ksize;

    HashIntoType * buf = new HashIntoType[tagset_size];

    outfile.write(SAVED_SIGNATURE, 4);
    unsigned char version = SAVED_FORMAT_VERSION;
    outfile.write((const char *) &version, 1);

    unsigned char ht_type = SAVED_TAGS;
    outfile.write((const char *) &ht_type, 1);

    outfile.write((const char *) &save_ksize, sizeof(save_ksize));
    outfile.write((const char *) &tagset_size, sizeof(tagset_size));
    outfile.write((const char *) &_tag_density, sizeof(_tag_density));

    unsigned int i = 0;
    for (SeenSet::iterator pi = all_tags.begin(); pi != all_tags.end();
            ++pi, i++) {
        buf[i] = *pi;
    }

    outfile.write((const char *) buf, sizeof(HashIntoType) * tagset_size);
    if (outfile.fail()) {
        delete[] buf;
        throw khmer_file_exception(strerror(errno));
    }
    outfile.close();

    delete[] buf;
}

void Hashgraph::load_tagset(std::string infilename, bool clear_tags)
{
    ifstream infile;

    // configure ifstream to raise exceptions for everything.
    infile.exceptions(std::ifstream::failbit | std::ifstream::badbit |
                      std::ifstream::eofbit);

    try {
        infile.open(infilename.c_str(), ios::binary);
    } catch (std::ifstream::failure &e) {
        std::string err;
        if (!(infile.is_open())) {
            err = "Cannot open tagset file: " + infilename;
        } else {
            err = "Unknown error in opening file: " + infilename;
        }
        throw khmer_file_exception(err);
    } catch (const std::exception &e) {
        // Catching std::exception is a stopgap for
        // https://gcc.gnu.org/bugzilla/show_bug.cgi?id=66145
        std::string err = "Unknown error opening file: " + infilename + " "
                          + strerror(errno);
        throw khmer_file_exception(err);
    }

    if (clear_tags) {
        all_tags.clear();
    }

    unsigned char version, ht_type;
    unsigned int save_ksize = 0;

    size_t tagset_size = 0;
    HashIntoType * buf = NULL;

    try {
        char signature[4];
        infile.read(signature, 4);
        infile.read((char *) &version, 1);
        infile.read((char *) &ht_type, 1);
        if (!(std::string(signature, 4) == SAVED_SIGNATURE)) {
            std::ostringstream err;
            err << "Incorrect file signature 0x";
            for(size_t i=0; i < 4; ++i) {
                err << std::hex << (int) signature[i];
            }
            err << " while reading tagset from " << infilename
                << "; should be " << SAVED_SIGNATURE;
            throw khmer_file_exception(err.str());
        } else if (!(version == SAVED_FORMAT_VERSION)) {
            std::ostringstream err;
            err << "Incorrect file format version " << (int) version
                << " while reading tagset from " << infilename
                << "; should be " << (int) SAVED_FORMAT_VERSION;
            throw khmer_file_exception(err.str());
        } else if (!(ht_type == SAVED_TAGS)) {
            std::ostringstream err;
            err << "Incorrect file format type " << (int) ht_type
                << " while reading tagset from " << infilename;
            throw khmer_file_exception(err.str());
        }

        infile.read((char *) &save_ksize, sizeof(save_ksize));
        if (!(save_ksize == _ksize)) {
            std::ostringstream err;
            err << "Incorrect k-mer size " << save_ksize
                << " while reading tagset from " << infilename;
            throw khmer_file_exception(err.str());
        }

        infile.read((char *) &tagset_size, sizeof(tagset_size));
        infile.read((char *) &_tag_density, sizeof(_tag_density));

        buf = new HashIntoType[tagset_size];

        infile.read((char *) buf, sizeof(HashIntoType) * tagset_size);

        for (unsigned int i = 0; i < tagset_size; i++) {
            all_tags.insert(buf[i]);
        }

        delete[] buf;
    } catch (std::ifstream::failure &e) {
        std::string err = "Error reading data from: " + infilename;
        if (buf != NULL) {
            delete[] buf;
        }
        throw khmer_file_exception(err);
        /* Yes, this is boneheaded. Unfortunately, there is a bug in gcc > 5
         * regarding the basic_ios::failure that makes it impossible to catch
         * with more specificty. So, we catch *all* exceptions after trying to
         * get the ifstream::failure, and assume it must have been the buggy one.
         * Unfortunately, this would also cause us to catch the
         * khmer_file_exceptions thrown above, so we catch them again first and
         * rethrow them :) If this is understandably irritating to you, please
         * bother the gcc devs at:
         *     https://gcc.gnu.org/bugzilla/show_bug.cgi?id=66145
         *
         * See also: http://media4.giphy.com/media/3o6UBpHgaXFDNAuttm/giphy.gif
         */
    } catch (khmer_file_exception &e) {
        throw e;
    } catch (const std::exception &e) {
        std::string err = "Unknown error opening file: " + infilename + " "
                          + strerror(errno);
        throw khmer_file_exception(err);
    }
}

void Hashgraph::consume_sequence_and_tag(const std::string& seq,
        unsigned long long& n_consumed,
        SeenSet * found_tags)
{
    bool kmer_tagged;

    KmerIterator kmers(seq.c_str(), _ksize);
    HashIntoType kmer;

    unsigned int since = _tag_density / 2 + 1;

    while(!kmers.done()) {
        kmer = kmers.next();
        bool is_new_kmer;

        // Set the bits for the kmer in the various hashtables,
        // and report on whether or not they had already been set.
        // This is probably better than first testing and then setting the bits,
        // as a failed test essentially results in doing the same amount of work
        // twice.
        if ((is_new_kmer = test_and_set_bits( kmer ))) {
            ++n_consumed;
        }

#if (1)
        if (is_new_kmer) {
            ++since;
        } else {
            ACQUIRE_ALL_TAGS_SPIN_LOCK
            kmer_tagged = set_contains(all_tags, kmer);
            RELEASE_ALL_TAGS_SPIN_LOCK
            if (kmer_tagged) {
                since = 1;
                if (found_tags) {
                    found_tags->insert(kmer);
                }
            } else {
                ++since;
            }
        }
#else
        if (!is_new_kmer && set_contains(all_tags, kmer)) {
            since = 1;
            if (found_tags) {
                found_tags->insert(kmer);
            }
        } else {
            since++;
        }
#endif

        if (since >= _tag_density) {
            ACQUIRE_ALL_TAGS_SPIN_LOCK
            all_tags.insert(kmer);
            RELEASE_ALL_TAGS_SPIN_LOCK
            if (found_tags) {
                found_tags->insert(kmer);
            }
            since = 1;
        }

    } // iteration over kmers

    if (since >= _tag_density/2 - 1) {
        ACQUIRE_ALL_TAGS_SPIN_LOCK
        all_tags.insert(kmer);	// insert the last k-mer, too.
        RELEASE_ALL_TAGS_SPIN_LOCK
        if (found_tags) {
            found_tags->insert(kmer);
        }
    }
}

//
// consume_fasta_and_tag: consume a FASTA file of reads, tagging reads every
//     so often.
//

// TODO? Inline in header.
void
Hashgraph::
consume_fasta_and_tag(
    std:: string const  &filename,
    unsigned int	      &total_reads, unsigned long long	&n_consumed
)
{
    IParser *	  parser =
        IParser::get_parser( filename );

    consume_fasta_and_tag(
        parser,
        total_reads, n_consumed
    );

    delete parser;
}

void
Hashgraph::
consume_fasta_and_tag(
    read_parsers:: IParser *  parser,
    unsigned int		    &total_reads,   unsigned long long	&n_consumed
)
{
    Read			  read;

    // TODO? Delete the following assignments.
    total_reads = 0;
    n_consumed = 0;

    // Iterate through the reads and consume their k-mers.
    while (!parser->is_complete( )) {

        try {
            read = parser->get_next_read( );
        } catch (NoMoreReadsAvailable &e) {
            // Bail out if this error is raised
            break;
        }

        if (check_and_normalize_read( read.sequence )) {
            unsigned long long this_n_consumed = 0;
            consume_sequence_and_tag( read.sequence, this_n_consumed );

            __sync_add_and_fetch( &n_consumed, this_n_consumed );
            __sync_add_and_fetch( &total_reads, 1 );
        }
    } // while reads left for parser

}

//
// divide_tags_into_subsets - take all of the tags in 'all_tags', and
//   divide them into subsets (based on starting tag) of <= given size.
//

void Hashgraph::divide_tags_into_subsets(unsigned int subset_size,
        SeenSet& divvy)
{
    unsigned int i = 0;

    for (SeenSet::const_iterator si = all_tags.begin(); si != all_tags.end();
            ++si) {
        if (i % subset_size == 0) {
            divvy.insert(*si);
            i = 0;
        }
        i++;
    }
}

//
// consume_partitioned_fasta: consume a FASTA file of reads
//

void Hashgraph::consume_partitioned_fasta(const std::string &filename,
        unsigned int &total_reads,
        unsigned long long &n_consumed)
{
    total_reads = 0;
    n_consumed = 0;

    IParser* parser = IParser::get_parser(filename.c_str());
    Read read;

    string seq = "";

    // reset the master subset partition
    delete partition;
    partition = new SubsetPartition(this);

    //
    // iterate through the FASTA file & consume the reads.
    //

    while(!parser->is_complete())  {
        try {
            read = parser->get_next_read();
        } catch (NoMoreReadsAvailable &exc) {
            break;
        }
        seq = read.sequence;

        if (check_and_normalize_read(seq)) {
            // First, figure out what the partition is (if non-zero), and save that.
            PartitionID p = _parse_partition_id(read.name);

            // Then consume the sequence
            n_consumed += consume_string(seq); // @CTB why are we doing this?

            // Next, compute the tag & set the partition, if nonzero
            HashIntoType kmer = hash_dna(seq.c_str());
            all_tags.insert(kmer);
            if (p > 0) {
                partition->set_partition_id(kmer, p);
            }
        }

        // reset the sequence info, increment read number
        total_reads++;
    }

    delete parser;
}

//////////////////////////////////////////////////////////////////////
// graph stuff

void Hashgraph::calc_connected_graph_size(Kmer start,
        unsigned long long& count,
        KmerSet& keeper,
        const unsigned long long threshold,
        bool break_on_circum)
const
{
    const BoundedCounterType val = get_count(start);

    if (val == 0) {
        return;
    }

    Traverser traverser(this);
    KmerQueue node_q;
    node_q.push(start);

    // Avoid high-circumference k-mers
    auto filter = [&] (Kmer& n) {
        return !(break_on_circum &&
                 traverser.degree(n) > 4);
    };

    while(!node_q.empty()) {
        Kmer node = node_q.front();
        node_q.pop();

        // have we already seen me? don't count; exit.
        if (set_contains(keeper, node)) {
            continue;
        }

        // is this in stop_tags?
        if (set_contains(stop_tags, node)) {
            continue;
        }

        // keep track of both seen kmers, and counts.
        keeper.insert(node);

        count += 1;

        // are we past the threshold? truncate search.
        if (threshold && count >= threshold) {
            return;
        }

        // otherwise, explore in all directions.
        traverser.traverse_right(node, node_q, filter);
        traverser.traverse_left(node, node_q, filter);
    }
}

unsigned int Hashgraph::kmer_degree(HashIntoType kmer_f, HashIntoType kmer_r)
{
    Traverser traverser(this);
    Kmer node = build_kmer(kmer_f, kmer_r);
    return traverser.degree(node);
}

unsigned int Hashgraph::kmer_degree(const char * kmer_s)
{
    Traverser traverser(this);
    Kmer node = build_kmer(kmer_s);
    return traverser.degree(node);
}

size_t Hashgraph::trim_on_stoptags(std::string seq) const
{
    if (!check_and_normalize_read(seq)) {
        return 0;
    }

    KmerIterator kmers(seq.c_str(), _ksize);

    size_t i = _ksize - 2;
    while (!kmers.done()) {
        HashIntoType kmer = kmers.next();
        if (set_contains(stop_tags, kmer)) {
            return i;
        }
        i++;
    }

    return seq.length();
}

unsigned int Hashgraph::traverse_from_kmer(Kmer start,
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

        if (set_contains(stop_tags, node)) {
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

void Hashgraph::load_stop_tags(std::string infilename, bool clear_tags)
{
    ifstream infile;

    // configure ifstream to raise exceptions for everything.
    infile.exceptions(std::ifstream::failbit | std::ifstream::badbit |
                      std::ifstream::eofbit);

    try {
        infile.open(infilename.c_str(), ios::binary);
    } catch (std::ifstream::failure &e) {
        std::string err;
        if (!(infile.is_open())) {
            err = "Cannot open stoptags file: " + infilename;
        } else {
            err = "Unknown error in opening file: " + infilename;
        }
        throw khmer_file_exception(err);
    } catch (const std::exception &e) {
        // Catching std::exception is a stopgap for
        // https://gcc.gnu.org/bugzilla/show_bug.cgi?id=66145
        std::string err = "Unknown error opening file: " + infilename + " "
                          + strerror(errno);
        throw khmer_file_exception(err);
    }

    if (clear_tags) {
        stop_tags.clear();
    }

    unsigned char version, ht_type;
    unsigned int save_ksize = 0;

    size_t tagset_size = 0;

    try {
        char signature[4];
        infile.read(signature, 4);
        infile.read((char *) &version, 1);
        infile.read((char *) &ht_type, 1);
        if (!(std::string(signature, 4) == SAVED_SIGNATURE)) {
            std::ostringstream err;
            err << "Incorrect file signature 0x";
            for(size_t i=0; i < 4; ++i) {
                err << std::hex << (int) signature[i];
            }
            err << " while reading stoptags from " << infilename
                << "; should be " << SAVED_SIGNATURE;
            throw khmer_file_exception(err.str());
        } else if (!(version == SAVED_FORMAT_VERSION)) {
            std::ostringstream err;
            err << "Incorrect file format version " << (int) version
                << " while reading stoptags from " << infilename
                << "; should be " << (int) SAVED_FORMAT_VERSION;
            throw khmer_file_exception(err.str());
        } else if (!(ht_type == SAVED_STOPTAGS)) {
            std::ostringstream err;
            err << "Incorrect file format type " << (int) ht_type
                << " while reading stoptags from " << infilename;
            throw khmer_file_exception(err.str());
        }

        infile.read((char *) &save_ksize, sizeof(save_ksize));
        if (!(save_ksize == _ksize)) {
            std::ostringstream err;
            err << "Incorrect k-mer size " << save_ksize
                << " while reading stoptags from " << infilename;
            throw khmer_file_exception(err.str());
        }
        infile.read((char *) &tagset_size, sizeof(tagset_size));

        HashIntoType * buf = new HashIntoType[tagset_size];

        infile.read((char *) buf, sizeof(HashIntoType) * tagset_size);

        for (unsigned int i = 0; i < tagset_size; i++) {
            stop_tags.insert(buf[i]);
        }
        delete[] buf;
    } catch (std::ifstream::failure &e) {
        std::string err = "Error reading stoptags from: " + infilename;
        throw khmer_file_exception(err);
    } catch (khmer_file_exception &e) {
        throw e;
    } catch (const std::exception &e) {
        // Catching std::exception is a stopgap for
        // https://gcc.gnu.org/bugzilla/show_bug.cgi?id=66145
        std::string err = "Unknown error opening file: " + infilename + " "
                          + strerror(errno);
        throw khmer_file_exception(err);
    }
}

void Hashgraph::save_stop_tags(std::string outfilename)
{
    ofstream outfile(outfilename.c_str(), ios::binary);
    size_t tagset_size = stop_tags.size();

    HashIntoType * buf = new HashIntoType[tagset_size];

    outfile.write(SAVED_SIGNATURE, 4);
    unsigned char version = SAVED_FORMAT_VERSION;
    outfile.write((const char *) &version, 1);

    unsigned char ht_type = SAVED_STOPTAGS;
    outfile.write((const char *) &ht_type, 1);

    unsigned int save_ksize = _ksize;
    outfile.write((const char *) &save_ksize, sizeof(save_ksize));
    outfile.write((const char *) &tagset_size, sizeof(tagset_size));

    unsigned int i = 0;
    for (SeenSet::iterator pi = stop_tags.begin(); pi != stop_tags.end();
            ++pi, i++) {
        buf[i] = *pi;
    }

    outfile.write((const char *) buf, sizeof(HashIntoType) * tagset_size);
    outfile.close();

    delete[] buf;
}

void Hashgraph::print_stop_tags(std::string infilename)
{
    ofstream printfile(infilename.c_str());

    unsigned int i = 0;
    for (SeenSet::iterator pi = stop_tags.begin(); pi != stop_tags.end();
            ++pi, i++) {
        std::string kmer = _revhash(*pi, _ksize);
        printfile << kmer << "\n";
    }

    printfile.close();
}

void Hashgraph::print_tagset(std::string infilename)
{
    ofstream printfile(infilename.c_str());

    unsigned int i = 0;
    for (SeenSet::iterator pi = all_tags.begin(); pi != all_tags.end();
            ++pi, i++) {
        std::string kmer = _revhash(*pi, _ksize);
        printfile << kmer << "\n";
    }

    printfile.close();
}

void Hashgraph::extract_unique_paths(std::string seq,
                                     unsigned int min_length,
                                     float min_unique_f,
                                     std::vector<std::string> &results)
{
    if (seq.size() < min_length) {
        return;
    }

    float max_seen = 1.0 - min_unique_f;

    min_length = min_length - _ksize + 1; // adjust for k-mer size.

    KmerIterator kmers(seq.c_str(), _ksize);

    std::deque<bool> seen_queue;
    unsigned int n_already_seen = 0;
    unsigned int n_kmers = 0;

    // first, put together an array for presence/absence of the k-mer
    // at each given position.
    while (!kmers.done()) {
        HashIntoType kmer = kmers.next();

        if (get_count(kmer)) {
            seen_queue.push_back(true);
            n_already_seen++;
        } else {
            seen_queue.push_back(false);
        }
        n_kmers++;
    }

    // next, run through this array with 'i'.

    unsigned int i = 0;
    while (i < n_kmers - min_length) {
        unsigned int seen_counter, j;

        // For each starting 'i', count the number of 'seen' k-mers in the
        // given window.

        // yes, inefficient n^2 algorithm.  sue me.
        for (seen_counter = 0, j = 0; j < min_length; j++) {
            if (seen_queue[i + j]) {
                seen_counter++;
            }
        }

        // If the fraction seen is small enough to be interesting, suggesting
        // that this, in fact, a "new" window -- extend until it isn't, and
        // then extract.

        if (!(j == min_length)) {
            throw khmer_exception();
        }
        if ( ((float)seen_counter / (float) j) <= max_seen) {
            unsigned int start = i;

            // extend the window until the end of the sequence...
            while ((start + min_length) < n_kmers) {
                if (seen_queue[start]) {
                    seen_counter--;
                }
                if (seen_queue[start + min_length]) {
                    seen_counter++;
                }
                start++;

                // ...or until we've seen too many of the k-mers.
                if (((float)seen_counter / (float) min_length) > max_seen) {
                    break;
                }
            }

            // adjust for ending point.
            if (start + min_length == n_kmers) {	// potentially decrement twice at end
                if (((float)seen_counter / (float) min_length) > max_seen) {
                    start--;
                }
                start--;
            } else {
                start -= 2;
            }

            // ...and now extract the relevant portion of the sequence, and adjust
            // starting pos'n.
            results.push_back(seq.substr(i, start + min_length + _ksize - i));

            i = start + min_length + 1;
        } else {
            i++;
        }
    }
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
    KmerIterator kmers(s.c_str(), _ksize);

    while(!kmers.done()) {
        HashIntoType kmer = kmers.next();
        kmers_vec.push_back(kmer);
    }
}


void Hashtable::get_kmer_hashes_as_hashset(const std::string &s,
        SeenSet& hashes) const
{
    KmerIterator kmers(s.c_str(), _ksize);

    while(!kmers.done()) {
        HashIntoType kmer = kmers.next();
        hashes.insert(kmer);
    }
}


void Hashtable::get_kmer_counts(const std::string &s,
                                std::vector<BoundedCounterType> &counts) const
{
    KmerIterator kmers(s.c_str(), _ksize);

    while(!kmers.done()) {
        HashIntoType kmer = kmers.next();
        BoundedCounterType c = this->get_count(kmer);
        counts.push_back(c);
    }
}

void Hashgraph::find_high_degree_nodes(const char * s,
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

unsigned int Hashgraph::traverse_linear_path(const Kmer seed_kmer,
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

// Starting from the given seed k-mer, assemble the maximal linear path in
// both directions.
//
// No guarantees on direction, of course - this may return the reverse
// complement of the input sequence.
//
// Note: as written, will ignore branches to the left and continue
// past them; this probably needs to be fixed.  For now, this means
// that assembling from two different directions may yield different
// results.

std::string Hashgraph::assemble_linear_path(const Kmer seed_kmer,
        const Hashtable * stop_bf)
const
{
    if (get_count(seed_kmer) == 0) {
        // If the seed k-mer is not in the de Bruijn graph, stop trying to make
        // something happen. It's not going to happen!
        return "";
    }

    std::string start_kmer = seed_kmer.get_string_rep(_ksize);
    std::string right = _assemble_right(start_kmer.c_str(), stop_bf);

    start_kmer = _revcomp(start_kmer);
    std::string left = _assemble_right(start_kmer.c_str(), stop_bf);

    left = left.substr(_ksize);
    return _revcomp(left) + right;
}

std::string Hashgraph::_assemble_right(const char * start_kmer,
                                       const Hashtable * stop_bf)
const
{
    const char bases[] = "ACGT";
    std::string kmer = start_kmer;
    std::string contig = kmer;

    // This loop extends the starting k-mer to the right as long as it can
    // do so unambiguously (or not at all).  This involves checking each
    // possible nucleotide suffix for presence; extension is continued until
    // either more than one such k-mer is present ('found2' is true), or no
    // such k-mer is present ('found' is false).

    while (1) {
        const char * base = &bases[0];
        bool found = false;
        char found_base;
        bool found2 = false;

        // check all four suffixes for presence.
        while(*base != 0) {
            std::string try_kmer = kmer.substr(1) + (char) *base;

            // a hit!
            if (this->get_count(try_kmer.c_str()) &&
                    (!stop_bf || !stop_bf->get_count(try_kmer.c_str()))) {
                if (found) {
                    found2 = true;
                    break;
                }
                found_base = (char) *base;
                found = true;
            }
            base++;
        }

        // exit condition: no suffix k-mer, or more than one.
        if (!found or found2) {
            break;
        }

        // extend assembly!
        contig += found_base;
        kmer = kmer.substr(1) + found_base;
        found = true;
    }
    return contig;
}

// vim: set sts=2 sw=2:
