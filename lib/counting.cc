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

void CountingHash::save(std::string outfilename)
{
    CountingHashFile::save(outfilename, *this);
}

void CountingHash::load(std::string infilename)
{
    CountingHashFile::load(infilename, *this);
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


void CountingHashFile::load(
    const std::string   &infilename,
    CountingHash    &ht)
{
    std::string filename(infilename);
    size_t found = filename.find_last_of(".");
    std::string type = filename.substr(found + 1);

    if (type == "gz") {
        CountingHashGzFileReader(filename, ht);
    } else {
        CountingHashFileReader(filename, ht);
    }
}


void CountingHashFile::save(
    const std::string   &outfilename,
    const CountingHash  &ht)
{
    std::string filename(outfilename);
    size_t found = filename.find_last_of(".");
    std::string type = filename.substr(found + 1);

    if (type == "gz") {
        CountingHashGzFileWriter(filename, ht);
    } else {
        CountingHashFileWriter(filename, ht);
    }
}


CountingHashFileReader::CountingHashFileReader(
    const std::string   &infilename,
    CountingHash    &ht)
{
    ifstream infile;
    // configure ifstream to raise exceptions for everything.
    infile.exceptions(std::ifstream::failbit | std::ifstream::badbit |
                      std::ifstream::eofbit);

    try {
        infile.open(infilename.c_str(), ios::binary);
    } catch (std::ifstream::failure &e) {
        std::string err;
        if (!infile.is_open()) {
            err = "Cannot open k-mer count file: " + infilename;
        } else {
            err = "Unknown error in opening file: " + infilename;
        }
        throw khmer_file_exception(err + " " + strerror(errno));
    } catch (const std::exception &e) {
        std::string err = "Unknown error opening file: " + infilename + " "
                  + strerror(errno) + "; please poke the gcc devs at "
                  "https://gcc.gnu.org/bugzilla/show_bug.cgi?id=66145";
        throw khmer_file_exception(err);
    }

    if (ht._counts) {
        for (unsigned int i = 0; i < ht._n_tables; i++) {
            delete[] ht._counts[i];
            ht._counts[i] = NULL;
        }
        delete[] ht._counts;
        ht._counts = NULL;
    }
    ht._tablesizes.clear();

    try {
        unsigned int save_ksize = 0;
        unsigned char save_n_tables = 0;
        unsigned long long save_tablesize = 0;
        unsigned long long save_occupied_bins = 0;
        char signature [4];
        unsigned char version = 0, ht_type = 0, use_bigcount = 0;

        infile.read(signature, 4);
        infile.read((char *) &version, 1);
        infile.read((char *) &ht_type, 1);
        if (!(std::string(signature, 4) == SAVED_SIGNATURE)) {
            std::ostringstream err;
            err << "Does not start with signature for a khmer file: 0x";
            for(size_t i=0; i < 4; ++i) {
                err << std::hex << (int) signature[i];
            }
            err << " Should be: " << SAVED_SIGNATURE;
            throw khmer_file_exception(err.str());
        } else if (!(version == SAVED_FORMAT_VERSION)) {
            std::ostringstream err;
            err << "Incorrect file format version " << (int) version
                << " while reading k-mer count file from " << infilename
                << "; should be " << (int) SAVED_FORMAT_VERSION;
            throw khmer_file_exception(err.str());
        } else if (!(ht_type == SAVED_COUNTING_HT)) {
            std::ostringstream err;
            err << "Incorrect file format type " << (int) ht_type
                << " while reading k-mer count file from " << infilename;
            throw khmer_file_exception(err.str());
        }

        infile.read((char *) &use_bigcount, 1);
        infile.read((char *) &save_ksize, sizeof(save_ksize));
        infile.read((char *) &save_n_tables, sizeof(save_n_tables));
        infile.read((char *) &save_occupied_bins, sizeof(save_occupied_bins));

        ht._ksize = (WordLength) save_ksize;
        ht._n_tables = (unsigned int) save_n_tables;
        ht._occupied_bins = save_occupied_bins;
        ht._init_bitstuff();

        ht._use_bigcount = use_bigcount;

        ht._counts = new Byte*[ht._n_tables];
        for (unsigned int i = 0; i < ht._n_tables; i++) {
            ht._counts[i] = NULL;
        }

        for (unsigned int i = 0; i < ht._n_tables; i++) {
            uint64_t tablesize;

            infile.read((char *) &save_tablesize, sizeof(save_tablesize));

            tablesize = save_tablesize;
            ht._tablesizes.push_back(tablesize);

            ht._counts[i] = new Byte[tablesize];

            unsigned long long loaded = 0;
            while (loaded != tablesize) {
                infile.read((char *) ht._counts[i], tablesize - loaded);
                loaded += infile.gcount();
            }
        }

        uint64_t n_counts = 0;
        infile.read((char *) &n_counts, sizeof(n_counts));

        if (n_counts) {
            ht._bigcounts.clear();

            HashIntoType kmer;
            BoundedCounterType count;

            for (uint64_t n = 0; n < n_counts; n++) {
                infile.read((char *) &kmer, sizeof(kmer));
                infile.read((char *) &count, sizeof(count));
                ht._bigcounts[kmer] = count;
            }
        }

        infile.close();
    } catch (std::ifstream::failure &e) {
        std::string err;
        if (infile.eof()) {
            err = "Unexpected end of k-mer count file: " + infilename;
        } else {
            err = "Error reading from k-mer count file: " + infilename + " "
                  + strerror(errno);
        }
        throw khmer_file_exception(err);
    } catch (const std::exception &e) {
        std::string err = "Error reading from k-mer count file: " + infilename + " "
                  + strerror(errno);
        throw khmer_file_exception(err);
    }
}

CountingHashGzFileReader::CountingHashGzFileReader(
    const std::string   &infilename,
    CountingHash    &ht)
{
    gzFile infile = gzopen(infilename.c_str(), "rb");
    if (infile == Z_NULL) {
        std::string err = "Cannot open k-mer count file: " + infilename;
        throw khmer_file_exception(err);
    }

    if (ht._counts) {
        for (unsigned int i = 0; i < ht._n_tables; i++) {
            delete[] ht._counts[i];
            ht._counts[i] = NULL;
        }
        delete[] ht._counts;
        ht._counts = NULL;
    }
    ht._tablesizes.clear();

    unsigned int save_ksize = 0;
    unsigned char save_n_tables = 0;
    unsigned long long save_tablesize = 0;
    unsigned long long save_occupied_bins = 0;
    char signature [4];
    unsigned char version, ht_type, use_bigcount;

    int read_s = gzread(infile, signature, 4);
    int read_v = gzread(infile, (char *) &version, 1);
    int read_t = gzread(infile, (char *) &ht_type, 1);
    if (read_s <= 0 || read_v <= 0 || read_t <= 0) {
        std::string err = "K-mer count file read error: " + infilename + " "
                          + strerror(errno);
        gzclose(infile);
        throw khmer_file_exception(err);
    } else if (!(std::string(signature, 4) == SAVED_SIGNATURE)) {
        std::ostringstream err;
        err << "Does not start with signature for a khmer " <<
            "file: " << signature << " Should be: " <<
            SAVED_SIGNATURE;
        throw khmer_file_exception(err.str());
    } else if (!(version == SAVED_FORMAT_VERSION)
               || !(ht_type == SAVED_COUNTING_HT)) {
        if (!(version == SAVED_FORMAT_VERSION)) {
            std::ostringstream err;
            err << "Incorrect file format version " << (int) version
                << " while reading k-mer count file from " << infilename
                << "; should be " << (int) SAVED_FORMAT_VERSION;
            gzclose(infile);
            throw khmer_file_exception(err.str());
        } else if (!(ht_type == SAVED_COUNTING_HT)) {
            std::ostringstream err;
            err << "Incorrect file format type " << (int) ht_type
                << " while reading k-mer count file from " << infilename;
            gzclose(infile);
            throw khmer_file_exception(err.str());
        }
    }

    int read_b = gzread(infile, (char *) &use_bigcount, 1);
    int read_k = gzread(infile, (char *) &save_ksize, sizeof(save_ksize));
    int read_nt = gzread(infile, (char *) &save_n_tables,
                         sizeof(save_n_tables));
    int read_ob = gzread(infile, (char *) &save_occupied_bins,
                         sizeof(save_occupied_bins));

    if (read_b <= 0 || read_k <= 0 || read_nt <= 0 || read_ob <= 0) {
        std::string err = "K-mer count file header read error: " + infilename
                          + " " + strerror(errno);
        gzclose(infile);
        throw khmer_file_exception(err);
    }

    ht._ksize = (WordLength) save_ksize;
    ht._occupied_bins = save_occupied_bins;
    ht._n_tables = (unsigned int) save_n_tables;
    ht._init_bitstuff();

    ht._use_bigcount = use_bigcount;

    ht._counts = new Byte*[ht._n_tables];
    for (unsigned int i = 0; i < ht._n_tables; i++) {
        uint64_t tablesize;

        read_b = gzread(infile, (char *) &save_tablesize,
                        sizeof(save_tablesize));

        if (read_b <= 0) {
            std::string gzerr = gzerror(infile, &read_b);
            std::string err = "K-mer count file header read error: "
                              + infilename;
            if (read_b == Z_ERRNO) {
                err = err + " " + strerror(errno);
            } else {
                err = err + " " + gzerr;
            }
            gzclose(infile);
            throw khmer_file_exception(err);
        }

        tablesize = save_tablesize;
        ht._tablesizes.push_back(tablesize);

        ht._counts[i] = new Byte[tablesize];

        uint64_t loaded = 0;
        while (loaded != tablesize) {
            unsigned long long  to_read_ll = tablesize - loaded;
            unsigned int        to_read_int;
            // Zlib can only read chunks of at most INT_MAX bytes.
            if (to_read_ll > INT_MAX) {
                to_read_int = INT_MAX;
            } else {
                to_read_int = (unsigned int) to_read_ll;
            }
            read_b = gzread(infile, (char *) ht._counts[i], to_read_int);

            if (read_b <= 0) {
                std::string gzerr = gzerror(infile, &read_b);
                std::string err = "K-mer count file read error: " + infilename;
                if (read_b == Z_ERRNO) {
                    err = err + " " + strerror(errno);
                } else {
                    err = err + " " + gzerr;
                }
                gzclose(infile);
                throw khmer_file_exception(err);
            }

            loaded += read_b;
        }
    }

    uint64_t n_counts = 0;
    read_b = gzread(infile, (char *) &n_counts, sizeof(n_counts));
    if (read_b <= 0) {
        std::string gzerr = gzerror(infile, &read_b);
        std::string err = "K-mer count header read error: " + infilename;
        if (read_b == Z_ERRNO) {
            err = err + " " + strerror(errno);
        } else {
            err = err + " " + gzerr;
        }
        gzclose(infile);
        throw khmer_file_exception(err);
    }

    if (n_counts) {
        ht._bigcounts.clear();

        HashIntoType kmer;
        BoundedCounterType count;

        for (uint64_t n = 0; n < n_counts; n++) {
            int read_k = gzread(infile, (char *) &kmer, sizeof(kmer));
            int read_c = gzread(infile, (char *) &count, sizeof(count));

            if (read_k <= 0 || read_c <= 0) {
                std::string gzerr = gzerror(infile, &read_b);
                std::string err = "K-mer count read error: " + infilename;
                if (read_b == Z_ERRNO) {
                    err = err + " " + strerror(errno);
                } else {
                    err = err + " " + gzerr;
                }
                gzclose(infile);
                throw khmer_file_exception(err);
            }

            ht._bigcounts[kmer] = count;
        }
    }

    gzclose(infile);
}

CountingHashFileWriter::CountingHashFileWriter(
    const std::string   &outfilename,
    const CountingHash  &ht)
{
    if (!ht._counts[0]) {
        throw khmer_exception();
    }

    unsigned int save_ksize = ht._ksize;
    unsigned char save_n_tables = ht._n_tables;
    unsigned long long save_tablesize;
    unsigned long long save_occupied_bins = ht._occupied_bins;

    ofstream outfile(outfilename.c_str(), ios::binary);

    outfile.write(SAVED_SIGNATURE, 4);
    unsigned char version = SAVED_FORMAT_VERSION;
    outfile.write((const char *) &version, 1);

    unsigned char ht_type = SAVED_COUNTING_HT;
    outfile.write((const char *) &ht_type, 1);

    unsigned char use_bigcount = 0;
    if (ht._use_bigcount) {
        use_bigcount = 1;
    }
    outfile.write((const char *) &use_bigcount, 1);

    outfile.write((const char *) &save_ksize, sizeof(save_ksize));
    outfile.write((const char *) &save_n_tables, sizeof(save_n_tables));
    outfile.write((const char *) &save_occupied_bins,
                  sizeof(save_occupied_bins));

    for (unsigned int i = 0; i < save_n_tables; i++) {
        save_tablesize = ht._tablesizes[i];

        outfile.write((const char *) &save_tablesize, sizeof(save_tablesize));
        outfile.write((const char *) ht._counts[i], save_tablesize);
    }

    uint64_t n_counts = ht._bigcounts.size();
    outfile.write((const char *) &n_counts, sizeof(n_counts));

    if (n_counts) {
        KmerCountMap::const_iterator it = ht._bigcounts.begin();

        for (; it != ht._bigcounts.end(); ++it) {
            outfile.write((const char *) &it->first, sizeof(it->first));
            outfile.write((const char *) &it->second, sizeof(it->second));
        }
    }
    if (outfile.fail()) {
        throw khmer_file_exception(strerror(errno));
    }
    outfile.close();
}

CountingHashGzFileWriter::CountingHashGzFileWriter(
    const std::string   &outfilename,
    const CountingHash  &ht)
{
    if (!ht._counts[0]) {
        throw khmer_exception();
    }

    int errnum = 0;
    unsigned int save_ksize = ht._ksize;
    unsigned char save_n_tables = ht._n_tables;
    unsigned long long save_tablesize;
    unsigned long long save_occupied_bins = ht._occupied_bins;

    gzFile outfile = gzopen(outfilename.c_str(), "wb");
    if (outfile == NULL) {
        const char * error = gzerror(outfile, &errnum);
        if (errnum == Z_ERRNO) {
            throw khmer_file_exception(strerror(errno));
        } else {
            throw khmer_file_exception(error);
        }
    }

    gzwrite(outfile, SAVED_SIGNATURE, 4);
    unsigned char version = SAVED_FORMAT_VERSION;
    gzwrite(outfile, (const char *) &version, 1);

    unsigned char ht_type = SAVED_COUNTING_HT;
    gzwrite(outfile, (const char *) &ht_type, 1);

    unsigned char use_bigcount = 0;
    if (ht._use_bigcount) {
        use_bigcount = 1;
    }
    gzwrite(outfile, (const char *) &use_bigcount, 1);

    gzwrite(outfile, (const char *) &save_ksize, sizeof(save_ksize));
    gzwrite(outfile, (const char *) &save_n_tables, sizeof(save_n_tables));
    gzwrite(outfile, (const char *) &save_occupied_bins,
            sizeof(save_occupied_bins));

    for (unsigned int i = 0; i < save_n_tables; i++) {
        save_tablesize = ht._tablesizes[i];

        gzwrite(outfile, (const char *) &save_tablesize,
                sizeof(save_tablesize));
        unsigned long long written = 0;
        while (written != save_tablesize) {
            unsigned long long  to_write_ll = save_tablesize - written;
            unsigned int        to_write_int;
            int                 gz_result;
            // Zlib can only write chunks of at most INT_MAX bytes.
            if (to_write_ll > INT_MAX) {
                to_write_int = INT_MAX;
            } else {
                to_write_int = (unsigned int) to_write_ll;
            }
            gz_result = gzwrite(outfile, (const char *) ht._counts[i],
                                to_write_int);
            // Zlib returns 0 on error
            if (gz_result == 0) {
                int errcode = 0;
                const char *err_msg;
                std::ostringstream msg;

                msg << "gzwrite failed while writing counting hash: ";
                // Get zlib error
                err_msg = gzerror(outfile, &errcode);
                if (errcode != Z_ERRNO) {
                    // Zlib error, not stdlib
                    msg << err_msg;
                    gzclearerr(outfile);
                } else {
                    // stdlib error
                    msg << strerror(errno);
                }
                gzclose(outfile);
                throw khmer_file_exception(msg.str());
            }
            written += gz_result;
        }
    }

    uint64_t n_counts = ht._bigcounts.size();
    gzwrite(outfile, (const char *) &n_counts, sizeof(n_counts));

    if (n_counts) {
        KmerCountMap::const_iterator it = ht._bigcounts.begin();

        for (; it != ht._bigcounts.end(); ++it) {
            gzwrite(outfile, (const char *) &it->first, sizeof(it->first));
            gzwrite(outfile, (const char *) &it->second, sizeof(it->second));
        }
    }
    const char * error = gzerror(outfile, &errnum);
    if (errnum == Z_ERRNO) {
        throw khmer_file_exception(strerror(errno));
    } else if (errnum != Z_OK) {
        throw khmer_file_exception(error);
    }
    gzclose(outfile);
}

/* vim: set ft=cpp ts=8 sts=4 sw=4 et tw=79 */
