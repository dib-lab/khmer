//
// This file is part of khmer, https://github.com/dib-lab/khmer/, and is
// Copyright (C) Michigan State University, 2009-2015. It is licensed under
// the three-clause BSD license; see LICENSE.
// Contact: khmer-project@idyll.org
//

#include "hashtable.hh"
#include "counting.hh"
#include "hashbits.hh"
#include "read_parsers.hh"

#include "zlib.h"
#include <math.h>
#include <algorithm>
#include <sstream>
#include <errno.h>

using namespace std;
using namespace khmer;
using namespace khmer:: read_parsers;

///
/// output_fasta_kmer_pos_freq: outputs the kmer frequencies for each read
///

void CountingHash::output_fasta_kmer_pos_freq(
    const std::string &inputfile,
    const std::string &outputfile)
{
    IParser* parser = IParser::get_parser(inputfile.c_str());
    ofstream outfile;
    outfile.open(outputfile.c_str());
    string seq;
    Read read;

    while(!parser->is_complete()) {
        read = parser->get_next_read();
        seq = read.sequence;

        long numPos = seq.length() - _ksize + 1;

        for (long i = 0; i < numPos; i++)  {
            string kmer = seq.substr(i, _ksize);
            outfile << (int)get_count(kmer.c_str()) << " ";
        }
        outfile << endl;
    }

    delete parser;
    if (outfile.fail()) {
        throw khmer_file_exception(strerror(errno));
    }

    outfile.close();
}

const HashIntoType CountingHash::n_unique_kmers() const
{
    return _n_unique_kmers;
}

BoundedCounterType CountingHash::get_min_count(const std::string &s)
{
    KMerIterator kmers(s.c_str(), _ksize);

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
    KMerIterator kmers(s.c_str(), _ksize);

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

HashIntoType *
CountingHash::abundance_distribution(
    read_parsers::IParser * parser,
    Hashbits *          tracking)
{
    HashIntoType * dist = new HashIntoType[MAX_BIGCOUNT + 1];
    HashIntoType i;

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

    try {
        while(!parser->is_complete()) {
            read = parser->get_next_read();
            seq = read.sequence;

            if (check_and_normalize_read(seq)) {
                KMerIterator kmers(seq.c_str(), _ksize);

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
    } catch (NoMoreReadsAvailable) {
    }
    return dist;
}


HashIntoType * CountingHash::abundance_distribution(
    std::string filename,
    Hashbits *  tracking)
{
    IParser* parser = IParser::get_parser(filename.c_str());

    HashIntoType * distribution = abundance_distribution(parser, tracking);
    delete parser;
    return distribution;
}

HashIntoType * CountingHash::fasta_count_kmers_by_position(
    const std::string   &inputfile,
    const unsigned int  max_read_len,
    BoundedCounterType  limit_by_count,
    CallbackFn      callback,
    void *      callback_data)
{
    unsigned long long *counts = new unsigned long long[max_read_len];

    for (unsigned int i = 0; i < max_read_len; i++) {
        counts[i] = 0;
    }

    Read read;
    IParser* parser = IParser::get_parser(inputfile.c_str());
    string name;
    string seq;
    unsigned long long read_num = 0;

    while(!parser->is_complete()) {
        read = parser->get_next_read();

        seq = read.sequence;
        bool valid_read = check_and_normalize_read(seq);

        if (valid_read) {
            for (unsigned int i = 0; i < seq.length() - _ksize + 1; i++) {
                string kmer = seq.substr(i, i + _ksize);
                BoundedCounterType n = get_count(kmer.c_str());

                if (limit_by_count == 0 || n == limit_by_count) {
                    if (i < max_read_len) {
                        counts[i]++;
                    }
                }
            }
        }

        name.clear();
        seq.clear();

        read_num += 1;

        // run callback, if specified
        if (read_num % CALLBACK_PERIOD == 0 && callback) {
            try {
                callback("fasta_file_count_kmers_by_position", callback_data,
                         read_num, 0);
            } catch (...) {
                throw;
            }
        }
    } // while reads

    delete parser;

    return counts;
}

void CountingHash::fasta_dump_kmers_by_abundance(
    const std::string   &inputfile,
    BoundedCounterType  limit_by_count,
    CallbackFn      callback,
    void *      callback_data)
{
    Read read;
    IParser* parser = IParser::get_parser(inputfile.c_str());
    string name;
    string seq;
    unsigned long long read_num = 0;

    while(!parser->is_complete()) {
        read = parser->get_next_read();
        bool valid_read = check_and_normalize_read(seq);
        seq = read.sequence;

        if (valid_read) {
            for (unsigned int i = 0; i < seq.length() - _ksize + 1; i++) {
                string kmer = seq.substr(i, i + _ksize);
                BoundedCounterType n = get_count(kmer.c_str());
                char ss[_ksize + 1];
                strncpy(ss, kmer.c_str(), _ksize);
                ss[_ksize] = 0;

                if (n == limit_by_count) {
                    cout << ss << endl;
                }
            }
        }

        name.clear();
        seq.clear();

        read_num += 1;

        // run callback, if specified
        if (read_num % CALLBACK_PERIOD == 0 && callback) {
            try {
                callback("fasta_file_dump_kmers_by_abundance", callback_data,
                         read_num, 0);
            } catch (...) {
                throw;
            }
        }
    } // while reads

    delete parser;
}

void CountingHash::save(std::string outfilename)
{
    CountingHashFile::save(outfilename, *this);
}

void CountingHash::load(std::string infilename)
{
    CountingHashFile::load(infilename, *this);
}

void CountingHash::get_kadian_count(
    const std::string   &s,
    BoundedCounterType  &kadian,
    unsigned int    nk)
{
    std::vector<BoundedCounterType> counts;
    KMerIterator kmers(s.c_str(), _ksize);

    while(!kmers.done()) {
        HashIntoType kmer = kmers.next();
        BoundedCounterType count = this->get_count(kmer);
        counts.push_back(count);
    }

    if (!counts.size()) {
        throw khmer_exception();
    }
    unsigned int kpos = nk * _ksize;

    if (counts.size() < kpos) {
        kadian = 0;

        return;
    }

    sort(counts.begin(), counts.end());
    kadian = counts[kpos - 1];

#if 0
    std::cout << "k " << kpos << ": ";
    for (unsigned int i = 0; i < counts.size(); i++) {
        std::cout << i << "-" << counts[i] << " ";
    }
    std::cout << "\n";
#endif // 0
}

unsigned long CountingHash::trim_on_abundance(
    std::string     seq,
    BoundedCounterType  min_abund)
const
{
    if (!check_and_normalize_read(seq)) {
        return 0;
    }

    KMerIterator kmers(seq.c_str(), _ksize);

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

    KMerIterator kmers(seq.c_str(), _ksize);

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

    KMerIterator kmers(seq.c_str(), _ksize);

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
        unsigned char version = 0, ht_type = 0, use_bigcount = 0;

        infile.read((char *) &version, 1);
        infile.read((char *) &ht_type, 1);
        if (!(version == SAVED_FORMAT_VERSION)) {
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

        ht._ksize = (WordLength) save_ksize;
        ht._n_tables = (unsigned int) save_n_tables;
        ht._init_bitstuff();

        ht._use_bigcount = use_bigcount;

        ht._counts = new Byte*[ht._n_tables];
        for (unsigned int i = 0; i < ht._n_tables; i++) {
            ht._counts[i] = NULL;
        }

        for (unsigned int i = 0; i < ht._n_tables; i++) {
            HashIntoType tablesize;

            infile.read((char *) &save_tablesize, sizeof(save_tablesize));

            tablesize = (HashIntoType) save_tablesize;
            ht._tablesizes.push_back(tablesize);

            ht._counts[i] = new Byte[tablesize];

            unsigned long long loaded = 0;
            while (loaded != tablesize) {
                infile.read((char *) ht._counts[i], tablesize - loaded);
                loaded += infile.gcount();
            }
        }

        HashIntoType n_counts = 0;
        infile.read((char *) &n_counts, sizeof(n_counts));

        if (n_counts) {
            ht._bigcounts.clear();

            HashIntoType kmer;
            BoundedCounterType count;

            for (HashIntoType n = 0; n < n_counts; n++) {
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
    unsigned char version, ht_type, use_bigcount;

    int read_v = gzread(infile, (char *) &version, 1);
    int read_t = gzread(infile, (char *) &ht_type, 1);

    if (read_v <= 0 || read_t <= 0) {
        std::string err = "K-mer count file read error: " + infilename + " "
                          + strerror(errno);
        gzclose(infile);
        throw khmer_file_exception(err);
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

    if (read_b <= 0 || read_k <= 0 || read_nt <= 0) {
        std::string err = "K-mer count file header read error: " + infilename
                          + " " + strerror(errno);
        gzclose(infile);
        throw khmer_file_exception(err);
    }

    ht._ksize = (WordLength) save_ksize;
    ht._n_tables = (unsigned int) save_n_tables;
    ht._init_bitstuff();

    ht._use_bigcount = use_bigcount;

    ht._counts = new Byte*[ht._n_tables];
    for (unsigned int i = 0; i < ht._n_tables; i++) {
        HashIntoType tablesize;

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

        tablesize = (HashIntoType) save_tablesize;
        ht._tablesizes.push_back(tablesize);

        ht._counts[i] = new Byte[tablesize];

        HashIntoType loaded = 0;
        while (loaded != tablesize) {
            read_b = gzread(infile, (char *) ht._counts[i],
                            (unsigned) (tablesize - loaded));

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

    HashIntoType n_counts = 0;
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

        for (HashIntoType n = 0; n < n_counts; n++) {
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

    ofstream outfile(outfilename.c_str(), ios::binary);

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

    for (unsigned int i = 0; i < save_n_tables; i++) {
        save_tablesize = ht._tablesizes[i];

        outfile.write((const char *) &save_tablesize, sizeof(save_tablesize));
        outfile.write((const char *) ht._counts[i], save_tablesize);
    }

    HashIntoType n_counts = ht._bigcounts.size();
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

    gzFile outfile = gzopen(outfilename.c_str(), "wb");
    if (outfile == NULL) {
        const char * error = gzerror(outfile, &errnum);
        if (errnum == Z_ERRNO) {
            throw khmer_file_exception(strerror(errno));
        } else {
            throw khmer_file_exception(error);
        }
    }

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

    for (unsigned int i = 0; i < save_n_tables; i++) {
        save_tablesize = ht._tablesizes[i];

        gzwrite(outfile, (const char *) &save_tablesize,
                sizeof(save_tablesize));
        unsigned long long written = 0;
        while (written != save_tablesize) {
            written += gzwrite(outfile, (const char *) ht._counts[i],
                               (int) (save_tablesize - written));
        }
    }

    HashIntoType n_counts = ht._bigcounts.size();
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

void CountingHash::collect_high_abundance_kmers(
    const std::string   &filename,
    unsigned int    lower_count,
    unsigned int    upper_count,
    SeenSet&        found_kmers)
{
    unsigned long long total_reads = 0;

    IParser* parser = IParser::get_parser(filename.c_str());
    Read read;

    string currSeq = "";

    //
    // iterate through the FASTA file & consume the reads, until we hit
    // upper_count.
    //

    bool done = false;
    while(!parser->is_complete() && !done)  {
        read = parser->get_next_read();
        currSeq = read.sequence;

        // do we want to process it?
        if (check_and_normalize_read(currSeq)) {
            const char * sp = currSeq.c_str();

            KMerIterator kmers(sp, _ksize);

            while(!kmers.done()) {
                HashIntoType kmer = kmers.next();

                count(kmer);
                if (get_count(kmer) >= upper_count) {
                    done = true;
                }
            }
        }

        // increment read number
        total_reads++;

        if (total_reads % 100000 == 0) {
            std::cout << "..." << total_reads << "\n";
        }
    }

    delete parser;

    unsigned long long stop_at_read = total_reads;

    //
    // go back through the file again, and store all k-mers >= lower_count
    //

    parser = IParser::get_parser(filename.c_str());

    total_reads = 0;
    while(!parser->is_complete() && total_reads != stop_at_read)  {
        read = parser->get_next_read();
        currSeq = read.sequence;

        // do we want to process it?
        if (check_and_normalize_read(currSeq)) {
            const char * sp = currSeq.c_str();

            KMerIterator kmers(sp, _ksize);

            while(!kmers.done()) {
                HashIntoType kmer = kmers.next();

                if (get_count(kmer) >= lower_count) {
                    found_kmers.insert(kmer);
                }
            }
        }

        // increment read number
        total_reads++;

        if (total_reads % 100000 == 0) {
            std::cout << "... x 2 " << total_reads << "\n";
        }
    }
    delete parser;
    parser = NULL;
}

/* vim: set ft=cpp ts=8 sts=4 sw=4 et tw=79 */
