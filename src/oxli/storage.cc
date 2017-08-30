/*
This file is part of khmer, https://github.com/dib-lab/khmer/, and is
Copyright (C) 2016, The Regents of the University of California.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are
met:

    * Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.

    * Redistributions in binary form must reproduce the above
      copyright notice, this list of conditions and the following
      disclaimer in the documentation and/or other materials provided
      with the distribution.

    * Neither the name of the University of California nor the names
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
#include <sstream> // IWYU pragma: keep
#include <fstream>
#include <iostream>

#include "oxli/oxli_exception.hh"
#include "oxli/hashtable.hh"
#include "zlib.h"

using namespace oxli;
using namespace std;

void Storage::set_use_bigcount(bool b)
{
    if (!_supports_bigcount) {
        throw oxli_exception("bigcount is not supported for this storage.");
    }
    _use_bigcount = b;
}

bool Storage::get_use_bigcount()
{
    return _use_bigcount;
}

void BitStorage::update_from(const BitStorage& other)
{
    if (_tablesizes != other._tablesizes) {
        throw oxli_exception("both nodegraphs must have same table sizes");
    }

    Byte tmp = 0;
    for (unsigned int table_num = 0; table_num < _n_tables; table_num++) {
        Byte * me = _counts[table_num];
        Byte * ot = other._counts[table_num];
        uint64_t tablesize = _tablesizes[table_num];
        uint64_t tablebytes = tablesize / 8 + 1;

        for (uint64_t index = 0; index < tablebytes; index++) {
            // Bloom filters can be unioned with bitwise OR.
            // First, get the new value
            tmp = me[index] | ot[index];
            if (table_num == 0) {
                // We'd like for the merged filter to have an accurate
                // count of occupied bins.  First, observe that
                // HammingDistance(x,y) is equivalent to
                // HammingWeight(x^y).  Then, observe that the number
                // of additional occupied bins from the update is the
                // hamming distance between the original bin and the
                // OR'd bin. Thus, we can use the builtin popcountll
                // function, which calls a hardware instruction for
                // hamming weight, with the original and merged bin,
                // to find the number of additional occupied bins.
                _occupied_bins += __builtin_popcountll(me[index] ^ tmp);
            }
            me[index] = tmp;
        }
    }
}


void BitStorage::save(std::string outfilename, WordLength ksize)
{
    if (!_counts[0]) {
        throw oxli_exception();
    }

    unsigned int save_ksize = ksize;
    unsigned char save_n_tables = _n_tables;
    unsigned long long save_tablesize;
    unsigned long long save_occupied_bins = _occupied_bins;

    ofstream outfile(outfilename.c_str(), ios::binary);

    outfile.write(SAVED_SIGNATURE, 4);
    unsigned char version = SAVED_FORMAT_VERSION;
    outfile.write((const char *) &version, 1);

    unsigned char ht_type = SAVED_HASHBITS;
    outfile.write((const char *) &ht_type, 1);

    outfile.write((const char *) &save_ksize, sizeof(save_ksize));
    outfile.write((const char *) &save_n_tables, sizeof(save_n_tables));
    outfile.write((const char *) &save_occupied_bins,
                  sizeof(save_occupied_bins));

    for (unsigned int i = 0; i < _n_tables; i++) {
        save_tablesize = _tablesizes[i];
        unsigned long long tablebytes = save_tablesize / 8 + 1;

        outfile.write((const char *) &save_tablesize, sizeof(save_tablesize));

        outfile.write((const char *) _counts[i], tablebytes);
    }
    if (outfile.fail()) {
        throw oxli_file_exception(strerror(errno));
    }
    outfile.close();
}

/**
 * Loads @param infilename into BitStorage, with error checking on
 * file type and file version.  Populates _counts internally.
 */
void BitStorage::load(std::string infilename, WordLength &ksize)
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
            err = "Cannot open k-mer graph file: " + infilename;
        } else {
            err = "Unknown error in opening file: " + infilename;
        }
        throw oxli_file_exception(err);
    } catch (const std::exception &e) {
        // Catching std::exception is a stopgap for
        // https://gcc.gnu.org/bugzilla/show_bug.cgi?id=66145
        std::string err = "Unknown error opening file: " + infilename + " "
                          + strerror(errno);
        throw oxli_file_exception(err);
    }

    if (_counts) {
        for (unsigned int i = 0; i < _n_tables; i++) {
            delete[] _counts[i];
            _counts[i] = NULL;
        }
        delete[] _counts;
        _counts = NULL;
    }
    _tablesizes.clear();

    try {
        unsigned int save_ksize = 0;
        unsigned char save_n_tables = 0;
        unsigned long long save_tablesize = 0;
        unsigned long long save_occupied_bins = 0;
        char signature[4];
        unsigned char version, ht_type;

        infile.read(signature, 4);
        infile.read((char *) &version, 1);
        infile.read((char *) &ht_type, 1);
        if (!(std::string(signature, 4) == SAVED_SIGNATURE)) {
            std::ostringstream err;
            err << "Does not start with signature for a oxli file: 0x";
            for(size_t i=0; i < 4; ++i) {
                err << std::hex << (int) signature[i];
            }
            err << " Should be: " << SAVED_SIGNATURE;
            throw oxli_file_exception(err.str());
        } else if (!(version == SAVED_FORMAT_VERSION)) {
            std::ostringstream err;
            err << "Incorrect file format version " << (int) version
                << " while reading k-mer graph from " << infilename
                << "; should be " << (int) SAVED_FORMAT_VERSION;
            throw oxli_file_exception(err.str());
        } else if (!(ht_type == SAVED_HASHBITS)) {
            std::ostringstream err;
            err << "Incorrect file format type " << (int) ht_type
                << " while reading k-mer graph from " << infilename;
            throw oxli_file_exception(err.str());
        }

        infile.read((char *) &save_ksize, sizeof(save_ksize));
        infile.read((char *) &save_n_tables, sizeof(save_n_tables));
        infile.read((char *) &save_occupied_bins, sizeof(save_occupied_bins));

        ksize = (WordLength) save_ksize;
        _n_tables = (unsigned int) save_n_tables;
        _occupied_bins = save_occupied_bins;

        _counts = new Byte*[_n_tables];
        for (unsigned int i = 0; i < _n_tables; i++) {
            uint64_t tablesize;
            unsigned long long tablebytes;

            infile.read((char *) &save_tablesize, sizeof(save_tablesize));

            tablesize = save_tablesize;
            _tablesizes.push_back(tablesize);

            tablebytes = tablesize / 8 + 1;
            _counts[i] = new Byte[tablebytes];

            unsigned long long loaded = 0;
            while (loaded != tablebytes) {
                infile.read((char *) _counts[i], tablebytes - loaded);
                loaded += infile.gcount();
            }
        }
        infile.close();
    } catch (std::ifstream::failure &e) {
        std::string err;
        if (infile.eof()) {
            err = "Unexpected end of k-mer graph file: " + infilename;
        } else {
            err = "Error reading from k-mer graph file: " + infilename;
        }
        throw oxli_file_exception(err);
    } catch (const std::exception &e) {
        // Catching std::exception is a stopgap for
        // https://gcc.gnu.org/bugzilla/show_bug.cgi?id=66145
        std::string err = "Unknown error opening file: " + infilename + " "
                          + strerror(errno);
        throw oxli_file_exception(err);
    }
}

void ByteStorageFile::save(
    const std::string   &outfilename,
    WordLength ksize,
    const ByteStorage &store)
{
    std::string filename(outfilename);
    size_t found = filename.find_last_of(".");
    std::string type = filename.substr(found + 1);

    if (type == "gz") {
        ByteStorageGzFileWriter(filename, ksize, store);
    } else {
        ByteStorageFileWriter(filename, ksize, store);
    }
}


ByteStorageFileReader::ByteStorageFileReader(
    const std::string   &infilename,
    WordLength& ksize,
    ByteStorage &store)
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
        throw oxli_file_exception(err + " " + strerror(errno));
    } catch (const std::exception &e) {
        std::string err = "Unknown error opening file: " + infilename + " "
                          + strerror(errno);
        throw oxli_file_exception(err);
    }

    if (store._counts) {
        for (unsigned int i = 0; i < store._n_tables; i++) {
            delete[] store._counts[i];
            store._counts[i] = NULL;
        }
        delete[] store._counts;
        store._counts = NULL;
    }
    store._tablesizes.clear();

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
            err << "Does not start with signature for a oxli file: 0x";
            for(size_t i=0; i < 4; ++i) {
                err << std::hex << (int) signature[i];
            }
            err << " Should be: " << SAVED_SIGNATURE;
            throw oxli_file_exception(err.str());
        } else if (!(version == SAVED_FORMAT_VERSION)) {
            std::ostringstream err;
            err << "Incorrect file format version " << (int) version
                << " while reading k-mer count file from " << infilename
                << "; should be " << (int) SAVED_FORMAT_VERSION;
            throw oxli_file_exception(err.str());
        } else if (!(ht_type == SAVED_COUNTING_HT)) {
            std::ostringstream err;
            err << "Incorrect file format type " << (int) ht_type
                << " while reading k-mer count file from " << infilename;
            throw oxli_file_exception(err.str());
        }

        infile.read((char *) &use_bigcount, 1);
        infile.read((char *) &save_ksize, sizeof(save_ksize));
        infile.read((char *) &save_n_tables, sizeof(save_n_tables));
        infile.read((char *) &save_occupied_bins, sizeof(save_occupied_bins));

        ksize = (WordLength) save_ksize;
        store._n_tables = (unsigned int) save_n_tables;
        store._occupied_bins = save_occupied_bins;

        store._use_bigcount = use_bigcount;

        store._counts = new Byte*[store._n_tables];
        for (unsigned int i = 0; i < store._n_tables; i++) {
            store._counts[i] = NULL;
        }

        for (unsigned int i = 0; i < store._n_tables; i++) {
            uint64_t tablesize;

            infile.read((char *) &save_tablesize, sizeof(save_tablesize));

            tablesize = save_tablesize;
            store._tablesizes.push_back(tablesize);

            store._counts[i] = new Byte[tablesize];

            unsigned long long loaded = 0;
            while (loaded != tablesize) {
                infile.read((char *) store._counts[i], tablesize - loaded);
                loaded += infile.gcount();
            }
        }

        uint64_t n_counts = 0;
        infile.read((char *) &n_counts, sizeof(n_counts));

        if (n_counts) {
            store._bigcounts.clear();

            HashIntoType kmer;
            BoundedCounterType count;

            for (uint64_t n = 0; n < n_counts; n++) {
                infile.read((char *) &kmer, sizeof(kmer));
                infile.read((char *) &count, sizeof(count));
                store._bigcounts[kmer] = count;
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
        throw oxli_file_exception(err);
    } catch (const std::exception &e) {
        std::string err = "Error reading from k-mer count file: " + infilename + " "
                          + strerror(errno);
        throw oxli_file_exception(err);
    }
}

ByteStorageGzFileReader::ByteStorageGzFileReader(
    const std::string   &infilename,
    WordLength &ksize,
    ByteStorage    &store)
{
    gzFile infile = gzopen(infilename.c_str(), "rb");
    if (infile == Z_NULL) {
        std::string err = "Cannot open k-mer count file: " + infilename;
        throw oxli_file_exception(err);
    }

    if (store._counts) {
        for (unsigned int i = 0; i < store._n_tables; i++) {
            delete[] store._counts[i];
            store._counts[i] = NULL;
        }
        delete[] store._counts;
        store._counts = NULL;
    }
    store._tablesizes.clear();

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
        throw oxli_file_exception(err);
    } else if (!(std::string(signature, 4) == SAVED_SIGNATURE)) {
        std::ostringstream err;
        err << "Does not start with signature for a oxli " <<
            "file: " << signature << " Should be: " <<
            SAVED_SIGNATURE;
        throw oxli_file_exception(err.str());
    } else if (!(version == SAVED_FORMAT_VERSION)
               || !(ht_type == SAVED_COUNTING_HT)) {
        if (!(version == SAVED_FORMAT_VERSION)) {
            std::ostringstream err;
            err << "Incorrect file format version " << (int) version
                << " while reading k-mer count file from " << infilename
                << "; should be " << (int) SAVED_FORMAT_VERSION;
            gzclose(infile);
            throw oxli_file_exception(err.str());
        } else if (!(ht_type == SAVED_COUNTING_HT)) {
            std::ostringstream err;
            err << "Incorrect file format type " << (int) ht_type
                << " while reading k-mer count file from " << infilename;
            gzclose(infile);
            throw oxli_file_exception(err.str());
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
        throw oxli_file_exception(err);
    }

    ksize = (WordLength) save_ksize;
    store._occupied_bins = save_occupied_bins;
    store._n_tables = (unsigned int) save_n_tables;

    store._use_bigcount = use_bigcount;

    store._counts = new Byte*[store._n_tables];
    for (unsigned int i = 0; i < store._n_tables; i++) {
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
            throw oxli_file_exception(err);
        }

        tablesize = save_tablesize;
        store._tablesizes.push_back(tablesize);

        store._counts[i] = new Byte[tablesize];

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
            read_b = gzread(infile, (char *) store._counts[i], to_read_int);

            if (read_b <= 0) {
                std::string gzerr = gzerror(infile, &read_b);
                std::string err = "K-mer count file read error: " + infilename;
                if (read_b == Z_ERRNO) {
                    err = err + " " + strerror(errno);
                } else {
                    err = err + " " + gzerr;
                }
                gzclose(infile);
                throw oxli_file_exception(err);
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
        throw oxli_file_exception(err);
    }

    if (n_counts) {
        store._bigcounts.clear();

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
                throw oxli_file_exception(err);
            }

            store._bigcounts[kmer] = count;
        }
    }

    gzclose(infile);
}

ByteStorageFileWriter::ByteStorageFileWriter(
    const std::string   &outfilename,
    const WordLength ksize,
    const ByteStorage& store)
{
    if (!store._counts[0]) {
        throw oxli_exception();
    }

    unsigned int save_ksize = ksize;
    unsigned char save_n_tables = store._n_tables;
    unsigned long long save_tablesize;
    unsigned long long save_occupied_bins = store._occupied_bins;

    ofstream outfile(outfilename.c_str(), ios::binary);

    outfile.write(SAVED_SIGNATURE, 4);
    unsigned char version = SAVED_FORMAT_VERSION;
    outfile.write((const char *) &version, 1);

    unsigned char ht_type = SAVED_COUNTING_HT;
    outfile.write((const char *) &ht_type, 1);

    unsigned char use_bigcount = 0;
    if (store._use_bigcount) {
        use_bigcount = 1;
    }
    outfile.write((const char *) &use_bigcount, 1);

    outfile.write((const char *) &save_ksize, sizeof(save_ksize));
    outfile.write((const char *) &save_n_tables, sizeof(save_n_tables));
    outfile.write((const char *) &save_occupied_bins,
                  sizeof(save_occupied_bins));

    for (unsigned int i = 0; i < save_n_tables; i++) {
        save_tablesize = store._tablesizes[i];

        outfile.write((const char *) &save_tablesize, sizeof(save_tablesize));
        outfile.write((const char *) store._counts[i], save_tablesize);
    }

    uint64_t n_counts = store._bigcounts.size();
    outfile.write((const char *) &n_counts, sizeof(n_counts));

    if (n_counts) {
        KmerCountMap::const_iterator it = store._bigcounts.begin();

        for (; it != store._bigcounts.end(); ++it) {
            outfile.write((const char *) &it->first, sizeof(it->first));
            outfile.write((const char *) &it->second, sizeof(it->second));
        }
    }
    if (outfile.fail()) {
        throw oxli_file_exception(strerror(errno));
    }
    outfile.close();
}

ByteStorageGzFileWriter::ByteStorageGzFileWriter(
    const std::string   &outfilename,
    const WordLength ksize,
    const ByteStorage &store)
{
    if (!store._counts[0]) {
        throw oxli_exception();
    }

    int errnum = 0;
    unsigned int save_ksize = ksize;
    unsigned char save_n_tables = store._n_tables;
    unsigned long long save_tablesize;
    unsigned long long save_occupied_bins = store._occupied_bins;

    gzFile outfile = gzopen(outfilename.c_str(), "wb");
    if (outfile == NULL) {
        const char * error = gzerror(outfile, &errnum);
        if (errnum == Z_ERRNO) {
            throw oxli_file_exception(strerror(errno));
        } else {
            throw oxli_file_exception(error);
        }
    }

    gzwrite(outfile, SAVED_SIGNATURE, 4);
    unsigned char version = SAVED_FORMAT_VERSION;
    gzwrite(outfile, (const char *) &version, 1);

    unsigned char ht_type = SAVED_COUNTING_HT;
    gzwrite(outfile, (const char *) &ht_type, 1);

    unsigned char use_bigcount = 0;
    if (store._use_bigcount) {
        use_bigcount = 1;
    }
    gzwrite(outfile, (const char *) &use_bigcount, 1);

    gzwrite(outfile, (const char *) &save_ksize, sizeof(save_ksize));
    gzwrite(outfile, (const char *) &save_n_tables, sizeof(save_n_tables));
    gzwrite(outfile, (const char *) &save_occupied_bins,
            sizeof(save_occupied_bins));

    for (unsigned int i = 0; i < save_n_tables; i++) {
        save_tablesize = store._tablesizes[i];

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
            gz_result = gzwrite(outfile, (const char *) store._counts[i],
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
                throw oxli_file_exception(msg.str());
            }
            written += gz_result;
        }
    }

    uint64_t n_counts = store._bigcounts.size();
    gzwrite(outfile, (const char *) &n_counts, sizeof(n_counts));

    if (n_counts) {
        KmerCountMap::const_iterator it = store._bigcounts.begin();

        for (; it != store._bigcounts.end(); ++it) {
            gzwrite(outfile, (const char *) &it->first, sizeof(it->first));
            gzwrite(outfile, (const char *) &it->second, sizeof(it->second));
        }
    }
    const char * error = gzerror(outfile, &errnum);
    if (errnum == Z_ERRNO) {
        throw oxli_file_exception(strerror(errno));
    } else if (errnum != Z_OK) {
        throw oxli_file_exception(error);
    }
    gzclose(outfile);
}

void ByteStorageFile::load(
    const std::string   &infilename,
    WordLength &ksize,
    ByteStorage &store)
{
    std::string filename(infilename);
    size_t found = filename.find_last_of(".");
    std::string type = filename.substr(found + 1);

    if (type == "gz") {
        ByteStorageGzFileReader(filename, ksize, store);
    } else {
        ByteStorageFileReader(filename, ksize, store);
    }
}

void ByteStorage::save(std::string outfilename, WordLength ksize)
{
    ByteStorageFile::save(outfilename, ksize, *this);
}

void ByteStorage::load(std::string infilename, WordLength& ksize)
{
    ByteStorageFile::load(infilename, ksize, *this);
}


void NibbleStorage::save(std::string outfilename, WordLength ksize)
{
    if (!_counts[0]) {
        throw oxli_exception();
    }

    unsigned int save_ksize = ksize;
    unsigned char save_n_tables = _n_tables;
    unsigned long long save_tablesize;
    unsigned long long save_occupied_bins = _occupied_bins;

    ofstream outfile(outfilename.c_str(), ios::binary);

    outfile.write(SAVED_SIGNATURE, 4);
    unsigned char version = SAVED_FORMAT_VERSION;
    outfile.write((const char *) &version, 1);

    unsigned char ht_type = SAVED_SMALLCOUNT;
    outfile.write((const char *) &ht_type, 1);

    outfile.write((const char *) &save_ksize, sizeof(save_ksize));
    outfile.write((const char *) &save_n_tables, sizeof(save_n_tables));
    outfile.write((const char *) &save_occupied_bins,
                  sizeof(save_occupied_bins));

    for (unsigned int i = 0; i < save_n_tables; i++) {
        save_tablesize = _tablesizes[i];

        outfile.write((const char *) &save_tablesize, sizeof(save_tablesize));
        outfile.write((const char *) _counts[i], save_tablesize / 2 + 1);
    }
}

void NibbleStorage::load(std::string infilename, WordLength& ksize)
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
        throw oxli_file_exception(err + " " + strerror(errno));
    } catch (const std::exception &e) {
        std::string err = "Unknown error opening file: " + infilename + " "
                          + strerror(errno);
        throw oxli_file_exception(err);
    }

    if (_counts) {
        for (unsigned int i = 0; i < _n_tables; i++) {
            delete[] _counts[i];
            _counts[i] = NULL;
        }
        delete[] _counts;
        _counts = NULL;
    }
    _tablesizes.clear();

    try {
        unsigned int save_ksize = 0;
        unsigned char save_n_tables = 0;
        unsigned long long save_tablesize = 0;
        unsigned long long save_occupied_bins = 0;
        char signature [4];
        unsigned char version = 0, ht_type = 0;

        infile.read(signature, 4);
        infile.read((char *) &version, 1);
        infile.read((char *) &ht_type, 1);
        if (!(std::string(signature, 4) == SAVED_SIGNATURE)) {
            std::ostringstream err;
            err << "Does not start with signature for a oxli file: 0x";
            for(size_t i=0; i < 4; ++i) {
                err << std::hex << (int) signature[i];
            }
            err << " Should be: " << SAVED_SIGNATURE;
            throw oxli_file_exception(err.str());
        } else if (!(version == SAVED_FORMAT_VERSION)) {
            std::ostringstream err;
            err << "Incorrect file format version " << (int) version
                << " while reading k-mer count file from " << infilename
                << "; should be " << (int) SAVED_FORMAT_VERSION;
            throw oxli_file_exception(err.str());
        } else if (!(ht_type == SAVED_SMALLCOUNT)) {
            std::ostringstream err;
            err << "Incorrect file format type " << (int) ht_type
                << " while reading k-mer count file from " << infilename;
            throw oxli_file_exception(err.str());
        }

        infile.read((char *) &save_ksize, sizeof(save_ksize));
        infile.read((char *) &save_n_tables, sizeof(save_n_tables));
        infile.read((char *) &save_occupied_bins, sizeof(save_occupied_bins));

        ksize = (WordLength) save_ksize;
        _n_tables = (unsigned int) save_n_tables;
        _occupied_bins = save_occupied_bins;

        _counts = new Byte*[_n_tables];
        for (unsigned int i = 0; i < _n_tables; i++) {
            _counts[i] = NULL;
        }

        for (unsigned int i = 0; i < _n_tables; i++) {
            uint64_t tablesize;
            uint64_t tablebytes;

            infile.read((char *) &save_tablesize, sizeof(save_tablesize));

            tablebytes = save_tablesize / 2 + 1;
            tablesize = save_tablesize;
            _tablesizes.push_back(tablesize);

            _counts[i] = new Byte[tablebytes];

            unsigned long long loaded = 0;
            while (loaded != tablebytes) {
                infile.read((char *) _counts[i], tablebytes - loaded);
                loaded += infile.gcount();
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
        throw oxli_file_exception(err);
    } catch (const std::exception &e) {
        std::string err = "Error reading from k-mer count file: " + infilename + " "
                          + strerror(errno);
        throw oxli_file_exception(err);
    }
}


void QFStorage::save(std::string outfilename, WordLength ksize)
{
    ofstream outfile(outfilename.c_str(), ios::binary);

    unsigned char version = SAVED_FORMAT_VERSION;
    unsigned char ht_type = SAVED_QFCOUNT;

    outfile.write(SAVED_SIGNATURE, 4);
    outfile.write((const char *) &version, 1);
    outfile.write((const char *) &ht_type, 1);
    outfile.write((const char *) &ksize, sizeof(ksize));

    /* just a hack to handle __uint128_t value. Don't know a better to handle it
     * right now */
    uint64_t tmp_range;
    tmp_range = cf.range;

    outfile.write((const char *) &cf.nslots, sizeof(cf.nslots));
    outfile.write((const char *) &cf.xnslots, sizeof(cf.xnslots));
    outfile.write((const char *) &cf.key_bits, sizeof(cf.key_bits));
    outfile.write((const char *) &cf.value_bits, sizeof(cf.value_bits));
    outfile.write((const char *) &cf.key_remainder_bits, sizeof(cf.key_remainder_bits));
    outfile.write((const char *) &cf.bits_per_slot, sizeof(cf.bits_per_slot));
    outfile.write((const char *) &tmp_range, sizeof(tmp_range));
    outfile.write((const char *) &cf.nblocks, sizeof(cf.nblocks));
    outfile.write((const char *) &cf.nelts, sizeof(cf.nelts));
    outfile.write((const char *) &cf.ndistinct_elts, sizeof(cf.ndistinct_elts));
    outfile.write((const char *) &cf.noccupied_slots, sizeof(cf.noccupied_slots));

    #if BITS_PER_SLOT == 8 || BITS_PER_SLOT == 16 || BITS_PER_SLOT == 32 || BITS_PER_SLOT == 64
        outfile.write((const char *) cf.blocks, sizeof(qfblock) * cf.nblocks);
    #else
        outfile.write((const char *) cf.blocks,
                      (sizeof(qfblock) + SLOTS_PER_BLOCK * cf.bits_per_slot / 8) * cf.nblocks);
    #endif
    outfile.close();
}


void QFStorage::load(std::string infilename, WordLength &ksize)
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
        throw oxli_file_exception(err + " " + strerror(errno));
    } catch (const std::exception &e) {
        std::string err = "Unknown error opening file: " + infilename + " "
                          + strerror(errno);
        throw oxli_file_exception(err);
    }
    WordLength save_ksize = 0;
    char signature [4];
    unsigned char version = 0, ht_type = 0;
    uint64_t tmp_range;

    infile.read(signature, 4);
    infile.read((char *) &version, 1);
    infile.read((char *) &ht_type, 1);

    if (!(std::string(signature, 4) == SAVED_SIGNATURE)) {
        std::ostringstream err;
        err << "Does not start with signature for a oxli file: 0x";
        for(size_t i=0; i < 4; ++i) {
            err << std::hex << (int) signature[i];
        }
        err << " Should be: " << SAVED_SIGNATURE;
        throw oxli_file_exception(err.str());
    } else if (!(version == SAVED_FORMAT_VERSION)) {
        std::ostringstream err;
        err << "Incorrect file format version " << (int) version
            << " while reading k-mer count file from " << infilename
            << "; should be " << (int) SAVED_FORMAT_VERSION;
        throw oxli_file_exception(err.str());
    } else if (!(ht_type == SAVED_QFCOUNT)) {
        std::ostringstream err;
        err << "Incorrect file format type " << (int) ht_type
            << " expected " << (int) SAVED_QFCOUNT
            << " while reading k-mer count file from " << infilename;
        throw oxli_file_exception(err.str());
    }

    infile.read((char *) &save_ksize, sizeof(save_ksize));
    ksize = save_ksize;

    infile.read((char *) &cf.nslots, sizeof(cf.nslots));
    infile.read((char *) &cf.xnslots, sizeof(cf.xnslots));
    infile.read((char *) &cf.key_bits, sizeof(cf.key_bits));
    infile.read((char *) &cf.value_bits, sizeof(cf.value_bits));
    infile.read((char *) &cf.key_remainder_bits, sizeof(cf.key_remainder_bits));
    infile.read((char *) &cf.bits_per_slot, sizeof(cf.bits_per_slot));
    infile.read((char *) &tmp_range, sizeof(tmp_range));

    infile.read((char *) &cf.nblocks, sizeof(cf.nblocks));
    infile.read((char *) &cf.nelts, sizeof(cf.nelts));
    infile.read((char *) &cf.ndistinct_elts, sizeof(cf.ndistinct_elts));
    infile.read((char *) &cf.noccupied_slots, sizeof(cf.noccupied_slots));
    /* just a hack to handle __uint128_t value. Don't know a better to handle it
     * right now */
    cf.range = tmp_range;
    // deallocate previously allocated blocks
    free(cf.blocks);
    /* allocate the space for the actual qf blocks */
    #if BITS_PER_SLOT == 8 || BITS_PER_SLOT == 16 || BITS_PER_SLOT == 32 || BITS_PER_SLOT == 64
        cf.blocks = (qfblock *)calloc(cf.nblocks, sizeof(qfblock));
    #else
        cf.blocks = (qfblock *)calloc(cf.nblocks, sizeof(qfblock) + SLOTS_PER_BLOCK * cf.bits_per_slot / 8);
    #endif
    #if BITS_PER_SLOT == 8 || BITS_PER_SLOT == 16 || BITS_PER_SLOT == 32 || BITS_PER_SLOT == 64
        infile.read((char *) cf.blocks, sizeof(qfblock) * cf.nblocks);
    #else
        infile.read((char *) cf.blocks,
                    (sizeof(qfblock) + SLOTS_PER_BLOCK * cf.bits_per_slot / 8) * cf.nblocks);
    #endif
    infile.close();
}
