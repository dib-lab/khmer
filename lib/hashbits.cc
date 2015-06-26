//
// This file is part of khmer, https://github.com/dib-lab/khmer/, and is
// Copyright (C) Michigan State University, 2009-2015. It is licensed under
// the three-clause BSD license; see LICENSE.
// Contact: khmer-project@idyll.org
//

#include <iostream>
#include "hashtable.hh"
#include "hashbits.hh"
#include "read_parsers.hh"

#include <sstream>
#include <errno.h>

using namespace std;
using namespace khmer;
using namespace khmer:: read_parsers;

void Hashbits::save(std::string outfilename)
{
    if (!_counts[0]) {
        throw khmer_exception();
    }

    unsigned int save_ksize = _ksize;
    unsigned char save_n_tables = _n_tables;
    unsigned long long save_tablesize;

    ofstream outfile(outfilename.c_str(), ios::binary);

    outfile.write(SAVED_SIGNATURE, 4);
    unsigned char version = SAVED_FORMAT_VERSION;
    outfile.write((const char *) &version, 1);

    unsigned char ht_type = SAVED_HASHBITS;
    outfile.write((const char *) &ht_type, 1);

    outfile.write((const char *) &save_ksize, sizeof(save_ksize));
    outfile.write((const char *) &save_n_tables, sizeof(save_n_tables));

    for (unsigned int i = 0; i < _n_tables; i++) {
        save_tablesize = _tablesizes[i];
        unsigned long long tablebytes = save_tablesize / 8 + 1;

        outfile.write((const char *) &save_tablesize, sizeof(save_tablesize));

        outfile.write((const char *) _counts[i], tablebytes);
    }
    if (outfile.fail()) {
        throw khmer_file_exception(strerror(errno));
    }
    outfile.close();
}

/**
 * Loads @param infilename into Hashbits, with error checking on
 * file type and file version.  Populates _counts internally.
 */
void Hashbits::load(std::string infilename)
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
        throw khmer_file_exception(err);
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
        char signature[4];
        unsigned char version, ht_type;

        infile.read(signature, 4);
        infile.read((char *) &version, 1);
        infile.read((char *) &ht_type, 1);
        if (!(std::string(signature, 4) == SAVED_SIGNATURE)) {
            std::ostringstream err;
            err << "Does not start with signature for a khmer " <<
                "file: " << signature << " Should be: " <<
                SAVED_SIGNATURE;
            throw khmer_file_exception(err.str());
        } else if (!(version == SAVED_FORMAT_VERSION)) {
            std::ostringstream err;
            err << "Incorrect file format version " << (int) version
                << " while reading k-mer graph from " << infilename
                << "; should be " << (int) SAVED_FORMAT_VERSION;
            throw khmer_file_exception(err.str());
        } else if (!(ht_type == SAVED_HASHBITS)) {
            std::ostringstream err;
            err << "Incorrect file format type " << (int) ht_type
                << " while reading k-mer graph from " << infilename;
            throw khmer_file_exception(err.str());
        }

        infile.read((char *) &save_ksize, sizeof(save_ksize));
        infile.read((char *) &save_n_tables, sizeof(save_n_tables));

        _ksize = (WordLength) save_ksize;
        _n_tables = (unsigned int) save_n_tables;
        _init_bitstuff();

        _counts = new Byte*[_n_tables];
        for (unsigned int i = 0; i < _n_tables; i++) {
            HashIntoType tablesize;
            unsigned long long tablebytes;

            infile.read((char *) &save_tablesize, sizeof(save_tablesize));

            tablesize = (HashIntoType) save_tablesize;
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
        throw khmer_file_exception(err);
    }
}

/**
 * Checks for non-ACGT characters before consuming read.
 * This is specifically for counting overlap k-mers.
 */
unsigned int Hashbits::check_and_process_read_overlap(std::string &read,
        bool &is_valid,
        Hashbits &ht2)
{
    is_valid = check_and_normalize_read(read);

    if (!is_valid) {
        return 0;
    }

    return consume_string_overlap(read, ht2);
}

/**
 * Consume a FASTA file of reads.
 */
void Hashbits::consume_fasta_overlap(const std::string &filename,
                                     HashIntoType curve[2][100],Hashbits &ht2,
                                     unsigned int &total_reads,
                                     unsigned long long &n_consumed)
{
    total_reads = 0;
    n_consumed = 0;
    Read read;

//get total number of reads in dataset

    IParser* parser = IParser::get_parser(filename.c_str());
    while(!parser->is_complete())  {
        read = parser->get_next_read();
        total_reads++;
    }
//block size for curve
    int block_size = total_reads/100;

// reads number <100, block size =1
    if (block_size == 0) {
        block_size = 1;
    }
// set the remaining as 0
    for (int n=total_reads; n<100; n++) {
        curve[0][n] = 0;
        curve[1][n] = 0;
    }

    total_reads = 0;

    delete parser;
    parser = IParser::get_parser(filename.c_str());



    string currSeq = "";

    //
    // iterate through the FASTA file & consume the reads.
    //

    while(!parser->is_complete())  {
        read = parser->get_next_read();
        currSeq = read.sequence;

        unsigned int this_n_consumed;
        bool is_valid;

        this_n_consumed = check_and_process_read_overlap(currSeq,
                          is_valid, ht2);

        n_consumed += this_n_consumed;

        // reset the sequence info, increment read number

        total_reads++;

        if (total_reads%block_size == 0) {
            curve[0][total_reads/block_size-1] = n_overlap_kmers();
            curve[1][total_reads/block_size-1] = n_unique_kmers();
        }
    } // while

    delete parser;
}

/**
 * Run through every k-mer in the given string, & hash it.
 */
unsigned int Hashbits::consume_string_overlap(const std::string &s,
        Hashbits &ht2)
{
    const char * sp = s.c_str();
    unsigned int n_consumed = 0;

    KMerIterator kmers(sp, _ksize);

    while(!kmers.done()) {
        HashIntoType kmer = kmers.next();

        count_overlap(kmer,ht2);
        n_consumed++;
    }

    return n_consumed;
}

void Hashbits::update_from(const Hashbits &other)
{
    if (_ksize != other._ksize) {
        throw khmer_exception("both nodegraphs must have same k size");
    }
    if (_tablesizes != other._tablesizes) {
        throw khmer_exception("both nodegraphs must have same table sizes");
    }
    for (unsigned int table_num = 0; table_num < _n_tables; table_num++) {
        Byte * me = _counts[table_num];
        Byte * ot = other._counts[table_num];
        HashIntoType tablesize = _tablesizes[table_num];
        HashIntoType tablebytes = tablesize / 8 + 1;

        for (HashIntoType index = 0; index < tablebytes; index++) {
            me[index] |= ot[index];     // bitwise or
        }
    }
}

// vim: set sts=2 sw=2:
