//
// This file is part of khmer, http://github.com/ged-lab/khmer/, and is
// Copyright (C) Michigan State University, 2009-2013. It is licensed under
// the three-clause BSD license; see doc/LICENSE.txt. 
// Contact: khmer-project@idyll.org
//

#include <iostream>
#include "hashtable.hh"
#include "hashbits.hh"
#include "read_parsers.hh"

using namespace std;
using namespace khmer;
using namespace khmer:: read_parsers;

void Hashbits::save(std::string outfilename)
{
    assert(_counts[0]);

    unsigned int save_ksize = _ksize;
    unsigned char save_n_tables = _n_tables;
    unsigned long long save_tablesize;

    ofstream outfile(outfilename.c_str(), ios::binary);

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
    outfile.close();
}

void Hashbits::load(std::string infilename)
{
    if (_counts) {
        for (unsigned int i = 0; i < _n_tables; i++) {
            delete _counts[i];
            _counts[i] = NULL;
        }
        delete _counts;
        _counts = NULL;
    }
    _tablesizes.clear();

    unsigned int save_ksize = 0;
    unsigned char save_n_tables = 0;
    unsigned long long save_tablesize = 0;
    unsigned char version, ht_type;

    ifstream infile(infilename.c_str(), ios::binary);
    assert(infile.is_open());

    infile.read((char *) &version, 1);
    infile.read((char *) &ht_type, 1);
    assert(version == SAVED_FORMAT_VERSION);
    assert(ht_type == SAVED_HASHBITS);

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
            loaded += infile.gcount();	// do I need to do this loop?
        }
    }
    infile.close();
}

// for counting overlap k-mers specifically!!

//
// check_and_process_read: checks for non-ACGT characters before consuming
//

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

//
// consume_fasta: consume a FASTA file of reads
//

void Hashbits::consume_fasta_overlap(const std::string &filename,
                                     HashIntoType curve[2][100],Hashbits &ht2,
                                     unsigned int &total_reads,
                                     unsigned long long &n_consumed,
                                     CallbackFn callback,
                                     void * callback_data)
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
    for (int n=total_reads; n<100; n++){
            curve[0][n] = 0;
            curve[1][n] = 0;
    }

    total_reads = 0;
    HashIntoType start = 0, stop = 0;

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
            curve[0][total_reads/block_size-1] = n_overlap_kmers(start,stop);
            curve[1][total_reads/block_size-1] = n_kmers(start,stop);
        }
        // run callback, if specified
        if (total_reads % CALLBACK_PERIOD == 0 && callback) {
            try {
                callback("consume_fasta", callback_data, total_reads, n_consumed);
            } catch (...) {
                throw;
            }
        }

    } // while

    delete parser;
}

//
// consume_string: run through every k-mer in the given string, & hash it.
//

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

// vim: set sts=2 sw=2:
