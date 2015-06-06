//
// This file is part of khmer, https://github.com/dib-lab/khmer/, and is
// Copyright (C) Michigan State University, 2009-2015. It is licensed under
// the three-clause BSD license; see LICENSE.
// Contact: khmer-project@idyll.org
//

// Simple C++ implementation of the 'load-graph' Python script.
// Author: Eric A. McDonald


#include <cstring>
#include <cstdio>
#include <cerrno>
#include <cstdlib>
#include <unistd.h>
#include <getopt.h>

#include <omp.h>

#define HASH_TYPE_TO_TEST   1 // Counting Hash
//#define HASH_TYPE_TO_TEST   2 // Bit Hash

// #define OUTPUT_HASHTABLE


#include "error.hh"
#include "read_parsers.hh"
#if HASH_TYPE_TO_TEST == 1
#  include "counting.hh"
#elif HASH_TYPE_TO_TEST == 2
#  include "hashbits.hh"
#else
#  error "No HASH_TYPE_TO_TEST macro defined."
#endif
#include "primes.hh"

using namespace std;
using namespace khmer;
using namespace khmer:: read_parsers;


static const char *	    SHORT_OPTS		= "k:N:x:s:";


int main( int argc, char * argv[ ] )
{
    unsigned long	kmer_length	    = 32;
    float		ht_size_FP	    = 1.0E6;
    unsigned long	ht_count	    = 4;
    uint64_t		cache_size	    = 4L * 1024 * 1024 * 1024;

    int			rc		    = 0;
    int			opt		    = -1;
    char *		conv_residue	    = NULL;
    string		ofile_name;
    string		ifile_name;
    // FILE *		ofile		    = NULL;

    while (-1 != (opt = getopt( argc, argv, SHORT_OPTS ))) {

        switch (opt) {

        case 'k':
            kmer_length = strtoul( optarg, &conv_residue, 10 );
            if (!strcmp( optarg, conv_residue )) {
                error( EINVAL, EINVAL, "Invalid kmer length" );
            }
            break;

        case 'N':
            ht_count = strtoul( optarg, &conv_residue, 10 );
            if (!strcmp( optarg, conv_residue )) {
                error( EINVAL, EINVAL, "Invalid number of hashtables" );
            }
            break;

        case 'x':
            ht_size_FP = strtof( optarg, &conv_residue );
            if (!strcmp( optarg, conv_residue )) {
                error( EINVAL, EINVAL, "Invalid hashtable size" );
            }
            break;

        case 's':
            cache_size = strtoull( optarg, &conv_residue, 10 );
            if (!strcmp( optarg, conv_residue )) {
                error( EINVAL, EINVAL, "Invalid cache size" );
            }
            break;

        default:
            error( 0, 0, "Skipping unknown arg, '%c'", optopt );
        }

    }

    if (optind < argc) {
        ofile_name = string( argv[ optind++ ] );
    } else {
        error( EINVAL, 0, "Output file name required" );
    }

    if (optind < argc) {
        ifile_name = string( argv[ optind++ ] );
    } else {
        error( EINVAL, 0, "Input file name required" );
    }

    HashIntoType	    ht_size		= (HashIntoType)ht_size_FP;
    Primes primetab( ht_size );
    vector<HashIntoType> ht_sizes;
    for ( unsigned int i = 0; i < ht_count; ++i ) {
        ht_sizes.push_back( primetab.get_next_prime( ) );
    }

    unsigned int	    reads_total		= 0;
    unsigned long long int  n_consumed		= 0;

    Config		    &the_config		= get_active_config( );
    the_config.set_number_of_threads( omp_get_max_threads( ) );

#if HASH_TYPE_TO_TEST == 1
    CountingHash ht( kmer_length, ht_sizes );
    IParser * parser = IParser:: get_parser(
                           ifile_name, the_config.get_number_of_threads( ), cache_size
                       );
    #pragma omp parallel shared( reads_total, n_consumed )
    {
        ht.consume_fasta( parser, reads_total, n_consumed );
    }
#elif HASH_TYPE_TO_TEST == 2
    Hashbits ht( kmer_length, ht_sizes );
    ht.consume_fasta_and_tag( ifile_name, reads_total, n_consumed );
#endif

#ifdef OUTPUT_HASHTABLE
#if	HASH_TYPE_TO_TEST == 1
    ht.save( ofile_name + ".ht_count" );
#elif	HASH_TYPE_TO_TEST == 2
    ht.save( ofile_name + ".ht_bits" );
    ht.save_tagset( ofile_name + ".tagset" );
#endif
#endif

    return rc;
}


// vim: set sts=4 sw=4 tw=80:
