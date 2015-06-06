//
// This file is part of khmer, https://github.com/dib-lab/khmer/, and is
// Copyright (C) Michigan State University, 2009-2015. It is licensed under
// the three-clause BSD license; see LICENSE.
// Contact: khmer-project@idyll.org
//

// Simple C++ implementation of a diff between counting hashes.
// Author: Eric A. McDonald

// You can learn which hash bins have differing values with a simple 'cmp'
// between two hash files (if yo uaccount for file header length),
// but this program actually loads the tables into memory
// and checks things such as number of hash tables and hash table sizes.
// Also, any differences in count or bigcount values are reported in
// human-readable form.

#if (__cplusplus >= 201103L)
#   include <cstdint>
#else
extern "C"
{
#   include <stdint.h>
}
#endif
#include <cstring>
#include <cstdio>
#include <cerrno>
#include <cstdlib>

#include <string>

#include "khmer.hh"
#include "error.hh"
#include "counting.hh"

using namespace std;
using namespace khmer;


static const char *	    SHORT_OPTS		= "C:R";


int main( int argc, char * argv[ ] )
{
    int			rc		    = 0;
    int			opt		    = -1;
    char *		conv_residue	    = NULL;
    uint32_t		max_count	    = MAX_KCOUNT;
    bool		report_all	    = false;
    string		ifile_name_1;
    string		ifile_name_2;

    while (-1 != (opt = getopt( argc, argv, SHORT_OPTS ))) {

        switch (opt) {
        case 'C':
            max_count = (uint32_t)strtoul( optarg, &conv_residue, 10 );
            if (!strcmp( optarg, conv_residue )) {
                error( EINVAL, EINVAL, "Invalid count threshold" );
            }
            break;
        case 'R':
            report_all = true;
            break;
        default:
            error( 0, 0, "Skipping unknown arg, '%c'", optopt );
        }

    }

    if (optind < argc) {
        ifile_name_1 = string( argv[ optind++ ] );
    } else {
        error( EINVAL, 0, "Name of first hash table file required" );
    }

    if (optind < argc) {
        ifile_name_2 = string( argv[ optind++ ] );
    } else {
        error( EINVAL, 0, "Name of second hash table file required" );
    }

    CountingHash ht1( 20, 1 );
    CountingHash ht2( 20, 1 );
    printf( "Loading hash tables into memory....\n" );
    ht1.load( ifile_name_1 );
    ht2.load( ifile_name_2 );

    HashIntoType i = 0, max_ht_size = 0;
    std:: vector<HashIntoType> ht1_sizes = ht1.get_tablesizes( );
    std:: vector<HashIntoType> ht2_sizes = ht2.get_tablesizes( );

    // Compare number of tables.
    if (ht1_sizes.size( ) != ht2_sizes.size( )) {
        fprintf(
            stderr, "Unequal number of hashtables (%lu and %lu).\n",
            (unsigned long int)ht1_sizes.size( ),
            (unsigned long int)ht2_sizes.size( )
        );
        exit( 1 );
    } else
        printf(
            "Number of Hash Tables: %lu\n",
            (unsigned long int)ht1_sizes.size( )
        );

    // Compare sizes of each table.
    for (i = 0; i < ht1_sizes.size( ); ++i) {
        if (ht1_sizes[ i ] != ht2_sizes[ i ]) {
            fprintf(
                stderr, "Hash table %lu has mismatched sizes of %llu and %llu.\n",
                (unsigned long int)i, ht1_sizes[ i ], ht2_sizes[ i ]
            );
            exit( 1 );
        } else {
            printf(
                "Size of Hash Table %lu: %llu bins\n",
                (unsigned long int)i, ht1_sizes[ i ]
            );
            if (max_ht_size < ht1_sizes[ i ]) {
                max_ht_size = ht1_sizes[ i ];
            }
        }
    }

    printf( "Scanning hash key space....\n" );
    for (i = 0; i < max_ht_size; ++i) {
        // Truncate counts at specified saturation threshold.
        // (This accounts for the sloppy counting used for >1 threads.)
        uint32_t count1 = MIN( ht1.get_count( i ), max_count );
        uint32_t count2 = MIN( ht2.get_count( i ), max_count );
        if (count1 != count2) {
            fprintf(
                stderr, "Hash key %llu has mismatched counts of %u and %u.\n",
                i, ht1.get_count( i ), ht2.get_count( i )
            );
            if (!report_all) {
                exit( 1 );
            }
        }
    }
    // TODO: Implement bigcount checking.

    return rc;

}

// vim: set sts=4 sw=4 tw=80:
