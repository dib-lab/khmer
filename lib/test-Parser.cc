//
// This file is part of khmer, https://github.com/dib-lab/khmer/, and is
// Copyright (C) Michigan State University, 2009-2015. It is licensed under
// the three-clause BSD license; see LICENSE.
// Contact: khmer-project@idyll.org
//

// Test driver for the Parser class.
// Author: Eric McDonald


#include <cerrno>
#include <cstring>
#include <cstdio>
#include <cstdlib>
#include <fcntl.h>
#include <getopt.h>

#include <omp.h>

#include "error.hh"
#include "read_parsers.hh"


// #define WRITE_SUMMARY


using namespace khmer;
using namespace khmer:: read_parsers;


// s: Cache Size
static char const *	    SHORT_OPTS	    = "s:";


int main( int argc, char * argv[ ] )
{
    int		    rc		    = 0;
    Config	    &the_config	    = get_active_config( );
    uint64_t	    cache_size	    =
        the_config.get_reads_input_buffer_size( );
    char *	    ifile_name	    = NULL;

    int		    opt		    = -1;
    char *	    conv_residue    = NULL;
    while (-1 != (opt = getopt( argc, argv, SHORT_OPTS ))) {

        switch (opt) {

        case 's':
            cache_size = strtoull( optarg, &conv_residue, 10 );
            if (!strcmp( optarg, conv_residue )) {
                error( EINVAL, EINVAL, "Invalid cache size" );
            }
            break;

        default:
            error( 0, 0, "Skipping unknown arg, '%c'", optopt );

        } // option switch

    } // getopt loop

    if (optind < argc) {
        ifile_name = argv[ optind++ ];
    } else {
        error( EINVAL, 0, "Input file name required" );
    }
    std:: string    ifile_name_STRING( ifile_name );

    the_config.set_input_buffer_trace_level( TraceLogger:: TLVL_ALL );
    uint32_t	    number_of_threads	    = omp_get_max_threads( );
    IParser *	    parser		    = IParser:: get_parser(
                                          ifile_name_STRING, number_of_threads, cache_size,
                                          TraceLogger:: TLVL_DEBUG6
                                      );

    #pragma omp parallel default( shared )
    {
        uint32_t	thread_id	    = (uint32_t)omp_get_thread_num( );
        Read		the_read;
        uint64_t	seq_len;
        char		ofile_name[ FILENAME_MAX + 1 ];
        FILE *		ofile_handle	    = NULL;

#ifdef WRITE_SUMMARY
        ofile_name[ FILENAME_MAX ] = '\0';
#endif

        fprintf(
            stderr,
            "OMP thread %lu reporting for duty.\n",
            (unsigned long int)thread_id
        );

#ifdef WRITE_SUMMARY
        snprintf(
            ofile_name, FILENAME_MAX, "summary-%lu.log",
            (unsigned long int)thread_id
        );
        ofile_handle = fopen( ofile_name, "w" );
        if (NULL == ofile_handle)
            // TODO: Report an error.
            ;
#endif

        for (uint64_t readnum = 0; !parser->is_complete( ); ++readnum) {

            if (0 == readnum % 100000)
                fprintf(
                    stderr,
                    "OMP thread %lu is on read number %llu.\n",
                    (unsigned long int)thread_id,
                    (unsigned long long int)readnum
                );

            the_read = parser->get_next_read( );

#if (1)
            printf(
                "@%s\n%s\n+\n%s\n",
                the_read.name.c_str( ),
                the_read.sequence.c_str( ),
                the_read.accuracy.c_str( )
            );
#endif

#ifdef WRITE_SUMMARY
            fflush( ofile_handle );
#endif

        }

#ifdef WRITE_SUMMARY
        fclose( ofile_handle );
#endif

    } // parallel region

    delete parser;
    return rc;
}


// vim: set ft=cpp sts=4 sw=4 tw=80:
