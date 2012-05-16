// Test driver for the Parser class.
// Author: Eric McDonald


#include <cerrno>
#include <climits>
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
#include <cstdlib>
#include <fcntl.h>
#include <error.h>
#include <getopt.h>

#include <omp.h>

#include "read_parsers.hh"


using namespace khmer;
using namespace khmer:: read_parsers;


// t: Stream Type
// s: Cache Size
static char const *	    SHORT_OPTS	    = "t:s:";


int main( int argc, char * argv[ ] )
{
    int		    rc		    = 0;
    char *	    ifile_type	    = (char *)"raw";
    uint64_t	    cache_size	    = 4L * 1024 * 1024 * 1024;
    char *	    ifile_name	    = NULL;

    int		    opt		    = -1;
    char *	    conv_residue    = NULL;
    while (-1 != (opt = getopt( argc, argv, SHORT_OPTS )))
    {
	
	switch (opt)
	{

	case 't':
	    if (    strcmp( optarg, "raw" ) 
		&&  strcmp( optarg, "gz" )
		&&  strcmp( optarg, "bz2" )
	       )
		error( EINVAL, EINVAL, "Invalid file type" );
	    ifile_type = new char[ strlen( optarg ) ];
	    strcpy( ifile_type, optarg );
	    break;
	
	case 's':
	    cache_size = strtoull( optarg, &conv_residue, 10 );
	    if (!strcmp( optarg, conv_residue ))
		error( EINVAL, EINVAL, "Invalid cache size" );
	    break;
	
	default:
	    error( 0, 0, "Skipping unknown arg, '%c'", optopt );
	    
	} // option switch

    } // getopt loop

    if (optind < argc) ifile_name = argv[ optind++ ];
    else error( EINVAL, 0, "Input file name required" );
    std:: string    ifile_name_STRING( ifile_name );

    uint32_t	    number_of_threads	    = omp_get_max_threads( );
    IParser *	    parser		    = IParser:: get_parser(
	ifile_name_STRING, number_of_threads, cache_size,
	TraceLogger:: TLVL_DEBUG8
    );

#pragma omp parallel default( shared )
    {
	uint32_t	thread_id	    = (uint32_t)omp_get_thread_num( );
	Read		the_read;
	uint64_t	seq_len;
	char		ofile_name[ FILENAME_MAX + 1 ];
	FILE *		ofile_handle	    = NULL;

	ofile_name[ FILENAME_MAX ] = '\0';

	fprintf(
	    stderr,
	    "OMP thread %lu reporting for duty.\n",
	    (unsigned long int)thread_id
	);

	snprintf(
	    ofile_name, FILENAME_MAX, "summary-%lu.log",
	    (unsigned long int)thread_id
	);
	ofile_handle = fopen( ofile_name, "w" );
	if (NULL == ofile_handle)
	    // TODO: Report an error.
	    ;

	for (uint64_t readnum = 0; !parser->is_complete( ); ++readnum)
	{

	    //if (0 == readnum % 100000) 
		fprintf(
		    stderr,
		    "OMP thread %lu is on read number %llu.\n",
		    (unsigned long int)thread_id,
		    (unsigned long long int)readnum
		);

	    the_read = parser->get_next_read( );

	    fprintf( ofile_handle, ">%s\n", the_read.name.c_str( ) );
	    seq_len = the_read.seq.length( );
	    fprintf(
		ofile_handle, "%s[%llu]\n",
		the_read.seq.c_str( ), (unsigned long long int)seq_len
	    );

	}

	fclose( ofile_handle );

    } // parallel region

    delete parser;
    return rc;
}


// vim: set ft=cpp sts=4 sw=4 tw=80:
