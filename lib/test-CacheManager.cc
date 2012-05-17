// Test driver for the CacheManager class.
// Author: Eric McDonald


#include <cerrno>
#include <climits>
extern "C"
{
#include <stdint.h>
}
#include <cstring>
#include <cstdio>
#include <fcntl.h>
#include <sys/stat.h>
#include <cstdlib>
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
    int		    ifd		    = -1;
    IStreamReader * sr		    = NULL;

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

    // TODO: Handle stdin.
    // TODO: Play with O_DIRECT.
    if (-1 == (ifd = open( ifile_name, O_RDONLY )))
	error( errno, errno, "Failed to open input file" );

    try
    {
	if	(!strcmp( "raw", ifile_type ))
	    sr = new RawStreamReader( ifd );
	else if (!strcmp( "gz", ifile_type ))
	    sr = new GzStreamReader( ifd );
	else if (!strcmp( "bz2", ifile_type ))
	    sr = new Bz2StreamReader( ifd );
    }
    catch (InvalidStreamBuffer & exc)
    {
	error( EBADF, EBADF, "Failed to initialize stream reader" );
    }

    uint32_t	    number_of_threads	    = omp_get_max_threads( );
    CacheManager *  cmgr		    = new CacheManager(
	*sr, number_of_threads, cache_size, 3
    );

#pragma omp parallel default( shared )
    {
	uint32_t	thread_id	    = (uint32_t)omp_get_thread_num( );
	drand48_data    rng_state;
	long int	randnum		    = 0;
	uint8_t		buffer[ 127 ];
	uint64_t	segment_cut_pos	    = 0;
	uint64_t	nbread		    = 0;
	uint64_t	nbread_total	    = 0;
	timespec	sleep_duration;
	timespec	sleep_duration_rem;

	fprintf(
	    stderr,
	    "OMP thread %lu reporting for duty.\n",
	    (unsigned long int)thread_id
	);

	srand48_r( (long int)thread_id, &rng_state ); 
	for (uint64_t i = 1; cmgr->has_more_data( ); ++i)
	{

	    if (0 == i % 1000000) 
		fprintf(
		    stderr,
"OMP thread %lu is on data processing iteration %llu.\n",
		    (unsigned long int)thread_id,
		    (unsigned long long int)i
		);
	    
	    lrand48_r( &rng_state, &randnum );
	    randnum %= 128;
	    nbread  = cmgr->get_bytes(
		(uint8_t * const)buffer, (uint64_t)randnum, segment_cut_pos
	    ); 
	    nbread_total += nbread;

	    // Pretend to work for a random duration.
	    lrand48_r( &rng_state, &randnum );
	    //sleep_duration_rem.tv_sec	    = randnum / 1000000000;
	    //sleep_duration_rem.tv_nsec    = randnum % 1000000000;
	    sleep_duration_rem.tv_sec	= 0;
	    sleep_duration_rem.tv_nsec	= randnum % 1000000;
	    while ( sleep_duration_rem.tv_sec && sleep_duration_rem.tv_nsec )
	    {
		sleep_duration.tv_sec	= sleep_duration_rem.tv_sec;
		sleep_duration.tv_nsec	= sleep_duration_rem.tv_nsec;
		nanosleep( &sleep_duration, &sleep_duration_rem );
	    }

	    // Occasionally split setaside buffer,
	    // when opportunity exists.
	    lrand48_r( &rng_state, &randnum );
	    if (    (0 == randnum % 1024)
		&&  (!cmgr->_sa_buffer_avail( ))
		&&  (!sr->is_at_end_of_stream( )))
	    {
		lrand48_r( &rng_state, &randnum );
		randnum %= 1024;
		cmgr->split_at( (uint64_t)randnum );
	    }

	} // work simulator loop

	fprintf(
	    stderr,
	    "OMP thread %lu finished work.\n",
	    (unsigned long int)thread_id
	);

    } // parallel block

    delete cmgr;
    delete sr;
    return rc;
}


// vim: set ft=cpp sts=4 sw=4 tw=80:
