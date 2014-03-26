//
// This file is part of khmer, http://github.com/ged-lab/khmer/, and is
// Copyright (C) Michigan State University, 2009-2013. It is licensed under
// the three-clause BSD license; see doc/LICENSE.txt. 
// Contact: khmer-project@idyll.org
//

// Test driver for the StreamReader classes.
// Author: Eric McDonald


#include <cerrno>
#include <climits>
extern "C"
{
#include <stdint.h>
}
//#define SSIZE_MAX	(SIZE_MAX / 2)
#include <cstring>
#include <cstdio>
#include <fcntl.h>
#include <sys/stat.h>
#include <cstdlib>
#include <getopt.h>

#include "error.hh"
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
    uint8_t *	    cache	    = NULL;
    char *	    ifile_name	    = NULL;
    char *	    ofile_name	    = NULL;
    int		    ifd		    = -1;
    int		    ofd		    = -1;
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

    if (optind < argc) ofile_name = argv[ optind++ ];
    else error( EINVAL, 0, "Output file name required" );

    // TODO: Handle stdin.
    // TODO: Play with O_DIRECT.
    if (-1 == (ifd = open( ifile_name, O_RDONLY )))
	error( errno, errno, "Failed to open input file" );

    // TODO: Handle stdout.
    if (-1 == (ofd = creat( ofile_name, 0644 )))
	error( errno, errno, "Failed to open output file" );

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

    try
    {
	cache = new uint8_t[ cache_size ];
    }
    catch (std:: bad_alloc & exc)
    {
	error( ENOMEM, ENOMEM, "Failed to allocate cache" );
    }

    uint64_t	    nbread	    = 0;
    ssize_t	    nbwrote	    = 0;
    uint64_t	    nbread_total    = 0;
    uint64_t	    nbwrote_total   = 0;
    while (!sr->is_at_end_of_stream( ))
    {
	uint64_t    nbwrote_subtotal	= 0;
	try
	{
	    nbread = sr->read_into_cache( cache, cache_size );
	    nbread_total += nbread;
	    for ( uint64_t nbrem = nbread; 0 < nbrem; nbrem -= nbwrote )
	    {
		nbwrote = 
		write( ofd, 
		       cache + nbwrote_subtotal, 
		       (nbrem > SSIZE_MAX ? SSIZE_MAX : nbrem)
		);
		if (-1 == nbwrote)
		    error( EIO, EIO, "Error during write of output stream" );
		nbwrote_subtotal += nbwrote;
	    }
	    nbwrote_total += nbwrote_subtotal;
	}
	catch (StreamReadError & exc)
	{
	    error( EIO, EIO, "Error during read of input stream" );
	}
	catch (...)
	{
	    throw;
	}
	fprintf( stdout, 
		 "Read %llu bytes from disk.\n", 
		 (long long unsigned int)nbread
	);
	fprintf( stdout, 
		 "Wrote %llu bytes to disk.\n", 
		 (long long unsigned int)nbwrote_subtotal
	);
    } // stream reader read loop
    fprintf( stdout,
	     "Read %llu bytes in total from disk.\n",
	     (long long unsigned int)nbread_total
    );
    fprintf( stdout,
	     "Wrote %llu bytes in total to disk.\n",
	     (long long unsigned int)nbwrote_total
    );

    close( ofd );

    return rc;
}


// vim: set ft=cpp sts=4 sw=4 tw=80:
