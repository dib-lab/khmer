#include <climits>
extern "C"
{
#include <stdint.h>
}
#ifndef SSIZE_MAX
#   define SSIZE_MAX	((ssize_t)(SIZE_MAX / 2))
#endif
#include <cstring>

#include <cstdio>
#include <cassert>

#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#ifdef __linux__
#   include <sys/syscall.h>
#endif

#include <ctime>

#include "read_parsers.hh"


// Convenience function for getting the difference, in nanoseconds, bewteen 
// two timespec structures.
inline
static
uint64_t
_timespec_diff_in_nsecs( timespec & tstart, timespec & tstop )
{
    return
	    ((tstop.tv_sec * 1000000000U) + (uint64_t)tstop.tv_nsec)
	-   ((tstart.tv_sec * 1000000000U) + (uint64_t)tstart.tv_nsec);
}


namespace khmer
{


ThreadIDMap::
ThreadIDMap( uint32_t number_of_threads )
:   _number_of_threads( number_of_threads ),
    _thread_counter( 0 ),
    _tid_map_spin_lock( 0 )
{
    if (0 == number_of_threads) throw InvalidNumberOfThreadsRequested( );
}


ThreadIDMap::
~ThreadIDMap( )
{
    _thread_id_map.clear( );
}


uint32_t const
ThreadIDMap::
get_thread_id( )
{
#ifdef __linux__
    // Note: No error handling because this call always succeeds, allegedly.
    pid_t native_thread_id = syscall( SYS_gettid ); 
    std:: map< pid_t, uint32_t > :: iterator match;

#else
    // TODO: Maybe try something with pthread_self for the general case.
#endif

    match = _thread_id_map.find( native_thread_id );
    if (match == _thread_id_map.end( ))
    {
	uint32_t thread_id;
	if (_number_of_threads <= _thread_id_map.size( ))
	    throw TooManyThreads( );
	thread_id = __sync_fetch_and_add( &_thread_counter, 1 );
	while ( !__sync_bool_compare_and_swap( &_tid_map_spin_lock, 0, 1 ) );
	_thread_id_map[ native_thread_id ] = thread_id;
	__sync_bool_compare_and_swap( &_tid_map_spin_lock, 1, 0 );
	return thread_id;
    }

    return (*match).second; // Return the found value.
}


namespace read_parsers
{


IStreamReader::
IStreamReader( )
{
    
    at_eos	    = false;

}


RawStreamReader::
RawStreamReader( int const fd )
{

    if (0 > fd) throw InvalidStreamBuffer( );
    _stream_handle    = fd;

}


GzStreamReader::
GzStreamReader( int const fd )
{

    if (0 > fd) throw InvalidStreamBuffer( );
    _stream_handle    = gzdopen( fd, "rb" );
    if (NULL == _stream_handle) throw InvalidStreamBuffer( );

}


Bz2StreamReader::
Bz2StreamReader( int const fd )
{

    if (0 > fd) throw InvalidStreamBuffer( );
    if (NULL == (_stream_handle = fdopen( fd, "r" )))
	throw InvalidStreamBuffer( );

    _block_handle = NULL;

}


IStreamReader::
~IStreamReader( )
{ }


RawStreamReader::
~RawStreamReader( )
{
    
    if (0 <= _stream_handle) close( _stream_handle );
    _stream_handle = -1;

}


GzStreamReader::
~GzStreamReader( )
{
    
    if (NULL != _stream_handle) gzclose( _stream_handle );
    _stream_handle = NULL;

}


Bz2StreamReader::
~Bz2StreamReader( )
{
    int		bz2_error	= BZ_OK;
    
    if (NULL != _block_handle) BZ2_bzReadClose( &bz2_error, _block_handle );
    _block_handle = NULL;
    if (NULL != _stream_handle) fclose( _stream_handle );
    _stream_handle = NULL;

}


bool const
IStreamReader::
is_at_end_of_stream( ) const
{
    return at_eos;
}


uint64_t const
RawStreamReader::
read_into_cache( uint8_t * const cache, uint64_t const cache_size )
{
    ssize_t	nbread		    = 0;
    uint64_t	nbread_total	    = 0;

    assert (NULL != cache);
    if (0 == cache_size) return 0;

    for ( uint64_t nbrem = cache_size; (0 < nbrem) && !at_eos; nbrem -= nbread )
    {
	nbread =
	read(
	    _stream_handle, 
	    cache + nbread_total,
	    (size_t)( nbrem > SSIZE_MAX ? SSIZE_MAX : nbrem )
	);
	if (-1 == nbread) throw StreamReadError( );
	at_eos = !nbread;
	nbread_total += nbread;
    }

    return nbread_total;
}


uint64_t const
GzStreamReader::
read_into_cache( uint8_t * const cache, uint64_t const cache_size )
{
    int		nbread		    = 0;
    uint64_t	nbread_total	    = 0;

    assert (NULL != cache);
    if (0 == cache_size) return 0;

    for ( uint64_t nbrem = cache_size; (0 < nbrem) && !at_eos; nbrem -= nbread )
    {
	nbread =
	gzread(
	    _stream_handle, 
	    cache + nbread_total,
	    (unsigned int)( nbrem > INT_MAX ? INT_MAX : nbrem )
	);
	if (-1 == nbread) throw StreamReadError( );
	at_eos = !nbread;
	nbread_total += nbread;
    }

    return nbread_total;
}


uint64_t const
Bz2StreamReader::
read_into_cache( uint8_t * const cache, uint64_t const cache_size )
{
    int		bz2_error	    = BZ_OK;
    bool	block_complete	    = false;
    uint8_t	bz2_unused[ BZ_MAX_UNUSED ];
    uint8_t *	bz2_unused_temp	    = NULL;
    int		bz2_unused_nbread   = 0;
    int		nbread		    = 0;
    uint64_t	nbread_total	    = 0;

    assert (NULL != cache);
    if (0 == cache_size) return 0;

    for ( uint64_t nbrem = cache_size; (0 < nbrem) && !at_eos; nbrem -= nbread )
    {

	if (NULL == _block_handle)
	{
	    _block_handle = 
	    BZ2_bzReadOpen( 
		&bz2_error, 
		_stream_handle, 
		0, 0, 
		bz2_unused, bz2_unused_nbread
	    );
	    if (BZ_OK != bz2_error) throw InvalidStreamBuffer( );
	}

	nbread =
	BZ2_bzRead(
	    &bz2_error, 
	    _block_handle, 
	    cache + nbread_total, 
	    (int)( nbrem > INT_MAX ? INT_MAX : nbrem )
	);
	switch (bz2_error)
	{
	    
	    case BZ_STREAM_END: block_complete = true;
	    case BZ_OK:
		nbread_total += nbread;
		break;

	    // TODO: Inject BZ2 error code or error string into exception.
	    default: throw StreamReadError( );

	}

	if (block_complete)
	{
	    BZ2_bzReadGetUnused(
		&bz2_error,
		_block_handle,
		(void **)&bz2_unused_temp, &bz2_unused_nbread
	    );
	    // TODO: Inject BZ2 error code or error string into exception.
	    if (BZ_OK != bz2_error) throw StreamReadError( );
	    for (int i = 0; i < bz2_unused_nbread; ++i)
		bz2_unused[ i ] = bz2_unused_temp[ i ];

	    BZ2_bzReadClose( &bz2_error, _block_handle );
	    _block_handle = NULL;
	    if (feof( _stream_handle )) at_eos = true;
	    block_complete = false;
	}

    } // loop to fill cache from disk

    return nbread_total;
}


CacheManager::
CacheManager(
    IStreamReader & stream_reader,
    uint32_t const  number_of_threads,
    uint64_t const  cache_size,
    uint8_t const   trace_level
)
:   _trace_level( trace_level),
    _stream_reader( stream_reader ),
    _number_of_threads( number_of_threads ),
    _thread_id_map( ThreadIDMap( number_of_threads ) ),
    _segment_ref_count( 0 ),
    _segment_to_fill( 0 ),
    _fill_counter( 0 )
{

    if (cache_size < number_of_threads)	throw InvalidCacheSizeRequested( );
    _segment_size	= cache_size / number_of_threads;
    _segments		= new CacheSegment *[ number_of_threads ];
    for (uint32_t i = 0; i < number_of_threads; ++i) _segments[ i ] = NULL;

}


CacheManager::
~CacheManager( )
{

    for (uint32_t i = 0; i < _number_of_threads; ++i)
    {
	if (NULL != _segments[ i ])
	{
	    delete _segments[ i ];
	    _segments[ i ]	= NULL;
	}
    }
    delete [ ] _segments;
    _segments		= NULL;

}


CacheManager:: CacheSegment::
CacheSegment(
    uint32_t const  thread_id,
    uint64_t const  size,
    uint8_t const   trace_level
)
:   thread_id( thread_id ),
    size( size ),
    memory( new uint8_t[ size ] ),
    cursor( 0 ),
    cursor_in_sa_buffer( false ),
    fill_id( 0 ),
    _trace_level( trace_level )
{

    _sa_buffer_avail	= false;
    sa_buffer_size	= 0;

    if (trace_level)
    {
	char	    tfile_name[ 1024 ];
	sprintf(
	    tfile_name,
	    "cmgr_trace-thread_%lu",
	    (unsigned long int)thread_id
	);
	trace_file_handle  = fopen( tfile_name, "w" );
	if (NULL == trace_file_handle)
	    // TODO: Throw exception.
	    ;
	fprintf(
	    trace_file_handle,
	    "Trace of thread %lu started.\n",
	    (unsigned long int)thread_id
	);
	fflush( trace_file_handle );
    }
    else trace_file_handle = NULL;

    _nsecs_waiting_to_set_sa_buffer	= 0;
    _nsecs_waiting_to_acquire_sa_buffer	= 0;
    _nsecs_waiting_to_fill_from_stream	= 0;
    _nsecs_reading_from_stream		= 0;
    _nsecs_in_final_sync_barrier	= 0;

    avail		= true;
}


CacheManager:: CacheSegment::
~CacheSegment( )
{

    avail		= false;
    _sa_buffer_avail	= false;
    sa_buffer_size	= 0;
    size		= 0;

    delete [ ] memory;
    memory		= NULL;

    if (_trace_level && (NULL != trace_file_handle))
    {
	fprintf(
	    trace_file_handle,
	    "Waited to set my setaside buffer for %.9f seconds.\n" \
	    "Waited to get higher setaside buffer for %.9f seconds.\n" \
	    "Waited to fill from stream for %.9f seconds.\n" \
	    "Read from stream for %.9f seconds.\n" \
	    "Waited in final synchronization barrier for %.9f seconds.\n",
	    ((double)_nsecs_waiting_to_set_sa_buffer / 1000000000),
	    ((double)_nsecs_waiting_to_acquire_sa_buffer / 1000000000),
	    ((double)_nsecs_waiting_to_fill_from_stream / 1000000000),
	    ((double)_nsecs_reading_from_stream / 1000000000),
	    ((double)_nsecs_in_final_sync_barrier / 1000000000)
	);
    }

    if (_trace_level && (NULL != trace_file_handle))
	fclose( trace_file_handle );

}


bool const
CacheManager::
has_more_data( )
{
    CacheSegment &	segment		= _get_segment( );
    timespec		tspec_start;
    timespec		tspec_stop;

    // Return true immediately, if segment can provide more data.
    if (segment.avail) return true;

    if (1 < _trace_level)
    {
	fprintf(
	    segment.trace_file_handle,
	    "Before 'has_more_data' synchronization barrier.\n"
	);
	fflush( segment.trace_file_handle );
    }

    // Block indefinitely, if some other segment can provide more data.
    // (This is a synchronization barrier.)
    clock_gettime( CLOCK_THREAD_CPUTIME_ID, &tspec_start );
sync_barrier:
    for (uint64_t i = 0; _segment_ref_count; ++i)
    {
	if ((2 < _trace_level) && (0 == i % 100000000))
	{
	    fprintf(
		segment.trace_file_handle,
		"Waited in synchronization barrier for %llu iterations.\n",
		(unsigned long long int)i
	    );
	    fflush( segment.trace_file_handle );
	}
    }

    // Return false, if no segment can provide more data.
    if (!_get_segment_ref_count_ATOMIC( ))
    {
	clock_gettime( CLOCK_THREAD_CPUTIME_ID, &tspec_stop );
	segment._nsecs_in_final_sync_barrier =
	_timespec_diff_in_nsecs( tspec_start, tspec_stop );
	if (1 < _trace_level)
	{
	    fprintf(
		segment.trace_file_handle,
		"After 'has_more_data' synchronization barrier.\n"
	    );
	    fflush( segment.trace_file_handle );
	}
	return false;
    }
    // If we somehow got here and there are still active segments,
    // then go back to waiting.
    else goto sync_barrier;
}


uint8_t const
CacheManager::
get_byte( )
{
    CacheSegment &	segment		= _get_segment( );
    uint8_t *		memory		= NULL;

    if (!segment.avail) throw CacheSegmentUnavailable( );

    _perform_segment_maintenance( segment );
    if (segment.cursor_in_sa_buffer)
    {
	CacheSegment &
	hsegment    = _get_segment( true );
	memory	    = hsegment.memory;
    }
    else memory	    = segment.memory;

    return memory[ segment.cursor++ ];
}


uint64_t const
CacheManager::
get_bytes( uint8_t * const buffer, uint64_t buffer_len )
{
    CacheSegment &	segment		= _get_segment( );
    uint64_t		nbcopied	= 0;
    uint64_t		nbcopied_total	= 0;
    uint8_t *		memory		= NULL;
    uint64_t		size		= 0;

    if (!segment.avail) throw CacheSegmentUnavailable( );

    for (uint64_t nbrem = buffer_len; (nbrem > 0); nbrem -= nbcopied)
    {
	_perform_segment_maintenance( segment );
	if (segment.cursor_in_sa_buffer)
	{
	    CacheSegment &
	    hsegment	    = _get_segment( true );
	    memory	    = hsegment.memory;
	    size	    = hsegment.sa_buffer_size;
	}
	else
	{
	    memory	    = segment.memory;
	    size	    = segment.size;
	}

	nbcopied = MIN( nbrem, size - segment.cursor );
	if (0 == nbcopied) break;
	memcpy( buffer, memory + segment.cursor, nbcopied );
	segment.cursor += nbcopied;

	nbcopied_total += nbcopied;
    }

    return nbcopied_total;
}


uint64_t const
CacheManager::
whereis_cursor( )
{ return _get_segment( ).cursor; }


void
CacheManager::
split_at( uint64_t const pos )
{

    CacheSegment &	segment		= _get_segment( );
    timespec		tspec_start;
    timespec		tspec_stop;

    if (1 < _trace_level)
    {
	fprintf(
	    segment.trace_file_handle,
	    "Splitting off setaside buffer at byte %llu.\n",
	    (unsigned long long int)pos
	);
	fflush( segment.trace_file_handle );
    }

    // Wait until the lower segment has consumed the setaside buffer.
    clock_gettime( CLOCK_THREAD_CPUTIME_ID, &tspec_start );
wait_for_sa_buffer:
    for (uint64_t i = 0; segment.get_sa_buffer_avail( ); ++i)
    {
	if ((2 < _trace_level) && (0 == i % 100000000))
	{
	    fprintf(
		segment.trace_file_handle,
		"Waited to set setaside buffer for %llu iterations.\n",
		(unsigned long long int)i
	    );
	    fflush( segment.trace_file_handle );
	}
    }

    // If we get here but are not ready to proceed,
    // then go back and wait some more.
    if (segment.get_sa_buffer_avail_ATOMIC( ))
	goto wait_for_sa_buffer;

    clock_gettime( CLOCK_THREAD_CPUTIME_ID, &tspec_stop );
    segment._nsecs_waiting_to_set_sa_buffer +=
    _timespec_diff_in_nsecs( tspec_start, tspec_stop );

    // Setup the setaside buffer.
    segment.sa_buffer_size = pos;
    segment.set_sa_buffer_avail_ATOMIC( true );

    if (1 < _trace_level)
    {
	fprintf(
	    segment.trace_file_handle,
	    "Finished 'split_at'.\n"
	);
	fflush( segment.trace_file_handle );
    }

}


uint64_t const
CacheManager::
get_fill_id( )
{ return _get_segment( ).fill_id; }


bool const
CacheManager::
_in_sa_buffer( )
{ return _get_segment( ).cursor_in_sa_buffer; }


bool const
CacheManager::
_sa_buffer_avail( )
{ return _get_segment( ).get_sa_buffer_avail_ATOMIC( ); }


void
CacheManager::
_perform_segment_maintenance( CacheSegment & segment )
{

    assert( segment.avail );

    CacheSegment &	hsegment	    = _get_segment( true );
    timespec		tspec_start;
    timespec		tspec_stop;

    // If at end of segment and not already in setaside buffer, 
    // then jump into setaside buffer from higher segment.
    if (!segment.cursor_in_sa_buffer && (segment.cursor == segment.size))
    {

	// If there is only 1 thread, 
	// then force setaside buffer to be available.
	if (segment.thread_id == hsegment.thread_id)
	    hsegment.set_sa_buffer_avail_ATOMIC( true );

	// Wait while higher segment is available 
	// and its setaside buffer is not ready for consumption.
	clock_gettime( CLOCK_THREAD_CPUTIME_ID, &tspec_start );
wait_for_sa_buffer:
	for (	uint64_t i = 0;
		hsegment.avail && !hsegment.get_sa_buffer_avail( );
		++i)
	{
	    if ((2 < _trace_level) && (0 == i % 100000000))
	    {
		fprintf(
		    segment.trace_file_handle,
		    "Waited to get setaside buffer for %llu iterations.\n",
		    (unsigned long long int)i
		);
		fflush( segment.trace_file_handle );
	    }
	}

	// Atomically test that the setaside buffer is available.
	// If so, then jump into it.
	if (hsegment.get_sa_buffer_avail_ATOMIC( ))
	{
	    clock_gettime( CLOCK_THREAD_CPUTIME_ID, &tspec_stop );
	    segment._nsecs_waiting_to_acquire_sa_buffer +=
	    _timespec_diff_in_nsecs( tspec_start, tspec_stop );
	    segment.cursor_in_sa_buffer	    = true;
	    segment.cursor		    = 0;
	    if (1 < _trace_level)
	    {
		fprintf(
		    segment.trace_file_handle,
		    "Jumped into setaside buffer.\n"
		);
		fflush( segment.trace_file_handle );
	    }
	}

	// If the higher segment is no longer available,
	// and its setaside buffer was never set,
	// then jump into a dummy buffer with no bytes remaining.
	else if (!hsegment.avail)
	{
	    segment.cursor_in_sa_buffer	    = true;
	    segment.cursor		    = hsegment.sa_buffer_size;
	}

	// If we somehow got here and shouldn't have, 
	// then go back and wait some more.
	else goto wait_for_sa_buffer;

    } // jump into setaside buffer

    // If at end of setaside buffer...
    if (    segment.cursor_in_sa_buffer
	&&  (segment.cursor == hsegment.sa_buffer_size))
    {
	
	// Jump out of setaside buffer and reset it.
	segment.cursor_in_sa_buffer	= false;
	segment.cursor			= 0;
	hsegment.sa_buffer_size		= 0;
	hsegment.set_sa_buffer_avail_ATOMIC( false );
	if (1 < _trace_level)
	{
	    fprintf(
		segment.trace_file_handle,
		"Jumped out of setaside buffer.\n"
	    );
	    fflush( segment.trace_file_handle );
	}

	// Jump past end of setaside buffer
	// so as not to clobber what the lower segment will want to use.
	if (segment.get_sa_buffer_avail_ATOMIC( ))
	    segment.cursor		= segment.sa_buffer_size;
	
	_fill_segment_from_stream( segment );
	
    } // refill or mark unavailable

}


inline
bool const
CacheManager::
_check_segment_to_fill_ATOMIC( uint32_t const thread_id )
{
    uint32_t	segment_idx	= 
    __sync_and_and_fetch( &_segment_to_fill, (uint32_t)0xffffffff );
    return (thread_id == segment_idx);
}


inline
void
CacheManager::
_select_segment_to_fill_ATOMIC( )
{
    uint32_t	segment_idx =
    __sync_add_and_fetch( &_segment_to_fill, 1 );
    if (_number_of_threads == segment_idx)
	__sync_bool_compare_and_swap(
	    &_segment_to_fill, _number_of_threads, 0
	);
}


inline
CacheManager:: CacheSegment &
CacheManager::
_get_segment( bool const higher )
{
    uint32_t	    thread_id		= _thread_id_map.get_thread_id( );
    CacheSegment *  segment_PTR		= NULL;

    assert( NULL != _segments );

    // If referring to a segment to snoop,
    // then index is for the thread with the next higher ID.
    if (higher) thread_id = ((thread_id + 1) % _number_of_threads);

    segment_PTR	    = _segments[ thread_id ];
    // TODO: Protect with a mutex in case another thread 
    //	     is trying to create the segment at the same time.
    if (NULL == segment_PTR)
    {
	_segments[ thread_id ]	    = new CacheSegment(
	    thread_id, _segment_size, _trace_level
	);
	segment_PTR		    = _segments[ thread_id ];
	_increment_segment_ref_count_ATOMIC( );
	_fill_segment_from_stream( *segment_PTR );
    }

    return *segment_PTR;
}


inline
void
CacheManager::
_fill_segment_from_stream( CacheSegment & segment )
{
    timespec tspec_start;
    timespec tspec_stop;

    // Wait while segment not selected and not end of stream.
    clock_gettime( CLOCK_THREAD_CPUTIME_ID, &tspec_start );
wait_to_fill:
    for (   uint64_t i = 0;
	    !_stream_reader.is_at_end_of_stream( )
	&&  (_segment_to_fill != segment.thread_id);
	    ++i)
    {
	if ((2 < _trace_level) && (0 == i % 100000000))
	{
	    fprintf(
		segment.trace_file_handle,
		"Waited to fill segment for %llu iterations.\n",
		(unsigned long long int)i
	    );
	    fflush( segment.trace_file_handle );
	}
    }
    clock_gettime( CLOCK_THREAD_CPUTIME_ID, &tspec_stop );
    segment._nsecs_waiting_to_fill_from_stream +=
    _timespec_diff_in_nsecs( tspec_start, tspec_stop );

    // If at end of stream, then mark segment unavailable.
    if (_stream_reader.is_at_end_of_stream( ))
    {
	if (1 < _trace_level)
	{
	    fprintf(
		segment.trace_file_handle,
		"At end of input stream.\n"
	    );
	    fflush( segment.trace_file_handle );
	}
	segment.avail		= false;
	_decrement_segment_ref_count_ATOMIC( );
    }

    // Else, refill the segment.
    else if (_check_segment_to_fill_ATOMIC( segment.thread_id ))
    {
	clock_gettime( CLOCK_THREAD_CPUTIME_ID, &tspec_start );
	segment.size		=
	    segment.cursor
	+   _stream_reader.read_into_cache(
		segment.memory + segment.cursor,
		_segment_size - segment.cursor
	    );
	clock_gettime( CLOCK_THREAD_CPUTIME_ID, &tspec_stop );
	segment.fill_id = _fill_counter++;
	_select_segment_to_fill_ATOMIC( );
	segment._nsecs_reading_from_stream +=
	_timespec_diff_in_nsecs( tspec_start, tspec_stop );
	if (2 < _trace_level)
	{
	    fprintf(
		segment.trace_file_handle,
		"Read %llu bytes into segment.\n",
		(unsigned long long int)segment.size
	    );
	    fflush( segment.trace_file_handle );
	}
    }

    // If we somehow get here, then go back and wait some more.
    else goto wait_to_fill;

}


inline
void
CacheManager::
_increment_segment_ref_count_ATOMIC( )
{
    __sync_add_and_fetch( &_segment_ref_count, 1 );
}


inline
void
CacheManager::
_decrement_segment_ref_count_ATOMIC( )
{
    __sync_sub_and_fetch( &_segment_ref_count, 1 );
}

inline
uint32_t const
CacheManager::
_get_segment_ref_count_ATOMIC( )
{
    return __sync_and_and_fetch( &_segment_ref_count, (uint32_t)0xffffffff );
}


inline
bool
CacheManager:: CacheSegment::
get_sa_buffer_avail( ) const
{
    return _sa_buffer_avail;
}


inline
bool
CacheManager:: CacheSegment::
get_sa_buffer_avail_ATOMIC( )
{
    return __sync_and_and_fetch( &_sa_buffer_avail, 1 );
}


inline
void
CacheManager:: CacheSegment::
set_sa_buffer_avail_ATOMIC( bool const avail )
{
    __sync_bool_compare_and_swap( &_sa_buffer_avail, !avail, avail );
}


IParser:: IParser * const
IParser::
get_parser(
    std:: string const &    ifile_name,
    uint32_t const	    number_of_threads,
    uint64_t const	    cache_size,
    uint8_t const	    trace_level
)
{
    // TODO: Replace file extension detection with header magic detection.

    IStreamReader * stream_reader   = NULL;
    IParser *	    parser	    = NULL;

    int		    ifile_handle    = open( ifile_name.c_str( ), O_RDONLY );
    if (-1 == ifile_handle)
	// TODO: Raise exception.
	;
    std:: string    ext	    = "";
    std:: string    ifile_name_chopped( ifile_name );
    size_t	    ext_pos = ifile_name.find_last_of( "." );
    bool	    rechop  = false;

    if (0 < ext_pos)
    {
	ext		    = ifile_name.substr( ext_pos + 1 );
	ifile_name_chopped  = ifile_name.substr( 0, ext_pos );
    }

    if	    ("gz" == ext)
    {
	stream_reader	= new GzStreamReader( ifile_handle );
	rechop		= true;
    }
    else if ("bz2" == ext)
    {
	stream_reader	= new Bz2StreamReader( ifile_handle );
	rechop		= true;
    }
    else
	stream_reader	= new RawStreamReader( ifile_handle );

    if (rechop)
    {
	ext_pos		    = ifile_name_chopped.find_last_of( "." );
	ext		    = ifile_name_chopped.substr( ext_pos + 1 );
	ifile_name_chopped  = ifile_name_chopped.substr( 0, ext_pos );
    }

    if (("fq" == ext) || ("fastq" == ext))
	parser =
	new FastqParser(
	    *stream_reader,
	    number_of_threads,
	    cache_size,
	    trace_level
	);
    else
	parser =
	new FastaParser(
	    *stream_reader,
	    number_of_threads,
	    cache_size,
	    trace_level
	);

    return parser;
}


IParser::
IParser(
    IStreamReader &  stream_reader,
    uint32_t const  number_of_threads,
    uint64_t const  cache_size,
    uint8_t const   trace_level
)
:   _trace_level( trace_level ),
    _cache_manager(
	CacheManager(
	    stream_reader, number_of_threads, cache_size, trace_level
	)
    ),
    _thread_id_map( ThreadIDMap( number_of_threads ) )
{
    
    _states  = new ParserState *[ number_of_threads ];
    for (uint32_t i = 0; i < number_of_threads; ++i) _states[ i ] = NULL;

}


IParser::
~IParser( )
{ }


bool
IParser::
is_complete( )
{ return !_cache_manager.has_more_data( ); }


inline
void
IParser::
_copy_line( )
{
    ParserState &   state	= _get_state( );
    uint8_t (&	    buffer)[ ParserState:: BUFFER_SIZE + 1 ]
				= state.buffer;
    uint64_t &	    pos		= state.buffer_pos;
    uint64_t &	    rem		= state.buffer_rem;
    std:: string &  line	= state.line;
    uint64_t	    i		= 0;
    bool	    hit		= false;

    line.clear( );

    while( true )
    {

	for (i = 0; (i < rem) && ('\n' != buffer[ pos + i ]); i++);
	if (0 < rem)
	{
	    if (i < rem)
	    {
		buffer[ pos + i ]   = '\0';
		hit		    = true;
	    }
	    line += (char const *)(buffer + pos);
	    rem -= (i + 1); pos += (i + 1);
	    if (hit) break;
	}
	
	if (_cache_manager.has_more_data( ))
	{
	    rem = _cache_manager.get_bytes( buffer, ParserState:: BUFFER_SIZE );
	    pos = 0;
	}
	else break;

    }

}


inline
IParser:: ParserState &
IParser::
_get_state( )
{
    uint32_t	    thread_id		= _thread_id_map.get_thread_id( );
    ParserState *   state_PTR		= NULL;

    assert( NULL != _states );

    state_PTR	    = _states[ thread_id ];
    if (NULL == state_PTR)
    {
	_states[ thread_id ]	    = new ParserState( );
	state_PTR		    = _states[ thread_id ];
    }

    return *state_PTR;
}


IParser:: ParserState::
ParserState( )
:   at_start( true ),
    need_new_line( true ),
    buffer_pos( 0 ),
    buffer_rem( 0 )
{ memset( buffer, 0, BUFFER_SIZE + 1 ); }


IParser:: ParserState::
~ParserState( )
{ }


FastaParser::
FastaParser(
    IStreamReader &  stream_reader,
    uint32_t const  number_of_threads,
    uint64_t const  cache_size,
    uint8_t const   trace_level
)
: IParser( stream_reader, number_of_threads, cache_size, trace_level )
{ }


FastqParser::
FastqParser(
    IStreamReader &  stream_reader,
    uint32_t const  number_of_threads,
    uint64_t const  cache_size,
    uint8_t const   trace_level
)
: IParser( stream_reader, number_of_threads, cache_size, trace_level )
{ }


FastaParser::
~FastaParser( )
{ }


FastqParser::
~FastqParser( )
{ }


Read
FastaParser::
get_next_read( )
{
    
    ParserState &   parser_state    = _get_state( );
    Read	    the_read;
    
    while (_cache_manager.has_more_data( ))
    {
	the_read.name	= "";
	the_read.seq	= "";

	if (parser_state.need_new_line) _copy_line( );
	parser_state.need_new_line = true;

	// If a cache segment fill occurred during line retrieval, 
	// then flag the event.
	if (!parser_state.at_start)
	    parser_state.at_start =
	    (parser_state.fill_id != _cache_manager.get_fill_id( ));

	parser_state.fill_id = _cache_manager.get_fill_id( );

	// If at start of file, then error on garbage.
	// Else, skip forward to next read boundary.
	if ('>' != parser_state.line[ 0 ])
	{
	    if (parser_state.at_start && (0 == parser_state.fill_id))
		// TODO: Throw exception.
		;
	    continue;
	}

	// Skip over an unmatched second part of a paired read,
	// when at the beginning of a new fill.
	if (	parser_state.at_start
	    &&	(    (parser_state.line.length( ) - 2)
		 ==  parser_state.line.rfind( "/2" )))
	{
	    while (_cache_manager.has_more_data( ))
	    {
		_copy_line( );
		if ('>' == parser_state.line[ 0 ]) break;
	    }
	}

	// Split off skipped over data into a setaside buffer.
	_cache_manager.split_at(
		_cache_manager.whereis_cursor( )
	    -   parser_state.line.length( )
	);

	parser_state.at_start = false;

	// Parse read.
	the_read.name		    = parser_state.line.substr( 1 );
	parser_state.need_new_line  = false;
	while (_cache_manager.has_more_data( ))
	{
	    _copy_line( );
	    if ('>' == parser_state.line[ 0 ]) break;
	    the_read.seq += parser_state.line;
	}

	// Discard invalid read.
	if (std:: string:: npos != the_read.seq.find_first_of( "Nn" )) continue;

	break;
    } // while invalid read

    return the_read;
}


Read
FastqParser::
get_next_read( )
{
    // TODO: Implement.
}


} // namespace read_parsers


} // namespace khmer

// vim: set ft=cpp sts=4 sw=4 tw=80:
