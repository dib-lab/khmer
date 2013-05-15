#include <cstdlib>
#include <cstring>
#include <cstdio>

#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#ifdef __linux__
#   include <sys/ioctl.h>
#   include <linux/fs.h>
#endif

#include "read_parsers.hh"


namespace khmer
{


namespace read_parsers
{


StreamReaderPerformanceMetrics::
StreamReaderPerformanceMetrics( )
:   IPerformanceMetrics( ),
    numbytes_read( 0 ),
    clock_nsecs_reading( 0 ),
    cpu_nsecs_reading( 0 )
{ }


StreamReaderPerformanceMetrics::
~StreamReaderPerformanceMetrics( )
{ }


void
StreamReaderPerformanceMetrics::
accumulate_timer_deltas( uint32_t metrics_key )
{

    switch (metrics_key)
    {
    case MKEY_TIME_READING:
	clock_nsecs_reading +=
	_timespec_diff_in_nsecs( _temp_clock_start, _temp_clock_stop );
	cpu_nsecs_reading   +=
	_timespec_diff_in_nsecs( _temp_cpu_start, _temp_cpu_stop );
	break;
    default: throw InvalidPerformanceMetricsKey( );
    }

}


IStreamReader::
IStreamReader( )
:   pmetrics( StreamReaderPerformanceMetrics( ) ),
    _alignment( 0 ),
    _max_aligned( SSIZE_MAX ),
    _at_eos( false )
{ }


RawStreamReader::
RawStreamReader( int const fd, size_t const alignment )
: IStreamReader( )
{

    if (0 > fd) throw InvalidStreamBuffer( );

#ifdef __linux__
    if (alignment)
    {
	_alignment	= alignment;
	_max_aligned	= _alignment * (SSIZE_MAX / _alignment);
    }
#endif
	
    _stream_handle    = fd;

}


GzStreamReader::
GzStreamReader( int const fd )
: IStreamReader( )
{

    if (0 > fd) throw InvalidStreamBuffer( );
    _stream_handle    = gzdopen( fd, "rb" );
    if (NULL == _stream_handle) throw InvalidStreamBuffer( );

}


Bz2StreamReader::
Bz2StreamReader( int const fd )
: IStreamReader( )
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
{ return _at_eos; }


size_t const
IStreamReader::
get_memory_alignment( ) const
{ return _alignment; }


uint64_t const
RawStreamReader::
read_into_cache( uint8_t * const cache, uint64_t const cache_size )
{
    uint64_t	cache_size_adjusted = cache_size;
    ssize_t	nbread		    = 0;
    uint64_t	nbread_total	    = 0;

    assert (NULL != cache);
    if (0 == cache_size) return 0;

#ifdef __linux__
    if (_alignment && (0 != (cache_size % _alignment)))
	cache_size_adjusted = _alignment * ((cache_size / _alignment) + 1);
#endif
    for (uint64_t nbrem = cache_size_adjusted;
	 (0 < nbrem) && !_at_eos;
	 nbrem -= nbread)
    {
#ifdef WITH_INTERNAL_METRICS
	pmetrics.start_timers( );
#endif
	nbread =
	read(
	    _stream_handle, 
	    cache + nbread_total,
	    (size_t)(nbrem > _max_aligned ? _max_aligned : nbrem )
	);
#ifdef WITH_INTERNAL_METRICS
	pmetrics.stop_timers( );
	pmetrics.accumulate_timer_deltas(
	    (uint32_t)StreamReaderPerformanceMetrics:: MKEY_TIME_READING
	);
#endif
	if (-1 == nbread) throw StreamReadError( );
	_at_eos = !nbread;
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

    for (uint64_t nbrem = cache_size; (0 < nbrem) && !_at_eos; nbrem -= nbread)
    {
#ifdef WITH_INTERNAL_METRICS
	pmetrics.start_timers( );
#endif
	nbread =
	gzread(
	    _stream_handle, 
	    cache + nbread_total,
	    (unsigned int)( nbrem > INT_MAX ? INT_MAX : nbrem )
	);
#ifdef WITH_INTERNAL_METRICS
	pmetrics.stop_timers( );
	pmetrics.accumulate_timer_deltas(
	    (uint32_t)StreamReaderPerformanceMetrics:: MKEY_TIME_READING
	);
#endif
	if (-1 == nbread) throw StreamReadError( );
	_at_eos = !nbread;
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

    for (uint64_t nbrem = cache_size; (0 < nbrem) && !_at_eos; nbrem -= nbread)
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

#ifdef WITH_INTERNAL_METRICS
	pmetrics.start_timers( );
#endif
	nbread =
	BZ2_bzRead(
	    &bz2_error, 
	    _block_handle, 
	    cache + nbread_total, 
	    (int)( nbrem > INT_MAX ? INT_MAX : nbrem )
	);
#ifdef WITH_INTERNAL_METRICS
	pmetrics.stop_timers( );
	pmetrics.accumulate_timer_deltas(
	    (uint32_t)StreamReaderPerformanceMetrics:: MKEY_TIME_READING
	);
#endif
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
	    if (feof( _stream_handle )) _at_eos = true;
	    block_complete = false;
	}

    } // loop to fill cache from disk

    return nbread_total;
}


CacheSegmentPerformanceMetrics::
CacheSegmentPerformanceMetrics( )
:   IPerformanceMetrics( ),
    numbytes_filled_from_stream( 0 ),
    numbytes_copied_from_sa_buffer( 0 ),
    numbytes_reserved_as_sa_buffer( 0 ),
    numbytes_copied_to_caller_buffer( 0 ),
    clock_nsecs_waiting_to_set_sa_buffer( 0 ),
    cpu_nsecs_waiting_to_set_sa_buffer( 0 ),
    clock_nsecs_waiting_to_get_sa_buffer( 0 ),
    cpu_nsecs_waiting_to_get_sa_buffer( 0 ),
    clock_nsecs_waiting_to_fill_from_stream( 0 ),
    cpu_nsecs_waiting_to_fill_from_stream( 0 ),
    clock_nsecs_filling_from_stream( 0 ),
    cpu_nsecs_filling_from_stream( 0 ),
    clock_nsecs_in_sync_barrier( 0 ),
    cpu_nsecs_in_sync_barrier( 0 ),
    _accumulated_count( 1 )
{ }


CacheSegmentPerformanceMetrics::
~CacheSegmentPerformanceMetrics( )
{ }


void
CacheSegmentPerformanceMetrics::
accumulate_timer_deltas( uint32_t metrics_key )
{
    
    switch (metrics_key)
    {
    case MKEY_TIME_WAITING_TO_SET_SA_BUFFER:
	clock_nsecs_waiting_to_set_sa_buffer +=
	_timespec_diff_in_nsecs( _temp_clock_start, _temp_clock_stop );
	cpu_nsecs_waiting_to_set_sa_buffer   +=
	_timespec_diff_in_nsecs( _temp_cpu_start, _temp_cpu_stop );
	break;
    case MKEY_TIME_WAITING_TO_GET_SA_BUFFER:
	clock_nsecs_waiting_to_get_sa_buffer +=
	_timespec_diff_in_nsecs( _temp_clock_start, _temp_clock_stop );
	cpu_nsecs_waiting_to_get_sa_buffer   +=
	_timespec_diff_in_nsecs( _temp_cpu_start, _temp_cpu_stop );
	break;
    case MKEY_TIME_WAITING_TO_FILL_FROM_STREAM:
	clock_nsecs_waiting_to_fill_from_stream +=
	_timespec_diff_in_nsecs( _temp_clock_start, _temp_clock_stop );
	cpu_nsecs_waiting_to_fill_from_stream   +=
	_timespec_diff_in_nsecs( _temp_cpu_start, _temp_cpu_stop );
	break;
    case MKEY_TIME_FILLING_FROM_STREAM:
	clock_nsecs_filling_from_stream +=
	_timespec_diff_in_nsecs( _temp_clock_start, _temp_clock_stop );
	cpu_nsecs_filling_from_stream   +=
	_timespec_diff_in_nsecs( _temp_cpu_start, _temp_cpu_stop );
	break;
    case MKEY_TIME_IN_SYNC_BARRIER:
	clock_nsecs_in_sync_barrier +=
	_timespec_diff_in_nsecs( _temp_clock_start, _temp_clock_stop );
	cpu_nsecs_in_sync_barrier   +=
	_timespec_diff_in_nsecs( _temp_cpu_start, _temp_cpu_stop );
	break;
    default: throw InvalidPerformanceMetricsKey( );
    }

}


void
CacheSegmentPerformanceMetrics::
accumulate_metrics( CacheSegmentPerformanceMetrics &source )
{

    numbytes_filled_from_stream		    +=
    source.numbytes_filled_from_stream;
    numbytes_copied_from_sa_buffer	    +=
    source.numbytes_copied_from_sa_buffer;
    numbytes_reserved_as_sa_buffer	    +=
    source.numbytes_reserved_as_sa_buffer;
    numbytes_copied_to_caller_buffer	    +=
    source.numbytes_copied_to_caller_buffer;
    clock_nsecs_waiting_to_set_sa_buffer    +=
    source.clock_nsecs_waiting_to_set_sa_buffer;
    cpu_nsecs_waiting_to_set_sa_buffer	    +=
    source.cpu_nsecs_waiting_to_set_sa_buffer;
    clock_nsecs_waiting_to_get_sa_buffer    +=
    source.clock_nsecs_waiting_to_get_sa_buffer;
    cpu_nsecs_waiting_to_get_sa_buffer	    +=
    source.cpu_nsecs_waiting_to_get_sa_buffer;
    clock_nsecs_waiting_to_fill_from_stream +=
    source.clock_nsecs_waiting_to_fill_from_stream;
    cpu_nsecs_waiting_to_fill_from_stream   +=
    source.cpu_nsecs_waiting_to_fill_from_stream;
    clock_nsecs_filling_from_stream	    +=
    source.clock_nsecs_filling_from_stream;
    cpu_nsecs_filling_from_stream	    +=
    source.cpu_nsecs_filling_from_stream;
    clock_nsecs_in_sync_barrier		    +=
    source.clock_nsecs_in_sync_barrier;
    cpu_nsecs_in_sync_barrier		    +=
    source.cpu_nsecs_in_sync_barrier;
    _accumulated_count			    +=
    source._accumulated_count;

}


CacheManager::
CacheManager(
    IStreamReader & stream_reader,
    uint32_t const  number_of_threads,
    uint64_t const  cache_size,
    uint8_t const   trace_level
)
:   _trace_level( trace_level ),
    _stream_reader( stream_reader ),
    _number_of_threads( number_of_threads ),
    _thread_id_map( ThreadIDMap( number_of_threads ) ),
    _segment_ref_count( 0 ),
    _segment_to_fill( 0 ),
    _fill_counter( 0 ),
    _ca_spin_lock( 0 )
{

    if (cache_size < number_of_threads)	throw InvalidCacheSizeRequested( );
    _segment_size	= cache_size / number_of_threads;
    _alignment		= stream_reader.get_memory_alignment( );
#ifdef __linux__
    if (_alignment)
    {
	_segment_size = _alignment * (_segment_size / _alignment);
	if (!_segment_size) _segment_size = _alignment;
    }
#endif
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
    size_t const    alignment,
    uint8_t const   trace_level
)
:   thread_id( thread_id ),
    size( size ),
    alignment( alignment ),
    cursor( 0 ),
    cursor_in_ca_buffer( false ),
    fill_id( 0 ),
    pmetrics( CacheSegmentPerformanceMetrics( ) ),
    trace_logger(
	TraceLogger(
	    trace_level, "cmgr-%lu.log", (unsigned long int)thread_id
	)
    )
{
    
#ifdef __linux__
    if (alignment)
    {
	if (!posix_memalign(
	    (void **)&memory, alignment, size * sizeof( uint8_t )
	))
	    throw std:: bad_alloc( );
    }
    else
#endif
	memory = new uint8_t[ size ];

    trace_logger(
	TraceLogger:: TLVL_INFO0, 
	"Trace of thread %lu started.\n", (unsigned long int)thread_id
    );

    avail		= true;

}


CacheManager:: CacheSegment::
~CacheSegment( )
{

    avail		= false;
    size		= 0;
    ca_buffer.clear( );

#ifdef __linux__
    if (alignment) free( memory );
    else
#endif
	delete [ ] memory;
    memory		= NULL;

}


bool const
CacheManager::
has_more_data( )
{
    CacheSegment &	segment		= _get_segment( );

    // Return true immediately, if segment can provide more data.
    if (segment.avail) return true;

    segment.trace_logger(
	TraceLogger:: TLVL_DEBUG1,
	"Before 'has_more_data' synchronization barrier.\n"
    );

    // Block indefinitely, if some other segment can provide more data.
    // (This is a synchronization barrier.)
#ifdef WITH_INTERNAL_METRICS
    segment.pmetrics.start_timers( );
#endif
sync_barrier:
    for (uint64_t i = 0; _segment_ref_count; ++i)
    {
	// TODO: Determine optimal period. (Probably arch-dependent.)
	if (0 == i % 100000)
	{
	    if (0 == i % 100000000)
		segment.trace_logger(
		    TraceLogger:: TLVL_DEBUG3,
		    "Waited in synchronization barrier for %llu iterations.\n",
		    (unsigned long long int)i
		);
	    // HACK: Occasionally issue a memory barrier to break things up.
	    _get_segment_ref_count_ATOMIC( );
	}
    }

    // Return false, if no segment can provide more data.
    if (!_get_segment_ref_count_ATOMIC( ))
    {
#ifdef WITH_INTERNAL_METRICS
	segment.pmetrics.stop_timers( );
	segment.pmetrics.accumulate_timer_deltas(
	    CacheSegmentPerformanceMetrics:: MKEY_TIME_IN_SYNC_BARRIER
	);
#endif
	segment.trace_logger(
	    TraceLogger:: TLVL_DEBUG1,
	    "After 'has_more_data' synchronization barrier.\n"
	);

	return false;
    }
    // If we somehow got here and there are still active segments,
    // then go back to waiting.
    else goto sync_barrier;
}


uint64_t const
CacheManager::
get_bytes( uint8_t * const buffer, uint64_t buffer_len )
{
    CacheSegment	&segment	= _get_segment( );
    uint64_t		nbcopied	= 0;
    uint64_t		nbcopied_total	= 0;
    uint8_t *		memory		= NULL;
    uint64_t		size		= 0;
    bool		in_sa_buffer	= false;
    TraceLogger		&trace_logger	= segment.trace_logger;

    if (!segment.avail) throw CacheSegmentUnavailable( );

    for (uint64_t nbrem = buffer_len; (0 < nbrem); nbrem -= nbcopied)
    {

	_perform_segment_maintenance( segment );

	if (segment.cursor_in_ca_buffer)
	{
	    memory	    = (uint8_t *)segment.ca_buffer.c_str( );
	    size	    = (uint64_t)segment.ca_buffer.length( );
	    in_sa_buffer    = true;
	}
	else
	{
	    memory	    = segment.memory;
	    size	    = segment.size;
	    if (!segment.avail) break;
	    if (in_sa_buffer) in_sa_buffer = false;
	}

	nbcopied = MIN( nbrem, size - segment.cursor );
	memcpy( buffer + nbcopied_total, memory + segment.cursor, nbcopied );
	segment.cursor += nbcopied;

	trace_logger(
	    TraceLogger:: TLVL_DEBUG8,
	    "get_bytes: Copied %llu bytes from %s.\n",
	    (unsigned long long int)nbcopied,
	    in_sa_buffer ? "setaside buffer" : "cache segment"
	);

#ifdef WITH_INTERNAL_METRICS
	segment.pmetrics.numbytes_copied_to_caller_buffer += nbcopied;
	if (in_sa_buffer)
	    segment.pmetrics.numbytes_copied_from_sa_buffer += nbcopied;
#endif
	nbcopied_total += nbcopied;
    }

    return nbcopied_total;
}


uint64_t const
CacheManager::
whereis_cursor( void )
{ return _get_segment( ).cursor; }


bool const
CacheManager::
is_cursor_in_ca_buffer( void )
{ return _get_segment( ).cursor_in_ca_buffer; }


// TODO? Change type of 'pos' to 'size_t'.
void
CacheManager::
split_at( uint64_t const pos )
{

    CacheSegment &	segment		= _get_segment( );

    if (1 == _number_of_threads) return;

    segment.trace_logger(
	TraceLogger:: TLVL_DEBUG2,
	"Creating copyaside buffer for fill ID %llu up to byte %llu....\n",
	(unsigned long long int)segment.fill_id, (unsigned long long int)pos
    );

#ifdef WITH_INTERNAL_METRICS
    segment.pmetrics.start_timers( );
#endif
    // Acquire copyaside buffers spinlock.
    for (   uint64_t i = 0;
	    !__sync_bool_compare_and_swap( &_ca_spin_lock, 0, 1 );
	    ++i
	)
    {
	if (0 == i % 100000000)
	{
	    segment.trace_logger(
		TraceLogger:: TLVL_DEBUG3,
		"Waited to acquire copyaside buffers spinlock " \
		"for %llu iterations. [write buffer]\n",
		(unsigned long long int)i
	    );
	    segment.trace_logger(
		TraceLogger:: TLVL_DEBUG4,
		"\tSpinlock is probably %s.\n",
		_ca_spin_lock ? "set" : "unset"
	    );
	}
    }
    // Create and register copyaside buffer,
    // keyed to segment's current fill ID.
    _ca_buffers[ segment.fill_id ].append(
	(char *)segment.memory, (size_t)pos
    );
    segment.trace_logger(
	TraceLogger:: TLVL_DEBUG3,
	"Contents of created copyaside buffer: %s\n",
	_ca_buffers[ segment.fill_id ].c_str( )
    );
    // Release copyaside buffers spinlock.
    __sync_bool_compare_and_swap( &_ca_spin_lock, 1, 0 );
#ifdef WITH_INTERNAL_METRICS
    segment.pmetrics.stop_timers( );
    segment.pmetrics.accumulate_timer_deltas(
	CacheSegmentPerformanceMetrics:: MKEY_TIME_WAITING_TO_SET_SA_BUFFER
    );
    segment.pmetrics.numbytes_reserved_as_sa_buffer += pos;
#endif

}


uint64_t const
CacheManager::
get_fill_id( )
{ return _get_segment( ).fill_id; }


void
CacheManager::
_perform_segment_maintenance( CacheSegment &segment )
{

    assert( segment.avail );

    segment.trace_logger(
	TraceLogger:: TLVL_DEBUG3,
	"Performing segment maintenance....\n"
    );
    if (segment.cursor_in_ca_buffer)
	segment.trace_logger(
	    TraceLogger:: TLVL_DEBUG3,
	    "\tCursor at byte %llu in copyaside buffer for fill %llu.\n",
	    (unsigned long long int)segment.cursor,
	    (unsigned long long int)(segment.fill_id + 1)
	);
    else
	segment.trace_logger(
	    TraceLogger:: TLVL_DEBUG3,
	    "\tCursor at byte %llu in fill %llu.\n",
	    (unsigned long long int)segment.cursor,
	    (unsigned long long int)segment.fill_id
	);

    // If at end of segment and not already in copyaside buffer, 
    // then jump into copyaside buffer from next fill.
    if (!segment.cursor_in_ca_buffer && (segment.cursor == segment.size))
    {
	if (1 == _number_of_threads)
	{
	    segment.ca_buffer.clear( );
	    segment.cursor_in_ca_buffer = true;
	    segment.cursor		= 0;
	}

	else // multi-threaded
	{

	    std:: map< uint64_t, std:: string >:: iterator ca_buffers_ITER;

	    // Loop while copyaside buffer from next fill does not exist.
	    // (TODO: And there is a next fill and EOS not reached.)
	    for (uint64_t j = 0; !segment.cursor_in_ca_buffer; ++j)
	    {
		
		// If there are no more fills and EOS has been reached,
		// then we know that no copyaside buffer will be found.
		if (    (_fill_counter == (segment.fill_id + 1))
		    &&  _stream_reader.is_at_end_of_stream( ))
		{
		    segment.ca_buffer.clear( );
		    segment.cursor_in_ca_buffer = true;
		    segment.cursor		= 0;
		    break;
		}

    #ifdef WITH_INTERNAL_METRICS
		segment.pmetrics.start_timers( );
    #endif

		// Acquire copyaside buffers spinlock.
		for (   uint64_t i = 0;
			!__sync_bool_compare_and_swap( &_ca_spin_lock, 0, 1 );
			++i
		    )
		{
		    if (0 == i % 100000000)
		    {
			segment.trace_logger(
			    TraceLogger:: TLVL_DEBUG3,
			    "Waited to acquire copyaside buffers spinlock " \
			    "for %llu iterations. [read buffer]\n",
			    (unsigned long long int)i
			);
			segment.trace_logger(
			    TraceLogger:: TLVL_DEBUG4,
			    "\tSpinlock is probably %s.\n",
			    _ca_spin_lock ? "set" : "unset"
			);
		    }
		}

		// Test for existence of copyaside buffer from next fill.
		// If copyaside buffer exists, then copy it local.
		ca_buffers_ITER = _ca_buffers.find( segment.fill_id + 1 );
		if (ca_buffers_ITER != _ca_buffers.end( ))
		{
		    segment.cursor_in_ca_buffer = true;
		    segment.ca_buffer = _ca_buffers[ segment.fill_id + 1 ];
		    _ca_buffers.erase( ca_buffers_ITER );
		}

		// Release copyaside buffers spinlock.
		__sync_bool_compare_and_swap( &_ca_spin_lock, 1, 0 );

    #ifdef WITH_INTERNAL_METRICS
		segment.pmetrics.stop_timers( );
		segment.pmetrics.accumulate_timer_deltas(
		    CacheSegmentPerformanceMetrics:: MKEY_TIME_WAITING_TO_SET_SA_BUFFER
		);
    #endif

		if (0 == j % 100000000)
		    segment.trace_logger(
			TraceLogger:: TLVL_DEBUG3,
			"Waited for copyaside buffer from next fill " \
			"for %llu iterations.\n",
			(unsigned long long int)j
		    );

	    } // loop until copyaside buffer created

	    if (segment.cursor_in_ca_buffer)
		segment.cursor = 0;
		segment.trace_logger(
		    TraceLogger:: TLVL_DEBUG2, "Jumped into copyaside buffer.\n"
		);
		segment.trace_logger(
		    TraceLogger:: TLVL_DEBUG3,
		    "Contents of copyaside buffer in use: %s\n",
		    segment.ca_buffer.c_str( )
		);

	} // if multi-threaded

    } // end of segment

    // If at end of copyaside buffer...
    if (    segment.cursor_in_ca_buffer
	&&  (segment.cursor == segment.ca_buffer.length( )))
    {
	segment.cursor_in_ca_buffer	= false;
	segment.cursor			= 0;
	segment.trace_logger(
	    TraceLogger:: TLVL_DEBUG2, "Jumped out of copyaside buffer.\n"
	);
	
	_fill_segment_from_stream( segment );
    } // end of copyaside buffer

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
	    thread_id, _segment_size, _alignment, _trace_level
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
    uint64_t	nbfilled    = 0;

    // Wait while segment not selected and not end of stream.
#ifdef WITH_INTERNAL_METRICS
    segment.pmetrics.start_timers( );
#endif
wait_to_fill:
    for (   uint64_t i = 0;
	    !_stream_reader.is_at_end_of_stream( )
	&&  (_segment_to_fill != segment.thread_id);
	    ++i)
    {
	if (0 == i % 100000000)
	    segment.trace_logger(
		TraceLogger:: TLVL_DEBUG3,
		"Waited to fill segment for %llu iterations.\n",
		(unsigned long long int)i
	    );
    }
#ifdef WITH_INTERNAL_METRICS
    segment.pmetrics.stop_timers( );
    segment.pmetrics.accumulate_timer_deltas(
	CacheSegmentPerformanceMetrics::
	MKEY_TIME_WAITING_TO_FILL_FROM_STREAM
    );
#endif

    // If at end of stream, then mark segment unavailable.
    if (_stream_reader.is_at_end_of_stream( ))
    {
	segment.trace_logger(
	    TraceLogger:: TLVL_DEBUG1, "At end of input stream.\n"
	);
	segment.size	= 0;
	segment.avail	= false;
	_decrement_segment_ref_count_ATOMIC( );
    }

    // Else, refill the segment.
    else if (_check_segment_to_fill_ATOMIC( segment.thread_id ))
    {
#ifdef WITH_INTERNAL_METRICS
	segment.pmetrics.start_timers( );
#endif
	segment.size =
	    segment.cursor
	+   (	nbfilled =
		_stream_reader.read_into_cache(
		    segment.memory + segment.cursor,
		    _segment_size - segment.cursor
		));
#ifdef WITH_INTERNAL_METRICS
	segment.pmetrics.stop_timers( );
#endif
	segment.fill_id	= _fill_counter++;
	_select_segment_to_fill_ATOMIC( );
#ifdef WITH_INTERNAL_METRICS
	segment.pmetrics.numbytes_filled_from_stream += nbfilled;
	segment.pmetrics.accumulate_timer_deltas(
	    CacheSegmentPerformanceMetrics::
	    MKEY_TIME_FILLING_FROM_STREAM
	);
#endif
	segment.trace_logger(
	    TraceLogger:: TLVL_DEBUG2,
	    "Read %llu bytes into segment.\n",
	    (unsigned long long int)segment.size
	);
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


ParserPerformanceMetrics::
ParserPerformanceMetrics( )
:   numlines_copied( 0 ),
    numreads_parsed_total( 0 ),
    numreads_parsed_valid( 0 )
{ }


ParserPerformanceMetrics::
~ParserPerformanceMetrics( )
{ }


void
ParserPerformanceMetrics::
accumulate_timer_deltas( uint32_t metrics_key )
{ }


IParser * const
IParser::
get_parser(
    std:: string const	    &ifile_name,
    uint32_t const	    number_of_threads,
    uint64_t const	    cache_size,
    uint8_t const	    trace_level
)
{
    // TODO: Replace file extension detection with header magic detection.

    IStreamReader * stream_reader   = NULL;
    IParser *	    parser	    = NULL;

    std:: string    ext	    = "";
    std:: string    ifile_name_chopped( ifile_name );
    size_t	    ext_pos = ifile_name.find_last_of( "." );
    bool	    rechop  = false;

    int		    ifile_handle    = -1;
    int		    ifile_flags	    = O_RDONLY;

    if (0 < ext_pos)
    {
	ext		    = ifile_name.substr( ext_pos + 1 );
	ifile_name_chopped  = ifile_name.substr( 0, ext_pos );
    }

    if	    ("gz" == ext)
    {
	ifile_handle    = open( ifile_name.c_str( ), ifile_flags );
	if (-1 == ifile_handle) throw InvalidStreamHandle( );
#ifdef __linux__
	posix_fadvise(
	    ifile_handle, 0, 0, POSIX_FADV_SEQUENTIAL | POSIX_FADV_WILLNEED
	);
#endif
	stream_reader	= new GzStreamReader( ifile_handle );
	rechop		= true;
    } // gz
    else if ("bz2" == ext)
    {
	ifile_handle    = open( ifile_name.c_str( ), ifile_flags );
	if (-1 == ifile_handle) throw InvalidStreamHandle( );
#ifdef __linux__
	posix_fadvise(
	    ifile_handle, 0, 0, POSIX_FADV_SEQUENTIAL | POSIX_FADV_WILLNEED
	);
#endif
	stream_reader	= new Bz2StreamReader( ifile_handle );
	rechop		= true;
    } // bz2
    else // Uncompressed file.
    {
	size_t	alignment   = 0;	// 512 bytes is Chaotic Good?

#ifdef __linux__
	ifile_handle	= open( ifile_name.c_str( ), ifile_flags | O_DIRECT );
	if (-1 != ifile_handle)
	{
	    if (0 > ioctl( ifile_handle, BLKSSZGET, &alignment ))
	    {
		close( ifile_handle );
		ifile_handle = -1;
		alignment    = 0;
	    }
	}
#else
	ifile_handle    = -1;
#endif
	if (-1 == ifile_handle)
	{
	    ifile_handle    = open( ifile_name.c_str( ), ifile_flags );
	    if (-1 == ifile_handle) throw InvalidStreamHandle( );
	}
#ifdef __linux__
	if (!alignment) // Lawful Evil
	    posix_fadvise(
		ifile_handle, 0, 0, POSIX_FADV_SEQUENTIAL | POSIX_FADV_WILLNEED
	    );
#endif
	stream_reader	= new RawStreamReader( ifile_handle, alignment );
    } // uncompressed

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
    IStreamReader   &stream_reader,
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
    _number_of_threads( number_of_threads ),
    _thread_id_map( ThreadIDMap( number_of_threads ) ),
    _unithreaded( 1 == number_of_threads ),
    _states( new ParserState *[ number_of_threads ] )
{
    _uuid = rand( );

    for (uint32_t i = 0; i < number_of_threads; ++i) _states[ i ] = NULL;

    int regex_rc = -1;
    regex_rc = 
    regcomp(
	&_re_read_2_nosub,
	// ".+(/2| 2:[YN]:[[:digit:]]+:[[:alpha:]]+)$",
	"^.+(/2| 2:[YN]:[[:digit:]]+:[[:alpha:]]+).{0}",
	REG_EXTENDED | REG_NOSUB
    );
    assert( !regex_rc );
    regex_rc =
    regcomp(
	&_re_read_1,
	"^.+(/1| 1:[YN]:[[:digit:]]+:[[:alpha:]]+).{0}", REG_EXTENDED
    );
    assert( !regex_rc );
    regex_rc = 
    regcomp(
	&_re_read_2,
	"^.+(/2| 2:[YN]:[[:digit:]]+:[[:alpha:]]+).{0}", REG_EXTENDED
    );
    assert( !regex_rc );
}


IParser::
~IParser( )
{
    for (uint32_t i = 0; i < _number_of_threads; ++i)
    {
	if (_states[ i ])
	{
	    delete _states[ i ];
	    _states[ i ] = NULL;
	}
    }

    regfree( &_re_read_2_nosub );
    regfree( &_re_read_1 ); regfree( &_re_read_2 ); 
}


IParser:: ParserState::
ParserState( uint32_t const thread_id, uint8_t const trace_level )
:   at_start( true ),
    need_new_line( true ),
    buffer_pos( 0 ),
    buffer_rem( 0 ), 
    pmetrics( ParserPerformanceMetrics( ) ),
    trace_logger(
	TraceLogger(
	    trace_level, "parser-%lu.log", (unsigned long int)thread_id
	)
    )
{ memset( buffer, 0, BUFFER_SIZE + 1 ); }


IParser:: ParserState::
~ParserState( )
{ }


inline
void
IParser::
_copy_line( ParserState &state )
{
    TraceLogger	    &trace_logger   = state.trace_logger;
    uint8_t	    (&buffer)[ ParserState:: BUFFER_SIZE + 1 ]
				    = state.buffer;
    uint64_t	    &pos	    = state.buffer_pos;
    uint64_t	    &rem	    = state.buffer_rem;
    std:: string    &line	    = state.line;
    uint64_t	    i		    = 0;
    bool	    hit		    = false;

    line.clear( );

    while (true)
    {

	for (i = 0; (i < rem) && ('\n' != buffer[ pos + i ]); i++);
	if (i < rem)
	{
	    buffer[ pos + i ]   = '\0';
	    hit			= true;
	}

	trace_logger(
	    TraceLogger:: TLVL_DEBUG8,
	    "_copy_line: Detected line fragment: \"%s\"[%llu]\n",
	    (char const *)(buffer + pos), (unsigned long long int)i
	);
	line.append( (char const *)(buffer + pos), i );

	if (hit)
	{
	    rem -= (i + 1); pos += (i + 1);
	    break;
	}
	else { rem = 0; pos += i; }
	
	if (_cache_manager.has_more_data( ))
	{
	    rem = _cache_manager.get_bytes( buffer, ParserState:: BUFFER_SIZE );
	    pos = 0;
	    trace_logger(
		TraceLogger:: TLVL_DEBUG8,
		"_copy_line: Copied %llu bytes into parser buffer.\n",
		(unsigned long long int)rem
	    );
	}
	else break;

    } // while true

#ifdef WITH_INTERNAL_METRICS
    state.pmetrics.numlines_copied++;
#endif

}


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


inline
void
FastaParser::
_parse_read( ParserState &state, Read &the_read )
{
    std:: string    &line	    = state.line;

    // Validate and consume the 'name' field.
    state.trace_logger(
	TraceLogger:: TLVL_DEBUG9,
	"_parse_read: Read Name: %s\n", line.c_str( )
    );
    the_read.bytes_consumed += (line.length( ) + 1);
    if ('>' != line[ 0 ]) throw InvalidFASTAFileFormat( );    
    the_read.name = line.substr( 1 );

    // Grab sequence lines until exit conditions are met.
    while (!is_complete( ))
    {
	_copy_line( state );

	state.trace_logger(
	    TraceLogger:: TLVL_DEBUG9,
	    "_parse_read: Read Sequence (candidate): %s\n", line.c_str( )
	);

	// If a new record is detected, then existing one is complete.
	if ('>' == line[ 0 ]) break;

	// TODO? Uppercase and validate entire sequence here.
	the_read.sequence += line;
	the_read.bytes_consumed += (line.length( ) + 1);
    }

    state.trace_logger(
	TraceLogger:: TLVL_DEBUG8,
	"_parse_read: Successfully parsed FASTA record of %llu bytes.\n",
	the_read.bytes_consumed
    );
}


/* WARNING!
 * Under the following extremely rare condition, 
 * this method will incorrectly parse a FASTQ record:
 *  * A thread begins parsing in the quality scores section.
 *  * There are four or more lines of quality scores remaining.
 *  * The first line starts with '@'.
 *  * The second line starts with an alphabetic character.
 *  * The third line starts with '+' or '#'.
 *  * The fourth line is the same length as the second line.
 * Even rarer cases may occur if the third line starts with an alphabetic 
 * character and the fifth and sixth lines are the same length in total as 
 * the second and third lines, etc....
 * NOTE: This potential problem can be made even rarer by validating sequences 
 * at parse time.
 * This potential bug cannot occur if only one thread is being used to parse.
 */
inline
void
FastqParser::
_parse_read( ParserState &state, Read &the_read )
{
    std:: string    &line	    = state.line;

    // Validate and consume the 'name' field.
    state.trace_logger(
	TraceLogger:: TLVL_DEBUG9,
	"_parse_read: Read Name: %s\n", line.c_str( )
    );
    the_read.bytes_consumed += (line.length( ) + 1);
    if ('@' != line[ 0 ]) throw InvalidFASTQFileFormat( );    
    the_read.name = line.substr( 1 );

    // Grab sequence lines until exit conditions are met.
    while (!is_complete( ))
    {
	_copy_line( state );

	state.trace_logger(
	    TraceLogger:: TLVL_DEBUG9,
	    "_parse_read: Read Sequence (candidate): %s\n", line.c_str( )
	);

	// If separator line is detected, then assume sequence is complete.
	if (('+' == line[ 0 ]) || ('#' == line[ 0 ])) break;
	// TODO? Uppercase and validate entire sequence here.
	// If line starts with non-alphabetic character, 
	// then record is corrupt.
	if (	!(('A' <= line[ 0 ]) && ('Z' >= line[ 0 ]))
	    &&	!(('a' <= line[ 0 ]) && ('z' >= line[ 0 ])))
	    throw InvalidFASTQFileFormat( );

	the_read.sequence += line;
	the_read.bytes_consumed += (line.length( ) + 1);
    }

    // Ignore repeated name field.
    the_read.bytes_consumed += (line.length( ) + 1);

    // Grab quality score lines until exit conditions are met.
    while (	!is_complete( )
	    &&	(the_read.accuracy.length( ) < the_read.sequence.length( )))
    {
	_copy_line( state );

	state.trace_logger(
	    TraceLogger:: TLVL_DEBUG9,
	    "_parse_read: Read Quality Scores (candidate): %s\n",
	    line.c_str( )
	);

	the_read.accuracy += line;
	the_read.bytes_consumed += (line.length( ) + 1);
    }

    // Validate quality score lines versus sequence lines.
    if (the_read.accuracy.length( ) != the_read.sequence.length( ))
	throw InvalidFASTQFileFormat( );

    state.trace_logger(
	TraceLogger:: TLVL_DEBUG8,
	"_parse_read: Successfully parsed FASTQ record of %llu bytes.\n",
	the_read.bytes_consumed
    );

    // Prefetch next line. (Needed to be consistent with FASTA logic.)
    _copy_line( state );
}


void
IParser::
imprint_next_read( Read &the_read )
{
    
    ParserState	    &state	    = _get_state( );
    uint64_t	    &fill_id	    = state.fill_id;
    bool	    &at_start	    = state.at_start;
    bool	    &need_new_line  = state.need_new_line;
    TraceLogger	    &trace_logger   = state.trace_logger;
    uint64_t	    split_pos	    = 0;
    bool	    skip_read	    = false;
    
    while (!is_complete( ))
    {
	the_read.reset( );

	if (need_new_line) _copy_line( state );
	need_new_line = true;

	if (!at_start)
	    at_start =
		!_unithreaded
	    &&  (fill_id != _cache_manager.get_fill_id( ))
	    &&  (state.buffer_rem <= _cache_manager.whereis_cursor( ));
	if (at_start) fill_id = _cache_manager.get_fill_id( );
	trace_logger(
	    TraceLogger:: TLVL_DEBUG8,
	    "_get_next_read: fill_id = %llu, at_start = %d\n",
	    (unsigned long long int)fill_id, at_start
	);

	// Attempt to parse a read.
	// If at start of file, then error on garbage.
	// Else, skip forward to next read boundary.
	try { _parse_read( state, the_read ); }
	catch (InvalidReadFileFormat &exc)
	{
	    trace_logger(
		TraceLogger:: TLVL_DEBUG7,
		"get_next_read: Parse error on line: %s\n" \
		"\t(fill_id = %llu, at_start = %d)\n",
		state.line.c_str( ),
		(unsigned long long int)fill_id,
		at_start
	    );

	    if (!at_start || (at_start && (0 == fill_id))) throw;

	    split_pos += the_read.bytes_consumed;
	    continue;
	}

	// Skip over an unmatched second part of a paired read,
	// when at the beginning of a new fill.
	skip_read =
		at_start && (0 != fill_id)
	    &&	!regexec(
		    &_re_read_2_nosub, the_read.name.c_str( ), 0, NULL, 0
		);
	if (skip_read)
	{
	    trace_logger(
		TraceLogger:: TLVL_DEBUG7,
		"get_next_read: Skipped a split pair, using %llu bytes, " \
		"looking for next read.\n",
		(unsigned long long int)the_read.bytes_consumed
	    );
	    split_pos += the_read.bytes_consumed;
	}

	// Copy skipped over data into copyaside buffer.
	if (at_start)
	{

	    trace_logger(
		TraceLogger:: TLVL_DEBUG7,
		"get_next_read: Memory cursor is at byte %llu " \
		"in segment (fill %llu).\n",
		(unsigned long long int)_cache_manager.whereis_cursor( ),
		(unsigned long long int)_cache_manager.get_fill_id( )
	    );
	    trace_logger(
		TraceLogger:: TLVL_DEBUG7,
		"get_next_read: Parser buffer has %llu bytes remaining.\n", 
		(unsigned long long int)state.buffer_rem
	    );

	    _cache_manager.split_at( split_pos );

	    trace_logger(
		TraceLogger:: TLVL_DEBUG6,
		"get_next_read: Skipped %llu bytes of data total " \
		"at segment start.\n",
		(unsigned long long int)split_pos
	    );

	}

	split_pos	= 0;
	at_start	= false;
	need_new_line	= false;
	if (skip_read) continue;

#ifdef WITH_INTERNAL_METRICS
	state.pmetrics.numreads_parsed_total++;
#endif

	// Discard invalid read.
	if (std:: string:: npos != the_read.sequence.find_first_of( "Nn" ))
	{
	    trace_logger(
		TraceLogger:: TLVL_DEBUG6,
		"get_next_read: Discarded read \"%s\" (length %lu).\n",
		the_read.name.c_str( ),
		(unsigned long int)the_read.sequence.length( )
	    );
	    continue;
	}
	else
	    trace_logger(
		TraceLogger:: TLVL_DEBUG6,
		"get_next_read: Accepted read \"%s\" (length %lu).\n",
		the_read.name.c_str( ),
		(unsigned long int)the_read.sequence.length( )
	    );

#ifdef WITH_INTERNAL_METRICS
	state.pmetrics.numreads_parsed_valid++;
#endif
	break;
    } // while invalid read

    if (is_complete( ) && (0 == the_read.name.length( )))
	throw NoMoreReadsAvailable( );
} // imprint_next_read


void
IParser::
imprint_next_read_pair( ReadPair &the_read_pair, uint8_t mode )
{
    switch (mode)
    {
#if (0)
    case IParser:: PAIR_MODE_ALLOW_UNPAIRED:
	_imprint_next_read_pair_in_allow_mode( the_read_pair );
	break;
#endif
    case IParser:: PAIR_MODE_IGNORE_UNPAIRED:
	_imprint_next_read_pair_in_ignore_mode( the_read_pair );
	break;
    case IParser:: PAIR_MODE_ERROR_ON_UNPAIRED:
	_imprint_next_read_pair_in_error_mode( the_read_pair );
	break;
    default:
	throw UnknownPairReadingMode( );
    }
}


#if (0)
void
IParser::
_imprint_next_read_pair_in_allow_mode( ReadPair &the_read_pair )
{
    // TODO: Implement.
    //	     Probably need caching of reads between invocations 
    //	     and the ability to return pairs which are half empty. 
}
#endif


void
IParser::
_imprint_next_read_pair_in_ignore_mode( ReadPair &the_read_pair )
{
    Read	    &read_1		= the_read_pair.first;
    Read	    &read_2		= the_read_pair.second;
    regmatch_t	    match_1, match_2;

    // Hunt for a read pair until one is found or end of reads is reached.
    while (true)
    {

	// Toss out all reads which are not marked as first of a pair.
	// Note: We let any exception, which flies out of the following,
	//	 pass through unhandled.
	while (true)
	{
	    imprint_next_read( read_1 );
	    if (!regexec(
		&_re_read_1, read_1.name.c_str( ), 1, &match_1, 0
	    )) break;
	}

	// If first read of a pair was found, then insist upon second read.
	// If not found, then restart search for pair.
	// If found, then validate match.
	// If invalid pair, then restart search for pair.
	imprint_next_read( read_2 );
	if (!regexec(
	    &_re_read_2, read_2.name.c_str( ), 1, &match_2, 0
	))
	{
	    if (_is_valid_read_pair( the_read_pair, match_1, match_2 ))
		break;
	}

    } // while pair not found

} // _imprint_next_read_pair_in_ignore_mode


void
IParser::
_imprint_next_read_pair_in_error_mode( ReadPair &the_read_pair )
{
    Read	    &read_1		= the_read_pair.first;
    Read	    &read_2		= the_read_pair.second;
    regmatch_t	    match_1, match_2;

    // Note: We let any exception, which flies out of the following,
    //	     pass through unhandled.
    imprint_next_read( read_1 );
    imprint_next_read( read_2 );

    // Is the first read really the first member of a pair?
    if (REG_NOMATCH == regexec(
	&_re_read_1, read_1.name.c_str( ), 1, &match_1, 0
    )) throw InvalidReadPair( );
    // Is the second read really the second member of a pair?
    if (REG_NOMATCH == regexec(
	&_re_read_2, read_2.name.c_str( ), 1, &match_2, 0
    )) throw InvalidReadPair( );

    // Is the pair valid?
    if (!_is_valid_read_pair( the_read_pair, match_1, match_2 ))
	throw InvalidReadPair( );

} // _imprint_next_read_pair_in_error_mode


bool
IParser::
_is_valid_read_pair(
    ReadPair &the_read_pair, regmatch_t &match_1, regmatch_t &match_2
)
{
    return	(match_1.rm_so == match_2.rm_so)
	    &&	(match_1.rm_eo == match_2.rm_eo)
	    &&	(	the_read_pair.first.name.substr( 0, match_1.rm_so )
		    ==	the_read_pair.second.name.substr( 0, match_1.rm_so ));
}


} // namespace read_parsers


} // namespace khmer

// vim: set ft=cpp sts=4 sw=4 tw=80:
