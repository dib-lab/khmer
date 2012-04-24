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

#ifdef __linux__
#   include <sys/types.h>
#   include <sys/syscall.h>
#endif

#include "read_parsers.hh"


namespace khmer
{


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
    IStreamReader * stream_reader_PTR,
    uint32_t const  number_of_threads,
    uint64_t const  cache_size
)
{
    
    _stream_reader_PTR	= stream_reader_PTR;

    // TODO: Verify that the number of threads > 0.
    _number_of_threads	= number_of_threads;
    // Renormalize cache size by way of integer division.
    _segment_size	= cache_size / number_of_threads;
    _cache_size		= _segment_size * number_of_threads;

    // Allocate cache and segment information.
    _cache		= new uint8_t[ _cache_size ];
    _psegment_infos	= new SegmentInfo *[ number_of_threads ];
    _ssegment_infos	= new SegmentInfo *[ number_of_threads ];

    _thread_counter	= 0;
    _tid_map_spin_lock	= 0;

    _fill_cursor	= 0;

}


CacheManager::
~CacheManager( )
{
    
    for (uint64_t i = 0; i < _number_of_threads ; ++i)
    {
	delete _psegment_infos[ i ];
	delete _ssegment_infos[ i ];
    }
    delete [ ] _psegment_infos;
    _psegment_infos	= NULL;
    delete [ ] _ssegment_infos;
    _ssegment_infos	= NULL;
    delete [ ] _cache;
    _cache		= NULL;

}


CacheManager:: SegmentInfo::
SegmentInfo( uint64_t const start, uint64_t const end )
{
    cursor	    = end;
    segment_start   = start;
    segment_end	    = end;
    _available	    = false;
}


CacheManager:: SegmentInfo::
~SegmentInfo( )
{
    // TODO: Implement.
}


inline
bool const
CacheManager:: SegmentInfo::
get_availability( ) const
{
    return _available;
}


inline
bool const
CacheManager:: SegmentInfo::
get_availability_ATOMIC( )
{
    return __sync_and_and_fetch( &_available, 1 );
}


inline
void
CacheManager:: SegmentInfo::
set_availability_ATOMIC( bool const available )
{
    __sync_bool_compare_and_swap( &_available, !available, available );
}


uint8_t const
CacheManager::
get_byte( )
{
    SegmentInfo &	segment_info	    = _get_psegment( );   

    _perform_segment_maintenance( segment_info );

    return _cache[ segment_info.cursor++ ];
}


inline
uint64_t const
CacheManager::
get_bytes( uint8_t * const buffer, uint64_t buffer_len )
{
    SegmentInfo &	segment_info	    = _get_psegment( );   
    uint64_t		nbcopied	    = 0;
    uint64_t		nbcopied_total	    = 0;

    for (uint64_t nbrem = buffer_len; (nbrem > 0); nbrem -= nbcopied)
    {
	_perform_segment_maintenance( segment_info );

	nbcopied = MIN( nbrem, segment_info.segment_end - segment_info.cursor );
	memcpy( buffer, _cache + segment_info.cursor, nbcopied );
	segment_info.cursor += nbcopied;

	nbcopied_total += nbcopied;
    }

    return nbcopied_total;
}


void
CacheManager::
split_at( uint64_t const pos )
{
    // TODO: Implement.
    //	     Place start address of current segment into previous segment.
    //	     Set size of previous segment.to pos.
    //	     Set start address of current segment to pos.
    //	     Set size of current segment to size - pos.
    //	     Set offset of current segment to offset - pos.
}


inline
uint64_t const
CacheManager::
tell( )
{
    return _get_psegment( ).cursor;
}


inline
void
CacheManager::
seek( uint64_t const offset, uint8_t const whence )
{
    SegmentInfo &	segment_info	    = _get_psegment( );
    uint64_t		cursor;

    switch (whence)
    {

    case SEEK_SET:
	cursor = segment_info.segment_start + offset;
	break;

    case SEEK_CUR:
	cursor = segment_info.cursor + offset;
	break;

    default:
	// TODO: UnknownSegmentSeekTypeError.
	// TEMP: Do nothing.
	return;

    }

    if (    (cursor < segment_info.segment_start)
	||  (cursor > segment_info.segment_end))
	// TODO: Throw SegmentBoundaryError.
	;
    segment_info.cursor = cursor;
}

void
CacheManager::
_perform_segment_maintenance( SegmentInfo & psegment_info )
{
    SegmentInfo &	    ssegment_info   = _get_ssegment( false );

    // Attempt to merge ssegment into psegment.
    if (psegment_info.cursor == psegment_info.segment_end)
    {
	SegmentInfo &	    npsegment_info  = _get_psegment( true );
	
	// Wait while next psegment is not in final state,
	// and our ssegment has not been made available.
	while (	(npsegment_info.segment_start != npsegment_info.segment_end)
	    &&	!ssegment_info.get_availability( ));
	// Atomically test whether our ssegment is available.
	// If the test succeeds, then proceed with merge.
	if (	ssegment_info.get_availability_ATOMIC( )
	    &&	(ssegment_info.segment_start != ssegment_info.segment_end))
	{
	    psegment_info.segment_end = ssegment_info.segment_end;
	}
    } // Segment Merge

    // If at EOS, then mark psegment off-line.
    if (    (psegment_info.cursor == psegment_info.segment_end)
	&&  _stream_reader_PTR->is_at_end_of_stream( ))
    {
	psegment_info.set_availability_ATOMIC( false );
	psegment_info.segment_end   = 0;
	psegment_info.segment_start = 0;
	psegment_info.cursor	    = 0;
	// TODO? ssegment bookkeeping.
	return;
    }

    // Attempt to refill psegment.
    if (psegment_info.cursor == psegment_info.segment_end)
    {
	// Wait while fill cursor is not at psegment start.
	while (	(_fill_cursor != psegment_info.segment_start)
	    &&	!_stream_reader_PTR->is_at_end_of_stream( ));
	// Atomically test whether fill cursor is at our psegment.
	// If test succeeds, then fill the segment.
	if (	_check_fill_cursor_ATOMIC( psegment_info.segment_start )
	    &&	!_stream_reader_PTR->is_at_end_of_stream( ))
	{
	    ssegment_info.set_availability_ATOMIC( false );
	    psegment_info.cursor	= psegment_info.segment_start;
	    psegment_info.segment_end	= 
		psegment_info.segment_start
	    +	_stream_reader_PTR->read_into_cache(
		    _cache + psegment_info.segment_start, 
		    psegment_info.segment_end - psegment_info.segment_start
		);
	    _set_fill_cursor_ATOMIC( psegment_info.segment_end );
	}
    }
}


inline
bool const
CacheManager::
_check_fill_cursor_ATOMIC( uint64_t const pos )
{
    return (pos == __sync_val_compare_and_swap( &_fill_cursor, pos, pos ));
}


inline
void
CacheManager::
_set_fill_cursor_ATOMIC( uint64_t const pos )
{
    __sync_bool_compare_and_swap( &_fill_cursor, _fill_cursor, pos );
}


inline
CacheManager:: SegmentInfo &
CacheManager::
_get_psegment( bool const higher )
{
    uint32_t	    thread_id		= _get_thread_id( );
    SegmentInfo *   segment_info_PTR	= NULL;

    assert( NULL != _psegment_infos );

    // If referring to a segment to snoop,
    // then index is for the thread with the next higher ID.
    if (higher) thread_id = ((thread_id + 1) % _number_of_threads);

    segment_info_PTR = _psegment_infos[ thread_id ];
    if (NULL == segment_info_PTR)
    {
	_psegment_infos[ thread_id ]	= 
	new SegmentInfo(
	    thread_id * _segment_size, thread_id * (_segment_size + 1)
	);
	segment_info_PTR		= _psegment_infos[ thread_id ];
    }

    return *segment_info_PTR;
}


inline
CacheManager:: SegmentInfo &
CacheManager::
_get_ssegment( bool const lower )
{
    uint32_t	    thread_id		= _get_thread_id( );
    SegmentInfo *   segment_info_PTR	= NULL;

    assert( NULL != _ssegment_infos );

    // If referring to a split-off segment, 
    // then index is for the thread with the next lower ID.
    if (lower)
    {
	if (0 == thread_id) thread_id = _number_of_threads - 1;
	else		    thread_id--;
    }

    segment_info_PTR = _ssegment_infos[ thread_id ];
    if (NULL == segment_info_PTR)
    {
	_ssegment_infos[ thread_id ]	= new SegmentInfo( 0, 0 );
	segment_info_PTR		= _psegment_infos[ thread_id ];
    }

    return *segment_info_PTR;
}


uint32_t const
CacheManager::
_get_thread_id( )
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
	uint32_t	    thread_id;
	if (_number_of_threads <= _thread_id_map.size( ))
	    // TODO: Implement. Throw exception.
	    ;
	thread_id = __sync_fetch_and_add( &_thread_counter, 1 );
	while ( !__sync_bool_compare_and_swap( &_tid_map_spin_lock, 0, 1 ) );
	_thread_id_map[ native_thread_id ] = thread_id;
	__sync_bool_compare_and_swap( &_tid_map_spin_lock, 1, 0 );
	return thread_id;
    }

    return (*match).second; // Return the found value.
}


IParser::
IParser( uint32_t const number_of_threads, uint64_t const cache_size )
{
    // TODO: Implement.
}


IParser::
~IParser( )
{
    // TODO: Implement.
}


} // namespace read_parsers


} // namespace khmer

// vim: set ft=cpp sts=4 sw=4 tw=80:
