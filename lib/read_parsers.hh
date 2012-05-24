#ifndef READ_PARSERS_HH
#define READ_PARSERS_HH


#include <cassert>
#include <cstdarg>

#include <string>
#include <map>

#ifdef __linux__
#   include <sys/types.h>
#else
#   error "Your current operating system is not supported by this software."
#endif

#include "zlib/zlib.h"
#include "bzip2/bzlib.h"

#include "khmer_config.hh"


#define MIN( a, b )	(((a) > (b)) ? (b) : (a))


namespace khmer
{


struct InvalidStreamHandle : public std:: exception
{ };

struct InvalidStreamBuffer : public std:: exception
{ };

struct StreamReadError : public std:: exception
{ };

struct InvalidPerformanceMetricsKey : public std:: exception
{ };

struct InvalidNumberOfThreadsRequested : public std:: exception
{ };

struct TooManyThreads : public std:: exception
{ };


struct TraceLogger
{
    
    enum
    {
	TLVL_ALL	= 0,
	TLVL_DEBUG9, TLVL_DEBUG8, TLVL_DEBUG7, TLVL_DEBUG6, TLVL_DEBUG5,
	TLVL_DEBUG4, TLVL_DEBUG3, TLVL_DEBUG2, TLVL_DEBUG1, TLVL_DEBUG0,
	TLVL_INFO9, TLVL_INFO8, TLVL_INFO7, TLVL_INFO6, TLVL_INFO5,
	TLVL_INFO4, TLVL_INFO3, TLVL_INFO2, TLVL_INFO1, TLVL_INFO0,
	TLVL_WARNING	= 30,
	TLVL_ERROR	= 40,
	TLVL_CRITICAL	= 50,
	TLVL_NONE
    };
    
    TraceLogger( uint8_t const level, FILE * stream_handle );
    TraceLogger(
	uint8_t const level, char const * const file_name_format, ...
    );
    ~TraceLogger( );

    inline void	    operator( )(
	uint8_t const level, char const * const format, ...
    ) const
    {
	va_list varargs;
	
	if (_level <= level)
	{
	    va_start( varargs, format );
	    vfprintf( _stream_handle, format, varargs );
	    va_end( varargs );
	    fflush( _stream_handle );
	}

    }

private:
    
    uint8_t	    _level;
    bool	    _shared_stream;
    FILE *	    _stream_handle;

};


struct IPerformanceMetrics
{

	    IPerformanceMetrics( );
    virtual ~IPerformanceMetrics( );
    
    inline void	    start_timers( )
    {
	clock_gettime( CLOCK_REALTIME, &_temp_clock_start );
	clock_gettime( CLOCK_THREAD_CPUTIME_ID, &_temp_cpu_start );
    }
    inline void	    stop_timers( )
    {
	clock_gettime( CLOCK_THREAD_CPUTIME_ID, &_temp_cpu_stop );
	clock_gettime( CLOCK_REALTIME, &_temp_clock_stop );
    }
    virtual void    accumulate_timer_deltas( uint32_t metrics_key )	= 0;
    
    // TODO: Add a printing or log file feature.

protected:
    
    timespec	_temp_cpu_start;
    timespec	_temp_cpu_stop;
    timespec	_temp_clock_start;
    timespec	_temp_clock_stop;

    uint64_t const  _timespec_diff_in_nsecs(
	timespec const &start, timespec const &stop
    );

};


struct ThreadIDMap
{

    ThreadIDMap( uint32_t number_of_threads );
    ~ThreadIDMap( );

    uint32_t const get_thread_id( );

private:

    uint32_t			    _number_of_threads;
    uint32_t			    _thread_counter;
#ifdef __linux__
    std:: map< pid_t, uint32_t >    _thread_id_map;
#else
    // TODO: Maybe try something with pthreads for the general case.
#endif
    uint32_t			    _tid_map_spin_lock;

};


namespace read_parsers
{


struct InvalidFASTAFileFormat: public std:: exception
{ };

struct InvalidFASTQFileFormat: public std:: exception
{ };

struct CacheSegmentUnavailable : public std:: exception
{ };

struct CacheSegmentBoundaryViolation : public std:: exception
{ };

struct InvalidCacheSizeRequested : public std:: exception
{ };


struct StreamReaderPerformanceMetrics : public IPerformanceMetrics
{
    
    enum
    {
	MKEY_TIME_READING
    };
    
    uint64_t	    numbytes_read;
    uint64_t	    clock_nsecs_reading;
    uint64_t	    cpu_nsecs_reading;
    
	    StreamReaderPerformanceMetrics( );
    virtual ~StreamReaderPerformanceMetrics( );

    virtual void    accumulate_timer_deltas( uint32_t metrics_key );

};


struct IStreamReader
{

    StreamReaderPerformanceMetrics  pmetrics;
    
	    IStreamReader( );
    virtual ~IStreamReader( );

    bool const			    is_at_end_of_stream( ) const;

    virtual uint64_t const	    read_into_cache(
	uint8_t * const cache, uint64_t const cache_size
    ) = 0;

protected:
    
    bool			    _at_eos;

};


struct RawStreamReader : public IStreamReader
{
    
	    RawStreamReader( int const fd );
    virtual ~RawStreamReader( );

    virtual uint64_t const  read_into_cache(
	uint8_t * const cache, uint64_t const cache_size
    );

protected:
    
    int			    _stream_handle;

};


struct GzStreamReader : public IStreamReader
{
    
	    GzStreamReader( int const fd );
    virtual ~GzStreamReader( );

    virtual uint64_t const  read_into_cache(
	uint8_t * const cache, uint64_t const cache_size
    );

private:
    
    gzFile		    _stream_handle;

};


struct Bz2StreamReader : public IStreamReader
{
    
	    Bz2StreamReader( int const fd );
    virtual ~Bz2StreamReader( );

    virtual uint64_t const  read_into_cache(
	uint8_t * const cache, uint64_t const cache_size
    );

private:
    
    FILE *		    _stream_handle;
    BZFILE *		    _block_handle;

};


struct CacheSegmentPerformanceMetrics : public IPerformanceMetrics
{
    
    enum
    {
	MKEY_TIME_WAITING_TO_SET_SA_BUFFER,
	MKEY_TIME_WAITING_TO_GET_SA_BUFFER,
	MKEY_TIME_WAITING_TO_FILL_FROM_STREAM,
	MKEY_TIME_FILLING_FROM_STREAM,
	MKEY_TIME_IN_SYNC_BARRIER
    };
    
    uint64_t	    numbytes_filled_from_stream;
    uint64_t	    numbytes_copied_from_sa_buffer;
    uint64_t	    numbytes_reserved_as_sa_buffer;
    uint64_t	    numbytes_copied_to_caller_buffer;
    uint64_t	    clock_nsecs_waiting_to_set_sa_buffer;
    uint64_t	    cpu_nsecs_waiting_to_set_sa_buffer;
    uint64_t	    clock_nsecs_waiting_to_get_sa_buffer;
    uint64_t	    cpu_nsecs_waiting_to_get_sa_buffer;
    uint64_t	    clock_nsecs_waiting_to_fill_from_stream;
    uint64_t	    cpu_nsecs_waiting_to_fill_from_stream;
    uint64_t	    clock_nsecs_filling_from_stream;
    uint64_t	    cpu_nsecs_filling_from_stream;
    uint64_t	    clock_nsecs_in_sync_barrier;
    uint64_t	    cpu_nsecs_in_sync_barrier;

	    CacheSegmentPerformanceMetrics( );
    virtual ~CacheSegmentPerformanceMetrics( );

    virtual void    accumulate_timer_deltas( uint32_t metrics_key );

    virtual void    accumulate_metrics(
	CacheSegmentPerformanceMetrics &source
    );

protected:
    
    uint32_t	    _accumulated_count;

};


struct CacheManager
{
    
    CacheManager(
	IStreamReader &	stream_reader,
	uint32_t const	number_of_threads,
	uint64_t const	cache_size,
	uint8_t const	trace_level = 0
    );
    ~CacheManager( );

    // Returns true, if current thread has more bytes to consume.
    // Blocks, if current thread has no more bytes to consume, 
    //   but other threads still do. (Synchronization barrier.)
    // Returns false, if no threads have more bytes to consume.
    bool const		has_more_data( );

    uint64_t const	get_bytes(
	uint8_t * const buffer, uint64_t buffer_len
    );

    uint64_t const	whereis_cursor( );
    void		split_at( uint64_t const pos );

    uint64_t const	get_fill_id( );

    // NOTE: The following methods should not be needed in "real world"
    //	     situations; they exist to help testing and tracing.
    bool const		_in_sa_buffer( );
    bool const		_sa_buffer_avail( );
    
private:
    
    struct CacheSegment
    {

	bool				avail;
	uint32_t			thread_id;
	uint64_t			size;
	uint8_t *			memory;
	uint64_t			cursor;
	bool				cursor_in_sa_buffer;
	uint64_t			sa_buffer_size;
	uint64_t			fill_id;
	CacheSegmentPerformanceMetrics	pmetrics;
	TraceLogger			trace_logger;
	
	CacheSegment(
	    uint32_t const  thread_id,
	    uint64_t const  size,
	    uint8_t const   trace_level = 0
	);
	~CacheSegment( );

	bool				get_sa_buffer_avail( ) const;
	bool				get_sa_buffer_avail_ATOMIC( );
	void				set_sa_buffer_avail_ATOMIC( bool const avail );

    private:
	
	bool				_sa_buffer_avail;

    }; // struct CacheSegment

    uint8_t		_trace_level;

    IStreamReader &	_stream_reader;

    uint32_t		_number_of_threads;
    ThreadIDMap		_thread_id_map;

    uint64_t		_segment_size;
    CacheSegment **	_segments;
    uint32_t		_segment_ref_count;
    uint32_t		_segment_to_fill;
    uint64_t		_fill_counter;

    // Extends or refills segment for current thread, as needed.
    void		_perform_segment_maintenance(
	CacheSegment & segment
    );

    bool const		_check_segment_to_fill_ATOMIC(
	uint32_t const thread_id
    );
    void		_select_segment_to_fill_ATOMIC( );
    CacheSegment &	_get_segment( bool const higher = false );
    void		_fill_segment_from_stream(
	CacheSegment & segment
    );
    void		_increment_segment_ref_count_ATOMIC( );
    void		_decrement_segment_ref_count_ATOMIC( );
    uint32_t const	_get_segment_ref_count_ATOMIC( );
    
}; // struct CacheManager


struct Read
{
    std:: string name;
    std:: string seq;
};


struct ParserPerformanceMetrics: public IPerformanceMetrics
{
    
    uint64_t	    numlines_copied;
    uint64_t	    numreads_parsed_total;
    uint64_t	    numreads_parsed_valid;

	    ParserPerformanceMetrics( );
    virtual ~ParserPerformanceMetrics( );

    virtual void    accumulate_timer_deltas( uint32_t metrics_key );

};


struct IParser
{

    static IParser * const  get_parser(
	std:: string const 	&ifile_name,
	uint32_t const		number_of_threads   =
	khmer:: get_active_config( ).get_number_of_threads( ),
	uint64_t const		cache_size	    =
	khmer:: get_active_config( ).get_reads_file_chunk_size( ),
	uint8_t const		trace_level	    = 0
    );
    
	    IParser(
	IStreamReader	&stream_reader,
	uint32_t const	number_of_threads   =
	khmer:: get_active_config( ).get_number_of_threads( ),
	uint64_t const	cache_size	    =
	khmer:: get_active_config( ).get_reads_file_chunk_size( ),
	uint8_t const	trace_level	    = 0
    );
    virtual ~IParser( );

    inline bool		is_complete( )
    { return !_cache_manager.has_more_data( ) && !_get_state( ).buffer_rem; }

    virtual Read	get_next_read( )    = 0;

protected:
    
    struct ParserState
    {

	// TODO: Set buffer size from Config.
	static uint64_t const	    BUFFER_SIZE		= 127;

	uint32_t		    thread_id;
	
	bool			    at_start;
	uint64_t		    fill_id;

	std:: string		    line;
	bool			    need_new_line;

	uint8_t			    buffer[ BUFFER_SIZE + 1 ];
	uint64_t		    buffer_pos;
	uint64_t		    buffer_rem;

	ParserPerformanceMetrics    pmetrics;
	TraceLogger		    trace_logger;
	
	ParserState( uint32_t const thread_id, uint8_t const trace_level );
	~ParserState( );

    }; // struct ParserState
    
    uint8_t		_trace_level;

    CacheManager	_cache_manager;

    ThreadIDMap		_thread_id_map;
    bool		_unithreaded;

    ParserState **	_states;

    void		_copy_line( ParserState &state );

    inline ParserState	&_get_state( )
    {
	uint32_t	thread_id	= _thread_id_map.get_thread_id( );
	ParserState *	state_PTR	= NULL;

	assert( NULL != _states );

	state_PTR = _states[ thread_id ];
	if (NULL == state_PTR)
	{
	    _states[ thread_id ]    =
	    new ParserState( thread_id, _trace_level );
	    state_PTR		    = _states[ thread_id ];
	}

	return *state_PTR;
    }

}; // struct IParser


struct FastaParser : public IParser
{
    
	    FastaParser(
	IStreamReader &	stream_reader,
	uint32_t const	number_of_threads,
	uint64_t const	cache_size,
	uint8_t const	trace_level
    );
    virtual ~FastaParser( );

    virtual Read    get_next_read( );

};


struct FastqParser : public IParser
{

	    FastqParser(
	IStreamReader &	stream_reader,
	uint32_t const	number_of_threads,
	uint64_t const	cache_size,
	uint8_t const	trace_level
    );
    virtual ~FastqParser( );

    virtual Read    get_next_read( );

};


} // namespace read_parsers


} // namespace khmer


#endif // READ_PARSERS_HH

// vim: set ft=cpp sts=4 sw=4 tw=80:
