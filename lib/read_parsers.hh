#ifndef READ_PARSERS_HH
#define READ_PARSERS_HH


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


#define MIN( a, b )	((a > b) ? (b) : (a))


namespace khmer
{


struct InvalidNumberOfThreadsRequested : public std:: exception
{ };

struct TooManyThreads : public std:: exception
{ };


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


struct InvalidStreamBuffer : public std:: exception
{ };

struct StreamReadError : public std:: exception
{ };

struct CacheSegmentUnavailable : public std:: exception
{ };

struct CacheSegmentBoundaryViolation : public std:: exception
{ };

struct InvalidCacheSizeRequested : public std:: exception
{ };


struct Read
{
    std:: string name;
    std:: string seq;
};


struct IStreamReader
{
    
	    IStreamReader( );
    virtual ~IStreamReader( );

    bool const		    is_at_end_of_stream( ) const;

    virtual uint64_t const  read_into_cache(
	uint8_t * const cache, uint64_t const cache_size
    ) = 0;

protected:
    
    bool		    at_eos;

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

    uint8_t const	get_byte( );
    uint64_t const	get_bytes(
	uint8_t * const buffer, uint64_t buffer_len
    );

    uint64_t const	whereis_cursor( );
    void		split_at( uint64_t const pos );

    uint64_t const	get_fill_id( );

    // NOTE: The following methods should not be needed in "real world"
    //	     sitatuions. They exist to help the test harness perform some more
    //	     intelligent testing.
    bool const		_in_sa_buffer( );
    bool const		_sa_buffer_avail( );
    
private:
    
    struct CacheSegment
    {

	bool		avail;
	uint32_t	thread_id;
	uint64_t	size;
	uint8_t *	memory;
	uint64_t	cursor;
	bool		cursor_in_sa_buffer;
	uint64_t	sa_buffer_size;
	uint64_t	fill_id;
	
	CacheSegment(
	    uint32_t const  thread_id,
	    uint64_t const  size,
	    uint8_t const   trace_level = 0
	);
	~CacheSegment( );

	bool		get_sa_buffer_avail( ) const;
	bool		get_sa_buffer_avail_ATOMIC( );
	void		set_sa_buffer_avail_ATOMIC( bool const avail );

	FILE *		trace_file_handle;

	uint64_t	_nsecs_waiting_to_set_sa_buffer;
	uint64_t	_nsecs_waiting_to_acquire_sa_buffer;
	uint64_t	_nsecs_waiting_to_fill_from_stream;
	uint64_t	_nsecs_reading_from_stream;
	uint64_t	_nsecs_in_final_sync_barrier;

    private:
	
	uint8_t		_trace_level;
	bool		_sa_buffer_avail;

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


struct IParser
{
    
    static IParser * const  get_parser(
	std:: string const &	ifile_name,
	uint32_t const		number_of_threads   =
	khmer:: get_active_config( ).get_number_of_threads( ),
	uint64_t const		cache_size	    =
	khmer:: get_active_config( ).get_reads_file_chunk_size( ),
	uint8_t const		trace_level	    = 0
    );
    
	    IParser(
	IStreamReader &	stream_reader,
	uint32_t const	number_of_threads   =
	khmer:: get_active_config( ).get_number_of_threads( ),
	uint64_t const	cache_size	    =
	khmer:: get_active_config( ).get_reads_file_chunk_size( ),
	uint8_t const	trace_level	    = 0
    );
    virtual ~IParser( );

	    bool	is_complete( );
    virtual Read	get_next_read( )	    = 0;

protected:
    
    struct ParserState
    {

	// TODO: Set buffer size from Config.
	static uint64_t const	BUFFER_SIZE	    = 127;
	
	bool		at_start;
	uint64_t	fill_id;

	std:: string	line;
	bool		need_new_line;

	uint8_t		buffer[ BUFFER_SIZE + 1 ];
	uint64_t	buffer_pos;
	uint64_t	buffer_rem;
	
	ParserState( );
	~ParserState( );

    };
    
    uint8_t		_trace_level;

    CacheManager	_cache_manager;

    ThreadIDMap		_thread_id_map;

    ParserState **	_states;

    void		_copy_line( );

    ParserState &	_get_state( );

};


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
