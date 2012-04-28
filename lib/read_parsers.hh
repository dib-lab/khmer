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


#define MIN( a, b )	((a > b) ? (b) : (a))


namespace khmer
{


namespace read_parsers
{


struct InvalidStreamBuffer : public std:: exception
{ };

struct StreamReadError : public std:: exception
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
	IStreamReader *	stream_reader,
	uint32_t const	number_of_threads,
	uint64_t const	cache_size
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

    void		split_at( uint64_t const pos );

    uint64_t const	tell( );
    void		seek( uint64_t const offset, uint8_t const whence );
    
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
	
	CacheSegment( uint32_t const thread_id, uint64_t const size );
	~CacheSegment( );

	bool		get_sa_buffer_avail( ) const;
	bool		get_sa_buffer_avail_ATOMIC( );
	void		set_sa_buffer_avail_ATOMIC( bool const avail );

	uint64_t	_nsecs_waiting_to_set_sa_buffer;
	uint64_t	_nsecs_waiting_to_acquire_sa_buffer;
	uint64_t	_nsecs_waiting_to_fill_from_stream;
	uint64_t	_nsecs_reading_from_stream;
	uint64_t	_nsecs_in_final_sync_barrier;

    private:
	
	bool		_sa_buffer_avail;

    }; // struct CacheSegment

    IStreamReader *	_stream_reader_PTR;

    uint32_t		_number_of_threads;
    uint32_t		_thread_counter;
#ifdef __linux__
    std:: map< pid_t, uint32_t >    _thread_id_map;
#else
    // TODO: Maybe try something with pthreads for the general case.
#endif
    uint32_t		_tid_map_spin_lock;

    uint64_t		_segment_size;
    CacheSegment **	_segments;
    uint32_t		_segment_ref_count;
    uint32_t		_segment_to_fill;

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
    uint32_t const	_get_thread_id( );
    
}; // struct CacheManager


struct IParser
{
    
    // TODO: Get defaults from Config interface.
	    IParser(
		uint32_t const number_of_threads,
		uint64_t const cache_size
	    );
    virtual ~IParser( );

    virtual bool	is_complete( )		    = 0;
    virtual Read	get_next_read( uint32_t )   = 0;

protected:
    
    CacheManager *	_cache_manager;

};


struct FastaParser : public IParser
{

};


struct FastqParser : public IParser
{

};


} // namespace read_parsers


} // namespace khmer


#endif // READ_PARSERS_HH

// vim: set ft=cpp sts=4 sw=4 tw=80:
