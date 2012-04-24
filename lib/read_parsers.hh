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

/* Theory of Operation
 *  The greatest I/O efficiency can be achieved by pulling a large chunk of data
 *  from disk at one time and then parsing it in memory.
 *  The in-memory chunk can be split up among threads, each section being
 *  accessed via thread key supplied to 'get_next_read'. An internal index table
 *  of thread keys and section offsets will need to be maintained.
 *  This approach potentially allows I/O bandwidth saturation and also mostly 
 *  makes use of Titus' existing interfaces without introducing a parser 
 *  factory, such as Tim's approach does.
 *  Another advantage is that compressed streams do not need to be treated
 *  specially (i.e., they do not need to be worked with in a single-threaded
 *  manner due to an inability to seek on compressed streams).
 *  Note: Using mmap'd I/O may be one way to pull up different sections within a
 *  chunk for different threads. However, in some cases threads need to read
 *  into one another's windows. If they can do this with overlapping memory
 *  maps, then all is well. If not, then things may be too complicated.
 */

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

    uint8_t const	get_byte( );
    uint64_t const	get_bytes(
	uint8_t * const buffer, uint64_t buffer_len
    );

    void		split_at( uint64_t const pos );

    uint64_t const	tell( );
    void		seek( uint64_t const offset, uint8_t const whence );

    // Returns false, if current thread has more bytes to consume.
    // Blocks, if current thread has no more bytes to consume, 
    //   but other threads still do. (Synchronization barrier.)
    // Returns true, if no threads have more bytes to consume.
    bool const		has_more_data( );
    
private:
    
    struct SegmentInfo
    {

	uint64_t	    segment_start;
	uint64_t	    segment_end;
	uint64_t	    cursor;
	
	SegmentInfo(
	    uint64_t const segment_start, uint64_t const segment_end
	);
	~SegmentInfo( );

	bool const	    get_availability( ) const;
	bool const	    get_availability_ATOMIC( );
	void		    set_availability_ATOMIC( bool const available );

    private:
	
	bool		    _available;

    }; // struct SegmentInfo

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
    uint64_t		_cache_size;
    uint8_t *		_cache;
    SegmentInfo **	_psegment_infos;
    SegmentInfo **	_ssegment_infos;

    // Offset from start of cache.
    // Coincides with next segment to refill.
    uint64_t		_fill_cursor;

    // Extends or refills segment for current thread, as needed.
    void		_perform_segment_maintenance(
	SegmentInfo & psegment_info
    );

    // Atmoically checks the fill cursor against the given position.
    bool const		_check_fill_cursor_ATOMIC( uint64_t const pos );
    // Atomically sets the fill cursor to the given position.
    void		_set_fill_cursor_ATOMIC( uint64_t const pos );

    // Retrieves psegment for current thread.
    SegmentInfo &	_get_psegment( bool const higher = false );
    // Retrieves ssegment for current thread or the one with the next lower ID.
    SegmentInfo &	_get_ssegment( bool const lower = false );
    // Used to identify which segments to use.
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
