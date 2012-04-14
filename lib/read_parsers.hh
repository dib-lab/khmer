#ifndef READ_PARSERS_HH
#define READ_PARSERS_HH


#include <string>

#include "zlib/zlib.h"
#include "bzip2/bzlib.h"


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

    bool const is_at_end_of_stream( ) const;

    virtual uint64_t const
    read_into_cache( uint8_t * const cache, uint64_t const cache_size ) = 0;

protected:
    
    bool at_eos;

};


struct RawStreamReader : public IStreamReader
{
    
	    RawStreamReader( int const fd );
    virtual ~RawStreamReader( );

    virtual uint64_t const
    read_into_cache( uint8_t * const cache, uint64_t const cache_size );

protected:
    
    int			_stream_handle;

};


struct GzStreamReader : public IStreamReader
{
    
	    GzStreamReader( int const fd );
    virtual ~GzStreamReader( );

    virtual uint64_t const
    read_into_cache( uint8_t * const cache, uint64_t const cache_size );

private:
    
    gzFile	    _stream_handle;

};


struct Bz2StreamReader : public IStreamReader
{
    
	    Bz2StreamReader( int const fd );
    virtual ~Bz2StreamReader( );

    virtual uint64_t const
    read_into_cache( uint8_t * const cache, uint64_t const cache_size );

private:
    
    FILE *	    _stream_handle;
    BZFILE *	    _block_handle;

};


struct IParser
{


    struct ReadsCacheManager
    {
	
	uint8_t *   extend_segment( uint32_t );
	uint32_t    get_next_segment( );

    private:
	
	struct SegmentInfo
	{
	    
	    SegmentInfo( uint8_t *, uint64_t );
	    ~SegmentInfo( );

	private:
	    
	    uint8_t *	_start_address;
	    uint8_t *	_end_address;
	    // TODO? Add lock.
	    // TODO? Add stream pointer.
	    // TODO? Add thread ID.

	};
	
	//IStreamReader &	_stream_reader;
	//uint32_t	_segment_counter;

	bool		_read_into_segment( uint32_t );

    };
    
	    IParser( uint32_t const number_of_threads = 1 );
    virtual ~IParser( );

    virtual bool    is_complete( )		= 0;
    virtual Read    get_next_read( uint32_t )	= 0;

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
