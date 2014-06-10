//
// This file is part of khmer, http://github.com/ged-lab/khmer/, and is
// Copyright (C) Michigan State University, 2009-2013. It is licensed under
// the three-clause BSD license; see doc/LICENSE.txt.
// Contact: khmer-project@idyll.org
//

#ifndef READ_PARSERS_HH
#define READ_PARSERS_HH


#include <cassert>
#include <cstdarg>
#include <iostream>
#include <string>
#include <utility>
#include <stdlib.h>

extern "C"
{
#include <regex.h>
}

#include "zlib.h"
#include "bzlib.h"

#include "khmer.hh"
#include "khmer_config.hh"
#include "thread_id_map.hh"
#include "trace_logger.hh"
#include "perf_metrics.hh"

namespace khmer
{



namespace read_parsers
{


class InvalidReadFileFormat: public khmer:: khmer_file_exception  {

public:

    InvalidReadFileFormat(
        char const * exc_name,
        char const * reason	= NULL,
        char const * evidence	= NULL
    );

    virtual char const *    what( ) const throw( );

protected:

    char		    _reason[ CHAR_MAX ];

};

class InvalidFASTAFileFormat: public InvalidReadFileFormat {

public:
    InvalidFASTAFileFormat(
        char const * reason	= NULL,
        char const * evidence	= NULL
    );

};

class InvalidFASTQFileFormat: public InvalidReadFileFormat {

public:
    InvalidFASTQFileFormat(
        char const * reason	= NULL,
        char const * evidence	= NULL
    );

};

struct CacheSegmentUnavailable : public std:: exception {
};

struct CacheSegmentBoundaryViolation : public std:: exception {
};

struct InvalidCacheSizeRequested : public std:: exception {
};

struct NoMoreReadsAvailable : public std:: exception {
};

struct UnknownPairReadingMode : public std:: exception {
};

struct InvalidReadPair : public std:: exception {
};

#ifdef WITH_INTERNAL_METRICS
struct StreamReaderPerformanceMetrics : public IPerformanceMetrics {

    enum {
        MKEY_TIME_READING
    };

    uint64_t	    numbytes_read;
    uint64_t	    clock_nsecs_reading;
    uint64_t	    cpu_nsecs_reading;

    StreamReaderPerformanceMetrics( );
    virtual ~StreamReaderPerformanceMetrics( );

    virtual void    accumulate_timer_deltas( uint32_t metrics_key );

};
#endif

struct IStreamReader {
#ifdef WITH_INTERNAL_METRICS
    StreamReaderPerformanceMetrics  pmetrics;
#endif
    IStreamReader( );
    virtual ~IStreamReader( );

    size_t const		    get_memory_alignment( ) const;

    bool const			    is_at_EOS_ATOMIC( );

    virtual uint64_t const	    read_into_cache(
        uint8_t * const cache, uint64_t const cache_size
    ) = 0;

protected:

    size_t			    _alignment;
    size_t			    _max_aligned;

    bool			    _at_eos;

    void			    _set_EOS_ATOMIC( );

};


struct RawStreamReader : public IStreamReader {

    RawStreamReader( int const fd, size_t const alignment = 0 );
    virtual ~RawStreamReader( );

    virtual uint64_t const  read_into_cache(
        uint8_t * const cache, uint64_t const cache_size
    );

protected:

    int			    _stream_handle;

};


struct GzStreamReader : public IStreamReader {

    GzStreamReader( int const fd );
    virtual ~GzStreamReader( );

    virtual uint64_t const  read_into_cache(
        uint8_t * const cache, uint64_t const cache_size
    );

private:

    gzFile		    _stream_handle;

};


struct Bz2StreamReader : public IStreamReader {

    Bz2StreamReader( int const fd );
    virtual ~Bz2StreamReader( );

    virtual uint64_t const  read_into_cache(
        uint8_t * const cache, uint64_t const cache_size
    );

private:

    FILE *		    _stream_handle;
    BZFILE *		    _block_handle;

};

#ifdef WITH_INTERNAL_METRICS
struct CacheSegmentPerformanceMetrics : public IPerformanceMetrics {

    enum {
        MKEY_TIME_WAITING_TO_SET_SA_BUFFER,
        MKEY_TIME_WAITING_TO_GET_SA_BUFFER,
        MKEY_TIME_WAITING_TO_FILL_FROM_STREAM,
        MKEY_TIME_FILLING_FROM_STREAM,
        MKEY_TIME_IN_SYNC_BARRIER
    };

    uint64_t	    numbytes_filled_from_stream;
    uint64_t	    numbytes_copied_from_ca_buffer;
    uint64_t	    numbytes_reserved_as_ca_buffer;
    uint64_t	    numbytes_copied_to_caller_buffer;
    uint64_t	    clock_nsecs_waiting_to_set_ca_buffer;
    uint64_t	    cpu_nsecs_waiting_to_set_ca_buffer;
    uint64_t	    clock_nsecs_waiting_to_get_ca_buffer;
    uint64_t	    cpu_nsecs_waiting_to_get_ca_buffer;
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
#endif

struct CacheManager {

    CacheManager(
        IStreamReader &	stream_reader,
        uint32_t const	number_of_threads,
        uint64_t const	cache_size,
        uint8_t const	trace_level =
            khmer:: get_active_config( ).get_input_buffer_trace_level( )
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

    uint64_t const	whereis_cursor( void );
    bool const		is_cursor_in_ca_buffer( void );
    void		split_at( uint64_t const pos );

    uint64_t const	get_fill_id( );

private:

    struct CacheSegment {

        bool				avail;
        uint32_t			thread_id;
        uint64_t			size;
        size_t				alignment;
        uint8_t *			memory;
        uint64_t			cursor;
        bool				cursor_in_ca_buffer;
        std:: string			ca_buffer;
        uint64_t			fill_id;
        bool				found_EOS;
#ifdef WITH_INTERNAL_METRICS
        CacheSegmentPerformanceMetrics	pmetrics;
#endif
        TraceLogger			trace_logger;

        CacheSegment(
            uint32_t const  thread_id,
            uint64_t const  size,
            size_t const    alignment = 0,
            uint8_t const   trace_level =
                khmer:: get_active_config( ).get_input_buffer_trace_level( )
        );
        ~CacheSegment( );

    }; // struct CacheSegment

    uint8_t		_trace_level;

    IStreamReader &	_stream_reader;

    uint32_t		_number_of_threads;
    ThreadIDMap		_thread_id_map;

    size_t		_alignment;
    uint64_t		_segment_size;
    CacheSegment **	_segments;
    uint32_t		_segment_ref_count;
    uint32_t		_segment_to_fill;
    uint64_t		_fill_counter;

    // Copyaside buffers.
    std:: map< uint64_t, std:: string >	_ca_buffers;
    uint32_t				_ca_spin_lock;

    // Extends or refills segment for current thread, as needed.
    void		_perform_segment_maintenance(
        CacheSegment & segment
    );

    uint64_t		_get_fill_counter_ATOMIC( );
    void		_increment_fill_counter_ATOMIC( );
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


struct Read {
    std:: string    name;
    std:: string    annotations;
    std:: string    sequence;
    std:: string    accuracy;
    // TODO? Add description field.
    uint64_t	    bytes_consumed;

    inline void reset ( ) {
        name.clear( );
        annotations.clear( );
        sequence.clear( );
        accuracy.clear( );
        bytes_consumed = 0;
    }
};

typedef std:: pair< Read, Read >	ReadPair;

#ifdef WITH_INTERNAL_METRICS
struct ParserPerformanceMetrics: public IPerformanceMetrics {

    uint64_t	    numlines_copied;
    uint64_t	    numreads_parsed_total;
    uint64_t	    numreads_parsed_valid;

    ParserPerformanceMetrics( );
    virtual ~ParserPerformanceMetrics( );

    virtual void    accumulate_timer_deltas( uint32_t metrics_key );

};
#endif


struct IParser {

    enum {
        PAIR_MODE_ALLOW_UNPAIRED = 0,
        PAIR_MODE_IGNORE_UNPAIRED,
        PAIR_MODE_ERROR_ON_UNPAIRED
    };

    static IParser * const  get_parser(
        std:: string const 	&ifile_name,
        uint32_t const		number_of_threads   =
            khmer:: get_active_config( ).get_number_of_threads( ),
        uint64_t const		cache_size	    =
            khmer:: get_active_config( ).get_reads_input_buffer_size( ),
        uint8_t const		trace_level	    =
            khmer:: get_active_config( ).get_reads_parser_trace_level( )
    );

    IParser(
        IStreamReader	&stream_reader,
        uint32_t const	number_of_threads   =
            khmer:: get_active_config( ).get_number_of_threads( ),
        uint64_t const	cache_size	    =
            khmer:: get_active_config( ).get_reads_input_buffer_size( ),
        uint8_t const	trace_level	    =
            khmer:: get_active_config( ).get_reads_parser_trace_level( )
    );
    virtual ~IParser( );

    inline int		uuid( ) {
        return _uuid;
    }

    inline bool		is_complete( ) {
        return !_cache_manager.has_more_data( ) && !_get_state( ).buffer_rem;
    }

    // Note: 'get_next_read' exists for legacy reasons.
    //	     In the long term, it should be eliminated in favor of direct use of
    //	     'imprint_next_read'. A potentially costly copy-by-value happens
    //	     upon return.
    // TODO: Eliminate all calls to 'get_next_read'.
    inline Read		get_next_read( ) {
        Read the_read;
        imprint_next_read( the_read );
        return the_read;
    }
    virtual void	imprint_next_read( Read &the_read );

    virtual void	imprint_next_read_pair(
        ReadPair &the_read_pair,
        uint8_t mode = PAIR_MODE_ERROR_ON_UNPAIRED
    );

protected:

    struct ParserState {

        // TODO: Set buffer size from Config.
        static uint64_t const	    BUFFER_SIZE		= 127;

        bool			    at_start;
        uint64_t		    fill_id;

        std:: string		    line;
        bool			    need_new_line;

        uint8_t			    buffer[ BUFFER_SIZE + 1 ];
        uint64_t		    buffer_pos;
        uint64_t		    buffer_rem;
#ifdef WITH_INTERNAL_METRICS
        ParserPerformanceMetrics    pmetrics;
#endif
        TraceLogger		    trace_logger;

        ParserState( uint32_t const thread_id, uint8_t const trace_level );
        ~ParserState( );

    }; // struct ParserState

    // TODO: Use a 16-octet IETF RFC 4122 UUID or equivalent.
    int			_uuid;

    uint8_t		_trace_level;

    CacheManager	_cache_manager;

    uint32_t		_number_of_threads;
    ThreadIDMap		_thread_id_map;
    bool		_unithreaded;

    ParserState **	_states;

    regex_t		_re_read_2_nosub;
    regex_t		_re_read_1;
    regex_t		_re_read_2;

    void		_copy_line( ParserState &state );

    virtual void	_parse_read( ParserState &, Read & )	    = 0;

#if (0)
    void		_imprint_next_read_pair_in_allow_mode(
        ReadPair &the_read_pair
    );
#endif
    void		_imprint_next_read_pair_in_ignore_mode(
        ReadPair &the_read_pair
    );
    void		_imprint_next_read_pair_in_error_mode(
        ReadPair &the_read_pair
    );
    bool		_is_valid_read_pair(
        ReadPair &the_read_pair, regmatch_t &match_1, regmatch_t &match_2
    );

    inline ParserState	&_get_state( ) {
        uint32_t	thread_id	= _thread_id_map.get_thread_id( );
        ParserState *	state_PTR	= NULL;

        assert( NULL != _states );

        state_PTR = _states[ thread_id ];
        if (NULL == state_PTR) {
            _states[ thread_id ]    =
                new ParserState( thread_id, _trace_level );
            state_PTR		    = _states[ thread_id ];
        }

        return *state_PTR;
    }

}; // struct IParser


struct FastaParser : public IParser {

    FastaParser(
        IStreamReader &	stream_reader,
        uint32_t const	number_of_threads,
        uint64_t const	cache_size,
        uint8_t const	trace_level
    );
    virtual ~FastaParser( );

    virtual void    _parse_read( ParserState &, Read &);

};


struct FastqParser : public IParser {

    FastqParser(
        IStreamReader &	stream_reader,
        uint32_t const	number_of_threads,
        uint64_t const	cache_size,
        uint8_t const	trace_level
    );
    virtual ~FastqParser( );

    virtual void    _parse_read( ParserState &, Read &);

};

inline PartitionID _parse_partition_id(std::string name)
{
    PartitionID p = 0;
    const char * s = name.c_str() + name.length() - 1;
    assert(*(s + 1) == (unsigned int) NULL);

    while(*s != '\t' && s >= name.c_str()) {
        s--;
    }

    if (*s == '\t') {
        p = (PartitionID) atoi(s + 1);
    } else {
        std::cerr << "consume_partitioned_fasta barfed on read "  << name << "\n";
        assert(0);
    }

    return p;
}



} // namespace read_parsers


} // namespace khmer


#endif // READ_PARSERS_HH

// vim: set ft=cpp sts=4 sw=4 tw=80:
