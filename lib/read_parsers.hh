//
// This file is part of khmer, http://github.com/ged-lab/khmer/, and is
// Copyright (C) Michigan State University, 2009-2013. It is licensed under
// the three-clause BSD license; see doc/LICENSE.txt.
// Contact: khmer-project@idyll.org
//

#ifndef READ_PARSERS_HH
#define READ_PARSERS_HH

#include <iostream>
#include <seqan/sequence.h>
#include <seqan/seq_io.h>
#include <seqan/stream.h>
#include <cstdarg>
#include <string>
#include <utility>
#include <stdlib.h>
#include <mutex>

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


class InvalidReadFileFormat: public khmer:: khmer_file_exception
{

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

class InvalidFASTAFileFormat: public InvalidReadFileFormat
{

public:

public:
    InvalidFASTAFileFormat(
        char const * reason	= NULL,
        char const * evidence	= NULL
    );

};

class InvalidFASTQFileFormat: public InvalidReadFileFormat
{

public:

public:
    InvalidFASTQFileFormat(
        char const * reason	= NULL,
        char const * evidence	= NULL
    );

};

struct CacheSegmentUnavailable : public  khmer_exception {
};

struct CacheSegmentBoundaryViolation : public  khmer_exception {
};

struct InvalidCacheSizeRequested : public  khmer_exception {
};

struct NoMoreReadsAvailable : public  khmer_exception {
};

struct UnknownPairReadingMode : public  khmer_exception {
};

struct InvalidReadPair : public  khmer_exception {
};

struct Read {
    std:: string    name;
    std:: string    annotations;
    std:: string    sequence;
    std:: string    accuracy;
    // TODO? Add description field.

    inline void reset ( )
    {
        name.clear( );
        annotations.clear( );
        sequence.clear( );
        accuracy.clear( );
    }
};

typedef std:: pair< Read, Read >	ReadPair;

struct IParser {

    enum {
        PAIR_MODE_ALLOW_UNPAIRED = 0,
        PAIR_MODE_IGNORE_UNPAIRED,
        PAIR_MODE_ERROR_ON_UNPAIRED
    };

    static IParser * const  get_parser(
        std:: string const 	&ifile_name
    );

    IParser( );
    virtual ~IParser( );

    virtual bool		is_complete( ) = 0;

    // Note: 'get_next_read' exists for legacy reasons.
    //	     In the long term, it should be eliminated in favor of direct use of
    //	     'imprint_next_read'. A potentially costly copy-by-value happens
    //	     upon return.
    // TODO: Eliminate all calls to 'get_next_read'.
    // Or switch to C++11 w/ move constructors
    inline Read		get_next_read( )
    {
        Read the_read;
        imprint_next_read( the_read );
        return the_read;
    }
    virtual void	imprint_next_read( Read &the_read ) = 0;

    virtual void	imprint_next_read_pair(
        ReadPair &the_read_pair,
        uint8_t mode = PAIR_MODE_ERROR_ON_UNPAIRED
    );

protected:

    regex_t		_re_read_2_nosub;
    regex_t		_re_read_1;
    regex_t		_re_read_2;

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

}; // struct IParser

class SeqAnParser : public IParser
{

public:
    SeqAnParser( const char * filename );
    ~SeqAnParser( );

    bool is_complete( );
    void imprint_next_read(Read &the_read);

private:
    seqan::SequenceStream _stream;
    std::mutex _imprint_mutex;

};

inline PartitionID _parse_partition_id(std::string name)
{
    PartitionID p = 0;
    const char * s = name.c_str() + name.length() - 1;
    if (!(*(s + 1) == (unsigned int) NULL)) {
        throw khmer_exception();
    }

    while(*s != '\t' && s >= name.c_str()) {
        s--;
    }

    if (*s == '\t') {
        p = (PartitionID) atoi(s + 1);
    } else {
        std::cerr << "consume_partitioned_fasta barfed on read "  << name << "\n";
        throw khmer_exception();
    }

    return p;
}



} // namespace read_parsers


} // namespace khmer


#endif // READ_PARSERS_HH

// vim: set ft=cpp sts=4 sw=4 tw=80:
