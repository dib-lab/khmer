//
// This file is part of khmer, https://github.com/dib-lab/khmer/, and is
// Copyright (C) Michigan State University, 2009-2015. It is licensed under
// the three-clause BSD license; see LICENSE.
// Contact: khmer-project@idyll.org
//

#ifndef READ_PARSERS_HH
#define READ_PARSERS_HH

#include <regex.h>
#include <iostream>
#include <cstdlib>
#include "khmer.hh"

namespace khmer
{



namespace read_parsers
{

struct NoMoreReadsAvailable : public  khmer_exception {
    explicit NoMoreReadsAvailable(const char *msg) :
        khmer_exception(msg) {}
    NoMoreReadsAvailable() :
        khmer_exception("No more reads available in this stream.") {}
};

struct InvalidRead : public  khmer_exception {
    explicit InvalidRead(const char *msg) :
        khmer_exception(msg) {}
    InvalidRead() :
        khmer_exception("Invalid read") {}
};

struct UnknownPairReadingMode : public  khmer_exception {
    explicit UnknownPairReadingMode(const char *msg) :
        khmer_exception(msg) {}
    UnknownPairReadingMode() :
        khmer_exception("Unknown pair reading mode supplied.") {}
};

struct InvalidReadPair : public  khmer_exception {
    explicit InvalidReadPair(const char *msg) :
        khmer_exception(msg) {}
    InvalidReadPair() :
        khmer_exception("Invalid read pair detected.") {}
};

struct Read {
    std:: string    name;
    std:: string    annotations;
    std:: string    sequence;
    std:: string    quality;
    // TODO? Add description field.

    inline void reset ( )
    {
        name.clear( );
        annotations.clear( );
        sequence.clear( );
        quality.clear( );
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

    size_t		    get_num_reads()
    {
        return _num_reads;
    }

protected:

    size_t		_num_reads;
    bool        _have_qualities;
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
    struct Handle;
    Handle* _private;

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
