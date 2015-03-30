//
// This file is part of khmer, http://github.com/ged-lab/khmer/, and is
// Copyright (C) Michigan State University, 2009-2013. It is licensed under
// the three-clause BSD license; see LICENSE.
// Contact: khmer-project@idyll.org
//

#include "read_parsers.hh"

#include <cstring>
#include <string.h>
#include "khmer_exception.hh"
#include <seqan/sequence.h>
#include <seqan/seq_io.h>
#include <seqan/stream.h>
#include <pthread.h>
#include "seqan/basic.h"

namespace khmer
{


namespace read_parsers
{

struct SeqAnParser::Handle {
    seqan::SeqFileIn file;
    uint32_t seqan_spin_lock;
};

SeqAnParser::SeqAnParser( char const * filename ) : IParser( )
{
    _private = new SeqAnParser::Handle();
    bool ret = false;
    if (strlen(filename) == 1 && strcmp(filename, "-") == 0) {
	ret = seqan::open(_private->file, std::cin);
    } else {
	ret = seqan::open(_private->file, filename);
    }
    if (!ret) {
	std::string message = "Could not open ";
        message = message + filename + " for reading.";
        throw InvalidStreamHandle(message.c_str());
    } else if (seqan::atEnd(_private->file)) {
        std::string message = "File ";
        message = message + filename + " does not contain any sequences!";
        throw InvalidStreamHandle(message.c_str());
    }
    __asm__ __volatile__ ("" ::: "memory");
    _private->seqan_spin_lock = 0;
}

bool SeqAnParser::is_complete()
{
    while (!__sync_bool_compare_and_swap(& _private->seqan_spin_lock, 0, 1));
    bool end = seqan::atEnd(_private->file);
    __asm__ __volatile__ ("" ::: "memory");
    _private->seqan_spin_lock = 0;
    return end;
}

void SeqAnParser::imprint_next_read(Read &the_read)
{
    the_read.reset();
    while (!__sync_bool_compare_and_swap(& _private->seqan_spin_lock, 0, 1));
    try {
	seqan::readRecord(the_read.name, the_read.sequence, the_read.quality,
		_private->file);
    } catch (std::exception &e) {
	__asm__ __volatile__ ("" ::: "memory");
	_private->seqan_spin_lock = 0;
	if (seqan::atEnd(_private->file)) {
	    throw NoMoreReadsAvailable();
	} else {
	    throw StreamReadError(e.what());
	}
    }
    __asm__ __volatile__ ("" ::: "memory");
    _private->seqan_spin_lock = 0;
}

SeqAnParser::~SeqAnParser()
{
    while (!__sync_bool_compare_and_swap(& _private->seqan_spin_lock, 0, 1));
    seqan::close(_private->file);
    __asm__ __volatile__ ("" ::: "memory");
    _private->seqan_spin_lock = 0;
    delete _private;
}

IParser * const
IParser::
get_parser(
    std:: string const	    &ifile_name
)
{

    return new SeqAnParser(ifile_name.c_str());
}


IParser::
IParser(
)
{
    int regex_rc =
        regcomp(
            &_re_read_2_nosub,
            // ".+(/2| 2:[YN]:[[:digit:]]+:[[:alpha:]]+)$",
            "^.+(/2| 2:[YN]:[[:digit:]]+:[[:alpha:]]+).{0}",
            REG_EXTENDED | REG_NOSUB
        );
    if (regex_rc) {
        throw khmer_exception();
    }
    regex_rc =
        regcomp(
            &_re_read_1,
            "^.+(/1| 1:[YN]:[[:digit:]]+:[[:alpha:]]+).{0}", REG_EXTENDED
        );
    if (regex_rc) {
        throw khmer_exception();
    }
    regex_rc =
        regcomp(
            &_re_read_2,
            "^.+(/2| 2:[YN]:[[:digit:]]+:[[:alpha:]]+).{0}", REG_EXTENDED
        );
    if (regex_rc) {
        throw khmer_exception();
    }
    _num_reads = 0;
    _have_qualities = false;
}

IParser::
~IParser( )
{
    regfree( &_re_read_2_nosub );
    regfree( &_re_read_1 );
    regfree( &_re_read_2 );
}

void
IParser::
imprint_next_read_pair( ReadPair &the_read_pair, uint8_t mode )
{
    switch (mode) {
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
    while (true) {

        // Toss out all reads which are not marked as first of a pair.
        // Note: We let any exception, which flies out of the following,
        //	 pass through unhandled.
        while (true) {
            imprint_next_read( read_1 );
            if (!regexec(
                        &_re_read_1, read_1.name.c_str( ), 1, &match_1, 0
                    )) {
                break;
            }
        }

        // If first read of a pair was found, then insist upon second read.
        // If not found, then restart search for pair.
        // If found, then validate match.
        // If invalid pair, then restart search for pair.
        imprint_next_read( read_2 );
        if (!regexec(
                    &_re_read_2, read_2.name.c_str( ), 1, &match_2, 0
                )) {
            if (_is_valid_read_pair( the_read_pair, match_1, match_2 )) {
                break;
            }
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
            )) {
        throw InvalidReadPair( );
    }
    // Is the second read really the second member of a pair?
    if (REG_NOMATCH == regexec(
                &_re_read_2, read_2.name.c_str( ), 1, &match_2, 0
            )) {
        throw InvalidReadPair( );
    }

    // Is the pair valid?
    if (!_is_valid_read_pair( the_read_pair, match_1, match_2 )) {
        throw InvalidReadPair( );
    }

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
