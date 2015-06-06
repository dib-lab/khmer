//
// This file is part of khmer, https://github.com/dib-lab/khmer/, and is
// Copyright (C) Michigan State University, 2009-2015. It is licensed under
// the three-clause BSD license; see LICENSE.
// Contact: khmer-project@idyll.org
//

#include "read_parsers.hh"

#include <cstring>
#include "khmer_exception.hh"
#include <seqan/sequence.h>
#include <seqan/seq_io.h>
#include <seqan/stream.h>
#include <pthread.h>

namespace khmer
{


namespace read_parsers
{

struct SeqAnParser::Handle {
    seqan::SequenceStream stream;
    uint32_t seqan_spin_lock;
};

SeqAnParser::SeqAnParser( char const * filename ) : IParser( )
{
    _private = new SeqAnParser::Handle();
    seqan::open(_private->stream, filename);
    if (!seqan::isGood(_private->stream)) {
        std::string message = "Could not open ";
        message = message + filename + " for reading.";
        throw InvalidStreamHandle(message.c_str());
    } else if (seqan::atEnd(_private->stream)) {
        std::string message = "File ";
        message = message + filename + " does not contain any sequences!";
        throw InvalidStreamHandle(message.c_str());
    }
    __asm__ __volatile__ ("" ::: "memory");
    _private->seqan_spin_lock = 0;
}

bool SeqAnParser::is_complete()
{
    return !seqan::isGood(_private->stream) || seqan::atEnd(_private->stream);
}

void SeqAnParser::imprint_next_read(Read &the_read)
{
    the_read.reset();
    int ret = -1;
    const char *invalid_read_exc = NULL;
    while (!__sync_bool_compare_and_swap(& _private->seqan_spin_lock, 0, 1));
    bool atEnd = seqan::atEnd(_private->stream);
    if (!atEnd) {
        ret = seqan::readRecord(the_read.name, the_read.sequence,
                                the_read.quality, _private->stream);
        if (ret == 0) {
            // Detect if we're parsing something w/ qualities on the first read
            // only
            if (_num_reads == 0 && the_read.quality.length() != 0) {
                _have_qualities = true;
            }

            // Handle error cases, or increment number of reads on success
            if (the_read.sequence.length() == 0) {
                invalid_read_exc = "Sequence is empty";
            } else if (_have_qualities && (the_read.sequence.length() != \
                                           the_read.quality.length())) {
                invalid_read_exc = "Sequence and quality lengths differ";
            } else {
                _num_reads++;
            }
        }
    }
    __asm__ __volatile__ ("" ::: "memory");
    _private->seqan_spin_lock = 0;
    // Throw any error in the read, even if we're at the end
    if (invalid_read_exc != NULL) {
        throw InvalidRead(invalid_read_exc);
    }
    // Throw NoMoreReadsAvailable if none of the above errors were raised, even
    // if ret == 0
    if (atEnd) {
        throw NoMoreReadsAvailable();
    }
    // Catch-all error in readRecord that isn't one of the above
    if (ret != 0) {
        throw StreamReadError();
    }
}

SeqAnParser::~SeqAnParser()
{
    seqan::close(_private->stream);
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
