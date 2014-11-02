//
// This file is part of khmer, http://github.com/ged-lab/khmer/, and is
// Copyright (C) Michigan State University, 2009-2013. It is licensed under
// the three-clause BSD license; see doc/LICENSE.txt.
// Contact: khmer-project@idyll.org
//

#include "read_parsers.hh"

#include <cstring>
#include "khmer_exception.hh"

namespace khmer
{


namespace read_parsers
{

SeqAnParser::SeqAnParser( char const * filename ) : IParser( )
{
    seqan::open(_stream, filename);
    if (!seqan::isGood(_stream) || seqan::atEnd(_stream)) {
        throw InvalidStreamHandle();
    }
    pthread_mutex_init(&_imprint_mutex, NULL);
}

bool SeqAnParser::is_complete()
{
    return !seqan::isGood(_stream) || seqan::atEnd(_stream);
}

void SeqAnParser::imprint_next_read(Read &the_read)
{
    the_read.reset();
    pthread_mutex_lock(&_imprint_mutex);
    int ret = seqan::readRecord(the_read.name, the_read.sequence,
                                the_read.accuracy, _stream);
    pthread_mutex_unlock(&_imprint_mutex);
    if (ret != 0) {
        throw NoMoreReadsAvailable();
    }
}

SeqAnParser::~SeqAnParser()
{
    pthread_mutex_destroy(&_imprint_mutex);
    seqan::close(_stream);
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
