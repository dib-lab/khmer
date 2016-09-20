/*
This file is part of khmer, https://github.com/dib-lab/khmer/, and is
Copyright (C) 2012-2015, Michigan State University.
Copyright (C) 2015-2016, The Regents of the University of California.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are
met:

    * Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.

    * Redistributions in binary form must reproduce the above
      copyright notice, this list of conditions and the following
      disclaimer in the documentation and/or other materials provided
      with the distribution.

    * Neither the name of the Michigan State University nor the names
      of its contributors may be used to endorse or promote products
      derived from this software without specific prior written
      permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
"AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
LICENSE (END)

Contact: khmer-project@idyll.org
*/

#include <fstream>
#include "khmer_exception.hh"
#include "read_parsers.hh"

namespace khmer
{

namespace read_parsers
{

void ReadParser::_init()
{
    int regex_rc =
        regcomp(
            &_re_read_2_nosub,
            // ".+(/2| 2:[YN]:[[:digit:]]+:[[:alpha:]]+)$",
            "^.+(/2| 2:[YN]:[[:digit:]]+:[[:alpha:]]+).{0}",
            REG_EXTENDED | REG_NOSUB
        );
    if (regex_rc) {
        throw khmer_exception("Could not compile R2 nosub regex");
    }
    regex_rc =
        regcomp(
            &_re_read_1,
            "^.+(/1| 1:[YN]:[[:digit:]]+:[[:alpha:]]+).{0}", REG_EXTENDED
        );
    if (regex_rc) {
        throw khmer_exception("Could not compile R1 regex");
    }
    regex_rc =
        regcomp(
            &_re_read_2,
            "^.+(/2| 2:[YN]:[[:digit:]]+:[[:alpha:]]+).{0}", REG_EXTENDED
        );
    if (regex_rc) {
        throw khmer_exception("Could not compile R2 regex");
    }
}

template<typename ParseFunctor>
ReadParser::ReadParser(ParseFunctor pf)
        : _parser(pf), _num_reads(0), _have_qualities(false)
{
    _init();
}

ReadParser::ReadParser(ReadParser& other)
        : _parser(other._parser),
          _num_reads(other._num_reads),
          _have_qualities(other._have_qualities)
{
    _init();
}

ReadParser::~ReadParser()
{
    regfree(&_re_read_2_nosub);
    regfree(&_re_read_1);
    regfree(&_re_read_2);
}

void ReadParser::imprint_next_read(Read &read)
{
    parser(read);
}

void ReadParser::imprint_next_read_pair(ReadPair &pair, uint8_t mode)
{
    if (mode == IParser::PAIR_MODE_IGNORE_UNPAIRED) {
        _imprint_next_read_pair_in_ignore_mode(pair);
    }
    else if (mode == IParser::PAIR_MODE_ERROR_ON_UNPAIRED) {
        _imprint_next_read_pair_in_error_mode(pair);
    }
#if (0)
    else if (mode == IParser::PAIR_MODE_ALLOW_UNPAIRED) {
        _imprint_next_read_pair_in_allow_mode(pair);
    }
#endif
    else {
        std::ostringstream oss;
        oss << "Unknown pair reading mode: " << mode;
        throw UnknownPairReadingMode(oss.str());
    }
}

size_t ReadParser::get_num_reads()
{
    return _num_reads;
}

bool ReadParser::is_complete()
{
    return parser.is_complete();
}

#if (0)
void
ReadParser::_imprint_next_read_pair_in_allow_mode(ReadPair& pair)
{
    // TODO: Implement.
    //	     Probably need caching of reads between invocations
    //	     and the ability to return pairs which are half empty.
}
#endif

void ReadParser::_imprint_next_read_pair_in_ignore_mode(ReadPair& pair)
{
    Read& read_1 = pair.first;
    Read& read_2 = pair.second;
    regmatch_t match_1, match_2;

    // Hunt for a read pair until one is found or end of reads is reached.
    while (true) {

        // Toss out all reads which are not marked as first of a pair.
        // Note: We let any exception, which flies out of the following,
        //	 pass through unhandled.
        while (true) {
            imprint_next_read(read_1);
            if (!regexec(&_re_read_1, read_1.name.c_str(), 1, &match_1, 0)) {
                break;
            }
        }

        // If first read of a pair was found, then insist upon second read.
        // If not found, then restart search for pair.
        // If found, then validate match.
        // If invalid pair, then restart search for pair.
        imprint_next_read(read_2);
        if (!regexec(&_re_read_2, read_2.name.c_str(), 1, &match_2, 0)) {
            if (_is_valid_read_pair(pair, match_1, match_2)) {
                break;
            }
        }

    } // while pair not found

} // _imprint_next_read_pair_in_ignore_mode

void ReadParser:: _imprint_next_read_pair_in_error_mode(ReadPair& pair)
{
    Read & read_1 = the_read_pair.first;
    Read & read_2 = the_read_pair.second;
    regmatch_t match_1, match_2;

    // Note: We let any exception, which flies out of the following,
    //	     pass through unhandled.
    imprint_next_read(read_1);
    imprint_next_read(read_2);

    // Is the first read really the first member of a pair?
    if (REG_NOMATCH == regexec(
                &_re_read_1, read_1.name.c_str(), 1, &match_1, 0
            )) {
        throw InvalidReadPair();
    }
    // Is the second read really the second member of a pair?
    if (REG_NOMATCH == regexec(
                &_re_read_2, read_2.name.c_str(), 1, &match_2, 0
            )) {
        throw InvalidReadPair();
    }

    // Is the pair valid?
    if (!_is_valid_read_pair( the_read_pair, match_1, match_2 )) {
        throw InvalidReadPair();
    }

} // _imprint_next_read_pair_in_error_mode

bool ReadParser::_is_valid_read_pair(
    ReadPair& pair,
    regmatch_t &match_1,
    regmatch_t &match_2
)
{
    return (match_1.rm_so == match_2.rm_so)
            && (match_1.rm_eo == match_2.rm_eo)
            && (the_read_pair.first.name.substr(0, match_1.rm_so)
                    ==	the_read_pair.second.name.substr(0, match_1.rm_so));
}


void FastxParser::_init()
{
    seqan::open(_stream, filename);
    if (!seqan::isGood(_private->stream)) {
        std::string message = "Could not open ";
        message = message + filename + " for reading.";
        throw InvalidStream(message);
    } else if (seqan::atEnd(_private->stream)) {
        std::string message = "File ";
        message = message + filename + " does not contain any sequences!";
        throw InvalidStream(message);
    }
    __asm__ __volatile__ ("" ::: "memory");
}

FastxParser::FastxParser() : _filename("-"), _spin_lock(0)
{
    _init();
}

FastxParser::FastxParser(std::string& infile) : _filename(infile), _spin_lock(0)
{
    _init();
}

FastxParser::~FastxParser()
{
    seqan::close(_stream);
}

bool FastxParser::is_complete()
{
    return !seqan::isGood(_private->stream) || seqan::atEnd(_private->stream);
}

void FastxParser::imprint_next_read(Read& read)
{
    read.reset();
    int ret = -1;
    const char *invalid_read_exc = NULL;
    while (!__sync_bool_compare_and_swap(&_spin_lock, 0, 1));
    bool atEnd = seqan::atEnd(_stream);
    if (!atEnd) {
        ret = seqan::readRecord(read.name, read.sequence, read.quality, _stream);
        if (ret == 0) {
            // Detect if we're parsing something w/ qualities on the first read
            // only
            if (_num_reads == 0 && read.quality.length() != 0) {
                _have_qualities = true;
            }

            // Handle error cases, or increment number of reads on success
            if (read.sequence.length() == 0) {
                invalid_read_exc = "Sequence is empty";
            } else if (_have_qualities && (read.sequence.length() != \
                                           read.quality.length())) {
                invalid_read_exc = "Sequence and quality lengths differ";
            } else {
                _num_reads++;
            }
        }
    }
    __asm__ __volatile__ ("" ::: "memory");
    _spin_lock = 0;
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

} // namespace read_parsers

} // namespace khmer
