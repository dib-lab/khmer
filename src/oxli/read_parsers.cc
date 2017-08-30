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
#include <utility>  
#include "seqan/seq_io.h" // IWYU pragma: keep
#include "seqan/sequence.h" // IWYU pragma: keep
#include "seqan/stream.h" // IWYU pragma: keep
#include "oxli/oxli_exception.hh"
#include "oxli/read_parsers.hh"


namespace oxli
{

namespace read_parsers
{

unsigned char _to_valid_dna(const unsigned char c)
{
    switch(c) {
    case 'A':
    case 'C':
    case 'G':
    case 'T':
        return c;
    case 'a':
    case 'c':
    case 'g':
    case 't':
        return toupper(c);
    default:
        return 'A';
    }
}

template<typename SeqIO>
void ReadParser<SeqIO>::_init()
{
    int regex_rc =
        regcomp(
            &_re_read_2_nosub,
            // ".+(/2| 2:[YN]:[[:digit:]]+:[[:alpha:]]+)$",
            "^.+(/2| 2:[YN]:[[:digit:]]+:[[:alpha:]]+).{0}",
            REG_EXTENDED | REG_NOSUB
        );
    if (regex_rc) {
        throw oxli_exception("Could not compile R2 nosub regex");
    }
    regex_rc =
        regcomp(
            &_re_read_1,
            "^.+(/1| 1:[YN]:[[:digit:]]+:[[:alpha:]]+).{0}", REG_EXTENDED
        );
    if (regex_rc) {
        throw oxli_exception("Could not compile R1 regex");
    }
    regex_rc =
        regcomp(
            &_re_read_2,
            "^.+(/2| 2:[YN]:[[:digit:]]+:[[:alpha:]]+).{0}", REG_EXTENDED
        );
    if (regex_rc) {
        throw oxli_exception("Could not compile R2 regex");
    }
}

template<typename SeqIO>
ReadPair ReadParser<SeqIO>::_get_next_read_pair_in_ignore_mode()
{
    ReadPair pair;
    regmatch_t match_1, match_2;

    // Hunt for a read pair until one is found or end of reads is reached.
    while (true) {

        // Toss out all reads which are not marked as first of a pair.
        // Note: We let any exception, which flies out of the following,
        //	 pass through unhandled.
        while (true) {
            pair.first = get_next_read();
            if (!regexec(
                        &_re_read_1, pair.first.name.c_str( ), 1, &match_1, 0
                    )) {
                break;
            }
        }

        // If first read of a pair was found, then insist upon second read.
        // If not found, then restart search for pair.
        // If found, then validate match.
        // If invalid pair, then restart search for pair.
        pair.second = get_next_read();
        if (!regexec(
                    &_re_read_2, pair.second.name.c_str( ), 1, &match_2, 0
                )) {
            if (_is_valid_read_pair(pair, match_1, match_2)) {
                break;
            }
        }

    } // while pair not found

    return pair;
} // _get_next_read_pair_in_ignore_mode

template<typename SeqIO>
ReadPair ReadParser<SeqIO>::_get_next_read_pair_in_error_mode()
{
    ReadPair pair;
    regmatch_t match_1, match_2;

    // Note: We let any exception, which flies out of the following,
    //	     pass through unhandled.
    pair.first = get_next_read();
    pair.second = get_next_read();

    // Is the first read really the first member of a pair?
    if (REG_NOMATCH == regexec(
                &_re_read_1, pair.first.name.c_str( ), 1, &match_1, 0
            )) {
        throw InvalidReadPair( );
    }
    // Is the second read really the second member of a pair?
    if (REG_NOMATCH == regexec(
                &_re_read_2, pair.second.name.c_str( ), 1, &match_2, 0
            )) {
        throw InvalidReadPair( );
    }

    // Is the pair valid?
    if (!_is_valid_read_pair(pair, match_1, match_2)) {
        throw InvalidReadPair( );
    }

    return pair;
} // _get_next_read_pair_in_error_mode

template<typename SeqIO>
bool ReadParser<SeqIO>::_is_valid_read_pair(
    ReadPair &the_read_pair, regmatch_t &match_1, regmatch_t &match_2
)
{
    return	(match_1.rm_so == match_2.rm_so)
            &&	(match_1.rm_eo == match_2.rm_eo)
            &&	(	the_read_pair.first.name.substr( 0, match_1.rm_so )
                    ==	the_read_pair.second.name.substr( 0, match_1.rm_so ));
}


template<typename SeqIO>
ReadParser<SeqIO>::ReadParser(std::unique_ptr<SeqIO> pf)
{
    _parser = std::move(pf);
    _init();
}


template<typename SeqIO>
ReadParser<SeqIO>::ReadParser(ReadParser<SeqIO>& other)
{
    _parser = std::move(other._parser);
    _init();
}

template<typename SeqIO>
ReadParser<SeqIO>&
ReadParser<SeqIO>::operator=(ReadParser<SeqIO>& other) {
    _parser = std::move(other._parser);
    return *this;
}

template<typename SeqIO>
ReadParser<SeqIO>::ReadParser(ReadParser<SeqIO>&&) noexcept {}


template<typename SeqIO>
ReadParser<SeqIO>::~ReadParser()
{
    regfree(&_re_read_2_nosub);
    regfree(&_re_read_1);
    regfree(&_re_read_2);
}

template<typename SeqIO>
Read ReadParser<SeqIO>::get_next_read()
{
    return _parser->get_next_read();
}

template<typename SeqIO>
ReadPair ReadParser<SeqIO>::get_next_read_pair(uint8_t mode)
{
    if (mode == ReadParser<SeqIO>::PAIR_MODE_IGNORE_UNPAIRED) {
        return _get_next_read_pair_in_ignore_mode();
    } else if (mode == ReadParser<SeqIO>::PAIR_MODE_ERROR_ON_UNPAIRED) {
        return _get_next_read_pair_in_error_mode();
    } else {
        std::ostringstream oss;
        oss << "Unknown pair reading mode: " << mode;
        throw UnknownPairReadingMode(oss.str());
    }
}

template<typename SeqIO>
size_t ReadParser<SeqIO>::get_num_reads()
{
    return _parser->get_num_reads();
}

template<typename SeqIO>
bool ReadParser<SeqIO>::is_complete()
{
    return _parser->is_complete();
}

template<typename SeqIO>
void ReadParser<SeqIO>::close()
{
    _parser->close();
}

void FastxReader::_init()
{
    _stream = std::unique_ptr<seqan::SequenceStream>(new seqan::SequenceStream());
    seqan::open(*_stream, _filename.c_str());
    if (!seqan::isGood(*_stream)) {
        std::string message = "File ";
        message = message + _filename + " contains badly formatted sequence";
        message = message + " or does not exist.";
        throw InvalidStream(message);
    } else if (seqan::atEnd(*_stream)) {
        std::string message = "File ";
        message = message + _filename + " does not contain any sequences!";
        throw InvalidStream(message);
    }
    __asm__ __volatile__ ("" ::: "memory");
}

FastxReader::FastxReader()
    : _filename("-"), _spin_lock(0), _num_reads(0), _have_qualities(false)
{
    _init();
}

FastxReader::FastxReader(const std::string& infile)
    : _filename(infile),
      _spin_lock(0),
      _num_reads(0),
      _have_qualities(false)
{
    _init();
}

FastxReader::FastxReader(FastxReader& other) {
    _filename = other._filename;
    _spin_lock = other._spin_lock;
    _num_reads = other._num_reads;
    _have_qualities = other._have_qualities;
    _stream = std::move(other._stream);
}

//FastxReader::FastxReader(const FastxReader& other) = delete;
FastxReader& FastxReader::operator= (FastxReader& other) {
    _filename = other._filename;
    _spin_lock = other._spin_lock;
    _num_reads = other._num_reads;
    _have_qualities = other._have_qualities;
    _stream = std::move(other._stream);
    return *this;
}

FastxReader::FastxReader(FastxReader&&) noexcept {}

FastxReader::~FastxReader()
{
    seqan::close(*_stream);
}

bool FastxReader::is_complete()
{
    return !seqan::isGood(*_stream) || seqan::atEnd(*_stream);
}

size_t FastxReader::get_num_reads()
{
    return _num_reads;
}

void FastxReader::close()
{
    seqan::close(*_stream);
}

Read FastxReader::get_next_read()
{
    Read read;
    int ret = -1;
    const char *invalid_read_exc = NULL;
    while (!__sync_bool_compare_and_swap(&_spin_lock, 0, 1));
    bool atEnd = seqan::atEnd(*_stream);
    if (!atEnd) {
        ret = seqan::readRecord(read.name, read.sequence, read.quality, *_stream);
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
    return read;
}

template<typename SeqIO>
ReadParserPtr<SeqIO> get_parser(const std::string& filename)
{
    return ReadParserPtr<SeqIO>(
               new ReadParser<SeqIO>(
                   std::unique_ptr<SeqIO>(new SeqIO(filename))
               )
           );
}

// All template instantiations used in the codebase must be declared here.
template class ReadParser<FastxReader>;
template FastxParserPtr get_parser<FastxReader>(const std::string& filename);

} // namespace read_parsers

} // namespace oxli

// vim: set ft=cpp sts=4 sw=4 tw=80:

