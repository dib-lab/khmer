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
#ifndef READ_PARSERS_HH
#define READ_PARSERS_HH

#include <regex.h>
#include <stddef.h>
#include <stdint.h>
#include <cstdlib>
#include <iostream>
#include <string>
#include <utility>

#include "khmer.hh"
#include "khmer_exception.hh"

namespace khmer
{

namespace read_parsers
{

struct NoMoreReadsAvailable : public  khmer_file_exception {
    explicit NoMoreReadsAvailable(const std::string& msg) :
        khmer_file_exception(msg) {}
    NoMoreReadsAvailable() :
        khmer_file_exception("No more reads available in this stream.") {}
};

struct InvalidRead : public  khmer_value_exception {
    explicit InvalidRead(const std::string& msg) :
        khmer_value_exception(msg) {}
    InvalidRead() :
        khmer_value_exception("Invalid FASTA/Q read") {}
};

struct UnknownPairReadingMode : public  khmer_value_exception {
    explicit UnknownPairReadingMode(const std::string& msg) :
        khmer_value_exception(msg) {}
    UnknownPairReadingMode() :
        khmer_value_exception("Unknown pair reading mode supplied.") {}
};

struct InvalidReadPair : public  khmer_value_exception {
    explicit InvalidReadPair(const std::string& msg) :
        khmer_value_exception(msg) {}
    InvalidReadPair() :
        khmer_value_exception("Invalid read pair detected.") {}
};

class Read
{
public:
    std::string name;
    std::string annotations;
    std::string sequence;
    std::string quality;

    inline void reset()
    {
        name.clear();
        annotations.clear();
        sequence.clear();
        quality.clear();
    }

    inline void write_fastx(std::ostream& output)
    {
        if (quality.length() != 0) {
            output << "@" << name << '\n'
                   << sequence << '\n'
                   << "+" << '\n'
                   << quality << '\n';
        } else {
            output << ">" << name << '\n'
                   << sequence << '\n';
        }
    }
};
typedef std::pair<Read, Read> ReadPair;

template <typename ParseFunctor>
class ReadParser
{
protected:
    ParseFunctor parser;
    size_t _num_reads;
    bool _have_qualities;
    regex_t _re_read_2_nosub;
    regex_t _re_read_1;
    regex_t _re_read_2;

#if (0)
    void  _imprint_next_read_pair_in_allow_mode(ReadPair& pair);
#endif
    void  _imprint_next_read_pair_in_ignore_mode(ReadPair& pair);
    void  _imprint_next_read_pair_in_error_mode(ReadPair& pair);
    bool  _is_valid_read_pair(
        ReadPair& pair,
        regmatch_t &match_1,
        regmatch_t &match_2
    );

public:
    enum {
        PAIR_MODE_ALLOW_UNPAIRED = 0,
        PAIR_MODE_IGNORE_UNPAIRED,
        PAIR_MODE_ERROR_ON_UNPAIRED
    };

    // Constructors / destructors
    explicit ReadParser(ParseFunctor parser);
    explicit ReadParser(ReadParser& other);
    ReadParser& operator=(const ReadParser& other);
    virtual ~ReadParser();

    // Class methods
    void imprint_next_read(Read &read);
    void imprint_next_read_pair(
        ReadPair &pair,
        uint8_t mode = PAIR_MODE_ERROR_ON_UNPAIRED
    );
    size_t get_num_reads();
    virtual bool is_complete();

    // Note: 'get_next_read' exists for legacy reasons.
    //       In the long term, it should be eliminated in favor of direct use of
    //       'imprint_next_read'. A potentially costly copy-by-value happens
    //       upon return.
    // TODO: Eliminate all calls to 'get_next_read'.
    // Or switch to C++11 w/ move constructors
    inline Read get_next_read()
    {
        Read read;
        imprint_next_read(read);
        return read;
    }
}; // struct ReadParser

class FastxParser
{
private:
    std::string filename;
    seqan::SequenceStream stream;
    uint32_t seqan_spin_lock;

public:
    FastxParser();
    FastxParser(std::string& infile);
    FastxParser(FastxParser& other);
    FastxParser& operator=(const FastxParser& other);

    void operator()(Read &read);
    bool is_complete();
}

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
        std::cerr << "consume_partitioned_fasta barfed on read "
                  << name << "\n";
        throw khmer_exception();
    }

    return p;
}

} // namespace read_parsers

} // namespace khmer

#endif // READ_PARSERS_HH
