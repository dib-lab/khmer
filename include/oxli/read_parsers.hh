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
#include <algorithm>
#include <cstdlib>
#include <iostream>
#include <string>
#include <utility>
#include <memory>

#include "oxli.hh"
#include "oxli_exception.hh"


namespace seqan
{
    class SequenceStream; // forward dec seqan dep
}

namespace oxli
{

namespace read_parsers
{

struct NoMoreReadsAvailable : public  oxli_file_exception {
    explicit NoMoreReadsAvailable(const std::string& msg) :
        oxli_file_exception(msg) {}
    NoMoreReadsAvailable() :
        oxli_file_exception("No more reads available in this stream.") {}
};

struct InvalidRead : public  oxli_value_exception {
    explicit InvalidRead(const std::string& msg) :
        oxli_value_exception(msg) {}
    InvalidRead() :
        oxli_value_exception("Invalid FASTA/Q read") {}
};

struct UnknownPairReadingMode : public  oxli_value_exception {
    explicit UnknownPairReadingMode(const std::string& msg) :
        oxli_value_exception(msg) {}
    UnknownPairReadingMode() :
        oxli_value_exception("Unknown pair reading mode supplied.") {}
};

struct InvalidReadPair : public  oxli_value_exception {
    explicit InvalidReadPair(const std::string& msg) :
        oxli_value_exception(msg) {}
    InvalidReadPair() :
        oxli_value_exception("Invalid read pair detected.") {}
};


unsigned char _to_valid_dna(const unsigned char c);


struct Read {
    std::string name;
    std::string description;
    std::string sequence;
    std::string quality;
    std::string cleaned_seq;

    inline void reset()
    {
        name.clear();
        description.clear();
        sequence.clear();
        quality.clear();
        cleaned_seq.clear();
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

    // Compute cleaned_seq from sequence. Call this after changing sequence.
    inline void set_clean_seq()
    {
        cleaned_seq = std::string(sequence.size(), 0);
        std::transform(sequence.begin(), sequence.end(), cleaned_seq.begin(),
                       _to_valid_dna);
    }
};
typedef std::pair<Read, Read> ReadPair;


template<typename SeqIO>
class ReadParser
{
protected:
    std::unique_ptr<SeqIO> _parser;
    regex_t _re_read_2_nosub;
    regex_t _re_read_1;
    regex_t _re_read_2;
    void _init();

    ReadPair _get_next_read_pair_in_ignore_mode();
    ReadPair _get_next_read_pair_in_error_mode();
    bool _is_valid_read_pair(
        ReadPair &the_read_pair,
        regmatch_t &match_1,
        regmatch_t &match_2
    );

public:
    enum {
        PAIR_MODE_IGNORE_UNPAIRED,
        PAIR_MODE_ERROR_ON_UNPAIRED
    };

    ReadParser(std::unique_ptr<SeqIO> pf);

    ReadParser(ReadParser& other);
    ReadParser& operator=(ReadParser& other);

    ReadParser(ReadParser&&) noexcept;
    ReadParser& operator=(ReadParser&&) noexcept;

    virtual ~ReadParser();

    Read get_next_read();
    ReadPair get_next_read_pair(uint8_t mode = PAIR_MODE_ERROR_ON_UNPAIRED);

    size_t get_num_reads();
    bool is_complete();
    void close();
}; // class ReadParser


class FastxReader
{
private:
    std::string _filename;
    std::unique_ptr<seqan::SequenceStream> _stream;
    uint32_t _spin_lock;
    size_t _num_reads;
    bool _have_qualities;
    void _init();

public:
    FastxReader();
    FastxReader(const std::string& infile);

    FastxReader(FastxReader& other);
    FastxReader& operator=(FastxReader& other);

    FastxReader(FastxReader&&) noexcept;
    FastxReader& operator=(FastxReader&&) noexcept;

    ~FastxReader();

    Read get_next_read();
    bool is_complete();
    size_t get_num_reads();
    void close();
}; // class FastxReader


inline PartitionID _parse_partition_id(std::string name)
{
    PartitionID p = 0;
    const char * s = name.c_str() + name.length() - 1;
    if (!(*(s + 1) == (unsigned int) NULL)) {
        throw oxli_exception();
    }

    while(*s != '\t' && s >= name.c_str()) {
        s--;
    }

    if (*s == '\t') {
        p = (PartitionID) atoi(s + 1);
    } else {
        std::string err;
        err = "consume_partitioned_fasta cannot find partition ID for read ";
        err += name;

        throw oxli_value_exception(err);
    }

    return p;
}

// Alias for generic/templated ReadParser pointer
template<typename T> using ReadParserPtr = std::shared_ptr<ReadParser<T>>;
template<typename T> using WeakReadParserPtr = std::weak_ptr<ReadParser<T>>;

// Convenience function
template<typename SeqIO>
ReadParserPtr<SeqIO> get_parser(const std::string& filename);

// Alias for instantiated ReadParsers
typedef std::shared_ptr<ReadParser<FastxReader>> FastxParserPtr;
typedef std::weak_ptr<ReadParser<FastxReader>> WeakFastxParserPtr;

} // namespace read_parsers

} // namespace oxli

#endif // READ_PARSERS_HH
