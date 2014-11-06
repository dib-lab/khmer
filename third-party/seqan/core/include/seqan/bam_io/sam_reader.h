// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2013, Knut Reinert, FU Berlin
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are met:
//
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//     * Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in the
//       documentation and/or other materials provided with the distribution.
//     * Neither the name of Knut Reinert or the FU Berlin nor the names of
//       its contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL KNUT REINERT OR THE FU BERLIN BE LIABLE
// FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
// DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
// SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
// CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
// LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY
// OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
// DAMAGE.
//
// ==========================================================================
// Author: Manuel Holtgrewe <manuel.holtgrewe@fu-berlin.de>
// ==========================================================================

#ifndef CORE_INCLUDE_SEQAN_BAM_IO_SAM_READER_H_
#define CORE_INCLUDE_SEQAN_BAM_IO_SAM_READER_H_

#include <memory>

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// Implementation of reading BAM files.

class SamReader_ :
    public XamReader_
{
public:
    // The size of the file in bytes.
    __int64 _fileSize;

    // Pointer to the input stream to read from.  Used such that we can also read from stdin.
    std::ifstream _fstream;
    // The file stream to read from if we do not read from stdout.  Pointed to by _stream.
    std::istream * _stream;

    // Record reader to use for parsing.
#if __cplusplus <= 199711L
    // C++98
    std::auto_ptr<RecordReader<std::istream, SinglePass<> > > _reader;
#else  // #if __cplusplus <= 199711L
    // C++11
    std::unique_ptr<RecordReader<std::istream, SinglePass<> > > _reader;
#endif  // #if __cplusplus <= 199711L

    SamReader_() :
        XamReader_(), _fileSize(0), _stream(0),
#if __cplusplus <= 199711L
        // C++98
        _reader(0)
#else  // #if __cplusplus <= 199711L
        // C++11
        _reader()
#endif  // #if __cplusplus <= 199711L
    {}

    SamReader_(CharString const & filename);

    // XamReader_ interface.

    virtual int open(CharString const & filename);
    virtual bool isGood();
    virtual bool atEnd();
    virtual int readHeader(BamHeader & header, BamIOContext<StringSet<CharString> > & context);
    virtual int readRecord(BamAlignmentRecord & record, BamIOContext<StringSet<CharString> > & context);
    virtual int close();

    virtual __int64 fileSize() const;
    virtual __int64 positionInFile() const;
};

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Member Function SamReader_::SamReader_()
// ----------------------------------------------------------------------------

inline SamReader_::SamReader_(CharString const & filename) :
    XamReader_(filename), _stream(0), _reader()
{
    this->open(_filename);
}

// ----------------------------------------------------------------------------
// Member Function SamReader_::open()
// ----------------------------------------------------------------------------

inline int SamReader_::open(CharString const & filename)
{
    // Open file.
    if (filename == "-")
    {
        this->_stream = &std::cin;
    }
    else
    {
        this->_stream = &this->_fstream;
        this->_fstream.open(toCString(filename), std::ios::binary | std::ios::in);
        if (!this->_fstream.good())
            return 1;

        // Determine file size.
        this->_stream->seekg(0, std::ios::end);
        this->_fileSize = this->_stream->tellg();
        this->_stream->seekg(0, std::ios::beg);
    }
    this->_reader.reset(new RecordReader<std::istream, SinglePass<> >(*this->_stream));

    return 0;
}

// ----------------------------------------------------------------------------
// Member Function SamReader_::isGood()
// ----------------------------------------------------------------------------

inline bool SamReader_::isGood()
{
    return this->_stream->good();
}

// ----------------------------------------------------------------------------
// Member Function SamReader_::atEnd()
// ----------------------------------------------------------------------------

inline bool SamReader_::atEnd()
{
    return seqan::atEnd(*this->_reader);
}

// ----------------------------------------------------------------------------
// Member Function SamReader_::readHeader()
// ----------------------------------------------------------------------------

inline int SamReader_::readHeader(BamHeader & header, BamIOContext<StringSet<CharString> > & context)
{
    return seqan::readRecord(header, context, *this->_reader, Sam());
}

// ----------------------------------------------------------------------------
// Member Function SamReader_::readRecord()
// ----------------------------------------------------------------------------

inline int SamReader_::readRecord(BamAlignmentRecord & record, BamIOContext<StringSet<CharString> > & context)
{
    return seqan::readRecord(record, context, *this->_reader, Sam());
}

// ----------------------------------------------------------------------------
// Member Function SamReader_::readRecord()
// ----------------------------------------------------------------------------

inline int SamReader_::close()
{
    if (this->_stream == &this->_fstream)
        this->_fstream.close();
    return 0;
}

// ----------------------------------------------------------------------------
// Member Function SamReader_::fileSize()
// ----------------------------------------------------------------------------

inline __int64 SamReader_::fileSize() const
{
    return this->_fileSize;
}

// ----------------------------------------------------------------------------
// Member Function SamReader_::positionInFile()
// ----------------------------------------------------------------------------

// TODO(holtgrew): This does not work yet.

inline __int64 SamReader_::positionInFile() const
{
    return 0;
}

}  // namespace seqan

#endif  // #ifndef CORE_INCLUDE_SEQAN_BAM_IO_SAM_READER_H_
