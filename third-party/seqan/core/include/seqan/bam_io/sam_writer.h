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

#ifndef CORE_INCLUDE_SEQAN_BAM_IO_SAM_WRITER_H_
#define CORE_INCLUDE_SEQAN_BAM_IO_SAM_WRITER_H_

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

class SamWriter_ :
    public XamWriter_
{
public:
    // Pointer to the output stream to write to.  Used such that we can also write SAM to stdout.
    std::ostream * _stream;
    // The file stream to write to if we do not write to stdout.  Pointed to by _stream.
    std::ofstream _fstream;

    SamWriter_() :
        XamWriter_(), _stream(0)
    {}

    SamWriter_(CharString const & filename);

    // XamWriter_ interface.

    virtual int open(CharString const & filename);
    virtual bool isGood();
    virtual int writeHeader(BamHeader const & header, BamIOContext<StringSet<CharString> > const & context);
    virtual int writeRecord(BamAlignmentRecord const & record, BamIOContext<StringSet<CharString> > const & context);
    virtual int flush();
    virtual int close();
};

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Member Function SamWriter_::SamWriter_()
// ----------------------------------------------------------------------------

inline SamWriter_::SamWriter_(CharString const & filename) :
    XamWriter_(filename), _stream(0)
{
    this->open(filename);
}

// ----------------------------------------------------------------------------
// Member Function SamWriter_::open()
// ----------------------------------------------------------------------------

inline int SamWriter_::open(CharString const & filename)
{
    if (filename == "-")
    {
        this->_stream = &std::cout;
    }
    else
    {
        this->_stream = &this->_fstream;
        this->_fstream.open(toCString(filename), std::ios::binary | std::ios::out);
        if (!this->_fstream.good())
            return 1;
    }

    return 0;
}

// ----------------------------------------------------------------------------
// Member Function SamWriter_::isGood()
// ----------------------------------------------------------------------------

inline bool SamWriter_::isGood()
{
    return this->_stream->good();
}

// ----------------------------------------------------------------------------
// Member Function SamWriter_::writeHeader()
// ----------------------------------------------------------------------------

inline int SamWriter_::writeHeader(BamHeader const & header, BamIOContext<StringSet<CharString> > const & context)
{
    return seqan::write2(*this->_stream, header, context, Sam());
}

// ----------------------------------------------------------------------------
// Member Function SamWriter_::writeRecord()
// ----------------------------------------------------------------------------

inline int SamWriter_::writeRecord(BamAlignmentRecord const & record, BamIOContext<StringSet<CharString> > const & context)
{
    return seqan::write2(*this->_stream, record, context, Sam());
}

// ----------------------------------------------------------------------------
// Member Function SamWriter_::flush()
// ----------------------------------------------------------------------------

inline int SamWriter_::flush()
{
    return streamFlush(*this->_stream);
}

// ----------------------------------------------------------------------------
// Member Function SamWriter_::close()
// ----------------------------------------------------------------------------

inline int SamWriter_::close()
{
    if (this->_stream == &this->_fstream)  // not to stdout
        this->_fstream.close();
    return 0;
}

}  // namespace seqan

#endif  // #ifndef CORE_INCLUDE_SEQAN_BAM_IO_SAM_WRITER_H_
