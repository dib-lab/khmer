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

#ifndef CORE_INCLUDE_SEQAN_BAM_IO_BAM_WRITER_H_
#define CORE_INCLUDE_SEQAN_BAM_IO_BAM_WRITER_H_

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// TODO(holtgrew): Allow writing BAM to stdout? Extend Stream<Bgzf>?

class BamWriter_ :
    public XamWriter_
{
public:
    // The BGZF stream to write to.
    Stream<Bgzf> _stream;
    // Flag indicating whether there was an error or not.
    // TODO(holtgrew): Could we also use streamError()?
    bool _isGood;

    BamWriter_() :
        XamWriter_()
    {}

    BamWriter_(CharString const & filename);

    // XamWriter_ interface.

    virtual int open(CharString const & filename);
    virtual bool isGood();
    virtual int writeHeader(BamHeader const & header,
                            BamIOContext<StringSet<CharString> > const & context);
    virtual int writeRecord(BamAlignmentRecord const & record,
                            BamIOContext<StringSet<CharString> > const & context);
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
// Member Function BamWriter_::BamWriter_()
// ----------------------------------------------------------------------------

inline BamWriter_::BamWriter_(CharString const & filename) :
    XamWriter_(filename)
{
    this->open(filename);
}

// ----------------------------------------------------------------------------
// Member Function BamWriter_::open()
// ----------------------------------------------------------------------------

inline int BamWriter_::open(CharString const & filename)
{
    if (!seqan::open(this->_stream, toCString(filename), "w"))
    {
        _isGood = false;
        return 1;
    }

    return 0;
}

// ----------------------------------------------------------------------------
// Member Function BamWriter_::isGood()
// ----------------------------------------------------------------------------

inline bool BamWriter_::isGood()
{
    return this->_isGood;
}

// ----------------------------------------------------------------------------
// Member Function BamWriter_::writeHeader()
// ----------------------------------------------------------------------------

inline int BamWriter_::writeHeader(BamHeader const & header,
                            BamIOContext<StringSet<CharString> > const & context)
{
    return write2(this->_stream, header, context, Bam());
}

// ----------------------------------------------------------------------------
// Member Function BamWriter_::writeRecord()
// ----------------------------------------------------------------------------

inline int BamWriter_::writeRecord(BamAlignmentRecord const & record,
                                   BamIOContext<StringSet<CharString> > const & context)
{
    return write2(this->_stream, record, context, Bam());
}

// ----------------------------------------------------------------------------
// Member Function BamWriter_::flush()
// ----------------------------------------------------------------------------

inline int BamWriter_::flush()
{
    return streamFlush(this->_stream);
}

// ----------------------------------------------------------------------------
// Member Function BamWriter_::close()
// ----------------------------------------------------------------------------

inline int BamWriter_::close()
{
    seqan::close(this->_stream);
    return 0;
}

}  // namespace seqan

#endif  // #ifndef CORE_INCLUDE_SEQAN_BAM_IO_BAM_WRITER_H_
