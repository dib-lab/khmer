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

#ifndef CORE_INCLUDE_SEQAN_BAM_IO_XAM_READER_H_
#define CORE_INCLUDE_SEQAN_BAM_IO_XAM_READER_H_

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// Base class implementing reading of SAM and BAM files.

class XamReader_
{
public:
    CharString _filename;

    XamReader_()
    {}

    XamReader_(CharString const & filename) :
        _filename(filename)
    {}

    virtual ~XamReader_()
    {}

    // Those of the following functions that return integers return status codes.  As usual, a status code of 0 means
    // "OK", one != 0 means that there was an error.

    // Open file with given name.
    virtual int open(CharString const & filename) = 0;

    // Return true if there is no error, false otherwise.
    virtual bool isGood() = 0;

    // Return true if read is at end.
    virtual bool atEnd() = 0;

    // Read the BAM header from the wrapped file.
    virtual int readHeader(BamHeader & header, BamIOContext<StringSet<CharString> > & context) = 0;

    // Read the next BAM record from the wrapped file.
    virtual int readRecord(BamAlignmentRecord & record, BamIOContext<StringSet<CharString> > & context) = 0;

    // Close the file again.
    virtual int close() = 0;

    // Return file size in bytes.
    virtual __int64 fileSize() const = 0;

    // Return approximate position in file.  For BAM, position of previous chunk is returned.
    virtual __int64 positionInFile() const = 0;
};

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

}  // namespace seqan

#endif  // #ifndef CORE_INCLUDE_SEQAN_BAM_IO_XAM_READER_H_
