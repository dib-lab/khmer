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
// This header provides the readFasta() function to be used in demos.
// ==========================================================================

#ifndef CORE_INCLUDE_SEQAN_SEQ_IO_SIMPLE_READ_FASTA_H_
#define CORE_INCLUDE_SEQAN_SEQ_IO_SIMPLE_READ_FASTA_H_

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function readFasta()
// ----------------------------------------------------------------------------

/**
.Function.readFasta
..cat:Input/Output
..signature:int readFasta(seq, filename)
..summary:Read first sequence from a FASTA file.
..description:
This function is meant for demo programs and allows to read the first sequence from a FASTA file very easily.
..param.seq:@Class.String@ to store the result.
...type:Class.String
..param.filename:Path to the file to read from.
...type:nolink:$char const *$
..returns:$int$ with status code, $0$ for success, $1$ for error.
..example.text:Read the file given by the first parameter to the program and print it to the screen.
..example.code:
#include <iostream>
#include <seqan/sequence.h>
#include <seqan/seq_io.h>

int main(int argc, char const ** argv)
{
    if (argc != 2)
        return 1;

    seqan::Dna5String seq;
    if (readFasta(seq, argv[1]) != 0)
        return 1;

    std::cout << seq << '\n';
    return 0;
}
..include:seqan/seq_io.h
*/

template <typename TValue, typename TSpec>
int readFasta(String<TValue, TSpec> & seq, char const * filename)
{
    SequenceStream seqStream(filename);
    if (!isGood(seqStream))
        return 1;
    CharString id;
    return readRecord(id, seq, seqStream);
}

}  // namespace seqan

#endif  // #ifndef CORE_INCLUDE_SEQAN_SEQ_IO_SIMPLE_READ_FASTA_H_
