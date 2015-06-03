// ==========================================================================
//                 SeqAn - The Library for Sequence Analysis
// ==========================================================================
// Copyright (c) 2006-2015, Knut Reinert, FU Berlin
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
// Example program used in the documentation of the Segment class.
// ==========================================================================

#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/stream.h>  // for I/O

using namespace seqan;

int main()
{
//![basic operations]
    typedef Prefix<String<char> >::Type TPrefix;
    typedef Infix<String<char> >::Type  TInfix;
    typedef Suffix<String<char> >::Type TSuffix;

    String<char> text = "This is a text!";

    TPrefix preA(text, 4);
    TInfix infA(text, 10, 14);
    TSuffix sufA(text, 10);
    std::cout << preA << " " << infA << " " << sufA << "\n";  // => "This text text!"

    String<char> str;
    append(str, preA);
    append(str, infA);
    append(str, sufA);
    std::cout << str << "\n";  // => "This text text!"

    std::cout << preA[0] << " " << infA[0] << " " << sufA[0] << "\n";  // => "T t t"

    preA[0] = 'X';
    infA[0] = 'X';
    sufA[1] = 'X';
    std::cout << text << "\n";  // => "Xhis is a XXxt!"

    typedef Iterator<TInfix, Standard>::Type TIter;
    TIter it = begin(preA, Standard());
    it += 2;
    *it = 'Y';
    std::cout << text << "\n";  // => "XhYs is a XXxt!"
//![basic operations]

//![metafunction examples]
    typedef Infix<TInfix>::Type  TInfix2;  // == TInfix
    typedef Prefix<TInfix>::Type TInfix3;  // == TInfix
    typedef Suffix<TInfix>::Type TInfix4;  // == TInfix

    typedef Infix<TPrefix>::Type  TInfix5;   // == TInfix
    typedef Prefix<TPrefix>::Type TPrefix2;  // == TPrefix
    typedef Suffix<TPrefix>::Type TInfix6;   // == TInfix

    typedef Infix<TSuffix>::Type  TInfix7;   // == TInfix
    typedef Prefix<TSuffix>::Type TInfix8;   // == TPrefix
    typedef Suffix<TSuffix>::Type TSuffix2;  // == TSuffix
//![metafunction examples]

//![explicit segment]
    typedef Segment<TSuffix, PrefixSegment> TExplicitPrefix;
    TExplicitPrefix preB(sufA, 3);
    std::cout << preB << "\n";  // => "XXx"
//![explicit segment]

    return 0;
}
