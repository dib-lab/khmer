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
// Example program used in the documentation of the String class that shows
// the basic usage of this class.
// ==========================================================================

#include <iostream>

#include <seqan/sequence.h>
#include <seqan/stream.h>  // for output

using namespace seqan;

int main()
{
    // Examples for constructing strings.
//![initializing strings]
    String<char> strA;        // default construction
    String<char> strB(strA);  // copy construction
    String<char> strC("copy from other sequence");
//![initializing strings]

    // Common string operations.
//![usual operations]
    // Assignment of sequence with the same alphabet and of another string.
    strA = "Hello World!";
    strB = strA;
    std::cout << strA << "\n"   // => "Hello World!"
              << strB << "\n";  // => "Hello World!"

    // Appending of values (characters) and whole strings.
    appendValue(strA, ' ');
    append(strA, strB);
    std::cout << strA << "\n";  // => "Hello World! Hello World!"

    // Element-wise access and replacing.
    std::cout << strB[3] << "\n";  // => "l";
    strB[3] = 'g';
    std::cout << strB[3] << "\n";  // => "g"

    replace(strB, 5, 12, "land");
    std::cout << strB << "\n";  // => "Helgoland"

    // Removal of elements and strings.
    erase(strA, 5, 18);
    erase(strA, length(strA) - 1);
    std::cout << strA << "\n";  // => "Hello World"
//![usual operations]

    // The swap function is implemented for the String class for efficient
    // swapping of string contents.
//![swap function]
    swap(strA, strB);
    std::cout << strA << "\n"   // => "Helgoland"
              << strB << "\n";  // => "Hello World"
//![swap function]

    // A demonstration of memory allocation related function.
//![clear and resize]
    std::cout << "length(strA) = " << length(strA) << "\n"       // "length(strA) == 9"
              << "capacity(strA) = " << capacity(strA) << "\n";  // "capacity(strA) == 32"
    clear(strA);
    std::cout << "length(strA) = " << length(strA) << "\n"       // "length(strA) == 0"
              << "capacity(strA) = " << capacity(strA) << "\n";  // "capacity(strA) == 32"
    shrinkToFit(strA);
    std::cout << "length(strA) = " << length(strA) << "\n"       // "length(strA) == 0"
              << "capacity(strA) = " << capacity(strA) << "\n";  // "capacity(strA) == 0"
//![clear and resize]

    return 0;
}
