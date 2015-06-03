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
// Demo for the Switch metaprogramming statement.
// ==========================================================================

#include <iostream>
#include <seqan/basic.h>

using namespace seqan;

//![switch demo]
int switchTest(Nothing const &)
{
    return -1;
}

int switchTest(False const &)
{
    return 0;
}

int switchTest(True const &)
{
    return 1;
}

int switchTest(NilCase const &)
{
    return 2;
}

template <int X>
struct SwitchTest
{
    typedef typename Switch<
        X,
        Case<-1, Nothing,
             Case<0, False,
                  Case<1, True
                       > > > >::Type Type;
};

int main()
{
    typedef SwitchTest<-1>::Type T1;
    typedef SwitchTest<0>::Type T2;
    typedef SwitchTest<1>::Type T3;
    typedef SwitchTest<2>::Type T4;

    std::cout << switchTest(T1()) << "\n"   // => "-1"
              << switchTest(T2()) << "\n"   // => "0"
              << switchTest(T3()) << "\n"   // => "1"
              << switchTest(T4()) << "\n";  // => "2"

    return 0;
}

//![switch demo]
