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
// Example program used in the documentation of the enable-if metaprogramming
// utilities.
// ==========================================================================

#include <string>

#include <seqan/basic.h>

using namespace seqan;

#if !defined(_MSC_VER)  // Currently, there are some issues with MSVC and concepts.
class EnableIfExample
{
public:
    int num;

//![enable if example constructor]
    template <typename T>
    EnableIfExample(T const & n, SEQAN_CTOR_ENABLE_IF(Is<IntegerConcept<T> >)) :
        num(0)
    {
        ignoreUnusedVariableWarning(dummy);
    }

//![enable if example constructor]

//![disable if example constructor]
    template <typename T>
    EnableIfExample(T const & n, SEQAN_CTOR_DISABLE_IF(Is<IntegerConcept<T> >)) :
        num(0)
    {
        ignoreUnusedVariableWarning(dummy);
    }

//![disable if example constructor]

//![enable if example function]
    template <typename T>
    SEQAN_FUNC_ENABLE_IF(Is<IntegerConcept<T> >)
    f(T /* x */)
    { /* ... */ }
//![enable if example function]

//![disable if example function]
    template <typename T>
    SEQAN_FUNC_DISABLE_IF(Is<IntegerConcept<T> >)
    f(T /* x */)
    { /* ... */}
//![disable if example function]
};
#endif  // #if !defined(_MSC_VER)

int main()
{
#if !defined(_MSC_VER)
    EnableIfExample ex1(1);
    (void)ex1;
    EnableIfExample ex2("asdf");
    (void)ex2;
#endif  // #if !defined(_MSC_VER)

    return 0;
}
