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

#ifndef SEQAN_HEADER_PRIORITY_TYPE_BASE_H
#define SEQAN_HEADER_PRIORITY_TYPE_BASE_H

#include <functional>

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////

struct PriorityHeap;

//////////////////////////////////////////////////////////////////////////////


/**
.Class.PriorityType:
..cat:Miscellaneous
..summary:Stores items in such a way that the item with the highest priority is at the top.
..signature:PriorityType<TValue, TLess, TSpec>
..param.TValue:The value type that is stored.
...default:int
..param.TLess:The comparator type that is used for sorting the items stored.
...default:std::less<TValue>
..param.TSpec:The specializing type.
...default:@Spec.PriorityHeap@
..include:seqan/misc.h
*/
template <typename TValue = int, typename TLess = ::std::less<TValue>, typename TSpec = PriorityHeap>
class PriorityType;

//////////////////////////////////////////////////////////////////////////////

template <typename TValue, typename TLess, typename TSpec>
struct Value< PriorityType<TValue, TLess, TSpec> >
{
	typedef TValue Type;
};

template <typename TValue, typename TLess, typename TSpec>
struct Size< PriorityType<TValue, TLess, TSpec> >
{
	typedef typename Size<TValue>::Type Type;
};

//////////////////////////////////////////////////////////////////////////////

}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
