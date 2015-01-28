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
// Author: David Weese <david.weese@fu-berlin.de>
// ==========================================================================

#ifndef SEQAN_HEADER_STORE_MATEPAIR_H
#define SEQAN_HEADER_STORE_MATEPAIR_H

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////
// Mate Store
//////////////////////////////////////////////////////////////////////////////

/**
.Class.MatePairStoreElement
..summary:Represents a mate-pair.
..cat:Fragment Store
..signature:MatePairStoreElement<>
..signature:MatePairStoreElement<TSpec>
..param.TSpec:The specialization type.
...default:$void$
..remarks:A mate-pair consists of two reads sequenced from opposite ends and strands of the same fragment.
The insert size of a mate-pair is the size of the fragment.
..remarks:Value type of the @Memvar.FragmentStore#matePairStore@ string.

.Memfunc.MatePairStoreElement#MatePairStoreElement
..summary:Constructor
..signature:MatePairStoreElement<> ()
..signature:MatePairStoreElement<TSpec> ()
..remarks:Sets $readId[0]$, $readId[1]$ and $libId$ to $INVALID_ID$.
..class:Class.MatePairStoreElement
.Memvar.MatePairStoreElement#readId[2]
..summary:Refers to two paired reads in the @Memvar.FragmentStore#readStore@ or contains $INVALID_ID$ values.
..type:Metafunction.Id
..class:Class.MatePairStoreElement
.Memvar.MatePairStoreElement#libId
..summary:Refers to a library in the @Memvar.FragmentStore#libraryStore@ or is $INVALID_ID$ if the mate-pair has no library.
..type:Metafunction.Id
..class:Class.MatePairStoreElement
.Memvar.MatePairStoreElement#INVALID_ID
..summary:Constant to represent an invalid id.
..type:Metafunction.Id
..class:Class.MatePairStoreElement
..include:seqan/store.h
*/

template <typename TSpec = void>
struct MatePairStoreElement
{
	typedef typename Id<MatePairStoreElement>::Type TId;

	static const TId INVALID_ID;
	
	TId		readId[2];	// refers to the two reads of a mate-pair, INVALID_ID if this is a singleton fragment (e.g. in afg: reads refer to fragments (mate pairs) and these refer to libraries, singletons refer to an empty fragment)
	TId		libId;

	MatePairStoreElement() : libId(INVALID_ID) 
	{
		readId[0] = INVALID_ID;
		readId[1] = INVALID_ID;
	}

    inline bool operator==(MatePairStoreElement const & other) const
    {
        return readId[0] == other.readId[0] &&
                readId[1] == other.readId[1] &&
                libId == other.libId;
    }
};

//////////////////////////////////////////////////////////////////////////////

template <typename TSpec>
const typename Id<MatePairStoreElement<TSpec> >::Type
MatePairStoreElement<TSpec>::INVALID_ID = MaxValue<typename Id<MatePairStoreElement<TSpec> >::Type>::VALUE;

//////////////////////////////////////////////////////////////////////////////

}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
