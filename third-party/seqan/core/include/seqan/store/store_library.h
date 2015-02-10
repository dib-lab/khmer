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

#ifndef SEQAN_HEADER_STORE_LIBRARY_H
#define SEQAN_HEADER_STORE_LIBRARY_H

namespace SEQAN_NAMESPACE_MAIN
{

//////////////////////////////////////////////////////////////////////////////
// Library Store
//////////////////////////////////////////////////////////////////////////////

/**
.Class.LibraryStoreElement
..summary:Represents a fragment library.
..cat:Fragment Store
..signature:LibraryStoreElement<>
..signature:LibraryStoreElement<TMean[, TStd[, TSpec]]>
..param.TMean:The type to represent the library size mean.
...default:$double$
..param.TStd:The type to represent the library size standard deviation.
...default:$double$
..param.TSpec:The specialization type.
...default:$void$
..remarks:A fragment library is a set of mate-pairs having a certain distribution of insert sizes.
..remarks:Value type of the @Memvar.FragmentStore#libraryStore@ string.

.Memfunc.LibraryStoreElement#LibraryStoreElement
..summary:Constructor
..signature:LibraryStoreElement<>()
..signature:LibraryStoreElement<TMean[, TStd[, TSpec]]> ()
..remarks:Sets $mean$ and $std$ to $0$.
..class:Class.LibraryStoreElement
.Memvar.LibraryStoreElement#mean
..summary:The library size mean.
..class:Class.LibraryStoreElement
.Memvar.LibraryStoreElement#std
..summary:The library size standard deviation.
..class:Class.LibraryStoreElement
..include:seqan/store.h
*/

template <typename TMean = double, typename TStd = double, typename TSpec = void>
struct LibraryStoreElement
{
	TMean		mean;		// mean library size in bps
	TStd		std;	// library size variance

	LibraryStoreElement() : mean(0), std(0) {}

    inline bool operator==(LibraryStoreElement const & other) const
    {
        return mean == other.mean && std == other.std;
    }
};

//////////////////////////////////////////////////////////////////////////////


}// namespace SEQAN_NAMESPACE_MAIN

#endif //#ifndef SEQAN_HEADER_...
