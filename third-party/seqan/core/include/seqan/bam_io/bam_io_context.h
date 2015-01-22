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
// Class BamIOContext, accessor functions.
// ==========================================================================

#ifndef CORE_INCLUDE_SEQAN_BAM_IO_BAM_IO_CONTEXT_H_
#define CORE_INCLUDE_SEQAN_BAM_IO_BAM_IO_CONTEXT_H_

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

/**
.Class.BamIOContext
..cat:BAM I/O
..signature:BamIOContext<TNameStore[, TNameStoreCache]>
..summary:The I/O context to use for BAM I/O.
..param.TNameStore:The name store class.
..param.TNameStoreCache:The name store cache class.
...default:@Class.NameStoreCache@<TNameStore>
..include:bam_io.h
..example.text:Creating a @Class.BamIOContext@ for a raw @Class.StringSet@ of @Shortcut.CharString@.
..example.code:
StringSet<CharString> nameStore;
NameStoreCache<StringSet<CharString> > nameStoreCache(nameStore);
BamIOContext<StringSet<CharString> > bamIOContext(nameStore, nameStoreCache);
// ...
..example.text:Using a @Class.BamIOContext@ with a @Class.FragmentStore@.
..example.code:
typedef FragmentStore<>::TContigNameStore         TNameStore;
typedef NameStoreCache<TNameStore>                TNameStoreCache;
FragmentStore<> store;
// Optionally, do something with store.
typedef BamIOContext<TNameStore, TNameStoreCache> TBamIOContext;
TBamIOContext bamIOContext(store.contigNameStore, store.contigNameStoreCache);
// ...

.Memfunc.BamIOContext#BamIOContext
..class:Class.BamIOContext
..signature:BamIOContext()
..summary:Constructor.
..remarks:Only the default constructor is provided.

.Typedef.BamIOContext#TNameStore
..class:Class.BamIOContext
..summary:The name store class.

.Typedef.BamIOContext#TNameStoreCache
..class:Class.BamIOContext
..summary:The name store cache class.
*/

template <typename TNameStore_, typename TNameStoreCache_ = NameStoreCache<TNameStore_> >
class BamIOContext
{
public:
    typedef TNameStore_ TNameStore;
    typedef TNameStoreCache_ TNameStoreCache;

    TNameStore * _nameStore;
    TNameStoreCache * _nameStoreCache;
    CharString buffer;

    BamIOContext() : _nameStore(0), _nameStoreCache(0)
    {}

    BamIOContext(TNameStore & nameStore, TNameStoreCache & nameStoreCache) :
            _nameStore(&nameStore), _nameStoreCache(&nameStoreCache)
    {}
};

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function nameStore()
// ----------------------------------------------------------------------------

/**
.Function.BamIOContext#nameStore
..class:Class.BamIOContext
..cat:BAM I/O
..summary:Return reference to name store from @Class.BamIOContext@.
..signature:nameStore(context)
..param.context:The @Class.BamIOContext@ to query.
...type:Class.BamIOContext
..see:Typedef.BamIOContext#TNameStore
..include:seqan/bam_io.h
*/

// TODO(holtgrew): Rename to referenceNameStore
template <typename TNameStore, typename TNameStoreCache>
TNameStore &
nameStore(BamIOContext<TNameStore, TNameStoreCache> & context)
{
    SEQAN_ASSERT(context._nameStore != 0);
    return *context._nameStore;
}

template <typename TNameStore, typename TNameStoreCache>
TNameStore const &
nameStore(BamIOContext<TNameStore, TNameStoreCache> const & context)
{
    SEQAN_ASSERT(context._nameStore != 0);
    return *context._nameStore;
}

// ----------------------------------------------------------------------------
// Function nameStoreCache()
// ----------------------------------------------------------------------------

/**
.Function.BamIOContext#nameStoreCache
..class:Class.BamIOContext
..cat:BAM I/O
..summary:Return reference to name store cache from @Class.BamIOContext@.
..signature:nameStoreCache(context)
..param.context:The @Class.BamIOContext@ to query.
...type:Class.BamIOContext
..see:Typedef.BamIOContext#TNameStoreCache
..include:seqan/bam_io.h
..see:Function.BamIOContext#nameStore
*/

// TODO(holtgrew): Rename to referenceNameStoreCache
template <typename TNameStore, typename TNameStoreCache>
TNameStoreCache &
nameStoreCache(BamIOContext<TNameStore, TNameStoreCache> & context)
{
    SEQAN_ASSERT(context._nameStoreCache != 0);
    return *context._nameStoreCache;
}

template <typename TNameStore, typename TNameStoreCache>
TNameStoreCache const &
nameStoreCache(BamIOContext<TNameStore, TNameStoreCache> const & context)
{
    SEQAN_ASSERT(context._nameStoreCache != 0);
    return *context._nameStoreCache;
}

}  // namespace seqan

#endif  // #ifndef CORE_INCLUDE_SEQAN_BAM_IO_BAM_IO_CONTEXT_H_
