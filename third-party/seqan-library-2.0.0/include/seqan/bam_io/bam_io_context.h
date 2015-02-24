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
// Class BamIOContext, accessor functions.
// ==========================================================================

#ifndef INCLUDE_SEQAN_BAM_IO_BAM_IO_CONTEXT_H_
#define INCLUDE_SEQAN_BAM_IO_BAM_IO_CONTEXT_H_

namespace seqan {

// ============================================================================
// Tags
// ============================================================================

// ----------------------------------------------------------------------------
// Tag NameStoreMember
// ----------------------------------------------------------------------------

/*!
 * @defgroup BamIOContextMemberTag BamIOContext Member Tags
 * @brief Defines standard tags used to get the type of the members of the @link BamIOContext @endlink using the @link Member @endlink metafunction.
 */

/*!
 * @tag BamIOContextMemberTag#NameStoreMember
 * @brief Tag used to get the type for the <tt>NameStore</tt>.
 * @headerfile <seqan/stream.h>
 *
 * @signature typedef Tag<NameStoreMember_> NameStoreMember;
 */

struct NameStoreMember_;
typedef Tag<NameStoreMember_> NameStoreMember;

// ----------------------------------------------------------------------------
// Tag NameStoreCacheMember
// ----------------------------------------------------------------------------

/*!
 * @tag BamIOContextMemberTag#NameStoreCacheMember
 * @brief Tag used to get the type for the <tt>NameStoreCache</tt>.
 * @headerfile <seqan/stream.h>
 *
 * @signature typedef Tag<NameStoreCacheMember_> NameStoreCacheMember;
 */

struct NameStoreCacheMember_;
typedef Tag<NameStoreCacheMember_> NameStoreCacheMember;

// ----------------------------------------------------------------------------
// Tag LengthStoreMember
// ----------------------------------------------------------------------------

/*!
 * @tag BamIOContextMemberTag#LengthStoreMember
 * @brief Tag used to get the type for the <tt>LengthStore</tt>.
 * @headerfile <seqan/stream.h>
 *
 * @signature typedef Tag<NameStoreMember_> LengthStoreMember;
 */

struct LengthStoreMember_;
typedef Tag<LengthStoreMember_> LengthStoreMember;

// ============================================================================
// Classes
// ============================================================================

/*!
 * @class BamIOContext
 * @headerfile <seqan/bam_io.h>
 * @brief The I/O context to use for BAM I/O.
 *
 * @signature template <typename TNameStore[, typename TNameStoreCache]>
 *            class BamIOContext;
 *
 * @tparam TNameStore      The type used to represent the names.
 * @tparam TNameStoreCache The type used to cache the names. Defaults to @link NameStoreCache @endlink &lt;TNameStore&gtl;.
 *
 * BamIOContext objects store the names of (and provide a cache for) reference contig names.
 *
 * @section Examples
 *
 * Creating a @link BamIOContext @endlink for a raw @link StringSet @endlink of @link CharString @endlink.
 *
 * @code{.cpp}
 * StringSet<CharString> contigNames;
 * NameStoreCache<StringSet<CharString> > contigNamesCache(contigNames);
 * BamIOContext<StringSet<CharString> > bamIOContext(contigNames, contigNamesCache);
 * // ...
 * @endcode
 *
 * Using a @link BamIOContext @endlink with a @link FragmentStore @endlink.
 *
 * @code{.cpp}
 * typedef FragmentStore<>::TContigNameStore         TNameStore;
 * typedef NameStoreCache<TNameStore>                TNameStoreCache;
 * FragmentStore<> store;
 * // Optionally, do something with store.
 * typedef BamIOContext<TNameStore, TNameStoreCache> TBamIOContext;
 * TBamIOContext bamIOContext(store.contigNameStore, store.contigNameStoreCache);
 * // ...
 * @endcode
 */

/*!
 * @fn BamIOContext::BamIOContext
 * @headerfile <seqan/bam_io.h>
 * @brief Constructor.
 *
 * @signature BamIOContext::BamIOContext();
 * @signature BamIOContext::BamIOContext(contigNameStore, contigNamesStoreCache);
 *
 * Default constructor or construction with references to sequence and sample names.
 */

template <typename TNameStore_        = StringSet<CharString>,
          typename TNameStoreCache_   = NameStoreCache<TNameStore_>,
          typename TStorageSpec       = void>
class BamIOContext
{
public:
    typedef typename Member<BamIOContext, NameStoreMember>::Type            TNameStore;
    typedef typename Member<BamIOContext, NameStoreCacheMember>::Type       TNameStoreCache;
    typedef typename Member<BamIOContext, LengthStoreMember>::Type          TLengthStore;

    typedef typename If<IsSameType<TStorageSpec, void>,
                        Dependent<>, TStorageSpec>::Type                    TNSStorageSpec;
    typedef typename If<IsSameType<TStorageSpec, void>,
                        Owner<>, TStorageSpec>::Type                        TSLStorageSpec;

    typedef typename StorageSwitch<TNameStore, TNSStorageSpec>::Type        TNameStoreMember;
    typedef typename StorageSwitch<TNameStoreCache, TNSStorageSpec>::Type   TNameStoreCacheMember;
    typedef typename StorageSwitch<TLengthStore, TSLStorageSpec>::Type      TLengthStoreMember;

    TNameStoreMember        _contigNames;
    TNameStoreCacheMember   _contigNamesCache;
    TLengthStoreMember      _contigLengths;
    CharString              buffer;
    String<CharString>      buffers;
    String<unsigned>        translateFile2GlobalRefId;

    BamIOContext() :
        _contigNames(TNameStoreMember()),
        _contigNamesCache(ifSwitch(typename IsPointer<TNameStoreCacheMember>::Type(),
                                 (TNameStoreCache*)NULL,
                                 _contigNames)),
        _contigLengths(TLengthStoreMember())
    {}

    BamIOContext(TNameStore & contigNames_, TNameStoreCache & contigNamesCache_) :
        _contigNames(_referenceCast<typename Parameter_<TNameStoreMember>::Type>(contigNames_)),
        _contigNamesCache(ifSwitch(typename IsPointer<TNameStoreCacheMember>::Type(),
                                 &contigNamesCache_,
                                 _contigNames)),
        _contigLengths(TLengthStoreMember())
    {}

    template <typename TOtherStorageSpec>
    BamIOContext(BamIOContext<TNameStore, TNameStoreCache, TOtherStorageSpec> & other) :
        _contigNames(_referenceCast<typename Parameter_<TNameStoreMember>::Type>(contigNames(other))),
        _contigNamesCache(ifSwitch(typename IsPointer<TNameStoreCacheMember>::Type(),
                                 &contigNamesCache(other),
                                 _contigNames)),
        _contigLengths(_referenceCast<typename Parameter_<TLengthStoreMember>::Type>(contigLengths(other)))
    {}
};

// ============================================================================
// Metafunctions
// ============================================================================

// ----------------------------------------------------------------------------
// Metafunction NameStore
// ----------------------------------------------------------------------------

template <typename TNameStore, typename TNameStoreCache, typename TStorageSpec>
struct Member<BamIOContext<TNameStore, TNameStoreCache, TStorageSpec>,
              NameStoreMember>
{
    typedef TNameStore Type;
};

// ----------------------------------------------------------------------------
// Metafunction NameStoreCache
// ----------------------------------------------------------------------------

template <typename TNameStore, typename TNameStoreCache, typename TStorageSpec>
struct Member<BamIOContext<TNameStore, TNameStoreCache, TStorageSpec>,
              NameStoreCacheMember>
{
    typedef TNameStoreCache Type;
};

// ----------------------------------------------------------------------------
// Metafunction LengthStore
// ----------------------------------------------------------------------------

template <typename TNameStore, typename TNameStoreCache, typename TStorageSpec>
struct Member<BamIOContext<TNameStore, TNameStoreCache, TStorageSpec>,
              LengthStoreMember>
{
    typedef String<__int32> Type;
};

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// Function contigNames()
// ----------------------------------------------------------------------------

/*!
 * @fn BamIOContext#contigNames
 * @brief Return reference to contig names from @link BamIOContext @endlink.
 *
 * @signature TNameStoreRef contigNames(context);
 *
 * @param[in] context The @link BamIOContext @endlink to query.
 *
 * @return TNameStoreRef A reference to the <tt>TNameStore</tt> of the context.
 */

template <typename TNameStore, typename TNameStoreCache, typename TStorageSpec>
inline TNameStore &
contigNames(BamIOContext<TNameStore, TNameStoreCache, TStorageSpec> & context)
{
    return _referenceCast<TNameStore &>(context._contigNames);
}

template <typename TNameStore, typename TNameStoreCache, typename TStorageSpec>
inline TNameStore const &
contigNames(BamIOContext<TNameStore, TNameStoreCache, TStorageSpec> const & context)
{
    return _referenceCast<TNameStore const &>(context._contigNames);
}

template <typename TNameStore, typename TNameStoreCache>
inline void
setContigNames(BamIOContext<TNameStore, TNameStoreCache, Dependent<> > & context, TNameStore & contigNames)
{
    context._contigNames = &contigNames;
}

// ----------------------------------------------------------------------------
// Function contigLengths()
// ----------------------------------------------------------------------------

/*!
 * @fn BamIOContext#contigLengths
 * @brief Return reference to contig lengths from @link BamIOContext @endlink.
 *
 * @signature TLengthStoreRef contigLengths(context);
 *
 * @param[in] context The @link BamIOContext @endlink to query.
 *
 * @return TLengthStoreRef A reference to the <tt>TLengthStore</tt> of the context.
 */

template <typename TNameStore, typename TNameStoreCache, typename TStorageSpec>
inline typename BamIOContext<TNameStore, TNameStoreCache, TStorageSpec>::TLengthStore &
contigLengths(BamIOContext<TNameStore, TNameStoreCache, TStorageSpec> & context)
{
    typedef typename BamIOContext<TNameStore, TNameStoreCache, TStorageSpec>::TLengthStore TLengthStore;
    return _referenceCast<TLengthStore &>(context._contigLengths);
}

template <typename TNameStore, typename TNameStoreCache, typename TStorageSpec>
inline typename BamIOContext<TNameStore, TNameStoreCache, TStorageSpec>::TLengthStore const &
contigLengths(BamIOContext<TNameStore, TNameStoreCache, TStorageSpec> const & context)
{
    typedef typename BamIOContext<TNameStore, TNameStoreCache, TStorageSpec>::TLengthStore TLengthStore;
    return _referenceCast<TLengthStore const &>(context._contigLengths);
}

template <typename TNameStore, typename TNameStoreCache, typename TLengthStore>
inline void
setContigLengths(BamIOContext<TNameStore, TNameStoreCache, Dependent<> > & context, TLengthStore & contigLengths)
{
    context._contigLengths = &contigLengths;
}

// ----------------------------------------------------------------------------
// Function contigNamesCache()
// ----------------------------------------------------------------------------

/*!
 * @fn BamIOContext#contigNamesCache
 * @brief Return reference to contig names cache from @link BamIOContext @endlink.
 *
 * @signature TNameStoreCacheRef contigNamesCache(context);
 *
 * @param[in] context The @link BamIOContext @endlink to query.
 *
 * @return TNameStoreCacheRef A reference to the <tt>TNameStoreCache</tt> of the context.
 */

template <typename TNameStore, typename TNameStoreCache, typename TStorageSpec>
inline TNameStoreCache &
contigNamesCache(BamIOContext<TNameStore, TNameStoreCache, TStorageSpec> & context)
{
    return _referenceCast<TNameStoreCache &>(context._contigNamesCache);
}

template <typename TNameStore, typename TNameStoreCache, typename TStorageSpec>
inline TNameStoreCache const &
contigNamesCache(BamIOContext<TNameStore, TNameStoreCache, TStorageSpec> const & context)
{
    return _referenceCast<TNameStoreCache const &>(context._contigNamesCache);
}

template <typename TNameStore, typename TNameStoreCache>
inline void
setContigNamesCache(BamIOContext<TNameStore, TNameStoreCache, Dependent<> > & context, TNameStoreCache & contigNamesCache)
{
    context._contigNamesCache = &contigNamesCache;
}

}  // namespace seqan

#endif  // #ifndef INCLUDE_SEQAN_BAM_IO_BAM_IO_CONTEXT_H_
