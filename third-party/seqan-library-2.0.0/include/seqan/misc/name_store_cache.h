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
// Author: David Weese <david.weese@fu-berlin.de>
// ==========================================================================

//TODO(weese): Make the name stores thread-safe, i.e. provide an atomic
//             getIdByName() that returns the id of an existing name or append
//             the new one with its id


#ifndef SEQAN_HEADER_MISC_NAME_STORE_CACHE_H
#define SEQAN_HEADER_MISC_NAME_STORE_CACHE_H

namespace seqan {

// ============================================================================
// Forwards
// ============================================================================

// ============================================================================
// Tags, Classes, Enums
// ============================================================================

// ----------------------------------------------------------------------------
// struct NameStoreLess_
// ----------------------------------------------------------------------------

template <typename TNameStore, typename TName>
struct NameStoreLess_
{
    typedef typename Position<TNameStore>::Type TId;

    TNameStore *nameStore;
    TName *name;

    NameStoreLess_() {}

    NameStoreLess_(TNameStore &_nameStore, TName &_name):
        nameStore(&_nameStore),
        name(&_name) {}

    template <typename TId>
    inline bool operator() (TId a, TId b) const
    {
        if (a != maxValue(a))
        {
            if (b != maxValue(b))
                return (*nameStore)[a] < (*nameStore)[b];
            else
                return (*nameStore)[a] < *name;
        } else
        {
            if (b != maxValue(b))
                return *name < (*nameStore)[b];
            else
                return false;
        }
    }
};

// ----------------------------------------------------------------------------
// class NameStoreCache
// ----------------------------------------------------------------------------

/*!
 * @class NameStoreCache
 * @headerfile <seqan/misc/name_store_cache.h>
 * @brief Fast mapping from string names to numeric ids.
 *
 * NameStore objects store a binary search tree (using <tt>std::map&lt;&gt;</tt>) on a @link StringSet @endlink (the
 * name store).  They store a pointer to this name store but not the name store itself.
 *
 * When adding values to the name store using @link NameStoreCache#appendName @endlink, the cache (i.e, the binary
 * search tree) is automatically updated.  When modifying the name store when not using @link NameStoreCache#appendName
 * @endlink, you have to use @link NameStoreCache#refresh @endlink to update the cache after modifying and before
 * querying.
 *
 * The fast lookups can be performed using @link NameStoreCache#getIdByName @endlink and @link NameStoreCache#nameToId
 * @endlink.  The query function @link NameStoreCache#nameToId @endlink, the cache can also be modified (and thus
 * updated).
 *
 * @signature template <typename TNameStore[, typename TName]>
 *            class NameStoreCache;
 *
 * @tparam TNameStore The type to use for the name store.  Usually a @link StringSet @endlink of
 *                    @link CharString @endlink.
 * @tparam TName      The type to use for the names, defaults to <tt>Value&lt;TNameStore&gt;::Type</tt>.
 *
 * @section Example
 *
 * The demo below shows how to initialize a NameStoreCache with an existing name store, lookup existing names, add new
 * names, and add names during lookup.
 *
 * @include demos/misc/name_store_cache.cpp
 *
 * Here is the output:
 *
 * @include demos/misc/name_store_cache.cpp.stdout
 */

/*!
 * @fn NameStoreCache::NameStoreCache
 * @brief Constructors.
 *
 * NameStore cache offers the default constructor, copy constructor, and construction using an existing name store.
 *
 * @signature NameStoreCache::NameStoreCache();
 * @signature NameStoreCache::NameStoreCache(other);
 * @signature NameStoreCache::NameStoreCache(nameStore);
 *
 * @param[in] other     The other NameStoreCache to copy from.
 * @param[in] nameStore A NameStore for which a pointer is stored.
 */

template <typename TNameStore, typename TName = String<typename Value<typename Value<TNameStore>::Type>::Type> >
class NameStoreCache
{
public:
    typedef typename Position<TNameStore>::Type TId;
    typedef NameStoreLess_<TNameStore, TName> TLess;
    typedef std::set<TId, TLess> TSet;

    TSet nameSet;
    // TODO(holtgrew): Mutable here necessary for conceptual const-ness.  However, we would rather have a thread-safe interface!
    TName mutable name;

    NameStoreCache()
    {}

    NameStoreCache(TNameStore & nameStore):
        nameSet(TLess(nameStore, name))
    {
        for (unsigned i = 0; i < length(nameStore); ++i)
            nameSet.insert(i);
    }

    NameStoreCache(NameStoreCache const & other):
        nameSet(TLess(host(other), name))
    {
        for (unsigned i = 0; i < length(host(other)); ++i)
            nameSet.insert(i);
    }
};

// ============================================================================
// Metafunctions
// ============================================================================

// ============================================================================
// Functions
// ============================================================================

// ----------------------------------------------------------------------------
// refresh()
// ----------------------------------------------------------------------------

template <typename TNameStore, typename TName>
inline TNameStore &
host(NameStoreCache<TNameStore, TName> & cache)
{
    return *cache.nameSet.key_comp().nameStore;
}

template <typename TNameStore, typename TName>
inline TNameStore &
host(NameStoreCache<TNameStore, TName> const & cache)
{
    return *cache.nameSet.key_comp().nameStore;
}

// ----------------------------------------------------------------------------
// Function clear()
// ----------------------------------------------------------------------------

/*!
 * @fn NameStoreCache#clear
 * @brief Reset the NameStoreCache (not the name store).
 *
 * @signature void clear(cache);
 *
 * @param[in,out] cache The NameStoreCache to clear.
 */

template <typename TNameStore, typename TName>
inline void
clear(NameStoreCache<TNameStore, TName> &cache)
{
    cache.nameSet.clear();
}

// ----------------------------------------------------------------------------
// Function empty()
// ----------------------------------------------------------------------------

/*!
 * @fn NameStoreCache#empty
 * @brief Query whether there are any entries in the cache (not the name store).
 *
 * @signature bool empty(cache);
 *
 * @param[in,out] cache The NameStoreCache to clear.
 *
 * @return bool <tt>true</tt> if the NameStoreCache is empty.
 */

template <typename TNameStore, typename TName>
inline bool
empty(NameStoreCache<TNameStore, TName> const &cache)
{
    return cache.nameSet.empty();
}

// ----------------------------------------------------------------------------
// Function refresh()
// ----------------------------------------------------------------------------

/*!
 * @fn NameStoreCache#refresh
 * @brief Rebuild the name store cache.
 *
 * Use this after modifying the underlying NameStore before querying.
 *
 * @signature void refresh(cache);
 *
 * @param[in,out] nameStore The NameStoreCache to rebuild.
 */

template <typename TNameStore, typename TName>
inline void
refresh(NameStoreCache<TNameStore, TName> &cache)
{
    clear(cache);
    for (unsigned i = 0; i < length(*cache.nameSet.key_comp().nameStore); ++i)
        cache.nameSet.insert(i);
}

// ----------------------------------------------------------------------------
// Function appendName()
// ----------------------------------------------------------------------------

/*!
 * @fn NameStoreCache#appendName
 * @brief Append a name to a name store and register it in the cache.
 *
 * The NameStoreCache only registers update to the name store when performed by this function.
 *
 * @signature void appendName(cache, name);
 *
 * @param[in,out] cache     The NameStoreCache to use for faster access.
 * @param[in]     name      The name to append to the store (@link ContainerConcept#Value @endlink of
 *                          <tt>TNameStore</tt>).
 */

template <typename TNameStore, typename TName>
void appendName(TNameStore & nameStore, TName const & name)
{
    appendValue(nameStore, name, Generous());
}

template <typename TCNameStore, typename TCName, typename TName>
void appendName(NameStoreCache<TCNameStore, TCName> & cache, TName const & name)
{
    appendValue(host(cache), name, Generous());
    cache.nameSet.insert(length(host(cache)) - 1);
}

// TODO(holtgrew): Add deprecation annotation for compiler warnings.

// deprecated.
// In the future we want to use only one argument either nameStore or nameStoreCache (has a reference to the nameStore)
template <typename TNameStore, typename TName, typename TContext>
void appendName(TNameStore &nameStore, TName const & name, TContext &)
{
    appendName(nameStore, name);
}

// deprecated.
template <typename TNameStore, typename TName, typename TCNameStore, typename TCName>
void appendName(TNameStore &nameStore, TName const & name, NameStoreCache<TCNameStore, TCName> &context)
{
    appendValue(nameStore, name, Generous());
    context.nameSet.insert(length(nameStore) - 1);
}

// ----------------------------------------------------------------------------
// Function getIdByName()
// ----------------------------------------------------------------------------

/*!
 * @fn NameStoreCache#getIdByName
 * @brief Get id/index of a string in a name store using a NameStoreCache.
 *
 * @signature bool getIdByName(idx, cache, name);
 *
 * @param[out]    idx       The variable to store the index in the store of (@link IntegerConcept @endlink).
 * @param[in]     cache     The NameStoreCache to use for speeding up the lookup.
 * @param[in]     name      The name to search in the name store (@link ContainerConcept#Value @endlink of
 *                          <tt>TNameStore</tt>).
 *
 * @return bool <tt>true</tt> if the name could be found and <tt>false</tt> otherwise.
 */

template <typename TPos, typename TNameStore, typename TName>
bool getIdByName(TPos & pos, TNameStore const & nameStore, TName const & name)
{
    typedef typename Iterator<TNameStore const, Standard>::Type TNameStoreIter;

    // Iterator over read names
    for (TNameStoreIter iter = begin(nameStore, Standard()); iter != end(nameStore, Standard()); ++iter)
    {
        // if the element was found
        if (name == getValue(iter))
        {
            // set the ID
            pos = iter - begin(nameStore, Standard());
            // and break the loop
            return true;
        }
    }
    return false;
}

template <typename TCNameStore, typename TCName, typename TName, typename TPos>
inline bool
getIdByName(TPos & pos, NameStoreCache<TCNameStore, TCName> const & context, TName const & name)
{
    typedef typename Position<TCNameStore const>::Type TId;
    typedef NameStoreCache<TCNameStore, TCName> const TNameStoreCache;
    typedef typename TNameStoreCache::TSet TSet;

    TSet const &set = context.nameSet;
    typename TSet::const_iterator it;

    // (weese:)
    // To avoid local variables to copy the name into we use a member in context.
    // However, changing the NameStoreCache per query is not thread-safe and the user might not notice it.
    // To avoid pitfalls, we should introduce a critical section.
    //SEQAN_OMP_PRAGMA(critical (nameStoreFind))
    {

        context.name = name;
        it = set.find(maxValue<TId>());
    }

    if (it != set.end())
    {
        pos = *it;
        return true;
    }
    return false;
}

// deprecated.
template <typename TNameStore, typename TName, typename TPos, typename TContext>
inline bool
getIdByName(TNameStore const & nameStore, TName const & name, TPos & pos, TContext const & /*not a cache*/)
{
    return getIdByName(pos, nameStore, name);
}

// deprecated.
template<typename TNameStore, typename TName, typename TPos, typename TCNameStore, typename TCName>
inline bool
getIdByName(TNameStore const & /*nameStore*/, TName const & name, TPos & pos, NameStoreCache<TCNameStore, TCName> const & context)
{
    return getIdByName(pos, context, name);
}

// ----------------------------------------------------------------------------
// Function nametoId()
// ----------------------------------------------------------------------------

/*!
 * @fn NameStoreCache#nameToId
 * @brief Translate a name to a numeric id, adding the name to the store and cache if new.
 *
 * @signature TPos nameToId(cache, name);
 *
 * @param[in,out] cache The NameStoreCache use for translating the name to a numeric id.
 * @param[in]     name  The name to add (@link ContainerConcept#Value @endlink of <tt>TNameStore</tt>).
 *
 * @return TPos The numeric id of the name in the store and cache (@link ContainerConcept#Position Position @endlink of
 *              <tt>TNameStore</tt>).
 *
 * @note Since <tt>cache</tt> is modified if <tt>name</tt> is not known in cache, it is a <b>non-const</b> parameter
 *       for this function.
 *
 * If <tt>name</tt> is in <tt>cache</tt> then its numeric position/index/id in the name store is returned.  If it is not
 * in the name store then it is appended to the name store and registered with the NameStoreCache.
 */

// Append contig name to name store, if not known already.
template <typename TNameStore, typename TName, typename TName2>
typename Position<TNameStore>::Type
nameToId(NameStoreCache<TNameStore, TName> & cache, TName2 const & name)
{
    typename Size<TNameStore>::Type nameId = 0;
    if (!getIdByName(nameId, cache, name))
    {
        nameId = length(host(cache));
        appendName(cache, name);
    }
    return nameId;
}

}  // namespace seqan

#endif  // #ifndef SEQAN_HEADER_MISC_NAME_STORE_CACHE_H
