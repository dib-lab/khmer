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

/**
.Class.NameStoreCache
..summary:Stores a mapping from names to ids.
..cat:Fragment Store
..signature:FragmentStore<>
..signature:NameStoreCache<TNameStore[, TName]>
..param.TNameStore:The name store to be cached.
...see:Class.FragmentStore
..param.TName:The name type.
...default:$Value<TNameStore>::Type$
...type:Shortcut.CharString

.Memfunc.NameStoreCache#NameStoreCache
..summary:Constructor
..signature:NameStoreCache<TNameStore, TName> (nameStore)
..param.nameStore:A name store, e.g. @Memvar.FragmentStore#readNameStore@
...see:Class.FragmentStore
..class:Class.NameStoreCache
..include:seqan/store.h
*/
	
template <typename TNameStore, typename TName = typename Value<TNameStore>::Type>
class NameStoreCache
{
public:
    typedef typename Position<TNameStore>::Type TId;
    typedef NameStoreLess_<TNameStore, TName> TLess;
    typedef std::set<TId, TLess> TSet;
    
    TSet nameSet;
    TNameStore *nameStore;
    // TODO(holtgrew): Mutable here necessary for conceptual const-ness.  However, we would rather have a thread-safe interface!
    TName mutable name;

    NameStoreCache(TNameStore &_nameStore):
        nameSet(TLess(_nameStore, name)),
        nameStore(&_nameStore)
    {
        for (unsigned i = 0; i < length(*nameStore); ++i)
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

/**
.Function.refresh:
..class:Class.NameStoreCache
..summary:Recreate a name store cache.
..cat:Fragment Store
..signature:refresh(cache)
..param.cache:A @Class.NameStoreCache@ object.
...type:Class.NameStoreCache
..see:Function.getIdByName
..include:seqan/store.h
*/
    
template <typename TNameStore, typename TName>
inline void
refresh(NameStoreCache<TNameStore, TName> &cache)
{
    cache.nameSet.clear();
    for (unsigned i = 0; i < length(*cache.nameStore); ++i)
        cache.nameSet.insert(i);
}


// ----------------------------------------------------------------------------
// getIdByName()
// ----------------------------------------------------------------------------

/**
.Function.getIdByName:
..summary:Get the id/index of a string in a name store with a cache.
..cat:Fragment Store
..signature:getIdByName(nameStore, name, id[, cache])
..param.nameStore:A name store, e.g. @Memvar.FragmentStore#readNameStore@
...see:Class.FragmentStore
..param.name:The name to be searched.
...type:Shortcut.CharString
..param.id:The resulting id.
..param.cache:A structure to efficiently retrieve the id for a given name. If ommited a brute force method is used to search.
...default:Tag.Nothing
...type:Class.NameStoreCache
..returns:$true$ if the name was found and $false$ if not.
..see:Function.getIdByName
..include:seqan/store.h
*/

template <typename TNameStore, typename TName, typename TPos>
inline bool 
getIdByName(TNameStore const & nameStore, TName const & name, TPos & pos)
{
    typedef typename Iterator<TNameStore const, Standard>::Type TNameStoreIter;
    
    // Iterator over read names
    for (TNameStoreIter iter = begin(nameStore, Standard()); iter != end(nameStore, Standard()); ++iter)
    {
        // if the element was found
        if (name == getValue(iter))
        {
            // set the ID
            pos = position(iter);
            // and break the loop
            return true;
        }
    }
    return false;
}

template <typename TNameStore, typename TName, typename TPos, typename TContext>
inline bool
getIdByName(TNameStore const & nameStore, TName const & name, TPos & pos, TContext const & /*not a cache*/)
{
    return getIdByName(nameStore, name, pos);
}

template<typename TNameStore, typename TName, typename TPos, typename TCNameStore, typename TCName>
inline bool
getIdByName(TNameStore const & /*nameStore*/, TName const & name, TPos & pos, NameStoreCache<TCNameStore, TCName> const & context)
{
    typedef typename Position<TNameStore const>::Type TId;
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

// ----------------------------------------------------------------------------
// appendName()
// ----------------------------------------------------------------------------

/**
.Function.appendName:
..summary:Appends a name to a name store.
..cat:Fragment Store
..signature:appendName(nameStore, name[, cache])
..param.nameStore:A name store, e.g. @Memvar.FragmentStore#readNameStore@
...see:Class.FragmentStore
..param.name:The name to be appended.
...type:Shortcut.CharString
..param.cache:A structure to efficiently retrieve the id for a given name. See @Function.getIdByName@.
...default:Tag.Nothing
...type:Class.NameStoreCache
..see:Function.getIdByName
..include:seqan/store.h
*/

template <typename TNameStore, typename TName>
inline void
appendName(TNameStore &nameStore, TName const & name)
{
    appendValue(nameStore, name, Generous());
}

template <typename TNameStore, typename TName, typename TContext>
inline void
appendName(TNameStore &nameStore, TName const & name, TContext &)
{
    appendName(nameStore, name);
}

template <typename TNameStore, typename TName, typename TCNameStore, typename TCName>
inline void
appendName(TNameStore &nameStore, TName const & name, NameStoreCache<TCNameStore, TCName> &context)
{
    appendValue(nameStore, name, Generous());
    context.nameSet.insert(length(nameStore) - 1);
}

}  // namespace seqan

#endif  // #ifndef SEQAN_HEADER_MISC_NAME_STORE_CACHE_H
