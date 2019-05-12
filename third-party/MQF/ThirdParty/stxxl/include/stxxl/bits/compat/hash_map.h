/***************************************************************************
 *  include/stxxl/bits/compat/hash_map.h
 *
 *  compatibility interface to C++ standard extension hash_map
 *
 *  Part of the STXXL. See http://stxxl.sourceforge.net
 *
 *  Copyright (C) 2008, 2010, 2011 Andreas Beckmann <beckmann@cs.uni-frankfurt.de>
 *  Copyright (C) 2009, 2010 Johannes Singler <singler@kit.edu>
 *  Copyright (C) 2014 Timo Bingmann <tb@panthema.net>
 *
 *  Distributed under the Boost Software License, Version 1.0.
 *  (See accompanying file LICENSE_1_0.txt or copy at
 *  http://www.boost.org/LICENSE_1_0.txt)
 **************************************************************************/

#ifndef STXXL_COMPAT_HASH_MAP_HEADER
#define STXXL_COMPAT_HASH_MAP_HEADER

#include <stxxl/bits/config.h>
#include <stxxl/bits/namespace.h>

#if __cplusplus >= 201103L || (STXXL_MSVC && _MSC_VER >= 1900)
 #include <unordered_map>
#elif STXXL_MSVC
 #include <hash_map>
#elif defined(__GNUG__) && ((__GNUC__ * 10000 + __GNUC_MINOR__ * 100) >= 40200) && \
    (!defined(__ICC) || (__ICC > 1110))
 #include <tr1/unordered_map>
#else
 #include <ext/hash_map>
#endif

STXXL_BEGIN_NAMESPACE

template <class KeyType>
struct compat_hash {
#if __cplusplus >= 201103L || (STXXL_MSVC && _MSC_VER >= 1900)
    typedef std::hash<KeyType> result;
#elif STXXL_MSVC
    typedef stdext::hash_compare<KeyType> result;
#elif defined(__GNUG__) && ((__GNUC__ * 10000 + __GNUC_MINOR__ * 100) >= 40200) && \
    (!defined(__ICC) || (__ICC > 1110))
    typedef std::tr1::hash<KeyType> result;
#else
    typedef __gnu_cxx::hash<KeyType> result;
#endif
};

template <class KeyType, class MappedType,
          class HashType = typename compat_hash<KeyType>::result>
struct compat_hash_map {
#if __cplusplus >= 201103L || (STXXL_MSVC && _MSC_VER >= 1900)
    typedef std::unordered_map<KeyType, MappedType, HashType> result;
#elif STXXL_MSVC
    typedef stdext::hash_map<KeyType, MappedType, HashType> result;
#elif defined(__GNUG__) && ((__GNUC__ * 10000 + __GNUC_MINOR__ * 100) >= 40200) && \
    (!defined(__ICC) || (__ICC > 1110))
    typedef std::tr1::unordered_map<KeyType, MappedType, HashType> result;
#else
    typedef __gnu_cxx::hash_map<KeyType, MappedType, HashType> result;
#endif
};

template <class KeyType, class MappedType,
          class HashType = typename compat_hash<KeyType>::result>
struct compat_hash_multimap {
#if __cplusplus >= 201103L || (STXXL_MSVC && _MSC_VER >= 1900)
    typedef std::unordered_multimap<KeyType, MappedType, HashType> result;
#elif STXXL_MSVC
    typedef stdext::hash_multimap<KeyType, MappedType, HashType> result;
#elif defined(__GNUG__) && ((__GNUC__ * 10000 + __GNUC_MINOR__ * 100) >= 40200) && \
    (!defined(__ICC) || (__ICC > 1110))
    typedef std::tr1::unordered_multimap<KeyType, MappedType, HashType> result;
#else
    typedef __gnu_cxx::hash_multimap<KeyType, MappedType, HashType> result;
#endif
};

STXXL_END_NAMESPACE

#endif // !STXXL_COMPAT_HASH_MAP_HEADER
