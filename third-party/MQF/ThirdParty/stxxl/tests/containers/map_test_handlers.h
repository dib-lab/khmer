/***************************************************************************
 *  tests/containers/map_test_handlers.h
 *
 *  Part of the STXXL. See http://stxxl.sourceforge.net
 *
 *  Copyright (C) 2004, 2005 Thomas Nowak <t.nowak@imail.de>
 *  Copyright (C) 2006 Roman Dementiev <dementiev@ira.uka.de>
 *  Copyright (C) 2010 Andreas Beckmann <beckmann@cs.uni-frankfurt.de>
 *
 *  Distributed under the Boost Software License, Version 1.0.
 *  (See accompanying file LICENSE_1_0.txt or copy at
 *  http://www.boost.org/LICENSE_1_0.txt)
 **************************************************************************/

//! \file map_test_handlers.h
//! This file contains help functions for testing of stxxl::map.

#ifndef STXXL_TEST_HEADER__MAP_TEST_HANDLERS_H_
#define STXXL_TEST_HEADER__MAP_TEST_HANDLERS_H_

#include <stxxl/bits/namespace.h>

STXXL_BEGIN_NAMESPACE

// ***********************************************
// THERE
// ***********************************************

template <typename MAPTYPE>
bool there(const MAPTYPE& map_, const typename MAPTYPE::key_type& key, const typename MAPTYPE::mapped_type& data)
{
    typename MAPTYPE::const_iterator iter = map_.find(key);
    if (!(iter->second == data))
    {
        STXXL_VERBOSE2("iter=(" << (*iter).first << ":" << (*iter).second << ")");
        STXXL_VERBOSE2("key=" << key);
        STXXL_VERBOSE2("data=" << data);
        return false;
    }
    return true;
}

// ***********************************************
// IS EQUAL END
// ***********************************************

template <typename MAPTYPE>
bool is_equal_end(const MAPTYPE& map_, typename MAPTYPE::const_iterator& iter)
{
    return iter == map_.end();
}

// ***********************************************
// IS SAME
// ***********************************************

template <typename value_type>
bool is_same(value_type& v1, value_type& v2)
{
    return v1.first == v2.first && v1.second == v2.second;
}

template <typename value_type>
bool is_same(const value_type& v1, const value_type& v2)
{
    return v1.first == v2.first && v1.second == v2.second;
}

// ***********************************************
// NOT THERE
// ***********************************************

template <typename MAPTYPE>
bool not_there(const MAPTYPE& map_, const typename MAPTYPE::key_type& key)
{
    return map_.find(key) == map_.end();
}

// ***********************************************
// IS EMPTY
// ***********************************************

template <typename MAPTYPE>
bool is_empty(const MAPTYPE& map_)
{
    return map_.empty();
}

// ***********************************************
// IS END
// ***********************************************

template <typename MAPTYPE>
bool is_end(MAPTYPE& map_, typename MAPTYPE::iterator& iter)
{
    return iter == map_.end();
}

template <typename MAPTYPE>
bool is_end(const MAPTYPE& map_, typename MAPTYPE::const_iterator& iter)
{
    return iter == map_.end();
}

// ***********************************************
// IS SIZE
// ***********************************************

template <typename MAPTYPE>
bool is_size(const MAPTYPE& map_, const typename MAPTYPE::size_type size)
{
    return map_.size() == size;
}

STXXL_END_NAMESPACE

#endif // !STXXL_TEST_HEADER__MAP_TEST_HANDLERS_H_
