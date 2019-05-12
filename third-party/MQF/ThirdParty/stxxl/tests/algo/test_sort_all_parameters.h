/***************************************************************************
 *  tests/algo/test_sort_all_parameters.h
 *
 *  Part of the STXXL. See http://stxxl.sourceforge.net
 *
 *  Copyright (C) 2002 Roman Dementiev <dementiev@mpi-sb.mpg.de>
 *  Copyright (C) 2008, 2009 Andreas Beckmann <beckmann@cs.uni-frankfurt.de>
 *
 *  Distributed under the Boost Software License, Version 1.0.
 *  (See accompanying file LICENSE_1_0.txt or copy at
 *  http://www.boost.org/LICENSE_1_0.txt)
 **************************************************************************/

#include <limits>
#include <stxxl/bits/config.h>
#include <stxxl/types>

template <unsigned n>
struct bulk
{
    char m_data[n];

    bulk()
    {
#if STXXL_WITH_VALGRIND
        memset(m_data, 0, n);
#endif
    }
};

template <>
struct bulk<0>
{ };

template <typename KEY, unsigned SIZE>
struct my_type
{
    typedef KEY key_type;

    key_type m_key;
    bulk<SIZE - sizeof(key_type)> m_data;

    my_type() { }
    my_type(key_type k) : m_key(k) { }

#ifdef KEY_COMPARE
    key_type key() const
    {
        return m_key;
    }
#endif

    static my_type<KEY, SIZE> min_value()
    {
        return my_type<KEY, SIZE>(std::numeric_limits<key_type>::min());
    }
    static my_type<KEY, SIZE> max_value()
    {
        return my_type<KEY, SIZE>(std::numeric_limits<key_type>::max());
    }
};

template <typename KEY, unsigned SIZE>
std::ostream& operator << (std::ostream& o, const my_type<KEY, SIZE> obj)
{
    o << obj.m_key;
    return o;
}

#ifndef KEY_COMPARE

template <typename KEY, unsigned SIZE>
bool operator < (const my_type<KEY, SIZE>& a, const my_type<KEY, SIZE>& b)
{
    return a.m_key < b.m_key;
}

template <typename KEY, unsigned SIZE>
bool operator == (const my_type<KEY, SIZE>& a, const my_type<KEY, SIZE>& b)
{
    return a.m_key == b.m_key;
}

template <typename KEY, unsigned SIZE>
bool operator != (const my_type<KEY, SIZE>& a, const my_type<KEY, SIZE>& b)
{
    return a.m_key != b.m_key;
}

template <typename T>
struct Cmp : public std::less<T>
{
    bool operator () (const T& a, const T& b) const
    {
        return a.m_key < b.m_key;
    }

    static T min_value()
    {
        return T::min_value();
    }
    static T max_value()
    {
        return T::max_value();
    }
};

#else

template <typename KEY, unsigned SIZE>
bool operator < (const my_type<KEY, SIZE>& a, const my_type<KEY, SIZE>& b)
{
    return a.key() < b.key();
}

#endif

struct zero
{
    unsigned operator () ()
    {
        return 0;
    }
};
