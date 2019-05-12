/***************************************************************************
 *  tests/common/test_tuple.cpp
 *
 *  Part of the STXXL. See http://stxxl.sourceforge.net
 *
 *  Copyright (C) 2013 Timo Bingmann <tb@panthema.net>
 *
 *  Distributed under the Boost Software License, Version 1.0.
 *  (See accompanying file LICENSE_1_0.txt or copy at
 *  http://www.boost.org/LICENSE_1_0.txt)
 **************************************************************************/

#include <stxxl/bits/common/tuple.h>
#include <stxxl/bits/common/uint_types.h>
#include <stxxl/bits/verbose.h>

#include <sstream>

using stxxl::uint40;
using stxxl::uint64;

// force instantiation of some tuple types
template struct stxxl::tuple<int>;
template struct stxxl::tuple<int, int>;
template struct stxxl::tuple<int, int, int>;
template struct stxxl::tuple<int, int, int, int>;
template struct stxxl::tuple<int, int, int, int, int>;
template struct stxxl::tuple<int, int, int, int, int, int>;

// force instantiation of more tuple types
template struct stxxl::tuple<uint64>;
template struct stxxl::tuple<uint64, uint64>;
template struct stxxl::tuple<uint64, uint64, uint64>;
template struct stxxl::tuple<uint64, uint64, uint64, uint64>;
template struct stxxl::tuple<uint64, uint64, uint64, uint64, uint64>;
template struct stxxl::tuple<uint64, uint64, uint64, uint64, uint64, uint64>;

// force instantiation of even more tuple types
template struct stxxl::tuple<uint40>;
template struct stxxl::tuple<uint40, uint40>;
template struct stxxl::tuple<uint40, uint40, uint40>;
template struct stxxl::tuple<uint40, uint40, uint40, uint40>;
template struct stxxl::tuple<uint40, uint40, uint40, uint40, uint40>;
template struct stxxl::tuple<uint40, uint40, uint40, uint40, uint40, uint40>;

void test1()
{
    stxxl::tuple<int, int> pair;
    pair.first = 5;
    pair.second = 42;

    std::ostringstream oss1;
    oss1 << pair;
    STXXL_CHECK(oss1.str() == "(5,42)");

    pair = stxxl::tuple<int, int>(1, 2);

    std::ostringstream oss2;
    oss2 << pair;
    STXXL_CHECK(oss2.str() == "(1,2)");
}

int main(int, char**)
{
    test1();
    return 0;
}
