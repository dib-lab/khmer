/***************************************************************************
 *  tests/common/test_log2.cpp
 *
 *  Part of the STXXL. See http://stxxl.sourceforge.net
 *
 *  Copyright (C) 2008 Andreas Beckmann <beckmann@cs.uni-frankfurt.de>
 *  Copyright (C) 2013 Timo Bingmann <tb@panthema.net>
 *
 *  Distributed under the Boost Software License, Version 1.0.
 *  (See accompanying file LICENSE_1_0.txt or copy at
 *  http://www.boost.org/LICENSE_1_0.txt)
 **************************************************************************/

#include <iostream>
#include <cmath>
#include <cassert>
#include <stxxl/bits/config.h>
#include <stxxl/bits/common/tmeta.h>
#include <stxxl/bits/common/utils.h>
#include <stxxl/bits/verbose.h>

using stxxl::LOG2;
using stxxl::unsigned_type;

template <unsigned_type i>
void log_i(unsigned_type floorv, unsigned_type ceilv)
{
    std::cout << i << "\t" << (i < 1000000 ? "\t" : "")
              << stxxl::LOG2_floor<i>::value << "\t"
              << stxxl::LOG2<i>::floor << "\t"
              << stxxl::LOG2<i>::ceil << std::endl;

    std::cout << "\t\t"
              << stxxl::ilog2_floor(i) << "\t"
              << stxxl::ilog2_floor(i) << "\t"
              << stxxl::ilog2_ceil(i) << std::endl;

    STXXL_CHECK(stxxl::LOG2_floor<i>::value == stxxl::ilog2_floor(i));
    STXXL_CHECK(stxxl::LOG2<i>::floor == stxxl::ilog2_floor(i));
    STXXL_CHECK(stxxl::LOG2<i>::ceil == stxxl::ilog2_ceil(i));

    if (i <= 1)
    {
        STXXL_CHECK(stxxl::LOG2_floor<i>::value == 0);
        STXXL_CHECK(stxxl::LOG2<i>::floor == 0);
        STXXL_CHECK(stxxl::LOG2<i>::ceil == 0);
    }
    else if (i == 2)
    {
        STXXL_CHECK(stxxl::LOG2_floor<i>::value == 1);
        STXXL_CHECK(stxxl::LOG2<i>::floor == 1);
        STXXL_CHECK(stxxl::LOG2<i>::ceil == 1);
    }
    else
    {
        STXXL_CHECK(stxxl::LOG2_floor<i>::value == floorv);
        STXXL_CHECK(stxxl::LOG2<i>::floor == floorv);
        STXXL_CHECK(stxxl::LOG2<i>::ceil == ceilv);

#if 0                                        // not many compiler have log2l()
        if (i <= ((stxxl::uint64)(1) << 59)) // does not work for higher powers
        {
            STXXL_CHECK(stxxl::LOG2_floor<i>::value == (unsigned_type)floorl(log2l(i)));
            STXXL_CHECK(stxxl::LOG2<i>::floor == (unsigned_type)floorl(log2l(i)));
            STXXL_CHECK(stxxl::LOG2<i>::ceil == (unsigned_type)ceill(log2l(i)));
        }
#endif
    }

    std::cout << "\n";
}

template <unsigned_type i>
void log_ipm1(unsigned_type p)
{
    log_i<i - 1>(p - 1, p);
    log_i<i>(p, p);
    log_i<i + 1>(p, p + 1);
    std::cout << std::endl;
}

int main()
{
    std::cout << "i\t\tLOG2<i>\tLOG2<i>\tLOG2<i>" << std::endl;
    std::cout << "\t\t\tfloor\tceil" << std::endl;
    log_ipm1<1 << 0>(0);
    log_ipm1<1 << 1>(1);
    log_ipm1<1 << 2>(2);
    log_ipm1<1 << 3>(3);
    log_ipm1<1 << 4>(4);
    log_ipm1<1 << 5>(5);
    log_ipm1<1 << 6>(6);
    log_ipm1<1 << 7>(7);
    log_ipm1<1 << 8>(8);
    log_ipm1<1 << 9>(9);
    log_ipm1<1 << 10>(10);
    log_ipm1<1 << 11>(11);
    log_ipm1<1 << 12>(12);
    log_ipm1<1 << 16>(16);
    log_ipm1<1 << 24>(24);
    log_ipm1<1 << 30>(30);
    log_ipm1<1UL << 31>(31);
#if __WORDSIZE == 64
    log_ipm1<1UL << 32>(32);
    log_ipm1<1UL << 33>(33);
    log_ipm1<1UL << 48>(48);
    log_ipm1<1UL << 50>(50);
    log_ipm1<1UL << 55>(55);
    log_ipm1<1UL << 63>(63);
#endif
}
