/***************************************************************************
 *  tests/common/test_counting_ptr.cpp
 *
 *  Small test case for reference counting in stxxl::counting_ptr.
 *
 *  Part of the STXXL. See http://stxxl.sourceforge.net
 *
 *  Copyright (C) 2013 Timo Bingmann <tb@panthema.net>
 *
 *  Distributed under the Boost Software License, Version 1.0.
 *  (See accompanying file LICENSE_1_0.txt or copy at
 *  http://www.boost.org/LICENSE_1_0.txt)
 **************************************************************************/

#include <stxxl/bits/common/counting_ptr.h>
#include <stxxl/bits/verbose.h>

static unsigned int count_deletes = 0;

// derive from counted_object to include reference counter
struct my_int : public stxxl::counted_object
{
    int i;

    my_int(int _i) : i(_i) { }

    // count number of destructor calls
    ~my_int()
    { ++count_deletes; }
};

typedef stxxl::counting_ptr<my_int> int_ptr;
typedef stxxl::const_counting_ptr<my_int> int_cptr;

int_cptr run_test()
{
    // create object and pointer to it
    int_ptr i1 = new my_int(42);

    STXXL_CHECK(i1->i == 42);
    STXXL_CHECK((*i1).i == 42);
    STXXL_CHECK(i1.get()->i == 42);
    STXXL_CHECK(i1->unique());

    // make pointer sharing same object
    int_ptr i2 = i1;

    STXXL_CHECK(i2->i == 42);
    STXXL_CHECK(!i1->unique());
    STXXL_CHECK(i1 == i2);
    STXXL_CHECK(i1->get_reference_count() == 2);

    // make another pointer sharing the same object
    int_ptr i3 = i2;

    STXXL_CHECK(i3->i == 42);
    STXXL_CHECK(i3->get_reference_count() == 3);

    // replace object in i3 with new integer
    i3 = new my_int(5);
    STXXL_CHECK(i1 != i3);
    STXXL_CHECK(i1->get_reference_count() == 2);

    // create a const pointer
    int_cptr i4 = i3;

    // quitting the function will release the first object
    return i4;
}

int main()
{
    {
        // run tests and return a reference counted object
        int_cptr t1 = run_test();

        // check number of objects destructed
        STXXL_CHECK(count_deletes == 1);

        // quitting the block will release the ptr
    }

    STXXL_CHECK(count_deletes == 2);

    return 0;
}
