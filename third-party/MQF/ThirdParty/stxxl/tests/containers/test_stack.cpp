/***************************************************************************
 *  tests/containers/test_stack.cpp
 *
 *  Part of the STXXL. See http://stxxl.sourceforge.net
 *
 *  Copyright (C) 2003 Roman Dementiev <dementiev@mpi-sb.mpg.de>
 *  Copyright (C) 2010 Andreas Beckmann <beckmann@cs.uni-frankfurt.de>
 *
 *  Distributed under the Boost Software License, Version 1.0.
 *  (See accompanying file LICENSE_1_0.txt or copy at
 *  http://www.boost.org/LICENSE_1_0.txt)
 **************************************************************************/

//! \example containers/test_stack.cpp
//! This is an example of how to use \c stxxl::STACK_GENERATOR class
//! to generate an \b external stack type
//! with \c stxxl::grow_shrink_stack implementation, \b four blocks per page,
//! block size \b 4096 bytes

#include <stxxl/stack>

// forced instantiation
template class stxxl::STACK_GENERATOR<size_t, stxxl::external, stxxl::normal, 4, 4096>;
template class stxxl::STACK_GENERATOR<size_t, stxxl::migrating, stxxl::normal, 4, 4096>;
template class stxxl::STACK_GENERATOR<size_t, stxxl::external, stxxl::grow_shrink, 4, 4096>;
template class stxxl::STACK_GENERATOR<size_t, stxxl::external, stxxl::grow_shrink2, 1, 4096>;

template <typename stack_type>
void test_lvalue_correctness(stack_type& stack, int a, int b)
{
    int i;
    STXXL_CHECK(stack.empty());
    for (i = 0; i < a; ++i)
        stack.push(i);
    for (i = 0; i < b; ++i)
        stack.push(i);
    for (i = 0; i < b; ++i)
        stack.pop();
    stack.top() = 0xbeeff00d;
    for (i = 0; i < b; ++i)
        stack.push(i);
    for (i = 0; i < b; ++i)
        stack.pop();
    if ((stack.top() != (size_t)(0xbeeff00d))) {
        STXXL_ERRMSG("STACK MISMATCH AFTER top() LVALUE MODIFICATION (0x" << std::hex << stack.top() << " != 0xbeeff00d)");
        STXXL_CHECK(stack.top() == (size_t)(0xbeeff00d));
    }
    for (i = 0; i < a; ++i)
        stack.pop();
}

template <typename stack_type>
void simple_test(stack_type& my_stack, size_t test_size)
{
    for (size_t i = 0; i < test_size; i++)
    {
        my_stack.push(i);
        STXXL_CHECK(my_stack.top() == i);
        STXXL_CHECK(my_stack.size() == i + 1);
    }

    for (size_t i = test_size; i > 0; )
    {
        --i;
        STXXL_CHECK(my_stack.top() == i);
        my_stack.pop();
        STXXL_CHECK(my_stack.size() == i);
    }

    for (size_t i = 0; i < test_size; i++)
    {
        my_stack.push(i);
        STXXL_CHECK(my_stack.top() == i);
        STXXL_CHECK(my_stack.size() == i + 1);
    }

    // test swap
    stack_type s2;
    std::swap(s2, my_stack);
    std::swap(s2, my_stack);

    for (size_t i = test_size; i > 0; )
    {
        --i;
        STXXL_CHECK(my_stack.top() == i);
        my_stack.pop();
        STXXL_CHECK(my_stack.size() == i);
    }

    std::stack<size_t> int_stack;

    for (size_t i = 0; i < test_size; i++)
    {
        int_stack.push(i);
        STXXL_CHECK(int_stack.top() == i);
        STXXL_CHECK(int_stack.size() == i + 1);
    }

    stack_type my_stack1(int_stack);

    for (size_t i = test_size; i > 0; )
    {
        --i;
        STXXL_CHECK(my_stack1.top() == i);
        my_stack1.pop();
        STXXL_CHECK(my_stack1.size() == i);
    }

    STXXL_MSG("Test 1 passed.");

    test_lvalue_correctness(my_stack, 4 * 4096 / 4 * 2, 4 * 4096 / 4 * 2 * 20);
}

int main(int argc, char* argv[])
{
    typedef stxxl::STACK_GENERATOR<size_t, stxxl::external, stxxl::normal, 4, 4096>::result ext_normal_stack_type;
    typedef stxxl::STACK_GENERATOR<size_t, stxxl::migrating, stxxl::normal, 4, 4096>::result ext_migrating_stack_type;
    typedef stxxl::STACK_GENERATOR<size_t, stxxl::external, stxxl::grow_shrink, 4, 4096>::result ext_stack_type;
    typedef stxxl::STACK_GENERATOR<size_t, stxxl::external, stxxl::grow_shrink2, 1, 4096>::result ext_stack_type2;

    if (argc < 2)
    {
        STXXL_MSG("Usage: " << argv[0] << " test_size_in_pages");
        return -1;
    }
    {
        ext_normal_stack_type my_stack;
        simple_test(my_stack, atoi(argv[1]) * 4 * 4096 / sizeof(int));
    }
    {
        ext_migrating_stack_type my_stack;
        //simple_test(my_stack, atoi(argv[1]) * 4 * 4096 / sizeof(int));
    }
    {
        ext_stack_type my_stack;
        simple_test(my_stack, atoi(argv[1]) * 4 * 4096 / sizeof(int));
    }
    {
        // prefetch/write pool with 10 blocks prefetching and 10 block write cache (> D is recommended)
        stxxl::read_write_pool<ext_stack_type2::block_type> pool(10, 10);
        // create a stack that does not prefetch (level of prefetch aggressiveness 0)
        ext_stack_type2 my_stack(pool, 0);
        size_t test_size = atoi(argv[1]) * 4 * 4096 / sizeof(int);

        for (size_t i = 0; i < test_size; i++)
        {
            my_stack.push(i);
            STXXL_CHECK(my_stack.top() == i);
            STXXL_CHECK(my_stack.size() == i + 1);
        }
        my_stack.set_prefetch_aggr(10);
        for (size_t i = test_size; i > 0; )
        {
            --i;
            STXXL_CHECK(my_stack.top() == i);
            my_stack.pop();
            STXXL_CHECK(my_stack.size() == i);
        }

        for (size_t i = 0; i < test_size; i++)
        {
            my_stack.push(i);
            STXXL_CHECK(my_stack.top() == i);
            STXXL_CHECK(my_stack.size() == i + 1);
        }

        // test swap
        ext_stack_type2 s2(pool, 0);
        std::swap(s2, my_stack);
        std::swap(s2, my_stack);

        for (size_t i = test_size; i > 0; )
        {
            --i;
            STXXL_CHECK(my_stack.top() == i);
            my_stack.pop();
            STXXL_CHECK(my_stack.size() == i);
        }

        STXXL_MSG("Test 2 passed.");

        test_lvalue_correctness(my_stack, 4 * 4096 / 4 * 2, 4 * 4096 / 4 * 2 * 20);
    }

    return 0;
}

// vim: et:ts=4:sw=4
