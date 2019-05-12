/***************************************************************************
 *  tests/common/test_swap_vector.cpp
 *
 *  Part of the STXXL. See http://stxxl.sourceforge.net
 *
 *  Copyright (C) 2014 Thomas Keh <thomas.keh@student.kit.edu>
 *
 *  Distributed under the Boost Software License, Version 1.0.
 *  (See accompanying file LICENSE_1_0.txt or copy at
 *  http://www.boost.org/LICENSE_1_0.txt)
 **************************************************************************/

#include <stxxl/bits/common/swap_vector.h>
#include <stxxl/bits/verbose.h>

class Test
{
    int m_i;

public:
    Test() : m_i(-1)
    {
        STXXL_MSG("Construct Test " << m_i);
    }
    Test(int i) : m_i(i)
    {
        STXXL_MSG("Construct Test " << m_i);
    }
    // Copy constructor
    Test(const Test& o)
    {
        STXXL_ERRMSG("Copy Test " << o.m_i << " into Test " << m_i);
        abort();
        m_i = o.m_i;
    }
    // Copy assignment
    Test& operator = (const Test& o)
    {
        STXXL_ERRMSG("Copy-assign Test " << o.m_i << " to Test " << m_i);
        abort();
        m_i = o.m_i;
        return *this;
    }
    ~Test()
    {
        STXXL_MSG("Destruct Test " << m_i);
    }
    //! swap vector with another one
    void swap(Test& obj)
    {
        STXXL_MSG("Swap Test " << m_i << " with Test " << obj.m_i);
        std::swap(m_i, obj.m_i);
    }
    int get_i() const
    {
        return m_i;
    }
};

//! Unary operator which returns true if m_i==-1 or m_i==8.
struct test_eraser
{
    bool operator () (const Test& t) const
    {
        return (t.get_i() == -1 || t.get_i() == 8);
    }
};

namespace std {

template <>
void swap(Test& a,
          Test& b)
{
    a.swap(b);
}

} // namespace std

int main()
{
    // Create a vector with 3 default Test objects (-1) and reserved space for 6 objects.
    // Content afterwards: {-1,-1,-1}
    stxxl::swap_vector<Test> vec(3, 6);

    // Push back 10 values from 0 to 9. Internally a swap resize will happen
    // Content afterwards: {-1,-1,-1,0,1,2,3,4,5,6,7,8,9}
    for (unsigned i = 0; i < 10; ++i) {
        Test a(i);
        vec.swap_back(a);
    }

    // Delete the third and the fourth object.
    // Content afterwards: {-1,-1,1,2,3,4,5,6,7,8,9}
    vec.erase(vec.begin() + 2, vec.begin() + 4);

    // Delete the sixth object.
    // Content afterwards: {-1,-1,1,2,3,5,6,7,8,9}
    vec.erase(vec.begin() + 5);

    // Check the values.
    int expected_vals[10] = { -1, -1, 1, 2, 3, 5, 6, 7, 8, 9 };
    STXXL_CHECK_EQUAL(vec.size(), 10);
    for (unsigned i = 0; i < vec.size(); ++i) {
        STXXL_CHECK_EQUAL(vec[i].get_i(), expected_vals[i]);
    }

    // std::remove_if would fail because it makes use of copy assignment.
    // We test stxxl's swap_remove_if implementation instead.
    // STL: vec.erase(std::remove_if(vec.begin(), vec.end(), test_eraser()), vec.end());
    vec.erase(stxxl::swap_remove_if(vec.begin(), vec.end(), test_eraser()), vec.end());

    // Check the values.
    int expected_vals2[10] = { 1, 2, 3, 5, 6, 7, 9 };
    STXXL_CHECK_EQUAL(vec.size(), 7);
    for (unsigned i = 0; i < vec.size(); ++i) {
        STXXL_CHECK_EQUAL(vec[i].get_i(), expected_vals2[i]);
    }

    // Clear the vector.
    vec.clear();
    STXXL_CHECK(vec.empty());

    // Resize to 100 and overwrite the last value.
    // Content after resize and overwrite: {...,-1,100}
    // Note: clear() does not overwrite any values in the underlaying array.
    vec.resize(20);
    Test t(11);
    std::swap(vec[19], t);

    // Check the values.
    STXXL_CHECK_EQUAL(vec.size(), 20);
    STXXL_CHECK_EQUAL(vec[19].get_i(), 11);

    return EXIT_SUCCESS;
}
