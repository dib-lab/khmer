/***************************************************************************
 *  tools/benchmarks/benchmark_naive_matrix.cpp
 *
 *  Part of the STXXL. See http://stxxl.sourceforge.net
 *
 *  Copyright (C) 2010 Johannes Singler <singler@kit.edu>
 *
 *  Distributed under the Boost Software License, Version 1.0.
 *  (See accompanying file LICENSE_1_0.txt or copy at
 *  http://www.boost.org/LICENSE_1_0.txt)
 **************************************************************************/

//! \example containers/test_matrix.cpp
//! This is an example of how to represent a 2D matrix by using an \c stxxl::vector.
//! DISCLAIMER: The approach is very basic and the matrix multiplication is NOT I/O-EFFICIENT.

#include <iostream>
#include <stxxl/vector>

template <typename T>
class matrix2d
{
    stxxl::vector<T> v;
    stxxl::unsigned_type width, height;

public:
    matrix2d(stxxl::unsigned_type width, stxxl::unsigned_type height) : width(width), height(height)
    {
        v.resize(width * height);
    }

    stxxl::unsigned_type get_width() const
    {
        return width;
    }

    stxxl::unsigned_type get_height() const
    {
        return height;
    }

    T & element(stxxl::unsigned_type x, stxxl::unsigned_type y)
    {
        //row-major
        return v[y * width + x];
    }

    const T & const_element(stxxl::unsigned_type x, stxxl::unsigned_type y) const
    {
        //row-major
        return v[y * width + x];
    }
};

template <typename T>
void matrix_multiply(const matrix2d<T>& a, const matrix2d<T>& b, matrix2d<T>& c)
{
    assert(a.get_width() == b.get_height());
    assert(b.get_height() == c.get_height());
    assert(b.get_width() == c.get_width());

    for (stxxl::unsigned_type y = 0; y < c.get_height(); ++y)
        for (stxxl::unsigned_type x = 0; x < c.get_width(); ++x)
        {
            c.element(x, y) = 0;
            for (stxxl::unsigned_type t = 0; t < a.get_width(); ++t)
                c.element(x, y) += a.const_element(t, y) * b.const_element(x, t);
        }
}

int main()
{
    const stxxl::unsigned_type width = 416, height = 416;

    try
    {
        matrix2d<double> a(width, height), b(width, height), c(width, height);

        for (stxxl::unsigned_type y = 0; y < height; ++y)
            for (stxxl::unsigned_type x = 0; x < width; ++x)
            {
                a.element(x, y) = 1.0;
                b.element(x, y) = 1.0;
                c.element(x, y) = 1.0;
            }

        matrix_multiply(a, b, c);

        std::cout << c.const_element(0, 0) << std::endl;
        std::cout << c.const_element(width - 1, height - 1) << std::endl;
    }
    catch (const std::exception& ex)
    {
        STXXL_MSG("Caught exception: " << ex.what());
    }
    catch (...)
    {
        STXXL_MSG("Caught unknown exception.");
    }

    return 0;
}
