/***************************************************************************
 *  tests/algo/test_random_shuffle.cpp
 *
 *  Part of the STXXL. See http://stxxl.sourceforge.net
 *
 *  Copyright (C) 2007 Manuel Krings
 *  Copyright (C) 2007 Markus Westphal <mail@markuswestphal.de>
 *  Copyright (C) 2009 Andreas Beckmann <beckmann@cs.uni-frankfurt.de>
 *
 *  Distributed under the Boost Software License, Version 1.0.
 *  (See accompanying file LICENSE_1_0.txt or copy at
 *  http://www.boost.org/LICENSE_1_0.txt)
 **************************************************************************/

// TODO: test both vector and non-vector variant of random_shuffle
// TODO: test recursion, improve verboseness

//! \example algo/test_random_shuffle.cpp
//! Test \c stxxl::random_shuffle()

#include <stxxl/vector>
#include <stxxl/random_shuffle>

template <typename type>
struct counter
{
    type value;
    counter(type v = type(0)) : value(v) { }
    type operator () ()
    {
        type old_val = value;
        value++;
        return old_val;
    }
};

void long_test()
{
    typedef stxxl::vector<int> ext_vec_type;
    ext_vec_type STXXLVector(1024 * 1024 * 256 / sizeof(int));

    STXXL_MSG("Filling vector with increasing values...");
    stxxl::generate(STXXLVector.begin(), STXXLVector.end(),
                    counter<int>(), 4);

    stxxl::uint64 i;

    STXXL_MSG("Begin: ");
    for (i = 0; i < 10; i++)
        STXXL_MSG(STXXLVector[i]);

    STXXL_MSG("End: ");
    for (i = STXXLVector.size() - 10; i < STXXLVector.size(); i++)
        STXXL_MSG(STXXLVector[i]);

    STXXL_MSG("Permute randomly...");
    stxxl::random_shuffle(STXXLVector.begin(), STXXLVector.end(), 1024 * 1024 * 128);

    STXXL_MSG("Begin: ");
    for (i = 0; i < 10; i++)
        STXXL_MSG(STXXLVector[i]);

    STXXL_MSG("End: ");
    for (i = STXXLVector.size() - 10; i < STXXLVector.size(); i++)
        STXXL_MSG(STXXLVector[i]);
}

void short_test()
{
    STXXL_STATIC_ASSERT(sizeof(int) == 4);
    typedef stxxl::VECTOR_GENERATOR<int, 1, 2, 4096>::result vector_type;
    vector_type::size_type i;
    vector_type v(2048);
    for (i = 0; i < v.size(); ++i)
        v[i] = (int)(i / 1024 + 1);

    std::cout << v.size() << std::endl;
    std::cout << "before shuffle:" << std::endl;
    for (i = 0; i < v.size(); ++i) {
        std::cout << v[i] << " ";
        if ((i & 511) == 511)
            std::cout << std::endl;
    }
    std::cout << std::endl;
    v.flush();
    stxxl::random_shuffle(v.begin() + 512, v.begin() + 512 + 1024, 64 * 1024);
    std::cout << "after shuffle:" << std::endl;
    for (i = 0; i < v.size(); ++i) {
        std::cout << v[i] << " ";
        if ((i & 511) == 511)
            std::cout << std::endl;
    }
    std::cout << std::endl;
}

int main()
{
    short_test();
    long_test();
}
// vim: et:ts=4:sw=4
