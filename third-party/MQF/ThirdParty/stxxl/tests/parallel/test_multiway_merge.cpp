/***************************************************************************
 *  tests/parallel/test_multiway_merge.cpp
 *
 *  Part of the STXXL. See http://stxxl.sourceforge.net
 *
 *  Copyright (C) 2014-2015 Timo Bingmann <tb@panthema.net>
 *
 *  Distributed under the Boost Software License, Version 1.0.
 *  (See accompanying file LICENSE_1_0.txt or copy at
 *  http://www.boost.org/LICENSE_1_0.txt)
 **************************************************************************/

#include <stxxl/bits/parallel.h>
#include <stxxl/bits/verbose.h>
#include <stxxl/random>
#include <iostream>

struct Something
{
    unsigned int a, b;

    explicit Something(unsigned int x = 0)
        : a(x), b(x * x)
    { }

    bool operator < (const Something& other) const
    {
        return a < other.a;
    }

    bool operator == (const Something& other) const
    {
        return (a == other.a) && (b == other.b);
    }

    friend std::ostream& operator << (std::ostream& os, const Something& s)
    {
        return os << '(' << s.a << ',' << s.b << ')';
    }
};

template class stxxl::parallel::LoserTreeCopy<false, Something>;
template class stxxl::parallel::LoserTreeCopy<true, Something>;
template class stxxl::parallel::LoserTreeReference<Something>;
template class stxxl::parallel::LoserTreePointer<false, Something>;
template class stxxl::parallel::LoserTreePointer<true, Something>;
template class stxxl::parallel::LoserTreeCopyUnguarded<false, Something>;
template class stxxl::parallel::LoserTreeCopyUnguarded<true, Something>;
template class stxxl::parallel::LoserTreePointerUnguarded<false, Something>;
template class stxxl::parallel::LoserTreePointerUnguarded<true, Something>;

template <typename ValueType, bool Stable, bool Sentinels>
void test_vecs(unsigned int vecnum)
{
    static const bool debug = false;

    // construct many vectors of sorted random numbers

    stxxl::random_number32 rnd;
    std::vector<std::vector<ValueType> > vec(vecnum);

    size_t totalsize = 0;
    std::vector<ValueType> output, correct;

    for (size_t i = 0; i < vecnum; ++i)
    {
        // determine number of items in stream
        size_t inum = (rnd() % 128) + 64;
        if (i == 3) inum = 0; // add an empty sequence
        vec[i].resize(inum);
        totalsize += inum;

        // pick random items
        for (size_t j = 0; j < inum; ++j)
            vec[i][j] = ValueType(rnd() % (vecnum * 20));

        std::sort(vec[i].begin(), vec[i].end());

        // put items into correct vector as well
        correct.insert(correct.end(), vec[i].begin(), vec[i].end());
    }

    // append sentinels
    if (Sentinels)
    {
        for (size_t i = 0; i < vecnum; ++i)
            vec[i].push_back(ValueType(std::numeric_limits<unsigned int>::max()));
    }

    // prepare output and correct vector
    output.resize(totalsize);
    std::sort(correct.begin(), correct.end());

    if (debug)
    {
        for (size_t i = 0; i < vecnum; ++i)
        {
            std::cout << "vec[" << i << "]:";
            for (size_t j = 0; j < vec[i].size(); ++j)
                std::cout << " " << vec[i][j];
            std::cout << "\n";
        }
    }

    // construct vector of input iterator ranges
    typedef typename std::vector<ValueType>::iterator input_iterator;

    std::vector<std::pair<input_iterator, input_iterator> > sequences(vecnum);

    for (size_t i = 0; i < vecnum; ++i)
    {
        sequences[i] = std::make_pair(vec[i].begin(), vec[i].end() - (Sentinels ? 1 : 0));

        STXXL_CHECK(stxxl::is_sorted(vec[i].begin(), vec[i].end()));
    }

    if (!Sentinels) {
        if (!Stable)
            stxxl::potentially_parallel::multiway_merge(
                sequences.begin(), sequences.end(),
                output.begin(), totalsize, std::less<ValueType>());
        else
            stxxl::potentially_parallel::multiway_merge_stable(
                sequences.begin(), sequences.end(),
                output.begin(), totalsize, std::less<ValueType>());
    }
    else {
        if (!Stable)
            stxxl::potentially_parallel::multiway_merge_sentinels(
                sequences.begin(), sequences.end(),
                output.begin(), totalsize, std::less<ValueType>());
        else
            stxxl::potentially_parallel::multiway_merge_stable_sentinels(
                sequences.begin(), sequences.end(),
                output.begin(), totalsize, std::less<ValueType>());
    }

#ifdef TODO
    // sequential version's losertrees with sentinels does not work correctly!

    if (!Sentinels) {
        stxxl::parallel::sequential_multiway_merge<Stable, false>(
            sequences.begin(), sequences.end(),
            output.begin(), totalsize, std::less<ValueType>());
    }
    else {
        stxxl::parallel::sequential_multiway_merge<Stable, true>(
            sequences.begin(), sequences.end(),
            output.begin(), totalsize, std::less<ValueType>());
    }
#endif

    if (debug)
    {
        for (size_t i = 0; i < output.size(); ++i)
            std::cout << output[i] << " ";
    }

    STXXL_CHECK(output == correct);
}

void test_all()
{
    // run multiway merge tests for 0..256 sequences
    for (unsigned int n = 0; n <= 256; n += 1 + n / 64 + n / 128)
    {
        std::cout << "testing multiway_merge with " << n << " players\n";

        test_vecs<Something, false, false>(n);    // unstable, no-sentinels
        test_vecs<Something, true, false>(n);     // stable, no-sentinels
        test_vecs<Something, false, true>(n);     // unstable, sentinels
        test_vecs<Something, true, true>(n);      // stable, sentinels

        test_vecs<unsigned int, false, false>(n); // unstable, no-sentinels
        test_vecs<unsigned int, true, false>(n);  // stable, no-sentinels
        test_vecs<unsigned int, false, true>(n);  // unstable, sentinels
        test_vecs<unsigned int, true, true>(n);   // stable, sentinels
    }
}

int main()
{
    stxxl::parallel::SETTINGS::multiway_merge_splitting = stxxl::parallel::SETTINGS::EXACT;
    test_all();

    stxxl::parallel::SETTINGS::multiway_merge_splitting = stxxl::parallel::SETTINGS::SAMPLING;
    test_all();

    return 0;
}
