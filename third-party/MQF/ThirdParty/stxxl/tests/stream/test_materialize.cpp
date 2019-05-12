/***************************************************************************
 *  tests/stream/test_materialize.cpp
 *
 *  Part of the STXXL. See http://stxxl.sourceforge.net
 *
 *  Copyright (C) 2009 Andreas Beckmann <beckmann@cs.uni-frankfurt.de>
 *  Copyright (C) 2013 Timo Bingmann <tb@panthema.net>
 *
 *  Distributed under the Boost Software License, Version 1.0.
 *  (See accompanying file LICENSE_1_0.txt or copy at
 *  http://www.boost.org/LICENSE_1_0.txt)
 **************************************************************************/

#include <vector>
#include <stxxl/stream>
#include <stxxl/vector>
#include <stxxl/algorithm>

struct forty_two
{
    unsigned counter, length;

    forty_two(unsigned l) : counter(0), length(l) { }

    bool empty() const { return !(counter < length); }

    unsigned len() const { return length; }

    int operator * ()
    {
        STXXL_CHECK(!empty());
        return counter;
    }

    forty_two& operator ++ ()
    {
        STXXL_CHECK(!empty());
        ++counter;
        return *this;
    }

    forty_two & reset()
    {
        counter = 0;
        return *this;
    }
};

/*
    template <class OutputIterator_, class StreamAlgorithm_>
    OutputIterator_ materialize(StreamAlgorithm_ & in, OutputIterator_ out);

    template <class OutputIterator_, class StreamAlgorithm_>
    OutputIterator_ materialize(StreamAlgorithm_ & in, OutputIterator_ outbegin, OutputIterator_ outend);

    template <typename Tp_, typename AllocStr_, typename SzTp_, typename DiffTp_,
              unsigned BlkSize_, typename PgTp_, unsigned PgSz_, class StreamAlgorithm_>
    stxxl::vector_iterator<Tp_, AllocStr_, SzTp_, DiffTp_, BlkSize_, PgTp_, PgSz_>
    materialize(StreamAlgorithm_ & in,
                stxxl::vector_iterator<Tp_, AllocStr_, SzTp_, DiffTp_, BlkSize_, PgTp_, PgSz_> outbegin,
                stxxl::vector_iterator<Tp_, AllocStr_, SzTp_, DiffTp_, BlkSize_, PgTp_, PgSz_> outend,
                unsigned_type nbuffers = 0);

    template <typename Tp_, typename AllocStr_, typename SzTp_, typename DiffTp_,
              unsigned BlkSize_, typename PgTp_, unsigned PgSz_, class StreamAlgorithm_>
    stxxl::vector_iterator<Tp_, AllocStr_, SzTp_, DiffTp_, BlkSize_, PgTp_, PgSz_>
    materialize(StreamAlgorithm_ & in,
                stxxl::vector_iterator<Tp_, AllocStr_, SzTp_, DiffTp_, BlkSize_, PgTp_, PgSz_> out,
                unsigned_type nbuffers = 0);
*/

int generate_0()
{
    return 0;
}

template <typename VectorType>
void check_42_fill(VectorType& v, unsigned length)
{
    typename VectorType::const_iterator ci = v.begin();

    for (unsigned i = 0; i < length; ++i)
    {
        STXXL_CHECK(*ci == (int)i);
        ++ci;
    }

    for (unsigned i = length; i < v.size(); ++i)
    {
        STXXL_CHECK(*ci == 0);
        ++ci;
    }

    std::fill(v.begin(), v.end(), 0);
}

int main()
{
    stxxl::config::get_instance();

    {
        forty_two _42(42);

        // materialize into std vector
        std::vector<int> v(1000);
        std::generate(v.begin(), v.end(), generate_0);

        stxxl::stream::materialize(_42.reset(), v.begin());
        check_42_fill(v, _42.len());

        stxxl::stream::materialize(_42.reset(), v.begin(), v.end());
        check_42_fill(v, _42.len());
    }
    {
        forty_two _42(42);

        // materialize into stxxl vector
        stxxl::VECTOR_GENERATOR<int>::result v(1000);
        stxxl::generate(v.begin(), v.end(), generate_0, 42);

        stxxl::stream::materialize(_42.reset(), v.begin());
        check_42_fill(v, _42.len());

        stxxl::stream::materialize(_42.reset(), v.begin(), 42);
        check_42_fill(v, _42.len());

        stxxl::stream::materialize(_42.reset(), v.begin(), v.end());
        check_42_fill(v, _42.len());

        stxxl::stream::materialize(_42.reset(), v.begin(), v.end(), 42);
        check_42_fill(v, _42.len());
    }
    {
        forty_two _42mill(42 * 1000000);

        // materialize into larger stxxl vector (to cross block boundaries)
        stxxl::VECTOR_GENERATOR<int>::result v(60 * 1000000);
        stxxl::generate(v.begin(), v.end(), generate_0, 42);

        stxxl::stream::materialize(_42mill.reset(), v.begin());
        check_42_fill(v, _42mill.len());

        stxxl::stream::materialize(_42mill.reset(), v.begin(), 42);
        check_42_fill(v, _42mill.len());

        stxxl::stream::materialize(_42mill.reset(), v.begin(), v.end());
        check_42_fill(v, _42mill.len());

        stxxl::stream::materialize(_42mill.reset(), v.begin(), v.end(), 42);
        check_42_fill(v, _42mill.len());
    }
}
