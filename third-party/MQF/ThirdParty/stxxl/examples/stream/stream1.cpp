/***************************************************************************
 *  examples/stream/stream1.cpp
 *
 *  This file contains the example snippets from the stream tutorial.
 *
 *  Part of the STXXL. See http://stxxl.sourceforge.net
 *
 *  Copyright (C) 2012-2013 Timo Bingmann <tb@panthema.net>
 *
 *  Distributed under the Boost Software License, Version 1.0.
 *  (See accompanying file LICENSE_1_0.txt or copy at
 *  http://www.boost.org/LICENSE_1_0.txt)
 **************************************************************************/

#include <stxxl/stream>
#include <stxxl/vector>
#include <stxxl/sorter>

#include <vector>
#include <climits>

struct counter_object
{
    // This stream produces a sequence of integers.
    typedef int value_type;

private:
    // A class attribute to save the current value.
    int m_current_value;

public:
    // A constructor to set the initial value to 1.
    counter_object()
        : m_current_value(1)
    { }

    // The retrieve operator returning the current value.
    const value_type& operator * () const
    {
        return m_current_value;
    }

    // Increment operator advancing to the next integer.
    counter_object& operator ++ ()
    {
        ++m_current_value;
        return *this;
    }

    // Empty indicator, which in this case can check the current value.
    bool empty() const
    {
        return (m_current_value > 1000);
    }
};

template <typename InputStream>
struct squaring_object
{
    // This stream produces a sequence of integers.
    typedef int value_type;

private:
    // A reference to another stream of integers, which are our input.
    InputStream& m_input_stream;

    // A temporary value buffer to hold the current square in for retrieval.
    value_type m_current_value;

public:
    // A constructor taking another stream of integers as input.
    squaring_object(InputStream& input_stream)
        : m_input_stream(input_stream)
    {
        if (!m_input_stream.empty())
        {
            m_current_value = *m_input_stream;
            m_current_value = m_current_value * m_current_value;
        }
    }

    // The retrieve operator returning the square of the input stream.
    const value_type& operator * () const
    {
        return m_current_value;
    }

    // Increment operator: handled by incrementing the input stream.
    squaring_object& operator ++ ()
    {
        ++m_input_stream;
        if (!m_input_stream.empty())
        {
            m_current_value = *m_input_stream;
            m_current_value = m_current_value * m_current_value;
        }
        return *this;
    }

    // Empty indicator: this stream is empty when the input stream is.
    bool empty() const
    {
        return m_input_stream.empty();
    }
};

// define comparator class: compare right-most decimal and then absolute value
struct CompareMod10
{
    // comparison operator() returning true if (a < b)
    inline bool operator () (int a, int b) const
    {
        if ((a % 10) == (b % 10))
            return a < b;
        else
            return (a % 10) < (b % 10);
    }

    // smallest possible integer value
    int min_value() const { return INT_MIN; }
    // largest possible integer value
    int max_value() const { return INT_MAX; }
};

int main()
{
    {
        counter_object counter;

        while (!counter.empty())
        {
            std::cout << *counter << " ";
            ++counter;
        }
        std::cout << std::endl;
    }

    {
        for (counter_object cnt; !cnt.empty(); ++cnt)
        {
            std::cout << *cnt << " ";
        }
        std::cout << std::endl;
    }

    {
        counter_object counter;
        squaring_object<counter_object> squares(counter);

        while (!squares.empty())
        {
            std::cout << *squares << " ";
            ++squares;
        }
        std::cout << std::endl;
    }

    {
        std::vector<int> intvector;
        // (fill intvector)

        // define stream class iterating over an integer vector
        typedef stxxl::stream::iterator2stream<std::vector<int>::const_iterator> intstream_type;

        // instantiate the stream object, iterate from begin to end of intvector.
        intstream_type intstream(intvector.begin(), intvector.end());

        // plug in squaring object after vector iterator stream.
        squaring_object<intstream_type> squares(intstream);
    }

    {
        stxxl::vector<int> intvector;
        // (fill intvector)

        // define stream class iterating over an integer vector
        typedef stxxl::stream::vector_iterator2stream<stxxl::vector<int>::const_iterator> intstream_type;

        // instantiate the stream object, iterate from begin to end of intvector.
        intstream_type intstream(intvector.begin(), intvector.end());

        // plug in squaring object after vector iterator stream.
        squaring_object<intstream_type> squares(intstream);
    }

    {
        // construct the squared counter stream
        counter_object counter;
        squaring_object<counter_object> squares(counter);

        // allocate vector of 100 integers
        std::vector<int> intvector(100);

        // materialize 100 integers from stream and put into vector
        stxxl::stream::materialize(squares, intvector.begin(), intvector.end());
    }

    {
        // construct the squared counter stream
        counter_object counter;
        squaring_object<counter_object> squares(counter);

        // allocate STXXL vector of 100 integers
        stxxl::vector<int> intvector(100);

        // materialize 100 integers from stream and put into STXXL vector
        stxxl::stream::materialize(squares, intvector.begin(), intvector.end());
    }

    {
        static const int ram_use = 10 * 1024 * 1024; // amount of memory to use in runs creation

        counter_object counter;                      // the counter stream from first examples

        // define a runs sorter for the counter stream which order by CompareMod10 object.
        typedef stxxl::stream::runs_creator<counter_object, CompareMod10> rc_counter_type;

        // instance of CompareMod10 comparator class
        CompareMod10 comparemod10;

        // instance of runs_creator which reads the counter stream.
        rc_counter_type rc_counter(counter, comparemod10, ram_use);

        // define a runs merger for the sorted runs from rc_counter.
        typedef stxxl::stream::runs_merger<rc_counter_type::sorted_runs_type, CompareMod10> rm_counter_type;

        // instance of runs_merger which merges sorted runs from rc_counter.
        rm_counter_type rm_counter(rc_counter.result(), comparemod10, ram_use);

        // read sorted stream: runs_merger also conforms to the stream interface.
        while (!rm_counter.empty())
        {
            std::cout << *rm_counter << " ";
            ++rm_counter;
        }
        std::cout << std::endl;
    }

    {
        static const int ram_use = 10 * 1024 * 1024;   // amount of memory to use in runs creation

        // define a runs sorter which accepts imperative push()s and orders by CompareMod10 object.
        typedef stxxl::stream::runs_creator<stxxl::stream::use_push<int>, CompareMod10> rc_counter_type;

        // instance of CompareMod10 comparator class.
        CompareMod10 comparemod10;

        // instance of runs_creator which waits for input.
        rc_counter_type rc_counter(comparemod10, ram_use);

        // write sequence of integers into runs
        for (int i = 1; i <= 1000; ++i)
            rc_counter.push(i);

        // define a runs merger for the sorted runs from rc_counter.
        typedef stxxl::stream::runs_merger<rc_counter_type::sorted_runs_type, CompareMod10> rm_counter_type;

        // instance of runs_merger which merges sorted runs from rc_counter.
        rm_counter_type rm_counter(rc_counter.result(), comparemod10, ram_use);

        // read sorted stream: runs_merger also conforms to the stream interface.
        while (!rm_counter.empty())
        {
            std::cout << *rm_counter << " ";
            ++rm_counter;
        }
        std::cout << std::endl;
    }

    {
        static const int ram_use = 10 * 1024 * 1024;   // amount of memory to use in runs creation

        // define a runs sorter which accepts imperative push()s and orders by CompareMod10 object.
        typedef stxxl::sorter<int, CompareMod10> sr_counter_type;

        // instance of CompareMod10 comparator class.
        CompareMod10 comparemod10;

        // instance of sorter which waits for input.
        sr_counter_type sr_counter(comparemod10, ram_use);

        // write sequence of integers into sorter, which creates sorted runs
        for (int i = 1; i <= 1000; ++i)
            sr_counter.push(i);

        // signal sorter that the input stream is finished and switch to output mode.
        sr_counter.sort();

        // read sorted stream: sorter also conforms to the stream interface.
        while (!sr_counter.empty())
        {
            std::cout << *sr_counter << " ";
            ++sr_counter;
        }
        std::cout << std::endl;
    }
}
