/***************************************************************************
 *  examples/algo/phonebills.cpp
 *
 *  Part of the STXXL. See http://stxxl.sourceforge.net
 *
 *  Copyright (C) 2004 Roman Dementiev <dementiev@mpi-sb.mpg.de>
 *
 *  Distributed under the Boost Software License, Version 1.0.
 *  (See accompanying file LICENSE_1_0.txt or copy at
 *  http://www.boost.org/LICENSE_1_0.txt)
 **************************************************************************/

#define USE_STXXL

//! [prolog]
#include <stxxl.h>      // STXXL header
#include <algorithm>    // for STL std::sort
#include <vector>       // for STL std::vector
#include <fstream>      // for std::fstream
#include <ctime>        // for time_t type

#define CT_PER_MIN 2    // subscribers pay 2 cent per minute

struct LogEntry         // the event log data structure
{
    long long int from; // callers number (64 bit integer)
    long long int to;   // destination number (64 bit int)
    time_t timestamp;   // time of event
    int event;          // event type 1 - call started
                        //            2 - call ended
};

// input operator used for reading from the file
std::istream& operator >> (std::istream& i, LogEntry& entry)
{
    i >> entry.from;
    i >> entry.to;
    i >> entry.timestamp;
    i >> entry.event;
    return i;
}

// output operator used for writing to file
std::ostream& operator << (std::ostream& i, const LogEntry& entry)
{
    i << entry.from << " ";
    i << entry.to << " ";
    i << entry.timestamp << " ";
    i << entry.event;
    return i;
}

// unary function used for producing the bills
struct ProduceBill
{
    std::ostream& out;  // stream for outputting
    // the bills
    unsigned sum;       // current subscribers debit
    LogEntry last;      // the last record

    ProduceBill(std::ostream& o_) : out(o_), sum(0) { last.from = -1; }

    void operator () (const LogEntry& e)
    {
        if (last.from == e.from)
        {
            // either the last event was 'call started' and current event
            // is 'call ended' or the last event was 'call ended' and
            // current event is 'call started'
            assert((last.event == 1 && e.event == 2) ||
                   (last.event == 2 && e.event == 1));

            if (e.event == 2)            // call ended
                sum += CT_PER_MIN * (unsigned int)(e.timestamp - last.timestamp) / 60;
        }
        else if (last.from != -1)
        {
            assert(last.event == 2);            // must be 'call ended'
            assert(e.event == 1);               // must be 'call started'

            // output the total sum
            out << last.from << "; " << (sum / 100) << " EUR " << (sum % 100) <<
                " ct" << std::endl;

            sum = 0;             // reset the sum
        }

        last = e;
    }
};

struct SortByCaller
{
    // comparison function
    bool operator () (const LogEntry& a, const LogEntry& b) const
    {
        return a.from < b.from ||
               (a.from == b.from && a.timestamp < b.timestamp) ||
               (a.from == b.from && a.timestamp == b.timestamp &&
                a.event < b.event);
    }
    // min sentinel = value which is strictly smaller that any input element
    static LogEntry min_value()
    {
        LogEntry dummy;
        dummy.from = std::numeric_limits<long long int>::min();
        dummy.timestamp = std::numeric_limits<time_t>::min();
        dummy.event = std::numeric_limits<int>::min();
        return dummy;
    }
    // max sentinel = value which is strictly larger that any input element
    static LogEntry max_value()
    {
        LogEntry dummy;
        dummy.from = std::numeric_limits<long long int>::max();
        dummy.timestamp = std::numeric_limits<time_t>::max();
        dummy.event = std::numeric_limits<int>::max();
        return dummy;
    }
};

void print_usage(const char* program)
{
    std::cout << "Usage: " << program << " logfile main billfile" << std::endl;
    std::cout << " logfile  - file name of the input" << std::endl;
    std::cout << " main     - memory to use (in MiB)" << std::endl;
    std::cout << " billfile - file name of the output" << std::endl;
}
//! [prolog]

#ifdef USE_STXXL
typedef stxxl::vector<LogEntry> vector_type;
#else
typedef std::vector<LogEntry> vector_type;
#endif

int main(int argc, char* argv[])
{
    if (argc < 4)
    {
        print_usage(argv[0]);
        return 0;
    }

    // read from the file
    std::fstream in(argv[1], std::ios::in);
    vector_type v;

    std::copy(std::istream_iterator<LogEntry>(in),
              std::istream_iterator<LogEntry>(),
              std::back_inserter(v));

    // sort by callers
#ifndef USE_STXXL

    std::sort(v.begin(), v.end(), SortByCaller());
    std::fstream out(argv[3], std::ios::out);
    // output bills
    std::for_each(v.begin(), v.end(), ProduceBill(out));

#else
    const stxxl::internal_size_type M = atol(argv[2]) * 1024 * 1024;

    stxxl::sort(v.begin(), v.end(), SortByCaller(), M);
    std::fstream out(argv[3], std::ios::out);
    // output bills
    stxxl::for_each(v.begin(), v.end(), ProduceBill(out), M / vector_type::block_type::raw_size);

#endif

    return 0;
}
//! [example]
