#include <iostream>
#include <seqan/misc/interval_tree.h>

using namespace seqan;

int main()
{
    // fill a string of intervals and keys
    typedef IntervalAndCargo<int, CharString> TInterval;
    String<TInterval> intervals;
    appendValue(intervals, TInterval(100, 1000, "gene"));
    appendValue(intervals, TInterval(100, 300, "exon1"));
    appendValue(intervals, TInterval(150, 250, "coding1"));
    appendValue(intervals, TInterval(500, 800, "exon2"));
    appendValue(intervals, TInterval(600, 700, "coding2"));

    // create IntervalTree of that string
    IntervalTree<int, CharString> tree(intervals);

    // find intervals that overlap the query interval [550,900)
    String<CharString> results;
    findIntervals(results, tree, 550, 900);

    // output corresponding keys
    for (unsigned i = 0; i < length(results); ++i)
        std::cout << results[i] << std::endl;

    return 0;
}
