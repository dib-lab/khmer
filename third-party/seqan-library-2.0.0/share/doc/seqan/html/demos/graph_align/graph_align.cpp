#include <iostream>

#include <seqan/basic.h>
#include <seqan/stream.h>
#include <seqan/align.h>
#include <seqan/graph_align.h>

using namespace seqan;

int main()
{
    typedef StringSet<DnaString, Dependent<> > TStringSet;
    typedef Graph<Alignment<TStringSet, void> > TAlignmentGraph;

    // Define scorings scheme.
    //
    // In this case, affine gap costs with match = 0, mismatch = -1,
    // gapextend = -1 and gapopen = -2.
    Score<int> scoringScheme(0, -1, -1, -2);

    // Define the two sequence to be allocated.
    DnaString seq1 = "atcgaatgcgga";
    DnaString seq2 = "actcgttgca";

    // Create StringSet with these two sequences and construct an
    // AlignmentGraph for them.
    TStringSet stringSet;
    appendValue(stringSet, seq1);
    appendValue(stringSet, seq2);
    TAlignmentGraph alignmentGraph(stringSet);

    // Compute global alignment of seq1 and seq1 using Gotoh's algorithm and
    // print score and alignment graph.
    int score = globalAlignment(alignmentGraph, scoringScheme, Gotoh());
    std::cout << "Score = " << score << "\n"
              << alignmentGraph << std::endl;
    return 0;
}
