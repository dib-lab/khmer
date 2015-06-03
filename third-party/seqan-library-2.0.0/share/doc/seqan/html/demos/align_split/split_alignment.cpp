#include <seqan/align.h>
#include <seqan/align_split.h>
#include <seqan/file.h>  // output of String
#include <seqan/sequence.h>
#include <seqan/score.h>

using namespace seqan;

int main()
{
    std::cout << "Situation\n"
              << "\n"
              << "REF  AGCATGTTAGATAAGATAG-----------CTGTGCTAGTAGGCAGTCAGCGCCAT\n"
              << "READ AGCATGTTAGATAAGATAGCCCCCCCCCCCCTGTGCTAGTAGGCAGTCAGCGCCAT\n"
              << "\n";

    // Demo for split alignment where the read contains an insertion with
    // respect to the reference.  The input of the function is the infix of
    // reference earlier identified.
    Dna5String ref =  "AGCATGTTAGATAAGATAG"         "CTGTGCTAGTAGGCAGTCAGCGCCAT";
    Dna5String read = "AGCATGTTAGATAAGATAGCCCCCCCCCCCCTGTGCTAGTAGGCAGTCAGCGCCAT";

    // Prepare Gaps objects.  We need one for the left part and one for the
    // right part of the alignment.
    Align<Dna5String> alignL;
    resize(rows(alignL), 2);
    assignSource(row(alignL, 0), ref);
    assignSource(row(alignL, 1), read);
    Align<Dna5String> alignR;
    resize(rows(alignR), 2);
    assignSource(row(alignR, 0), ref);
    assignSource(row(alignR, 1), read);

    // Define scoring scheme.
    Score<int, Simple> scoringScheme(1, -1, -1);

    // Call split alignment function.
    splitAlignment(alignL, alignR, scoringScheme);

    // Print resulting alignment to stdout.
    std::cout << "Resulting alignments\n"
              << "\n"
              << "Left\n"
              << alignL
              << "Right\n"
              << alignR
              << "\n";

    // Get relevant clipping positions.
    int refSplitPosition = toSourcePosition(row(alignL, 0), clippedEndPosition(row(alignL, 0)));
    SEQAN_ASSERT_EQ(refSplitPosition, toSourcePosition(row(alignR, 0), 0));
    int readSplitLPosition = toSourcePosition(row(alignL, 1), clippedEndPosition(row(alignL, 1)));
    int readSplitRPosition = toSourcePosition(row(alignR, 1), 0);

    std::cout << "refSplitPosition   == " << refSplitPosition << "\n"
              << "readSplitLPosition == " << readSplitLPosition << "\n"
              << "readSplitRPosition == " << readSplitRPosition << "\n\n";

    // Print sequence parts to stdout.
    std::cout << "Reference Left  " << prefix(ref, refSplitPosition) << "\n"
              << "Reference Right " << suffix(ref, refSplitPosition) << "\n"
              << "\n"
              << "Read Left       " << prefix(read, readSplitLPosition) << "\n"
              << "Read Center     " << infix(read, readSplitLPosition, readSplitRPosition) << "\n"
              << "Read Right      " << suffix(read, readSplitRPosition) << "\n";

    return 0;
}
