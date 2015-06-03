#include <seqan/align.h>
#include <seqan/basic.h>
#include <seqan/stream.h>  // for printing Seqan strings

using namespace seqan;

int main()
{
    Dna5String seq = "CGGTATC";

    // Construct the gaps object.
    Gaps<Dna5String> gaps(seq);

    // Insert gaps (view position, unclipped as of yet).  Note that inserting
    // the gaps changes the coordinates.
    insertGap(gaps, 3);
    insertGap(gaps, 5);
    insertGaps(gaps, 7, 2);

    // Clip the gaps object.  Note that the coordinates are in unclipped view position.
    setClippedBeginPosition(gaps, 1);
    setClippedEndPosition(gaps, 8);

    // Print the resulting gaps object and some coordinate transformation.
    std::cout << "Resulting gaps: " << gaps << "\n"
              << "toSourcePosition(gaps, 0) == " << toSourcePosition(gaps, 0) << "\n"
              << "toSourcePosition(gaps, 4) == " << toSourcePosition(gaps, 4) << "\n"
              << "toViewPosition(gaps, 0) == " << toViewPosition(gaps, 0) << "\n"
              << "toViewPosition(gaps, 5) == " << toViewPosition(gaps, 6) << "\n";

    return 0;
}
