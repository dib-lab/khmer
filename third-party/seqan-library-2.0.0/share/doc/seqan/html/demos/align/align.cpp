#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/stream.h>  // for I/O
#include <seqan/align.h>
#include <seqan/score.h>

using namespace seqan;

int main()
{
    StringSet<DnaString> stringSet;
    appendValue(stringSet, "AGTTTAATCA");
    appendValue(stringSet, "AGTATACGA");

    Align<DnaString> align(stringSet);                        // Initialize the Align object using a StringSet.

    int score = globalAlignment(align, EditDistanceScore());  // Compute a global alingment using the Align object.

    std::cout << "score = " << score << std::endl;
    std::cout << "align\n" << align << std::endl;

    return 0;
}
