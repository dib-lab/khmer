#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/stream.h>  // for I/O
#include <seqan/align.h>
#include <seqan/score.h>

using namespace seqan;

int main()
{
    StringSet<DnaString> stringSet;
    appendValue(stringSet, "ACGTGGATCGGTGACTTACGGACTG");
    appendValue(stringSet, "ACGTGTTGAGTGATTACGGACTG");

    Align<DnaString> align(stringSet);                // Initialize the Align object using a StringSet.

    Score<int> scoreScheme(5, -6, -2, -11);           // Initialize Score object with affine gap costs.
    int score = globalAlignment(align, scoreScheme);  // Computes a global alignment usign the defined scoring scheme.

    std::cout << "Score = " << score << std::endl;
    std::cout << "Alignment:" << std::endl;
    std::cout << align << std::endl;

    return 0;
}
