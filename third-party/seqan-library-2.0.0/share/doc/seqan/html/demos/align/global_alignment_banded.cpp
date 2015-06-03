#include <iostream>

#include <seqan/basic.h>
#include <seqan/align.h>
#include <seqan/sequence.h>
#include <seqan/stream.h>  // for printint strings

using namespace seqan;

int main()
{
    Dna5String seqH = "CGATT";
    Dna5String seqV = "CGAAATT";

    Align<Dna5String> align;
    resize(rows(align), 2);
    assignSource(row(align, 0), seqH);
    assignSource(row(align, 1), seqV);

    Score<int, Simple> scoringScheme(2, -1, -2);
    AlignConfig<> alignConfig;

    int lDiag = -2;
    int uDiag = 2;

    int result = globalAlignment(align, scoringScheme, alignConfig, lDiag, uDiag);

    std::cout << "Score: " << result << "\n";
    std::cout << "The resulting alignment is\n"
              << align << "\n";

    return 0;
}
