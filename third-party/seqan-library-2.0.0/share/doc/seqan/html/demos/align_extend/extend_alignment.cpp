#include <iostream>

#include <seqan/align.h>
#include <seqan/align_extend.h>
#include <seqan/sequence.h>

using namespace seqan;

int main()
{
    Score<int> sc(2, -1, -2);

    Align<Infix<CharString const>::Type> align;
    resize(rows(align), 2);

    // We create the following initial situation with subject/query and the
    // infixes thereof.
    //
    //                         infixes/seeds
    //                             <-->
    // subject  NNNNNNNNNNTTCCGGGACGGTACACACACGGGGGGGGGG
    // query               CTCGGGACGGTACAGGCACGGTTTTTTTT
    CharString subject = "NNNNNNNNNNTTCCGGGACGGTACACACACGGGGGGGGGG";
    CharString query   = "CTCGGGACGGTACAGGCACGGTTTTTTTT";
    assignSource(row(align, 0), infix(subject, 19, 23));
    assignSource(row(align, 1), infix(query, 8, 12));
    int score = globalAlignment(align, sc);

    std::cout << "Initial alignment of infixes (score == " << score << ")\n\n"
              << align;

    // The alignment starts at diagonal (23 - 19) = 4.  A band of 4 in each direction has
    // the following diagonals.
    int lDiag = 0, uDiag = 8;
    // Set the x-Drop value to 5.
    int xDrop = 5;
    Tuple<unsigned, 4> positions = { {19u, 8u, 23u, 12u} };
    score = extendAlignment(align, score, subject, query, positions, EXTEND_BOTH,
                            lDiag, uDiag, xDrop, sc);

    std::cout << "Resulting alignment (score == " << score << ")\n\n"
              << align;

    std::cout << "source(row(align, 0)) == " << source(row(align, 0)) << " (full sequence)\n"
              << "source(row(align, 1)) == " << source(row(align, 1)) << " (full sequence)\n"
              << "\n"
              << "clipping positions of row 0: " << clippedBeginPosition(row(align, 0))
              << ", " << clippedEndPosition(row(align, 0)) << "\n"
              << "clipping positions of row 1: " << clippedBeginPosition(row(align, 1))
              << ", " << clippedEndPosition(row(align, 1)) << "\n";


    return 0;
}
