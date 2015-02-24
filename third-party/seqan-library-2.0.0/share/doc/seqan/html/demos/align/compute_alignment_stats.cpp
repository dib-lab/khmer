#include <iostream>

#include <seqan/align.h>
#include <seqan/sequence.h>

using namespace seqan;

int main()
{
    // Create an alignment between subject and query.
    Peptide subject =
        "MGLSDGEWQLVLNVWGKVEADIPGHGQEVLIRLFKGHPETLEKFDKFKHLKSEDEMKASE"
        "DLKKHGATVLTALGGILKKKGHHEAEIKPLAQSHATKHKIPVKYLEFISECIIQVLQSKH"
        "PGDFGADAQGAMNKALELFRKDMASNYK";
    Peptide query =
        "MSLTKTERTIIVSMWAKISTQADTIGTETLERLFLSHPQTKTYFPHFDLHPGSA"
        "QLRAHGSKVVAAVGDAVKSIDDIGGALSKLSELHAYILRVDPVNFKLLSHCLLVTLAARF"
        "PADFTAEAHAAWDKFLSVTEKYR";

    Align<Peptide> align;
    resize(rows(align), 2);
    setSource(row(align, 0), subject);
    setSource(row(align, 1), query);

    Blosum62 scoringScheme(-1, -12);
    globalAlignment(align, scoringScheme);

    // Compute the statistics of the alignment.
    AlignmentStats stats;
    int scoreVal = computeAlignmentStats(stats, align, scoringScheme);
    SEQAN_ASSERT_EQ(scoreVal, stats.alignmentScore);
    std::cout << align
              << "gap opens:           " << stats.numGapOpens << "\n"
              << "gap extensions:      " << stats.numGapExtensions << "\n"
              << "num insertions:      " << stats.numInsertions << "\n"
              << "num deletions:       " << stats.numDeletions << "\n"
              << "num matches:         " << stats.numMatches << "\n"
              << "num mismatches:      " << stats.numMismatches << "\n"
              << "num positive scores: " << stats.numPositiveScores << "\n"
              << "num negative scores: " << stats.numNegativeScores << "\n"
              << "percent similarity:  " << stats.alignmentSimilarity << "\n"
              << "percent identity:    " << stats.alignmentIdentity << "\n\n\n";

    // Clip alignment rows and compute score of this view.
    setClippedEndPosition(row(align, 0), 100);
    setClippedEndPosition(row(align, 1), 100);
    setClippedBeginPosition(row(align, 0), 5);
    setClippedBeginPosition(row(align, 1), 5);

    scoreVal = computeAlignmentStats(stats, align, scoringScheme);
    SEQAN_ASSERT_EQ(scoreVal, stats.alignmentScore);
    std::cout << "Clipping alignment to (5, 100)\n"
              << align
              << "gap opens:           " << stats.numGapOpens << "\n"
              << "gap extensions:      " << stats.numGapExtensions << "\n"
              << "num insertions:      " << stats.numInsertions << "\n"
              << "num deletions:       " << stats.numDeletions << "\n"
              << "num matches:         " << stats.numMatches << "\n"
              << "num mismatches:      " << stats.numMismatches << "\n"
              << "num positive scores: " << stats.numPositiveScores << "\n"
              << "num negative scores: " << stats.numNegativeScores << "\n"
              << "percent similarity:  " << stats.alignmentSimilarity << "\n"
              << "percent identity:    " << stats.alignmentIdentity << "\n";

    return 0;
}
