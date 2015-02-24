#include <iostream>
#include <fstream>

#include <seqan/index.h>
#include <seqan/alignment_free.h>
#include <seqan/statistics.h>
#include <seqan/seq_io.h>

using namespace seqan;

int main()
{
    // Build path to serialized MarkovModel.
    CharString mmPath = SEQAN_PATH_TO_ROOT();
    append(mmPath, "/demos/statistics/zscore_example_mm.3");

    // Open the file.
    FILE * mmFile = fopen(toCString(mmPath), "rb");
    if (!mmFile)
    {
        std::cerr << "ERROR: Could not open " << mmPath << "\n";
        return 1;
    }

    // Create MarkovModel of order 3 and load it from the file.
    MarkovModel<Dna> mm(3);
    read(mmFile, mm);
    fclose(mmFile);  // close file again

    // Build set of words that we want to compute the zscore of.
    DnaString word = "CCCAAAGC";

    // Compute variance.
    double variance = 0;
    int n = 10000;  // assumed text length
    calculateVariance(variance, word, mm, n);
    std::cout << "variance: " << variance << "\n";

    return 0;
}
