#include <seqan/basic.h>
#include <seqan/seq_io.h>
#include <seqan/sequence.h>

using namespace seqan;

int main()
{
    CharString path = SEQAN_PATH_TO_ROOT();
    append(path, "/demos/seq_io/example.fa");

    FaiIndex faiIndex;

    // Try to read the FAI index.
    if (!open(faiIndex, toCString(path)))
    {
        std::cerr << "Could not read the FAI index.  Not fatal, we can just build it.\n";
        return 1;
    }

    // Try to build the FAI index (in memory) if reading was unsuccessful.  If
    // building into memory succeeded, we try to write it out.
    if (!build(faiIndex, toCString(path)))
    {
        std::cerr << "FATAL: Could not build FAI index.\n";
        return 1;
    }

    if (!save(faiIndex))
    {
        std::cerr << "FATAL: Could not write out FAI index after building.\n";
        return 1;
    }

    // Now, read the first 1000 characters of chr1.
    unsigned idx = 0;
    if (!getIdByName(idx, faiIndex, "chr"))
    {
        std::cerr << "FATAL: chr1 not found in FAI index.\n";
        return 1;
    }
    CharString seq;
    readRegion(seq, faiIndex, idx, 0, 100);

    // Now print the first 100 characters we just read.
    std::cout << "chr:1-100 = " << seq << "\n";

    return 0;
}
