#include <iostream>

#include <seqan/consensus.h>
#include <seqan/sequence.h>
#include <seqan/store.h>

using namespace seqan;

int main()
{
    // Reference to simulate reads from.
    Dna5String ref =
        "AATGGATGGCAAAATAGTTGTTCCATGAATACATCTCTAAAGAGCTTTGATGCTAATTTAGTCAAATTTT"
        "CAATACTGTACAATCTTCTCTAGAGCAGAGCAAAAGAATAAAAGCACTTCTAGCTAATATTATGTGGCAT";

    // Read length and step width for reads.
    int const READ_LENGTH = 50;
    int const STEP = 5;

    // Compute reads and append to FragmentStore.
    FragmentStore<> store;
    for (unsigned pos = 0, i = 0; pos + READ_LENGTH < length(ref); pos += STEP, ++i)
    {
        // Append a new read sequence.
        unsigned readID = appendRead(store, infix(ref, pos, pos + READ_LENGTH));
        // Create small perturbation of the position but not left of position 0.
        int pos2 = std::max(0, (int)pos + ((int)i % 6 - 3));
        // Append a new read alignment for the just added read.
        appendAlignedRead(store, readID, 0, pos2, pos2 + READ_LENGTH);
    }

    // Add contig and contig name for printing.
    resize(store.contigStore, 1);
    store.contigStore[0].seq = ref;
    resize(store.contigNameStore, 1);

    // Print initial perturbated alignment.
    std::cout << "Initial Alignment\n\n";
    AlignedReadLayout layout;
    layoutAlignment(layout, store);
    printAlignment(std::cout, layout, store, /*contigID=*/ 0,
                   /*beginPos=*/ 0, /*endPos=*/ (int)length(ref), 0, 30);

    // Compute consensus alignment.
    ConsensusAlignmentOptions options;
    consensusAlignment(store, options);

    // Print final consensus alignment.
    std::cout << "\n\nFinal Consensus Alignment\n\n";
    layoutAlignment(layout, store);
    printAlignment(std::cout, layout, store, /*contigID=*/ 0,
                   /*beginPos=*/ 0, /*endPos=*/ (int)length(ref), 0, 30);

    return 0;
}
