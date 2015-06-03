#include <iostream>
#include <seqan/store.h>

using namespace seqan;

int main()
{
    // instantiate emtpy FragmentStore and set file paths
    FragmentStore<> store;
    std::string pathGenome = std::string(SEQAN_PATH_TO_ROOT()) + "/demos/tutorial/store/ex1.fa";
    std::string pathSAM    = std::string(SEQAN_PATH_TO_ROOT()) + "/demos/tutorial/store/ex1.sam";

    // load example genome and example reads and alignments
    loadContigs(store, pathGenome.c_str());
    BamFileIn bamFile(pathSAM.c_str());
    readRecords(store, bamFile);

    // compute staircase read layout and print from position 30..129
    AlignedReadLayout layout;
    layoutAlignment(layout, store);
    printAlignment(std::cout, layout, store, 1, 30, 130, 0, 36);

    return 0;
}
