#include <fstream>

#include <seqan/basic.h>
#include <seqan/index.h>
#include <seqan/seq_io.h>

using namespace seqan;

int main()
{
    // Get path to file to search for repeats in.
    std::string path = (std::string)SEQAN_PATH_TO_ROOT() + "/demos/index/ref.fa";

    // Load first sequence from file.
    CharString id;
    Dna5String seq;
    SeqFileIn file(path.c_str());
    readRecord(id, seq, file);

    // Find repeats and print them.
    String<Repeat<unsigned, unsigned> > repeats;
    findRepeats(repeats, seq, 3);

    std::cerr << "# of repeats: " << length(repeats) << "\n";
    for (unsigned i = 0; i < length(repeats); ++i)
        std::cerr << "i == " << i << ", beginPosition = " << repeats[i].beginPosition
                  << ", endPosition = " << repeats[i].endPosition
                  << ", period = " << repeats[i].period << "\n";

    return 0;
}
