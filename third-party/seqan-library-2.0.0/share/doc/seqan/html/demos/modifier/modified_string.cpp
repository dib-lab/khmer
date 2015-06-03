#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/stream.h>  // for I/O
#include <seqan/modifier.h>

using namespace seqan;

int main()
{
    String<Dna> seq = "TATACGCGAAAA";
    ModifiedString<String<Dna>, ModReverse> myModifier(seq);

    std::cout << seq << std::endl;
    std::cout << myModifier << std::endl;
    replace(seq, 8, 12, "TTTT");
    std::cout << std::endl;
    std::cout << seq << std::endl;
    std::cout << myModifier << std::endl;

    return 0;
}
