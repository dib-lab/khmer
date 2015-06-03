#include <seqan/sequence.h>
#include <seqan/index.h>

using namespace seqan;

int main()
{
    DnaString text = "AAAACACAGTTTGA";
    Shape<Dna, UngappedShape<3> > myShape;

    // loop using hash() and hashNext() starts at position 1
    std::cout << hash(myShape, begin(text)) << '\t';
    for (unsigned i = 1; i < length(text) - length(myShape) + 1; ++i)
        std::cout << hashNext(myShape, begin(text) + i) << '\t';
    std::cout << std::endl;

    // loop using hashInit() and hashNext() starts at position 0
    hashInit(myShape, begin(text));
    for (unsigned i = 0; i < length(text) - length(myShape) + 1; ++i)
        std::cout << hashNext(myShape, begin(text) + i) << '\t';
    std::cout << std::endl;

    return 0;
}
