#include <iostream>
#include <seqan/modifier.h>
#include <seqan/stream.h>

using namespace seqan;

int main()
{
    // A nested modifier.
    typedef ModifiedString<String<Dna>, ModComplementDna>   TMyComplement;
    typedef ModifiedString<TMyComplement, ModReverse>       TMyReverseComplement;

    // The original string.
    String<Dna> myString = "attacgg";
    // A reverse complemented string.
    TMyReverseComplement myReverseComplement(myString);
    std::cout << myString << "\n"
              << myReverseComplement << "\n";

    replace(myString, 1, 1, "cgt");
    std::cout << myString << "\n"
              << myReverseComplement << "\n"
              << DnaStringReverseComplement(myString) << "\n";

    return 0;
}
