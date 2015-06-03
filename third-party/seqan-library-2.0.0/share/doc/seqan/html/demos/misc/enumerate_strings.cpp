#include <iostream>

#include <seqan/basic.h>
#include <seqan/sequence.h>
#include <seqan/stream.h>      // For printing SeqAn Strings.
#include <seqan/misc/edit_environment.h>

using namespace seqan;

int main()
{
    Dna5String original = "CGAT";

    // Enumerate neighbourhood using Hamming distance.
    typedef StringEnumerator<Dna5String, EditEnvironment<HammingDistance, 2> > THammingEnumerator;
    typedef Iterator<THammingEnumerator>::Type THammingIterator;
    std::cout << "Enumerating Hamming distance environment of " << original << " of distance 2\n";
    THammingEnumerator hammingEnumerator(original);
    for (THammingIterator itH = begin(hammingEnumerator); !atEnd(itH); goNext(itH))
        std::cout << *itH << '\n';

    // Enumerate neighbourhood using edit distance.
    typedef StringEnumerator<Dna5String, EditEnvironment<LevenshteinDistance, 2> > TEditEnumerator;
    typedef Iterator<TEditEnumerator>::Type TEditIterator;
    std::cout << "\nEnumerating edit distance environment of " << original << " of distance 1-2\n";
    TEditEnumerator editEnumerator(original);
    for (TEditIterator itE = begin(editEnumerator); !atEnd(itE); goNext(itE))
        std::cout << *itE << '\n';

    return 0;
}
