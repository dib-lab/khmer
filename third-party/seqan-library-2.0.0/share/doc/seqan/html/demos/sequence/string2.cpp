#include <seqan/sequence.h>
#include <seqan/stream.h>

using namespace seqan;
int main()
{
    // Creating text
    String<char> text = "to be";
    std::cout << text << std::endl;
    appendValue(text, ' ');
    std::cout << "Last sign is whitespace? " << endsWith(text, " ") << std::endl;
    // Erasing whitespaces in text
    eraseBack(text);
    erase(text, 2);
    // Appending another string
    append(text, "ornottobe");
    std::cout << text << std::endl;

    // Pattern
    String<char> pattern = "be";
    // Number of consecutive matching characters per position
    String<int> score;
    resize(score, length(text) - length(pattern) + 1);

    // Brute force pattern matching for every position
    for (unsigned i = 0; i < length(text) - length(pattern) + 1; ++i)
    {
        int localScore = 0;
        for (unsigned j = 0; j < length(pattern); ++j)
            if (text[i + j] == pattern[j])
                ++localScore;
        score[i] = localScore;
    }

    std::cout << "hit at ";
    for (unsigned i = 0; i < length(score); ++i)
        if (score[i] == (int)length(pattern))
            std::cout << i << " ";
    std::cout << std::endl;

    return 0;
}
