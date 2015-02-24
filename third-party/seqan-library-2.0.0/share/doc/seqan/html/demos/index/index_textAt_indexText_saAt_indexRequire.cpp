#include <seqan/index.h>

using namespace seqan;

int main()
{
    Index<String<char> > index("MISSISSIPPI");

    // Because indices are build on demand we force the index creation here.
    indexRequire(index, FibreSA());

    std::cout << "BWT\tSuffices" << std::endl;

    for (unsigned i = 0; i < length(indexSA(index)); ++i)
    {
        unsigned textPos = (saAt(i, index) == 0) ? length(index) - 1 : saAt(i, index) - 1;
        std::cout << textAt(textPos, index) << "\t" << suffix(indexText(index), textPos) << std::endl;
    }
    return 0;
}
