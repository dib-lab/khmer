#include <iostream>
#include <fstream>
#include <seqan/index.h>
#include <seqan/seq_io.h>

using namespace seqan;


template <typename TIndex>
void findMUMs(TIndex & esa, unsigned minLen)
{
    typename Iterator<TIndex, Mums>::Type it(esa, minLen);  // set minimum MUM length
    String<typename SAValue<TIndex>::Type> occs;            // temp. string storing the hit positions

    std::cout << std::resetiosflags(std::ios::left);
    while (!atEnd(it))
    {
        occs = getOccurrences(it);                          // gives hit positions (seqNo,seqOfs)
        orderOccurrences(occs);                             // order them by seqNo

        for (unsigned i = 0; i < length(occs); ++i)
            std::cout << std::setw(8)
                      << getValueI2(occs[i]) + 1            // output them in MUMmer's output format
                      << "  ";

        std::cout << std::setw(8)
                  << repLength(it)
                  << "\n";

        ++it;
    }
    std::cout << std::setiosflags(std::ios::left) << "\n";
}

template <typename TSpec>
int runMummy(int argc, const char * argv[], unsigned seqCount, unsigned minLen)
{
    typedef String<Dna5, TSpec> TString;
    typedef StringSet<TString>  TStringSet;
    typedef Index<TStringSet>   TIndex;

    TIndex  index;

    // import sequences
    StringSet<CharString> meta;
    for (int arg = 1, seq = 0; arg < argc; ++arg)
    {
        // skip two argument options
        if (strcmp(argv[arg], "-p") == 0 || strcmp(argv[arg], "--profile") == 0 ||
            strcmp(argv[arg], "-l") == 0 || strcmp(argv[arg], "--minlen") == 0)
        {
            ++arg;
            continue;
        }

        if (argv[arg][0] != '-')
        {
            SeqFileIn file;
            if (!open(file, argv[arg]))
            {
                std::cout << "Import of sequence " << argv[arg] << " failed.\n";
                return 1;
            }
            readRecords(meta, indexText(index), file);
            clear(meta);
            close(file);
            ++seq;
        }
    }
    std::cout << lengthSum(indexText(index)) << " bps sequence imported.\n";

    findMUMs(index, minLen);

    return 0;
}

void printHelp(int, const char *[], bool longHelp = false)
{
    std::cout << "***************************************\n";
    std::cout << "***        Simple MUM finder        ***\n";
    std::cout << "*** written by David Weese (c) 2007 ***\n";
    std::cout << "***************************************\n\n";
    std::cout << "Usage: mummy [OPTION]... <SEQUENCE FILE> ... <SEQUENCE FILE>\n";
    if (longHelp)
    {
        std::cout << "\nOptions:\n";
        std::cout << "  -e, --extern          \tuse external memory (for large datasets)\n";
        std::cout << "  -l, --minlen          \tset minimum MUM length\n";
        std::cout << "                        \tif not set, default value is 20\n";
        std::cout << "  -h, --help            \tprint this help\n";
    }
    else
    {
        std::cout << "Try 'mummy --help' for more information.\n";
    }
}

int main(int argc, const char * argv[])
{
    if (argc < 2)
    {
        printHelp(argc, argv);
        return 0;
    }

    unsigned optMinLen = 20;
    bool     optExtern = false;

    unsigned seqCount = 0;
    for (int arg = 1; arg < argc; ++arg)
    {
        if (argv[arg][0] == '-')
        {
            // parse option

            if (strcmp(argv[arg], "-e") == 0 || strcmp(argv[arg], "--extern") == 0)
            {
                // use external memory algorithms
                optExtern = true;
                continue;
            }
            if (strcmp(argv[arg], "-l") == 0 || strcmp(argv[arg], "--minlen") == 0)
            {
                // minimum match length
                if (arg + 1 == argc)
                {
                    printHelp(argc, argv);
                    return 0;
                }
                ++arg;
                optMinLen = atoi(argv[arg]);
                continue;
            }
            if (strcmp(argv[arg], "-h") == 0 || strcmp(argv[arg], "--help") == 0)
            {
                // print help
                printHelp(argc, argv, true);
                return 0;
            }
        }
        else
        {
            // parse sequence file
            ++seqCount;
        }
    }

    if (optExtern)
        return runMummy<External<> >(argc, argv, seqCount, optMinLen);
    else
        return runMummy<Alloc<> >(argc, argv, seqCount, optMinLen);
}
