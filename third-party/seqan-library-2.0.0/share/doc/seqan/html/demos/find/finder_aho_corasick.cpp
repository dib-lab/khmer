#include <seqan/find.h>

using namespace seqan;

int main()
{
    typedef String<AminoAcid> AminoAcidString;

    // A simple amino acid database.
    StringSet<AminoAcidString> dbs;
    appendValue(dbs, "MARDPLY");
    appendValue(dbs, "AVGGGGAAA");
    // We put some words of the database into the queries.
    String<AminoAcidString> queries;
    appendValue(queries, "MARD");
    appendValue(queries, "AAA");
    appendValue(queries, "DPLY");
    appendValue(queries, "VGGGG");

    // Define the Aho-Corasick pattern over the queries with the preprocessing
    // data structure.
    Pattern<String<AminoAcidString>, AhoCorasick> pattern(queries);

    // Search for the queries in the databases.  We have to search database
    // sequence by database sequence.
    std::cout << "DB\tPOS\tENDPOS\tTEXT\n";
    for (unsigned i = 0; i < length(dbs); ++i)
    {
        Finder<AminoAcidString> finder(dbs[i]);  // new finder for each seq
        while (find(finder, pattern))
            std::cout << i << "\t" << position(finder) << "\t"
                      << endPosition(finder) << "\t"
                      << infix(finder) << "\n";
    }

    return 0;
}
