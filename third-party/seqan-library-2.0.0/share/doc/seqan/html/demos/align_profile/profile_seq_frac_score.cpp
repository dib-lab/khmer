#include <iostream>
#include <seqan/align_profile.h>

using namespace seqan;

int main()
{
    typedef ProfileChar<Dna, int> TDnaProfile;
    typedef String<TDnaProfile> TProfileString;

    TProfileString profile = "CGAT";
    DnaString seq = "CGGAAT";

    Gaps<TProfileString> gapsH(profile);
    Gaps<DnaString> gapsV(seq);

    Score<int, ProfileSeqFracScore> sScheme(profile);

    int val = globalAlignment(gapsH, gapsV, sScheme, NeedlemanWunsch());
    std::cout << "score value = " << val << "\n";

    std::cout << "gaps in profile/sequence\n"
              << "pos\tG\tS\n";
    for (unsigned i = 0; i < length(gapsH); ++i)
        std::cerr << i << "\t" << isGap(gapsH, i) << "\t" << isGap(gapsV, i) << "\n";

    return 0;
}
