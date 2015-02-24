#include <iostream>

#include <seqan/align.h>
#include <seqan/sequence.h>

using namespace seqan;

int main()
{
    // We will create the following alignment between the subject and the
    // query (with clipping in subject).
    //
    //                    1    1    2    2
    //          0    5    0    5    0    5
    //          :    .    :    .    :    .
    // subject  CPISRTW-SIFRCWALLKAMEAALL
    //             |||| | |||||
    // query       SRTWAS-FRCWA
    //
    //             | clipping |
    //             '----------'

    Peptide subject = "CPISRTWSIFRCWALLKAMEAALL";
    Peptide query   = "SRTWASFRCWA";

    std::cout << "subject: " << subject << "\n"
              << "query:   " << query << "\n\n";

    // Build initial alignment.  Note that the resulting alignment does not
    // make sense yet since the row size is different.
    Align<Peptide> align;
    resize(rows(align), 2);
    setSource(row(align, 0), subject);
    setSource(row(align, 1), query);

    // Build the alignment that is to be integrated into align.
    Infix<Peptide>::Type subjectInfix = infix(subject, 3, 9);
    Infix<Peptide>::Type queryInfix = infix(query, 0, 6);

    std::cout << "subject infix: " << subjectInfix << "\n"
              << "query infix:   " << queryInfix << "\n\n";

    Align<Infix<Peptide>::Type> infixAlign;
    resize(rows(infixAlign), 2);
    setSource(row(infixAlign, 0), subjectInfix);
    setSource(row(infixAlign, 1), queryInfix);

    globalAlignment(infixAlign, Blosum62());

    std::cout << "infix alignment\n"
              << infixAlign;

    // Now integrate infixAlign into align.  Note that the alignment itself
    // does not make sense yet either: the whole query is aligned to a part of
    // the subject and we have to limit/clip the subject row of align.
    integrateAlign(align, infixAlign);

    // For the clipping, we have to transform the alignment begin position in
    // sequence space of subject into the view space of row(align, 0) since
    // there are (and generally can be) gaps in row(align, 0);
    int beginSourcePos = 3;
    int endSourcePos = 14;
    int beginViewPos = toViewPosition(row(align, 0), beginSourcePos);
    int endViewPos = toViewPosition(row(align, 0), endSourcePos);

    std::cout << "clipping row(align, 0) to source range (" << beginSourcePos
              << ", " << endSourcePos << ").  This is the view range ("
              << beginViewPos << ", " << endViewPos << ")\n";
    setClippedBeginPosition(row(align, 0), beginViewPos);
    setClippedEndPosition(row(align, 0), endViewPos);

    std::cout << "align with clipping\n"
              << align;

    return 0;
}
