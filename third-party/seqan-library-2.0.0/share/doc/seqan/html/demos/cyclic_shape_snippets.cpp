// ==========================================================================
// Author: Sascha Meiers <meiers@inf.fu-berlin.de>
// ==========================================================================
// Snippets for CyclicShape demonstrations
// ==========================================================================


#include <seqan/sequence.h>
#include <seqan/stream.h>      // For printing SeqAn Strings.
#include <seqan/modifier.h>

using namespace seqan;
int main(int argc, char const ** argv)
{

    {
        //![Define FixedCyclicShape]
        typedef GappedShape<HardwiredShape<1, 1, 3, 2> >       TInnerShape; // 11100101
        // ^--You can also use predefied Shapes here, e.g. Patternhunter

        typedef CyclicShape<FixedShape<2, TInnerShape, 1> > TShape;     // 00111001010
        TShape shape;
        //![Define FixedCyclicShape]


        //![Define FixedCyclicShape Modified String]
        typedef ModifiedString<CharString, ModCyclicShape<TShape> > TModString;

        CharString str = "das ist das Haus vom Nikolaus";
        TModString modStr(str);

        std::cout << str << " => " << modStr << std::endl;
        //![Define FixedCyclicShape Modified String]
    }


    {
        //![Define GenericCyclicShape]
        typedef CyclicShape<GenericShape> TShape;
        TShape shape;
        stringToCyclicShape(shape, "00111001010");
        //![Define GenericCyclicShape]


        //![Define GenericCyclicShape Modified String]
        typedef ModifiedString<CharString, ModCyclicShape<TShape> > TModString;

        CharString str = "das ist das Haus vom Nikolaus";
        TModString modStr(str, shape);

        std::cout << str << " => " << modStr << std::endl;
        //![Define GenericCyclicShape Modified String]


        //![Define CyclicShape Modified Iterator]
        typedef ModifiedString<CharString, ModCyclicShape<TShape> > TModString;
        typedef Iterator<TModString>::Type TModIter;
        TModIter it = begin(modStr);
        TModIter itBeg = it;
        TModIter itEnd = end(modStr);

        // output iter position and host iter position
        for (; it != itEnd; ++it)
            std::cout << (it - itBeg) << "/" << (host(it) - begin(str)) << ", ";

        //  0/2, 1/3, 2/4, 3/7, 4/9, 5/13, 6/14, 7/15, 8/18, 9/20, 10/24, 11/25, 12/26,
        //![Define CyclicShape Modified Iterator]


        //![CyclicShape Care Positions]
        std::cout << std::endl << "relative care positions: ";
        for (unsigned i = 0; i < weight(shape); ++i)
            std::cout << (int)shape.diffs[i] << ",";    // output: 1,1,3,2,4,

        std::cout << std::endl << "absolute care positions: ";
        String<int> carePos;
        carePositions(carePos, shape);

        for (unsigned i = 0; i < weight(shape); ++i)
            std::cout << carePos[i] << ",";             // output: 2,3,4,7,9,
        std::cout << std::endl;
        //![CyclicShape Care Positions]
    }

    return 0;
}
