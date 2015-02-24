#include <seqan/modifier.h>
#include <seqan/stream.h>

using namespace seqan;

int main(int argc, char const ** argv)
{
    // create a (generic) CyclicShape
    typedef CyclicShape<GenericShape> TShape;
    TShape shape;
    stringToCyclicShape(shape, "1110010");

    // print cyclic Shape
    CharString out;
    cyclicShapeToString(out, shape);
    std::cout << "shape: " << out << std::endl;

    // determine weight and span
    std::cout << "weight: " << weight(shape);
    std::cout << ", span: " << shape.span << std::endl;

    // modify a text to leave out characters
    CharString str = "this is an original string";
    ModifiedString<CharString, ModCyclicShape<TShape> > modStr(str, shape);

    // modStr can be used like a normal String
    std::cout << str << " => " << modStr << std::endl;
    std::cout << "length: " << length(str) << " => " << length(modStr) << std::endl;

    return 0;
}
