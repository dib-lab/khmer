#include <iostream>
#include <string>

#include <seqan/sequence.h>
#include <seqan/stream.h>

using namespace seqan;

int main()
{
    bool success;

    // interpret string literal as int
    int x = 0;
    // successful cast
    std::cout << "success=" << lexicalCast(x, "123") << ", x=" << x << "\n";
    // unsuccessful cast
    std::cout << "success=" << lexicalCast(x, "42a") << ", x=" << x << "\n";

    // interpret std::string as int
    std::cout << "success=" << lexicalCast(x, std::string("234")) << ", x=" << x << "\n";
    // interpret CharString as int
    std::cout << "success=" << lexicalCast(x, CharString("345")) << ", x=" << x << "\n";
    // interpret infix as int
    CharString str = "123";
    std::cout << "success=" << lexicalCast(x, infix(str, 1, 2)) << ", x=" << x << "\n";

    // interpret literal as float and double
    float f = -1;
    double d = -1;
    std::cout << "success=" << lexicalCast(f, "3.1") << ", f=" << f << "\n"
              << "success=" << lexicalCast(d, "4.2") << ", d=" << d << "\n";

    // demonstrate throwing of exceptions with lexicalCast()
    try
    {
        x = lexicalCast<int>("a");
    }
    catch (BadLexicalCast const & badCast)
    {
        std::cout << "error: " << badCast.what() << "\n";
    }

    // demonstrate throwing of exceptions with lexicalCastWithException
    try
    {
        lexicalCastWithException(x, "b");
    }
    catch (BadLexicalCast const & badCast)
    {
        std::cout << "error: " << badCast.what() << "\n";
    }

    return 0;
}
