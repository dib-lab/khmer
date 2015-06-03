#include <iostream>
#include <seqan/stream.h>
#include <seqan/modifier.h>

using namespace seqan;

// A user-defined modifier that transforms all characters to upper case.
struct MyFunctor :
    public std::unary_function<char, char>
{
    inline char operator()(char x) const
    {
        if (('a' <= x) && (x <= 'z'))
            return x + ('A' - 'a');

        return x;
    }

};


int main()
{
    // Construct a String and a ModifiedString over it.
    String<char> myString = "A man, a plan, a canal-Panama";
    ModifiedString<String<char>, ModView<MyFunctor> > myModifier(myString);

    // Print the result and demonstrate that the changes to myString are also
    // visible through myModifier.
    std::cout << myString << "\n"
              << myModifier << "\n";

    replace(myString, 9, 9, "master ");

    std::cout << myString << "\n"
              << myModifier << "\n";

    return 0;
}
