#include <iostream>

#include <seqan/misc/terminal.h>

using namespace seqan;

int main()
{
    // Get terminal size and print it to stdout.
    unsigned cols = 0, rows = 0;
    getTerminalSize(cols, rows);
    std::cout << "cols == " << cols << ", rows == " << rows << "\n";

    return 0;
}
