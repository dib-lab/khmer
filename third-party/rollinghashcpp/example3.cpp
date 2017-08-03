#include <string>
#include <vector>
#include <memory>
#include <iostream>


#include "cyclichash.h"

// given hash value of "BCD", can I have value of
// "ABC"quicky?

int demo1() {
    CyclicHash<> hf(3, 32);
    string input = "ABCD";
    hf.eat(input[1]);    //B
    hf.eat(input[2]);    //C
    hf.eat(input[3]);    //D
    cout << "Hash value of BCD is " << hf.hashvalue << endl;
    // we check the answer going the long way...
    const std::vector<unsigned char> charvectslice(input.begin() + 1,
            input.begin() + 4);
    uint32_t trueanswerslice = hf.hash(charvectslice);
    if (trueanswerslice != hf.hashvalue)
        throw runtime_error("bug");
    // we continue
    hf.reverse_update(input[0], input[3]);    //remove D, prepend A
    cout << "Hash value of ABC is " << hf.hashvalue << endl;
    // we check the answer going the long way
    const std::vector<unsigned char> charvect(input.begin(), input.begin() + 3);
    uint32_t trueanswer = hf.hash(charvect);
    if (trueanswer != hf.hashvalue)
        throw runtime_error("bug");
    return 0;
}


int main(int argc, char * argv[])
{
    demo1();
}
