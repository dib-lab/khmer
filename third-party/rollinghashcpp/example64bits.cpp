#include <string>
#include <vector>
#include <memory>
#include <iostream>

// Example of 64-bit hashing

#include "cyclichash.h"


int main(int argc, char * argv[])
{
    CyclicHash<uint64>  hf(5,64);
    string input = "ABCDE";
    hf.eat(input[0]);//A
    hf.eat(input[1]);//B
    hf.eat(input[2]);//C
    hf.eat(input[3]);//D
    cout<<"Hash value of ABCD is " << hf.hashvalue << endl;
    // we check the answer going the long way...
    const std::vector<unsigned char> charvectslice(input.begin(), input.begin()+4);
    uint64_t trueanswerslice  = hf.hash(charvectslice);
    if(trueanswerslice != hf.hashvalue ) throw runtime_error("bug");
    // we continue
    hf.eat(input[4]);//E
    cout<<"Hash value of ABCDE is " << hf.hashvalue << endl;
    // we check the answer going the long way
    const std::vector<unsigned char> charvect(input.begin(), input.end());
    uint64_t trueanswer  = hf.hash(charvect);
    if(trueanswer != hf.hashvalue ) throw runtime_error("bug");
    return 0;

}
