#include <string>
#include <vector>
#include <memory>
#include <iostream>

// given hash value of "ABCD", can I have value of
// "ABCDE", without computing the whole hash value?

#include "cyclichash.h"


int main(int argc, char * argv[])
{
    CyclicHash<>  hf(5,19);
    string input = "ABCDE";
    hf.eat(input[0]);//A
    hf.eat(input[1]);//B
    hf.eat(input[2]);//C
    hf.eat(input[3]);//D
    cout<<"Hash value of ABCD is " << hf.hashvalue << endl;
    // we check the answer going the long way...
    const std::vector<unsigned char> charvectslice(input.begin(), input.begin()+4);
    uint32_t trueanswerslice  = hf.hash(charvectslice);
    if(trueanswerslice != hf.hashvalue ) throw runtime_error("bug");
    // we continue
    hf.eat(input[4]);//E
    cout<<"Hash value of ABCDE is " << hf.hashvalue << endl;
    // we check the answer going the long way
    const std::vector<unsigned char> charvect(input.begin(), input.end());
    uint32_t trueanswer  = hf.hash(charvect);
    if(trueanswer != hf.hashvalue ) throw runtime_error("bug");
    return 0;

}
