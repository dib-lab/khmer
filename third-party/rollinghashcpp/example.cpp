#include <string>
#include <vector>
#include <memory>
#include <iostream>

#include "rabinkarphash.h"

int main(int argc, char * argv[])
{
    size_t q = 3;
    size_t k = 4;
    typedef KarpRabinHash<>  HashFunction;
    std::vector<std::unique_ptr<HashFunction> > hashPtr(q);
    for(size_t z = 0; z < hashPtr.size(); ++z)
    {
        std::unique_ptr<HashFunction> & ptr  = hashPtr[z];
        ptr.reset(new HashFunction(k, 12));
    }

    std::string str = "ACGTAACGT";
    for (size_t j = 0; j < k; j++)
    {
        for(size_t z = 0; z < hashPtr.size(); ++z)
        {
            std::unique_ptr<HashFunction> & ptr  = hashPtr[z];
            ptr->eat(str[j]);
        }

    }

    for (size_t i = 0;; i++)
    {
        std::cout << std::string(str.begin() + i, str.begin() + i + k);
        for(size_t z = 0; z < hashPtr.size(); ++z)
        {
            std::unique_ptr<HashFunction> & ptr  = hashPtr[z];
            std::cout << ' ' << ptr->hashvalue;
        }

        std::cout << std::endl;
        if (i + k < str.size())
        {
            for(size_t z = 0; z < hashPtr.size(); ++z)
            {
                std::unique_ptr<HashFunction> & ptr  = hashPtr[z];
                ptr->update(str[i], str[i + k]);
            }
        }
        else
        {
            break;
        }
    }

    return 0;
}
