#ifndef CHARACTERHASH
#define CHARACTERHASH

#include <cassert>
#include <iostream>
#include <stdexcept>
#include <random>

using namespace std;


class mersenneRNG {
public:
    mersenneRNG(uint32_t maxval, uint32_t seed) : dist(0, maxval), n(maxval) {
        mtr = std::mt19937(seed);
    }
    uint32_t operator()() {
        return dist(mtr);
    }
    void seed(uint32_t seedval) {
        mtr.seed(seedval);
    }
    void seed() {
        mtr.seed();
    }
    uint32_t rand_max() {
        return n;
    }
private:
    std::mt19937 mtr;
    std::uniform_int_distribution<uint32_t> dist;
    uint32_t n;
};

template <typename hashvaluetype>
hashvaluetype maskfnc(int bits) {
    assert(bits>0);
    assert(bits<=sizeof(hashvaluetype)*8);
    hashvaluetype x = static_cast<hashvaluetype>(1) << (bits - 1);
    return x ^ (x - 1);
}

template <typename hashvaluetype = uint32_t, typename chartype =  unsigned char>
class CharacterHash {
public:
    CharacterHash(hashvaluetype maxval, hashvaluetype seed) {
        if(sizeof(hashvaluetype) <=4) {
            mersenneRNG randomgenerator(maxval, seed);
            for(size_t k =0; k<nbrofchars; ++k)
                hashvalues[k] = static_cast<hashvaluetype>(randomgenerator());
        } else if (sizeof(hashvaluetype) == 8) {
            mersenneRNG randomgenerator(maxval>>32, seed);
            mersenneRNG randomgeneratorbase(
              (maxval>>32) ==0 ? maxval : 0xFFFFFFFFU, seed+1);
            for(size_t k =0; k<nbrofchars; ++k)
                hashvalues[k] = static_cast<hashvaluetype>(randomgeneratorbase())
                                | (static_cast<hashvaluetype>(randomgenerator()) << 32);
        } else throw runtime_error("unsupported hash value type");
    }

    enum {nbrofchars = 1 << ( sizeof(chartype)*8 )};

    hashvaluetype hashvalues[1 << ( sizeof(chartype)*8 )];
};

#endif
