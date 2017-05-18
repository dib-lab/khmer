#ifndef GENERALHASH
#define GENERALHASH

#include <iostream>
#include <vector>

#include "characterhash.h"

using namespace std;

enum {NOPRECOMP,FULLPRECOMP};

/**
* Each instance is a rolling hash function meant to hash streams of characters.
* Each new instance of this class comes with new random keys.
*
* Recommended usage to get L-bit hash values over n-grams:
*        GeneralHash<> hf(n,L );
*        for(uint32 k = 0; k<n;++k) {
*                  unsigned char c = ... ; // grab some character
*                  hf.eat(c); // feed it to the hasher
*        }
*        while(...) { // go over your string
*           hf.hashvalue; // at all times, this contains the hash value
*           unsigned char c = ... ;// points to the next character
*           unsigned char out = ...; // character we want to forget
*           hf.update(out,c); // update hash value
*        }
*/
template <int precomputationtype=NOPRECOMP,typename hashvaluetype = uint32, typename chartype =  unsigned char>
class GeneralHash {
public:

    // myn is the length of the sequences, e.g., 3 means that you want to hash sequences of 3 characters
    // mywordsize is the number of bits you which to receive as hash values, e.g., 19 means that the hash values are 19-bit integers
    GeneralHash(int myn, int mywordsize = 19):
        hashvalue(0),
        wordsize(mywordsize),
        n(myn),
        irreduciblepoly(0),
        hasher(maskfnc<hashvaluetype>(wordsize)),
        lastbit(static_cast<hashvaluetype>(1)<<wordsize),
        precomputedshift(precomputationtype==FULLPRECOMP ? (1<<n) : 0) {
        if(wordsize == 19) {
            irreduciblepoly = 1 + (1<<1) + (1<<2) + (1<<5)
                               + (1<<19);
        } else if (wordsize == 9) {
            irreduciblepoly = 1+(1<<2)+(1<<3)+(1<<5)+(1<<9);
        } else {
            cerr << "unsupported wordsize "<<wordsize << " bits, try 19 or 9"<< endl;
        }
        // in case the precomp is activated at the template level
        if(precomputationtype==FULLPRECOMP) {
            for(hashvaluetype x = 0; x<precomputedshift.size(); ++x) {
                hashvaluetype leftover = x << (wordsize-n);
                fastleftshift(leftover, n);
                precomputedshift[x]=leftover;
            }
        }
    }


    void fastleftshift(hashvaluetype & x, int r) const {
        for (int i = 0; i < r; ++i) {
            x  <<= 1;
            if(( x & lastbit) == lastbit)
                x ^= irreduciblepoly;
        }
    }

    void fastleftshiftn(hashvaluetype & x) const {
        x=
            // take the last n bits and look-up the result
            precomputedshift[(x >> (wordsize-n))]
            ^
            // then just shift the first L-n bits
            ((x << n) & (lastbit -1 ));
    }

    // add inchar as an input and remove outchar, the hashvalue is updated
    // this function can be used to update the hash value from the hash value of [outchar]ABC to the hash value of ABC[inchar]
    void update(chartype outchar, chartype inchar) {
        hashvalue <<= 1;
        if(( hashvalue & lastbit) == lastbit)
            hashvalue ^= irreduciblepoly;
        //
        hashvaluetype z (hasher.hashvalues[outchar]);
        // the compiler should optimize away the next if/else
        if(precomputationtype==FULLPRECOMP) {
            fastleftshiftn(z);
            hashvalue ^= z ^ hasher.hashvalues[inchar];
        } else {
            fastleftshift(z,n);
            hashvalue ^= z ^ hasher.hashvalues[inchar];
        }
    }



    // add inchar as an input, this is used typically only at the start
    // the hash value is updated to that of a longer string (one where inchar was appended)
    void eat(chartype inchar) {
        fastleftshift(hashvalue,1);
        hashvalue ^=  hasher.hashvalues[inchar];
    }

    // this is a convenience function, use eat,update and .hashvalue to use as a rolling hash function
    template<class container>
    hashvaluetype  hash(container & c) const {
        hashvaluetype answer(0);
        for(uint k = 0; k<c.size(); ++k) {
            fastleftshift(answer, 1) ;
            answer ^= hasher.hashvalues[c[k]];
        }
        return answer;
    }

    hashvaluetype hashvalue;
    const int wordsize;
    int n;
    hashvaluetype irreduciblepoly;
    CharacterHash<hashvaluetype,chartype> hasher;
    const hashvaluetype lastbit;
    vector<hashvaluetype> precomputedshift;

};



#endif

