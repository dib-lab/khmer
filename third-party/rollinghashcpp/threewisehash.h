#ifndef THREEWISEHASH
#define THREEWISEHASH

#include <deque>
#include <vector>
#include "characterhash.h"

using namespace std;


/**
* Each instance is a rolling hash function meant to hash streams of characters.
* Each new instance of this class comes with new random keys.
*
* Recommended usage to get L-bit hash values over n-grams:
*        ThreeWiseHash<> hf(n,L );
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
template <typename hashvaluetype = uint32, typename chartype =  unsigned char>
class ThreeWiseHash {

public:

    // myn is the length of the sequences, e.g., 3 means that you want to hash sequences of 3 characters
    // mywordsize is the number of bits you which to receive as hash values, e.g., 19 means that the hash values are 19-bit integers
    ThreeWiseHash(int myn, int mywordsize=19) : n(myn), wordsize(mywordsize),
        hashers(),hasher(0) {
        if(static_cast<uint>(wordsize) > 8*sizeof(hashvaluetype)) {
            cerr<<"Can't create "<<wordsize<<"-bit hash values"<<endl;
            throw "abord";
        }
        for (int i=0; i < n; ++i) {
            CharacterHash<hashvaluetype,chartype> ch(maskfnc<hashvaluetype>(wordsize));
            hashers.push_back(ch);
        }
    }

    // add inchar as an input, this is used typically only at the start
    // the hash value is updated to that of a longer string (one where inchar was appended)
    void eat(chartype inchar) {
        ngram.push_back(inchar);
        __updateHashValue();
    }

    // add inchar as an input and remove outchar, the hashvalue is updated
    // this function can be used to update the hash value from the hash value of [outchar]ABC to the hash value of ABC[inchar]
    void update(chartype outchar, chartype inchar) {
        ngram.push_back(inchar);
        ngram.pop_front();
        __updateHashValue();
    }



    void __updateHashValue() {
        hashvalue = 0;
        for(size_t k = 0; k<ngram.size(); ++k) {
            hashvalue ^= hashers[k].hashvalues[ngram[k]];
        }
    }

    // this is a convenience function, use eat,update and .hashvalue to use as a rolling hash function
    template<class container>
    hashvaluetype hash(container & c) {
        hashvaluetype answer(0);
        for(size_t k = 0; k<c.size(); ++k) {
            answer ^= hashers[k].hashvalues[c[k]];
        }
        return answer;
    }

    hashvaluetype hashvalue;
    int n;
    const int wordsize;
    deque<chartype> ngram;
    vector<CharacterHash<hashvaluetype,chartype> > hashers;
    CharacterHash<hashvaluetype,chartype> hasher;//placeholder

};


#endif

