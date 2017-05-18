#ifndef CYCLICHASH
#define CYCLICHASH

#include "characterhash.h"

/**
* Each instance is a rolling hash function meant to hash streams of characters.
* Each new instance of this class comes with new random keys.
*
* Recommended usage to get L-bit hash values over n-grams:
*        CyclicHash<> hf(n,L );
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
class CyclicHash {

public:
    // myn is the length of the sequences, e.g., 3 means that you want to hash sequences of 3 characters
    // mywordsize is the number of bits you which to receive as hash values, e.g., 19 means that the hash values are 19-bit integers
    CyclicHash(int myn, int mywordsize=19) : hashvalue(0),
        n(myn), wordsize(mywordsize),
        hasher(maskfnc<hashvaluetype>(wordsize)),
        mask1(maskfnc<hashvaluetype>(wordsize-1)),
        myr(n%wordsize),
        maskn(maskfnc<hashvaluetype>(wordsize-myr))
    {
        if(static_cast<uint>(wordsize) > 8*sizeof(hashvaluetype)) {
            cerr<<"Can't create "<<wordsize<<"-bit hash values"<<endl;
            throw "abord";
        }
    }

    void fastleftshiftn(hashvaluetype & x) const {
        x =  ((x & maskn) << myr ) | (x >> (wordsize-myr)) ;
    }

    void fastleftshift1(hashvaluetype & x) const {
        x =  ((x & mask1) << 1 ) | (x >> (wordsize-1)) ;
    }

    void fastrightshift1(hashvaluetype & x) const {
        x =  (x  >> 1 ) | ((x & 1)<< (wordsize-1)) ;
    }


    hashvaluetype getfastleftshift1(hashvaluetype  x) const {
        return ((x & mask1) << 1 ) | (x >> (wordsize-1)) ;
    }


    hashvaluetype getfastrightshift1(hashvaluetype  x) const {
        return (x  >> 1 ) | ((x & 1)<< (wordsize-1)) ;
    }

    // this is a convenience function, use eat,update and .hashvalue to use as a rolling hash function
    template<class container>
    hashvaluetype  hash(container & c) {
        hashvaluetype answer(0);
        for(uint k = 0; k<c.size(); ++k) {
            fastleftshift1(answer);
            answer ^= hasher.hashvalues[static_cast<unsigned int>(c[k])];
        }
        return answer;
    }

    hashvaluetype  hashz(chartype outchar,uint n) {
        hashvaluetype answer = hasher.hashvalues[static_cast<unsigned int>(outchar)];
        for(uint k = 0; k<n; ++k) {
            fastleftshift1(answer);
        }
        return answer;
    }


    // add inchar as an input and remove outchar, the hashvalue is updated
    // this function can be used to update the hash value from the hash value of [outchar]ABC to the hash value of ABC[inchar]
    void update(chartype outchar, chartype inchar) {
        hashvaluetype z (hasher.hashvalues[outchar]);
        fastleftshiftn(z);
        hashvalue =  getfastleftshift1(hashvalue)
                     ^ z
                     ^ hasher.hashvalues[inchar];

    }

    // this is the reverse of the update function.
    // this function can be used to update the hash value from the hash value of ABC[inchar] to the hash value of [outchar]ABC
    void reverse_update(chartype outchar, chartype inchar) {
        hashvaluetype z (hasher.hashvalues[outchar]);
        fastleftshiftn(z);
        hashvalue ^=   z ^ hasher.hashvalues[inchar];
        hashvalue = getfastrightshift1(hashvalue);
    }

    // add inchar as an input, this is used typically only at the start
    // the hash value is updated to that of a longer string (one where inchar was appended)
    void eat(chartype inchar) {
        fastleftshift1(hashvalue);
        hashvalue ^= hasher.hashvalues[inchar];
    }

    //for an n-gram X it returns hash value of (n + 1)-gram XY without changing the object X. For example, if X = "ABC", then X.hash_extend("D") returns value of "ABCD" without changing the state of X
    hashvaluetype hash_extend(chartype Y) {
        return getfastleftshift1(hashvalue) ^ hasher.hashvalues[Y];
    }

    //  same as hash_extend, but with prepending the n-gram with character Y. If X = "ABC", then X.hash_prepend("D") returns value of "DABC" without changing the state of X
    hashvaluetype hash_prepend(chartype Y) {
        hashvaluetype z (hasher.hashvalues[Y]);
        fastleftshiftn(z);
        return z ^ hashvalue;
    }


    hashvaluetype hashvalue;
    int n;
    const int wordsize;
    CharacterHash<hashvaluetype,chartype> hasher;
    const hashvaluetype mask1;
    const int myr;
    const hashvaluetype maskn;

};



#endif
