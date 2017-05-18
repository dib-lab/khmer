/*
 *
 * nthash.hpp
 * Author: Hamid Mohamadi
 * Genome Sciences Centre,
 * British Columbia Cancer Agency
 */

#ifndef NT_HASH_H
#define NT_HASH_H

#include <stdint.h>

// offset for the complement base in the random seeds table
const uint8_t cpOff = 0x07;

// shift for gerenerating multiple hash values
const int multiShift = 27;

// seed for gerenerating multiple hash values
static const uint64_t multiSeed = 0x90b45d39fb6da1fa;

// 64-bit random seeds corresponding to bases and their complements
static const uint64_t seedA = 0x3c8bfbb395c60474;
static const uint64_t seedC = 0x3193c18562a02b4c;
static const uint64_t seedG = 0x20323ed082572324;
static const uint64_t seedT = 0x295549f54be24456;
static const uint64_t seedN = 0x0000000000000000;

static const uint64_t vecA[64] = {
    0x3c8bfbb395c60474,0x7917f7672b8c08e8,0xf22feece571811d0,0xe45fdd9cae3023a1,0xc8bfbb395c604743,0x917f7672b8c08e87,0x22feece571811d0f,0x45fdd9cae3023a1e,
    0x8bfbb395c604743c,0x17f7672b8c08e879,0x2feece571811d0f2,0x5fdd9cae3023a1e4,0xbfbb395c604743c8,0x7f7672b8c08e8791,0xfeece571811d0f22,0xfdd9cae3023a1e45,
    0xfbb395c604743c8b,0xf7672b8c08e87917,0xeece571811d0f22f,0xdd9cae3023a1e45f,0xbb395c604743c8bf,0x7672b8c08e87917f,0xece571811d0f22fe,0xd9cae3023a1e45fd,
    0xb395c604743c8bfb,0x672b8c08e87917f7,0xce571811d0f22fee,0x9cae3023a1e45fdd,0x395c604743c8bfbb,0x72b8c08e87917f76,0xe571811d0f22feec,0xcae3023a1e45fdd9,
    0x95c604743c8bfbb3,0x2b8c08e87917f767,0x571811d0f22feece,0xae3023a1e45fdd9c,0x5c604743c8bfbb39,0xb8c08e87917f7672,0x71811d0f22feece5,0xe3023a1e45fdd9ca,
    0xc604743c8bfbb395,0x8c08e87917f7672b,0x1811d0f22feece57,0x3023a1e45fdd9cae,0x604743c8bfbb395c,0xc08e87917f7672b8,0x811d0f22feece571,0x023a1e45fdd9cae3,
    0x04743c8bfbb395c6,0x08e87917f7672b8c,0x11d0f22feece5718,0x23a1e45fdd9cae30,0x4743c8bfbb395c60,0x8e87917f7672b8c0,0x1d0f22feece57181,0x3a1e45fdd9cae302,
    0x743c8bfbb395c604,0xe87917f7672b8c08,0xd0f22feece571811,0xa1e45fdd9cae3023,0x43c8bfbb395c6047,0x87917f7672b8c08e,0x0f22feece571811d,0x1e45fdd9cae3023a
};

static const uint64_t vecC[64] = {
    0x3193c18562a02b4c,0x6327830ac5405698,0xc64f06158a80ad30,0x8c9e0c2b15015a61,0x193c18562a02b4c3,0x327830ac54056986,0x64f06158a80ad30c,0xc9e0c2b15015a618,
    0x93c18562a02b4c31,0x27830ac540569863,0x4f06158a80ad30c6,0x9e0c2b15015a618c,0x3c18562a02b4c319,0x7830ac5405698632,0xf06158a80ad30c64,0xe0c2b15015a618c9,
    0xc18562a02b4c3193,0x830ac54056986327,0x06158a80ad30c64f,0x0c2b15015a618c9e,0x18562a02b4c3193c,0x30ac540569863278,0x6158a80ad30c64f0,0xc2b15015a618c9e0,
    0x8562a02b4c3193c1,0x0ac5405698632783,0x158a80ad30c64f06,0x2b15015a618c9e0c,0x562a02b4c3193c18,0xac54056986327830,0x58a80ad30c64f061,0xb15015a618c9e0c2,
    0x62a02b4c3193c185,0xc54056986327830a,0x8a80ad30c64f0615,0x15015a618c9e0c2b,0x2a02b4c3193c1856,0x54056986327830ac,0xa80ad30c64f06158,0x5015a618c9e0c2b1,
    0xa02b4c3193c18562,0x4056986327830ac5,0x80ad30c64f06158a,0x015a618c9e0c2b15,0x02b4c3193c18562a,0x056986327830ac54,0x0ad30c64f06158a8,0x15a618c9e0c2b150,
    0x2b4c3193c18562a0,0x56986327830ac540,0xad30c64f06158a80,0x5a618c9e0c2b1501,0xb4c3193c18562a02,0x6986327830ac5405,0xd30c64f06158a80a,0xa618c9e0c2b15015,
    0x4c3193c18562a02b,0x986327830ac54056,0x30c64f06158a80ad,0x618c9e0c2b15015a,0xc3193c18562a02b4,0x86327830ac540569,0x0c64f06158a80ad3,0x18c9e0c2b15015a6
};

static const uint64_t vecG[64] = {
    0x20323ed082572324,0x40647da104ae4648,0x80c8fb42095c8c90,0x0191f68412b91921,0x0323ed0825723242,0x0647da104ae46484,0x0c8fb42095c8c908,0x191f68412b919210,
    0x323ed08257232420,0x647da104ae464840,0xc8fb42095c8c9080,0x91f68412b9192101,0x23ed082572324203,0x47da104ae4648406,0x8fb42095c8c9080c,0x1f68412b91921019,
    0x3ed0825723242032,0x7da104ae46484064,0xfb42095c8c9080c8,0xf68412b919210191,0xed08257232420323,0xda104ae464840647,0xb42095c8c9080c8f,0x68412b919210191f,
    0xd08257232420323e,0xa104ae464840647d,0x42095c8c9080c8fb,0x8412b919210191f6,0x08257232420323ed,0x104ae464840647da,0x2095c8c9080c8fb4,0x412b919210191f68,
    0x8257232420323ed0,0x04ae464840647da1,0x095c8c9080c8fb42,0x12b919210191f684,0x257232420323ed08,0x4ae464840647da10,0x95c8c9080c8fb420,0x2b919210191f6841,
    0x57232420323ed082,0xae464840647da104,0x5c8c9080c8fb4209,0xb919210191f68412,0x7232420323ed0825,0xe464840647da104a,0xc8c9080c8fb42095,0x919210191f68412b,
    0x232420323ed08257,0x464840647da104ae,0x8c9080c8fb42095c,0x19210191f68412b9,0x32420323ed082572,0x64840647da104ae4,0xc9080c8fb42095c8,0x9210191f68412b91,
    0x2420323ed0825723,0x4840647da104ae46,0x9080c8fb42095c8c,0x210191f68412b919,0x420323ed08257232,0x840647da104ae464,0x080c8fb42095c8c9,0x10191f68412b9192
};

static const uint64_t vecT[64] = {
    0x295549f54be24456,0x52aa93ea97c488ac,0xa55527d52f891158,0x4aaa4faa5f1222b1,0x95549f54be244562,0x2aa93ea97c488ac5,0x55527d52f891158a,0xaaa4faa5f1222b14,
    0x5549f54be2445629,0xaa93ea97c488ac52,0x5527d52f891158a5,0xaa4faa5f1222b14a,0x549f54be24456295,0xa93ea97c488ac52a,0x527d52f891158a55,0xa4faa5f1222b14aa,
    0x49f54be244562955,0x93ea97c488ac52aa,0x27d52f891158a555,0x4faa5f1222b14aaa,0x9f54be2445629554,0x3ea97c488ac52aa9,0x7d52f891158a5552,0xfaa5f1222b14aaa4,
    0xf54be24456295549,0xea97c488ac52aa93,0xd52f891158a55527,0xaa5f1222b14aaa4f,0x54be24456295549f,0xa97c488ac52aa93e,0x52f891158a55527d,0xa5f1222b14aaa4fa,
    0x4be24456295549f5,0x97c488ac52aa93ea,0x2f891158a55527d5,0x5f1222b14aaa4faa,0xbe24456295549f54,0x7c488ac52aa93ea9,0xf891158a55527d52,0xf1222b14aaa4faa5,
    0xe24456295549f54b,0xc488ac52aa93ea97,0x891158a55527d52f,0x1222b14aaa4faa5f,0x24456295549f54be,0x488ac52aa93ea97c,0x91158a55527d52f8,0x222b14aaa4faa5f1,
    0x4456295549f54be2,0x88ac52aa93ea97c4,0x1158a55527d52f89,0x22b14aaa4faa5f12,0x456295549f54be24,0x8ac52aa93ea97c48,0x158a55527d52f891,0x2b14aaa4faa5f122,
    0x56295549f54be244,0xac52aa93ea97c488,0x58a55527d52f8911,0xb14aaa4faa5f1222,0x6295549f54be2445,0xc52aa93ea97c488a,0x8a55527d52f89115,0x14aaa4faa5f1222b
};

static const uint64_t vecN[64] = {
    seedN,seedN,seedN,seedN,seedN,seedN,seedN,seedN,
    seedN,seedN,seedN,seedN,seedN,seedN,seedN,seedN,
    seedN,seedN,seedN,seedN,seedN,seedN,seedN,seedN,
    seedN,seedN,seedN,seedN,seedN,seedN,seedN,seedN,
    seedN,seedN,seedN,seedN,seedN,seedN,seedN,seedN,
    seedN,seedN,seedN,seedN,seedN,seedN,seedN,seedN,
    seedN,seedN,seedN,seedN,seedN,seedN,seedN,seedN,
    seedN,seedN,seedN,seedN,seedN,seedN,seedN,seedN
};

static const uint64_t *msTab[256] = {
    vecN, vecT, vecN, vecG, vecA, vecN, vecN, vecC, // 0..7
    vecN, vecN, vecN, vecN, vecN, vecN, vecN, vecN, // 8..15
    vecN, vecN, vecN, vecN, vecN, vecN, vecN, vecN, // 16..23
    vecN, vecN, vecN, vecN, vecN, vecN, vecN, vecN, // 24..31
    vecN, vecN, vecN, vecN, vecN, vecN, vecN, vecN, // 32..39
    vecN, vecN, vecN, vecN, vecN, vecN, vecN, vecN, // 40..47
    vecN, vecN, vecN, vecN, vecN, vecN, vecN, vecN, // 48..55
    vecN, vecN, vecN, vecN, vecN, vecN, vecN, vecN, // 56..63
    vecN, vecA, vecN, vecC, vecN, vecN, vecN, vecG, // 64..71
    vecN, vecN, vecN, vecN, vecN, vecN, vecN, vecN, // 72..79
    vecN, vecN, vecN, vecN, vecT, vecN, vecN, vecN, // 80..87
    vecN, vecN, vecN, vecN, vecN, vecN, vecN, vecN, // 88..95
    vecN, vecA, vecN, vecC, vecN, vecN, vecN, vecG, // 96..103
    vecN, vecN, vecN, vecN, vecN, vecN, vecN, vecN, // 104..111
    vecN, vecN, vecN, vecN, vecT, vecN, vecN, vecN, // 112..119
    vecN, vecN, vecN, vecN, vecN, vecN, vecN, vecN, // 120..127
    vecN, vecN, vecN, vecN, vecN, vecN, vecN, vecN, // 128..135
    vecN, vecN, vecN, vecN, vecN, vecN, vecN, vecN, // 136..143
    vecN, vecN, vecN, vecN, vecN, vecN, vecN, vecN, // 144..151
    vecN, vecN, vecN, vecN, vecN, vecN, vecN, vecN, // 152..159
    vecN, vecN, vecN, vecN, vecN, vecN, vecN, vecN, // 160..167
    vecN, vecN, vecN, vecN, vecN, vecN, vecN, vecN, // 168..175
    vecN, vecN, vecN, vecN, vecN, vecN, vecN, vecN, // 176..183
    vecN, vecN, vecN, vecN, vecN, vecN, vecN, vecN, // 184..191
    vecN, vecN, vecN, vecN, vecN, vecN, vecN, vecN, // 192..199
    vecN, vecN, vecN, vecN, vecN, vecN, vecN, vecN, // 200..207
    vecN, vecN, vecN, vecN, vecN, vecN, vecN, vecN, // 208..215
    vecN, vecN, vecN, vecN, vecN, vecN, vecN, vecN, // 216..223
    vecN, vecN, vecN, vecN, vecN, vecN, vecN, vecN, // 224..231
    vecN, vecN, vecN, vecN, vecN, vecN, vecN, vecN, // 232..239
    vecN, vecN, vecN, vecN, vecN, vecN, vecN, vecN, // 240..247
    vecN, vecN, vecN, vecN, vecN, vecN, vecN, vecN  // 248..255
};

static const uint64_t seedTab[256] = {
    seedN, seedT, seedN, seedG, seedA, seedN, seedN, seedC, // 0..7
    seedN, seedN, seedN, seedN, seedN, seedN, seedN, seedN, // 8..15
    seedN, seedN, seedN, seedN, seedN, seedN, seedN, seedN, // 16..23
    seedN, seedN, seedN, seedN, seedN, seedN, seedN, seedN, // 24..31
    seedN, seedN, seedN, seedN, seedN, seedN, seedN, seedN, // 32..39
    seedN, seedN, seedN, seedN, seedN, seedN, seedN, seedN, // 40..47
    seedN, seedN, seedN, seedN, seedN, seedN, seedN, seedN, // 48..55
    seedN, seedN, seedN, seedN, seedN, seedN, seedN, seedN, // 56..63
    seedN, seedA, seedN, seedC, seedN, seedN, seedN, seedG, // 64..71
    seedN, seedN, seedN, seedN, seedN, seedN, seedN, seedN, // 72..79
    seedN, seedN, seedN, seedN, seedT, seedN, seedN, seedN, // 80..87
    seedN, seedN, seedN, seedN, seedN, seedN, seedN, seedN, // 88..95
    seedN, seedA, seedN, seedC, seedN, seedN, seedN, seedG, // 96..103
    seedN, seedN, seedN, seedN, seedN, seedN, seedN, seedN, // 104..111
    seedN, seedN, seedN, seedN, seedT, seedN, seedN, seedN, // 112..119
    seedN, seedN, seedN, seedN, seedN, seedN, seedN, seedN, // 120..127
    seedN, seedN, seedN, seedN, seedN, seedN, seedN, seedN, // 128..135
    seedN, seedN, seedN, seedN, seedN, seedN, seedN, seedN, // 136..143
    seedN, seedN, seedN, seedN, seedN, seedN, seedN, seedN, // 144..151
    seedN, seedN, seedN, seedN, seedN, seedN, seedN, seedN, // 152..159
    seedN, seedN, seedN, seedN, seedN, seedN, seedN, seedN, // 160..167
    seedN, seedN, seedN, seedN, seedN, seedN, seedN, seedN, // 168..175
    seedN, seedN, seedN, seedN, seedN, seedN, seedN, seedN, // 176..183
    seedN, seedN, seedN, seedN, seedN, seedN, seedN, seedN, // 184..191
    seedN, seedN, seedN, seedN, seedN, seedN, seedN, seedN, // 192..199
    seedN, seedN, seedN, seedN, seedN, seedN, seedN, seedN, // 200..207
    seedN, seedN, seedN, seedN, seedN, seedN, seedN, seedN, // 208..215
    seedN, seedN, seedN, seedN, seedN, seedN, seedN, seedN, // 216..223
    seedN, seedN, seedN, seedN, seedN, seedN, seedN, seedN, // 224..231
    seedN, seedN, seedN, seedN, seedN, seedN, seedN, seedN, // 232..239
    seedN, seedN, seedN, seedN, seedN, seedN, seedN, seedN, // 240..247
    seedN, seedN, seedN, seedN, seedN, seedN, seedN, seedN  // 248..255
};

/*// assembly rol
inline uint64_t rol (uint64_t v, size_t n) {
    asm("rol %b1, %0":"+r,r"(v):"i,c"(n));
    return (v);
}

// assembly ror
inline uint64_t ror (uint64_t v, size_t n) {
    asm("ror %b1, %0":"+r,r"(v):"i,c"(n));
    return (v);
}*/

// rotate "v" to the left by "s" positions
inline uint64_t rol(const uint64_t v, const int s) {
    return (v << s) | (v >> (64 - s));
}

// rotate "v" to the right by "s" positions
inline uint64_t ror(const uint64_t v, const int s) {
    return (v >> s) | (v << (64 - s));
}

// forward-strand hash value of the base kmer, i.e. fhval(kmer_0)
inline uint64_t getFhval(const char * kmerSeq, const unsigned k) {
    uint64_t hVal=0;
    for(unsigned i=0; i<k; i++)
        hVal ^= rol(seedTab[(unsigned char)kmerSeq[i]], k-1-i);
    return hVal;
}

// reverse-strand hash value of the base kmer, i.e. rhval(kmer_0)
inline uint64_t getRhval(const char * kmerSeq, const unsigned k) {
    uint64_t hVal=0;
    for(unsigned i=0; i<k; i++)
        hVal ^= rol(seedTab[(unsigned char)kmerSeq[i]&cpOff], i);
    return hVal;
}

// ntHash basic function, i.e. ntBase, using rotate ops
inline uint64_t NT64(const char * kmerSeq, const unsigned k) {
    return getFhval(kmerSeq, k);
}

// ntHash for sliding k-mers, using rotate ops
inline uint64_t NT64(const uint64_t fhVal, const unsigned char charOut, const unsigned char charIn, const unsigned k) {
    return(rol(fhVal, 1) ^ rol(seedTab[charOut], k) ^ seedTab[charIn]);
}

// ntHash with seeding option, using rotate ops
inline uint64_t NT64(const char * kmerSeq, const unsigned k, const unsigned seed) {
    uint64_t hVal=0;
    for(unsigned i=0; i<k; i++)
        hVal ^= rol(seedTab[(unsigned char)kmerSeq[i]], k-1-i);
    hVal *= seed ^ k * multiSeed;
    hVal ^= hVal >> multiShift;
    return hVal;
}

// canonical ntBase, using rotate ops
inline uint64_t NTC64(const char * kmerSeq, const unsigned k) {
    uint64_t fhVal=0, rhVal=0;
    fhVal=getFhval(kmerSeq, k);
    rhVal=getRhval(kmerSeq, k);
    return (rhVal<fhVal)? rhVal : fhVal;
}

// canonical ntHash, using rotate ops
inline uint64_t NTC64(const char * kmerSeq, const unsigned k, uint64_t& fhVal, uint64_t& rhVal) {
    fhVal = getFhval(kmerSeq, k);
    rhVal = getRhval(kmerSeq, k);
    return (rhVal<fhVal)? rhVal : fhVal;
}

// canonical ntHash for sliding k-mers, using rotate ops
inline uint64_t NTC64(const unsigned char charOut, const unsigned char charIn, const unsigned k, uint64_t& fhVal, uint64_t& rhVal) {
    fhVal = rol(fhVal, 1) ^ rol(seedTab[charOut], k) ^ seedTab[charIn];
    rhVal = ror(rhVal, 1) ^ ror(seedTab[charOut&cpOff], 1) ^ rol(seedTab[charIn&cpOff], k-1);
    return (rhVal<fhVal)? rhVal : fhVal;
}

// canonical ntBase with seeding option, using rotate ops
inline uint64_t NTC64(const char * kmerSeq, const unsigned k, const unsigned seed) {
    uint64_t hVal = NTC64(kmerSeq,k);
    hVal *= seed ^ k * multiSeed;
    hVal ^= hVal >> multiShift;
    return hVal;
}

/*
 * Using pre-computed seed matrix msTab instead of rotate opts
*/

// ntHash
inline uint64_t NTP64(const char * kmerSeq, const unsigned k) {
    uint64_t hVal=0;
    for(unsigned i=0; i<k; i++)
        hVal ^= msTab[(unsigned char)kmerSeq[i]][(k-1-i)%64];
    return hVal;
}

// ntHash for sliding k-mers
inline uint64_t NTP64(const uint64_t fhVal, const unsigned char charOut, const unsigned char charIn, const unsigned k) {
    return(rol(fhVal, 1) ^ msTab[charOut][k%64] ^ msTab[charIn][0]);
}

// ntBase with seeding option
inline uint64_t NTP64(const char * kmerSeq, const unsigned k, const unsigned seed) {
    uint64_t hVal=0;
    for(unsigned i=0; i<k; i++)
        hVal ^= msTab[(unsigned char)kmerSeq[i]][(k-1-i)%64];
    hVal *= seed ^ k * multiSeed;
    hVal ^= hVal >> multiShift;
    return hVal;
}

// canonical ntBase
inline uint64_t NTPC64(const char * kmerSeq, const unsigned k) {
    uint64_t fhVal=0, rhVal=0;
    for(unsigned i=0; i<k; i++) {
        fhVal ^= msTab[(unsigned char)kmerSeq[i]][(k-1-i)%64];
        rhVal ^= msTab[(unsigned char)kmerSeq[i]&cpOff][i%64];
    }
    return (rhVal<fhVal)? rhVal : fhVal;
}

// canonical ntHash
inline uint64_t NTPC64(const char * kmerSeq, const unsigned k, uint64_t& fhVal, uint64_t& rhVal) {
    fhVal=0, rhVal=0;
    for(unsigned i=0; i<k; i++) {
        fhVal ^= msTab[(unsigned char)kmerSeq[i]][(k-1-i)%64];
        rhVal ^= msTab[(unsigned char)kmerSeq[i]&cpOff][i%64];
    }
    return (rhVal<fhVal)? rhVal : fhVal;
}

// canonical ntHash for sliding k-mers
inline uint64_t NTPC64(const unsigned char charOut, const unsigned char charIn, const unsigned k, uint64_t& fhVal, uint64_t& rhVal) {
    fhVal = rol(fhVal, 1) ^ msTab[charOut][k%64] ^ msTab[charIn][0];
    rhVal = ror(rhVal, 1) ^ msTab[charOut&cpOff][63] ^ msTab[charIn&cpOff][(k-1)%64];
    return (rhVal<fhVal)? rhVal : fhVal;
}

// canonical ntBase with seeding option
inline uint64_t NTPC64(const char * kmerSeq, const unsigned k, const unsigned seed) {
    uint64_t hVal = NTPC64(kmerSeq,k);
    hVal *= seed ^ k * multiSeed;
    hVal ^= hVal >> multiShift;
    return hVal;
}

// multihash ntHash, ntBase
inline void NTM64(const char * kmerSeq, const unsigned k, const unsigned m, uint64_t *hVal) {
    uint64_t bVal=0, tVal=0;
    bVal = NTP64(kmerSeq, k);
    hVal[0] = bVal;
    for(unsigned i=1; i<m; i++) {
        tVal = bVal * (i ^ k * multiSeed);
        tVal ^= tVal >> multiShift;
        hVal[i] =  tVal;
    }
}

// one extra hash for given base hash
inline uint64_t NTE64(const uint64_t hVal, const unsigned k, const unsigned i) {
    uint64_t tVal = hVal;
    tVal *= (i ^ k * multiSeed);
    tVal ^= tVal >> multiShift;
    return tVal;
}

// multihash ntHash for sliding k-mers
inline void NTM64(const unsigned char charOut, const unsigned char charIn, const unsigned k, const unsigned m, uint64_t *hVal) {
    uint64_t bVal=0, tVal=0;
    bVal = rol(hVal[0], 1) ^ msTab[charOut][k%64] ^ msTab[charIn][0];
    hVal[0] = bVal;
    for(unsigned i=1; i<m; i++) {
        tVal = bVal * (i ^ k * multiSeed);
        tVal ^= tVal >> multiShift;
        hVal[i] =  tVal;
    }
}

// canonical multihash ntBase
inline void NTMC64(const char * kmerSeq, const unsigned k, const unsigned m, uint64_t *hVal) {
    uint64_t bVal=0, tVal=0;
    bVal = NTPC64(kmerSeq, k);
    hVal[0] = bVal;
    for(unsigned i=1; i<m; i++) {
        tVal = bVal * (i ^ k * multiSeed);
        tVal ^= tVal >> multiShift;
        hVal[i] =  tVal;
    }
}

// canonical multihash ntHash
inline void NTMC64(const char * kmerSeq, const unsigned k, const unsigned m, uint64_t& fhVal, uint64_t& rhVal, uint64_t *hVal) {
    uint64_t bVal=0, tVal=0;
    bVal = NTPC64(kmerSeq, k, fhVal, rhVal);
    hVal[0] = bVal;
    for(unsigned i=1; i<m; i++) {
        tVal = bVal * (i ^ k * multiSeed);
        tVal ^= tVal >> multiShift;
        hVal[i] =  tVal;
    }
}

// canonical multihash ntHash for sliding k-mers
inline void NTMC64(const unsigned char charOut, const unsigned char charIn, const unsigned k, const unsigned m, uint64_t& fhVal, uint64_t& rhVal, uint64_t *hVal) {
    uint64_t bVal=0, tVal=0;
    bVal = NTPC64(charOut, charIn, k, fhVal, rhVal);
    hVal[0] = bVal;
    for(unsigned i=1; i<m; i++) {
        tVal = bVal * (i ^ k * multiSeed);
        tVal ^= tVal >> multiShift;
        hVal[i] =  tVal;
    }
}

/*
 * ignoring k-mers containing nonACGT using ntHash function
*/

// canonical ntBase
inline bool NTPC64(const char *kmerSeq, const unsigned k, uint64_t& hVal, unsigned& locN) {
    hVal=0;
    locN=0;
    uint64_t fhVal=0,rhVal=0;
    for(int i=k-1; i>=0; i--) {
        if(msTab[(unsigned char)kmerSeq[i]][(k-1-i)%64]==seedN) {
            locN=i;
            return false;
        }
        fhVal ^= msTab[(unsigned char)kmerSeq[i]][(k-1-i)%64];
        rhVal ^= msTab[(unsigned char)kmerSeq[i]&cpOff][i%64];
    }
    hVal = (rhVal<fhVal)? rhVal : fhVal;
    return true;
}

// canonical multihash ntBase
inline bool NTMC64(const char *kmerSeq, const unsigned k, const unsigned m, unsigned& locN, uint64_t* hVal) {
    uint64_t bVal=0, tVal=0, fhVal=0, rhVal=0;
    locN=0;
    for(int i=k-1; i>=0; i--) {
        if(msTab[(unsigned char)kmerSeq[i]][(k-1-i)%64]==seedN) {
            locN=i;
            return false;
        }
        fhVal ^= msTab[(unsigned char)kmerSeq[i]][(k-1-i)%64];
        rhVal ^= msTab[(unsigned char)kmerSeq[i]&cpOff][i%64];
    }
    bVal = (rhVal<fhVal)? rhVal : fhVal;
    hVal[0] = bVal;
    for(unsigned i=1; i<m; i++) {
        tVal = bVal * (i ^ k * multiSeed);
        tVal ^= tVal >> multiShift;
        hVal[i] =  tVal;
    }
    return true;
}

// canonical ntHash
inline bool NTPC64(const char *kmerSeq, const unsigned k, uint64_t& fhVal, uint64_t& rhVal, uint64_t& hVal, unsigned& locN) {
    hVal=fhVal=rhVal=0;
    locN=0;
    for(int i=k-1; i>=0; i--) {
        if(msTab[(unsigned char)kmerSeq[i]][(k-1-i)%64]==seedN) {
            locN=i;
            return false;
        }
        fhVal ^= msTab[(unsigned char)kmerSeq[i]][(k-1-i)%64];
        rhVal ^= msTab[(unsigned char)kmerSeq[i]&cpOff][i%64];
    }
    hVal = (rhVal<fhVal)? rhVal : fhVal;
    return true;
}

// canonical multihash ntHash
inline bool NTMC64(const char *kmerSeq, const unsigned k, const unsigned m, uint64_t& fhVal, uint64_t& rhVal, unsigned& locN, uint64_t* hVal) {
    fhVal=rhVal=0;
    uint64_t bVal=0, tVal=0;
    locN=0;
    for(int i=k-1; i>=0; i--) {
        if(msTab[(unsigned char)kmerSeq[i]][(k-1-i)%64]==seedN) {
            locN=i;
            return false;
        }
        fhVal ^= msTab[(unsigned char)kmerSeq[i]][(k-1-i)%64];
        rhVal ^= msTab[(unsigned char)kmerSeq[i]&cpOff][i%64];
    }
    bVal = (rhVal<fhVal)? rhVal : fhVal;
    hVal[0] = bVal;
    for(unsigned i=1; i<m; i++) {
        tVal = bVal * (i ^ k * multiSeed);
        tVal ^= tVal >> multiShift;
        hVal[i] =  tVal;
    }
    return true;
}

#endif
