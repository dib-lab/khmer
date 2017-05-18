
/**
* High performance random generator.
* Mersenne Twister

@article{matsumoto1998mtd,
  title={{Mersenne Twister: A 623-Dimensionally Equidistributed Uniform Pseudo-Random Number Generator}},
  author={MATSUMOTO, M. and NISHIMURA, T.},
  journal={ACM Transactions on Modeling and Computer Simulation},
  volume={8},
  number={1},
  pages={3-30},
  year={1998}
}
*/
// MersenneTwister.h
// Mersenne Twister random number generator -- a C++ class MTRand
// Based on code by Makoto Matsumoto, Takuji Nishimura, and Shawn Cokus
// Richard J. Wagner  v1.0  15 May 2003  rjwagner@writeme.com

// The Mersenne Twister is an algorithm for generating random numbers.  It
// was designed with consideration of the flaws in various other generators.
// The period, 2^19937-1, and the order of equidistribution, 623 dimensions,
// are far greater.  The generator is also fast; it avoids multiplication and
// division, and it benefits from caches and pipelines.  For more information
// see the inventors' web page at http://www.math.keio.ac.jp/~matumoto/emt.html

// Reference
// M. Matsumoto and T. Nishimura, "Mersenne Twister: A 623-Dimensionally
// Equidistributed Uniform Pseudo-Random Number Generator", ACM Transactions on
// Modeling and Computer Simulation, Vol. 8, No. 1, January 1998, pp 3-30.

// Copyright (C) 1997 - 2002, Makoto Matsumoto and Takuji Nishimura,
// Copyright (C) 2000 - 2003, Richard J. Wagner
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions
// are met:
//
//   1. Redistributions of source code must retain the above copyright
//      notice, this list of conditions and the following disclaimer.
//
//   2. Redistributions in binary form must reproduce the above copyright
//      notice, this list of conditions and the following disclaimer in the
//      documentation and/or other materials provided with the distribution.
//
//   3. The names of its contributors may not be used to endorse or promote
//      products derived from this software without specific prior written
//      permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT OWNER OR
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

// The original code included the following notice:
//
//     When you use this, send an email to: matumoto@math.keio.ac.jp
//     with an appropriate reference to your work.
//
// It would be nice to CC: rjwagner@writeme.com and Cokus@math.washington.edu
// when you write.

#ifndef ROLLINGHASH_MERSENNETWISTER_H
#define ROLLINGHASH_MERSENNETWISTER_H

// Not thread safe (unless auto-initialization is avoided and each thread has
// its own MTRand object)

#include <iostream>
#include <limits.h>
#include <stdio.h>
#include <time.h>
#include <math.h>

class MTRand {
// Data
public:
    typedef unsigned long uint32;  // unsigned integer type, at least 32 bits

    enum { N = 624 };       // length of state vector
    enum { SAVE = N + 1 };  // length of array for save()

protected:
    enum { M = 397 };  // period parameter

    uint32 state[N];   // internal state
    uint32 *pNext;     // next value to get from state
    int left;          // number of values left before reload needed


//Methods
public:
    MTRand( const uint32& oneSeed );  // initialize with a simple uint32
    MTRand( uint32 *const bigSeed, uint32 const seedLength = N );  // or an array
    MTRand();  // auto-initialize with /dev/urandom or time() and clock()

    // Do NOT use for CRYPTOGRAPHY without securely hashing several returned
    // values together, otherwise the generator state can be learned after
    // reading 624 consecutive values.

    // Access to 32-bit random numbers
    double rand();                          // real number in [0,1]
    double rand( const double& n );         // real number in [0,n]
    double randExc();                       // real number in [0,1)
    double randExc( const double& n );      // real number in [0,n)
    double randDblExc();                    // real number in (0,1)
    double randDblExc( const double& n );   // real number in (0,n)
    uint32 randInt();                       // integer in [0,2^32-1]
    uint32 randInt( const uint32& n );      // integer in [0,n] for n < 2^32
    double operator()() {
        return rand();    // same as rand()
    }

    // Access to 53-bit random numbers (capacity of IEEE double precision)
    double rand53();  // real number in [0,1)

    // Access to nonuniform random number distributions
    double randNorm( const double& mean = 0.0, const double& variance = 0.0 );

    // Re-seeding functions with same behavior as initializers
    void seed( const uint32 oneSeed );
    void seed( uint32 *const bigSeed, const uint32 seedLength = N );
    void seed();

    // Saving and loading generator state
    void save( uint32* saveArray ) const;  // to array of size SAVE
    void load( uint32 *const loadArray );  // from such array
    friend std::ostream& operator<<( std::ostream& os, const MTRand& mtrand );
    friend std::istream& operator>>( std::istream& is, MTRand& mtrand );

protected:
    void initialize( const uint32 oneSeed );
    void reload();
    uint32 hiBit( const uint32& u ) const {
        return u & 0x80000000UL;
    }
    uint32 loBit( const uint32& u ) const {
        return u & 0x00000001UL;
    }
    uint32 loBits( const uint32& u ) const {
        return u & 0x7fffffffUL;
    }
    uint32 mixBits( const uint32& u, const uint32& v ) const
    {
        return hiBit(u) | loBits(v);
    }
    uint32 twist( const uint32& m, const uint32& s0, const uint32& s1 ) const
    {
        return m ^ (mixBits(s0,s1)>>1) ^ (-static_cast<long>(loBit(s1)) & 0x9908b0dfUL);
    }
    static uint32 hash( time_t t, clock_t c );
};


#endif  // MERSENNETWISTER_H

// Change log:
//
// v0.1 - First release on 15 May 2000
//      - Based on code by Makoto Matsumoto, Takuji Nishimura, and Shawn Cokus
//      - Translated from C to C++
//      - Made completely ANSI compliant
//      - Designed convenient interface for initialization, seeding, and
//        obtaining numbers in default or user-defined ranges
//      - Added automatic seeding from /dev/urandom or time() and clock()
//      - Provided functions for saving and loading generator state
//
// v0.2 - Fixed bug which reloaded generator one step too late
//
// v0.3 - Switched to clearer, faster reload() code from Matthew Bellew
//
// v0.4 - Removed trailing newline in saved generator format to be consistent
//        with output format of built-in types
//
// v0.5 - Improved portability by replacing static const int's with enum's and
//        clarifying return values in seed(); suggested by Eric Heimburg
//      - Removed MAXINT constant; use 0xffffffffUL instead
//
// v0.6 - Eliminated seed overflow when uint32 is larger than 32 bits
//      - Changed integer [0,n] generator to give better uniformity
//
// v0.7 - Fixed operator precedence ambiguity in reload()
//      - Added access for real numbers in (0,1) and (0,n)
//
// v0.8 - Included time.h header to properly support time_t and clock_t
//
// v1.0 - Revised seeding to match 26 Jan 2002 update of Nishimura and Matsumoto
//      - Allowed for seeding with arrays of any length
//      - Added access for real numbers in [0,1) with 53-bit resolution
//      - Added access for real numbers from normal (Gaussian) distributions
//      - Increased overall speed by optimizing twist()
//      - Doubled speed of integer [0,n] generation
//      - Fixed out-of-range number generation on 64-bit machines
//      - Improved portability by substituting literal constants for long enum's
//      - Changed license from GNU LGPL to BSD
