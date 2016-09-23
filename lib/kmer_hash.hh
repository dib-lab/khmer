/*
This file is part of khmer, https://github.com/dib-lab/khmer/, and is
Copyright (C) 2010-2015, Michigan State University.
Copyright (C) 2015-2016, The Regents of the University of California.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are
met:

    * Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.

    * Redistributions in binary form must reproduce the above
      copyright notice, this list of conditions and the following
      disclaimer in the documentation and/or other materials provided
      with the distribution.

    * Neither the name of the Michigan State University nor the names
      of its contributors may be used to endorse or promote products
      derived from this software without specific prior written
      permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
"AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
LICENSE (END)

Contact: khmer-project@idyll.org
*/
#ifndef KMER_HASH_HH
#define KMER_HASH_HH

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <string>


#include "MurmurHash3.h"
#include "khmer.hh"
#include "khmer_exception.hh"

// test validity
#ifdef KHMER_EXTRA_SANITY_CHECKS
#   define is_valid_dna(ch) ((toupper(ch)) == 'A' || (toupper(ch)) == 'C' || \
			    (toupper(ch)) == 'G' || (toupper(ch)) == 'T')
#else
// NOTE: Assumes data is already sanitized as it should be by parsers.
//	     This assumption eliminates 4 function calls.
#   define is_valid_dna(ch) ((ch) == 'A' || (ch) == 'C' || \
			     (ch) == 'G' || (ch) == 'T')
#endif

// bit representation of A/T/C/G.
#ifdef KHMER_EXTRA_SANITY_CHECKS
#   define twobit_repr(ch) ((toupper(ch)) == 'A' ? 0LL : \
			    (toupper(ch)) == 'T' ? 1LL : \
			    (toupper(ch)) == 'C' ? 2LL : 3LL)
#else
// NOTE: Assumes data is already sanitized as it should be by parsers.
//	     This assumption eliminates 4 function calls.
#   define twobit_repr(ch) ((ch) == 'A' ? 0LL : \
			    (ch) == 'T' ? 1LL : \
			    (ch) == 'C' ? 2LL : 3LL)
#endif

#define revtwobit_repr(n) ((n) == 0 ? 'A' : \
                           (n) == 1 ? 'T' : \
                           (n) == 2 ? 'C' : 'G')

#ifdef KHMER_EXTRA_SANITY_CHECKS
#   define twobit_comp(ch) ((toupper(ch)) == 'A' ? 1LL : \
			    (toupper(ch)) == 'T' ? 0LL : \
			    (toupper(ch)) == 'C' ? 3LL : 2LL)
#else
// NOTE: Assumes data is already sanitized as it should be by parsers.
//	     This assumption eliminates 4 function calls.
#   define twobit_comp(ch) ((ch) == 'A' ? 1LL : \
			    (ch) == 'T' ? 0LL : \
			    (ch) == 'C' ? 3LL : 2LL)
#endif

// choose wisely between forward and rev comp.
#ifndef NO_UNIQUE_RC
#define uniqify_rc(f, r) ((f) < (r) ? (f) : (r))
#else
#define uniqify_rc(f,r)(f)
#endif


namespace khmer
{


std::string _revcomp(const std::string& kmer)
{
    std::string out = kmer;
    size_t ksize = out.size();

    for (size_t i=0; i < ksize; ++i) {
        char complement;

        switch(kmer[i]) {
        case 'A':
            complement = 'T';
            break;
        case 'C':
            complement = 'G';
            break;
        case 'G':
            complement = 'C';
            break;
        case 'T':
            complement = 'A';
            break;
        default:
            throw khmer::khmer_exception("Invalid base in read");
            break;
        }
        out[ksize - i - 1] = complement;
    }
    return out;
}


/**
 * \class Kmer
 *
 * \brief Hold the hash values corresponding to a single k-mer.
 *
 * This class stores the forward, reverse complement, and
 * uniqified hash values for a given k-mer. It also defines
 * some basic operators and a utility function for getting
 * the string representation of the sequence. This is meant
 * to replace the original inelegant macros used for hashing.
 *
 * \author Camille Scott
 *
 * Contact: camille.scott.w@gmail.com
 *
 */
template <typename HashType>
class Kmer
{

public:

    /// The forward hash
    HashType kmer_f;
    /// The reverse (complement) hash
    HashType kmer_r;
    /// The uniqified hash
    HashType kmer_u;

    /** @param[in]   f forward hash.
     *  @param[in]   r reverse (complement) hash.
     *  @param[in]   u uniqified hash.
     */
    Kmer(HashType f, HashType r, HashType u)
    {
        kmer_f = f;
        kmer_r = r;
        kmer_u = u;
    }

    Kmer(HashType f, HashType r)
    {
        kmer_u = uniqify_rc(f, r);
        kmer_f = f;
        kmer_r = r;
    }

    /// @warning The default constructor builds an invalid k-mer.
    Kmer()
    {
        kmer_f = kmer_r = kmer_u = 0;
    }

    /// Allows complete backwards compatibility
    operator HashType() const
    {
        return kmer_u;
    }

    bool operator< (const Kmer &other) const
    {
        return kmer_u < other.kmer_u;
    }

};


class BitRepFunctor
{

public:

    WordLength K;
    HashIntoType bitmask;
    unsigned int rc_left_shift;
    unsigned int _nbits_sub_1;

    explicit BitRepFunctor(WordLength K): 
        K(K), bitmask(0) {

        for (unsigned char i = 0; i < K; i++) {
            bitmask = (bitmask << 2) | 3;
        }
        rc_left_shift = K * 2 - 2;
    }

    Kmer<HashIntoType> operator()(const char * sequence) const {
        // sizeof(HashIntoType) * 8 bits / 2 bits/base
        if (!(K <= sizeof(HashIntoType)*4) || !(strlen(sequence) >= K)) {
            throw khmer_exception("Supplied kmer string doesn't match the underlying k-size.");
        }

        HashIntoType h = 0, r = 0;

        h |= twobit_repr(sequence[0]);
        r |= twobit_comp(sequence[K-1]);

        for (WordLength i = 1, j = K - 2; i < K; i++, j--) {
            h = h << 2;
            r = r << 2;

            h |= twobit_repr(sequence[i]);
            r |= twobit_comp(sequence[j]);
        }

        return Kmer<HashIntoType>(h, r, uniqify_rc(h, r));
    }

    std::string inverse(Kmer<HashIntoType>& kmer) const {
        std::string s = "";
        HashIntoType hash = (HashIntoType) kmer;

        unsigned int val = hash & 3;
        s += revtwobit_repr(val);

        for (WordLength i = 1; i < K; i++) {
            hash = hash >> 2;
            val = hash & 3;
            s += revtwobit_repr(val);
        }

        reverse(s.begin(), s.end());

        return s;
    }

    Kmer<HashIntoType> operator()(Kmer<HashIntoType>& kmer, const char * seq) {
        HashIntoType kmer_f, kmer_r;
        kmer_f = (((kmer.kmer_f) << 2) & bitmask) | (twobit_repr(seq[K-1]));
        kmer_r = ((kmer.kmer_r) >> 2) | (twobit_comp(seq[K-1]) << rc_left_shift);

        return Kmer<HashIntoType>(kmer_f, kmer_r);
    }

    Kmer<HashIntoType> operator()(const char * seq, Kmer<HashIntoType>& kmer) {
        HashIntoType kmer_f, kmer_r;
        kmer_f = ((kmer.kmer_f) >> 2 | twobit_repr(seq[0]) << rc_left_shift);
        kmer_r = (((kmer.kmer_r) << 2) & bitmask) | (twobit_comp(seq[0]));
        
        return Kmer<HashIntoType>(kmer_f, kmer_r);
    }

};


class MurmurFunctor
{
public:
    Kmer operator()(const std::string& sequence) {
        HashIntoType out[2];
        uint32_t seed = 0;
        MurmurHash3_x64_128((void *)sequence.c_str(), sequence.size(), seed, &out);
        HashIntoType h = out[0];

        std::string rev = khmer::_revcomp(sequence);
        MurmurHash3_x64_128((void *)rev.c_str(), rev.size(), seed, &out);
        HashIntoType r = out[0];

        return Kmer(h, r, h ^ r);
    }
};


/**
 * \class KmerIterator
 *
 * \brief Emit Kmer objects generated from the given sequence.
 *
 * Given a string \f$S\f$ and a length \f$K > 0\f$, we define
 * the k-mers of \f$S\f$ as the set \f$S_{i..i+K} \forall i \in \{0..|S|-K+1\}\f$,
 * where \f$|S|\f$ is the length and \f$S_{j..k}\f$ is the half-open
 * substring starting at \f$j\f$ and terminating at \f$k\f$.
 *
 * KmerIterator mimics a python-style generator function which
 * emits the k-mers of the given sequence, in order, as Kmer objects.
 *
 * @warning This is not actually a valid C++ iterator, though it is close.
 *
 * \author Camille Scott
 *
 * Contact: camille.scott.w@gmail.com
 *
 */
template <class HashFunctorType>
class KmerIterator
{

protected:
    const char * _seq;
    WordLength _ksize;

    Kmer cur_kmer;
    unsigned int index;
    size_t length;
    bool initialized;
    HashFunctorType  _hash;
    
public:
    KmerIterator(const char * seq,
                 unsigned char ksize,
                 HashFunctorType hash_function)
                _seq(seq), _ksize(ksize), bitmask(0),
                _hash(hash_function)

    {
        index = _ksize - 1;
        length = strlen(_seq);
        initialized = false;
    }

    Kmer first()
    {
        index = 1;
        cur_kmer = _hash(_seq);

        return cur_kmer;
    }

    Kmer next()
    {
        if (done()) {
            throw khmer_exception();
        }

        if (!initialized) {
            initialized = true;
            return first();
        }

        const unsigned char * ch = _seq + index;
        index++;
        if (!(index <= length)) {
            throw khmer_exception();
        }

        cur_kmer = _hash(cur_kmer, ch);
        return cur_kmer;   
    }

    /// @return Whether or not the iterator has completed.
    bool done()
    {
        return index >= length;
    }

    unsigned int get_start_pos() const
    {
        return index - _ksize;
    }

    unsigned int get_end_pos() const
    {
        return index;
    }
}; // class KmerIterator

}

#endif // KMER_HASH_HH
