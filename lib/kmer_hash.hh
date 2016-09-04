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

#include "khmer.hh"

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
// two-way hash functions.
HashIntoType _hash(const char * kmer, const WordLength k);
HashIntoType _hash(const char * kmer, const WordLength k,
                   HashIntoType& h, HashIntoType& r);
HashIntoType _hash(const std::string kmer, const WordLength k);
HashIntoType _hash(const std::string kmer, const WordLength k,
                   HashIntoType& h, HashIntoType& r);
HashIntoType _hash_forward(const char * kmer, WordLength k);

std::string _revhash(HashIntoType hash, WordLength k);

// two-way hash functions, MurmurHash3.
HashIntoType _hash_murmur(const std::string& kmer);
HashIntoType _hash_murmur(const std::string& kmer,
                          HashIntoType& h, HashIntoType& r);
HashIntoType _hash_murmur_forward(const std::string& kmer);

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
class Kmer
{

public:

    std::string kmer_s;
    WordLength ksize;
    HashIntoType kmer_u, kmer_f, kmer_r;
    bool cached;

    /** @param[in]   s     DNA k-mer
        @param[in]   ksize k-mer size
     */
    Kmer(const std::string s, WordLength k)
    {
        kmer_s = s;
        ksize = k;

        kmer_u = _hash(s, k, kmer_f, kmer_r);
        cached = true;
    }

    Kmer(const std::string s, WordLength k,
         HashIntoType ku, HashIntoType kf, HashIntoType kr)
    {
        kmer_s = s;
        ksize = k;
        kmer_u = ku;
        kmer_f = kf;
        kmer_r = kr;
        cached = true;
    }

    Kmer()
    {
        kmer_s = "";
        ksize = 0;
        kmer_u = kmer_f = kmer_r = 0;
        cached = false;
    }

    /// Allows complete backwards compatibility
    operator HashIntoType() const
    {
        if (cached) {
            return kmer_u;
        }
        return _hash(kmer_s, ksize);
    }

    bool operator< (Kmer &other)
    {
        HashIntoType me = *this;
        HashIntoType other_val = other;

        return me < other_val;
    }

    std::string get_string_rep(WordLength K) const
    {
        return kmer_s;
    }

};

/**
 * \class KmerFactory
 *
 * \brief Build complete Kmer objects.
 *
 * The KmerFactory is a simple construct to emit complete
 * Kmer objects. The design decision leading to this class
 * stems from the issue of overloading the Kmer constructor
 * while also giving it a K size: you get ambiguous signatures
 * between the (kmer_u, K) and (kmer_f, kmer_r) cases. This
 * implementation also allows a logical architecture wherein
 * KmerIterator and Hashtable inherit directly from KmerFactory,
 * extending the latter's purpose of "create k-mers" to
 * "emit k-mers from a sequence" and "create and store k-mers".
 *
 * \author Camille Scott
 *
 * Contact: camille.scott.w@gmail.com
 *
 */
class KmerFactory
{
protected:
    WordLength _ksize;

public:

    explicit KmerFactory(WordLength K): _ksize(K) {}

    /** @param[in]  kmer_u Uniqified hash value.
     *  @return A complete Kmer object.
     */
    Kmer build_kmer(const char * kmer_s)
    {
        return Kmer(kmer_s, _ksize);
    }

    /** Hash the given sequence and call the uniqify function
     *  on its results to build a complete Kmer.
     *
     *  @param[in]  kmer_s String representation of a k-mer.
     *  @return A complete Kmer object hashed from the given string.
     */
    Kmer build_kmer(std::string kmer_s)
    {
        return Kmer(kmer_s.c_str(), _ksize);
    }
    
    Kmer build_kmer(std::string kmer_s,
                    HashIntoType kmer_f, HashIntoType kmer_r)
    {
        return Kmer(kmer_s.c_str(), _ksize, uniqify_rc(kmer_f, kmer_r),
                    kmer_f, kmer_r);
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
class KmerIterator: public KmerFactory
{
protected:
    const char * _seq;

    HashIntoType _kmer_f, _kmer_r;
    unsigned int index;
    size_t length;
    bool initialized;
    HashIntoType bitmask;
    unsigned int _nbits_sub_1;
public:
    KmerIterator(const char * seq, unsigned char k) :
        KmerFactory(k), _seq(seq)
    {
        index = _ksize;
        length = strlen(_seq);
        _kmer_f = 0;
        _kmer_r = 0;
        initialized = false;

        bitmask = 0;
        for (unsigned char i = 0; i < _ksize; i++) {
            bitmask = (bitmask << 2) | 3;
        }
        _nbits_sub_1 = (_ksize*2 - 2);
    }

    Kmer first()
    {
        initialized = true;

        _hash(_seq + index - _ksize, _ksize, _kmer_f, _kmer_r);
        Kmer k = build_kmer(_seq + index - _ksize, _kmer_f, _kmer_r);
        index++;
        return k;
    }

    Kmer next()
    {
        if (!initialized) {
            return first();
        }
        if (done()) {
            throw khmer_exception();
        }

        char ch = _seq[index];

        // left-shift the previous hash over
        _kmer_f = _kmer_f << 2;

        // 'or' in the current nt
        _kmer_f |= twobit_repr(ch);

        // mask off the 2 bits we shifted over.
        _kmer_f &= bitmask;

        // now handle reverse complement
        _kmer_r = _kmer_r >> 2;
        _kmer_r |= (twobit_comp(ch) << _nbits_sub_1);

        Kmer k = build_kmer(_seq + index - _ksize, _kmer_f, _kmer_r);
        index++;
        return k;
    }

    /// @return Whether or not the iterator has completed.
    bool done()
    {
        return index > length;
    }

    unsigned int get_start_pos() const
    {
        if (!initialized) {
            throw khmer_exception();
        }
        return index - _ksize - 1;
    }

    unsigned int get_end_pos() const
    {
        return index - 1;
    }
}; // class KmerIterator

}

#endif // KMER_HASH_HH
