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
#ifndef HASHING_HH 
#define HASHING_HH



#include "oxli/oxli.hh"
#include "oxli/alphabets.hh"
#include "oxli/kmer_hash.hh"
#include "oxli/oxli_exception.hh"

#include <string>
#include <utility>
#include <iostream>
#include <functional>


namespace oxli {


class HashedKmer
{

public:
    
    HashIntoType kmer_fw;
    HashIntoType kmer_rc;
    
    HashedKmer(HashIntoType kmer_fw, HashIntoType kmer_rc) : 
        kmer_fw(kmer_fw),
        kmer_rc(kmer_rc)
    {
    }

    HashedKmer() : kmer_fw(0), kmer_rc(0) {}

    HashedKmer(HashedKmer&& other)                 = default;
    HashedKmer(const HashedKmer& other)            = default;
    HashedKmer& operator=(const HashedKmer& other) = default;
    HashedKmer& operator=(HashedKmer&& other)      = default;

    operator HashIntoType() const
    {
       return uniqify_rc(kmer_fw, kmer_rc);
    }

    bool forward()  const
    {
        return kmer_fw == uniqify_rc(kmer_fw, kmer_rc);
    }

    bool operator< (const HashedKmer &other) const
    {
        return key() < other.key();
    }

    inline HashIntoType key() const {
        return uniqify_rc(kmer_fw, kmer_rc);
    }

};


class CompleteKmer: public HashedKmer
{

public:
    const std::string kmer_s;
    
    CompleteKmer(HashIntoType kmer_fw, 
                 HashIntoType kmer_rc,
                 std::string kmer_s) : 
        HashedKmer(kmer_fw, kmer_rc),
        kmer_s(kmer_s)
    {
    }

    CompleteKmer(const HashedKmer& kmer, std::string kmer_s) :
        HashedKmer(kmer),
        kmer_s(kmer_s)
    {
    }

    CompleteKmer(HashedKmer&& kmer, std::string kmer_s) :
        HashedKmer(std::move(kmer)),
        kmer_s(kmer_s)
    {
    }

    CompleteKmer(CompleteKmer&& other)                 = default;
    CompleteKmer(const CompleteKmer& other)            = default;
    CompleteKmer& operator=(const CompleteKmer& other) = default;
    CompleteKmer& operator=(CompleteKmer&& other)      = default;

    CompleteKmer revcomp() const
    {
        return CompleteKmer(kmer_fw, kmer_rc, _revcomp(kmer_s));
    }

};


class IHashedKmerIterator {
public:
    virtual WordLength ksize() const = 0;
    virtual HashedKmer next() = 0;
    virtual bool done() const = 0;
    virtual ~IHashedKmerIterator() { };
};


class ICompleteKmerIterator {
public:
    virtual WordLength ksize() const = 0;
    virtual CompleteKmer next() = 0;
    virtual bool done() const = 0;
    virtual ~ICompleteKmerIterator() { };
};


class KmerHasher;
class ShiftingKmerHasher;


class SequenceHasher: public IHashedKmerIterator
{
private:
    uint64_t index;
    const char * _sequence;
    const uint64_t length;
    HashedKmer previous;
    const KmerHasher * hasher;
    const WordLength _ksize;

    virtual HashedKmer _first(const KmerHasher * hasher);

    virtual HashedKmer _next(const KmerHasher * hasher);

    virtual HashedKmer _next(const ShiftingKmerHasher * hasher);

public:
    SequenceHasher(const std::string& sequence, const KmerHasher * hasher);

    WordLength ksize() const {
        return _ksize;
    }

    HashedKmer next()
    {
        if(done()) {
            throw oxli_exception("Iterator is done.");
        }

        if (index == 0) {
            index++;
            return _first(hasher);
        }

        return _next(hasher);
    }

        
    bool done() const
    {
        return index > (length - _ksize);
    }

    uint64_t get_start_pos() const
    {
        return index;
    }

    uint64_t get_end_pos() const
    {
        return index + _ksize - 1;
    }
};


class KmerHasher
{
public:

    virtual WordLength ksize() const = 0;
    virtual HashedKmer hash_dna(const char * kmer) const = 0;
    virtual HashedKmer hash_dna(const std::string& kmer) const = 0;

    virtual HashIntoType hash_dna_fw(const std::string& kmer) const = 0;
    virtual HashIntoType hash_dna_fw(const char * kmer) const = 0;

    /*
    virtual CompleteKmer hash_dna_complete(const std::string& kmer) const
    {
        return CompleteKmer(std::move(hash_dna(kmer)), kmer);
    }
    virtual CompleteKmer hash_dna_complete(const char * kmer) const
    {
        return CompleteKmer(std::move(hash_dna(kmer)), std::string(kmer));
    }
    */

    virtual SequenceHasher hash_sequence(const std::string& sequence) const = 0;
    //virtual CompleteKmerIterator hash_sequence_complete(const std::string& sequence) const = 0;
};


class ReversibleKmerHasher: public virtual KmerHasher
{
public:

    virtual std::string unhash_dna(HashedKmer& kmer) const = 0;

};


class ShiftingKmerHasher: public virtual KmerHasher
{
public:

    virtual HashedKmer _get_left(const HashedKmer& right,
                                 const char base) const = 0;

    //virtual CompleteKmer _get_left_complete(const CompleteKmer& right,
    //                                        const char base) const = 0;

    virtual HashedKmer _get_right(const HashedKmer& left,
                                  const char base) const = 0;

    //virtual CompleteKmer _get_right_complete(const CompleteKmer& left,
    //                                         const char base) const = 0;

};

class TwoBitKmerHasher: public ReversibleKmerHasher,
                        public ShiftingKmerHasher
{
protected:
    
    const WordLength _ksize;
    HashIntoType bitmask;
    unsigned int rc_left_shift;

public:

    TwoBitKmerHasher(WordLength K);

    ~TwoBitKmerHasher() {}

    virtual WordLength ksize() const {
        return _ksize;
    }

    virtual HashedKmer hash_dna(const std::string& kmer) const
    {
        return hash_dna(kmer.data());
    }

    virtual HashedKmer hash_dna(const char * kmer) const
    {
        HashIntoType kmer_fw, kmer_rc;
        _hash(kmer, _ksize, kmer_fw, kmer_rc);
        return HashedKmer(kmer_fw, kmer_rc);
    }

    virtual HashIntoType hash_dna_fw(const std::string& kmer) const
    {
        return _hash_forward(kmer.data(), _ksize);
    }

    virtual HashIntoType hash_dna_fw(const char * kmer) const
    {
        return _hash_forward(kmer, _ksize);
    }

    virtual std::string unhash_dna(HashedKmer& kmer) const
    {
        return unhash_dna(kmer.key());
    }

    virtual std::string unhash_dna(HashIntoType kmer) const
    {
        return _revhash(kmer, _ksize);
    }

    virtual HashedKmer _get_left(const HashedKmer& right,
                                 const char base) const;

    virtual HashedKmer _get_right(const HashedKmer& left,
                                  const char base) const;

    virtual CompleteKmer _get_left_complete(const CompleteKmer& right,
                                            const char base) const
    {
        std::string left = base + right.kmer_s.substr(1, _ksize - 1);
        return CompleteKmer(_get_left(right, base), left);
    }

    virtual CompleteKmer _get_right_complete(const CompleteKmer& left,
                                                 const char base) const
    {
        std::string right = left.kmer_s.substr(0, _ksize - 1) + base;
        return CompleteKmer(_get_right(left, base), right);
    }


    virtual SequenceHasher hash_sequence(const std::string& sequence) const
    {
        return SequenceHasher(sequence, this);
    }
};


class Murmur3KmerHasher: public KmerHasher
{
protected:
    const WordLength _ksize;

public:
    Murmur3KmerHasher(WordLength K): KmerHasher(), _ksize(K) {}

    virtual WordLength ksize() const {
        return _ksize;
    }

    virtual HashedKmer hash_dna(const std::string& kmer) const 
    {
        HashIntoType kmer_fw, kmer_rc;
        _hash_murmur(kmer, _ksize, kmer_fw, kmer_rc);
        return HashedKmer(kmer_fw, kmer_rc);
    }

    virtual HashedKmer hash_dna(const char * kmer) const
    {
        HashIntoType kmer_fw, kmer_rc;
        std::string _kmer(kmer);
        _hash_murmur(_kmer, _ksize, kmer_fw, kmer_rc);
        return HashedKmer(kmer_fw, kmer_rc);
    }

    virtual HashedKmer hash_dna(const char * kmer,
                                const char * revkmer) const
    {
        HashIntoType kmer_fw, kmer_rc;
        kmer_fw = hash_dna_fw(kmer);
        kmer_rc = hash_dna_fw(revkmer);
        return HashedKmer(kmer_fw, kmer_rc);
    }

    virtual HashIntoType hash_dna_fw(const std::string& kmer) const
    {
        return _hash_murmur_forward(kmer, _ksize);
    }

    virtual HashIntoType hash_dna_fw(const char * kmer) const;

};

}

#endif
