#include "MurmurHash3.h"
#include "oxli/hashing.hh"
#include <iostream>

namespace oxli
{



TwoBitKmerHasher::TwoBitKmerHasher(WordLength K) : 
    KmerHasher(), _ksize(K)
{
    std::cout << "TwoBitKmerHasher(" << static_cast<unsigned>(_ksize) << ")" << std::endl;
    bitmask = 0;
    for (unsigned int i = 0; i < _ksize; i++) {
        bitmask = (bitmask << 2) | 3;
    }
    rc_left_shift = _ksize * 2 - 2;
}


HashedKmer TwoBitKmerHasher::_get_left(const HashedKmer& right,
                                               const char base) const
{
    HashIntoType kmer_fw, kmer_rc;
    kmer_fw = ((right.kmer_fw) >> 2 | twobit_repr(base) << rc_left_shift);
    kmer_rc = (((right.kmer_rc) << 2) & bitmask) | (twobit_comp(base));
    return HashedKmer(kmer_fw, kmer_rc);
}


HashedKmer TwoBitKmerHasher::_get_right(const HashedKmer& left,
                                                const char base) const
{
    HashIntoType kmer_fw, kmer_rc;
    kmer_fw = (((left.kmer_fw) << 2) & bitmask) | (twobit_repr(base));
    kmer_rc = ((left.kmer_rc) >> 2) | (twobit_comp(base) << rc_left_shift);
    return HashedKmer(kmer_fw, kmer_rc);
}


SequenceHasher::SequenceHasher(const std::string& sequence, const KmerHasher * hasher) :
    hasher(hasher),
    index(0),
    _sequence(sequence.data()),
    _ksize(hasher->ksize()),
    length(sequence.length())
{
    std::cout << "SequenceHasher(" << static_cast<unsigned>(_ksize) << ")" << std::endl;
    std::cout << "_sequence = " << _sequence << std::endl;
}


HashedKmer SequenceHasher::_first(const KmerHasher * hasher)
{
    if (_sequence == NULL) {
        throw oxli_exception("Invalid sequence pointer.");
    }
    previous = hasher->hash_dna(_sequence);
    return previous;
}

HashedKmer SequenceHasher::_next(const KmerHasher * hasher)
{
    //if (!(index <= (length)) {
    //    throw oxli_exception("KmerIterator index <= length; should have finished.");
    //}

    const char * kmer = _sequence + index;
    HashedKmer next_kmer = hasher->hash_dna(kmer);
    previous = next_kmer;
    index++;
    
    return next_kmer;
}


HashedKmer SequenceHasher::_next(const ShiftingKmerHasher * hasher)
{
    //if (!(index <= length)) {
    //    throw oxli_exception("KmerIterator index <= length; should have finished.");
    //}

    unsigned char next_base = _sequence[index+_ksize];
    HashedKmer next_kmer = hasher->_get_right(previous, next_base);
    previous = next_kmer;
    index++;

    return next_kmer;
}


HashIntoType Murmur3KmerHasher::hash_dna_fw(const char * kmer) const
{
    uint64_t out[2];
    uint32_t seed = 0;
    MurmurHash3_x64_128((void *)kmer, _ksize, seed, &out);
    return out[0];
}

}
