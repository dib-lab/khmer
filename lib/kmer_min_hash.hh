#ifndef KMER_MIN_HASH_HH
#define KMER_MIN_HASH_HH

#include <set>
#include <map>
#include <exception>
#include <string>

#include "MurmurHash3.h"
#include "khmer.hh"
#include "kmer_hash.hh"

uint64_t _hash_murmur64(const std::string& kmer);

namespace khmer {

////

typedef std::set<HashIntoType> CMinHashType;

class minhash_exception : public std::exception
{
public:
    explicit minhash_exception(const std::string& msg = "Generic minhash exception")
        : _msg(msg) { }

    virtual ~minhash_exception() throw() { }
    virtual const char* what() const throw ()
    {
        return _msg.c_str();
    }

protected:
    const std::string _msg;
};


class KmerMinHash
{
public:
    const unsigned int num;
    const unsigned int ksize;
    const bool is_protein;
    CMinHashType mins;

    KmerMinHash(unsigned int n, unsigned int k, bool prot) :
        num(n), ksize(k), is_protein(prot) { };

    void _shrink() {
        while (mins.size() > num) {
            CMinHashType::iterator mi = mins.end();
            mi--;
            mins.erase(mi);
        }
    }
    void add_hash(HashIntoType h) {
        mins.insert(h);
        _shrink();
    }
    void add_word(std::string word) {
        HashIntoType hash = _hash_murmur64(word);
        add_hash(hash);
    }
    void add_sequence(const char * sequence, bool force=false) {
        if (strlen(sequence) < ksize) {
            return;
        }
        const std::string seq = _checkdna(sequence, force);
        if (!is_protein) {
            for (unsigned int i = 0; i < seq.length() - ksize + 1; i++) {
                std::string kmer = seq.substr(i, ksize);
                std::string rc = _revcomp(kmer);
                if (kmer < rc) {
                    add_word(kmer);
                } else {
                    add_word(rc);
                }
            }
        } else {                      // protein
            for (unsigned int i = 0; i < seq.length() - ksize + 1; i ++) {
                std::string kmer = seq.substr(i, ksize);
                std::string aa = _dna_to_aa(kmer);

                add_word(aa);

                std::string rc = _revcomp(kmer);
                aa = _dna_to_aa(rc);

                add_word(aa);
            }
        }
    }

    std::string _dna_to_aa(const std::string& dna) {
        std::string aa;
        unsigned int dna_size = (dna.size() / 3) * 3; // floor it
        for (unsigned int j = 0; j < dna_size; j += 3) {
            std::string codon = dna.substr(j, 3);
            aa += (_codon_table)[codon];
        }
        return aa;
    }

    std::string _checkdna(const char * s, bool force=false) const {
        size_t seqsize = strlen(s);

        std::string seq = s;

        for (size_t i=0; i < seqsize; ++i) {
            switch(seq[i]) {
            case 'A':
            case 'C':
            case 'G':
            case 'T':
                break;
            default:
                if (force) {
                    seq[i] = 'N';
                } else {
                    std::string msg = "invalid DNA character in sequence: ";
                    msg += seq[i];
                    throw minhash_exception(msg);
                }
                break;
            }
        }
        return seq;
    }

    std::string _revcomp(const std::string& kmer) const {
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
            case 'N':
                complement = 'N';
                break;
            default:
                std::string msg = "invalid DNA character in sequence: ";
                msg += kmer[i];

                throw minhash_exception(msg);
            }
            out[ksize - i - 1] = complement;
        }
        return out;
    }

    void merge(const KmerMinHash& other) {
        CMinHashType::iterator mi;
        if (ksize != other.ksize) {
            throw minhash_exception("different ksizes cannot be merged");
        }
        if (is_protein != other.is_protein) {
            throw minhash_exception("DNA/prot minhashes cannot be merged");
        }
        for (mi = other.mins.begin(); mi != other.mins.end(); ++mi) {
            mins.insert(*mi);
        }
        _shrink();
    }
    unsigned int count_common(const KmerMinHash& other) {
        CMinHashType combined;

        if (ksize != other.ksize) {
            throw minhash_exception("different ksizes cannot be compared");
        }
        if (is_protein != other.is_protein) {
            throw minhash_exception("DNA/prot minhashes cannot be compared");
        }

        CMinHashType::iterator mi;
        for (mi = mins.begin(); mi != mins.end(); ++mi) {
            combined.insert(*mi);
        }
        for (mi = other.mins.begin(); mi != other.mins.end(); ++mi) {
            combined.insert(*mi);
        }
        return mins.size() + other.mins.size() - combined.size();
    }

private:
    std::map<std::string, std::string> _codon_table = {
        {"TTT", "F"}, {"TTC", "F"},
        {"TTA", "L"}, {"TTG", "L"},

        {"TCT", "S"}, {"TCC", "S"}, {"TCA", "S"}, {"TCG", "S"},

        {"TAT", "Y"}, {"TAC", "Y"},
        {"TAA", "*"}, {"TAG", "*"},

        {"TGT", "C"}, {"TGC", "C"},
        {"TGA", "*"},
        {"TGG", "W"},

        {"CTT", "L"}, {"CTC", "L"}, {"CTA", "L"}, {"CTG", "L"},

        {"CCT", "P"}, {"CCC", "P"}, {"CCA", "P"}, {"CCG", "P"},

        {"CAT", "H"}, {"CAC", "H"},
        {"CAA", "Q"}, {"CAG", "Q"},

        {"CGT", "R"}, {"CGC", "R"}, {"CGA", "R"}, {"CGG", "R"},

        {"ATT", "I"}, {"ATC", "I"}, {"ATA", "I"},
        {"ATG", "M"},

        {"ACT", "T"}, {"ACC", "T"}, {"ACA", "T"}, {"ACG", "T"},

        {"AAT", "N"}, {"AAC", "N"},
        {"AAA", "K"}, {"AAG", "K"},

        {"AGT", "S"}, {"AGC", "S"},
        {"AGA", "R"}, {"AGG", "R"},

        {"GTT", "V"}, {"GTC", "V"}, {"GTA", "V"}, {"GTG", "V"},

        {"GCT", "A"}, {"GCC", "A"}, {"GCA", "A"}, {"GCG", "A"},

        {"GAT", "D"}, {"GAC", "D"},
        {"GAA", "E"}, {"GAG", "E"},

        {"GGT", "G"}, {"GGC", "G"}, {"GGA", "G"}, {"GGG", "G"}
    };
};

typedef std::map<khmer::HashIntoType, khmer::TagSet> TagToTagSet;
typedef std::map<khmer::HashIntoType, khmer::HashIntoType> TagToHash;

class NeighborhoodMinHash
{
public:
    unsigned int ksize;
    bool is_protein;

    TagToTagSet tag_connections;
    TagToHash tag_to_hash;

    NeighborhoodMinHash(unsigned int _k,
                        bool _is_p=false) : ksize(_k),
                                            is_protein(_is_p)
    { ; }
};

class CombinedMinHash
{
public:
    TagSet tags;
    KmerMinHash * mh;

    CombinedMinHash() : mh(NULL) { }

    ~CombinedMinHash() {
        cleanup();
    }

    void cleanup() {
        delete mh;
        mh = NULL;
    }
};

}
#endif // KMER_MIN_HASH_HH
