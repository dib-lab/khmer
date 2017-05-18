#ifndef NTHASH__ITERATOR_H
#define NTHASH__ITERATOR_H 1

#include <string>
#include <limits>
#include "nthash.hpp"


/**
 * Iterate over hash values for k-mers in a
 * given DNA sequence.
 *
 * This implementation uses ntHash
 * function to efficiently calculate
 * hash values for successive k-mers.
 */

class ntHashIterator
{

public:

    /**
     * Default constructor. Creates an iterator pointing to
     * the end of the iterator range.
    */
    ntHashIterator():
        m_hVec(NULL),
        m_pos(std::numeric_limits<std::size_t>::max())
    {}

    /**
     * Constructor.
     * @param seq address of DNA sequence to be hashed
     * @param k k-mer size
     * @param h number of hashes
    */
    ntHashIterator(const std::string& seq, unsigned h, unsigned k):
        m_seq(seq), m_h(h), m_k(k), m_hVec(new uint64_t[h]), m_pos(0)
    {
        init();
    }

    /** Initialize internal state of iterator */
    void init()
    {
        if (m_k > m_seq.length()) {
            m_pos = std::numeric_limits<std::size_t>::max();
            return;
        }
        unsigned locN=0;
        while (m_pos<m_seq.length()-m_k+1 && !NTMC64(m_seq.data()+m_pos, m_k, m_h, m_fhVal, m_rhVal, locN, m_hVec))
            m_pos+=locN+1;
        if (m_pos >= m_seq.length()-m_k+1)
            m_pos = std::numeric_limits<std::size_t>::max();
    }

    /** Advance iterator right to the next valid k-mer */
    void next()
    {
        if (m_pos >= m_seq.length()-m_k+1) {
            m_pos = std::numeric_limits<std::size_t>::max();
            return;
        }
        if(seedTab[(unsigned char)(m_seq.at(m_pos+m_k-1))]==seedN) {
            m_pos+=m_k;
            init();
        }
        else
            NTMC64(m_seq.at(m_pos-1), m_seq.at(m_pos-1+m_k), m_k, m_h, m_fhVal, m_rhVal, m_hVec);
    }

    /** get pointer to hash values for current k-mer */
    const uint64_t* operator*() const
    {
        return m_hVec;
    }

    /** test equality with another iterator */
    bool operator==(const ntHashIterator& it) const
    {
        return m_pos == it.m_pos;
    }

    /** test inequality with another iterator */
    bool operator!=(const ntHashIterator& it) const
    {
        return !(*this == it);
    }

    /** pre-increment operator */
    ntHashIterator& operator++()
    {
        ++m_pos;
        next();
        return *this;
    }

    /** iterator pointing to one past last element */
    static const ntHashIterator end()
    {
        return ntHashIterator();
    }

    /** destructor */
    ~ntHashIterator() {
        if(m_hVec!=NULL)
            delete [] m_hVec;
    }

private:

    /** DNA sequence */
    std::string m_seq;

    /** number of hashes */
    unsigned m_h;

    /** k-mer size */
    unsigned m_k;

    /** hash values */
    uint64_t *m_hVec;

    /** position of current k-mer */
    size_t m_pos;

    /** forward-strand k-mer hash value */
    uint64_t m_fhVal;

    /** reverse-complement k-mer hash value */
    uint64_t m_rhVal;
};

#endif
