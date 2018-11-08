/***************************************************************************
 *  examples/applications/skew3-lcp.cpp
 *
 *  Implementation of the DC3-LCP aka skew3-lcp algorithm in external memory,
 *  based on the internal memory version described in Juha Kaerkkaeinen,
 *  Peter Sanders and Stefan Burkhardt. "Linear work suffix array construction".
 *  Journal of the ACM, Volume 53 Issue 6, November 2006, Pages 918-936.
 *
 *  Part of the STXXL. See http://stxxl.sourceforge.net
 *
 *  Copyright (C) 2004 Jens Mehnert <jmehnert@mpi-sb.mpg.de>
 *  Copyright (C) 2012-2015 Timo Bingmann <tb@panthema.net>
 *  Copyright (C) 2012-2015 Daniel Feist <daniel.feist@student.kit.edu>
 *
 *  Distributed under the Boost Software License, Version 1.0.
 *  (See accompanying file LICENSE_1_0.txt or copy at
 *  http://www.boost.org/LICENSE_1_0.txt)
 **************************************************************************/

#include <cmath>
#include <string>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>
#include <vector>
#include <limits>
#include <numeric>
#include <iterator>
#include <cassert>
#include <algorithm>

#include <stxxl/io>
#include <stxxl/timer>
#include <stxxl/stats>
#include <stxxl/timer>
#include <stxxl/random>
#include <stxxl/sorter>
#include <stxxl/stream>
#include <stxxl/vector>
#include <stxxl/cmdline>
#include <stxxl/sequence>
#include <stxxl/algorithm>
#include <stxxl/bits/common/uint_types.h>

using stxxl::internal_size_type;
using stxxl::external_size_type;

/// Global variables, types and helpers.

// limits max ram used by external data structures to 1.0 GiB the conservative
// way
internal_size_type ram_use = 0.8 * 1024 * 1024 * 1024;

// alphabet data type
typedef unsigned char alphabet_type;

// calculation data type
typedef external_size_type size_type;

/// Suffix array checker for correctness verification.

/**
 * Algorithm to check whether the suffix array is correct. Loosely based on the
 * ideas of Kaerkkaeinen und Burghardt, originally implemented in STXXL by Jens
 * Mehnert (2004), reimplemented using triples by Timo Bingmann (2012).
 *
 * \param InputT is the original text, from which the suffix array was build.
 * \param InputSA is the suffix array from InputT.
 * \return true iff InputSA is the correct suffix array of InputT, else false.
 *
 * Note: ISA := The inverse of SA.
 */
template <typename InputT, typename InputSA>
bool sacheck(InputT& inputT, InputSA& inputSA)
{
    typedef typename InputSA::value_type offset_type;
    typedef stxxl::tuple<offset_type, offset_type> pair_type;
    typedef stxxl::tuple<offset_type, offset_type, offset_type> triple_type;

    /// Pipeline Declaration

    // build tuples with index: (SA[i]) -> (i, SA[i])
    typedef stxxl::stream::counter<offset_type> index_counter_type;
    index_counter_type index_counter;

    typedef stxxl::stream::make_tuple<index_counter_type, InputSA> tuple_index_sa_type;
    tuple_index_sa_type tuple_index_sa(index_counter, inputSA);

    // take (i, SA[i]) and sort to (ISA[i], i)
    typedef stxxl::tuple_less2nd<pair_type> pair_less_type;
    typedef typename stxxl::stream::sort<tuple_index_sa_type, pair_less_type> build_isa_type;

    build_isa_type build_isa(tuple_index_sa, pair_less_type(), ram_use / 3);

    // build (ISA[i], T[i], ISA[i+1]) and sort to (i, T[SA[i]], ISA[SA[i]+1])
    typedef stxxl::tuple_less1st<triple_type> triple_less_type;             // comparison relation

    typedef typename stxxl::stream::use_push<triple_type> triple_push_type; // indicator use push()
    typedef typename stxxl::stream::runs_creator<triple_push_type, triple_less_type> triple_rc_type;
    typedef typename stxxl::stream::runs_merger<typename triple_rc_type::sorted_runs_type, triple_less_type> triple_rm_type;

    triple_rc_type triple_rc(triple_less_type(), ram_use / 3);

    /// Process

    // loop 1: read ISA and check for a permutation.
    // simultaneously create runs of triples by iterating ISA and T.
    size_type totalSize;
    {
        offset_type prev_isa = (*build_isa).first;
        offset_type counter = 0;
        while (!build_isa.empty()) {
            if ((*build_isa).second != counter) {
                std::cout << "Error: suffix array is not a permutation of 0..n-1." << std::endl;
                return false;
            }

            ++counter;
            ++build_isa; // ISA is one in front of T

            if (!build_isa.empty()) {
                triple_rc.push(triple_type(prev_isa, *inputT, (*build_isa).first));
                prev_isa = (*build_isa).first;
            }
            ++inputT;
        }

        totalSize = counter;
    }

    if (totalSize == 1) return true;

    // loop 2: read triples (i,T[SA[i]],ISA[SA[i]+1]) and check for correct ordering.
    triple_rm_type triple_rm(triple_rc.result(), triple_less_type(), ram_use / 3);

    {
        triple_type prev_triple = *triple_rm;
        size_type counter = 0;

        ++triple_rm;

        while (!triple_rm.empty()) {
            const triple_type& this_triple = *triple_rm;

            if (prev_triple.second > this_triple.second)
            {
                // simple check of first character of suffix
                std::cout << "Error: suffix array position "
                          << counter << " ordered incorrectly." << std::endl;
                return false;
            }
            else if (prev_triple.second == this_triple.second) {
                if (this_triple.third == (offset_type)totalSize) {
                    // last suffix of string must be first among those with same first character
                    std::cout << "Error: suffix array position "
                              << counter << " ordered incorrectly." << std::endl;
                    return false;
                }
                if (prev_triple.third != (offset_type)totalSize &&
                    prev_triple.third > this_triple.third) {
                    // positions SA[i] and SA[i-1] has same first character but their suffixes are ordered incorrectly:
                    // the suffix position of SA[i] is given by ISA[SA[i]]
                    std::cout << "Error: suffix array position "
                              << counter << " ordered incorrectly." << std::endl;
                    return false;
                }
            }

            prev_triple = this_triple;

            ++triple_rm;
            ++counter;
        }
    }
    return true;
}

/// Wrapper for suffix array checker class sachecker.

/**
 * Wrapper function for sachecker.
 *
 * \param InputT is the original text, from which the suffix array was built.
 * \param InputSA is the suffix array from InputT.
 * \return true iff InputSA is the correct suffix array of InputT, else false.
 */
template <typename InputT, typename InputSA>
bool sacheck_vectors(InputT& inputT, InputSA& inputSA)
{
    typename stxxl::stream::streamify_traits<typename InputT::iterator>::stream_type streamT
        = stxxl::stream::streamify(inputT.begin(), inputT.end());

    typename stxxl::stream::streamify_traits<typename InputSA::iterator>::stream_type streamSA
        = stxxl::stream::streamify(inputSA.begin(), inputSA.end());

    return sacheck(streamT, streamSA);
}

/// Kasai's semi external lcp array calculator.

template <typename InputStream>
struct sa_index_stream_type
{
    typedef typename InputStream::value_type offset_type;
    typedef stxxl::tuple<offset_type, offset_type, offset_type> value_type;

private:
    size_type m_counter;

    InputStream& m_input;

    value_type m_curr;

public:
    sa_index_stream_type(InputStream& input) : m_counter(0), m_input(input)
    {
        if (!m_input.empty())
            m_curr = value_type(*m_input, m_counter++, 0);
    }

    const value_type& operator * () const
    {
        return m_curr;
    }

    sa_index_stream_type& operator ++ ()
    {
        ++m_input;
        if (!m_input.empty())
            m_curr = value_type(*m_input, m_counter++, m_curr.first);
        return *this;
    }

    bool empty() const
    {
        return m_input.empty();
    }
};

/**
 * Algorithm to calculate the LCP array from a string and suffix array in linear
 * time.  Based on the ideas of Kasai et al. (2001), implemented by Timo
 * Bingmann (2012).
 *
 * \param &string is a reference to the input sequence container.
 * \param &SA is a reference to the suffix array container.
 * \param &lcp is a reference to the container which stores the lcp array
 * calculated by the algorithm.
 */
template <typename StringContainer, typename SAContainer, typename LCPContainer>
void lcparray_stxxl_kasai(const StringContainer& string,
                          const SAContainer& SA, LCPContainer& lcp)
{
    assert(string.size() == SA.size());

    /// Generate ISA.

    typedef typename StringContainer::value_type alphabet_type;
    typedef typename SAContainer::value_type offset_type;
    typedef typename SAContainer::size_type size_type;

    static const std::size_t block_size = sizeof(offset_type) * 1024 * 1024 / 2;

    typedef stxxl::tuple<offset_type, offset_type> offset_pair_type;
    typedef stxxl::tuple<offset_type, offset_type, offset_type> offset_triple_type;

    // define stream iterating over a STXXL vector
    typedef stxxl::stream::iterator2stream<typename SAContainer::const_iterator> sa_stream_type;
    sa_stream_type sa_stream(SA.begin(), SA.end());
    sa_index_stream_type<sa_stream_type> sa_index_stream(sa_stream);

    // comparator for ISA sort: 1st component ascending
    typedef stxxl::tuple_less1st<offset_triple_type> offset_triple_less1st_type;
    offset_triple_less1st_type offset_triple_less1st;

    // runs creator for ISA stream
    typedef stxxl::stream::sort<sa_index_stream_type<sa_stream_type>, offset_triple_less1st_type, block_size> isa_sort_type;
    isa_sort_type isa_sort(sa_index_stream, offset_triple_less1st, ram_use / 8, ram_use / 8);

    /// Output sorter.

    typedef stxxl::tuple_less1st_less2nd<offset_pair_type> offset_pair_less1st_type;
    offset_pair_less1st_type offset_pair_less1st;

    typedef stxxl::sorter<offset_pair_type, offset_pair_less1st_type, block_size> lcp_sorter_type;
    lcp_sorter_type lcp_sorter(offset_pair_less1st, ram_use / 8);

    /// Kasai algorithm: iterate over ISA.

    {
        size_type h = 0; // current height
        size_type i = 0; // ISA index counter;

        lcp_sorter.push(offset_pair_type(0, 0));

        std::vector<alphabet_type> stringRAM(string.begin(), string.end());

        while (!isa_sort.empty()) {
            offset_type k = isa_sort->second;  // k = ISA[i]

            if (k > offset_type(0)) {
                size_type j = isa_sort->third; // j = SA[k-1];

                while (i + h < stringRAM.size() && j + h < stringRAM.size() &&
                       stringRAM[i + h] == stringRAM[j + h])
                    h++;

                lcp_sorter.push(offset_pair_type(k, h));
            }
            if (h > 0) h--;

            ++isa_sort, ++i;
        }
    }

    /// Collect output.

    stxxl::stream::choose<lcp_sorter_type, 2> pair_2nd_stream(lcp_sorter);

    lcp.resize(string.size());
    lcp_sorter.sort();

    stxxl::stream::materialize(pair_2nd_stream, lcp.begin(), lcp.end());
}

namespace algo {

/*
 * Class which implements the DC3-LCP aka skew3-lcp algorithm.
 * T := input string. The recursion works as follows:
 * Step 1: a) convert T into stream of quadruples (i,T[i],T[i+1],T[i+2],T[i+3]) (-> make_quads class)
 *         b) pick all mod1/mod2 triples T[i],T[i+1],T[i+2] at position i mod 3 != 0 (-> extract_mod12 class)
 *         b) sort mod1/mod2 triples lexicographically (-> build_sa class)
 *         c) give mod1/mod2 triples lexicographical ascending names and
 *            save overlap of neighbored triples in LCPn (-> naming class)
 *         d) check lexicographical names for uniqueness (-> naming class)
 *            If yes: proceed with Step 2, If no: set T' := "lexicographical names" and recurse
              (i.e. run Step 1 again) which returns SA12 := SA(T') and LCP12 := LCPn(T')
 * Step 2: a) by sorting the lexicographical names n we receive ranks r which represent ISA.
 *         b) construct mod0-quints, mod1-quads and mod2-quints  (-> build_sa class)
 *         c) prepare for merging by:
 *            sort mod0-quints by 2 components, sort mod1-quads / mod2-quints by one component (-> build_sa class)
 *         c) merge {mod0-quints, mod1-quads, mod2-quints} and
 *            compare first characters of tuples,
 *            perform an RMQ on LCP12,
 *            perform an RMQ on LCPn (-> merge_sa_lcp class)
 * Step 3: a) return Suffix and LCP array of T
 *
 * \param offset_type later suffix array data type
 */
template <typename offset_type>
class skew
{
public:
    /// Types and comparators needed by the skew algorithm.

    struct intPair
    {
        offset_type id, i, j;

        intPair() { }

        intPair(offset_type _id, offset_type _i, offset_type _j)
            : id(_id), i(_i), j(_j) { }
    };

    struct leftRmq
    {
        offset_type id, i, j, k;

        leftRmq() { }

        leftRmq(offset_type _id, offset_type _i, offset_type _j, offset_type _k)
            : id(_id), i(_i), j(_j), k(_k) { }
    };

    struct leftRmqComparator
    {
        bool operator () (const leftRmq& a, const leftRmq& b) const
        {
            return a.k < b.k;
        }

        leftRmq min_value() const
        {
            return leftRmq(std::numeric_limits<offset_type>::min(),
                           std::numeric_limits<offset_type>::min(),
                           std::numeric_limits<offset_type>::min(),
                           std::numeric_limits<offset_type>::min());
        }

        leftRmq max_value() const
        {
            return leftRmq(std::numeric_limits<offset_type>::max(),
                           std::numeric_limits<offset_type>::max(),
                           std::numeric_limits<offset_type>::max(),
                           std::numeric_limits<offset_type>::max());
        }
    };

    struct rightRmq
    {
        offset_type id, i, j, k, tmpK;

        rightRmq() { }

        rightRmq(offset_type _id, offset_type _i, offset_type _j, offset_type _k, offset_type _tmpK)
            : id(_id), i(_i), j(_j), k(_k), tmpK(_tmpK) { }
    };

    struct rightRmqComparator
    {
        bool operator () (const rightRmq& a, const rightRmq& b) const
        {
            return a.k < b.k;
        }

        rightRmq min_value() const
        {
            return rightRmq(std::numeric_limits<offset_type>::min(),
                            std::numeric_limits<offset_type>::min(),
                            std::numeric_limits<offset_type>::min(),
                            std::numeric_limits<offset_type>::min(),
                            std::numeric_limits<offset_type>::min());
        }

        rightRmq max_value() const
        {
            return rightRmq(std::numeric_limits<offset_type>::max(),
                            std::numeric_limits<offset_type>::max(),
                            std::numeric_limits<offset_type>::max(),
                            std::numeric_limits<offset_type>::max(),
                            std::numeric_limits<offset_type>::max());
        }
    };

    struct finalPair
    {
        offset_type id, min;

        finalPair() { }

        finalPair(offset_type _id, offset_type _min)
            : id(_id), min(_min) { }
    };

    struct finalPairComparator
    {
        bool operator () (const finalPair& a, const finalPair& b) const
        {
            return a.id < b.id;
        }

        finalPair min_value() const
        {
            return finalPair(std::numeric_limits<offset_type>::min(),
                             std::numeric_limits<offset_type>::min());
        }

        finalPair max_value() const
        {
            return finalPair(std::numeric_limits<offset_type>::max(),
                             std::numeric_limits<offset_type>::max());
        }
    };

    struct two_tuple
    {
        offset_type first, second;

        two_tuple() { }

        two_tuple(offset_type _first, offset_type _second)
            : first(_first), second(_second) { }
    };

    struct l3Tuple
    {
        offset_type first;
        bool second;
        offset_type third, fourth;

        l3Tuple() { }

        l3Tuple(offset_type _first, bool _second, offset_type _third, offset_type _fourth)
            : first(_first), second(_second), third(_third), fourth(_fourth) { }
    };

    struct l3TupleComparator
    {
        bool operator () (const l3Tuple& a, const l3Tuple& b) const
        {
            return a.third < b.third;
        }

        l3Tuple min_value() const
        {
            return l3Tuple(std::numeric_limits<offset_type>::min(),
                           std::numeric_limits<bool>::min(),
                           std::numeric_limits<offset_type>::min(),
                           std::numeric_limits<offset_type>::min());
        }

        l3Tuple max_value() const
        {
            return l3Tuple(std::numeric_limits<offset_type>::max(),
                           std::numeric_limits<bool>::max(),
                           std::numeric_limits<offset_type>::max(),
                           std::numeric_limits<offset_type>::max());
        }
    };

    struct innerTuple
    {
        offset_type first;
        bool second;
        offset_type third;

        innerTuple() { }

        innerTuple(offset_type _first, bool _second, offset_type _third)
            : first(_first), second(_second), third(_third) { }
    };

    struct innerTupleComparator
    {
        bool operator () (const innerTuple& x, const innerTuple& y) const
        {
            return x.third < y.third;
        }

        innerTuple min_value() const
        {
            return innerTuple(std::numeric_limits<offset_type>::min(),
                              std::numeric_limits<bool>::min(),
                              std::numeric_limits<offset_type>::min());
        }

        innerTuple max_value() const
        {
            return innerTuple(std::numeric_limits<offset_type>::max(),
                              std::numeric_limits<bool>::max(),
                              std::numeric_limits<offset_type>::max());
        }
    };

    struct pos_pair
    {
        offset_type first;
        bool second;
        offset_type third;

        pos_pair() { }

        pos_pair(offset_type _first, bool _second, offset_type _third)
            : first(_first), second(_second), third(_third) { }
    };

    struct PairComparator
    {
        bool operator () (const pos_pair& x, const pos_pair& y) const
        {
            return x.first < y.first;
        }

        pos_pair min_value() const
        {
            return pos_pair(std::numeric_limits<offset_type>::min(),
                            std::numeric_limits<bool>::min(),
                            std::numeric_limits<offset_type>::min());
        }

        pos_pair max_value() const
        {
            return pos_pair(std::numeric_limits<offset_type>::max(),
                            std::numeric_limits<bool>::max(),
                            std::numeric_limits<offset_type>::max());
        }
    };

    // 2-tuple, 3-tuple, 4-tuple (=quads), 5-tuple(=quints) definition
    typedef stxxl::tuple<offset_type, offset_type> skew_pair_type;
    typedef stxxl::tuple<offset_type, offset_type, offset_type> skew_triple_type;
    typedef stxxl::tuple<offset_type, offset_type, offset_type, offset_type> skew_quad_type;
    typedef stxxl::tuple<offset_type, offset_type, offset_type, offset_type, offset_type> skew_quint_type;

    static const std::size_t block_size = sizeof(offset_type) * 1024 * 1024 / 2;

    typedef typename stxxl::VECTOR_GENERATOR<offset_type, 1, 2, block_size>::result offset_array_type;
    typedef stxxl::stream::vector_iterator2stream<typename offset_array_type::iterator> offset_array_it_rg;

    typedef stxxl::sequence<offset_type, block_size> simpleDeq;
    typedef stxxl::sequence<intPair, block_size> pairDeq;
    typedef stxxl::sequence<two_tuple, block_size> twoDeq;

    /**
     * Comparison function for the mod0 tuples.
     */
    struct less_mod0
    {
        typedef skew_quint_type value_type;

        bool operator () (const value_type& a, const value_type& b) const
        {
            if (a.second == b.second)
                return a.fourth < b.fourth;
            else
                return a.second < b.second;
        }

        static value_type min_value() { return value_type::min_value(); }
        static value_type max_value() { return value_type::max_value(); }
    };

    /**
     * Counter for creating tuple indexes for example.
     */
    template <class ValueType>
    struct counter
    {
        typedef ValueType value_type;

        value_type cnt;

        counter() : cnt(0) { }

        const value_type& operator * () const
        {
            return cnt;
        }

        counter& operator ++ ()
        {
            ++cnt;
            return *this;
        }

        bool empty() const
        {
            return false;
        }
    };

    /**
     * Compares quadruples with respect to the first component.
     */
    struct less_quad_2nd
    {
        typedef skew_quad_type value_type;

        bool operator () (const skew_quad_type& a, const skew_quad_type& b) const
        {
            return a.second < b.second;
        }

        static value_type min_value() { return value_type::min_value(); }
        static value_type max_value() { return value_type::max_value(); }
    };

    /**
     * Compares five-tuples with respect to the second component.
     */
    struct less_quint_2nd
    {
        typedef skew_quint_type value_type;

        bool operator () (const value_type& a, const value_type& b) const
        {
            return a.second < b.second;
        }

        static value_type min_value() { return value_type::min_value(); }
        static value_type max_value() { return value_type::max_value(); }
    };

    /**
     * Put the (0 mod 2) [which are the 1,2 mod 3 tuples] tuples at the begin.
     */
    struct less_skew
    {
        typedef skew_pair_type value_type;

        bool operator () (const value_type& a, const value_type& b) const
        {
            if ((a.first & 1) == (b.first & 1))
                return a.first < b.first;
            else
                return (a.first & 1) < (b.first & 1);
        }

        static value_type min_value() { return value_type::min_value(); }
        static value_type max_value() { return value_type::max_value(); }
    };

    /**
     * Compare two pairs by their first component.
     */
    struct mod12Comparator
    {
        typedef skew_pair_type value_type;

        bool operator () (const value_type& a, const value_type& b) const
        {
            return a.first < b.first;
        }

        value_type min_value() const
        {
            return value_type(std::numeric_limits<offset_type>::min(), std::numeric_limits<offset_type>::min());
        }

        value_type max_value() const
        {
            return value_type(std::numeric_limits<offset_type>::max(), std::numeric_limits<offset_type>::max());
        }
    };

    /**
     * Concatenates two streams as streamA '+' streamB containing pairs.
     */
    template <class StreamType, class StreamA, class StreamB>
    class concat
    {
    public:
        typedef StreamType value_type;

    private:
        StreamA& A;
        StreamB& B;

    public:
        concat(StreamA& A_, StreamB& B_) : A(A_), B(B_)
        {
            assert(!A.empty());
            assert(!B.empty());
        }

        const value_type& operator * () const
        {
            if (empty())
                throw std::runtime_error("operator()* error");
            if (!A.empty())
                return *A;
            return *B;
        }

        concat& operator ++ ()
        {
            assert(!empty());

            if (!A.empty())
                ++A;
            else if (!B.empty())
                ++B;

            return *this;
        }

        bool empty() const
        {
            return ((A.empty()) && (B.empty()));
        }
    };

    /**
     * Sort skew_quad datatype.
     */
    template <typename alphabet_type>
    struct less_quad
    {
        typedef typename stxxl::tuple<offset_type, alphabet_type, alphabet_type, alphabet_type> value_type;

        bool operator () (const value_type& a, const value_type& b) const
        {
            if (a.second == b.second) {
                if (a.third == b.third)
                    return a.fourth < b.fourth;
                else
                    return a.third < b.third;
            }
            else {
                return a.second < b.second;
            }
        }

        static value_type min_value() { return value_type::min_value(); }
        static value_type max_value() { return value_type::max_value(); }
    };

    /// Skew naming.

    /**
     * Check, if last three components of two quads are equal.
     */
    template <class quad_type>
    static inline bool quad_eq(const quad_type& a, const quad_type& b)
    {
        return (a.second == b.second) && (a.third == b.third) && (a.fourth == b.fourth);
    }

    /**
     * Check overlapping of two triples in case of non-equality.
     */
    template <class quad_type>
    static inline char quad_neq(const quad_type& a, const quad_type& b)
    {
        if ((a.second == b.second) && (a.third == b.third))
            return 2;
        else if (a.second == b.second)
            return 1;
        else
            return 0;
    }

    /**
     * Naming pipe for the conventional skew algorithm without discarding.
     */
    template <class Input>
    class naming
    {
    public:
        typedef typename Input::value_type quad_type;
        typedef skew_pair_type value_type;

    private:
        Input& A;
        bool& unique;
        offset_type lexname;
        quad_type prev;
        skew_pair_type result;

        simpleDeq* LCPn;

    public:
        naming(Input& A_, bool& unique_, simpleDeq* lcp_names) : A(A_), unique(unique_), lexname(0)
        {
            assert(!A.empty());
            unique = true;

            LCPn = lcp_names;
            prev = *A;
            LCPn->push_back(0); // insert zero sentinel to hold LCPn[0] = 0

            result.first = prev.first;
            result.second = lexname;
        }

        const value_type& operator * () const
        {
            return result;
        }

        naming& operator ++ ()
        {
            assert(!A.empty());

            ++A;
            if (A.empty()) return *this;

            quad_type curr = *A;

            if (!quad_eq(prev, curr)) {
                ++lexname;
                LCPn->push_back(quad_neq(prev, curr));
            }
            else {
                if (!A.empty() && curr.second != offset_type(0)) {
                    LCPn->push_back(3);
                    unique = false;
                }
            }

            result.first = curr.first;
            result.second = lexname;

            prev = curr;
            return *this;
        }

        bool empty() const
        {
            return A.empty();
        }
    };

    /// Pipelining classes.

    /**
     * Create tuples until one of the input streams are empty.
     */
    template <class InputA, class InputB, const int add_alphabet = 0>
    class make_pairs
    {
    public:
        typedef skew_pair_type value_type;

    private:
        InputA& A;
        InputB& B;
        value_type result;

    public:
        make_pairs(InputA& a, InputB& b)
            : A(a), B(b)
        {
            assert(!A.empty());
            assert(!B.empty());
            if (!empty())
                result = value_type(*A, *B + add_alphabet);
        }

        const value_type& operator * () const
        {
            return result;
        }

        make_pairs& operator ++ ()
        {
            assert(!A.empty());
            assert(!B.empty());

            ++A;
            ++B;

            if (!A.empty() && !B.empty())
                result = value_type(*A, *B + add_alphabet);

            return *this;
        }

        bool empty() const
        {
            return (A.empty() || B.empty());
        }
    };

    /**
     * Collect three characters t_i, t_{i+1}, t_{i+2} beginning at the index i.
     * Since we need at least one unique endcaracter, we free the first characters,
     * i.e. we map (t_i) -> (i,t_i,t_{i+1},t_{i+2}).
     *
     * \param Input holds all characters t_i from input string t.
     * \param alphabet_type
     * \param add_alphabet
     */
    template <class Input, typename alphabet_type, const int add_alphabet = 0>
    class make_quads
    {
    public:
        typedef stxxl::tuple<offset_type, alphabet_type, alphabet_type, alphabet_type> value_type;

    private:
        Input& A;
        value_type current;
        offset_type counter;
        unsigned int z3z; // = counter mod 3, ("+",Z/3Z) is cheaper than '%'
        bool finished;

        offset_array_type& backup;

    public:
        make_quads(Input& data_in_, offset_array_type& backup_)
            : A(data_in_),
              current(0, 0, 0, 0),
              counter(0),
              z3z(0),
              finished(false),
              backup(backup_)
        {
            assert(!A.empty());

            current.first = counter;
            current.second = (*A).second + add_alphabet;
            ++A;

            if (!A.empty()) {
                current.third = (*A).second + add_alphabet;
                ++A;
            }
            else {
                current.third = 0;
                current.fourth = 0;
            }

            if (!A.empty())
                current.fourth = (*A).second + add_alphabet;
            else
                current.fourth = 0;
        }

        const value_type& operator * () const
        {
            return current;
        }

        make_quads& operator ++ ()
        {
            assert(!A.empty() || !finished);

            if (current.second != offset_type(0))
                backup.push_back(current.second);

            if (++z3z == 3) z3z = 0;

            current.first = ++counter;
            current.second = current.third;
            current.third = current.fourth;

            if (!A.empty())
                ++A;

            if (!A.empty())
                current.fourth = (*A).second + add_alphabet;
            else
                current.fourth = 0;

            if ((current.second == offset_type(0)) && (z3z != 1))
                finished = true;

            return *this;
        }

        bool empty() const
        {
            return (A.empty() && finished);
        }
    };

    /**
     * Drop 1/3 of the input, more exactly the offsets at positions (0 mod 3).
     * Index begins with 0.
     */
    template <class Input>
    class extract_mod12
    {
    public:
        typedef typename Input::value_type value_type;

    private:
        Input& A;
        offset_type counter;
        offset_type output_counter;
        value_type result;

    public:
        extract_mod12(Input& A_) : A(A_), counter(0), output_counter(0)
        {
            assert(!A.empty());
            ++A, ++counter; // skip 0 = mod0 offset
            if (!A.empty()) {
                result = *A;
                result.first = output_counter;
            }
        }

        const value_type& operator * () const
        {
            return result;
        }

        extract_mod12& operator ++ ()
        {
            assert(!A.empty());

            ++A, ++counter, ++output_counter;

            if (!A.empty() && (counter % 3) == 0)
                ++A, ++counter; // skip mod0 offsets

            if (!A.empty()) {
                result = *A;
                result.first = output_counter;
            }

            return *this;
        }

        bool empty() const
        {
            return A.empty();
        }
    };

    /// Basis for (Batched) Range Minimum Queries (RMQ) construction.

    /**
     * Sparse Table construction based on "Bender, M. A., and Farach-Colton, M.
     * The LCA problem revisited (2000)", that creates a table of size n log n for
     * every minimum in 'power-of-two'-ranges and allows to answer an RMQ(i,j)
     * in O(1) by simple table lookups afterwards.
     */
    class sparseTable
    {
    private:
        std::vector<std::vector<offset_type> > QTable;

    public:
        sparseTable(std::vector<offset_type>& fieldVector)
        {
            std::size_t dimI = fieldVector.size();    // y-Dimension of QTable
            std::size_t dimK = floor(log2(dimI)) + 1; // x-Dimension of QTable
            createQTable(fieldVector, dimI, dimK);
        }

        void createQTable(std::vector<offset_type>& field, std::size_t dimI, std::size_t dimK)
        {
            assert(dimI > 0);
            assert(dimK > 0);

            std::size_t iLimit = dimI - 1;
            std::size_t kLimit = dimK - 1;
            std::size_t magicCounter = 0;

            QTable.push_back(std::vector<offset_type>());
            std::vector<offset_type>& tempVector = QTable[QTable.size() - 1];
            tempVector.reserve(iLimit);

            for (std::size_t i = 0; i <= iLimit; ++i)
                tempVector.push_back(field[i]);

            for (std::size_t k = 0; k < kLimit; ++k) {
                QTable.push_back(std::vector<offset_type>());
                std::vector<offset_type>& tempVector = QTable[QTable.size() - 1];
                tempVector.reserve(dimI - (0 + (1 << k)));

                magicCounter += (1 << k);
                for (std::size_t i = 0; i + (1 << k) < dimI; ++i) {
                    if (QTable[k][i] <= QTable[k][i + (1 << k)])
                        tempVector.push_back(QTable[k][i]);
                    else
                        tempVector.push_back(QTable[k][i + (1 << k)]);

                    if (i == (iLimit - magicCounter))
                        break;
                }
            }
        }

        offset_type query(offset_type I, offset_type J)
        {
            std::size_t i = I;
            std::size_t j = J;
            std::size_t k = floor(log2(j - i + 1));
            offset_type r = QTable[k][i];
            offset_type s = QTable[k][j - (1 << k) + 1];

            if (r <= s) return r;
            else return s;
        }

        void clear()
        {
            QTable.clear();
        }
    };

    /**
     * This class uses the sparseTable construction to answer batched RMQs.
     * First: calculate number_of_pages and page borders.
     * Second: split up sequence of RMQs (id, i, j) into two parts based
     * on page borders and sort splitted RMQs.
     * Third: build up minBorderArray.
     * Fourth: sort and answer RMQs with respect to their id.
     *
     * \param Stream1 holds a reference to a sequence of values the RMQs are performed on.
     * \param Stream2 holds a reference to a sequence of RMQs.
     */
    template <class Stream1, class Stream2>
    class sparseTableAlgo
    {
        typedef leftRmqComparator left_cmp;
        typedef rightRmqComparator right_cmp;
        typedef finalPairComparator final_cmp;

        typedef stxxl::sorter<leftRmq, left_cmp, block_size> left_sorter_type;
        typedef stxxl::sorter<rightRmq, right_cmp, block_size> right_sorter_type;
        typedef stxxl::sorter<finalPair, final_cmp, block_size> final_sorter_type;

        left_sorter_type left_sorter;
        right_sorter_type right_sorter;
        final_sorter_type collective_sorter;

    private:
        std::vector<offset_type> minBorderArray;
        offset_type span, end_part, numberOfPages;

        // Compare min of page k with last minBorderArray entries
        offset_type getArrayMin(offset_type pageMin, offset_type k, offset_type tmpK)
        {
            offset_type min = pageMin;
            offset_type lastPages = k - 1 - tmpK;

            if (lastPages > offset_type(0)) {
                for (offset_type p = 1; p <= lastPages; ++p) {
                    if (min > minBorderArray[k - p])
                        min = minBorderArray[k - p];
                }
            }
            return min;
        }

    public:
        sparseTableAlgo()
            : left_sorter(left_cmp(), ram_use / 6),
              right_sorter(right_cmp(), ram_use / 6),
              collective_sorter(final_cmp(), ram_use / 6) { }

        // compute numberOfPages and their page-size (span)
        void initialize(offset_type lengthOfVector)
        {
            offset_type M = (offset_type(ram_use / 2) / sizeof(offset_type));
            span = end_part = 0;

            assert(lengthOfVector > offset_type(0));

            if (lengthOfVector >= offset_type((ram_use / 2) / sizeof(offset_type)))
                span = 10000;

            do {
                // space inequation: M - c >= v + v*log_2(v) + n/v
                if ((M - (span + 1) - ((span + 1) * log2(span + 1)) - (lengthOfVector / (span + 1))) > 0)
                    ++span;
                else
                    break;
            } while (span < lengthOfVector);

            assert(span > offset_type(0));

            offset_type tmp = (offset_type)lengthOfVector / span;
            if ((lengthOfVector % span) != 0) {
                end_part = lengthOfVector - (tmp * span);
                numberOfPages = tmp + 1;
            }
            else {
                numberOfPages = tmp;
            }
        }

        // split up RMQs (id, i, j)) into two parts
        void splitAndSortRMQ(Stream2& pair_stream)
        {
            intPair temp;
            offset_type left, leftBorder, right, rightBorder, id, tmpK;
            bool finishedI;

            while (!pair_stream.empty()) {
                temp = *pair_stream;
                id = temp.id;
                left = temp.i;
                right = temp.j;

                finishedI = 0;
                tmpK = 0;
                leftBorder = 0;
                rightBorder = span - 1;

                // i := left, j := right, | := left or right border, *|* := currently considered border

                for (offset_type k = 0; k < numberOfPages; ++k) {
                    // find beginning of RMQ (id, i, j)
                    if (left >= leftBorder && left <= rightBorder && finishedI == 0) {
                        // (trivial) case: ...*|*..i..j..*|*...
                        // => push leftRMQ(id, relative_i, relative_j, k) into leftRMQSorter
                        if (right >= leftBorder && right <= rightBorder) {
                            left_sorter.push(leftRmq(id, (left - leftBorder), (right - leftBorder), k));
                            finishedI = 1;
                            break;
                        }
                        else {   // case: ...|..i..*|*... => push leftRMQ(id, relative_i, relative_j, k) into leftRMQSorter
                            left_sorter.push(leftRmq(id, (left - leftBorder), (rightBorder - leftBorder), k));
                            finishedI = 1;
                            tmpK = k;

                            assert(numberOfPages > offset_type(1));

                            leftBorder = (rightBorder + 1);
                            if ((k == offset_type(numberOfPages - 2)) && (end_part > offset_type(0)))
                                rightBorder += end_part;
                            else
                                rightBorder += span;

                            continue;
                        }
                    }

                    // case: *|*.....*|*..j..|...
                    // => "empty" page, do nothing and keep on walking

                    // case: *|*..j..|.....|.....|
                    // => push rightRMQ(id, relative_i, relative_j, k, tmpK) into rightRMQSorter
                    if (right >= leftBorder && right <= rightBorder && finishedI == 1) {
                        right_sorter.push(rightRmq(id, 0, (right - leftBorder), k, tmpK));
                        break;
                    }

                    assert(numberOfPages > offset_type(1));

                    leftBorder = (rightBorder + 1);
                    if ((k == offset_type(numberOfPages - 2)) && (end_part > offset_type(0)))
                        rightBorder += end_part;
                    else
                        rightBorder += span;
                }
                ++pair_stream;
            }
            left_sorter.sort();
            right_sorter.sort();
        }

        void answerTuple(Stream1& int_stream)
        {
            offset_type temp, min, left_min, right_min, leftBorder, rightBorder, page_min;
            std::vector<offset_type> tempVector;

            leftBorder = 0;
            rightBorder = span - 1;

            // build up minBorderArray
            for (offset_type n = 0; n < numberOfPages; ++n) {
                min = std::numeric_limits<offset_type>::max();

                for (offset_type i = leftBorder; i <= rightBorder; ++i) {
                    assert(!int_stream.empty());

                    temp = *int_stream;
                    tempVector.push_back(temp);

                    if (temp < min)
                        min = temp;

                    if (i <= rightBorder) {
                        assert(!int_stream.empty());
                        ++int_stream;
                    }
                }

                minBorderArray.push_back(min);

                // create sparseoffset_typeable of tempVector
                sparseTable myTable(tempVector);

                while (!left_sorter.empty() && left_sorter->k <= n) {
                    // getMin(i, j, k, tmpK, vector)
                    left_min = myTable.query(left_sorter->i, left_sorter->j);

                    collective_sorter.push(finalPair(left_sorter->id, left_min));
                    ++left_sorter;
                }

                while (!right_sorter.empty() && right_sorter->k <= n) {
                    // getMin(i, j, k, tmpK, vector)
                    page_min = myTable.query(right_sorter->i, right_sorter->j);

                    // mind previous minima in minBorderArray
                    right_min = getArrayMin(page_min, right_sorter->k, right_sorter->tmpK);

                    collective_sorter.push(finalPair(right_sorter->id, right_min));
                    ++right_sorter;
                }

                tempVector.clear();
                myTable.clear();

                leftBorder = (rightBorder + 1);
                if ((n == offset_type(numberOfPages - 2)) && (end_part > offset_type(0)))
                    rightBorder += end_part;
                else
                    rightBorder += span;
            }
            left_sorter.finish_clear();
            right_sorter.finish_clear();
        }

        void answerRMQ(simpleDeq& resultDeq)
        {
            finalPair tmp;
            offset_type id, min;

            collective_sorter.sort(); // sort by their id's

            while (!collective_sorter.empty()) {
                id = collective_sorter->id;
                min = collective_sorter->min;

                if (collective_sorter.size() > 1)
                    ++collective_sorter;

                if (id == collective_sorter->id) { // equal id's
                    if (min >= collective_sorter->min)
                        resultDeq.push_back(collective_sorter->min);
                    else
                        resultDeq.push_back(min);

                    if (collective_sorter.size() >= 1)
                        ++collective_sorter;
                }
                else {   // unequal id's
                    resultDeq.push_back(min);
                }
            }
            collective_sorter.finish_clear();
            minBorderArray.clear();
        }
    };

    /**
     * Construct the lcp array by:
     * \param lcpn defined as the overlapping of two lexicographical adjacent triples.
     * \param lcp12 defined as the lcp array of the last recursive step of dc3-lcp.
     * \param ISA1('old') defined as the inverse suffix array of the i%3==1 triples of the last recursive step of dc3.
     * \param ISA2('old') defined as the inverse suffix array of the i%3==2 triples of the last recursive step of dc3.
     * \param SA12('old') defined as the suffix array of the last recursive step of dc3.
     */
    class build_lcp
    {
    private:
        typedef l3Tuple l3_type;
        typedef innerTuple inner_type;
        typedef two_tuple two_tuple_type;
        typedef pos_pair pos_pair_type;

        typedef l3TupleComparator l3Tuple_cmp;
        typedef innerTupleComparator innerTuple_cmp;
        typedef PairComparator pair_cmp;

        typedef stxxl::sorter<l3Tuple, l3Tuple_cmp, block_size> l3_tuple_type;
        typedef stxxl::sorter<innerTuple, innerTuple_cmp, block_size> inner_tuple_type;
        typedef stxxl::sorter<pos_pair, pair_cmp, block_size> pair_type;

        typedef typename simpleDeq::stream simpleDeqStream;
        typedef typename pairDeq::stream pairDeqStream;
        typedef typename twoDeq::stream twoDeqStream;

        //typedef single_concat<simpleDeqStream, simpleDeqStream> s_concatenation;
        typedef concat<offset_type, simpleDeqStream, simpleDeqStream> s_concatenation;

        simpleDeq l1Deq, l2Deq, l3Deq;
        simpleDeq* lcpn_ptr, * isa1_ptr, * isa2_ptr, * sa12_ptr;
        simpleDeqStream* l1_ptr, * l2_ptr, * l3_ptr;

        pairDeq l2RMQ;
        twoDeq l3TempDeq;

        build_lcp* lcp12_ptr;

        bool finished_l2, blank;

        offset_type countl2;
        offset_type result;

    public:
        build_lcp(simpleDeq* lcpn, build_lcp* lcp12, simpleDeq* ISA1, simpleDeq* ISA2, simpleDeq* SA_12)
            : finished_l2(false), blank(false), countl2(0)
        {
            lcpn_ptr = lcpn;
            lcp12_ptr = lcp12;
            isa1_ptr = ISA1;
            isa2_ptr = ISA2;
            sa12_ptr = SA_12;

            l1_ptr = NULL;
        }

        void saveL1(char l)
        {
            l1Deq.push_back(l);
        }

        void saveRMQ(offset_type r_i, offset_type r_j)
        {
            if (lcp12_ptr != NULL) {
                if (r_j > offset_type(0))
                    r_j = r_j - offset_type(1);

                assert(lcpn_ptr->size() >= r_i);
                assert(lcpn_ptr->size() >= r_j);

                if (r_i <= r_j)
                    l2RMQ.push_back(intPair(countl2, r_i, r_j));
                else
                    l2RMQ.push_back(intPair(countl2, r_j, r_i));

                ++countl2;
            }
        }

        void finalize()
        {
            if (l1_ptr != NULL) return;

            l1_ptr = new simpleDeqStream(l1Deq);

            if (lcp12_ptr != NULL) {
                pairDeqStream rmqStream = l2RMQ.get_stream();

                sparseTableAlgo<build_lcp, pairDeqStream> sp_algo;
                sp_algo.initialize(lcpn_ptr->size()); // |lcpn| == |lcp12|
                sp_algo.splitAndSortRMQ(rmqStream);
                sp_algo.answerTuple(*lcp12_ptr);
                sp_algo.answerRMQ(l2Deq);
            }
            delete lcp12_ptr;
            lcp12_ptr = NULL;

            l2_ptr = new simpleDeqStream(l2Deq);

            finished_l2 = true;
            answerL3();
        }

        void preprocessL3(offset_type r_i, offset_type r_j)
        {
            l3TempDeq.push_back(two_tuple(r_i, r_j));
        }

        void answerL3()
        {
            assert(finished_l2);

            pairDeq rmqDeq;

            if (l2Deq.size() == 0) {
                offset_type id = 0;
                twoDeqStream l3Stream = l3TempDeq.get_stream();

                while (!l3Stream.empty()) {
                    two_tuple tmp = *l3Stream;
                    offset_type& r_i = tmp.first;
                    offset_type& r_j = tmp.second;

                    if ((r_i > offset_type(0)) && (r_j > offset_type(0))) {
                        r_j = r_j - offset_type(1);
                        if (r_i >= offset_type(lcpn_ptr->size()) || r_j >= offset_type(lcpn_ptr->size()))
                            rmqDeq.push_back(intPair(id, 0, 0));
                        else if (r_i <= r_j)
                            rmqDeq.push_back(intPair(id, r_i, r_j));
                        else
                            rmqDeq.push_back(intPair(id, r_j, r_i));
                    }
                    else {
                        rmqDeq.push_back(intPair(id, 0, 0));
                    }
                    ++l3Stream;
                    ++id;
                }

                pairDeqStream rmqStream = rmqDeq.get_stream();
                simpleDeqStream lcpnStream = lcpn_ptr->get_stream();

                sparseTableAlgo<simpleDeqStream, pairDeqStream> sp_algo;
                sp_algo.initialize(lcpn_ptr->size()); // |lcpn| == |lcp12|
                sp_algo.splitAndSortRMQ(rmqStream);
                sp_algo.answerTuple(lcpnStream);
                sp_algo.answerRMQ(l3Deq);
            }
            else {
                simpleDeqStream l2Stream = l2Deq.get_stream();
                twoDeqStream l3Stream = l3TempDeq.get_stream();

                l3_tuple_type l3_tuple_sorter(l3Tuple_cmp(), ram_use, ram_use / 2);

                offset_type id = 0;

                while (!l3Stream.empty()) {
                    two_tuple tmp = *l3Stream;

                    offset_type& r_i = tmp.first;
                    offset_type& r_j = tmp.second;

                    if (r_i > offset_type(0) && r_j > offset_type(0)) {
                        l3_tuple_sorter.push(l3_type(id, 0, r_i - 1, *l2Stream)); // (id, k, r_i - 1, l2)
                        l3_tuple_sorter.push(l3_type(id, 1, r_j - 1, *l2Stream)); // (id, k, r_j - 1, l2)
                    }
                    else {
                        rmqDeq.push_back(intPair(id, 0, 0));                      // trivial RMQ
                    }
                    ++l2Stream;
                    ++l3Stream;
                    ++id;
                }

                l3_tuple_sorter.sort();
                inner_tuple_type inner_tuple_sorter(innerTuple_cmp(), ram_use / 2, ram_use / 2);

                /// Compute SA12[...] part.

                {
                    simpleDeqStream sa12Stream = sa12_ptr->get_stream();
                    offset_type flag = 0;

                    while (!l3_tuple_sorter.empty()) {
                        l3_type tmp = *l3_tuple_sorter; // tmp = (id, k, r, l)

                        while (flag < tmp.third) {
                            if (!sa12Stream.empty())
                                ++sa12Stream;

                            ++flag;
                        }

                        inner_tuple_sorter.push(inner_type(tmp.first, tmp.second, (*sa12Stream + tmp.fourth)));
                        ++l3_tuple_sorter;
                    }

                    l3_tuple_sorter.finish_clear();
                    inner_tuple_sorter.sort();
                }

                /// Compute ISA[...] part.

                simpleDeqStream isa1Stream = isa1_ptr->get_stream();
                simpleDeqStream isa2Stream = isa2_ptr->get_stream();

                s_concatenation isa12Stream(isa1Stream, isa2Stream);
                pair_type isa_sorter(pair_cmp(), ram_use / 2, ram_use / 2);

                offset_type add = 0;

                while (!inner_tuple_sorter.empty()) {
                    inner_type tmp = *inner_tuple_sorter; // tmp (id, k, r)

                    while (add < tmp.third) {
                        if (!isa12Stream.empty())
                            ++isa12Stream;

                        ++add;
                    }

                    if (tmp.second == 1) {
                        if (tmp.third <= offset_type(sa12_ptr->size() - 1))
                            isa_sorter.push(pos_pair_type(tmp.first, tmp.second, (*isa12Stream - 1)));
                        else
                            isa_sorter.push(pos_pair_type(tmp.first, tmp.second, 0));
                    }
                    else {
                        if (tmp.third <= offset_type(sa12_ptr->size() - 1))
                            isa_sorter.push(pos_pair_type(tmp.first, tmp.second, *isa12Stream));
                        else
                            isa_sorter.push(pos_pair_type(tmp.first, tmp.second, 0));
                    }

                    ++inner_tuple_sorter;
                }

                inner_tuple_sorter.finish_clear();
                isa_sorter.sort();

                delete sa12_ptr;

                while (!isa_sorter.empty()) {
                    pos_pair_type tmp1 = *isa_sorter; // id, r_i, r_j
                    ++isa_sorter;
                    pos_pair_type tmp2 = *isa_sorter;

                    if (tmp1.third <= tmp2.third)
                        rmqDeq.push_back(intPair(tmp1.first, tmp1.third, tmp2.third));
                    else
                        rmqDeq.push_back(intPair(tmp1.first, tmp2.third, tmp1.third));

                    ++isa_sorter;
                }

                isa_sorter.finish_clear();

                pairDeqStream rmqStream = rmqDeq.get_stream();
                simpleDeqStream lcpnStream = lcpn_ptr->get_stream();

                sparseTableAlgo<simpleDeqStream, pairDeqStream> sp_algo;
                sp_algo.initialize(lcpn_ptr->size());
                sp_algo.splitAndSortRMQ(rmqStream);
                sp_algo.answerTuple(lcpnStream);
                sp_algo.answerRMQ(l3Deq);
            }   // else

            delete lcpn_ptr;

            l3_ptr = new simpleDeqStream(l3Deq);

            if (!l2_ptr->empty())
                result = *(*l1_ptr) + (3 * *(*l2_ptr)) + *(*l3_ptr);
            else
                result = *(*l1_ptr) + *(*l3_ptr);
        }   //answerL3()

        const offset_type& operator * () const
        {
            return result;
        }

        build_lcp& operator ++ ()
        {
            assert(!empty());
            assert(!l1_ptr->empty());
            assert(!l3_ptr->empty());

            if (!l2_ptr->empty())
                ++(*l2_ptr);

            ++(*l1_ptr);
            ++(*l3_ptr);
            if (!l1_ptr->empty() && !l3_ptr->empty()) {
                if (!l2_ptr->empty())
                    result = *(*l1_ptr) + (3 * *(*l2_ptr)) + *(*l3_ptr);
                else
                    result = *(*l1_ptr) + *(*l3_ptr);
            }
            else {
                blank = true;

                assert(l1_ptr != NULL && l2_ptr != NULL && l3_ptr != NULL);
                delete l1_ptr;
                delete l2_ptr;
                delete l3_ptr;
            }

            return *this;
        }

        bool empty() const
        {
            return blank;
        }
    };  //build_lcp()

    /**
     * Create the suffix array and the lcp array from the current sub problem
     * by simple comparison-based merging.  More precisely: compare characters(out of
     * text t) and ranks(out of ISA12) of the following parameter constellation:
     *
     * \param Mod0 5-tuple (quint): <i, t_i, t_{i+1}, ISA12[i+1], ISA12[i+2]>
     * \param Mod1 4-tuple (quad): <i, ISA12[i], t_i, ISA12[i+1]>
     * \param Mod2 5-tuple (quint): <i, ISA[i], t_i, t_{i+1}, ISA12[i+1]>
     */
    template <class Mod0, class Mod1, class Mod2>
    class merge_sa_lcp
    {
    public:
        typedef offset_type value_type;

    private:
        Mod0& A;
        Mod1& B;
        Mod2& C;

        skew_quint_type s0;
        skew_quad_type s1;
        skew_quint_type s2;

        int selected;

        bool exists[3];

        offset_type index;
        offset_type merge_result;

        offset_type last1_winner;
        offset_type last2_winner;

        skew_quint_type last_s0;
        skew_quad_type last_s1;
        skew_quint_type last_s2;

        offset_type last_type;
        offset_type bound;

        build_lcp* b_lcp;

        bool cmp_mod12()
        {
            return s1.second < s2.second;
        }

        bool cmp_mod02()
        {
            if (s0.second == s2.third) {
                if (s0.third == s2.fourth)
                    return s0.fifth < s2.fifth;
                else
                    return s0.third < s2.fourth;
            }
            else {
                return s0.second < s2.third;
            }
        }

        bool cmp_mod01()
        {
            if (s0.second == s1.third)
                return s0.fourth < s1.fourth;
            else
                return s0.second < s1.third;
        }

        /// Comparison part of merge.

        char l12_construction0(offset_type c1, offset_type c2)
        {
            if (last_type == offset_type(0)) {
                if (last1_winner == c1) {
                    if (last2_winner == c2) {
                        b_lcp->saveRMQ(last_s0.fifth, s0.fifth);
                        b_lcp->preprocessL3(last_s0.fifth, s0.fifth);
                        return 2;
                    }
                    else {
                        b_lcp->saveRMQ(0, 0);
                        b_lcp->preprocessL3(0, 0);
                        return 1;
                    }
                }
                else {
                    b_lcp->saveRMQ(0, 0);
                    b_lcp->preprocessL3(0, 0);
                    return 0;
                }
            }
            else if (last_type == offset_type(1)) {
                if (last1_winner == c1) {
                    b_lcp->saveRMQ(last_s1.fourth, s0.fourth);
                    b_lcp->preprocessL3(last_s1.fourth, s0.fourth);
                    return 1;
                }
                else {
                    b_lcp->saveRMQ(0, 0);
                    b_lcp->preprocessL3(0, 0);
                    return 0;
                }
            }
            else {
                if (last1_winner == c1) {
                    if (last2_winner == c2) {
                        b_lcp->saveRMQ(last_s2.fifth, s0.fifth);
                        b_lcp->preprocessL3(last_s2.fifth, s0.fifth);
                        return 2;
                    }
                    else {
                        b_lcp->saveRMQ(0, 0);
                        b_lcp->preprocessL3(0, 0);
                        return 1;
                    }
                }
                else {
                    b_lcp->saveRMQ(0, 0);
                    b_lcp->preprocessL3(0, 0);
                    return 0;
                }
            }
        }

        bool l12_construction1(offset_type c1)
        {
            if (last_type == offset_type(0)) {
                if (last1_winner == c1) {
                    b_lcp->saveRMQ(last_s0.fourth, s1.fourth);
                    b_lcp->preprocessL3(last_s0.fourth, s1.fourth);
                    return 1;
                }
                else {
                    b_lcp->saveRMQ(0, 0);
                    b_lcp->preprocessL3(0, 0);
                    return 0;
                }
            }
            else {   // i.e. S1-S1 or S2-S1
                if (last_type == offset_type(1)) {
                    b_lcp->saveRMQ(s1.second - 1, s1.second);
                    b_lcp->preprocessL3(last_s1.second, s1.second);
                }
                else {
                    b_lcp->saveRMQ(s1.second - 1, s1.second);
                    b_lcp->preprocessL3(last_s2.second, s1.second);
                }
                return 0;
            }
        }

        char l12_construction2(offset_type c1, offset_type c2)
        {
            if (last_type == offset_type(0)) {
                if (last1_winner == c1) {
                    if (last2_winner == c2) {
                        b_lcp->saveRMQ(last_s0.fifth, s2.fifth);
                        b_lcp->preprocessL3(last_s0.fifth, s2.fifth);
                        return 2;
                    }
                    else {
                        b_lcp->saveRMQ(0, 0);
                        b_lcp->preprocessL3(0, 0);
                        return 1;
                    }
                }
                else {
                    b_lcp->saveRMQ(0, 0);
                    b_lcp->preprocessL3(0, 0);
                    return 0;
                }
            }
            else {   // i.e. S1-S2 or S2-S2
                if (last_type == offset_type(1)) {
                    b_lcp->saveRMQ(s2.second - 1, s2.second);
                    b_lcp->preprocessL3(last_s1.second, s2.second);
                }
                else {
                    b_lcp->saveRMQ(s2.second - 1, s2.second);
                    b_lcp->preprocessL3(last_s2.second, s2.second);
                }
                return 0;
            }
        }

        void get012()
        {
            if (cmp_mod01()) {
                if (cmp_mod02()) {
                    selected = 0;
                    merge_result = s0.first;

                    b_lcp->saveL1(l12_construction0(s0.second, s0.third));
                    last1_winner = s0.second;
                    last2_winner = s0.third;

                    last_type = 0;
                    last_s0 = s0;
                }
                else {
                    selected = 2;
                    merge_result = s2.first;

                    b_lcp->saveL1(l12_construction2(s2.third, s2.fourth));
                    last1_winner = s2.third;
                    last2_winner = s2.fourth;

                    last_type = 2;
                    last_s2 = s2;
                }
            }
            else {
                if (cmp_mod12()) {
                    selected = 1;
                    merge_result = s1.first;

                    b_lcp->saveL1(l12_construction1(s1.third));
                    last1_winner = s1.third;

                    last_type = 1;
                    last_s1 = s1;
                }
                else {
                    selected = 2;
                    merge_result = s2.first;

                    b_lcp->saveL1(l12_construction2(s2.third, s2.fourth));
                    last1_winner = s2.third;
                    last2_winner = s2.fourth;

                    last_type = 2;
                    last_s2 = s2;
                }
            }
        }

        void get01()
        {
            if (cmp_mod01()) {
                selected = 0;
                merge_result = s0.first;

                b_lcp->saveL1(l12_construction0(s0.second, s0.third));
                last1_winner = s0.second;
                last2_winner = s0.third;

                last_type = 0;
                last_s0 = s0;
            }
            else {
                selected = 1;
                merge_result = s1.first;

                b_lcp->saveL1(l12_construction1(s1.third));
                last1_winner = s1.third;

                last_type = 1;
                last_s1 = s1;
            }
        }

        void get12()
        {
            if (cmp_mod12()) {
                selected = 1;
                merge_result = s1.first;

                b_lcp->saveL1(l12_construction1(s1.third));
                last1_winner = s1.third;

                last_type = 1;
                last_s1 = s1;
            }
            else {
                selected = 2;
                merge_result = s2.first;

                b_lcp->saveL1(l12_construction2(s2.third, s2.fourth));
                last1_winner = s2.third;
                last2_winner = s2.fourth;

                last_type = 2;
                last_s2 = s2;
            }
        }

        void get02()
        {
            if (cmp_mod02()) {
                selected = 0;
                merge_result = s0.first;

                b_lcp->saveL1(l12_construction0(s0.second, s0.third));
                last1_winner = s0.second;
                last2_winner = s0.third;

                last_type = 0;
                last_s0 = s0;
            }
            else {
                selected = 2;
                merge_result = s2.first;

                b_lcp->saveL1(l12_construction2(s2.third, s2.fourth));
                last1_winner = s2.third;
                last2_winner = s2.fourth;

                last_type = 2;
                last_s2 = s2;
            }
        }

        void solve()
        {
            if (exists[0]) {
                if (exists[1]) {
                    if (exists[2]) {
                        get012();
                    }
                    else {
                        get01();
                    }
                }
                else {
                    if (exists[2]) {
                        get02();
                    }
                    else {
                        selected = 0;
                        merge_result = s0.first;

                        b_lcp->saveL1(l12_construction0(s0.second, s0.third));
                        last1_winner = s0.second;
                        last2_winner = s0.third;

                        last_type = 0;
                        last_s0 = s0;
                    }
                }
            }
            else {
                if (exists[1]) {
                    if (exists[2]) {
                        get12();
                    }
                    else {
                        selected = 1;
                        merge_result = s1.first;

                        b_lcp->saveL1(l12_construction1(s1.third));
                        last1_winner = s1.third;

                        last_type = 1;
                        last_s1 = s1;
                    }
                }
                else {
                    if (exists[2]) {
                        selected = 2;
                        merge_result = s2.first;

                        b_lcp->saveL1(l12_construction2(s2.third, s2.fourth));
                        last1_winner = s2.third;
                        last2_winner = s2.fourth;

                        last_type = 2;
                        last_s2 = s2;
                    }
                    else {
                        assert(false);
                    }
                }
            }
        }

    public:
        bool empty() const
        {
            return (A.empty() && B.empty() && C.empty());
        }

        merge_sa_lcp(Mod0& x1, Mod1& x2, Mod2& x3, build_lcp* lcp_ptr)
            : A(x1), B(x2), C(x3), selected(-1), index(0)
        {
            assert(!A.empty());
            assert(!B.empty());
            assert(!C.empty());
            exists[0] = true;
            exists[1] = true;
            exists[2] = true;
            s0 = *A;
            s1 = *B;
            s2 = *C;

            last1_winner = std::numeric_limits<offset_type>::max();
            last2_winner = std::numeric_limits<offset_type>::max();

            last_type = 3; // last_type \in [0,2]

            last_s0 = s0;
            last_s1 = s1;
            last_s2 = s2;

            b_lcp = lcp_ptr;

            solve();
        }

        const value_type& operator * () const
        {
            return merge_result;
        }

        merge_sa_lcp& operator ++ ()
        {
            if (selected == 0) {
                assert(!A.empty());
                ++A;
                if (!A.empty())
                    s0 = *A;
                else
                    exists[0] = false;
            }
            else if (selected == 1) {
                assert(!B.empty());
                ++B;
                if (!B.empty())
                    s1 = *B;
                else
                    exists[1] = false;
            }
            else {
                assert(!C.empty());
                assert(selected == 2);
                ++C;
                if (!C.empty())
                    s2 = *C;
                else
                    exists[2] = false;
            }

            ++index;
            if (!empty())
                solve();

            return *this;
        }
    };

    /**
     * Helper function for computing the size of the 2/3 subproblem.
     */
    static inline offset_type subp_size(size_type n)
    {
        return offset_type((n / 3) * 2 + ((n % 3) == 2));
    }

    /**
     * Sort mod0-quints / mod1-quads / mod2-quints and
     * run merge_sa_lcp class to merge them together.
     * \param S input string pipe type.
     * \param Mod1 mod1 tuples input pipe type.
     * \param Mod2 mod2 tuples input pipe type.
     */
    template <class S, class Mod1, class Mod2>
    class build_sa
    {
    public:
        typedef offset_type value_type;
        static const unsigned int add_rank = 1; // free first rank to mark ranks beyond end of input

    private:
        // Mod1 types
        typedef typename stxxl::stream::use_push<skew_quad_type> mod1_push_type;
        typedef typename stxxl::stream::runs_creator<mod1_push_type, less_quad_2nd, block_size, stxxl::RC> mod1_runs_type;
        typedef typename mod1_runs_type::sorted_runs_type sorted_mod1_runs_type;
        typedef typename stxxl::stream::runs_merger<sorted_mod1_runs_type, less_quad_2nd> mod1_rm_type;

        // Mod2 types
        typedef typename stxxl::stream::use_push<skew_quint_type> mod2_push_type;
        typedef typename stxxl::stream::runs_creator<mod2_push_type, less_quint_2nd, block_size, stxxl::RC> mod2_runs_type;
        typedef typename mod2_runs_type::sorted_runs_type sorted_mod2_runs_type;
        typedef typename stxxl::stream::runs_merger<sorted_mod2_runs_type, less_quint_2nd> mod2_rm_type;

        // Mod0 types
        typedef typename stxxl::stream::use_push<skew_quint_type> mod0_push_type;
        typedef typename stxxl::stream::runs_creator<mod0_push_type, less_mod0, block_size, stxxl::RC> mod0_runs_type;
        typedef typename mod0_runs_type::sorted_runs_type sorted_mod0_runs_type;
        typedef typename stxxl::stream::runs_merger<sorted_mod0_runs_type, less_mod0> mod0_rm_type;

        // Merge type
        typedef merge_sa_lcp<mod0_rm_type, mod1_rm_type, mod2_rm_type> merge_sa_lcp_type;

        // Functions
        less_mod0 c0;
        less_quad_2nd c1;
        less_quint_2nd c2;

        // Runs merger
        mod1_rm_type* mod1_result;
        mod2_rm_type* mod2_result;
        mod0_rm_type* mod0_result;

        // Merger
        merge_sa_lcp_type* vmerge_sa_lcp;

        // Input
        S& source;
        Mod1& mod_1;
        Mod2& mod_2;

        // Tmp variables
        offset_type t[3];
        offset_type old_t2;
        offset_type old_mod2;
        bool exists[3];
        offset_type mod_one;
        offset_type mod_two;

        offset_type index;

        // Empty_flag
        bool ready;

        // Result
        value_type result;

        // lcp part
        simpleDeq ISA1;
        simpleDeq ISA2;
        build_lcp* b_lcp;

    public:
        build_sa(S& source_, Mod1& mod_1_, Mod2& mod_2_, size_type a_size, std::size_t memsize,
                 simpleDeq* lcp_names, build_lcp* lcp12, simpleDeq* SA_12)
            : source(source_), mod_1(mod_1_), mod_2(mod_2_), index(0), ready(false)
        {
            assert(!source_.empty());

            // Runs storage

            // input: ISA_1,2 from previous level
            mod0_runs_type mod0_runs(c0, memsize / 4);
            mod1_runs_type mod1_runs(c1, memsize / 4);
            mod2_runs_type mod2_runs(c2, memsize / 4);

            while (!source.empty())
            {
                exists[0] = false;
                exists[1] = false;
                exists[2] = false;

                if (!source.empty()) {
                    t[0] = *source;
                    ++source;
                    exists[0] = true;
                }

                if (!source.empty()) {
                    assert(!mod_1.empty());
                    t[1] = *source;
                    ++source;
                    mod_one = *mod_1 + add_rank;
                    ISA1.push_back(mod_one);
                    ++mod_1; // isa1
                    exists[1] = true;
                }

                if (!source.empty()) {
                    assert(!mod_2.empty());
                    t[2] = *source;
                    ++source;
                    mod_two = *mod_2 + add_rank;
                    ISA2.push_back(mod_two);
                    ++mod_2; // isa2
                    exists[2] = true;
                }

                // Check special cases in the middle of "source"
                // Cases are cx|xc cxx|cxx and cxxc|xxc

                assert(t[0] != offset_type(0));
                assert(t[1] != offset_type(0));
                assert(t[2] != offset_type(0));

                // modX = isaX
                // Mod 0 : (index0,char0,char1,mod1,mod2)
                // Mod 1 : (index1,mod1,char1,mod2)
                // Mod 2 : (index2,mod2)

                if (exists[2]) { // Nothing is missed
                    mod0_runs.push(skew_quint_type(index, t[0], t[1], mod_one, mod_two));
                    mod1_runs.push(skew_quad_type(index + 1, mod_one, t[1], mod_two));

                    if (index != offset_type(0))
                        mod2_runs.push(skew_quint_type((index - 1), old_mod2, old_t2, t[0], mod_one));
                }
                else if (exists[1]) { // Last element missed
                    mod0_runs.push(skew_quint_type(index, t[0], t[1], mod_one, 0));
                    mod1_runs.push(skew_quad_type(index + 1, mod_one, t[1], 0));

                    if (index != offset_type(0))
                        mod2_runs.push(skew_quint_type((index - 1), old_mod2, old_t2, t[0], mod_one));
                }
                else {  // Only one element left
                    assert(exists[0]);
                    mod0_runs.push(skew_quint_type(index, t[0], 0, 0, 0));

                    if (index != offset_type(0))
                        mod2_runs.push(skew_quint_type((index - 1), old_mod2, old_t2, t[0], 0));
                }

                old_mod2 = mod_two;
                old_t2 = t[2];
                index += 3;
            }

            if ((a_size % 3) == 0) { // changed
                if (index != offset_type(0))
                    mod2_runs.push(skew_quint_type((index - 1), old_mod2, old_t2, 0, 0));
            }

            if ((a_size % 3) == 1) {
                mod_one = *mod_1 + add_rank;
                ISA1.push_back(mod_one);
            }

            mod0_runs.deallocate();
            mod1_runs.deallocate();
            mod2_runs.deallocate();

            // Prepare for merging
            mod0_result = new mod0_rm_type(mod0_runs.result(), less_mod0(), memsize / 5);
            mod1_result = new mod1_rm_type(mod1_runs.result(), less_quad_2nd(), memsize / 5);
            mod2_result = new mod2_rm_type(mod2_runs.result(), less_quint_2nd(), memsize / 5);

            // output: ISA_1,2 for next level
            b_lcp = new build_lcp(lcp_names, lcp12, &ISA1, &ISA2, SA_12);
            vmerge_sa_lcp = new merge_sa_lcp_type(*mod0_result, *mod1_result, *mod2_result, b_lcp); // S_0, S_1, S_2, b_lcp

            result = *(*vmerge_sa_lcp);
        }

        const value_type& operator * () const
        {
            return result;
        }

        build_sa& operator ++ ()
        {
            assert(vmerge_sa_lcp != 0 && !vmerge_sa_lcp->empty());

            ++(*vmerge_sa_lcp);
            if (!vmerge_sa_lcp->empty()) {
                result = *(*vmerge_sa_lcp);
            }
            else {
                assert(vmerge_sa_lcp->empty());
                ready = true;

                assert(vmerge_sa_lcp != NULL);
                delete vmerge_sa_lcp;
                vmerge_sa_lcp = NULL;

                assert(mod0_result != NULL && mod1_result != NULL && mod2_result != NULL);
                delete mod0_result;
                mod0_result = NULL;
                delete mod1_result;
                mod1_result = NULL;
                delete mod2_result;
                mod2_result = NULL;
            }

            return *this;
        }

        bool empty() const
        {
            return ready;
        }

        build_lcp * finalize_lcp()
        {
            assert(ready); // merge *must* be finished here
            b_lcp->finalize();

            return b_lcp;
        }
    };

    /// The core of DC3-LCP aka skew3-lcp algorithm.

    /**
     * The algorithm class calls all sub-parts of DC3-LCP from this section and collect their returns.
     *
     * \param Input type of the input pipe.
     * \param block_size block size of an external sorter we use.
     */
    template <class Input>
    class algorithm
    {
    public:
        typedef offset_type value_type;
        typedef typename Input::value_type alphabet_type;

    private:
        typedef counter<offset_type> counter_stream_type;
        // (t_i) -> (i,t_i)
        typedef make_pairs<counter_stream_type, Input> make_pairs_input_type;

        bool finished;          // used for empty() check
        bool unique;            // is the current quad array unique?
        unsigned int rec_depth; // current recusion depth of the algorithm

        typedef mod12Comparator mod12cmp;
        typedef stxxl::sorter<stxxl::tuple<offset_type, offset_type>, mod12cmp, block_size> mod12_sorter_type;

        typedef stxxl::stream::choose<mod12_sorter_type, 2> isa_second_type;
        typedef build_sa<offset_array_it_rg, isa_second_type, isa_second_type> buildSA_type;
        typedef make_pairs<buildSA_type, counter_stream_type> precompute_isa_type;

        buildSA_type* out_sa; // points at final constructed suffix array

    private:
        /// Recursion of dc3-lcp.

        template <typename InputType1>
        buildSA_type * dc3_lcp(InputType1& p_Input)
        {
            // (t_i) -> (i,t_i,t_{i+1},t_{i+2})
            typedef make_quads<InputType1, offset_type, 1> make_quads_input_type;

            // (t_i) -> (i,t_i,t_{i+1},t_{i+2}) with i = 1,2 mod 3
            typedef extract_mod12<make_quads_input_type> mod12_quads_input_type;

            // sort (i,t_i,t_{i+1},t_{i+2}) by (t_i,t_{i+1},t_{i+2})
            typedef typename stxxl::stream::sort<mod12_quads_input_type, less_quad<offset_type>, block_size> sort_mod12_input_type;

            // name (i,t_i,t_{i+1},t_{i+2}) -> (i,n_i)
            typedef naming<sort_mod12_input_type> naming_input_type;

            unique = false;
            offset_type sticking_length = 0; // holds length of current s^12

            mod12_sorter_type m1_sorter(mod12cmp(), ram_use / 6);
            mod12_sorter_type m2_sorter(mod12cmp(), ram_use / 6);

            // sorted mod1 runs -concat- sorted mod2 runs
            typedef concat<stxxl::tuple<offset_type, offset_type>, mod12_sorter_type, mod12_sorter_type> concatenation;

            offset_array_type text;

            // (t_i) -> (i,t_i,t_{i+1},t_{i+2})
            make_quads_input_type quads_input(p_Input, text);

            // (t_i) -> (i,t_i,t_{i+1},t_{i+2}) with i = 1,2 mod 3
            mod12_quads_input_type mod12_quads_input(quads_input);

            // sort (i,t_i,t_{i+1},t_{i+2}) by (t_i,t_i+1},t_{i+2})
            sort_mod12_input_type sort_mod12_input(mod12_quads_input, less_quad<offset_type>(), ram_use / 6);

            // lcpn sequence
            simpleDeq* lcp_names = new simpleDeq;

            // name (i,t_i,t_{i+1},t_{i+2}) -> (i,"n_i")
            naming_input_type names_input(sort_mod12_input, unique, lcp_names); // fill lcp_names while lexnaming step

            // create (i, s^12[i])
            while (!names_input.empty()) {
                const skew_pair_type& tmp = *names_input;
                if (tmp.first & 1)
                    m2_sorter.push(tmp);
                else
                    m1_sorter.push(tmp);

                ++names_input;
                sticking_length = sticking_length + 1;
                //sticking_length++; //old
            }

            m1_sorter.sort();
            m2_sorter.sort();

            if (!unique) {
                ++rec_depth;
                (std::cout << "lexicographical names are not unique -> recusion depth=" << rec_depth << std::endl).flush();

                // compute s^12 := lexname[S[1 mod 3]] '+' lexname[S[2 mod 3]], (also known as reduced recursion string 'R')
                concatenation stick_mod1mod2(m1_sorter, m2_sorter);

                buildSA_type* recType = dc3_lcp(stick_mod1mod2);

                // deallocate runs_creator and runs_creator
                m1_sorter.finish_clear();
                m2_sorter.finish_clear();

                --rec_depth;
                (std::cout << "recusion depth=" << rec_depth << std::endl).flush();

                counter_stream_type sa_loop_index;
                precompute_isa_type sa_pairs(*recType, sa_loop_index); // => (SA12, i)

                // store beginning of mod2-tuples of s^12 in mod2_pos
                offset_type special = (sticking_length != subp_size(text.size()));
                offset_type mod2_pos = (subp_size(text.size()) >> 1) + (subp_size(text.size()) & 1) + special;

                mod12_sorter_type isa1_type(mod12cmp(), ram_use / 6);
                mod12_sorter_type isa2_type(mod12cmp(), ram_use / 6);

                simpleDeq* SA12 = new simpleDeq;

                while (!sa_pairs.empty()) {
                    const skew_pair_type& tmp = *sa_pairs; // tmp.first is exactly the SA12 array
                    SA12->push_back(tmp.first);            // save current SA^12 vector
                    if (tmp.first < mod2_pos)
                        isa1_type.push(tmp);
                    else
                        isa2_type.push(tmp);

                    ++sa_pairs;
                }

                isa1_type.finish();
                isa2_type.finish();

                build_lcp* lcp12 = recType->finalize_lcp();

                delete recType;

                offset_array_it_rg input(text.begin(), text.end());

                // => (i, ISA)
                isa1_type.sort(ram_use / 8);
                isa2_type.sort(ram_use / 8);

                // pick ISA of (i, ISA)
                isa_second_type isa1_(isa1_type);
                isa_second_type isa2_(isa2_type);

                return new buildSA_type(input, isa1_, isa2_, text.size(), ram_use, lcp_names, lcp12, SA12);
            }
            // else (if unique)
            (std::cout << "unique lexicographical names! " << std::endl).flush();

            isa_second_type isa1(m1_sorter);
            isa_second_type isa2(m2_sorter);

            offset_array_it_rg source(text.begin(), text.end());

            // build lcp12[i] = 0 forall i -> use NULL
            return new buildSA_type(source, isa1, isa2, text.size(), ram_use, lcp_names, NULL, NULL);
        }   //dc3-lcp()

    public:
        algorithm(Input& data_in) : finished(false), rec_depth(0)
        {
            // (t_i) -> (i,t_i)
            counter_stream_type dummy;
            make_pairs_input_type pairs_input(dummy, data_in);

            this->out_sa = dc3_lcp(pairs_input);
        }

        ~algorithm()
        {
            delete out_sa;
            out_sa = NULL;
        }

        const value_type& operator * () const
        {
            return *(*out_sa);
        }

        algorithm& operator ++ ()
        {
            assert(!out_sa->empty());

            ++(*out_sa);

            if ((out_sa != NULL) && (out_sa->empty()))
                finished = true;

            return *this;
        }

        bool empty() const
        {
            return finished;
        }

        build_lcp * finish_lcp()
        {
            return out_sa->finalize_lcp();
        }
    }; // algorithm class
};     // skew class

} // namespace algo

// Helper to print out readable characters.
template <typename alphabet_type>
static inline std::string dumpC(alphabet_type c)
{
    std::ostringstream oss;
    if (isalnum(c)) oss << '\'' << (char)c << '\'';
    else oss << (int)c;
    return oss.str();
}

/**
 * Helper class which cuts input off an input stream at a specified length.
 */
template <typename InputType>
class cut_stream
{
public:
    // same value type as input stream
    typedef typename InputType::value_type value_type;

protected:
    // instance of input stream
    InputType& m_input;

    // counter after which the stream ends
    size_type m_count;

public:
    cut_stream(InputType& input, size_type count)
        : m_input(input), m_count(count) { }

    const value_type& operator * () const
    {
        assert(m_count > 0);
        return *m_input;
    }

    cut_stream& operator ++ ()
    {
        assert(!empty());
        --m_count;
        ++m_input;
        return *this;
    }

    bool empty() const
    {
        return (m_count == 0) || m_input.empty();
    }
};

alphabet_type unary_generator()
{
    return 'a';
}

template <typename offset_type>
int process(const std::string& input_filename, const std::string& output_filename,
            size_type sizelimit, bool text_output_flag, bool check_flag, bool input_verbatim)
{
    static const std::size_t block_size = sizeof(offset_type) * 1024 * 1024 / 2;

    typedef typename stxxl::VECTOR_GENERATOR<alphabet_type, 1, 2>::result alphabet_vector_type;
    typedef typename stxxl::VECTOR_GENERATOR<offset_type, 1, 2, block_size>::result offset_vector_type;

    // input and output files (if supplied via command line)
    stxxl::syscall_file* input_file = NULL, * output_file = NULL;

    // input and output vectors for suffix array construction
    alphabet_vector_type input_vector;

    offset_vector_type output_vector;
    offset_vector_type lcparray;

    // to verify lcparray with kasai's semi external algorithm
    offset_vector_type checker_lcp;

    using stxxl::file;

    if (input_verbatim) {
        // copy input verbatim into vector
        input_vector.resize(input_filename.size());
        std::copy(input_filename.begin(), input_filename.end(), input_vector.begin());
    }
    else if (input_filename == "random") {
        if (sizelimit == std::numeric_limits<size_type>::max()) {
            std::cout << "You must provide -s <size> for generated inputs." << std::endl;
            return 1;
        }

        // fill input vector with random bytes
        input_vector.resize(sizelimit);
        stxxl::random_number8_r rand8;
        stxxl::generate(input_vector.begin(), input_vector.end(), rand8);
    }
    else if (input_filename == "unary") {
        if (sizelimit == std::numeric_limits<size_type>::max()) {
            std::cout << "You must provide -s <size> for generated inputs." << std::endl;
            return 1;
        }

        // fill input vector with a's
        input_vector.resize(sizelimit);
        stxxl::generate(input_vector.begin(), input_vector.end(), unary_generator);
    }
    else {
        // define input file object and map input_vector to input_file (no copying)
        input_file = new stxxl::syscall_file(input_filename, file::RDONLY | file::DIRECT);
        alphabet_vector_type file_input_vector(input_file);
        input_vector.swap(file_input_vector);
    }

    if (output_filename.size()) {
        // define output file object and map output_vector to output_file
        output_file = new stxxl::syscall_file(output_filename, file::RDWR | file::CREAT | file::DIRECT);
        offset_vector_type file_output_vector(output_file);
        output_vector.swap(file_output_vector);
    }

    // initialize and start I/O measurement
    stxxl::stats* Stats = stxxl::stats::get_instance();
    stxxl::stats_data stats_begin(*Stats);

    // construct skew class with bufreader input type
    typedef alphabet_vector_type::bufreader_type input_type;
    typedef cut_stream<input_type> cut_input_type;
    typedef typename algo::skew<offset_type>::template algorithm<cut_input_type> skew_type;

    size_type size = input_vector.size();
    if (size > sizelimit) size = sizelimit;

    std::cout << "input size = " << size << std::endl;

    if (size + 3 >= std::numeric_limits<offset_type>::max()) {
        std::cout << "error: input is too long for selected word size!" << std::endl;
        return -1;
    }

    input_type input(input_vector);
    cut_input_type cut_input(input, size);
    skew_type skew(cut_input);

    // make sure output vector has the right size
    output_vector.resize(size);
    lcparray.resize(size);
    checker_lcp.resize(size);

    // write suffix array stream into output vector
    stxxl::stream::materialize(skew, output_vector.begin(), output_vector.end());

    typename algo::skew<offset_type>::build_lcp * lcp = skew.finish_lcp();
    stxxl::stream::materialize(*lcp, lcparray.begin(), lcparray.end());

    std::cout << "output size = " << output_vector.size() << std::endl;
    std::cout << (stxxl::stats_data(*Stats) - stats_begin); // print i/o statistics

    if (text_output_flag) {
        std::cout << std::endl << "resulting suffix array:" << std::endl;

        for (size_type i = 0; i < output_vector.size(); i++) {
            std::cout << i << " : " << output_vector[i] << " : ";

            // We need a const reference because operator[] would write data
            const alphabet_vector_type& cinput = input_vector;

            for (size_type j = 0; output_vector[i] + j < cinput.size(); j++)
                std::cout << dumpC(cinput[output_vector[i] + j]) << " ";

            std::cout << std::endl;
        }

        std::cout << "resulting lcp array: \n";

        for (size_type k = 0; k < lcparray.size(); ++k)
            std::cout << k << " : " << dumpC(lcparray[k]) << std::endl;
    }

    int ret = 0;

    if (check_flag) {
        (std::cout << "checking suffix array... ").flush();

        if (!sacheck_vectors(input_vector, output_vector)) {
            std::cout << "failed!" << std::endl;
            ret = -1;
        }
        else
            std::cout << "ok." << std::endl;

        (std::cout << "checking lcp array... ").flush();

        lcparray_stxxl_kasai(input_vector, output_vector, checker_lcp);
        bool correct_lcp = true;

        for (unsigned int k = 0; k < lcparray.size(); ++k) {
            if (checker_lcp[k] != lcparray[k]) {
                std::cout << "failed!" << std::endl;
                std::cout << k << ", " << checker_lcp[k] << " - " << lcparray[k] << std::endl;
                correct_lcp = false;
                ret = -2;
                break;
            }
        }

        if (correct_lcp)
            std::cout << "ok." << std::endl;
    }

    // close file, but have to deallocate vector first!

    if (input_file) {
        input_vector = alphabet_vector_type();
        delete input_file;
    }
    if (output_file) {
        output_vector = offset_vector_type();
        delete output_file;
    }

    return ret;
}

int main(int argc, char* argv[])
{
    stxxl::cmdline_parser cp;

    cp.set_description(
        "DC3-LCP aka skew3-lcp algorithm for external memory suffix array and LCP array construction.");
    cp.set_author(
        "Jens Mehnert <jmehnert@mpi-sb.mpg.de>, \n"
        "Timo Bingmann <tb@panthema.net>, \n"
        "Daniel Feist <daniel.feist@student.kit.edu>");

    std::string input_filename, output_filename;
    size_type sizelimit = std::numeric_limits<size_type>::max();
    bool text_output_flag = false;
    bool check_flag = false;
    bool input_verbatim = false;
    unsigned wordsize = 32;

    cp.add_param_string("input", input_filename,
                        "Path to input file (or verbatim text).\n"
                        "  The special inputs 'random' and 'unary' generate "
                        "such text on-the-fly.");
    cp.add_flag('c', "check", check_flag,
                "Check suffix array for correctness.");
    cp.add_flag('t', "text", text_output_flag,
                "Print out suffix array in readable text.");
    cp.add_string('o', "output", output_filename,
                  "Output suffix array to given path.");
    cp.add_flag('v', "verbatim", input_verbatim,
                "Consider \"input\" as verbatim text to construct "
                "suffix array on.");
    cp.add_bytes('s', "size", sizelimit,
                 "Cut input text to given size,"
                 "units are [b,kb,mb,gb], e.g. 100kb.");
    cp.add_bytes('M', "memuse", ram_use,
                 "Amount of RAM to use, default: 1 GiB.");
    cp.add_uint('w', "wordsize", wordsize,
                "Set word size of suffix array to 32, 40 or 64 bit,"
                "default: 32-bit.");

    // process command line
    if (!cp.process(argc, argv))
        return -1;

    if (wordsize == 32)
        return process<stxxl::uint32>(
            input_filename, output_filename, sizelimit,
            text_output_flag, check_flag, input_verbatim);
    else if (wordsize == 40)
        return process<stxxl::uint40>(
            input_filename, output_filename, sizelimit,
            text_output_flag, check_flag, input_verbatim);
    else if (wordsize == 64)
        return process<stxxl::uint64>(
            input_filename, output_filename, sizelimit,
            text_output_flag, check_flag, input_verbatim);
    else
        std::cerr << "Invalid wordsize for suffix array: 32, 40 or 64 are allowed." << std::endl;

    return -1;
}
