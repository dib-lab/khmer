/***************************************************************************
 *  examples/applications/skew3.cpp
 *
 *  Implementation of the external memory suffix sorting algorithm DC3 aka
 *  skew3 as described in Roman Dementiev, Juha Kaerkkaeinen, Jens Mehnert and
 *  Peter Sanders. "Better External Memory Suffix Array Construction". Journal
 *  of Experimental Algorithmics (JEA), volume 12, 2008.
 *
 *  Part of the STXXL. See http://stxxl.sourceforge.net
 *
 *  Copyright (C) 2004 Jens Mehnert <jmehnert@mpi-sb.mpg.de>
 *  Copyright (C) 2012-2013 Timo Bingmann <tb@panthema.net>
 *  Copyright (C) 2012-2013 Daniel Feist <daniel.feist@student.kit.edu>
 *
 *  Distributed under the Boost Software License, Version 1.0.
 *  (See accompanying file LICENSE_1_0.txt or copy at
 *  http://www.boost.org/LICENSE_1_0.txt)
 **************************************************************************/

#include <algorithm>
#include <cassert>
#include <cctype>
#include <cstddef>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

#include <stxxl/algorithm>
#include <stxxl/cmdline>
#include <stxxl/io>
#include <stxxl/random>
#include <stxxl/sorter>
#include <stxxl/stats>
#include <stxxl/stream>
#include <stxxl/vector>
#include <stxxl/bits/common/uint_types.h>

using stxxl::uint64;
using stxxl::internal_size_type;
using stxxl::external_size_type;
namespace stream = stxxl::stream;

// 1 GiB ram used by external data structures / 1 MiB block size
internal_size_type ram_use = 1024 * 1024 * 1024;

// alphabet data type
typedef unsigned char alphabet_type;

// calculation data type
typedef external_size_type size_type;

/// Suffix Array checker for correctness verification

/**
 * Algorithm to check whether the suffix array is correct. Loosely based on the
 * ideas of Kaerkkaeinen und Burghardt, originally implemented in STXXL by Jens
 * Mehnert (2004), reimplemented using triples by Timo Bingmann (2012).
 *
 * \param InputT is the original text, from which the suffix array was build
 * \param InputSA is the suffix array from InputT
 *
 * Note: ISA := The inverse of SA
 */
template <typename InputT, typename InputSA>
bool sacheck(InputT& inputT, InputSA& inputSA)
{
    typedef typename InputSA::value_type offset_type;
    typedef stxxl::tuple<offset_type, offset_type> pair_type;
    typedef stxxl::tuple<offset_type, offset_type, offset_type> triple_type;

    // *** Pipeline Declaration ***

    // Build tuples with index: (SA[i]) -> (i, SA[i])
    typedef stxxl::stream::counter<offset_type> index_counter_type;
    index_counter_type index_counter;

    typedef stream::make_tuple<index_counter_type, InputSA> tuple_index_sa_type;
    tuple_index_sa_type tuple_index_sa(index_counter, inputSA);

    // take (i, SA[i]) and sort to (ISA[i], i)
    typedef stxxl::tuple_less2nd<pair_type> pair_less_type;
    typedef typename stream::sort<tuple_index_sa_type, pair_less_type> build_isa_type;

    build_isa_type build_isa(tuple_index_sa, pair_less_type(), ram_use / 3);

    // build (ISA[i], T[i], ISA[i+1]) and sort to (i, T[SA[i]], ISA[SA[i]+1])
    typedef stxxl::tuple_less1st<triple_type> triple_less_type;      // comparison relation

    typedef typename stream::use_push<triple_type> triple_push_type; // indicator use push()
    typedef typename stream::runs_creator<triple_push_type, triple_less_type> triple_rc_type;
    typedef typename stream::runs_merger<typename triple_rc_type::sorted_runs_type, triple_less_type> triple_rm_type;

    triple_rc_type triple_rc(triple_less_type(), ram_use / 3);

    // ************************* Process ******************************
    // loop 1: read ISA and check for a permutation. Simultaneously create runs
    // of triples by iterating ISA and T.

    size_type totalSize;
    {
        offset_type prev_isa = (*build_isa).first;
        offset_type counter = 0;
        while (!build_isa.empty())
        {
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

    // ************************************************************************
    // loop 2: read triples (i,T[SA[i]],ISA[SA[i]+1]) and check for correct
    // ordering.

    triple_rm_type triple_rm(triple_rc.result(), triple_less_type(), ram_use / 3);

    {
        triple_type prev_triple = *triple_rm;
        size_type counter = 0;

        ++triple_rm;

        while (!triple_rm.empty())
        {
            const triple_type& this_triple = *triple_rm;

            if (prev_triple.second > this_triple.second)
            {
                // simple check of first character of suffix
                std::cout << "Error: suffix array position " << counter << " ordered incorrectly." << std::endl;
                return false;
            }
            else if (prev_triple.second == this_triple.second)
            {
                if (this_triple.third == (offset_type)totalSize) {
                    // last suffix of string must be first among those with same
                    // first character
                    std::cout << "Error: suffix array position " << counter << " ordered incorrectly." << std::endl;
                    return false;
                }
                if (prev_triple.third != (offset_type)totalSize && prev_triple.third > this_triple.third) {
                    // positions SA[i] and SA[i-1] has same first character but
                    // their suffixes are ordered incorrectly: the suffix
                    // position of SA[i] is given by ISA[SA[i]]
                    std::cout << "Error: suffix array position " << counter << " ordered incorrectly." << std::endl;
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

template <typename InputT, typename InputSA>
bool sacheck_vectors(InputT& inputT, InputSA& inputSA)
{
    typename stream::streamify_traits<typename InputT::iterator>::stream_type streamT
        = stream::streamify(inputT.begin(), inputT.end());

    typename stream::streamify_traits<typename InputSA::iterator>::stream_type streamSA
        = stream::streamify(inputSA.begin(), inputSA.end());

    return sacheck(streamT, streamSA);
}

/// DC3 aka skew algorithm

/*
 * DC3 aka skew algorithm a short description. T := input string
 * The recursion works as follows:
 * Step 1: a) pick all mod1/mod2 triples (i.e. triples T[i,i+2] at position i mod 3 != 0) (-> extract_mod12 class)
 *         b) sort mod1/mod2 triples lexicographically (-> build_sa class)
 *         c) give mod1/mod2 triples lexicographical ascending names n (-> naming class)
 *         d) check lexicographical names for uniqueness (-> naming class)
 *            If yes: proceed to next Step, If no: set T := lexicographical names and run Step 1 again
 * Step 2: a) by sorting the lexicographical names n we receive ranks r
 *         b) construct mod0-quints, mod1-quads and mod2-quints  (-> build_sa class)
 *         c) prepare for merging by:
 *            sort mod0-quints by 2 components, sort mod1-quads / mod2-quints by one component (-> build_sa class)
 *         c) merge mod0-quints, mod1-quads and mod2-quints (-> merge_sa class)
 * Step 3: a) return Suffix Array of T
 *
 * \param offset_type later suffix array data type
 */
template <typename offset_type>
class skew
{
public:
    // 2-tuple, 3-tuple, 4-tuple (=quads), 5-tuple(=quints) definition
    typedef stxxl::tuple<offset_type, offset_type> skew_pair_type;
    typedef stxxl::tuple<offset_type, offset_type, offset_type> skew_triple_type;
    typedef stxxl::tuple<offset_type, offset_type, offset_type, offset_type> skew_quad_type;
    typedef stxxl::tuple<offset_type, offset_type, offset_type, offset_type, offset_type> skew_quint_type;

    typedef typename stxxl::VECTOR_GENERATOR<offset_type, 1, 2>::result offset_array_type;
    typedef stream::vector_iterator2stream<typename offset_array_type::iterator> offset_array_it_rg;

    /** Comparison function for the mod0 tuples. */
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

    typedef stxxl::tuple_less2nd<skew_quad_type> less_mod1;
    typedef stxxl::tuple_less2nd<skew_quint_type> less_mod2;

    /** Put the (0 mod 2) [which are the 1,2 mod 3 tuples] tuples at the begin. */
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

    /** Sort skew_quad datatype. */
    template <typename alphabet_type>
    struct less_quad
    {
        typedef stxxl::tuple<offset_type, alphabet_type, alphabet_type, alphabet_type> value_type;

        bool operator () (const value_type& a, const value_type& b) const
        {
            if (a.second == b.second) {
                if (a.third == b.third)
                    return a.fourth < b.fourth;
                else
                    return a.third < b.third;
            }
            else
                return a.second < b.second;
        }

        static value_type min_value() { return value_type::min_value(); }
        static value_type max_value() { return value_type::max_value(); }
    };

    /** Check, if last two components of tree quads are equal. */
    template <class quad_type>
    static inline bool quad_eq(const quad_type& a, const quad_type& b)
    {
        return (a.second == b.second) && (a.third == b.third) && (a.fourth == b.fourth);
    }

    /** Naming pipe for the conventional skew algorithm without discarding. */
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

    public:
        naming(Input& A_, bool& unique_)
            : A(A_), unique(unique_), lexname(0)
        {
            assert(!A.empty());
            unique = true;

            prev = *A;
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
            if (A.empty())
                return *this;

            quad_type curr = *A;
            if (!quad_eq(prev, curr)) {
                ++lexname;
            }
            else {
                if (!A.empty() && curr.second != offset_type(0)) {
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

    /** Create tuples of 2 components until one of the input streams are empty. */
    template <class InputA, class InputB, const int add_alphabet = 0>
    class make_pairs
    {
    public:
        typedef stxxl::tuple<typename InputA::value_type, offset_type> value_type;

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
            if (!empty()) {
                result = value_type(*A, *B + add_alphabet);
            }
        }

        const value_type& operator * () const
        { return result; }

        make_pairs& operator ++ ()
        {
            assert(!A.empty());
            assert(!B.empty());

            ++A;
            ++B;

            if (!A.empty() && !B.empty()) {
                result = value_type(*A, *B + add_alphabet);
            }

            return *this;
        }

        bool empty() const
        { return (A.empty() || B.empty()); }
    };

    /**
     * Collect three characters t_i, t_{i+1}, t_{i+2} beginning at the index
     * i. Since we need at least one unique endcaracter, we free the first
     * characters i.e. we map (t_i) -> (i,t_i,t_{i+1},t_{i+2})
     *
     * \param Input holds all characters t_i from input string t
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
        unsigned int z3z;  // = counter mod 3, ("+",Z/3Z) is cheaper than %
        bool finished;

        offset_array_type& text;

    public:
        make_quads(Input& data_in_, offset_array_type& text_)
            : A(data_in_),
              current(0, 0, 0, 0),
              counter(0),
              z3z(0),
              finished(false),
              text(text_)
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

            if (!A.empty()) {
                current.fourth = (*A).second + add_alphabet;
            }
            else {
                current.fourth = 0;
            }
        }

        const value_type& operator * () const
        { return current; }

        make_quads& operator ++ ()
        {
            assert(!A.empty() || !finished);

            if (current.second != offset_type(0)) {
                text.push_back(current.second);
            }

            // Calculate module
            if (++z3z == 3) z3z = 0;

            current.first = ++counter;
            current.second = current.third;
            current.third = current.fourth;

            if (!A.empty())
                ++A;

            if (!A.empty()) {
                current.fourth = (*A).second + add_alphabet;
            }
            else {
                current.fourth = 0;
            }

            // Inserts a dummy tuple for input sizes of n%3==1
            if ((current.second == offset_type(0)) && (z3z != 1)) {
                finished = true;
            }

            return *this;
        }

        bool empty() const
        { return (A.empty() && finished); }
    };

    /** Drop 1/3 of the input. More exactly the offsets at positions (0 mod
     * 3). Index begins with 0. */
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
        extract_mod12(Input& A_)
            : A(A_),
              counter(0),
              output_counter(0)
        {
            assert(!A.empty());
            ++A, ++counter;  // skip 0 = mod0 offset
            if (!A.empty()) {
                result = *A;
                result.first = output_counter;
            }
        }

        const value_type& operator * () const
        { return result; }

        extract_mod12& operator ++ ()
        {
            assert(!A.empty());

            ++A, ++counter, ++output_counter;

            if (!A.empty() && (counter % 3) == 0) {
                // skip mod0 offsets
                ++A, ++counter;
            }
            if (!A.empty()) {
                result = *A;
                result.first = output_counter;
            }

            return *this;
        }

        bool empty() const
        { return A.empty(); }
    };

    /** Create the suffix array from the current sub problem by simple
     *  comparison-based merging.  More precisely: compare characters(out of
     *  text t) and ranks(out of ISA12) of the following constellation:
     *  Input constellation:
     *  \param Mod0 5-tuple (quint): <i, t_i, t_{i+1}, ISA12[i+1], ISA12[i+2]>
     *  \param Mod1 4-tuple (quad): <i, ISA12[i], t_i, ISA12[i+1]>
     *  \param Mod2 5-tuple (quint): <i, ISA[i], t_i, t_{i+1}, ISA12[i+1]>
     */
    template <class Mod0, class Mod1, class Mod2>
    class merge_sa
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
        bool done[3];

        offset_type index;
        offset_type merge_result;

        bool cmp_mod1_less_mod2()
        {
            assert(!done[1] && !done[2]);

            return s1.second < s2.second;
        }

        bool cmp_mod0_less_mod2()
        {
            assert(!done[0] && !done[2]);

            if (s0.second == s2.third) {
                if (s0.third == s2.fourth)
                    return s0.fifth < s2.fifth;
                else
                    return s0.third < s2.fourth;
            }
            else
                return s0.second < s2.third;
        }

        bool cmp_mod0_less_mod1()
        {
            assert(!done[0] && !done[1]);

            if (s0.second == s1.third)
                return s0.fourth < s1.fourth;
            else
                return s0.second < s1.third;
        }

        void merge()
        {
            assert(!done[0] || !done[1] || !done[2]);

            if (done[0])
            {
                if (done[2] || (!done[1] && cmp_mod1_less_mod2()))
                {
                    selected = 1;
                    merge_result = s1.first;
                }
                else
                {
                    selected = 2;
                    merge_result = s2.first;
                }
            }
            else if (done[1] || cmp_mod0_less_mod1())
            {
                if (done[2] || cmp_mod0_less_mod2())
                {
                    selected = 0;
                    merge_result = s0.first;
                }
                else
                {
                    selected = 2;
                    merge_result = s2.first;
                }
            }
            else
            {
                if (done[2] || cmp_mod1_less_mod2())
                {
                    selected = 1;
                    merge_result = s1.first;
                }
                else
                {
                    selected = 2;
                    merge_result = s2.first;
                }
            }

            assert(!done[selected]);
        }

    public:
        bool empty() const
        {
            return (A.empty() && B.empty() && C.empty());
        }

        merge_sa(Mod0& x1, Mod1& x2, Mod2& x3)
            : A(x1), B(x2), C(x3), selected(-1), index(0)
        {
            assert(!A.empty());
            assert(!B.empty());
            assert(!C.empty());
            done[0] = false;
            done[1] = false;
            done[2] = false;
            s0 = *A;
            s1 = *B;
            s2 = *C;

            merge();
        }

        const value_type& operator * () const
        {
            return merge_result;
        }

        merge_sa& operator ++ ()
        {
            if (selected == 0) {
                assert(!A.empty());
                ++A;
                if (!A.empty())
                    s0 = *A;
                else
                    done[0] = true;
            }
            else if (selected == 1) {
                assert(!B.empty());
                ++B;
                if (!B.empty())
                    s1 = *B;
                else
                    done[1] = true;
            }
            else {
                assert(!C.empty());
                assert(selected == 2);
                ++C;
                if (!C.empty())
                    s2 = *C;
                else
                    done[2] = true;
            }

            ++index;
            if (!empty())
                merge();

            return *this;
        }
    };

    /** Helper function for computing the size of the 2/3 subproblem. */
    static inline size_type subp_size(size_type n)
    {
        return (n / 3) * 2 + ((n % 3) == 2);
    }

    /**
     * Sort mod0-quints / mod1-quads / mod2-quints and run merge_sa class to
     * merge them together.
     * \param S input string pipe type.
     * \param Mod1 mod1 tuples input pipe type.
     * \param Mod2 mod2 tuples input pipe type.
     */
    template <class S, class Mod1, class Mod2>
    class build_sa
    {
    public:
        typedef offset_type value_type;

        static const unsigned int add_rank = 1;  // free first rank to mark ranks beyond end of input

    private:
        // mod1 types
        typedef typename stream::use_push<skew_quad_type> mod1_push_type;
        typedef typename stream::runs_creator<mod1_push_type, less_mod1> mod1_runs_type;
        typedef typename mod1_runs_type::sorted_runs_type sorted_mod1_runs_type;
        typedef typename stream::runs_merger<sorted_mod1_runs_type, less_mod1> mod1_rm_type;

        // mod2 types
        typedef typename stream::use_push<skew_quint_type> mod2_push_type;
        typedef typename stream::runs_creator<mod2_push_type, less_mod2> mod2_runs_type;
        typedef typename mod2_runs_type::sorted_runs_type sorted_mod2_runs_type;
        typedef typename stream::runs_merger<sorted_mod2_runs_type, less_mod2> mod2_rm_type;

        // mod0 types
        typedef typename stream::use_push<skew_quint_type> mod0_push_type;
        typedef typename stream::runs_creator<mod0_push_type, less_mod0> mod0_runs_type;
        typedef typename mod0_runs_type::sorted_runs_type sorted_mod0_runs_type;
        typedef typename stream::runs_merger<sorted_mod0_runs_type, less_mod0> mod0_rm_type;

        // Merge type
        typedef merge_sa<mod0_rm_type, mod1_rm_type, mod2_rm_type> merge_sa_type;

        // Functions
        less_mod0 c0;
        less_mod1 c1;
        less_mod2 c2;

        // Runs merger
        mod1_rm_type* mod1_result;
        mod2_rm_type* mod2_result;
        mod0_rm_type* mod0_result;

        // Merger
        merge_sa_type* vmerge_sa;

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

    public:
        build_sa(S& source_, Mod1& mod_1_, Mod2& mod_2_, size_type a_size, size_t memsize)
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
                    ++mod_1;
                    exists[1] = true;
                }

                if (!source.empty()) {
                    assert(!mod_2.empty());
                    t[2] = *source;
                    ++source;
                    mod_two = *mod_2 + add_rank;
                    ++mod_2;
                    exists[2] = true;
                }

                // Check special cases in the middle of "source"
                // Cases are cx|xc cxx|cxx and cxxc|xxc

                assert(t[0] != offset_type(0));
                assert(t[1] != offset_type(0));
                assert(t[2] != offset_type(0));

                // Mod 0 : (index0,char0,char1,mod1,mod2)
                // Mod 1 : (index1,mod1,char1,mod2)
                // Mod 2 : (index2,mod2)

                if (exists[2]) { // Nothing is missed
                    mod0_runs.push(skew_quint_type(index, t[0], t[1], mod_one, mod_two));
                    mod1_runs.push(skew_quad_type(index + 1, mod_one, t[1], mod_two));

                    if (index != offset_type(0)) {
                        mod2_runs.push(skew_quint_type((index - 1), old_mod2, old_t2, t[0], mod_one));
                    }
                }
                else if (exists[1]) { // Last element missed
                    mod0_runs.push(skew_quint_type(index, t[0], t[1], mod_one, 0));
                    mod1_runs.push(skew_quad_type(index + 1, mod_one, t[1], 0));

                    if (index != offset_type(0)) {
                        mod2_runs.push(skew_quint_type((index - 1), old_mod2, old_t2, t[0], mod_one));
                    }
                }
                else { // Only one element left
                    assert(exists[0]);
                    mod0_runs.push(skew_quint_type(index, t[0], 0, 0, 0));

                    if (index != offset_type(0)) {
                        mod2_runs.push(skew_quint_type((index - 1), old_mod2, old_t2, t[0], 0));
                    }
                }

                old_mod2 = mod_two;
                old_t2 = t[2];
                index += 3;
            }

            if ((a_size % 3) == 0) { // changed
                if (index != offset_type(0)) {
                    mod2_runs.push(skew_quint_type((index - 1), old_mod2, old_t2, 0, 0));
                }
            }

            mod0_runs.deallocate();
            mod1_runs.deallocate();
            mod2_runs.deallocate();

            std::cout << "merging S0 = " << mod0_runs.size() << ", S1 = " << mod1_runs.size()
                      << ", S2 = " << mod2_runs.size() << " tuples" << std::endl;

            // Prepare for merging

            mod0_result = new mod0_rm_type(mod0_runs.result(), less_mod0(), memsize / 5);
            mod1_result = new mod1_rm_type(mod1_runs.result(), less_mod1(), memsize / 5);
            mod2_result = new mod2_rm_type(mod2_runs.result(), less_mod2(), memsize / 5);

            // output: ISA_1,2 for next level
            vmerge_sa = new merge_sa_type(*mod0_result, *mod1_result, *mod2_result);

            // read first suffix
            result = *(*vmerge_sa);
        }

        const value_type& operator * () const
        {
            return result;
        }

        build_sa& operator ++ ()
        {
            assert(vmerge_sa != 0 && !vmerge_sa->empty());

            ++(*vmerge_sa);
            if (!vmerge_sa->empty()) {
                result = *(*vmerge_sa);
            }
            else {  // cleaning up
                assert(vmerge_sa->empty());
                ready = true;

                assert(vmerge_sa != NULL);
                delete vmerge_sa, vmerge_sa = NULL;

                assert(mod0_result != NULL && mod1_result != NULL && mod2_result != NULL);
                delete mod0_result, mod0_result = NULL;
                delete mod1_result, mod1_result = NULL;
                delete mod2_result, mod2_result = NULL;
            }

            return *this;
        }

        ~build_sa()
        {
            if (vmerge_sa) delete vmerge_sa;

            if (mod0_result) delete mod0_result;
            if (mod1_result) delete mod1_result;
            if (mod2_result) delete mod2_result;
        }

        bool empty() const
        {
            return ready;
        }
    };

    /** The skew algorithm.
     *  \param Input type of the input pipe. */
    template <class Input>
    class algorithm
    {
    public:
        typedef offset_type value_type;
        typedef typename Input::value_type alphabet_type;

    protected:
        // finished reading final suffix array
        bool finished;

        // current recursion depth
        unsigned int rec_depth;

    protected:
        // generate (i) sequence
        typedef stxxl::stream::counter<offset_type> counter_stream_type;

        // Sorter
        typedef stxxl::tuple_less1st<skew_pair_type> mod12cmp;
        typedef stxxl::sorter<skew_pair_type, mod12cmp> mod12_sorter_type;

        // Additional streaming items
        typedef stream::choose<mod12_sorter_type, 2> isa_second_type;
        typedef build_sa<offset_array_it_rg, isa_second_type, isa_second_type> buildSA_type;
        typedef make_pairs<buildSA_type, counter_stream_type> precompute_isa_type;

        // Real recursive skew3 implementation
        // This part is the core of the skew algorithm and runs all class objects in their respective order
        template <typename RecInputType>
        buildSA_type * skew3(RecInputType& p_Input)
        {
            // (t_i) -> (i,t_i,t_{i+1},t_{i+2})
            typedef make_quads<RecInputType, offset_type, 1> make_quads_input_type;

            // (t_i) -> (i,t_i,t_{i+1},t_{i+2}) with i = 1,2 mod 3
            typedef extract_mod12<make_quads_input_type> mod12_quads_input_type;

            // sort (i,t_i,t_{i+1},t_{i+2}) by (t_i,t_{i+1},t_{i+2})
            typedef typename stream::sort<mod12_quads_input_type, less_quad<offset_type> > sort_mod12_input_type;

            // name (i,t_i,t_{i+1},t_{i+2}) -> (i,n_i)
            typedef naming<sort_mod12_input_type> naming_input_type;

            mod12_sorter_type m1_sorter(mod12cmp(), ram_use / 5);
            mod12_sorter_type m2_sorter(mod12cmp(), ram_use / 5);

            // sorted mod1 runs -concatenate- sorted mod2 runs
            typedef stxxl::stream::concatenate<mod12_sorter_type, mod12_sorter_type> concatenation_type;

            // (t_i) -> (i,t_i,t_{i+1},t_{i+2})
            offset_array_type text;
            make_quads_input_type quads_input(p_Input, text);

            // (t_i) -> (i,t_i,t_{i+1},t_{i+2}) with i = 1,2 mod 3
            mod12_quads_input_type mod12_quads_input(quads_input);

            // sort (i,t_i,t_{i+1},t_{i+2}) by (t_i,t_i+1},t_{i+2})
            sort_mod12_input_type sort_mod12_input(mod12_quads_input, less_quad<offset_type>(), ram_use / 5);

            // name (i,t_i,t_{i+1},t_{i+2}) -> (i,"n_i")
            bool unique = false;         // is the current quad array unique?
            naming_input_type names_input(sort_mod12_input, unique);

            // create (i, s^12[i])
            size_type concat_length = 0; // holds length of current S_12
            while (!names_input.empty()) {
                const skew_pair_type& tmp = *names_input;
                if (tmp.first & 1) {
                    m2_sorter.push(tmp); // sorter #2
                }
                else {
                    m1_sorter.push(tmp); // sorter #1
                }
                ++names_input;
                concat_length++;
            }

            std::cout << "recursion string length = " << concat_length << std::endl;

            m1_sorter.sort();
            m2_sorter.sort();

            if (!unique)
            {
                std::cout << "not unique -> next recursion level = " << ++rec_depth << std::endl;

                // compute s^12 := lexname[S[1 mod 3]] . lexname[S[2 mod 3]], (also known as reduced recursion string 'R')
                concatenation_type concat_mod1mod2(m1_sorter, m2_sorter);

                buildSA_type* recType = skew3(concat_mod1mod2);  // recursion with recursion string T' = concat_mod1mod2 lexnames

                std::cout << "exit recursion level = " << --rec_depth << std::endl;

                counter_stream_type isa_loop_index;
                precompute_isa_type isa_pairs(*recType, isa_loop_index); // add index as component => (SA12, i)

                // store beginning of mod2-tuples of s^12 in mod2_pos
                offset_type special = (concat_length != subp_size(text.size()));
                offset_type mod2_pos = offset_type((subp_size(text.size()) >> 1) + (subp_size(text.size()) & 1) + special);

                mod12_sorter_type isa1_pair(mod12cmp(), ram_use / 5);
                mod12_sorter_type isa2_pair(mod12cmp(), ram_use / 5);

                while (!isa_pairs.empty()) {
                    const skew_pair_type& tmp = *isa_pairs;
                    if (tmp.first < mod2_pos) {
                        if (tmp.first + special < mod2_pos) // else: special sentinel tuple is dropped
                            isa1_pair.push(tmp);            // sorter #1
                    }
                    else {
                        isa2_pair.push(tmp);                // sorter #2
                    }
                    ++isa_pairs;
                }

                delete recType;

                isa1_pair.finish();
                isa2_pair.finish();

                offset_array_it_rg input(text.begin(), text.end());

                // => (i, ISA)
                isa1_pair.sort(ram_use / 8);
                isa2_pair.sort(ram_use / 8);

                // pick ISA of (i, ISA)
                isa_second_type isa1(isa1_pair);
                isa_second_type isa2(isa2_pair);

                // prepare and run merger
                return new buildSA_type(input, isa1, isa2, text.size(), ram_use);
            }
            else // unique
            {
                std::cout << "unique names!" << std::endl;

                isa_second_type isa1(m1_sorter);
                isa_second_type isa2(m2_sorter);

                offset_array_it_rg source(text.begin(), text.end());

                // prepare and run merger
                return new buildSA_type(source, isa1, isa2, text.size(), ram_use);
            }
        } // end of skew3()

    protected:
        // Adapt (t_i) -> (i,t_i) for input to fit to recursive call
        typedef make_pairs<counter_stream_type, Input> make_pairs_input_type;

        // points to final constructed suffix array generator
        buildSA_type* out_sa;

    public:
        algorithm(Input& data_in)
            : finished(false), rec_depth(0)
        {
            // (t_i) -> (i,t_i)
            counter_stream_type dummy;
            make_pairs_input_type pairs_input(dummy, data_in);

            out_sa = skew3(pairs_input);
        }

        const value_type& operator * () const
        {
            return *(*out_sa);
        }

        algorithm& operator ++ ()
        {
            assert(out_sa);
            assert(!out_sa->empty());

            ++(*out_sa);

            if (out_sa->empty()) {
                finished = true;
                delete out_sa;
                out_sa = NULL;
            }
            return *this;
        }

        ~algorithm()
        {
            if (out_sa) delete out_sa;
        }

        bool empty() const
        {
            return finished;
        }
    }; // algorithm class
};     // skew class

//! helper to print out readable characters.
template <typename alphabet_type>
static inline std::string dumpC(alphabet_type c)
{
    std::ostringstream oss;
    if (isalnum(c)) oss << '\'' << (char)c << '\'';
    else oss << (int)c;
    return oss.str();
}

//! helper stream to cut input off a specified length.
template <typename InputType>
class cut_stream
{
public:
    //! same value type as input stream
    typedef typename InputType::value_type value_type;

protected:
    //! instance of input stream
    InputType& m_input;

    //! counter after which the stream ends
    size_type m_count;

public:
    cut_stream(InputType& input, size_type count)
        : m_input(input), m_count(count)
    { }

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
            size_type sizelimit,
            bool text_output_flag, bool check_flag, bool input_verbatim)
{
    static const size_t block_size = sizeof(offset_type) * 1024 * 1024 / 2;

    typedef typename stxxl::VECTOR_GENERATOR<alphabet_type, 1, 2>::result alphabet_vector_type;
    typedef typename stxxl::VECTOR_GENERATOR<offset_type, 1, 2, block_size>::result offset_vector_type;

    // input and output files (if supplied via command line)
    stxxl::syscall_file* input_file = NULL, * output_file = NULL;

    // input and output vectors for suffix array construction
    alphabet_vector_type input_vector;
    offset_vector_type output_vector;

    using stxxl::file;

    if (input_verbatim)
    {
        // copy input verbatim into vector
        input_vector.resize(input_filename.size());
        std::copy(input_filename.begin(), input_filename.end(), input_vector.begin());
    }
    else if (input_filename == "random")
    {
        if (sizelimit == std::numeric_limits<size_type>::max()) {
            std::cout << "You must provide -s <size> for generated inputs." << std::endl;
            return 1;
        }

        // fill input vector with random bytes
        input_vector.resize(sizelimit);
        stxxl::random_number8_r rand8;
        stxxl::generate(input_vector.begin(), input_vector.end(), rand8);
    }
    else if (input_filename == "unary")
    {
        if (sizelimit == std::numeric_limits<size_type>::max()) {
            std::cout << "You must provide -s <size> for generated inputs." << std::endl;
            return 1;
        }

        // fill input vector with random bytes
        input_vector.resize(sizelimit);
        stxxl::generate(input_vector.begin(), input_vector.end(), unary_generator);
    }
    else
    {
        // define input file object and map input_vector to input_file (no copying)
        input_file = new stxxl::syscall_file(input_filename, file::RDONLY | file::DIRECT);
        alphabet_vector_type file_input_vector(input_file);
        input_vector.swap(file_input_vector);
    }

    if (output_filename.size())
    {
        // define output file object and map output_vector to output_file
        output_file = new stxxl::syscall_file(output_filename, file::RDWR | file::CREAT | file::DIRECT);
        offset_vector_type file_output_vector(output_file);
        output_vector.swap(file_output_vector);
    }

    // I/O measurement
    stxxl::stats* Stats = stxxl::stats::get_instance();
    stxxl::stats_data stats_begin(*Stats);

    // construct skew class with bufreader input type
    typedef alphabet_vector_type::bufreader_type input_type;
    typedef cut_stream<input_type> cut_input_type;
    typedef typename skew<offset_type>::template algorithm<cut_input_type> skew_type;

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

    // write suffix array stream into output vector
    stream::materialize(skew, output_vector.begin(), output_vector.end());

    std::cout << "output size = " << output_vector.size() << std::endl;
    std::cout << (stxxl::stats_data(*Stats) - stats_begin); // print i/o statistics

    if (text_output_flag)
    {
        std::cout << std::endl << "resulting suffix array:" << std::endl;

        for (unsigned int i = 0; i < output_vector.size(); i++) {
            std::cout << i << " : " << output_vector[i] << " : ";

            // We need a const reference because operator[] would write data
            const alphabet_vector_type& cinput = input_vector;

            for (unsigned int j = 0; output_vector[i] + j < cinput.size(); j++) {
                std::cout << dumpC(cinput[output_vector[i] + j]) << " ";
            }

            std::cout << std::endl;
        }
    }

    int ret = 0;

    if (check_flag)
    {
        (std::cout << "checking suffix array... ").flush();

        if (!sacheck_vectors(input_vector, output_vector)) {
            std::cout << "failed!" << std::endl;
            ret = -1;
        }
        else
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
        "DC3 aka skew3 algorithm for external memory suffix array construction."
        );
    cp.set_author(
        "Jens Mehnert <jmehnert@mpi-sb.mpg.de>, "
        "Timo Bingmann <tb@panthema.net>, "
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
                 "Cut input text to given size, e.g. 2 GiB.");
    cp.add_bytes('M', "memuse", ram_use,
                 "Amount of RAM to use, default: 1 GiB.");
    cp.add_uint('w', "wordsize", wordsize,
                "Set word size of suffix array to 32, 40 or 64 bit, "
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
