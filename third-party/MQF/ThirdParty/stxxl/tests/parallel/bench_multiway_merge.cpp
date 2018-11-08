/***************************************************************************
 *  tests/parallel/bench_multiway_merge.cpp
 *
 *  Part of the STXXL. See http://stxxl.sourceforge.net
 *
 *  Copyright (C) 2014 Thomas Keh <thomas.keh@student.kit.edu>
 *  Copyright (C) 2014-2015 Timo Bingmann <tb@panthema.net>
 *
 *  Distributed under the Boost Software License, Version 1.0.
 *  (See accompanying file LICENSE_1_0.txt or copy at
 *  http://www.boost.org/LICENSE_1_0.txt)
 **************************************************************************/

#include <vector>
#include <utility>
#include <cstdlib>

#include <stxxl/bits/common/timer.h>
#include <stxxl/bits/common/rand.h>
#include <stxxl/bits/common/is_sorted.h>
#include <stxxl/bits/common/cmdline.h>
#include <stxxl/bits/parallel.h>
#include <stxxl/bits/parallel/multiway_merge.h>

using stxxl::uint64;

// number of repetitions of each benchmark
unsigned int g_outer_repeat = 3;

// number of inner repetitions of each benchmark
unsigned int g_inner_repeat = 1;

// factor to multiply while increasing number of arrays
unsigned int g_factor = 32;

// run only a few quick benchmark runs
bool g_quick = false;

struct DataStruct
{
    unsigned int key;
    char payload[32];

    explicit DataStruct(unsigned int k = 0)
        : key(k)
    { }

    bool operator < (const DataStruct& other) const
    {
        return key < other.key;
    }

    bool operator == (const DataStruct& other) const
    {
        return (key == other.key);
    }

    friend std::ostream& operator << (std::ostream& os, const DataStruct& s)
    {
        return os << '(' << s.key << ",...)";
    }
};

enum benchmark_type {
    SEQ_MWM_LT,
    SEQ_MWM_LT_STABLE,
    SEQ_MWM_LT_COMBINED,
    SEQ_GNU_MWM,

    PARA_MWM_EXACT_LT,
    PARA_MWM_EXACT_LT_STABLE,
    PARA_MWM_SAMPLING_LT,
    PARA_MWM_SAMPLING_LT_STABLE,
    PARA_GNU_MWM_EXACT,
    PARA_GNU_MWM_SAMPLING
};

template <typename ValueType, benchmark_type Method>
void test_multiway_merge(unsigned int seq_count, const size_t seq_size)
{
    // we allocate a list of blocks, each block being a sequence of items.
    static const size_t item_size = sizeof(ValueType);
    const size_t seq_items = seq_size / item_size;

    typedef std::vector<ValueType> sequence_type;

    size_t total_size = seq_count * seq_items;
    size_t total_bytes = seq_count * seq_size;

    std::less<ValueType> cmp;

    std::vector<sequence_type> seqs(seq_count);

    {
        stxxl::scoped_print_timer spt(
            "Filling sequences with random numbers in parallel", total_bytes);

        for (size_t i = 0; i < seq_count; ++i)
            seqs[i].resize(seq_items);

#if STXXL_PARALLEL
#pragma omp parallel
#endif
        {
#if STXXL_PARALLEL
            unsigned int seed = 1234 * omp_get_thread_num();
            stxxl::random_number32_r rnd(seed);
#pragma omp for
#else
            stxxl::random_number32_r rnd(1234);
#endif
            for (long i = 0; i < seq_count; ++i)
            {
                for (size_t j = 0; j < seq_items; ++j)
                    seqs[i][j] = ValueType(rnd());

                std::sort(seqs[i].begin(), seqs[i].end(), cmp);
            }
        }
    }

    std::vector<ValueType> out;

    {
        stxxl::scoped_print_timer spt("Allocating output buffer", total_bytes);
        out.resize(total_size);
    }

    {
        stxxl::scoped_print_timer spt("Merging", total_bytes);

        const char* method_name = NULL;

        typedef std::pair<typename sequence_type::iterator,
                          typename sequence_type::iterator>
            sequence_iterator_pair_type;

        std::vector<sequence_iterator_pair_type> iterpairs(seq_count);

        for (unsigned int r = 0; r < g_inner_repeat; ++r)
        {
            // (re-)set sequence iterators
            for (size_t i = 0; i < seq_count; ++i) {
                iterpairs[i] =
                    sequence_iterator_pair_type(seqs[i].begin(), seqs[i].end());
            }

            using stxxl::parallel::SETTINGS;

            switch (Method)
            {
            case SEQ_MWM_LT:
                method_name = "seq_mwm_lt";

                SETTINGS::multiway_merge_algorithm = SETTINGS::LOSER_TREE;

                stxxl::parallel::sequential_multiway_merge<false, false>(
                    iterpairs.begin(), iterpairs.end(),
                    out.begin(), total_size, cmp);
                break;

            case SEQ_MWM_LT_STABLE:
                method_name = "seq_mwm_lt_stable";

                SETTINGS::multiway_merge_algorithm = SETTINGS::LOSER_TREE;

                stxxl::parallel::sequential_multiway_merge<true, false>(
                    iterpairs.begin(), iterpairs.end(),
                    out.begin(), total_size, cmp);
                break;

            case SEQ_MWM_LT_COMBINED:
                method_name = "seq_mwm_lt_combined";

                SETTINGS::multiway_merge_algorithm = SETTINGS::LOSER_TREE_COMBINED;

                stxxl::parallel::sequential_multiway_merge<false, false>(
                    iterpairs.begin(), iterpairs.end(),
                    out.begin(), total_size, cmp);
                break;

#if STXXL_WITH_GNU_PARALLEL
            case SEQ_GNU_MWM:
                method_name = "seq_gnu_mwm";

                __gnu_parallel::multiway_merge(iterpairs.begin(), iterpairs.end(),
                                               out.begin(), total_size, cmp,
                                               __gnu_parallel::sequential_tag());
                break;
#endif          // STXXL_WITH_GNU_PARALLEL

#if STXXL_PARALLEL
            case PARA_MWM_EXACT_LT:
                method_name = "para_mwm_exact_lt";

                SETTINGS::multiway_merge_algorithm = SETTINGS::LOSER_TREE;
                SETTINGS::multiway_merge_splitting = SETTINGS::EXACT;

                stxxl::parallel::multiway_merge(
                    iterpairs.begin(), iterpairs.end(),
                    out.begin(), total_size, cmp);
                break;

            case PARA_MWM_EXACT_LT_STABLE:
                method_name = "para_mwm_exact_lt_stable";

                SETTINGS::multiway_merge_algorithm = SETTINGS::LOSER_TREE;
                SETTINGS::multiway_merge_splitting = SETTINGS::EXACT;

                stxxl::parallel::multiway_merge_stable(
                    iterpairs.begin(), iterpairs.end(),
                    out.begin(), total_size, cmp);
                break;

            case PARA_MWM_SAMPLING_LT:
                method_name = "para_mwm_sampling_lt";

                SETTINGS::multiway_merge_algorithm = SETTINGS::LOSER_TREE;
                SETTINGS::multiway_merge_splitting = SETTINGS::SAMPLING;

                stxxl::parallel::multiway_merge(
                    iterpairs.begin(), iterpairs.end(),
                    out.begin(), total_size, cmp);
                break;

            case PARA_MWM_SAMPLING_LT_STABLE:
                method_name = "para_mwm_sampling_lt_stable";

                SETTINGS::multiway_merge_algorithm = SETTINGS::LOSER_TREE;
                SETTINGS::multiway_merge_splitting = SETTINGS::SAMPLING;

                stxxl::parallel::multiway_merge_stable(
                    iterpairs.begin(), iterpairs.end(),
                    out.begin(), total_size, cmp);
                break;
#endif          // STXXL_PARALLEL

#if STXXL_WITH_GNU_PARALLEL
            case PARA_GNU_MWM_EXACT: {
                method_name = "para_gnu_mwm_exact";

                __gnu_parallel::_Settings s = __gnu_parallel::_Settings::get();
                s.multiway_merge_splitting = __gnu_parallel::EXACT;
                __gnu_parallel::_Settings::set(s);

                __gnu_parallel::multiway_merge(iterpairs.begin(), iterpairs.end(),
                                               out.begin(), total_size, cmp);
                break;
            }

            case PARA_GNU_MWM_SAMPLING: {
                method_name = "para_gnu_mwm_sampling";

                __gnu_parallel::_Settings s = __gnu_parallel::_Settings::get();
                s.multiway_merge_splitting = __gnu_parallel::SAMPLING;
                __gnu_parallel::_Settings::set(s);

                __gnu_parallel::multiway_merge(iterpairs.begin(), iterpairs.end(),
                                               out.begin(), total_size, cmp);
                break;
            }
#endif          // STXXL_WITH_GNU_PARALLEL

            default:
                STXXL_ERRMSG("Error: method " << Method << " is not available "
                             "in this compilation.");
                break;
            }
        }

        std::cout << "RESULT"
                  << " seq_count=" << seq_count
                  << " method=" << method_name
                  << " item_size=" << item_size
                  << " seq_items=" << seq_items
                  << " total_size=" << total_size
                  << " total_bytes=" << total_bytes
            #if STXXL_PARALLEL
            << " num_threads=" << omp_get_max_threads()
            #else
            << " num_threads=1"
            #endif
            << " time=" << spt.timer().seconds()
            << " inner_repeats=" << g_inner_repeat
            << " outer_repeats=" << g_outer_repeat
            << " time/item[ns]="
            << spt.timer().seconds() / g_inner_repeat / total_size * 1e9
            << std::endl;
    }

    STXXL_CHECK(stxxl::is_sorted(out.begin(), out.end(), cmp));
}

template <typename ValueType, benchmark_type Method>
void test_repeat(unsigned int seq_count, const size_t seq_size)
{
    for (unsigned int r = 0; r < g_outer_repeat; ++r)
        test_multiway_merge<ValueType, Method>(seq_count, seq_size);
}

template <typename ValueType, benchmark_type Method>
void test_seqnum(const size_t seq_size = 2* 1024* 1024)
{
    for (unsigned int s = 1; s < 64; s += 1 * g_factor)
        test_repeat<ValueType, Method>(s, seq_size);

    for (unsigned int s = 64; s < 128; s += 2 * g_factor)
        test_repeat<ValueType, Method>(s, seq_size);

    for (unsigned int s = 128; s < 256; s += 4 * g_factor)
        test_repeat<ValueType, Method>(s, seq_size);

    if (g_quick) return;

    for (unsigned int s = 256; s < 512; s += 16 * g_factor)
        test_repeat<ValueType, Method>(s, seq_size);

    for (unsigned int s = 512; s < 1024; s += 64 * g_factor)
        test_repeat<ValueType, Method>(s, seq_size);

    for (unsigned int s = 1024; s < 4 * 1024; s += 128 * g_factor)
        test_repeat<ValueType, Method>(s, seq_size);
}

template <typename ValueType>
void test_seqnum_sequential()
{
    test_seqnum<ValueType, SEQ_MWM_LT>();
    test_seqnum<ValueType, SEQ_MWM_LT_STABLE>();
    test_seqnum<ValueType, SEQ_MWM_LT_COMBINED>();
    test_seqnum<ValueType, SEQ_GNU_MWM>();
}

template <typename ValueType>
void test_seqnum_parallel()
{
    test_seqnum<ValueType, PARA_MWM_EXACT_LT>();
    // test_seqnum<ValueType, PARA_MWM_EXACT_LT_STABLE>();
    // test_seqnum<ValueType, PARA_MWM_SAMPLING_LT>();
    // test_seqnum<ValueType, PARA_MWM_SAMPLING_LT_STABLE>();
    // test_seqnum<ValueType, PARA_GNU_MWM_EXACT>();
    // test_seqnum<ValueType, PARA_GNU_MWM_SAMPLING>();
}

template <typename ValueType, benchmark_type Method>
void test_seqsize(const size_t seq_num = 64)
{
    for (size_t b = 1024; b < 1024 * 1024 / seq_num; b *= 2)
        test_repeat<ValueType, Method>(seq_num, b);

    for (size_t b = 1024 * 1024 / seq_num;
         b < 1024 * 1024 * (size_t)1024 / seq_num; b *= 2)
        test_repeat<ValueType, Method>(seq_num, b);
}

int main(int argc, char* argv[])
{
    std::string benchset;

    stxxl::cmdline_parser cp;
    cp.set_description("STXXL multiway_merge benchmark");

    cp.add_param_string("sequ/para/both", benchset,
                        "benchmark set: sequ(ential), para(llel) or both");

    cp.add_uint('r', "inner-repeat", g_inner_repeat,
                "number of inner repetitions within each benchmark");

    cp.add_uint('R', "outer-repeat", g_outer_repeat,
                "number of repetitions of each benchmark");

    cp.add_uint('f', "factor", g_factor,
                "factor to multiply while increasing number of arrays");

    cp.add_flag('q', "quick", g_quick,
                "run only a few quick benchmark runs");

    if (!cp.process(argc, argv))
        return EXIT_FAILURE;

    // run individually for debugging
    if (0)
    {
        test_repeat<uint64, PARA_GNU_MWM_EXACT>(2, 2 * 1024 * 1024);
        test_repeat<uint64, PARA_MWM_EXACT_LT>(2, 2 * 1024 * 1024);
        //test_repeat<uint64, PARA_MWM_EXACT_LT>(1024, 2 * 1024 * 1024);
        //test_repeat<uint64, SEQ_MWM_LT>(2, 2 * 1024 * 1024);
        //test_repeat<uint64, SEQ_GNU_MWM>(2, 2 * 1024 * 1024);
        return 0;
    }

    if (benchset == "sequ" || benchset == "sequential" || benchset == "both")
    {
        test_seqnum_sequential<uint64>();
        test_seqnum_sequential<DataStruct>();
    }
    if (benchset == "para" || benchset == "parallel" || benchset == "both")
    {
        test_seqnum_parallel<uint64>();
        test_seqnum_parallel<DataStruct>();
    }
    if (benchset == "vecsize")
    {
        test_seqsize<uint64, PARA_MWM_EXACT_LT>();
        test_seqsize<DataStruct, PARA_MWM_EXACT_LT>();
    }

    return 0;
}
