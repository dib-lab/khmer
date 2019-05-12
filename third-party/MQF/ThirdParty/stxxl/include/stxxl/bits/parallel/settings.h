/***************************************************************************
 *  include/stxxl/bits/parallel/settings.h
 *
 *  Settings and tuning parameters, heuristics to decide whether to use
 *  parallelized algorithms.
 *  Extracted from MCSTL http://algo2.iti.uni-karlsruhe.de/singler/mcstl/
 *
 *  Part of the STXXL. See http://stxxl.sourceforge.net
 *
 *  Copyright (C) 2007 Johannes Singler <singler@ira.uka.de>
 *  Copyright (C) 2014 Timo Bingmann <tb@panthema.net>
 *
 *  Distributed under the Boost Software License, Version 1.0.
 *  (See accompanying file LICENSE_1_0.txt or copy at
 *  http://www.boost.org/LICENSE_1_0.txt)
 **************************************************************************/

#ifndef STXXL_PARALLEL_SETTINGS_HEADER
#define STXXL_PARALLEL_SETTINGS_HEADER

#include <stxxl/bits/config.h>
#include <stxxl/bits/namespace.h>
#include <stxxl/bits/parallel/types.h>

#if STXXL_PARALLEL
  #include <omp.h>
#endif

STXXL_BEGIN_NAMESPACE

namespace parallel {

/** The extensible condition on whether the parallel variant of an algorithm sould be called.
 * \param c A condition that is overruled by stxxl::parallel::SETTINGS::force_parallel,
 * i. e. usually a decision based on the input size. */
#define STXXL_PARALLEL_CONDITION(c) (!(SETTINGS::force_sequential) && ((SETTINGS::num_threads > 1 && (c)) || SETTINGS::force_parallel))

/** Pseudo-integer that gets its initial value from \c omp_get_max_threads(). */
class NumberOfThreads
{
public:
    int num_threads;

    NumberOfThreads()
    {
#if STXXL_PARALLEL
        num_threads = omp_get_max_threads();
#else
        num_threads = 1;
#endif      // STXXL_PARALLEL
        if (num_threads < 1)
            num_threads = 1;
    }

    operator int ()
    {
        return num_threads;
    }

    int operator = (int nt)
    {
        num_threads = nt;
        return num_threads;
    }
};

/** Run-time settings for the MCSTL.
 * \param must_be_int This template parameter exists only to avoid having an
 * object file belonging to the library. It should always be set to int,
 * to ensure deterministic behavior. */
template <typename must_be_int = int>
struct Settings
{
public:
    /** Different parallel sorting algorithms to choose from: multi-way mergesort, quicksort, load-balanced quicksort. */
    enum SortAlgorithm { MWMS, QS, QS_BALANCED };
    /** Different merging algorithms: bubblesort-alike, loser-tree variants, enum sentinel */
    enum MultiwayMergeAlgorithm { BUBBLE, LOSER_TREE, LOSER_TREE_COMBINED, LOSER_TREE_SENTINEL, MWM_ALGORITHM_LAST };
    /** Different splitting strategies for sorting/merging: by sampling, exact */
    enum Splitting { SAMPLING, EXACT };

    /** Number of thread to be used.
     * Initialized to omp_get_max_threads(), but can be overridden by the user.
     */
    static NumberOfThreads num_threads;

    /** Force all algorithms to be executed sequentially.
     * This setting cannot be overwritten. */
    static volatile bool force_sequential;

    /** Force all algorithms to be executed in parallel.
     * This setting can be overriden by
     * mcstl::sequential_tag (compile-time), and
     * force_sequential (run-time). */
    static volatile bool force_parallel;

    /** Algorithm to use for sorting. */
    static volatile SortAlgorithm sort_algorithm;

    /** Strategy to use for splitting the input when sorting (MWMS). */
    static volatile Splitting sort_splitting;

    /** Minimal input size for parallel sorting. */
    static volatile sequence_index_t sort_minimal_n;
    /** Oversampling factor for parallel std::sort (MWMS). */
    static volatile unsigned int sort_mwms_oversampling;

    /** Oversampling factor for parallel std::merge.
     * Such many samples per thread are collected. */
    static volatile unsigned int merge_oversampling;

    /** Algorithm to use for parallel mcstl::multiway_merge. */
    static volatile MultiwayMergeAlgorithm multiway_merge_algorithm;
    /** Splitting strategy to use for parallel mcstl::multiway_merge. */
    static volatile Splitting multiway_merge_splitting;
    /** Oversampling factor for parallel mcstl::multiway_merge. */
    static volatile unsigned int multiway_merge_oversampling;
    /** Minimal input size for parallel mcstl::multiway_merge. */
    static volatile sequence_index_t multiway_merge_minimal_n;
    /** Oversampling factor for parallel mcstl::multiway_merge. */
    static volatile int multiway_merge_minimal_k;

//hardware dependent tuning parameters
    /** Size of the L1 cache in bytes (underestimation). */
    static volatile unsigned long long L1_cache_size;
    /** Size of the L2 cache in bytes (underestimation). */
    static volatile unsigned long long L2_cache_size;
    /** Size of the Translation Lookaside Buffer (underestimation). */
    static volatile unsigned int TLB_size;

    /** Overestimation of cache line size.
     * Used to avoid false sharing, i. e. elements of different threads are at least this amount apart. */
    static unsigned int cache_line_size;        //overestimation
};

/** Convenience typedef to avoid to have write \c Settings<>. */
typedef Settings<> SETTINGS;

template <typename must_be_int>
volatile bool Settings<must_be_int>::force_parallel = false;

template <typename must_be_int>
volatile bool Settings<must_be_int>::force_sequential = false;

template <typename must_be_int>
volatile typename Settings<must_be_int>::SortAlgorithm Settings<must_be_int>::sort_algorithm = Settings<must_be_int>::MWMS;

template <typename must_be_int>
volatile typename Settings<must_be_int>::Splitting Settings<must_be_int>::sort_splitting = Settings<must_be_int>::EXACT;

template <typename must_be_int>
volatile sequence_index_t Settings<must_be_int>::sort_minimal_n = 1000;

template <typename must_be_int>
volatile unsigned int Settings<must_be_int>::sort_mwms_oversampling = 10;

template <typename must_be_int>
volatile unsigned int Settings<must_be_int>::merge_oversampling = 10;

template <typename must_be_int>
volatile sequence_index_t Settings<must_be_int>::multiway_merge_minimal_n = 1000;

template <typename must_be_int>
volatile int Settings<must_be_int>::multiway_merge_minimal_k = 2;

template <typename must_be_int>
volatile typename Settings<must_be_int>::MultiwayMergeAlgorithm Settings<must_be_int>::multiway_merge_algorithm = Settings<must_be_int>::LOSER_TREE;

template <typename must_be_int>
volatile typename Settings<must_be_int>::Splitting Settings<must_be_int>::multiway_merge_splitting = Settings<must_be_int>::EXACT;

template <typename must_be_int>
volatile unsigned int Settings<must_be_int>::multiway_merge_oversampling = 10;

//hardware dependent tuning parameters
template <typename must_be_int>
volatile unsigned long long Settings<must_be_int>::L1_cache_size = 16 << 10;

template <typename must_be_int>
volatile unsigned long long Settings<must_be_int>::L2_cache_size = 256 << 10;

template <typename must_be_int>
volatile unsigned int Settings<must_be_int>::TLB_size = 128;

template <typename must_be_int>
unsigned int Settings<must_be_int>::cache_line_size = 64;       //overestimation

template <typename must_be_int>
NumberOfThreads Settings<must_be_int>::num_threads;

} // namespace parallel

STXXL_END_NAMESPACE

#endif // !STXXL_PARALLEL_SETTINGS_HEADER
