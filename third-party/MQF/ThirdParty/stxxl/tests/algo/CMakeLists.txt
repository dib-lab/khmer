############################################################################
#  tests/algo/CMakeLists.txt
#
#  Part of the STXXL. See http://stxxl.sourceforge.net
#
#  Copyright (C) 2013 Timo Bingmann <tb@panthema.net>
#
#  Distributed under the Boost Software License, Version 1.0.
#  (See accompanying file LICENSE_1_0.txt or copy at
#  http://www.boost.org/LICENSE_1_0.txt)
############################################################################

stxxl_build_test(test_asch)
stxxl_build_test(test_bad_cmp)
stxxl_build_test(test_ksort)
stxxl_build_test(test_random_shuffle)
stxxl_build_test(test_scan)
stxxl_build_test(test_sort)
stxxl_build_test(test_stable_ksort)

add_define(test_bad_cmp "STXXL_VERBOSE_LEVEL=0")
add_define(test_ksort "STXXL_VERBOSE_LEVEL=1" "STXXL_CHECK_ORDER_IN_SORTS")
add_define(test_random_shuffle "STXXL_VERBOSE_LEVEL=0")
add_define(test_sort "STXXL_VERBOSE_LEVEL=0")

stxxl_test(test_asch 3 100 1000 42)
stxxl_test(test_bad_cmp 16)
stxxl_test(test_ksort)
stxxl_test(test_random_shuffle)
stxxl_test(test_scan)
stxxl_test(test_sort)
stxxl_test(test_stable_ksort)

if(NOT CYGWIN AND NOT MINGW AND BUILD_EXTRAS) #-tb too big to build on cygwin

  stxxl_build_test(test_ksort_all_parameters)
  stxxl_build_test(test_stable_ksort_all_parameters)
  stxxl_build_test(test_sort_all_parameters)

  add_define(test_ksort_all_parameters "STXXL_VERBOSE_STABLE_KSORT=STXXL_VERBOSE0" "STXXL_VERBOSE_LEVEL=0" "STXXL_CHECK_ORDER_IN_SORTS")

  if(MSVC AND BUILD_TESTS)
    # requires /bigobj flag to build
    set_target_properties(test_ksort_all_parameters PROPERTIES COMPILE_FLAGS /bigobj)
    set_target_properties(test_sort_all_parameters PROPERTIES COMPILE_FLAGS /bigobj)
    set_target_properties(test_stable_ksort_all_parameters PROPERTIES COMPILE_FLAGS /bigobj)
  endif()

  ### extra sort tests:
  #test_sort_all_parameters         <MiB to sort> <MiB to use> <alloc_strategy [0..3]> <blk_size [0..14]> <seed>
  #test_ksort_all_parameters        <MiB to sort> <MiB to use> <alloc_strategy [0..3]> <blk_size [0..14]> <seed>
  #test_stable_ksort_all_parameters <MiB to sort> <MiB to use> <alloc_strategy [0..3]> <blk_size [0..14]> <seed>
  foreach(DATA 1024)
    foreach(RAM 384)
      foreach(STRATEGY 0 1 2 3)
        foreach(BLK_SIZE 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14)

          stxxl_extra_test(test_sort_all_parameters ${DATA} ${RAM} ${STRATEGY} ${BLK_SIZE} 42)
          stxxl_extra_test(test_ksort_all_parameters ${DATA} ${RAM} ${STRATEGY} ${BLK_SIZE} 42)
          #-tb not fully implemented
          #stxxl_extra_test(test_stable_ksort_all_parameters ${DATA} ${RAM} ${STRATEGY} ${BLK_SIZE} 42)

        endforeach(BLK_SIZE)
      endforeach(STRATEGY)
    endforeach(RAM)
  endforeach(DATA)

  if(NOT BUILD_EXTRAS)

    foreach(STRATEGY 0 1 2 3)
      stxxl_test(test_sort_all_parameters 512 192 ${STRATEGY} 10 42)
      stxxl_test(test_ksort_all_parameters 512 192 ${STRATEGY} 10 42)
      #-tb not fully implemented
      #stxxl_test(test_stable_ksort_all_parameters 512 192 ${STRATEGY} 10 42)
    endforeach(STRATEGY)

  endif(NOT BUILD_EXTRAS)

endif()

if(USE_GNU_PARALLEL)

  stxxl_build_test(test_parallel_sort)

  ### parallel sort tests
  #algo/test_parallel_sort [n in MiB] [p threads] [M in MiB] [sorting algorithm: m | q | qb | s] [merging algorithm: p | s | n]
  foreach(DATA 768)
    foreach(RAM 128)
      foreach(THREADS 1 2 3 4)
        foreach(SORT m q qb s)
          foreach(MERGE p s n)
            stxxl_extra_test(test_parallel_sort ${DATA} ${THREADS} ${RAM} ${SORT} ${MERGE})
          endforeach(MERGE)
        endforeach(SORT)
      endforeach(THREADS)
    endforeach(RAM)
  endforeach(DATA)

  if(NOT BUILD_EXTRAS)

    foreach(SORT m q qb s)
      foreach(MERGE p s n)
        stxxl_test(test_parallel_sort 256 4 96 ${SORT} ${MERGE})
      endforeach(MERGE)
    endforeach(SORT)

  endif(NOT BUILD_EXTRAS)

endif(USE_GNU_PARALLEL)
