############################################################################
#  tests/containers/CMakeLists.txt
#
#  Part of the STXXL. See http://stxxl.sourceforge.net
#
#  Copyright (C) 2013 Timo Bingmann <tb@panthema.net>
#
#  Distributed under the Boost Software License, Version 1.0.
#  (See accompanying file LICENSE_1_0.txt or copy at
#  http://www.boost.org/LICENSE_1_0.txt)
############################################################################

add_subdirectory(btree)
add_subdirectory(hash_map)
add_subdirectory(ppq)

stxxl_build_test(test_deque)
stxxl_build_test(test_ext_merger)
stxxl_build_test(test_ext_merger2)
stxxl_build_test(test_iterators)
stxxl_build_test(test_many_stacks)
stxxl_build_test(test_matrix)
stxxl_build_test(test_migr_stack)
stxxl_build_test(test_pqueue)
stxxl_build_test(test_queue)
stxxl_build_test(test_queue2)
stxxl_build_test(test_sequence)
stxxl_build_test(test_sorter)
stxxl_build_test(test_stack)
stxxl_build_test(test_vector)
stxxl_build_test(test_vector_buf)
stxxl_build_test(test_vector_export)
stxxl_build_test(test_vector_resize)
stxxl_build_test(test_vector_sizes)

stxxl_test(test_deque 3333333)
stxxl_test(test_ext_merger)
stxxl_test(test_ext_merger2)
stxxl_test(test_iterators)
stxxl_test(test_many_stacks 42)
stxxl_test(test_matrix)
stxxl_extra_test(test_matrix --rank 2000)
stxxl_test(test_migr_stack)
stxxl_test(test_pqueue)
stxxl_test(test_queue)
stxxl_test(test_queue2 200)
stxxl_test(test_sequence)
stxxl_test(test_sorter)
stxxl_test(test_stack 1024)
stxxl_test(test_vector)
stxxl_test(test_vector_buf)
stxxl_test(test_vector_export)
stxxl_test(test_vector_resize)
stxxl_test(test_vector_sizes "${STXXL_TMPDIR}/out" syscall)
if(STXXL_HAVE_MMAP_FILE)
  stxxl_test(test_vector_sizes "${STXXL_TMPDIR}/out" mmap)
endif(STXXL_HAVE_MMAP_FILE)
if(USE_BOOST)
  #-tb: boostfd does not work right with vector, because boostfd does not
  #correctly support set_size() (truncate and extending of files). FIXME
  #stxxl_test(test_vector_sizes "${STXXL_TMPDIR}/out" boostfd)
endif(USE_BOOST)

# TESTS_MAP
stxxl_build_test(test_map)
stxxl_build_test(test_map_random)

stxxl_test(test_map 8)
stxxl_test(test_map_random 2000)

#-tb longer test for map
stxxl_extra_test(test_map 16)
