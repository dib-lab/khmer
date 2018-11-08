############################################################################
#  tests/stream/CMakeLists.txt
#
#  Part of the STXXL. See http://stxxl.sourceforge.net
#
#  Copyright (C) 2013 Timo Bingmann <tb@panthema.net>
#
#  Distributed under the Boost Software License, Version 1.0.
#  (See accompanying file LICENSE_1_0.txt or copy at
#  http://www.boost.org/LICENSE_1_0.txt)
############################################################################

stxxl_build_test(test_loop)
stxxl_build_test(test_materialize)
stxxl_build_test(test_naive_transpose)
stxxl_build_test(test_push_sort)
stxxl_build_test(test_sorted_runs)
stxxl_build_test(test_stream)
stxxl_build_test(test_stream1)

add_define(test_stream1 "STXXL_VERBOSE_LEVEL=1")
add_define(test_push_sort "STXXL_VERBOSE_LEVEL=0")
add_define(test_stream "STXXL_VERBOSE_LEVEL=1")
add_define(test_sorted_runs "STXXL_VERBOSE_LEVEL=0")
add_define(test_materialize "STXXL_VERBOSE_LEVEL=0" "STXXL_VERBOSE_MATERIALIZE=STXXL_VERBOSE0")

stxxl_test(test_loop 100 -v)
stxxl_test(test_loop 1000000)
stxxl_test(test_materialize)
stxxl_test(test_naive_transpose)
stxxl_test(test_push_sort)
stxxl_test(test_sorted_runs)
stxxl_test(test_stream)
stxxl_test(test_stream1)
