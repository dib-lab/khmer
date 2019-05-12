############################################################################
#  tests/containers/ppq/CMakeLists.txt
#
#  Part of the STXXL. See http://stxxl.sourceforge.net
#
#  Copyright (C) 2014 Timo Bingmann <tb@panthema.net>
#
#  Distributed under the Boost Software License, Version 1.0.
#  (See accompanying file LICENSE_1_0.txt or copy at
#  http://www.boost.org/LICENSE_1_0.txt)
############################################################################

stxxl_build_test(benchmark_pqs)
stxxl_build_test(test_ppq)
stxxl_build_test(test_ppq_arrays_and_iterator)

stxxl_test(test_ppq)
stxxl_test(test_ppq_arrays_and_iterator)

if (MSVC)
 add_definitions(/bigobj)
endif()

add_define(test_ppq "STXXL_DEBUG_ASSERTIONS=1")
add_define(test_ppq_arrays_and_iterator "STXXL_DEBUG_ASSERTIONS=1")
