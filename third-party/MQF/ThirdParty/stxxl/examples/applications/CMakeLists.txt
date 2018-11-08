############################################################################
#  examples/applications/CMakeLists.txt
#
#  Part of the STXXL. See http://stxxl.sourceforge.net
#
#  Copyright (C) 2013 Timo Bingmann <tb@panthema.net>
#
#  Distributed under the Boost Software License, Version 1.0.
#  (See accompanying file LICENSE_1_0.txt or copy at
#  http://www.boost.org/LICENSE_1_0.txt)
############################################################################

if(NOT CYGWIN AND NOT MINGW) #-tb too big to build on cygwin

  stxxl_build_example(skew3)
  stxxl_build_example(skew3-lcp)

  if(MSVC AND BUILD_TESTS)
    # requires /bigobj flag to build
    set_target_properties(skew3 PROPERTIES COMPILE_FLAGS /bigobj)
  endif()

  stxxl_test(skew3 --size 1mib random --check)
  stxxl_test(skew3-lcp --size 1mib random --check)	

endif()