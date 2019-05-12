############################################################################
#  tools/benchmarks/CMakeLists.txt
#
#  Part of the STXXL. See http://stxxl.sourceforge.net
#
#  Copyright (C) 2013 Timo Bingmann <tb@panthema.net>
#
#  Distributed under the Boost Software License, Version 1.0.
#  (See accompanying file LICENSE_1_0.txt or copy at
#  http://www.boost.org/LICENSE_1_0.txt)
############################################################################

stxxl_build_test(benchmark_naive_matrix)
stxxl_build_test(matrix_benchmark)
stxxl_build_test(monotonic_pq)
stxxl_build_test(pq_benchmark)
stxxl_build_test(stack_benchmark)

add_define(benchmark_naive_matrix "STXXL_VERBOSE_LEVEL=0")

if(USE_TPIE AND BUILD_TESTS)
  file(GLOB TPIE_DIR_GLOB RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} "tpie_*")

  foreach(TPIE_TEST ${TPIE_DIR_GLOB})
    if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/${TPIE_TEST}/lib/libtpie.a")
      set(TPIE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/${TPIE_TEST}")
      break()
    endif()
  endforeach()

  if(EXISTS "${TPIE_DIR}/lib/libtpie.a")

    message(STATUS "Found compiled TPIE library at ${TPIE_DIR}")

    add_library(tpie STATIC IMPORTED)
    set_property(TARGET tpie PROPERTY IMPORTED_LOCATION "${TPIE_DIR}/lib/libtpie.a")

    include_directories("${TPIE_DIR}/include")

    stxxl_build_test(tpie_stack_benchmark)
    target_link_libraries(tpie_stack_benchmark ${STXXL_LIBRARIES} tpie)

    # try to locate a Berkeley DB version
    find_path(BERKELEYDB_INCLUDES db_cxx.h
      /usr/local/include/db4.8
      /usr/local/include/db4.7
      /usr/local/include/db4.6
      /usr/local/include/db4.5
      /usr/local/include/db4
      /usr/local/include/db
      /usr/include/db4.8
      /usr/include/db4.7
      /usr/include/db4.6
      /usr/include/db4.5
      /usr/include/db4
      /usr/include/db
      )

    find_library(BERKELEYDB_LIBRARIES NAMES db_cxx)

    include(FindPackageHandleStandardArgs)
    find_package_handle_standard_args(BERKELEYDB "Could not find Berkeley DB"
      BERKELEYDB_INCLUDES BERKELEYDB_LIBRARIES)
    mark_as_advanced(BERKELEYDB_INCLUDES BERKELEYDB_LIBRARIES)

    if(BERKELEYDB_FOUND)
      include_directories("${BERKELEYDB_INCLUDES}")

      stxxl_build_test(berkeley_db_benchmark)
      target_link_libraries(berkeley_db_benchmark ${STXXL_LIBRARIES} tpie ${BERKELEYDB_LIBRARIES})
    endif()
  else()
    message(SEND_ERROR "Could not find TPIE for extra benchmarks.")
  endif()
endif(USE_TPIE AND BUILD_TESTS)
