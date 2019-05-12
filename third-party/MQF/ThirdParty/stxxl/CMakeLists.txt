############################################################################
#  CMakeLists.txt
#
#  Base CMake file for building stxxl with different options.
#
#  Part of the STXXL. See http://stxxl.sourceforge.net
#
#  Copyright (C) 2013-2014 Timo Bingmann <tb@panthema.net>
#
#  Distributed under the Boost Software License, Version 1.0.
#  (See accompanying file LICENSE_1_0.txt or copy at
#  http://www.boost.org/LICENSE_1_0.txt)
############################################################################

set(CMAKE_LEGACY_CYGWIN_WIN32 0) # needed to silence warning on cygwin

# require cmake 2.6.2 (but please use 2.8.x)
cmake_minimum_required(VERSION 2.6.2 FATAL_ERROR)

# the STXXL project
project(stxxl)

# for additional cmake scripts
set(CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/misc/cmake)

# prohibit in-source builds
if("${PROJECT_SOURCE_DIR}" STREQUAL "${PROJECT_BINARY_DIR}")
  message(FATAL_ERROR "In-source builds are not allowed, use a separate build directory.")
endif()

# default to Debug building for single-config generators
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
  message("Defaulting CMAKE_BUILD_TYPE to Debug")
  set(CMAKE_BUILD_TYPE "Debug")
endif()

# STXXL version string
set(STXXL_VERSION_MAJOR "1")
set(STXXL_VERSION_MINOR "4")
set(STXXL_VERSION_PATCH "99")
set(STXXL_VERSION_STRING "${STXXL_VERSION_MAJOR}.${STXXL_VERSION_MINOR}.${STXXL_VERSION_PATCH}")
set(STXXL_VERSION_PHASE "prerelease/${CMAKE_BUILD_TYPE}")

# read git directory (if it exists) and find git sha
if(EXISTS ${PROJECT_SOURCE_DIR}/.git)
  find_package(Git)
  if(GIT_FOUND)
    execute_process(COMMAND ${GIT_EXECUTABLE} describe HEAD
      WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}"
      OUTPUT_VARIABLE "STXXL_VERSION_GIT_REFSPEC"
      ERROR_QUIET
      OUTPUT_STRIP_TRAILING_WHITESPACE)
    execute_process(COMMAND ${GIT_EXECUTABLE} rev-parse HEAD
      WORKING_DIRECTORY "${PROJECT_SOURCE_DIR}"
      OUTPUT_VARIABLE "STXXL_VERSION_GIT_SHA1"
      ERROR_QUIET
      OUTPUT_STRIP_TRAILING_WHITESPACE)
    message(STATUS "Detected git refspec ${STXXL_VERSION_GIT_REFSPEC} sha ${STXXL_VERSION_GIT_SHA1}")
  endif()
endif()

###############################################################################
# compilation options

# please document all build options in doc/install.doc:install_build_options

option(BUILD_EXAMPLES "Build all stxxl examples" OFF)

option(BUILD_TESTS "Build all stxxl test programs" OFF)
if(BUILD_TESTS)
  set(BUILD_EXAMPLES ON)
endif()

option(BUILD_EXTRAS "Build all extra stxxl tests and tool programs" OFF)

if(MSVC_VERSION LESS 1700)
  # require boost for Visual C++ versions older than VC11 = MSVS 2012
  set(USE_BOOST ON)
else()
  option(USE_BOOST "Use Boost libraries for threads,config,filesystem,random and date_time" OFF)
endif()

option(USE_OPENMP "Use OpenMP for multi-core parallelism" ON)

if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
  option(USE_GNU_PARALLEL "Use GNU parallel extensions for multi-core parallelism" ON)

  # GNU parallel mode but NOT OpenMP is impossible.
  if(USE_GNU_PARALLEL AND NOT USE_OPENMP)
    set(USE_GNU_PARALLEL OFF)
  endif()
endif()

option(TRY_COMPILE_HEADERS "Test stxxl header files for self-sufficiency: try to compile them." OFF)

if(NOT MSVC)
  option(USE_STD_THREADS "Force usage of C++ standard library thread support." OFF)
endif()

option(NO_CXX11 "Build without C++11 flags" OFF)

option(USE_VALGRIND "Run tests with valgrind, pre-initialize some memory buffers." OFF)

option(USE_MALLOC_COUNT "Compile all programs with included malloc_count heap profiler." OFF)

option(USE_GCOV "Compile and run tests with gcov for coverage analysis." OFF)

option(STXXL_DEBUG_ASSERTIONS "Enable more costly assertions for internal debugging." OFF)
if(STXXL_DEBUG_ASSERTIONS)
  set(STXXL_DEBUG_ASSERTIONS "1") # change from ON/OFF to 1/0
endif()

# see tools/benchmarks about older TPIE benchmarks.
option(USE_TPIE "Try to compile extra benchmarks from the 2007 S&PE paper with an old TPIE version." OFF)

### building shared and/or static libraries

# by default we currently only build a static library, since we do not aim to
# keep a stable binary interface.
option(BUILD_STATIC_LIBS "Build static library version of libstxxl" ON)
option(BUILD_SHARED_LIBS "Build shared library version of libstxxl" OFF)

### allow user to specify other installation paths

set(INSTALL_BIN_DIR "bin" CACHE PATH "Installation directory for executables")
set(INSTALL_LIB_DIR "lib" CACHE PATH "Installation directory for libraries")
set(INSTALL_INCLUDE_DIR "include" CACHE PATH "Installation directory for header files")
set(INSTALL_PKGCONFIG_DIR "lib/pkgconfig" CACHE PATH "Installation directory for pkg-config file")

if(WIN32 AND NOT CYGWIN)
  set(DEF_INSTALL_CMAKE_DIR "CMake")
else()
  set(DEF_INSTALL_CMAKE_DIR "lib/cmake/stxxl")
endif()
set(INSTALL_CMAKE_DIR "${DEF_INSTALL_CMAKE_DIR}" CACHE PATH "Installation directory for cmake files")

###############################################################################
# enable use of "make test"

enable_testing()

include(CTest)
set(CTEST_PROJECT_NAME "STXXL")

if(USE_VALGRIND)
  set(STXXL_WITH_VALGRIND 1)
  set(VALGRIND_OPTS --leak-check=full --error-exitcode=1 --suppressions=${PROJECT_SOURCE_DIR}/misc/valgrind.supp)
endif()

###############################################################################
# enable use of "make package"

# general package settings
set(CPACK_PACKAGE_VENDOR "STXXL Maintainer Team")
set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "Library of external memory algorithms")
set(CPACK_PACKAGE_DESCRIPTION_FILE "${CMAKE_CURRENT_SOURCE_DIR}/README")
set(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_CURRENT_SOURCE_DIR}/LICENSE_1_0.txt")
set(CPACK_PACKAGE_VERSION_MAJOR "${STXXL_VERSION_MAJOR}")
set(CPACK_PACKAGE_VERSION_MINOR "${STXXL_VERSION_MINOR}")
set(CPACK_PACKAGE_VERSION_PATCH "${STXXL_VERSION_PATCH}")

# source packaging
set(CPACK_SOURCE_GENERATOR "TGZ")
set(CPACK_SOURCE_PACKAGE_FILE_NAME "${CMAKE_PROJECT_NAME}-${CPACK_PACKAGE_VERSION_MAJOR}.${CPACK_PACKAGE_VERSION_MINOR}.${CPACK_PACKAGE_VERSION_PATCH}")
set(CPACK_SOURCE_IGNORE_FILES "/\\\\..*$;/build;${CPACK_SOURCE_IGNORE_FILES}")

# binary packaging
set(CPACK_GENERATOR "TGZ")

include(CPack)

###############################################################################
# check platform

if(WIN32 OR WIN64 OR MINGW)
  set(STXXL_WINDOWS "1")
endif()

if(MSVC)
  set(STXXL_MSVC ${MSVC_VERSION})
endif()

###############################################################################
# enable warnings

if(MSVC)
  # Force to always compile with W4
  if(CMAKE_CXX_FLAGS MATCHES "/W[0-4]")
    string(REGEX REPLACE "/W[0-4]" "/W4" CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")
  else()
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /W4")
  endif()
  ### disable verbose warnings:
  # warning C4127: conditional expression is constant
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /wd4127")
  # warning C4290: C++ exception specification ignored except to indicate a function is not __declspec(nothrow)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /wd4290")
  # warning C4250: '...' : inherits '...' via dominance
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /wd4250")
  # warning C4512: assignment operator could not be generated (contains const members)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /wd4512")
  # warning C4355: 'this' : used in base member initializer list
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /wd4355")
  # disable lots of warnings about "unsecure" C runtime function
  add_definitions(-D_CRT_SECURE_NO_WARNINGS)
  # disable lots of warnings about deprecated POSIX function names
  add_definitions(-D_CRT_NONSTDC_NO_DEPRECATE)
  # disable lots of warnings about "unsecure" STL functions
  add_definitions(-D_SCL_SECURE_NO_WARNINGS)
else()
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -W -Wall -pedantic -Wno-long-long")
  # detect -Wextra
  include(CheckCXXCompilerFlag)
  check_cxx_compiler_flag(-Wextra CXX_HAS_FLAGS_WEXTRA)
  if(CXX_HAS_FLAGS_WEXTRA)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wextra")
  endif()
endif()

###############################################################################
# enable C++11 and more compiler features

# more template depth for some tests and compilers
include(CheckCXXCompilerFlag)
check_cxx_compiler_flag(-ftemplate-depth=1024 CXX_HAS_TEMPLATE_DEPTH)
if(CXX_HAS_TEMPLATE_DEPTH)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -ftemplate-depth=1024")
endif()

# enable C++11
if(MSVC)
  # MSVC 11 or greater has C++11 automatically enabled
elseif(CYGWIN)
  # C++11 is very restrictive on Cygwin
elseif(NOT NO_CXX11)
  # try -std= flags to enable C++11
  check_cxx_compiler_flag(-std=c++11 CXX_HAS_STD_CXX11)
  if(CXX_HAS_STD_CXX11)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
  else()
    check_cxx_compiler_flag(-std=c++0x CXX_HAS_STD_CXX0X)
    if(CXX_HAS_STD_CXX0X)
      set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++0x")
    endif()
  endif()
  # on MacOSX with clang we need to use libc++ for C++11 headers
  if(APPLE)
    if (CMAKE_CXX_COMPILER MATCHES ".*clang[+][+]"
        OR CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
      check_cxx_compiler_flag(-stdlib=libc++ CXX_HAS_STDLIB_LIBCXX)
      if(CXX_HAS_STDLIB_LIBCXX)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -stdlib=libc++")
      else()
        message(SEND_ERROR "Compilation on MacOSX with clang requires libc++.")
      endif()
    endif()
  endif(APPLE)
endif()

# check C++ compiler for C++11 features

include(CheckCXXSourceCompiles)
check_cxx_source_compiles(
  "#include <vector>
  int main() { std::vector<int> v(42); for (auto i : v) { ++i; } return 0; }"
  STXXL_HAVE_CXX11_RANGE_FOR_LOOP)

###############################################################################
# enable gcov coverage analysis with gcc

if(USE_GCOV)
  # find programs
  find_program(GENHTML genhtml)
  find_program(LCOV lcov)

  if(NOT LCOV OR NOT GENHTML)
    message(SEND_ERROR "Coverage analysis requires lcov and genhtml programs.")
  endif()

  # add coverage anaylsis compile and link flags
  set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fprofile-arcs -ftest-coverage")
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fprofile-arcs -ftest-coverage")
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -lgcov")

  # add cached variable containing parameters for lcov/genhtml
  set(LCOV_FLAGS "" CACHE STRING "parameters for lcov")
  set(GENHTML_FLAGS --legend --no-branch-coverage
    CACHE STRING "parameters for genhtml")

  # custom target to run before tests
  add_custom_target(lcov-reset
    COMMAND ${LCOV} -q --directory ${CMAKE_BINARY_DIR} --zerocounters
    COMMENT "Resetting code coverage counters")

  # custom lcov target to run tests
  add_custom_target(lcov-runtests
    COMMAND ${CMAKE_CTEST_COMMAND} \${ARGS}
    DEPENDS ${ALL_TESTS} lcov-reset
    COMMENT "Running all unit tests")

  # get git version description
  execute_process(COMMAND git describe --tags
    WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
    OUTPUT_VARIABLE GITDESC
    OUTPUT_STRIP_TRAILING_WHITESPACE)

  # command sequence to gather, clean and generate HTML coverage report
  add_custom_target(lcov-html
    COMMAND ${LCOV} -q --directory . --capture --output-file lcov.info
    COMMAND ${LCOV} -q --remove lcov.info '/usr/*' ${LCOV_FLAGS} --output-file lcov-clean.info
    COMMAND ${GENHTML} -q -o coverage --title "STXXL ${GITDESC}" --prefix ${PROJECT_SOURCE_DIR} ${GENHTML_FLAGS} lcov-clean.info
    COMMENT "Capturing code coverage counters and create HTML coverage report"
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR})

  # top-level target to run tests and generate coverage report
  add_custom_target(test-coverage
    COMMENT "Generate HTML coverage report "
    DEPENDS lcov-runtests lcov-html)

endif(USE_GCOV)

###############################################################################
# find thread and random library

# for testing for c++ system include files
include(CheckIncludeFileCXX)

check_include_file_cxx(pthread.h HAVE_PTHREAD_H)

if(MINGW AND NOT HAVE_PTHREAD_H)
  set(USE_STD_THREADS ON)
endif()

if(MSVC OR USE_STD_THREADS)

  # check for std::mutex and std::threads avalability
  check_cxx_source_compiles(
    "#include <mutex>
    int main() { std::mutex mutex; mutex.lock(); return 0; }"
    HAVE_STD_MUTEX)

  check_cxx_source_compiles(
    "#include <thread>
    int main() { std::thread t; return 0; }"
    HAVE_STD_THREAD)

  if(HAVE_STD_THREAD AND HAVE_STD_MUTEX)
    set(STXXL_STD_THREADS "1")
  else()
    set(USE_BOOST ON)
    message("No std::thread and std::mutex found, trying to use Boost classes instead")
  endif()

  # using <thread> also requires -pthread on gcc
  if(CMAKE_COMPILER_IS_GNUCXX)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pthread")
  endif()

  # also check for <random>
  check_include_file_cxx(random STXXL_STD_RANDOM)

else()
  find_package(Threads REQUIRED)
  set(STXXL_POSIX_THREADS "1")
endif()

###############################################################################
# determine large file support
#     note: these flags need not be exported, all file operations must be
#     encapsulated within the library!

include(TestLargeFiles)
test_large_files(HAVE_LARGEFILES)

if(HAVE_LARGEFILES)
  add_definitions(-D_LARGEFILE_SOURCE -D_FILE_OFFSET_BITS=64 -D_LARGE_FILES)
else()
  message(FATAL_ERROR "Large file support was not detectable.")
endif()

###############################################################################
# check for O_DIRECT flag

if(CYGWIN)
  #-tb O_DIRECT messes up cygwin
  set(STXXL_DIRECT_IO_OFF 1)
elseif(MSVC OR MINGW)
  # have FILE_FLAG_NO_BUFFERING on Windows
  set(STXXL_DIRECT_IO_OFF 0)
elseif(APPLE)
  # we do not have to test anything, all is well.
  set(STXXL_DIRECT_IO_OFF 0)
else()

  include(CheckCXXSourceCompiles)
  check_cxx_source_compiles("
#include <unistd.h>
#include <fcntl.h>
int main() { return ((int)O_DIRECT) != 0; }
" STXXL_HAVE_O_DIRECT)

  if(STXXL_HAVE_O_DIRECT)
    set(STXXL_DIRECT_IO_OFF 0)
  else()
    set(STXXL_DIRECT_IO_OFF 1)
  endif()

endif()

if(STXXL_DIRECT_IO_OFF)
  message(WARNING
    "The open() flag O_DIRECT was disabled!\n"
    "This means that STXXL cannot bypass the system cache on your system, which may degrade performance.")
endif()

###############################################################################
# check for mmap() function in sys/mman.h

include(CheckSymbolExists)
check_symbol_exists(mmap "sys/mman.h" STXXL_HAVE_MMAP_FILE)

###############################################################################
# check for Linux aio syscalls

include(CheckCXXSourceCompiles)
check_cxx_source_compiles(
  "#include <unistd.h>
   #include <sys/syscall.h>
   #include <linux/aio_abi.h>
   int main() {
       aio_context_t context;
       long r = syscall(SYS_io_setup, 5, &context);
       return (r == 0) ? 0 : -1;
   }"
   STXXL_HAVE_LINUXAIO_FILE)

###############################################################################
# check for an atomic add-and-fetch intrinsic for counting_ptr

include(CheckCXXSourceCompiles)
check_cxx_source_compiles(
  "int main() { int x; __sync_add_and_fetch(&x, +1); return 0; }"
  STXXL_HAVE_SYNC_ADD_AND_FETCH)

###############################################################################
# optional Boost libraries

if(USE_BOOST)

  set(Boost_USE_MULTITHREADED ON)

  if(WIN32)
    set(Boost_USE_STATIC_LIBS ON)
  endif()

  # first try to find the version
  find_package(Boost 1.34.1 REQUIRED)

  if(Boost_VERSION GREATER 104699)
    #-tb Boost >= 1.47.0 found, we must also link with libboost_chrono for MSVC to work.
    find_package(Boost 1.47.0 REQUIRED COMPONENTS thread date_time chrono iostreams system filesystem)
  elseif(Boost_VERSION GREATER 103499)
    #-tb Boost >= 1.35.0 found, we must also link with libboost_system for MSVC to work.
    find_package(Boost 1.35.0 REQUIRED COMPONENTS thread date_time iostreams system filesystem)
  else()
    find_package(Boost 1.34.1 REQUIRED COMPONENTS thread date_time iostreams filesystem)
  endif()

  # force caching of variables to display them to the user
  set(BOOST_ROOT "${BOOST_ROOT}" CACHE PATH "Boost installation root." FORCE)

  if(Boost_FOUND)
    # globally add boost include directories
    include_directories(${Boost_INCLUDE_DIRS})
    # set defines in <stxxl/bits/config.h>
    set(STXXL_BOOST_CONFIG "1")
    set(STXXL_BOOST_FILESYSTEM "1")
    set(STXXL_BOOST_RANDOM "1")
    set(STXXL_BOOST_THREADS "1")
    set(STXXL_BOOST_TIMESTAMP "1")
  elseif(MSVC_VERSION LESS 1700)
    message(FATAL_ERROR "Boost libraries are required for MSVC < 2012.")
  else()
    message(FATAL_ERROR "Boost libraries not found. Try compilation without them.")
  endif()

endif(USE_BOOST)

###############################################################################
# optional OpenMP parallelization

if(USE_GNU_PARALLEL AND NOT USE_OPENMP)
  message(FATAL_ERROR "Configuration selects GNU parallel mode but not OpenMP."
    "This is impossible.")
endif()

if(USE_OPENMP OR USE_GNU_PARALLEL)

  include(FindOpenMP)

  if(NOT OPENMP_FOUND)
    message(STATUS "OpenMP not found. Continuing without parallel algorithm support.")
  else()
    set(STXXL_PARALLEL 1)
    message(STATUS "OpenMP found, enabling built-in parallel algorithms.")

    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
  endif()

else()

  message(STATUS "OpenMP disabled in STXXL (no parallelism is used).")
  set(STXXL_PARALLEL 0)

endif(USE_OPENMP OR USE_GNU_PARALLEL)

###############################################################################
# optional GNU parallel STL mode

if(USE_GNU_PARALLEL)

  if(NOT OPENMP_FOUND)
    message(FATAL_ERROR "OpenMP not found. Try compilation without GNU parallel mode.")
  else()
    check_include_file_cxx(parallel/algorithm HAVE_PARALLEL_ALGORITHM_H)

    if(NOT HAVE_PARALLEL_ALGORITHM_H)
      message(FATAL_ERROR "GNU parallel mode header not found. Try compilation without parallel mode.")
    else()
      set(STXXL_WITH_GNU_PARALLEL 1)
    endif()
  endif()

endif(USE_GNU_PARALLEL)

###############################################################################
# check if one thread and random library is available

if(STXXL_STD_THREADS)
  message(STATUS "Using std::thread and other C++11 library functions.")
  set(STXXL_BOOST_THREADS 0)
  set(STXXL_POSIX_THREADS 0)
elseif(STXXL_BOOST_THREADS)
  message(STATUS "Using boost::thread library functions.")
  set(STXXL_STD_THREADS 0)
  set(STXXL_POSIX_THREADS 0)
elseif(STXXL_POSIX_THREADS)
  message(STATUS "Using POSIX pthread library functions.")
  set(STXXL_STD_THREADS 0)
  set(STXXL_BOOST_THREADS 0)
else()
  message(SEND_ERROR "Could not detect a thread library. Check the compilation documentation.")
endif()

###############################################################################
# test for additional includes and features used by some stxxl_tool components

include(CheckSymbolExists)
check_symbol_exists(mallinfo "malloc.h" STXXL_HAVE_MALLINFO_PROTO)
check_symbol_exists(mlock "sys/mman.h" STXXL_HAVE_MLOCK_PROTO)

if(USE_MALLOC_COUNT)
  # malloc_count requires the linux dl (dynamic linker) library
  find_library(DL_LIBRARIES NAMES dl)
  set(STXXL_EXTRA_LIBRARIES ${STXXL_EXTRA_LIBRARIES} ${DL_LIBRARIES})
endif()

###############################################################################
# configure environment for building

# create config.h with define switches _inside binary dir_!
configure_file(${PROJECT_SOURCE_DIR}/include/stxxl/bits/config.h.in ${PROJECT_BINARY_DIR}/include/stxxl/bits/config.h)

# CXX_FLAGS required for compilation
set(STXXL_CXX_FLAGS ${OpenMP_CXX_FLAGS})

# globally adds top-level include directories
set(STXXL_INCLUDE_DIRS ${PROJECT_SOURCE_DIR}/include/ ${PROJECT_BINARY_DIR}/include/)
include_directories(${STXXL_INCLUDE_DIRS})

# for targets using stxxl library
set(STXXL_EXTRA_LIBRARIES ${STXXL_EXTRA_LIBRARIES} ${Boost_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT})
set(STXXL_LIBRARIES stxxl ${STXXL_EXTRA_LIBRARIES})

# export STXXL_INCLUDE_DIRS and STXXL_LIBRARIES to global CACHE
set(STXXL_CXX_FLAGS "${STXXL_CXX_FLAGS}" CACHE STRING "Compiler flags for STXXL")
set(STXXL_INCLUDE_DIRS "${STXXL_INCLUDE_DIRS}" CACHE STRING "Include paths for STXXL")
set(STXXL_LIBRARIES "${STXXL_LIBRARIES}" CACHE STRING "Libraries to link for STXXL")

# build libstxxl in /lib
add_subdirectory(lib)

###############################################################################
# macros for building stxxl programs and tests

# macro for building stxxl programs
macro(stxxl_build_tool PROGNAME)

  add_executable(${PROGNAME} ${PROGNAME}.cpp ${ARGN})
  target_link_libraries(${PROGNAME} ${STXXL_LIBRARIES})

endmacro(stxxl_build_tool)

# macro for building stxxl examples
macro(stxxl_build_example TESTNAME)

  if(BUILD_EXAMPLES)
    stxxl_build_tool(${TESTNAME} ${ARGN})
  endif(BUILD_EXAMPLES)

endmacro(stxxl_build_example)

# macro for building stxxl tests
macro(stxxl_build_test TESTNAME)

  if(BUILD_TESTS)
    stxxl_build_tool(${TESTNAME} ${ARGN})
  endif(BUILD_TESTS)

endmacro(stxxl_build_test)

# macro for registering stxxl tests
macro(stxxl_test TESTNAME)

  if(BUILD_TESTS)
    set(TESTFULLNAME ${TESTNAME} ${ARGN})
    string(REPLACE ";" "_" TESTFULLNAME "${TESTFULLNAME}") # stringify list

    if(USE_VALGRIND)
      # prepend valgrind call
      add_test(${TESTFULLNAME} /usr/bin/valgrind ${VALGRIND_OPTS} ./${TESTNAME} ${ARGN})
      set_tests_properties(${TESTFULLNAME} PROPERTIES TIMEOUT 7200)
    else()
      add_test(${TESTFULLNAME} ${TESTNAME} ${ARGN})
      set_tests_properties(${TESTFULLNAME} PROPERTIES TIMEOUT 3600)
    endif()

  endif(BUILD_TESTS)

endmacro(stxxl_test)

# macro for building stxxl extra program
macro(stxxl_build_extra_tool)

  if(BUILD_EXTRAS)
    stxxl_build_tool(${ARGN})
  endif(BUILD_EXTRAS)

endmacro(stxxl_build_extra_tool)

# macro for registering extra stxxl tests
macro(stxxl_extra_test)

  if(BUILD_EXTRAS)
    stxxl_test(${ARGN})
  endif(BUILD_EXTRAS)

endmacro(stxxl_extra_test)

# macro for setting additional defines for targets
macro(add_define PROGNAME)

  if(TARGET ${PROGNAME})
    set_property(TARGET ${PROGNAME} APPEND PROPERTY COMPILE_DEFINITIONS ${ARGN})
  else()
    if(BUILD_TESTS)
      message("Ignoring add_define(${PROGNAME} ${ARGN}) for unknown target ${PROGNAME}")
    endif(BUILD_TESTS)
  endif()

endmacro(add_define TESTNAME)

###############################################################################
# cmake script TRY_COMPILE all stxxl header files

if(TRY_COMPILE_HEADERS)

  include(CheckCXXSourceCompiles)
  set(CMAKE_REQUIRED_INCLUDES ${STXXL_INCLUDE_DIRS} ${Boost_INCLUDE_DIRS})
  set(CMAKE_REQUIRED_LIBRARIES ${Boost_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT})

  file(GLOB_RECURSE header_files FOLLOW_SYMLINKS "include/*")
  list(REMOVE_ITEM header_files "${PROJECT_SOURCE_DIR}/include/stxxl/bits/config.h.in")
  list(SORT header_files)

  foreach(file ${header_files})
    string(REPLACE "/" "_" compilename "${file}") # replace / to _ to fix warnings
    string(REPLACE "." "_" compilename "${compilename}")

    check_cxx_source_compiles(
      "#include \"${file}\"
      int main() { return 0; }" IsSelfContained${compilename})

    if(NOT IsSelfContained${compilename})
      message(SEND_ERROR "Compilation FAILED for ${file}\n\nCompiler output:\n${OUTPUT}")
    endif()
  endforeach()

endif(TRY_COMPILE_HEADERS)

###############################################################################
# figure out STXXL_TMPDIR for tests

if(BUILD_TESTS)

  if(DEFINED ENV{STXXL_TMPDIR})
    set(STXXL_TMPDIR "$ENV{STXXL_TMPDIR}")
  elseif(NOT DEFINED STXXL_TMPDIR)
    set(STXXL_TMPDIR ".")
  endif()

  message(STATUS "Using STXXL_TMPDIR ${STXXL_TMPDIR} for tests")

endif(BUILD_TESTS)

###############################################################################
# install header files

install(FILES include/stxxl.h DESTINATION ${INSTALL_INCLUDE_DIR})
install(DIRECTORY include/stxxl DESTINATION ${INSTALL_INCLUDE_DIR})

# also copy the config file with build options!
install(FILES ${PROJECT_BINARY_DIR}/include/stxxl/bits/config.h
  DESTINATION ${INSTALL_INCLUDE_DIR}/stxxl/bits/)

###############################################################################
# prepare pkg-config file

configure_file(misc/cmake/stxxl.pc
  "${PROJECT_BINARY_DIR}/${STXXL_LIBNAME}.pc" @ONLY)

# copy the stxxl.pc file into lib/pkgconfig
if(INSTALL_PKGCONFIG_DIR)
  install(FILES ${PROJECT_BINARY_DIR}/${STXXL_LIBNAME}.pc
    DESTINATION ${INSTALL_PKGCONFIG_DIR})
endif()

###############################################################################
# export targets to cmake project config file

if(${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION} GREATER 2.6)
  # register package for use from the global CMake-registry
  export(PACKAGE stxxl)
endif()

# add stxxl library targets to the build tree export set
export(TARGETS ${STXXL_EXPORTED_LIBS} FILE "${PROJECT_BINARY_DIR}/stxxl-targets.cmake")

# calculate absolute install paths
foreach(dir LIB BIN INCLUDE CMAKE)
  set(ABS_INSTALL_${dir}_DIR "${INSTALL_${dir}_DIR}")
  if(NOT IS_ABSOLUTE "${ABS_INSTALL_${dir}_DIR}")
    set(ABS_INSTALL_${dir}_DIR "${CMAKE_INSTALL_PREFIX}/${ABS_INSTALL_${dir}_DIR}")
  endif()
endforeach()

# calculate final installed include path relative to cmake install path
file(RELATIVE_PATH REL_INCLUDE_DIR "${ABS_INSTALL_CMAKE_DIR}" "${ABS_INSTALL_INCLUDE_DIR}")

# create common stxxl-config.cmake file
configure_file(misc/cmake/stxxl-version.cmake.in
  "${PROJECT_BINARY_DIR}/stxxl-version.cmake" @ONLY)

# create stxxl-version.cmake file for the build tree
set(CONF_INCLUDE_DIRS "${STXXL_INCLUDE_DIRS}")
configure_file(misc/cmake/stxxl-config.cmake.in
  "${PROJECT_BINARY_DIR}/stxxl-config.cmake" @ONLY)

# create stxxl-version.cmake file for the install tree
set(CONF_INCLUDE_DIRS "\${STXXL_CMAKE_DIR}/${REL_INCLUDE_DIR}")
configure_file(misc/cmake/stxxl-config.cmake.in
  "${PROJECT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/stxxl-config.cmake" @ONLY)

# install the stxxl-config.cmake and stxxl-version.cmake
install(FILES
  "${PROJECT_BINARY_DIR}/stxxl-version.cmake"
  "${PROJECT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}/stxxl-config.cmake"
  DESTINATION "${INSTALL_CMAKE_DIR}")

# Install the export set for use with the install-tree
install(EXPORT stxxl-targets DESTINATION "${INSTALL_CMAKE_DIR}")

###############################################################################
# build tests and tools

add_subdirectory(tools)
add_subdirectory(examples)
add_subdirectory(tests)
add_subdirectory(local)
