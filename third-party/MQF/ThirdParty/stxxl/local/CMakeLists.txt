############################################################################
#  local/CMakeLists.txt
#
#  CMake file which compiles and links all .cpp files in local/ with STXXL.
#
#  Part of the STXXL. See http://stxxl.sourceforge.net
#
#  Copyright (C) 2013 Timo Bingmann <tb@panthema.net>
#
#  Distributed under the Boost Software License, Version 1.0.
#  (See accompanying file LICENSE_1_0.txt or copy at
#  http://www.boost.org/LICENSE_1_0.txt)
############################################################################

# find all .cpp files in current directory
file(GLOB source_files RELATIVE ${CMAKE_CURRENT_SOURCE_DIR} "*.cpp")

# compile and link
foreach(file ${source_files})
  # remove .cpp
  string(REGEX REPLACE "\\.cpp$" "" file ${file})
  # call macro
  stxxl_build_tool(${file})
endforeach()
