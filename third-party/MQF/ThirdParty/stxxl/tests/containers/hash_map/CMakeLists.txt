############################################################################
#  tests/containers/hash_map/CMakeLists.txt
#
#  Part of the STXXL. See http://stxxl.sourceforge.net
#
#  Copyright (C) 2014 Timo Bingmann <tb@panthema.net>
#
#  Distributed under the Boost Software License, Version 1.0.
#  (See accompanying file LICENSE_1_0.txt or copy at
#  http://www.boost.org/LICENSE_1_0.txt)
############################################################################

stxxl_build_test(test_hash_map)
stxxl_build_test(test_hash_map_block_cache)
stxxl_build_test(test_hash_map_iterators)
stxxl_build_test(test_hash_map_reader_writer)

stxxl_test(test_hash_map)
stxxl_test(test_hash_map_block_cache)
stxxl_test(test_hash_map_iterators)
stxxl_test(test_hash_map_reader_writer)
