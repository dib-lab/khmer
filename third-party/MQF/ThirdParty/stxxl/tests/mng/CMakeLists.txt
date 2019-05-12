############################################################################
#  tests/mng/CMakeLists.txt
#
#  Part of the STXXL. See http://stxxl.sourceforge.net
#
#  Copyright (C) 2013 Timo Bingmann <tb@panthema.net>
#
#  Distributed under the Boost Software License, Version 1.0.
#  (See accompanying file LICENSE_1_0.txt or copy at
#  http://www.boost.org/LICENSE_1_0.txt)
############################################################################

stxxl_build_test(test_aligned)
stxxl_build_test(test_block_alloc_strategy)
stxxl_build_test(test_block_manager)
stxxl_build_test(test_block_manager1)
stxxl_build_test(test_block_manager2)
stxxl_build_test(test_block_scheduler)
stxxl_build_test(test_bmlayer)
stxxl_build_test(test_buf_streams)
stxxl_build_test(test_config)
stxxl_build_test(test_pool_pair)
stxxl_build_test(test_prefetch_pool)
stxxl_build_test(test_read_write_pool)
stxxl_build_test(test_write_pool)

stxxl_test(test_aligned)
stxxl_test(test_block_alloc_strategy)
stxxl_test(test_block_manager)
stxxl_test(test_block_manager1)
stxxl_test(test_block_manager2)
stxxl_test(test_block_scheduler)
stxxl_test(test_bmlayer)
stxxl_test(test_buf_streams)
stxxl_test(test_config)
stxxl_test(test_pool_pair)
stxxl_test(test_prefetch_pool)
stxxl_test(test_read_write_pool)
stxxl_test(test_write_pool)

add_define(test_block_manager "STXXL_VERBOSE_LEVEL=2")
add_define(test_block_manager1 "STXXL_VERBOSE_LEVEL=2")
add_define(test_block_manager2 "STXXL_VERBOSE_LEVEL=3")
add_define(test_prefetch_pool "STXXL_VERBOSE_LEVEL=2")
add_define(test_write_pool "STXXL_VERBOSE_LEVEL=1")
