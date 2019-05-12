############################################################################
#  tests/containers/btree/CMakeLists.txt
#
#  Part of the STXXL. See http://stxxl.sourceforge.net
#
#  Copyright (C) 2013-2014 Timo Bingmann <tb@panthema.net>
#
#  Distributed under the Boost Software License, Version 1.0.
#  (See accompanying file LICENSE_1_0.txt or copy at
#  http://www.boost.org/LICENSE_1_0.txt)
############################################################################

stxxl_build_test(test_btree)
stxxl_build_test(test_btree_const_scan)
stxxl_build_test(test_btree_insert_erase)
stxxl_build_test(test_btree_insert_find)
stxxl_build_test(test_btree_insert_scan)

stxxl_test(test_btree 10000)
stxxl_test(test_btree 100000)
stxxl_test(test_btree 1000000)
stxxl_test(test_btree_const_scan 10000)
stxxl_test(test_btree_const_scan 100000)
stxxl_test(test_btree_const_scan 1000000)
stxxl_test(test_btree_insert_erase 14)
stxxl_test(test_btree_insert_find 14)
stxxl_test(test_btree_insert_scan 14)
