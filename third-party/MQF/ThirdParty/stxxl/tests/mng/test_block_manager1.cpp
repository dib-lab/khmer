/***************************************************************************
 *  tests/mng/test_block_manager1.cpp
 *
 *  Part of the STXXL. See http://stxxl.sourceforge.net
 *
 *  Copyright (C) 2002 Roman Dementiev <dementiev@mpi-sb.mpg.de>
 *
 *  Distributed under the Boost Software License, Version 1.0.
 *  (See accompanying file LICENSE_1_0.txt or copy at
 *  http://www.boost.org/LICENSE_1_0.txt)
 **************************************************************************/

#include <stxxl/request>
#include <stxxl/mng>

int main()
{
    typedef stxxl::typed_block<128* 1024, double> block_type;
    std::vector<block_type::bid_type> bids;
    std::vector<stxxl::request_ptr> requests;
    stxxl::block_manager* bm = stxxl::block_manager::get_instance();
    bm->new_blocks<block_type>(32, stxxl::striping(), std::back_inserter(bids));
    std::vector<block_type, stxxl::new_alloc<block_type> > blocks(32);
    int vIndex;
    for (vIndex = 0; vIndex < 32; ++vIndex) {
        for (int vIndex2 = 0; vIndex2 < block_type::size; ++vIndex2) {
            blocks[vIndex][vIndex2] = vIndex2;
        }
    }
    for (vIndex = 0; vIndex < 32; ++vIndex) {
        requests.push_back(blocks[vIndex].write(bids[vIndex]));
    }
    stxxl::wait_all(requests.begin(), requests.end());
    bm->delete_blocks(bids.begin(), bids.end());
    return 0;
}
