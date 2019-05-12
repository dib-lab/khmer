/***************************************************************************
 *  tests/containers/test_vector_resize.cpp
 *
 *  Part of the STXXL. See http://stxxl.sourceforge.net
 *
 *  Copyright (C) 2015 Timo Bingmann <tb@panthema.net>
 *
 *  Distributed under the Boost Software License, Version 1.0.
 *  (See accompanying file LICENSE_1_0.txt or copy at
 *  http://www.boost.org/LICENSE_1_0.txt)
 **************************************************************************/

#include <stxxl/io>
#include <stxxl/vector>

typedef stxxl::VECTOR_GENERATOR<int, 4, 4>::result vector_type;

int main()
{
    stxxl::config* config = stxxl::config::get_instance();

    stxxl::disk_config disk1("/tmp/stxxl-###.tmp", 16 * 1024 * 1024,
                             "syscall autogrow=no");
    disk1.unlink_on_open = true;
    disk1.direct = stxxl::disk_config::DIRECT_OFF;
    config->add_disk(disk1);

    vector_type* vector = new vector_type();

    try {
        while (true)
            vector->push_back(0);
    }
    catch (std::exception& e) {
        STXXL_ERRMSG("Caught exception: " << e.what());
        delete vector; // here it will crash in block_manager::delete_block(bid)
    }

    STXXL_MSG("Delete done, all is well.");

    return 0;
}
