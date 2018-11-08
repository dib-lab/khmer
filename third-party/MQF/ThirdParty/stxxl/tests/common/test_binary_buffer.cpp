/***************************************************************************
 *  tests/common/test_binary_buffer.cpp
 *
 *  Part of the STXXL. See http://stxxl.sourceforge.net
 *
 *  Copyright (C) 2014 Timo Bingmann <tb@panthema.net>
 *
 *  Distributed under the Boost Software License, Version 1.0.
 *  (See accompanying file LICENSE_1_0.txt or copy at
 *  http://www.boost.org/LICENSE_1_0.txt)
 **************************************************************************/

#include <stxxl/bits/common/binary_buffer.h>
#include <stxxl/bits/verbose.h>

void test1()
{
//! [serialize]
    // construct a binary blob
    stxxl::binary_buffer bb;
    {
        bb.put<unsigned int>(1);
        bb.put_string("test");

        bb.put_varint(42);
        bb.put_varint(12345678);

        // add a sub block
        stxxl::binary_buffer sub;
        sub.put_string("sub block");
        sub.put_varint(6 * 9);

        bb.put_string(sub);
    }
//! [serialize]

    // read binary block and verify content

    stxxl::binary_buffer_ref bbr = bb;

    const unsigned char bb_data[] = {
        // bb.put<unsigned int>(1)
        0x01, 0x00, 0x00, 0x00,
        // bb.put_string("test")
        0x04, 0x74, 0x65, 0x73, 0x74,
        // bb.put_varint(42);
        0x2a,
        // bb.put_varint(12345678);
        0xce, 0xc2, 0xf1, 0x05,
        // begin sub block (length)
        0x0b,
        // sub.put_string("sub block");
        0x09, 0x73, 0x75, 0x62, 0x20, 0x62, 0x6c, 0x6f, 0x63, 0x6b,
        // sub.put_varint(6 * 9);
        0x36,
    };

    stxxl::binary_buffer_ref bb_verify(bb_data, sizeof(bb_data));

    if (bbr != bb_verify)
        std::cout << bbr.str();

    STXXL_CHECK(bbr == bb_verify);

//! [deserialize]
    // read binary block using binary_reader

    stxxl::binary_reader br(bb);

    STXXL_CHECK(br.get<unsigned int>() == 1);
    STXXL_CHECK(br.get_string() == "test");
    STXXL_CHECK(br.get_varint() == 42);
    STXXL_CHECK(br.get_varint() == 12345678);

    {
        stxxl::binary_reader sub_br = br.get_binary_buffer_ref();
        STXXL_CHECK(sub_br.get_string() == "sub block");
        STXXL_CHECK(sub_br.get_varint() == 6 * 9);
        STXXL_CHECK(sub_br.empty());
    }

    STXXL_CHECK(br.empty());
//! [deserialize]
}

int main(int, char**)
{
    test1();
    return 0;
}
