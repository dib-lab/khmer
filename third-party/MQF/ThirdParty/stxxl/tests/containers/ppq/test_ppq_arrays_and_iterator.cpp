/***************************************************************************
 *  tests/containers/ppq/test_ppq_arrays_and_iterator.cpp
 *
 *  Part of the STXXL. See http://stxxl.sourceforge.net
 *
 *  Copyright (C) 2014 Thomas Keh <thomas.keh@student.kit.edu>
 *
 *  Distributed under the Boost Software License, Version 1.0.
 *  (See accompanying file LICENSE_1_0.txt or copy at
 *  http://www.boost.org/LICENSE_1_0.txt)
 **************************************************************************/

//! \example containers/external_array.cpp
//! This is an example of how to use \c stxxl::external_array

#include <cstddef>
#include <stxxl/bits/common/cmdline.h>
#include <stxxl/bits/containers/parallel_priority_queue.h>
#include <stxxl/bits/mng/block_manager.h>
#include <stxxl/bits/common/utils.h>
#include <stxxl/bits/verbose.h>
#include <stxxl/timer>

typedef std::pair<size_t, size_t> value_type;
typedef stxxl::ppq_local::external_array<value_type> ea_type;
typedef stxxl::ppq_local::internal_array<value_type> ia_type;

//! Comparator for value_type / size_t-pairs
struct value_comp {
    bool operator () (value_type const& a, value_type const& b)
    {
        return a.first < b.first;
    }
};

//! Frequency of printing the progress
static const size_t printmod = 16 * 1024 * 1024;

//! Prints progress
static inline void progress(const char* text, size_t i, size_t nelements)
{
    if ((i % printmod) == 0) {
        STXXL_MSG(text << " " << i << " (" << std::setprecision(5) << (static_cast<double>(i) * 100. / static_cast<double>(nelements)) << " %)");
    }
}

//! Fills an external array with increasing numbers
void fill(ea_type& ea, size_t index, size_t n)
{
    STXXL_CHECK(n > 0);

    ea_type::writer_type ea_writer(ea);

    size_t i = index;

    STXXL_CHECK_EQUAL(ea_writer.end() - ea_writer.begin(), (ssize_t)n);

    for (ea_type::writer_type::iterator it = ea_writer.begin();
         it != ea_writer.end(); ++it)
    {
        STXXL_CHECK(i < index + n);
        progress("Inserting element", i, n);
        *it = std::make_pair(i + 1, i + 1);
        ++i;
    }

    STXXL_CHECK(i == index + n);
}

//! Run external array test.
//! Testing...
//!     - iterators
//!     - request_write_buffer()
//!     - random access
//!     - remove()
//!     - get_current_max()
//!     - request_further_block()
//!     - wait()
//!     - sizes
//!
//! \param volume Volume
//! \param numwbs Number of write buffer blocks for each external array.
int run_external_array_test(size_t volume)
{
    const size_t block_size = ea_type::block_size;
    const size_t N = volume / sizeof(value_type);
    const size_t num_blocks = stxxl::div_ceil(N, block_size);

    STXXL_VARDUMP(block_size);
    STXXL_VARDUMP(N);
    STXXL_VARDUMP(num_blocks);

    ea_type::pool_type rw_pool(6 * num_blocks, 4);
    ea_type ea(N, &rw_pool);

    {
        stxxl::scoped_print_timer timer("filling", volume);

        STXXL_CHECK(N >= 10 + 355 + 2 * block_size);

        fill(ea, 0, N);
    }

    {
        stxxl::scoped_print_timer timer("reading and checking", volume);

        STXXL_CHECK(N > 5 * block_size + 876);
        STXXL_CHECK(block_size > 34);

        STXXL_CHECK_EQUAL(ea.size(), N);

        // fetch first block
        ea.hint_next_block();
        ea.wait_next_blocks();

        STXXL_CHECK(ea.buffer_size() > 7);

        // Testing iterator...

        for (unsigned i = 0; i < 7; ++i) {
            STXXL_CHECK_EQUAL((ea.begin() + i)->first, i + 1);
        }

        STXXL_CHECK(ea.buffer_size() > 12);

        // Testing random access...

        for (unsigned i = 7; i < 12; ++i) {
            STXXL_CHECK_EQUAL(ea[i].first, i + 1);
        }

        // Testing remove...

        ea.remove_items(33);

        STXXL_CHECK_EQUAL(ea[33].first, 34);
        STXXL_CHECK_EQUAL(ea.begin().get_index(), 33);
        STXXL_CHECK_EQUAL(ea.begin()->first, 34);

        // Testing get_next_block_min()...

        unsigned maxround = 0;

        while (ea.get_next_block_min().first < 5 * block_size + 876) {
            switch (maxround) {
            case 0: STXXL_CHECK_EQUAL(ea.get_next_block_min().first, 131073);
                break;
            case 1: STXXL_CHECK_EQUAL(ea.get_next_block_min().first, 262145);
                break;
            case 2: STXXL_CHECK_EQUAL(ea.get_next_block_min().first, 393217);
                break;
            case 3: STXXL_CHECK_EQUAL(ea.get_next_block_min().first, 524289);
                break;
            case 4: STXXL_CHECK_EQUAL(ea.get_next_block_min().first, 655361);
                break;
            }

            ea.hint_next_block();
            ea.wait_next_blocks();
            ++maxround;
        }

        // Testing request_further_block() (called above)...

        STXXL_CHECK((ea.end() - 1)->first >= 5 * block_size + 876);

        size_t index = 33;

        for (ea_type::iterator it = ea.begin(); it != ea.end(); ++it) {
            progress("Extracting element", index, N);
            STXXL_CHECK_EQUAL(it.get_index(), index);
            STXXL_CHECK_EQUAL(it->first, index + 1);
            ++index;
        }

        ea.remove_items(ea.buffer_size());
        STXXL_CHECK(ea.buffer_size() == 0);

        // Extracting the rest...

        size_t round = 0;

        while (!ea.empty()) {
            if (ea.has_em_data() && round % 2 == 0) {
                ea.hint_next_block();
                ea.wait_next_blocks();
            }

            STXXL_CHECK((size_t)(ea.end() - ea.begin()) == ea.buffer_size());

            for (ea_type::iterator it = ea.begin(); it != ea.end(); ++it) {
                progress("Extracting element", index, N);
                STXXL_CHECK_EQUAL(it.get_index(), index);
                STXXL_CHECK_EQUAL(it->first, index + 1);
                ++index;
            }

            //std::cout << ea.buffer_size() << " - " << round << "\n";

            if (round % 3 == 0 && ea.buffer_size() > 0) {
                // remove all items but one
                ea.remove_items(ea.buffer_size() - 1);
                index--;
            }
            else {
                ea.remove_items(ea.buffer_size());
            }

            ++round;
        }
    }

    return EXIT_SUCCESS;
}

//! Run test for multiway_merge compatibility.
//! Fills NumEAs external arrays and one internal array.
//! Merges them into one external array and checks the content.
//!
//! \param volume Volume
//! \param numwbs Number of write buffer blocks for each external array.
template <int NumEAs>
int run_multiway_merge(size_t volume)
{
    const size_t block_size = ea_type::block_size;
    const size_t size = volume / sizeof(value_type);
    const size_t num_blocks = stxxl::div_ceil(size, block_size);

    STXXL_MSG("--- Running run_multiway_merge() test with " <<
              NumEAs << " external arrays");

    STXXL_VARDUMP(block_size);
    STXXL_VARDUMP(size);
    STXXL_VARDUMP(num_blocks);

    std::vector<value_type> v(size);

    {
        stxxl::scoped_print_timer timer("filling vector / internal array for multiway_merge test", volume);
        for (size_t i = 0; i < size; ++i) {
            v[i] = std::make_pair(i + 1, i + 1);
        }
    }

    ia_type ia(v);

    typedef std::vector<ea_type*> ea_vector_type;
    typedef ea_vector_type::iterator ea_iterator;

    ea_type::pool_type rw_pool1;
    ea_vector_type ealist;
    for (int i = 0; i < NumEAs; ++i)
        ealist.push_back(new ea_type(size, &rw_pool1));

    ea_type::pool_type rw_pool2;
    ea_type out(size * (NumEAs + 1), &rw_pool2);

    {
        stxxl::scoped_print_timer timer("filling external arrays for multiway_merge test", volume * NumEAs);

        for (ea_iterator ea = ealist.begin(); ea != ealist.end(); ++ea)
            fill(*(*ea), 0, size);
    }

    {
        stxxl::scoped_print_timer timer("loading input arrays into RAM and requesting output buffer", volume * NumEAs);

        if (ealist.size())
            rw_pool1.resize_prefetch(ealist.size() * ealist[0]->num_blocks());

        for (ea_iterator ea = ealist.begin(); ea != ealist.end(); ++ea)
        {
            while ((*ea)->has_unhinted_em_data())
                (*ea)->hint_next_block();

            while ((*ea)->num_hinted_blocks())
                (*ea)->wait_next_blocks();
        }

        for (ea_iterator ea = ealist.begin(); ea != ealist.end(); ++ea)
            STXXL_CHECK_EQUAL((*ea)->buffer_size(), size);
    }

    if (NumEAs > 0) // test ea ppq_iterators
    {
        stxxl::scoped_print_timer timer("test ea ppq_iterators", volume);

        ea_type::iterator begin = ealist[0]->begin(), end = ealist[0]->end();
        size_t index = 1;

        // prefix operator ++
        for (ea_type::iterator it = begin; it != end; ++it, ++index)
            STXXL_CHECK_EQUAL(it->first, index);

        // prefix operator --
        for (ea_type::iterator it = end; it != begin; )
        {
            --it, --index;
            STXXL_CHECK_EQUAL(it->first, index);
        }

        STXXL_CHECK_EQUAL(index, 1);

        // addition operator +
        for (ea_type::iterator it = begin; it != end; it = it + 1, ++index)
            STXXL_CHECK_EQUAL(it->first, index);

        // subtraction operator -
        for (ea_type::iterator it = end; it != begin; )
        {
            it = it - 1, --index;
            STXXL_CHECK_EQUAL(it->first, index);
        }
    }

    {
        stxxl::scoped_print_timer timer("running multiway merge", volume * (NumEAs + 1));

        std::vector<std::pair<ea_type::iterator, ea_type::iterator> > seqs;
        for (ea_iterator ea = ealist.begin(); ea != ealist.end(); ++ea)
            seqs.push_back(std::make_pair((*ea)->begin(), (*ea)->end()));

        seqs.push_back(std::make_pair(ia.begin(), ia.end()));

        ea_type::writer_type ea_writer(out);

        stxxl::potentially_parallel::multiway_merge(
            seqs.begin(), seqs.end(),
            ea_writer.begin(), (NumEAs + 1) * size, value_comp());

        // sequential:
        //__gnu_parallel::multiway_merge(seqs.begin(), seqs.end(), c.begin(),
        //    2*size, value_comp()); // , __gnu_parallel::sequential_tag()
    }

    {
        stxxl::scoped_print_timer timer("checking the order", volume * NumEAs);

        // each index is contained (NumEAs + 1) times
        size_t index = 1 * (NumEAs + 1);

        while (out.has_em_data())
        {
            out.hint_next_block();
            out.wait_next_blocks();

            for (ea_type::iterator it = out.begin(); it != out.end(); ++it)
            {
                progress("Extracting element", index, (NumEAs + 1) * size);
                STXXL_CHECK_EQUAL(it->first, index / (NumEAs + 1));

                ++index;
            }

            out.remove_items(out.buffer_size());
        }
    }

    return EXIT_SUCCESS;
}

//! Fill internal array, remove 3 elements, read and check result.
//! \param volume Volume
int run_internal_array_test(size_t volume)
{
    const size_t size = volume / sizeof(value_type);
    STXXL_VARDUMP(size);
    STXXL_CHECK(size > 3);

    std::vector<value_type> v(size);

    {
        stxxl::scoped_print_timer timer("filling vector / internal array", volume);
        for (size_t i = 0; i < size; ++i) {
            v[i] = std::make_pair(i + 1, i + 1);
        }
    }

    ia_type ia(v);

    {
        stxxl::scoped_print_timer timer("reading and checking the order", volume);

        ia.inc_min();
        ia.inc_min(2);

        typedef ia_type::iterator iterator;
        size_t index = 4;

        // test iterator
        for (iterator it = ia.begin(); it != ia.end(); ++it) {
            STXXL_CHECK_EQUAL(it->first, index++);
        }

        for (size_t i = 3; i < size; ++i) {
            STXXL_CHECK_EQUAL(ia.begin()[i - 3].first, i + 1);
            STXXL_CHECK_EQUAL(ia[i].first, i + 1);
        }

        STXXL_CHECK_EQUAL(index, size + 1);
        STXXL_CHECK_EQUAL(ia.size(), size - 3);
    }

    return EXIT_SUCCESS;
}

//! Testing upper_bound compatibility
//! \param volume Volume
int run_upper_bound_test(size_t volume)
{
    const size_t size = volume / sizeof(value_type);
    STXXL_CHECK(volume > ea_type::block_size);
    STXXL_VARDUMP(size);
    STXXL_CHECK(size > 2000);

    ea_type::pool_type rw_pool;

    std::vector<value_type> v(size);

    {
        stxxl::scoped_print_timer timer("filling vector / internal array", volume);
        for (size_t i = 0; i < size; ++i) {
            v[i] = std::make_pair(i + 1, i + 1);
        }
    }

    ia_type a(v);
    ea_type b(size - 10, &rw_pool);

    {
        stxxl::scoped_print_timer timer("filling external array", volume);

        ea_type::writer_type ea_writer(b, 1);

        ea_type::writer_type::iterator it = ea_writer.begin();

        for (size_t i = 0; i < size - 10; ++i, ++it) {
            *it = std::make_pair(2 * i, 2 * i);
        }

        it = ea_writer.begin();
        for (size_t i = 0; i < size - 10; ++i, ++it) {
            STXXL_CHECK_EQUAL(it->first, 2 * i);
        }
    }

    {
        stxxl::scoped_print_timer timer("running upper_bound", 2 * volume);

        b.hint_next_block();
        b.wait_next_blocks();

        value_type minmax = std::make_pair(2000, 2000);

        ea_type::iterator uba = std::upper_bound(a.begin(), a.end(), minmax, value_comp());
        ea_type::iterator ubb = std::upper_bound(b.begin(), b.end(), minmax, value_comp());

        STXXL_CHECK_EQUAL((uba - 1)->first, 2000);
        STXXL_CHECK_EQUAL((ubb - 1)->first, 2000);
        STXXL_CHECK_EQUAL(uba->first, 2001);
        STXXL_CHECK_EQUAL(ubb->first, 2002);
        STXXL_CHECK_EQUAL(uba.get_index(), 2000);
        STXXL_CHECK_EQUAL(ubb.get_index(), 1001);
    }

    return EXIT_SUCCESS;
}

int main(int argc, char** argv)
{
    stxxl::uint64 volume = 512 * 1024 * 1024;
    stxxl::uint64 mwmvolume = 50 * 1024 * 1024;
    stxxl::uint64 iavolume = 10 * 1024 * 1024;

    unsigned numpbs = 1;
    unsigned numwbs = 14;

    stxxl::cmdline_parser cp;
    cp.set_description("STXXL external array test");
    cp.set_author("Thomas Keh <thomas.keh@student.kit.edu>");
    cp.add_bytes('v', "volume", volume,
                 "Volume to fill into the external array, "
                 "default: 512 MiB, 0 = disable test");
    cp.add_bytes('m', "mwmvolume", mwmvolume,
                 "Testing multiway merge of two external arrays - "
                 "the volume of each input array, "
                 "default: 100 MiB, 0 = disable test");
    cp.add_bytes('i', "iavolume", iavolume,
                 "Volume to fill into the internal array, "
                 "default: 10 MiB, 0 = disable test");
    cp.add_uint('n', "num_prefetch_buffers", numpbs,
                "Number of prefetch buffer blocks, default: 1");
    cp.add_uint('w', "num_write_buffers", numwbs,
                "Number of write buffer blocks, default: 14");

    if (!cp.process(argc, argv))
        return EXIT_FAILURE;

    const size_t block_size = 2 * 1024 * 1024 / sizeof(value_type);
    if (volume > 0 && volume / sizeof(value_type) < 5 * block_size + 876) {
        STXXL_ERRMSG("The volume is too small for this test. It must be >= " <<
                     (5 * block_size + 876) * sizeof(value_type));
        return EXIT_FAILURE;
    }

    if (iavolume > 0 && iavolume / sizeof(value_type) < 3) {
        STXXL_ERRMSG("The internal array volume is too small for this test. "
                     "It must be > " << 3);
        return EXIT_FAILURE;
    }

    STXXL_MEMDUMP(volume);
    STXXL_MEMDUMP(mwmvolume);
    STXXL_MEMDUMP(iavolume);
    STXXL_MEMDUMP(sizeof(value_type));
    STXXL_VARDUMP(numpbs);
    STXXL_VARDUMP(numwbs);

    bool succ = EXIT_SUCCESS;

    if (volume > 0) {
        succ = run_external_array_test(volume) && succ;
    }

    if (iavolume > 0) {
        succ = run_internal_array_test(iavolume) && succ;
    }

    if (mwmvolume > 0) {
        succ = run_multiway_merge<0>(mwmvolume) && succ;
        succ = run_multiway_merge<1>(mwmvolume) && succ;
        succ = run_multiway_merge<2>(mwmvolume) && succ;
        succ = run_multiway_merge<3>(mwmvolume) && succ;
        succ = run_multiway_merge<4>(mwmvolume) && succ;
        succ = run_multiway_merge<5>(mwmvolume) && succ;
    }

    succ = run_upper_bound_test(3 * ea_type::block_size * sizeof(value_type)) && succ;

    STXXL_MSG("success = " << succ);

    return succ;
}
