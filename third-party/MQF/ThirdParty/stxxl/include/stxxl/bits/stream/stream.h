/***************************************************************************
 *  include/stxxl/bits/stream/stream.h
 *
 *  Part of the STXXL. See http://stxxl.sourceforge.net
 *
 *  Copyright (C) 2003-2005 Roman Dementiev <dementiev@mpi-sb.mpg.de>
 *  Copyright (C) 2009, 2010 Andreas Beckmann <beckmann@cs.uni-frankfurt.de>
 *  Copyright (C) 2010 Johannes Singler <singler@kit.edu>
 *
 *  Distributed under the Boost Software License, Version 1.0.
 *  (See accompanying file LICENSE_1_0.txt or copy at
 *  http://www.boost.org/LICENSE_1_0.txt)
 **************************************************************************/

#ifndef STXXL_STREAM_STREAM_HEADER
#define STXXL_STREAM_STREAM_HEADER

#include <stxxl/bits/namespace.h>
#include <stxxl/bits/mng/buf_istream.h>
#include <stxxl/bits/mng/buf_ostream.h>
#include <stxxl/bits/common/tuple.h>
#include <stxxl/bits/common/error_handling.h>
#include <stxxl/vector>
#include <stxxl/bits/compat/unique_ptr.h>

#ifndef STXXL_VERBOSE_MATERIALIZE
#define STXXL_VERBOSE_MATERIALIZE STXXL_VERBOSE3
#endif

STXXL_BEGIN_NAMESPACE

//! Stream package subnamespace.
namespace stream {

//! \addtogroup streampack
//! \{

////////////////////////////////////////////////////////////////////////
//     STREAMIFY                                                      //
////////////////////////////////////////////////////////////////////////

//! A model of stream that retrieves the data from an input iterator.
//! For convenience use \c streamify function instead of direct instantiation
//! of \c iterator2stream .
template <class InputIterator>
class iterator2stream
{
    InputIterator m_current, m_end;

public:
    //! Standard stream typedef.
    typedef typename std::iterator_traits<InputIterator>::value_type value_type;

    iterator2stream(InputIterator begin, InputIterator end)
        : m_current(begin), m_end(end)
    { }

    iterator2stream(const iterator2stream& a)
        : m_current(a.m_current), m_end(a.m_end)
    { }

    //! Standard stream method.
    const value_type& operator * () const
    {
        return *m_current;
    }

    const value_type* operator -> () const
    {
        return &(*m_current);
    }

    //! Standard stream method.
    iterator2stream& operator ++ ()
    {
        assert(m_end != m_current);
        ++m_current;
        return *this;
    }

    //! Standard stream method.
    bool empty() const
    {
        return (m_current == m_end);
    }
};

//! Input iterator range to stream converter.
//! \param begin iterator, pointing to the first value
//! \param end iterator, pointing to the last + 1 position, i.e. beyond the range
//! \return an instance of a stream object
template <class InputIterator>
iterator2stream<InputIterator> streamify(InputIterator begin, InputIterator end)
{
    return iterator2stream<InputIterator>(begin, end);
}

//! Traits class of \c streamify function.
template <class InputIterator>
struct streamify_traits
{
    //! return type (stream type) of \c streamify for \c InputIterator.
    typedef iterator2stream<InputIterator> stream_type;
};

//! A model of stream that retrieves data from an external \c stxxl::vector
//! iterator.  It is more efficient than generic \c iterator2stream thanks to
//! use of overlapping For convenience use \c streamify function instead of
//! direct instantiation of \c vector_iterator2stream .
template <class InputIterator>
class vector_iterator2stream
{
    InputIterator m_current, m_end;
    typedef buf_istream<typename InputIterator::block_type,
                        typename InputIterator::bids_container_iterator> buf_istream_type;

    typedef typename stxxl::compat_unique_ptr<buf_istream_type>::result buf_istream_unique_ptr_type;
    mutable buf_istream_unique_ptr_type in;

    void delete_stream()
    {
        in.reset();      // delete object
    }

public:
    typedef vector_iterator2stream<InputIterator> self_type;

    //! Standard stream typedef.
    typedef typename std::iterator_traits<InputIterator>::value_type value_type;

    vector_iterator2stream(InputIterator begin, InputIterator end,
                           unsigned_type nbuffers = 0)
        : m_current(begin), m_end(end),
          in(static_cast<buf_istream_type*>(NULL))
    {
        if (empty())
            return;

        begin.flush();         // flush container
        typename InputIterator::bids_container_iterator end_iter
            = end.bid() + ((end.block_offset()) ? 1 : 0);

        if (end_iter - begin.bid() > 0)
        {
            in.reset(new buf_istream_type(
                         begin.bid(), end_iter, nbuffers ? nbuffers :
                         (2 * config::get_instance()->disks_number())
                         )
                     );

            InputIterator cur = begin - begin.block_offset();

            // skip the beginning of the block
            for ( ; cur != begin; ++cur)
                ++(*in);
        }
    }

    vector_iterator2stream(const self_type& a)
        : m_current(a.m_current), m_end(a.m_end), in(a.in.release())
    { }

    //! Standard stream method.
    const value_type& operator * () const
    {
        return **in;
    }

    const value_type* operator -> () const
    {
        return &(**in);
    }

    //! Standard stream method.
    self_type& operator ++ ()
    {
        assert(m_end != m_current);
        ++m_current;
        ++(*in);
        if (UNLIKELY(empty()))
            delete_stream();

        return *this;
    }

    //! Standard stream method.
    bool empty() const
    {
        return (m_current == m_end);
    }
    virtual ~vector_iterator2stream()
    {
        delete_stream();          // not needed actually
    }
};

//! Input external \c stxxl::vector iterator range to stream converter.
//! It is more efficient than generic input iterator \c streamify thanks to use of overlapping
//! \param begin iterator, pointing to the first value
//! \param end iterator, pointing to the last + 1 position, i.e. beyond the range
//! \param nbuffers number of blocks used for overlapped reading (0 is default,
//! which equals to (2 * number_of_disks)
//! \return an instance of a stream object

template <typename ValueType, typename AllocStr, typename SizeType,
          typename DiffType, unsigned BlockSize, typename PagerType,
          unsigned PageSize>
vector_iterator2stream<
    stxxl::vector_iterator<ValueType, AllocStr, SizeType,
                           DiffType, BlockSize, PagerType, PageSize>
    >
streamify(
    stxxl::vector_iterator<ValueType, AllocStr, SizeType,
                           DiffType, BlockSize, PagerType, PageSize> begin,
    stxxl::vector_iterator<ValueType, AllocStr, SizeType,
                           DiffType, BlockSize, PagerType, PageSize> end,
    unsigned_type nbuffers = 0)
{
    STXXL_VERBOSE1("streamify for vector_iterator range is called");
    return vector_iterator2stream<
        stxxl::vector_iterator<ValueType, AllocStr, SizeType,
                               DiffType, BlockSize, PagerType, PageSize>
        >(begin, end, nbuffers);
}

template <typename ValueType, typename AllocStr, typename SizeType,
          typename DiffType, unsigned BlockSize, typename PagerType,
          unsigned PageSize>
struct streamify_traits<
    stxxl::vector_iterator<ValueType, AllocStr, SizeType,
                           DiffType, BlockSize, PagerType, PageSize>
    >
{
    typedef vector_iterator2stream<
            stxxl::vector_iterator<ValueType, AllocStr, SizeType,
                                   DiffType, BlockSize, PagerType, PageSize>
            > stream_type;
};

//! Input external \c stxxl::vector const iterator range to stream converter.
//! It is more efficient than generic input iterator \c streamify thanks to use of overlapping
//! \param begin const iterator, pointing to the first value
//! \param end const iterator, pointing to the last + 1 position, i.e. beyond the range
//! \param nbuffers number of blocks used for overlapped reading (0 is default,
//! which equals to (2 * number_of_disks)
//! \return an instance of a stream object

template <typename ValueType, typename AllocStr, typename SizeType,
          typename DiffType, unsigned BlockSize, typename PagerType,
          unsigned PageSize>
vector_iterator2stream<
    stxxl::const_vector_iterator<ValueType, AllocStr, SizeType,
                                 DiffType, BlockSize, PagerType, PageSize>
    >
streamify(
    stxxl::const_vector_iterator<ValueType, AllocStr, SizeType,
                                 DiffType, BlockSize, PagerType, PageSize> begin,
    stxxl::const_vector_iterator<ValueType, AllocStr, SizeType,
                                 DiffType, BlockSize, PagerType, PageSize> end,
    unsigned_type nbuffers = 0)
{
    STXXL_VERBOSE1("streamify for const_vector_iterator range is called");
    return vector_iterator2stream<
        stxxl::const_vector_iterator<ValueType, AllocStr, SizeType,
                                     DiffType, BlockSize, PagerType, PageSize>
        >(begin, end, nbuffers);
}

template <typename ValueType, typename AllocStr, typename SizeType,
          typename DiffType, unsigned BlockSize, typename PagerType,
          unsigned PageSize>
struct streamify_traits<
    stxxl::const_vector_iterator<ValueType, AllocStr, SizeType,
                                 DiffType, BlockSize, PagerType, PageSize>
    >
{
    typedef vector_iterator2stream<
            stxxl::const_vector_iterator<ValueType, AllocStr, SizeType,
                                         DiffType, BlockSize, PagerType, PageSize>
            > stream_type;
};

//! Version of  \c iterator2stream. Switches between \c vector_iterator2stream and \c iterator2stream .
//!
//! small range switches between
//! \c vector_iterator2stream and \c iterator2stream .
//! iterator2stream is chosen if the input iterator range
//! is small ( < B )
template <class InputIterator>
class vector_iterator2stream_sr
{
    vector_iterator2stream<InputIterator>* vec_it_stream;
    iterator2stream<InputIterator>* it_stream;

    typedef typename InputIterator::block_type block_type;

public:
    typedef vector_iterator2stream_sr<InputIterator> self_type;

    //! Standard stream typedef.
    typedef typename std::iterator_traits<InputIterator>::value_type value_type;

    vector_iterator2stream_sr(InputIterator begin, InputIterator end,
                              unsigned_type nbuffers = 0)
    {
        if (end - begin < block_type::size)
        {
            STXXL_VERBOSE1("vector_iterator2stream_sr::vector_iterator2stream_sr: Choosing iterator2stream<InputIterator>");
            it_stream = new iterator2stream<InputIterator>(begin, end);
            vec_it_stream = NULL;
        }
        else
        {
            STXXL_VERBOSE1("vector_iterator2stream_sr::vector_iterator2stream_sr: Choosing vector_iterator2stream<InputIterator>");
            it_stream = NULL;
            vec_it_stream = new vector_iterator2stream<InputIterator>(begin, end, nbuffers);
        }
    }

    vector_iterator2stream_sr(const self_type& a)
        : vec_it_stream(a.vec_it_stream), it_stream(a.it_stream)
    { }

    //! Standard stream method.
    const value_type& operator * () const
    {
        if (it_stream)
            return **it_stream;

        return **vec_it_stream;
    }

    const value_type* operator -> () const
    {
        if (it_stream)
            return &(**it_stream);

        return &(**vec_it_stream);
    }

    //! Standard stream method.
    self_type& operator ++ ()
    {
        if (it_stream)
            ++(*it_stream);

        else
            ++(*vec_it_stream);

        return *this;
    }

    //! Standard stream method.
    bool empty() const
    {
        if (it_stream)
            return it_stream->empty();

        return vec_it_stream->empty();
    }
    virtual ~vector_iterator2stream_sr()
    {
        if (it_stream)
            delete it_stream;

        else
            delete vec_it_stream;
    }
};

//! Version of \c streamify. Switches from \c vector_iterator2stream to \c
//! iterator2stream for small ranges.
template <typename ValueType, typename AllocStr, typename SizeType,
          typename DiffType, unsigned BlockSize, typename PagerType,
          unsigned PageSize>
vector_iterator2stream_sr<
    stxxl::vector_iterator<ValueType, AllocStr, SizeType,
                           DiffType, BlockSize, PagerType, PageSize>
    >
streamify_sr(
    stxxl::vector_iterator<ValueType, AllocStr, SizeType,
                           DiffType, BlockSize, PagerType, PageSize> begin,
    stxxl::vector_iterator<ValueType, AllocStr, SizeType,
                           DiffType, BlockSize, PagerType, PageSize> end,
    unsigned_type nbuffers = 0)
{
    STXXL_VERBOSE1("streamify_sr for vector_iterator range is called");
    return vector_iterator2stream_sr<
        stxxl::vector_iterator<ValueType, AllocStr, SizeType,
                               DiffType, BlockSize, PagerType, PageSize>
        >(begin, end, nbuffers);
}

//! Version of \c streamify. Switches from \c vector_iterator2stream to \c
//! iterator2stream for small ranges.
template <typename ValueType, typename AllocStr, typename SizeType,
          typename DiffType, unsigned BlockSize, typename PagerType,
          unsigned PageSize>
vector_iterator2stream_sr<
    stxxl::const_vector_iterator<ValueType, AllocStr, SizeType,
                                 DiffType, BlockSize, PagerType, PageSize>
    >
streamify_sr(
    stxxl::const_vector_iterator<ValueType, AllocStr, SizeType,
                                 DiffType, BlockSize, PagerType, PageSize> begin,
    stxxl::const_vector_iterator<ValueType, AllocStr, SizeType,
                                 DiffType, BlockSize, PagerType, PageSize> end,
    unsigned_type nbuffers = 0)
{
    STXXL_VERBOSE1("streamify_sr for const_vector_iterator range is called");
    return vector_iterator2stream_sr<
        stxxl::const_vector_iterator<ValueType, AllocStr, SizeType,
                                     DiffType, BlockSize, PagerType, PageSize>
        >(begin, end, nbuffers);
}

////////////////////////////////////////////////////////////////////////
//     MATERIALIZE                                                    //
////////////////////////////////////////////////////////////////////////

//! Stores consecutively stream content to an output iterator.
//! \param in stream to be stored used as source
//! \param out output iterator used as destination
//! \return value of the output iterator after all increments,
//! i.e. points to the first unwritten value
//! \pre Output (range) is large enough to hold the all elements in the input stream
template <class OutputIterator, class StreamAlgorithm>
OutputIterator materialize(StreamAlgorithm& in, OutputIterator out)
{
    STXXL_VERBOSE_MATERIALIZE(STXXL_PRETTY_FUNCTION_NAME);
    while (!in.empty())
    {
        *out = *in;
        ++out;
        ++in;
    }
    return out;
}

//! Stores consecutively stream content to an output iterator range \b until end of the stream or end of the iterator range is reached.
//! \param in stream to be stored used as source
//! \param outbegin output iterator used as destination
//! \param outend output end iterator, pointing beyond the output range
//! \return value of the output iterator after all increments,
//! i.e. points to the first unwritten value
//! \pre Output range is large enough to hold the all elements in the input stream
//!
//! This function is useful when you do not know the length of the stream beforehand.
template <class OutputIterator, class StreamAlgorithm>
OutputIterator materialize(StreamAlgorithm& in,
                           OutputIterator outbegin, OutputIterator outend)
{
    STXXL_VERBOSE_MATERIALIZE(STXXL_PRETTY_FUNCTION_NAME);
    while ((!in.empty()) && outend != outbegin)
    {
        *outbegin = *in;
        ++outbegin;
        ++in;
    }
    return outbegin;
}

//! Stores consecutively stream content to an output \c stxxl::vector iterator \b until end of the stream or end of the iterator range is reached.
//! \param in stream to be stored used as source
//! \param outbegin output \c stxxl::vector iterator used as destination
//! \param outend output end iterator, pointing beyond the output range
//! \param nbuffers number of blocks used for overlapped writing (0 is default,
//! which equals to (2 * number_of_disks)
//! \return value of the output iterator after all increments,
//! i.e. points to the first unwritten value
//! \pre Output range is large enough to hold the all elements in the input stream
//!
//! This function is useful when you do not know the length of the stream beforehand.
template <typename ValueType, typename AllocStr, typename SizeType,
          typename DiffType, unsigned BlockSize, typename PagerType,
          unsigned PageSize, class StreamAlgorithm>
stxxl::vector_iterator<ValueType, AllocStr, SizeType,
                       DiffType, BlockSize, PagerType, PageSize>
materialize(StreamAlgorithm& in,
            stxxl::vector_iterator<ValueType, AllocStr, SizeType,
                                   DiffType, BlockSize, PagerType, PageSize> outbegin,
            stxxl::vector_iterator<ValueType, AllocStr, SizeType,
                                   DiffType, BlockSize, PagerType, PageSize> outend,
            unsigned_type nbuffers = 0)
{
    STXXL_VERBOSE_MATERIALIZE(STXXL_PRETTY_FUNCTION_NAME);
    typedef stxxl::vector_iterator<ValueType, AllocStr, SizeType, DiffType, BlockSize, PagerType, PageSize> ExtIterator;
    typedef stxxl::const_vector_iterator<ValueType, AllocStr, SizeType, DiffType, BlockSize, PagerType, PageSize> ConstExtIterator;
    typedef buf_ostream<typename ExtIterator::block_type, typename ExtIterator::bids_container_iterator> buf_ostream_type;

    while (outbegin.block_offset())     //  go to the beginning of the block
    //  of the external vector
    {
        if (in.empty() || outbegin == outend)
            return outbegin;

        *outbegin = *in;
        ++outbegin;
        ++in;
    }

    if (nbuffers == 0)
        nbuffers = 2 * config::get_instance()->disks_number();

    outbegin.flush();     // flush container

    // create buffered write stream for blocks
    buf_ostream_type outstream(outbegin.bid(), nbuffers);

    assert(outbegin.block_offset() == 0);

    // delay calling block_externally_updated() until the block is
    // completely filled (and written out) in outstream
    ConstExtIterator prev_block = outbegin;

    while (!in.empty() && outend != outbegin)
    {
        if (outbegin.block_offset() == 0) {
            if (prev_block != outbegin) {
                prev_block.block_externally_updated();
                prev_block = outbegin;
            }
        }

        *outstream = *in;
        ++outbegin;
        ++outstream;
        ++in;
    }

    ConstExtIterator const_out = outbegin;

    while (const_out.block_offset())     // filling the rest of the block
    {
        *outstream = *const_out;
        ++const_out;
        ++outstream;
    }

    if (prev_block != outbegin)
        prev_block.block_externally_updated();

    outbegin.flush();

    return outbegin;
}

//! Stores consecutively stream content to an output \c stxxl::vector iterator.
//! \param in stream to be stored used as source
//! \param out output \c stxxl::vector iterator used as destination
//! \param nbuffers number of blocks used for overlapped writing (0 is default,
//! which equals to (2 * number_of_disks)
//! \return value of the output iterator after all increments,
//! i.e. points to the first unwritten value
//! \pre Output (range) is large enough to hold the all elements in the input stream
template <typename ValueType, typename AllocStr, typename SizeType,
          typename DiffType, unsigned BlockSize, typename PagerType,
          unsigned PageSize, class StreamAlgorithm>
stxxl::vector_iterator<ValueType, AllocStr, SizeType,
                       DiffType, BlockSize, PagerType, PageSize>
materialize(StreamAlgorithm& in,
            stxxl::vector_iterator<ValueType, AllocStr, SizeType,
                                   DiffType, BlockSize, PagerType, PageSize> out,
            unsigned_type nbuffers = 0)
{
    STXXL_VERBOSE_MATERIALIZE(STXXL_PRETTY_FUNCTION_NAME);
    typedef stxxl::vector_iterator<ValueType, AllocStr, SizeType, DiffType, BlockSize, PagerType, PageSize> ExtIterator;
    typedef stxxl::const_vector_iterator<ValueType, AllocStr, SizeType, DiffType, BlockSize, PagerType, PageSize> ConstExtIterator;
    typedef buf_ostream<typename ExtIterator::block_type, typename ExtIterator::bids_container_iterator> buf_ostream_type;

    // on the I/O complexity of "materialize":
    // crossing block boundary causes O(1) I/Os
    // if you stay in a block, then materialize function accesses only the cache of the
    // vector (only one block indeed), amortized complexity should apply here

    while (out.block_offset())     //  go to the beginning of the block
    //  of the external vector
    {
        if (in.empty())
            return out;

        *out = *in;
        ++out;
        ++in;
    }

    if (nbuffers == 0)
        nbuffers = 2 * config::get_instance()->disks_number();

    out.flush();     // flush container

    // create buffered write stream for blocks
    buf_ostream_type outstream(out.bid(), nbuffers);

    assert(out.block_offset() == 0);

    // delay calling block_externally_updated() until the block is
    // completely filled (and written out) in outstream
    ConstExtIterator prev_block = out;

    while (!in.empty())
    {
        if (out.block_offset() == 0) {
            if (prev_block != out) {
                prev_block.block_externally_updated();
                prev_block = out;
            }
        }

        // tells the vector that the block was modified
        *outstream = *in;
        ++out;
        ++outstream;
        ++in;
    }

    ConstExtIterator const_out = out;

    // copy over items remaining in block from vector.
    while (const_out.block_offset())
    {
        *outstream = *const_out;                 // might cause I/Os for loading the page that
        ++const_out;                             // contains data beyond out
        ++outstream;
    }

    if (prev_block != out)
        prev_block.block_externally_updated();

    out.flush();

    return out;
}

//! Reads stream content and discards it.
//! Useful where you do not need the processed stream anymore,
//! but are just interested in side effects, or just for debugging.
//! \param in input stream
template <class StreamAlgorithm>
void discard(StreamAlgorithm& in)
{
    while (!in.empty())
    {
        *in;
        ++in;
    }
}

////////////////////////////////////////////////////////////////////////
//     GENERATE                                                       //
////////////////////////////////////////////////////////////////////////

//! A model of stream that outputs data from an adaptable generator functor.
//! For convenience use \c streamify function instead of direct instantiation
//! of \c generator2stream .
template <class Generator, typename T = typename Generator::value_type>
class generator2stream
{
public:
    //! Standard stream typedef.
    typedef T value_type;

private:
    Generator gen_;
    value_type m_current;

public:
    generator2stream(Generator g)
        : gen_(g), m_current(gen_())
    { }

    generator2stream(const generator2stream& a) : gen_(a.gen_), m_current(a.m_current) { }

    //! Standard stream method.
    const value_type& operator * () const
    {
        return m_current;
    }

    const value_type* operator -> () const
    {
        return &m_current;
    }

    //! Standard stream method.
    generator2stream& operator ++ ()
    {
        m_current = gen_();
        return *this;
    }

    //! Standard stream method.
    bool empty() const
    {
        return false;
    }
};

//! Adaptable generator to stream converter.
//! \param gen_ generator object
//! \return an instance of a stream object
template <class Generator>
generator2stream<Generator> streamify(Generator gen_)
{
    return generator2stream<Generator>(gen_);
}

////////////////////////////////////////////////////////////////////////
//     TRANSFORM                                                      //
////////////////////////////////////////////////////////////////////////

struct Stopper { };

//! Processes (up to) 6 input streams using given operation functor.
//!
//! \tparam Operation type of the operation (type of an
//! adaptable functor that takes 6 parameters)
//! \tparam Input1 type of the 1st input
//! \tparam Input2 type of the 2nd input
//! \tparam Input3 type of the 3rd input
//! \tparam Input4 type of the 4th input
//! \tparam Input5 type of the 5th input
//! \tparam Input6 type of the 6th input
template <class Operation,
          class Input1,
          class Input2 = Stopper,
          class Input3 = Stopper,
          class Input4 = Stopper,
          class Input5 = Stopper,
          class Input6 = Stopper
          >
class transform
{
    Operation& op;
    Input1& i1;
    Input2& i2;
    Input3& i3;
    Input4& i4;
    Input5& i5;
    Input6& i6;

public:
    //! Standard stream typedef.
    typedef typename Operation::value_type value_type;

private:
    value_type current;

public:
    //! Construction.
    transform(Operation& o, Input1& i1_, Input2& i2_, Input3& i3_, Input4& i4_,
              Input5& i5_, Input5& i6_)
        : op(o), i1(i1_), i2(i2_), i3(i3_), i4(i4_), i5(i5_), i6(i6_)
    {
        if (!empty())
            current = op(*i1, *i2, *i3, *i4, *i5, *i6);
    }

    //! Standard stream method.
    const value_type& operator * () const
    {
        return current;
    }

    const value_type* operator -> () const
    {
        return &current;
    }

    //! Standard stream method.
    transform& operator ++ ()
    {
        ++i1;
        ++i2;
        ++i3;
        ++i4;
        ++i5;
        ++i6;
        if (!empty())
            current = op(*i1, *i2, *i3, *i4, *i5, *i6);

        return *this;
    }

    //! Standard stream method.
    bool empty() const
    {
        return i1.empty() || i2.empty() || i3.empty() ||
               i4.empty() || i5.empty() || i6.empty();
    }
};

// Specializations

////////////////////////////////////////////////////////////////////////
//     TRANSFORM (1 input stream)                                     //
////////////////////////////////////////////////////////////////////////

//! Processes an input stream using given operation functor.
//!
//! \tparam Operation type of the operation (type of an
//! adaptable functor that takes 1 parameter)
//! \tparam Input1 type of the input
//! \remark This is a specialization of \c transform .
template <class Operation,
          class Input1
          >
class transform<Operation, Input1, Stopper, Stopper, Stopper, Stopper, Stopper>
{
    Operation& op;
    Input1& i1;

public:
    //! Standard stream typedef.
    typedef typename Operation::value_type value_type;

private:
    value_type current;

public:
    //! Construction.
    transform(Operation& o, Input1& i1_) : op(o), i1(i1_)
    {
        if (!empty())
            current = op(*i1);
    }

    //! Standard stream method.
    const value_type& operator * () const
    {
        return current;
    }

    const value_type* operator -> () const
    {
        return &current;
    }

    //! Standard stream method.
    transform& operator ++ ()
    {
        ++i1;
        if (!empty())
            current = op(*i1);

        return *this;
    }

    //! Standard stream method.
    bool empty() const
    {
        return i1.empty();
    }
};

////////////////////////////////////////////////////////////////////////
//     TRANSFORM (2 input streams)                                    //
////////////////////////////////////////////////////////////////////////

//! Processes 2 input streams using given operation functor.
//!
//! \tparam Operation type of the operation (type of an
//! adaptable functor that takes 2 parameters)
//! \tparam Input1 type of the 1st input
//! \tparam Input2 type of the 2nd input
//! \remark This is a specialization of \c transform .
template <class Operation,
          class Input1,
          class Input2
          >
class transform<Operation, Input1, Input2, Stopper, Stopper, Stopper, Stopper>
{
    Operation& op;
    Input1& i1;
    Input2& i2;

public:
    //! Standard stream typedef.
    typedef typename Operation::value_type value_type;

private:
    value_type current;

public:
    //! Construction.
    transform(Operation& o, Input1& i1_, Input2& i2_) : op(o), i1(i1_), i2(i2_)
    {
        if (!empty())
            current = op(*i1, *i2);
    }

    //! Standard stream method.
    const value_type& operator * () const
    {
        return current;
    }

    const value_type* operator -> () const
    {
        return &current;
    }

    //! Standard stream method.
    transform& operator ++ ()
    {
        ++i1;
        ++i2;
        if (!empty())
            current = op(*i1, *i2);

        return *this;
    }

    //! Standard stream method.
    bool empty() const
    {
        return i1.empty() || i2.empty();
    }
};

////////////////////////////////////////////////////////////////////////
//     TRANSFORM (3 input streams)                                    //
////////////////////////////////////////////////////////////////////////

//! Processes 3 input streams using given operation functor.
//!
//! \tparam Operation type of the operation (type of an
//! adaptable functor that takes 3 parameters)
//! \tparam Input1 type of the 1st input
//! \tparam Input2 type of the 2nd input
//! \tparam Input3 type of the 3rd input
//! \remark This is a specialization of \c transform .
template <class Operation,
          class Input1,
          class Input2,
          class Input3
          >
class transform<Operation, Input1, Input2, Input3, Stopper, Stopper, Stopper>
{
    Operation& op;
    Input1& i1;
    Input2& i2;
    Input3& i3;

public:
    //! Standard stream typedef.
    typedef typename Operation::value_type value_type;

private:
    value_type current;

public:
    //! Construction.
    transform(Operation& o, Input1& i1_, Input2& i2_, Input3& i3_)
        : op(o), i1(i1_), i2(i2_), i3(i3_)
    {
        if (!empty())
            current = op(*i1, *i2, *i3);
    }

    //! Standard stream method.
    const value_type& operator * () const
    {
        return current;
    }

    const value_type* operator -> () const
    {
        return &current;
    }

    //! Standard stream method.
    transform& operator ++ ()
    {
        ++i1;
        ++i2;
        ++i3;
        if (!empty())
            current = op(*i1, *i2, *i3);

        return *this;
    }

    //! Standard stream method.
    bool empty() const
    {
        return i1.empty() || i2.empty() || i3.empty();
    }
};

////////////////////////////////////////////////////////////////////////
//     TRANSFORM (4 input streams)                                    //
////////////////////////////////////////////////////////////////////////

//! Processes 4 input streams using given operation functor.
//!
//! \tparam Operation type of the operation (type of an
//! adaptable functor that takes 4 parameters)
//! \tparam Input1 type of the 1st input
//! \tparam Input2 type of the 2nd input
//! \tparam Input3 type of the 3rd input
//! \tparam Input4 type of the 4th input
//! \remark This is a specialization of \c transform .
template <class Operation,
          class Input1,
          class Input2,
          class Input3,
          class Input4
          >
class transform<Operation, Input1, Input2, Input3, Input4, Stopper, Stopper>
{
    Operation& op;
    Input1& i1;
    Input2& i2;
    Input3& i3;
    Input4& i4;

public:
    //! Standard stream typedef.
    typedef typename Operation::value_type value_type;

private:
    value_type current;

public:
    //! Construction.
    transform(Operation& o, Input1& i1_, Input2& i2_, Input3& i3_, Input4& i4_)
        : op(o), i1(i1_), i2(i2_), i3(i3_), i4(i4_)
    {
        if (!empty())
            current = op(*i1, *i2, *i3, *i4);
    }

    //! Standard stream method.
    const value_type& operator * () const
    {
        return current;
    }

    const value_type* operator -> () const
    {
        return &current;
    }

    //! Standard stream method.
    transform& operator ++ ()
    {
        ++i1;
        ++i2;
        ++i3;
        ++i4;
        if (!empty())
            current = op(*i1, *i2, *i3, *i4);

        return *this;
    }

    //! Standard stream method.
    bool empty() const
    {
        return i1.empty() || i2.empty() || i3.empty() || i4.empty();
    }
};

////////////////////////////////////////////////////////////////////////
//     TRANSFORM (5 input streams)                                    //
////////////////////////////////////////////////////////////////////////

//! Processes 5 input streams using given operation functor.
//!
//! \tparam Operation type of the operation (type of an
//! adaptable functor that takes 5 parameters)
//! \tparam Input1 type of the 1st input
//! \tparam Input2 type of the 2nd input
//! \tparam Input3 type of the 3rd input
//! \tparam Input4 type of the 4th input
//! \tparam Input5 type of the 5th input
//! \remark This is a specialization of \c transform .
template <class Operation,
          class Input1,
          class Input2,
          class Input3,
          class Input4,
          class Input5
          >
class transform<Operation, Input1, Input2, Input3, Input4, Input5, Stopper>
{
    Operation& op;
    Input1& i1;
    Input2& i2;
    Input3& i3;
    Input4& i4;
    Input5& i5;

public:
    //! Standard stream typedef.
    typedef typename Operation::value_type value_type;

private:
    value_type current;

public:
    //! Construction.
    transform(Operation& o, Input1& i1_, Input2& i2_, Input3& i3_,
              Input4& i4_, Input5& i5_)
        : op(o), i1(i1_), i2(i2_), i3(i3_), i4(i4_), i5(i5_)
    {
        if (!empty())
            current = op(*i1, *i2, *i3, *i4, *i5);
    }

    //! Standard stream method.
    const value_type& operator * () const
    {
        return current;
    }

    const value_type* operator -> () const
    {
        return &current;
    }

    //! Standard stream method.
    transform& operator ++ ()
    {
        ++i1;
        ++i2;
        ++i3;
        ++i4;
        ++i5;
        if (!empty())
            current = op(*i1, *i2, *i3, *i4, *i5);

        return *this;
    }

    //! Standard stream method.
    bool empty() const
    {
        return i1.empty() || i2.empty() || i3.empty() || i4.empty() || i5.empty();
    }
};

////////////////////////////////////////////////////////////////////////
//     MAKE TUPLE                                                     //
////////////////////////////////////////////////////////////////////////

//! Creates stream of 6-tuples from 6 input streams.
//!
//! \tparam Input1 type of the 1st input
//! \tparam Input2 type of the 2nd input
//! \tparam Input3 type of the 3rd input
//! \tparam Input4 type of the 4th input
//! \tparam Input5 type of the 5th input
//! \tparam Input6 type of the 6th input
template <class Input1,
          class Input2,
          class Input3 = Stopper,
          class Input4 = Stopper,
          class Input5 = Stopper,
          class Input6 = Stopper
          >
class make_tuple
{
    Input1& i1;
    Input2& i2;
    Input3& i3;
    Input4& i4;
    Input5& i5;
    Input6& i6;

public:
    //! Standard stream typedef.
    typedef typename stxxl::tuple<
            typename Input1::value_type,
            typename Input2::value_type,
            typename Input3::value_type,
            typename Input4::value_type,
            typename Input5::value_type,
            typename Input6::value_type
            > value_type;

private:
    value_type current;

public:
    //! Construction.
    make_tuple(Input1& i1_,
               Input2& i2_,
               Input3& i3_,
               Input4& i4_,
               Input5& i5_,
               Input6& i6_)
        : i1(i1_), i2(i2_), i3(i3_), i4(i4_), i5(i5_), i6(i6_)
    {
        if (!empty())
            current = value_type(*i1, *i2, *i3, *i4, *i5, *i6);
    }

    //! Standard stream method.
    const value_type& operator * () const
    {
        return current;
    }

    const value_type* operator -> () const
    {
        return &current;
    }

    //! Standard stream method.
    make_tuple& operator ++ ()
    {
        ++i1;
        ++i2;
        ++i3;
        ++i4;
        ++i5;
        ++i6;

        if (!empty())
            current = value_type(*i1, *i2, *i3, *i4, *i5, *i6);

        return *this;
    }

    //! Standard stream method.
    bool empty() const
    {
        return i1.empty() || i2.empty() || i3.empty() ||
               i4.empty() || i5.empty() || i6.empty();
    }
};

//! Creates stream of 2-tuples (pairs) from 2 input streams.
//!
//! \tparam Input1 type of the 1st input
//! \tparam Input2 type of the 2nd input
//! \remark A specialization of \c make_tuple .
template <class Input1,
          class Input2
          >
class make_tuple<Input1, Input2, Stopper, Stopper, Stopper, Stopper>
{
    Input1& i1;
    Input2& i2;

public:
    //! Standard stream typedef.
    typedef typename stxxl::tuple<
            typename Input1::value_type,
            typename Input2::value_type
            > value_type;

private:
    value_type current;

public:
    //! Construction.
    make_tuple(Input1& i1_,
               Input2& i2_)
        : i1(i1_), i2(i2_)
    {
        if (!empty())
            current = value_type(*i1, *i2);
    }

    //! Standard stream method.
    const value_type& operator * () const
    {
        return current;
    }

    const value_type* operator -> () const
    {
        return &current;
    }

    //! Standard stream method.
    make_tuple& operator ++ ()
    {
        ++i1;
        ++i2;

        if (!empty())
            current = value_type(*i1, *i2);

        return *this;
    }

    //! Standard stream method.
    bool empty() const
    {
        return i1.empty() || i2.empty();
    }
};

//! Creates stream of 3-tuples from 3 input streams.
//!
//! \tparam Input1 type of the 1st input
//! \tparam Input2 type of the 2nd input
//! \tparam Input3 type of the 3rd input
//! \remark A specialization of \c make_tuple .
template <class Input1,
          class Input2,
          class Input3
          >
class make_tuple<Input1, Input2, Input3, Stopper, Stopper, Stopper>
{
    Input1& i1;
    Input2& i2;
    Input3& i3;

public:
    //! Standard stream typedef.
    typedef typename stxxl::tuple<
            typename Input1::value_type,
            typename Input2::value_type,
            typename Input3::value_type
            > value_type;

private:
    value_type current;

public:
    //! Construction.
    make_tuple(Input1& i1_,
               Input2& i2_,
               Input3& i3_)
        : i1(i1_), i2(i2_), i3(i3_)
    {
        if (!empty())
            current = value_type(*i1, *i2, *i3);
    }

    //! Standard stream method.
    const value_type& operator * () const
    {
        return current;
    }

    const value_type* operator -> () const
    {
        return &current;
    }

    //! Standard stream method.
    make_tuple& operator ++ ()
    {
        ++i1;
        ++i2;
        ++i3;

        if (!empty())
            current = value_type(*i1, *i2, *i3);

        return *this;
    }

    //! Standard stream method.
    bool empty() const
    {
        return i1.empty() || i2.empty() || i3.empty();
    }
};

//! Creates stream of 4-tuples from 4 input streams.
//!
//! \tparam Input1 type of the 1st input
//! \tparam Input2 type of the 2nd input
//! \tparam Input3 type of the 3rd input
//! \tparam Input4 type of the 4th input
//! \remark A specialization of \c make_tuple .
template <class Input1,
          class Input2,
          class Input3,
          class Input4
          >
class make_tuple<Input1, Input2, Input3, Input4, Stopper, Stopper>
{
    Input1& i1;
    Input2& i2;
    Input3& i3;
    Input4& i4;

public:
    //! Standard stream typedef.
    typedef typename stxxl::tuple<
            typename Input1::value_type,
            typename Input2::value_type,
            typename Input3::value_type,
            typename Input4::value_type
            > value_type;

private:
    value_type current;

public:
    //! Construction.
    make_tuple(Input1& i1_,
               Input2& i2_,
               Input3& i3_,
               Input4& i4_)
        : i1(i1_), i2(i2_), i3(i3_), i4(i4_)
    {
        if (!empty())
            current = value_type(*i1, *i2, *i3, *i4);
    }

    //! Standard stream method.
    const value_type& operator * () const
    {
        return current;
    }

    const value_type* operator -> () const
    {
        return &current;
    }

    //! Standard stream method.
    make_tuple& operator ++ ()
    {
        ++i1;
        ++i2;
        ++i3;
        ++i4;

        if (!empty())
            current = value_type(*i1, *i2, *i3, *i4);

        return *this;
    }

    //! Standard stream method.
    bool empty() const
    {
        return i1.empty() || i2.empty() || i3.empty() ||
               i4.empty();
    }
};

//! Creates stream of 5-tuples from 5 input streams.
//!
//! \tparam Input1 type of the 1st input
//! \tparam Input2 type of the 2nd input
//! \tparam Input3 type of the 3rd input
//! \tparam Input4 type of the 4th input
//! \tparam Input5 type of the 5th input
//! \remark A specialization of \c make_tuple .
template <
    class Input1,
    class Input2,
    class Input3,
    class Input4,
    class Input5
    >
class make_tuple<Input1, Input2, Input3, Input4, Input5, Stopper>
{
    Input1& i1;
    Input2& i2;
    Input3& i3;
    Input4& i4;
    Input5& i5;

public:
    //! Standard stream typedef.
    typedef typename stxxl::tuple<
            typename Input1::value_type,
            typename Input2::value_type,
            typename Input3::value_type,
            typename Input4::value_type,
            typename Input5::value_type
            > value_type;

private:
    value_type current;

public:
    //! Construction.
    make_tuple(Input1& i1_,
               Input2& i2_,
               Input3& i3_,
               Input4& i4_,
               Input5& i5_)
        : i1(i1_), i2(i2_), i3(i3_), i4(i4_), i5(i5_)
    {
        if (!empty())
            current = value_type(*i1, *i2, *i3, *i4, *i5);
    }

    //! Standard stream method.
    const value_type& operator * () const
    {
        return current;
    }

    const value_type* operator -> () const
    {
        return &current;
    }

    //! Standard stream method.
    make_tuple& operator ++ ()
    {
        ++i1;
        ++i2;
        ++i3;
        ++i4;
        ++i5;

        if (!empty())
            current = value_type(*i1, *i2, *i3, *i4, *i5);

        return *this;
    }

    //! Standard stream method.
    bool empty() const
    {
        return i1.empty() || i2.empty() || i3.empty() ||
               i4.empty() || i5.empty();
    }
};

//! \}

} // namespace stream

STXXL_END_NAMESPACE

#include <stxxl/bits/stream/choose.h>
#include <stxxl/bits/stream/unique.h>

#endif // !STXXL_STREAM_STREAM_HEADER
// vim: et:ts=4:sw=4
