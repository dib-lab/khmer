/***************************************************************************
 *  include/stxxl/bits/containers/vector.h
 *
 *  Part of the STXXL. See http://stxxl.sourceforge.net
 *
 *  Copyright (C) 2002-2008 Roman Dementiev <dementiev@mpi-sb.mpg.de>
 *  Copyright (C) 2007-2009 Johannes Singler <singler@ira.uka.de>
 *  Copyright (C) 2008-2010 Andreas Beckmann <beckmann@cs.uni-frankfurt.de>
 *  Copyright (C) 2013 Timo Bingmann <tb@panthema.net>
 *
 *  Distributed under the Boost Software License, Version 1.0.
 *  (See accompanying file LICENSE_1_0.txt or copy at
 *  http://www.boost.org/LICENSE_1_0.txt)
 **************************************************************************/

#ifndef STXXL_CONTAINERS_VECTOR_HEADER
#define STXXL_CONTAINERS_VECTOR_HEADER

#include <vector>
#include <queue>
#include <algorithm>

#include <stxxl/bits/deprecated.h>
#include <stxxl/bits/io/request_operations.h>
#include <stxxl/bits/mng/block_manager.h>
#include <stxxl/bits/mng/typed_block.h>
#include <stxxl/bits/common/tmeta.h>
#include <stxxl/bits/containers/pager.h>
#include <stxxl/bits/common/is_sorted.h>
#include <stxxl/bits/mng/buf_istream.h>
#include <stxxl/bits/mng/buf_istream_reverse.h>
#include <stxxl/bits/mng/buf_ostream.h>

STXXL_BEGIN_NAMESPACE

#define STXXL_VERBOSE_VECTOR(msg) STXXL_VERBOSE1("vector[" << static_cast<const void*>(this) << "]::" << msg)

//! \defgroup stlcont Containers
//! \ingroup stllayer
//! Containers with STL-compatible interface

//! \defgroup stlcont_vector vector
//! \ingroup stlcont
//! Vector and support classes
//! \{

template <typename SizeType, SizeType modulo2, SizeType modulo1>
class double_blocked_index
{
    typedef SizeType size_type;

    static const size_type modulo12 = modulo1 * modulo2;

    size_type pos;
    unsigned_type block1, block2, offset;

    //! \invariant block2 * modulo12 + block1 * modulo1 + offset == pos && 0 <= block1 &lt; modulo2 && 0 <= offset &lt; modulo1

    void set(size_type pos)
    {
        this->pos = pos;
        block2 = (int_type)(pos / modulo12);
        pos -= block2 * modulo12;
        block1 = (int_type)(pos / modulo1);
        offset = (int_type)(pos - block1 * modulo1);

        assert(block2 * modulo12 + block1 * modulo1 + offset == this->pos);
        assert(/* 0 <= block1 && */ block1 < modulo2);
        assert(/* 0 <= offset && */ offset < modulo1);
    }

public:
    double_blocked_index()
    {
        set(0);
    }

    double_blocked_index(size_type pos)
    {
        set(pos);
    }

    double_blocked_index(unsigned_type block2, unsigned_type block1, unsigned_type offset)
    {
        assert(/* 0 <= block1 && */ block1 < modulo2);
        assert(/* 0 <= offset && */ offset < modulo1);

        this->block2 = block2;
        this->block1 = block1;
        this->offset = offset;
        pos = block2 * modulo12 + block1 * modulo1 + offset;
    }

    double_blocked_index& operator = (size_type pos)
    {
        set(pos);
        return *this;
    }

    //pre-increment operator
    double_blocked_index& operator ++ ()
    {
        ++pos;
        ++offset;
        if (offset == modulo1)
        {
            offset = 0;
            ++block1;
            if (block1 == modulo2)
            {
                block1 = 0;
                ++block2;
            }
        }

        assert(block2 * modulo12 + block1 * modulo1 + offset == this->pos);
        assert(/* 0 <= block1 && */ block1 < modulo2);
        assert(/* 0 <= offset && */ offset < modulo1);

        return *this;
    }

    //post-increment operator
    double_blocked_index operator ++ (int)
    {
        double_blocked_index former(*this);
        operator ++ ();
        return former;
    }

    //pre-increment operator
    double_blocked_index& operator -- ()
    {
        --pos;
        if (offset == 0)
        {
            offset = modulo1;
            if (block1 == 0)
            {
                block1 = modulo2;
                --block2;
            }
            --block1;
        }
        --offset;

        assert(block2 * modulo12 + block1 * modulo1 + offset == this->pos);
        assert(/*0 <= block1 &&*/ block1 < modulo2);
        assert(/*0 <= offset &&*/ offset < modulo1);

        return *this;
    }

    //post-increment operator
    double_blocked_index operator -- (int)
    {
        double_blocked_index former(*this);
        operator -- ();
        return former;
    }

    double_blocked_index operator + (size_type addend) const
    {
        return double_blocked_index(pos + addend);
    }

    double_blocked_index& operator += (size_type addend)
    {
        set(pos + addend);
        return *this;
    }

    double_blocked_index operator - (size_type addend) const
    {
        return double_blocked_index(pos - addend);
    }

    size_type operator - (const double_blocked_index& dbi2) const
    {
        return pos - dbi2.pos;
    }

    double_blocked_index& operator -= (size_type subtrahend)
    {
        set(pos - subtrahend);
        return *this;
    }

    bool operator == (const double_blocked_index& dbi2) const
    {
        return pos == dbi2.pos;
    }

    bool operator != (const double_blocked_index& dbi2) const
    {
        return pos != dbi2.pos;
    }

    bool operator < (const double_blocked_index& dbi2) const
    {
        return pos < dbi2.pos;
    }

    bool operator <= (const double_blocked_index& dbi2) const
    {
        return pos <= dbi2.pos;
    }

    bool operator > (const double_blocked_index& dbi2) const
    {
        return pos > dbi2.pos;
    }

    bool operator >= (const double_blocked_index& dbi2) const
    {
        return pos >= dbi2.pos;
    }

    double_blocked_index& operator >>= (size_type shift)
    {
        set(pos >> shift);
        return *this;
    }

    size_type get_pos() const
    {
        return pos;
    }

    const unsigned_type & get_block2() const
    {
        return block2;
    }

    const unsigned_type & get_block1() const
    {
        return block1;
    }

    const unsigned_type & get_offset() const
    {
        return offset;
    }
};

////////////////////////////////////////////////////////////////////////////

template <
    typename ValueType,
    unsigned PageSize,
    typename PagerType,
    unsigned BlockSize,
    typename AllocStr,
    typename SizeType>
class vector;

template <typename ValueType, typename AllocStr, typename SizeType, typename DiffType,
          unsigned BlockSize, typename PagerType, unsigned PageSize>
class const_vector_iterator;

template <typename VectorIteratorType>
class vector_bufreader;

template <typename VectorIteratorType>
class vector_bufreader_reverse;

template <typename VectorIteratorType>
class vector_bufwriter;

////////////////////////////////////////////////////////////////////////////

//! External vector iterator, model of \c ext_random_access_iterator concept.
template <typename ValueType, typename AllocStr, typename SizeType, typename DiffType,
          unsigned BlockSize, typename PagerType, unsigned PageSize>
class vector_iterator
{
    typedef vector_iterator<ValueType, AllocStr, SizeType,
                            DiffType, BlockSize, PagerType, PageSize> self_type;

    typedef const_vector_iterator<ValueType, AllocStr, SizeType,
                                  DiffType, BlockSize, PagerType, PageSize> const_self_type;

    friend class const_vector_iterator<ValueType, AllocStr, SizeType,
                                       DiffType, BlockSize, PagerType, PageSize>;

public:
    //! \name Types
    //! \{

    typedef self_type iterator;
    typedef const_self_type const_iterator;

    typedef unsigned block_offset_type;
    typedef vector<ValueType, PageSize, PagerType, BlockSize, AllocStr, SizeType> vector_type;
    friend class vector<ValueType, PageSize, PagerType, BlockSize, AllocStr, SizeType>;
    typedef typename vector_type::bids_container_type bids_container_type;
    typedef typename bids_container_type::iterator bids_container_iterator;
    typedef typename bids_container_type::bid_type bid_type;
    typedef typename vector_type::block_type block_type;
    typedef typename vector_type::blocked_index_type blocked_index_type;

    typedef std::random_access_iterator_tag iterator_category;
    typedef typename vector_type::size_type size_type;
    typedef typename vector_type::difference_type difference_type;
    typedef typename vector_type::value_type value_type;
    typedef typename vector_type::reference reference;
    typedef typename vector_type::const_reference const_reference;
    typedef typename vector_type::pointer pointer;
    typedef typename vector_type::const_pointer const_pointer;

    //! \}

protected:
    blocked_index_type offset;
    vector_type* p_vector;

private:
    //! private constructor for initializing other iterators
    vector_iterator(vector_type* v, size_type o)
        : offset(o), p_vector(v)
    { }

public:
    //! constructs invalid iterator
    vector_iterator()
        : offset(0), p_vector(NULL)
    { }
    //! copy-constructor
    vector_iterator(const self_type& a)
        : offset(a.offset),
          p_vector(a.p_vector)
    { }

    //! \name Iterator Properties
    //! \{

    //! return pointer to vector containing iterator
    vector_type * parent_vector() const
    {
        return p_vector;
    }
    //! return block offset of current element
    block_offset_type block_offset() const
    {
        return static_cast<block_offset_type>(offset.get_offset());
    }
    //! return iterator to BID containg current element
    bids_container_iterator bid() const
    {
        return p_vector->bid(offset);
    }

    //! \}

    //! \name Access Operators
    //! \{

    //! return current element
    reference operator * ()
    {
        return p_vector->element(offset);
    }
    //! return pointer to current element
    pointer operator -> ()
    {
        return &(p_vector->element(offset));
    }
    //! return const reference to current element
    const_reference operator * () const
    {
        return p_vector->const_element(offset);
    }
    //! return const pointer to current element
    const_pointer operator -> () const
    {
        return &(p_vector->const_element(offset));
    }
    //! return mutable reference to element +i after the current element
    reference operator [] (size_type i)
    {
        return p_vector->element(offset.get_pos() + i);
    }

#ifdef _LIBCPP_VERSION
    //-tb 2013-11: libc++ defines std::reverse_iterator::operator[] in such a
    // way that it expects vector_iterator::operator[] to return a (mutable)
    // reference. Thus to remove confusion about the compiler error, we remove
    // the operator[] const for libc++. The const_reference actually violates
    // some version of the STL standard, but works well in gcc's libstdc++.
#else
    //! return const reference to element +i after the current element
    const_reference operator [] (size_type i) const
    {
        return p_vector->const_element(offset.get_pos() + i);
    }
#endif
    //! \}

    //! \name Relative Calculation of Iterators
    //! \{

    //! calculate different between two iterator
    difference_type operator - (const self_type& a) const
    {
        return offset - a.offset;
    }
    //! calculate different between two iterator
    difference_type operator - (const const_self_type& a) const
    {
        return offset - a.offset;
    }
    //! return iterator advanced -i positions in the vector
    self_type operator - (size_type i) const
    {
        return self_type(p_vector, offset.get_pos() - i);
    }
    //! return iterator advanced +i positions in the vector
    self_type operator + (size_type i) const
    {
        return self_type(p_vector, offset.get_pos() + i);
    }
    //! advance this iterator -i positions in the vector
    self_type& operator -= (size_type i)
    {
        offset -= i;
        return *this;
    }
    //! advance this iterator +i positions in the vector
    self_type& operator += (size_type i)
    {
        offset += i;
        return *this;
    }
    //! advance this iterator to next position in the vector
    self_type& operator ++ ()
    {
        offset++;
        return *this;
    }
    //! advance this iterator to next position in the vector
    self_type operator ++ (int)
    {
        self_type tmp = *this;
        offset++;
        return tmp;
    }
    //! advance this iterator to preceding position in the vector
    self_type& operator -- ()
    {
        offset--;
        return *this;
    }
    //! advance this iterator to preceding position in the vector
    self_type operator -- (int)
    {
        self_type tmp = *this;
        offset--;
        return tmp;
    }

    //! \}

    //! \name Comparison Operators
    //! \{

    bool operator == (const self_type& a) const
    {
        assert(p_vector == a.p_vector);
        return offset == a.offset;
    }
    bool operator != (const self_type& a) const
    {
        assert(p_vector == a.p_vector);
        return offset != a.offset;
    }
    bool operator < (const self_type& a) const
    {
        assert(p_vector == a.p_vector);
        return offset < a.offset;
    }
    bool operator <= (const self_type& a) const
    {
        assert(p_vector == a.p_vector);
        return offset <= a.offset;
    }
    bool operator > (const self_type& a) const
    {
        assert(p_vector == a.p_vector);
        return offset > a.offset;
    }
    bool operator >= (const self_type& a) const
    {
        assert(p_vector == a.p_vector);
        return offset >= a.offset;
    }

    bool operator == (const const_self_type& a) const
    {
        assert(p_vector == a.p_vector);
        return offset == a.offset;
    }
    bool operator != (const const_self_type& a) const
    {
        assert(p_vector == a.p_vector);
        return offset != a.offset;
    }
    bool operator < (const const_self_type& a) const
    {
        assert(p_vector == a.p_vector);
        return offset < a.offset;
    }
    bool operator <= (const const_self_type& a) const
    {
        assert(p_vector == a.p_vector);
        return offset <= a.offset;
    }
    bool operator > (const const_self_type& a) const
    {
        assert(p_vector == a.p_vector);
        return offset > a.offset;
    }
    bool operator >= (const const_self_type& a) const
    {
        assert(p_vector == a.p_vector);
        return offset >= a.offset;
    }

    //! \}

    //! \name Flushing Operation
    //! \{

    void block_externally_updated()
    {
        p_vector->block_externally_updated(offset);
    }

    void flush()
    {
        p_vector->flush();
    }

    //! \}
};

////////////////////////////////////////////////////////////////////////////

//! Const external vector iterator, model of \c ext_random_access_iterator concept.
template <typename ValueType, typename AllocStr, typename SizeType, typename DiffType,
          unsigned BlockSize, typename PagerType, unsigned PageSize>
class const_vector_iterator
{
    typedef const_vector_iterator<ValueType, AllocStr, SizeType, DiffType,
                                  BlockSize, PagerType, PageSize> self_type;

    typedef vector_iterator<ValueType, AllocStr, SizeType, DiffType,
                            BlockSize, PagerType, PageSize> mutable_self_type;

    friend class vector_iterator<ValueType, AllocStr, SizeType, DiffType,
                                 BlockSize, PagerType, PageSize>;

public:
    //! \name Types
    //! \{

    typedef self_type const_iterator;
    typedef mutable_self_type iterator;

    typedef unsigned block_offset_type;
    typedef vector<ValueType, PageSize, PagerType, BlockSize, AllocStr, SizeType> vector_type;
    friend class vector<ValueType, PageSize, PagerType, BlockSize, AllocStr, SizeType>;
    typedef typename vector_type::bids_container_type bids_container_type;
    typedef typename bids_container_type::iterator bids_container_iterator;
    typedef typename bids_container_type::bid_type bid_type;
    typedef typename vector_type::block_type block_type;
    typedef typename vector_type::blocked_index_type blocked_index_type;

    typedef std::random_access_iterator_tag iterator_category;
    typedef typename vector_type::size_type size_type;
    typedef typename vector_type::difference_type difference_type;
    typedef typename vector_type::value_type value_type;
    typedef typename vector_type::const_reference reference;
    typedef typename vector_type::const_reference const_reference;
    typedef typename vector_type::const_pointer pointer;
    typedef typename vector_type::const_pointer const_pointer;

    //! \}

protected:
    blocked_index_type offset;
    const vector_type* p_vector;

private:
    //! private constructor for initializing other iterators
    const_vector_iterator(const vector_type* v, size_type o)
        : offset(o), p_vector(v)
    { }

public:
    //! constructs invalid iterator
    const_vector_iterator()
        : offset(0), p_vector(NULL)
    { }
    //! copy-constructor
    const_vector_iterator(const self_type& a)
        : offset(a.offset), p_vector(a.p_vector)
    { }
    //! copy-constructor from mutable iterator
    const_vector_iterator(const mutable_self_type& a)
        : offset(a.offset), p_vector(a.p_vector)
    { }

    //! \name Iterator Properties
    //! \{

    //! return pointer to vector containing iterator
    const vector_type * parent_vector() const
    {
        return p_vector;
    }
    //! return block offset of current element
    block_offset_type block_offset() const
    {
        return static_cast<block_offset_type>(offset.get_offset());
    }
    //! return iterator to BID containg current element
    bids_container_iterator bid() const
    {
        return ((vector_type*)p_vector)->bid(offset);
    }

    //! \}

    //! \name Access Operators
    //! \{

    //! return current element
    const_reference operator * () const
    {
        return p_vector->const_element(offset);
    }
    //! return pointer to current element
    const_pointer operator -> () const
    {
        return &(p_vector->const_element(offset));
    }
    //! return const reference to element +i after the current element
    const_reference operator [] (size_type i) const
    {
        return p_vector->const_element(offset.get_pos() + i);
    }

    //! \}

    //! \name Relative Calculation of Iterators
    //! \{

    //! calculate different between two iterator
    difference_type operator - (const self_type& a) const
    {
        return offset - a.offset;
    }
    //! calculate different between two iterator
    difference_type operator - (const mutable_self_type& a) const
    {
        return offset - a.offset;
    }
    //! return iterator advanced -i positions in the vector
    self_type operator - (size_type i) const
    {
        return self_type(p_vector, offset.get_pos() - i);
    }
    //! return iterator advanced +i positions in the vector
    self_type operator + (size_type i) const
    {
        return self_type(p_vector, offset.get_pos() + i);
    }
    //! advance this iterator -i positions in the vector
    self_type& operator -= (size_type i)
    {
        offset -= i;
        return *this;
    }
    //! advance this iterator +i positions in the vector
    self_type& operator += (size_type i)
    {
        offset += i;
        return *this;
    }
    //! advance this iterator to next position in the vector
    self_type& operator ++ ()
    {
        offset++;
        return *this;
    }
    //! advance this iterator to next position in the vector
    self_type operator ++ (int)
    {
        self_type tmp_ = *this;
        offset++;
        return tmp_;
    }
    //! advance this iterator to preceding position in the vector
    self_type& operator -- ()
    {
        offset--;
        return *this;
    }
    //! advance this iterator to preceding position in the vector
    self_type operator -- (int)
    {
        self_type tmp = *this;
        offset--;
        return tmp;
    }

    //! \}

    //! \name Comparison Operators
    //! \{

    bool operator == (const self_type& a) const
    {
        assert(p_vector == a.p_vector);
        return offset == a.offset;
    }
    bool operator != (const self_type& a) const
    {
        assert(p_vector == a.p_vector);
        return offset != a.offset;
    }
    bool operator < (const self_type& a) const
    {
        assert(p_vector == a.p_vector);
        return offset < a.offset;
    }
    bool operator <= (const self_type& a) const
    {
        assert(p_vector == a.p_vector);
        return offset <= a.offset;
    }
    bool operator > (const self_type& a) const
    {
        assert(p_vector == a.p_vector);
        return offset > a.offset;
    }
    bool operator >= (const self_type& a) const
    {
        assert(p_vector == a.p_vector);
        return offset >= a.offset;
    }

    bool operator == (const mutable_self_type& a) const
    {
        assert(p_vector == a.p_vector);
        return offset == a.offset;
    }
    bool operator != (const mutable_self_type& a) const
    {
        assert(p_vector == a.p_vector);
        return offset != a.offset;
    }
    bool operator < (const mutable_self_type& a) const
    {
        assert(p_vector == a.p_vector);
        return offset < a.offset;
    }
    bool operator <= (const mutable_self_type& a) const
    {
        assert(p_vector == a.p_vector);
        return offset <= a.offset;
    }
    bool operator > (const mutable_self_type& a) const
    {
        assert(p_vector == a.p_vector);
        return offset > a.offset;
    }
    bool operator >= (const mutable_self_type& a) const
    {
        assert(p_vector == a.p_vector);
        return offset >= a.offset;
    }

    //! \}

    //! \name Flushing Operation
    //! \{

    void block_externally_updated()
    {
        p_vector->block_externally_updated(offset);
    }

    void flush()
    {
        p_vector->flush();
    }

    //! \}
};

////////////////////////////////////////////////////////////////////////////

//! External vector container. \n
//! <b>Introduction</b> to vector container: see \ref tutorial_vector tutorial. \n
//! <b>Design and Internals</b> of vector container: see \ref design_vector
//!
//! For semantics of the methods see documentation of the STL std::vector
//! \tparam ValueType type of contained objects (POD with no references to internal memory)
//! \tparam PageSize number of blocks in a page
//! \tparam PagerType pager type, \c random_pager<x> or \c lru_pager<x>, where x is the default number of pages,
//!  default is \c lru_pager<8>
//! \tparam BlockSize external block size in bytes, default is 2 MiB
//! \tparam AllocStr one of allocation strategies: \c striping , \c RC , \c SR , or \c FR
//!  default is RC
//!
//! Memory consumption: BlockSize*x*PageSize bytes
//! \warning Do not store references to the elements of an external vector. Such references
//! might be invalidated during any following access to elements of the vector
template <
    typename ValueType,
    unsigned PageSize = 4,
    typename PagerType = lru_pager<8>,
    unsigned BlockSize = STXXL_DEFAULT_BLOCK_SIZE(ValueType),
    typename AllocStr = STXXL_DEFAULT_ALLOC_STRATEGY,
    typename SizeType = stxxl::uint64       // will be deprecated soon
    >
class vector
{
public:
    //! \name Standard Types
    //! \{

    //! The type of elements stored in the vector.
    typedef ValueType value_type;
    //! reference to value_type
    typedef value_type& reference;
    //! constant reference to value_type
    typedef const value_type& const_reference;
    //! pointer to value_type
    typedef value_type* pointer;
    //! constant pointer to value_type
    typedef const value_type* const_pointer;
    //! an unsigned 64-bit integral type
    typedef SizeType size_type;
    typedef stxxl::int64 difference_type;

    typedef PagerType pager_type;
    typedef AllocStr alloc_strategy_type;

    enum constants {
        block_size = BlockSize,
        page_size = PageSize,
        on_disk = -1
    };

    //! iterator used to iterate through a vector, see \ref design_vector_notes.
    typedef vector_iterator<value_type, alloc_strategy_type, size_type,
                            difference_type, block_size, pager_type, page_size> iterator;
    friend class vector_iterator<value_type, alloc_strategy_type, size_type, difference_type, block_size, pager_type, page_size>;

    //! constant iterator used to iterate through a vector, see \ref design_vector_notes.
    typedef const_vector_iterator<value_type, alloc_strategy_type,
                                  size_type, difference_type, block_size, pager_type, page_size> const_iterator;
    friend class const_vector_iterator<value_type, alloc_strategy_type, size_type, difference_type, block_size, pager_type, page_size>;
    typedef std::reverse_iterator<iterator> reverse_iterator;
    typedef std::reverse_iterator<const_iterator> const_reverse_iterator;

    //! \}

    //! \name Extra Types
    //! \{

    //! vector_bufwriter compatible with this vector
    typedef vector_bufwriter<iterator> bufwriter_type;

    //! vector_bufreader compatible with this vector
    typedef vector_bufreader<const_iterator> bufreader_type;

    //! vector_bufreader compatible with this vector
    typedef vector_bufreader_reverse<const_iterator> bufreader_reverse_type;

    //! \internal
    class bid_vector : public std::vector<BID<block_size> >
    {
    public:
        typedef std::vector<BID<block_size> > super_type;
        typedef typename super_type::size_type size_type;
        typedef typename super_type::value_type bid_type;

        bid_vector(size_type sz) : super_type(sz)
        { }
    };

    typedef bid_vector bids_container_type;
    typedef typename bids_container_type::iterator bids_container_iterator;
    typedef typename bids_container_type::const_iterator const_bids_container_iterator;

    //! type of the block used in disk-memory transfers
    typedef typed_block<BlockSize, ValueType> block_type;
    //! double-index type to reference individual elements in a block
    typedef double_blocked_index<SizeType, PageSize, block_type::size> blocked_index_type;

    //! \}

private:
    alloc_strategy_type m_alloc_strategy;
    size_type m_size;
    bids_container_type m_bids;
    mutable pager_type m_pager;

    // enum specifying status of a page of the vector
    enum { valid_on_disk = 0, uninitialized = 1, dirty = 2 };
    //! status of each page (valid_on_disk, uninitialized or dirty)
    mutable std::vector<unsigned char> m_page_status;
    mutable std::vector<int_type> m_page_to_slot;
    mutable simple_vector<int_type> m_slot_to_page;
    mutable std::queue<int_type> m_free_slots;
    mutable simple_vector<block_type>* m_cache;
    file* m_from;
    block_manager* m_bm;
    bool m_exported;

    size_type size_from_file_length(stxxl::uint64 file_length) const
    {
        stxxl::uint64 blocks_fit = file_length / stxxl::uint64(block_type::raw_size);
        size_type cur_size = blocks_fit * stxxl::uint64(block_type::size);
        stxxl::uint64 rest = file_length - blocks_fit * stxxl::uint64(block_type::raw_size);
        return (cur_size + rest / stxxl::uint64(sizeof(value_type)));
    }

    stxxl::uint64 file_length() const
    {
        typedef stxxl::uint64 file_size_type;
        size_type cur_size = size();
        size_type num_full_blocks = cur_size / block_type::size;
        if (cur_size % block_type::size != 0)
        {
            size_type rest = cur_size - num_full_blocks * block_type::size;
            return file_size_type(num_full_blocks) * block_type::raw_size + rest * sizeof(value_type);
        }
        return file_size_type(num_full_blocks) * block_type::raw_size;
    }

public:
    //! \name Constructors/Destructors
    //! \{

    //! Constructs external vector with n elements.
    //!
    //! \param n Number of elements.
    //! \param npages Number of cached pages.
    vector(size_type n = 0, unsigned_type npages = pager_type().size())
        : m_size(n),
          m_bids((size_t)div_ceil(n, block_type::size)),
          m_pager(npages),
          m_page_status(div_ceil(m_bids.size(), page_size)),
          m_page_to_slot(div_ceil(m_bids.size(), page_size)),
          m_slot_to_page(npages),
          m_cache(NULL),
          m_from(NULL),
          m_exported(false)
    {
        m_bm = block_manager::get_instance();

        allocate_page_cache();

        for (size_t i = 0; i < m_page_status.size(); ++i)
        {
            m_page_status[i] = uninitialized;
            m_page_to_slot[i] = on_disk;
        }

        for (unsigned_type i = 0; i < numpages(); ++i)
            m_free_slots.push(i);

        m_bm->new_blocks(m_alloc_strategy, m_bids.begin(), m_bids.end(), 0);
    }

    //! \}

    //! \name Modifier
    //! \{

    //! swap content
    void swap(vector& obj)
    {
        std::swap(m_alloc_strategy, obj.m_alloc_strategy);
        std::swap(m_size, obj.m_size);
        std::swap(m_bids, obj.m_bids);
        std::swap(m_pager, obj.m_pager);
        std::swap(m_page_status, obj.m_page_status);
        std::swap(m_page_to_slot, obj.m_page_to_slot);
        std::swap(m_slot_to_page, obj.m_slot_to_page);
        std::swap(m_free_slots, obj.m_free_slots);
        std::swap(m_cache, obj.m_cache);
        std::swap(m_from, obj.m_from);
        std::swap(m_exported, obj.m_exported);
    }

    //! \}

    //! \name Miscellaneous
    //! \{

    //! Allocate page cache, must be called to allow access to elements.
    void allocate_page_cache() const
    {
        //  numpages() might be zero
        if (!m_cache && numpages() > 0)
            m_cache = new simple_vector<block_type>(numpages() * page_size);
    }

    //! allows to free the cache, but you may not access any element until call
    //! allocate_page_cache() again
    void deallocate_page_cache() const
    {
        flush();
        delete m_cache;
        m_cache = NULL;
    }

    //! \name Size and Capacity
    //! \{

    //! return the size of the vector.
    size_type size() const
    {
        return m_size;
    }
    //! true if the vector's size is zero.
    bool empty() const
    {
        return (!m_size);
    }

    //! Return the number of elelemtsn for which \a external memory has been
    //! allocated. capacity() is always greator than or equal to size().
    size_type capacity() const
    {
        return size_type(m_bids.size()) * block_type::size;
    }
    //! Returns the number of bytes that the vector has allocated on disks.
    size_type raw_capacity() const
    {
        return size_type(m_bids.size()) * block_type::raw_size;
    }

    /*! Reserves at least n elements in external memory.
     *
     * If n is less than or equal to capacity(), this call has no
     * effect. Otherwise, it is a request for allocation of additional \b
     * external memory. If the request is successful, then capacity() is
     * greater than or equal to n; otherwise capacity() is unchanged. In either
     * case, size() is unchanged.
     */
    void reserve(size_type n)
    {
        if (n <= capacity())
            return;

        unsigned_type old_bids_size = m_bids.size();
        unsigned_type new_bids_size = (unsigned_type)div_ceil(n, block_type::size);
        unsigned_type new_pages = div_ceil(new_bids_size, page_size);
        m_page_status.resize(new_pages, uninitialized);
        m_page_to_slot.resize(new_pages, on_disk);

        m_bids.resize(new_bids_size);
        if (m_from == NULL)
        {
            m_bm->new_blocks(m_alloc_strategy,
                             m_bids.begin() + old_bids_size, m_bids.end(),
                             old_bids_size);
        }
        else
        {
            size_type offset = size_type(old_bids_size) * size_type(block_type::raw_size);
            for (bids_container_iterator it = m_bids.begin() + old_bids_size;
                 it != m_bids.end(); ++it, offset += size_type(block_type::raw_size))
            {
                (*it).storage = m_from;
                (*it).offset = offset;
            }
            STXXL_VERBOSE_VECTOR("reserve(): Changing size of file " <<
                                 ((void*)m_from) << " to " << offset);
            m_from->set_size(offset);
        }
    }

    //! Resize vector contents to n items.
    //! \warning this will not call the constructor of objects in external memory!
    void resize(size_type n)
    {
        _resize(n);
    }

    //! Resize vector contents to n items, and allow the allocated external
    //! memory to shrink. Internal memory allocation remains unchanged.
    //! \warning this will not call the constructor of objects in external memory!
    void resize(size_type n, bool shrink_capacity)
    {
        if (shrink_capacity)
            _resize_shrink_capacity(n);
        else
            _resize(n);
    }

    //! \}

private:
    //! Resize vector, only allow capacity growth.
    void _resize(size_type n)
    {
        reserve(n);
        if (n < m_size) {
            // mark excess pages as uninitialized and evict them from cache
            unsigned_type first_page_to_evict = (unsigned_type)div_ceil(n, block_type::size * page_size);
            for (size_t i = first_page_to_evict; i < m_page_status.size(); ++i) {
                if (m_page_to_slot[i] != on_disk) {
                    m_free_slots.push(m_page_to_slot[i]);
                    m_page_to_slot[i] = on_disk;
                }
                m_page_status[i] = uninitialized;
            }
        }
        m_size = n;
    }

    //! Resize vector, also allow reduction of external memory capacity.
    void _resize_shrink_capacity(size_type n)
    {
        unsigned_type old_bids_size = m_bids.size();
        unsigned_type new_bids_size = (unsigned_type)div_ceil(n, block_type::size);

        if (new_bids_size > old_bids_size)
        {
            reserve(n);
        }
        else if (new_bids_size < old_bids_size)
        {
            unsigned_type new_pages_size = div_ceil(new_bids_size, page_size);

            STXXL_VERBOSE_VECTOR("shrinking from " << old_bids_size << " to " <<
                                 new_bids_size << " blocks = from " <<
                                 m_page_status.size() << " to " <<
                                 new_pages_size << " pages");

            // release blocks
            if (m_from != NULL)
                m_from->set_size(new_bids_size * block_type::raw_size);
            else
                m_bm->delete_blocks(m_bids.begin() + old_bids_size, m_bids.end());

            m_bids.resize(new_bids_size);

            // don't resize m_page_to_slot or m_page_status, because it is
            // still needed to check page status and match the mapping
            // m_slot_to_page

            // clear dirty flag, so these pages will be never written
            std::fill(m_page_status.begin() + new_pages_size,
                      m_page_status.end(), (unsigned char)valid_on_disk);
        }

        m_size = n;
    }

public:
    //! \name Modifiers
    //! \{

    //! Erases all of the elements and deallocates all external memory that is
    //! occupied.
    void clear()
    {
        m_size = 0;
        if (m_from == NULL)
            m_bm->delete_blocks(m_bids.begin(), m_bids.end());

        m_bids.clear();
        m_page_status.clear();
        m_page_to_slot.clear();
        while (!m_free_slots.empty())
            m_free_slots.pop();

        for (unsigned_type i = 0; i < numpages(); ++i)
            m_free_slots.push(i);
    }

    //! \name Front and Back Access
    //! \{

    //! Append a new element at the end.
    void push_back(const_reference obj)
    {
        size_type old_size = m_size;
        resize(old_size + 1);
        element(old_size) = obj;
    }
    //! Removes the last element (without returning it, see back()).
    void pop_back()
    {
        resize(m_size - 1);
    }

    //! \}

    //! \name Operators
    //! \{

    //! Returns a reference to the last element, see \ref design_vector_notes.
    reference back()
    {
        return element(m_size - 1);
    }
    //! Returns a reference to the first element, see \ref design_vector_notes.
    reference front()
    {
        return element(0);
    }
    //! Returns a constant reference to the last element, see \ref design_vector_notes.
    const_reference back() const
    {
        return const_element(m_size - 1);
    }
    //! Returns a constant reference to the first element, see \ref design_vector_notes.
    const_reference front() const
    {
        return const_element(0);
    }

    //! \}

    //! \name Constructors/Destructors
    //! \{

    //! Construct vector from a file.
    //! \param from file to be constructed from
    //! \param size Number of elements.
    //! \param npages Number of cached pages.
    //! \warning Only one \c vector can be assigned to a particular (physical) file.
    //! The block size of the vector must be a multiple of the element size
    //! \c sizeof(ValueType) and the page size (4096).
    vector(file* from, size_type size = size_type(-1), unsigned_type npages = pager_type().size())
        : m_size((size == size_type(-1)) ? size_from_file_length(from->size()) : size),
          m_bids((size_t)div_ceil(m_size, size_type(block_type::size))),
          m_pager(npages),
          m_page_status(div_ceil(m_bids.size(), page_size)),
          m_page_to_slot(div_ceil(m_bids.size(), page_size)),
          m_slot_to_page(npages),
          m_cache(NULL),
          m_from(from),
          m_exported(false)
    {
        // initialize from file
        if (!block_type::has_only_data)
        {
            std::ostringstream str;
            str << "The block size for a vector that is mapped to a file must be a multiple of the element size (" <<
                sizeof(value_type) << ") and the page size (4096).";
            throw std::runtime_error(str.str());
        }

        m_bm = block_manager::get_instance();

        allocate_page_cache();

        for (size_t i = 0; i < m_page_status.size(); ++i)
        {
            m_page_status[i] = valid_on_disk;
            m_page_to_slot[i] = on_disk;
        }

        for (unsigned_type i = 0; i < numpages(); ++i)
            m_free_slots.push(i);

        // allocate blocks equidistantly and in-order
        size_type offset = 0;
        for (bids_container_iterator it = m_bids.begin();
             it != m_bids.end(); ++it, offset += size_type(block_type::raw_size))
        {
            (*it).storage = from;
            (*it).offset = offset;
        }
        from->set_size(offset);
    }

    //! copy-constructor
    vector(const vector& obj)
        : m_size(obj.size()),
          m_bids((size_t)div_ceil(obj.size(), block_type::size)),
          m_pager(obj.numpages()),
          m_page_status(div_ceil(m_bids.size(), page_size)),
          m_page_to_slot(div_ceil(m_bids.size(), page_size)),
          m_slot_to_page(obj.numpages()),
          m_cache(NULL),
          m_from(NULL),
          m_exported(false)
    {
        assert(!obj.m_exported);
        m_bm = block_manager::get_instance();

        allocate_page_cache();

        for (size_t i = 0; i < m_page_status.size(); ++i)
        {
            m_page_status[i] = uninitialized;
            m_page_to_slot[i] = on_disk;
        }

        for (unsigned_type i = 0; i < numpages(); ++i)
            m_free_slots.push(i);

        m_bm->new_blocks(m_alloc_strategy, m_bids.begin(), m_bids.end(), 0);

        const_iterator inbegin = obj.begin();
        const_iterator inend = obj.end();
        std::copy(inbegin, inend, begin());
    }

    //! \}

    //! \name Operators
    //! \{

    //! assignment operator
    vector& operator = (const vector& obj)
    {
        if (&obj != this)
        {
            vector tmp(obj);
            this->swap(tmp);
        }
        return *this;
    }

    //! \}

    //! \name Iterator Construction
    //! \{

    //! returns an iterator pointing to the beginning of the vector, see \ref design_vector_notes.
    iterator begin()
    {
        return iterator(this, 0);
    }
    //! returns a const_iterator pointing to the beginning of the vector, see \ref design_vector_notes.
    const_iterator begin() const
    {
        return const_iterator(this, 0);
    }
    //! returns a const_iterator pointing to the beginning of the vector, see \ref design_vector_notes.
    const_iterator cbegin() const
    {
        return begin();
    }
    //! returns an iterator pointing beyond the end of the vector, see \ref design_vector_notes.
    iterator end()
    {
        return iterator(this, m_size);
    }
    //! returns a const_iterator pointing beyond the end of the vector, see \ref design_vector_notes.
    const_iterator end() const
    {
        return const_iterator(this, m_size);
    }
    //! returns a const_iterator pointing beyond the end of the vector, see \ref design_vector_notes.
    const_iterator cend() const
    {
        return end();
    }

    //! returns a reverse_iterator pointing to the end of the vector.
    reverse_iterator rbegin()
    {
        return reverse_iterator(end());
    }
    //! returns a reverse_iterator pointing to the end of the vector.
    const_reverse_iterator rbegin() const
    {
        return const_reverse_iterator(end());
    }
    //! returns a reverse_iterator pointing to the end of the vector.
    const_reverse_iterator crbegin() const
    {
        return const_reverse_iterator(end());
    }
    //! returns a reverse_iterator pointing beyond the beginning of the vector.
    reverse_iterator rend()
    {
        return reverse_iterator(begin());
    }
    //! returns a reverse_iterator pointing beyond the beginning of the vector.
    const_reverse_iterator rend() const
    {
        return const_reverse_iterator(begin());
    }
    //! returns a reverse_iterator pointing beyond the beginning of the vector.
    const_reverse_iterator crend() const
    {
        return const_reverse_iterator(begin());
    }

    //! \}

    //! \name Direct Element Access
    //! \{

    //! access the element at the given vector's offset
    reference operator [] (size_type offset)
    {
        return element(offset);
    }
    //! access the element at the given vector's offset
    const_reference operator [] (size_type offset) const
    {
        return const_element(offset);
    }

    //! access the element at the given vector's offset
    reference at(size_type offset)
    {
        assert(offset < (size_type)size());
        return element(offset);
    }
    //! access the element at the given vector's offset
    const_reference at(size_type offset) const
    {
        assert(offset < (size_type)size());
        return const_element(offset);
    }

    //! return true if the given vector offset is in cache
    bool is_element_cached(size_type offset) const
    {
        return is_page_cached(blocked_index_type(offset));
    }

    //! \}

    //! \name Modifiers
    //! \{

    //! Flushes the cache pages to the external memory.
    void flush() const
    {
        simple_vector<bool> non_free_slots(numpages());

        for (unsigned_type i = 0; i < numpages(); i++)
            non_free_slots[i] = true;

        while (!m_free_slots.empty())
        {
            non_free_slots[m_free_slots.front()] = false;
            m_free_slots.pop();
        }

        for (unsigned_type i = 0; i < numpages(); i++)
        {
            m_free_slots.push(i);
            int_type page_no = m_slot_to_page[i];
            if (non_free_slots[i])
            {
                STXXL_VERBOSE_VECTOR("flush(): flushing page " << i << " at address " <<
                                     (int64(page_no) * int64(block_type::size) * int64(page_size)));
                write_page(page_no, i);

                m_page_to_slot[page_no] = on_disk;
            }
        }
    }

    //! \}

    //! \name Constructors/Destructors
    //! \{

    ~vector()
    {
        STXXL_VERBOSE_VECTOR("~vector()");
        try
        {
            flush();
        }
        catch (io_error e)
        {
            STXXL_ERRMSG("io_error thrown in ~vector(): " << e.what());
        }
        catch (...)
        {
            STXXL_ERRMSG("Exception thrown in ~vector()");
        }

        if (!m_exported)
        {
            if (m_from == NULL) {
                m_bm->delete_blocks(m_bids.begin(), m_bids.end());
            }
            else // file must be truncated
            {
                STXXL_VERBOSE_VECTOR("~vector(): Changing size of file " <<
                                     ((void*)m_from) << " to " << file_length());
                STXXL_VERBOSE_VECTOR("~vector(): size of the vector is " << size());
                try
                {
                    m_from->set_size(file_length());
                }
                catch (...)
                {
                    STXXL_ERRMSG("Exception thrown in ~vector()...set_size()");
                }
            }
        }
        delete m_cache;
    }

    //! \}

    //! \name Miscellaneous
    //! \{

    //! Export data such that it is persistent on the file system. Resulting
    //! files will be numbered ascending.
    void export_files(std::string filename_prefix)
    {
        int64 no = 0;
        for (bids_container_iterator i = m_bids.begin(); i != m_bids.end(); ++i) {
            std::ostringstream number;
            number << std::setw(9) << std::setfill('0') << no;
            size_type current_block_size =
                ((i + 1) == m_bids.end() && m_size % block_type::size > 0) ?
                (m_size % block_type::size) * sizeof(value_type) :
                block_type::size * sizeof(value_type);
            (*i).storage->export_files((*i).offset, current_block_size, filename_prefix + number.str());
            ++no;
        }
        m_exported = true;
    }

    //! Get the file associated with this vector, or NULL.
    file * get_file() const
    {
        return m_from;
    }

    //! \}

    //! \name Capacity
    //! \{

    //! Set the blocks and the size of this container explicitly.
    //! The vector must be completely empty before.
    template <typename ForwardIterator>
    void set_content(const ForwardIterator& bid_begin, const ForwardIterator& bid_end, size_type n)
    {
        unsigned_type new_bids_size = div_ceil(n, block_type::size);
        m_bids.resize(new_bids_size);
        std::copy(bid_begin, bid_end, m_bids.begin());
        unsigned_type new_pages = div_ceil(new_bids_size, page_size);
        m_page_status.resize(new_pages, valid_on_disk);
        m_page_to_slot.resize(new_pages, on_disk);
        m_size = n;
    }

    //! Number of pages used by the pager.
    inline unsigned_type numpages() const
    {
        return m_pager.size();
    }

    //! \}

private:
    bids_container_iterator bid(const size_type& offset)
    {
        return (m_bids.begin() +
                static_cast<typename bids_container_type::size_type>
                (offset / block_type::size));
    }
    bids_container_iterator bid(const blocked_index_type& offset)
    {
        return (m_bids.begin() +
                static_cast<typename bids_container_type::size_type>
                (offset.get_block2() * PageSize + offset.get_block1()));
    }
    const_bids_container_iterator bid(const size_type& offset) const
    {
        return (m_bids.begin() +
                static_cast<typename bids_container_type::size_type>
                (offset / block_type::size));
    }
    const_bids_container_iterator bid(const blocked_index_type& offset) const
    {
        return (m_bids.begin() +
                static_cast<typename bids_container_type::size_type>
                (offset.get_block2() * PageSize + offset.get_block1()));
    }

    void read_page(int_type page_no, int_type cache_slot) const
    {
        assert(page_no < (int_type)m_page_status.size());
        if (m_page_status[page_no] == uninitialized)
            return;
        STXXL_VERBOSE_VECTOR("read_page(): page_no=" << page_no << " cache_slot=" << cache_slot);
        request_ptr* reqs = new request_ptr[page_size];
        int_type block_no = page_no * page_size;
        int_type last_block = STXXL_MIN(block_no + page_size, int_type(m_bids.size()));
        int_type i = cache_slot * page_size, j = 0;
        for ( ; block_no < last_block; ++block_no, ++i, ++j)
        {
            reqs[j] = (*m_cache)[i].read(m_bids[block_no]);
        }
        assert(last_block - page_no * page_size > 0);
        wait_all(reqs, last_block - page_no * page_size);
        delete[] reqs;
    }
    void write_page(int_type page_no, int_type cache_slot) const
    {
        assert(page_no < (int_type)m_page_status.size());
        if (!(m_page_status[page_no] & dirty))
            return;
        STXXL_VERBOSE_VECTOR("write_page(): page_no=" << page_no << " cache_slot=" << cache_slot);
        request_ptr* reqs = new request_ptr[page_size];
        int_type block_no = page_no * page_size;
        int_type last_block = STXXL_MIN(block_no + page_size, int_type(m_bids.size()));
        assert(block_no < last_block);
        int_type i = cache_slot * page_size, j = 0;
        for ( ; block_no < last_block; ++block_no, ++i, ++j)
        {
            reqs[j] = (*m_cache)[i].write(m_bids[block_no]);
        }
        m_page_status[page_no] = valid_on_disk;
        assert(last_block - page_no * page_size > 0);
        wait_all(reqs, last_block - page_no * page_size);
        delete[] reqs;
    }

    reference element(size_type offset)
    {
#ifdef STXXL_RANGE_CHECK
        assert(offset < (size_type)size());
#endif
        return element(blocked_index_type(offset));
    }

    reference element(const blocked_index_type& offset)
    {
#ifdef STXXL_RANGE_CHECK
        assert(offset.get_pos() < size());
#endif
        unsigned_type page_no = offset.get_block2();
        assert(page_no < m_page_to_slot.size());   // fails if offset is too large, out of bound access
        int_type cache_slot = m_page_to_slot[page_no];
        if (cache_slot < 0)                        // == on_disk
        {
            if (m_free_slots.empty())              // has to kick
            {
                int_type kicked_slot = m_pager.kick();
                m_pager.hit(kicked_slot);
                int_type old_page_no = m_slot_to_page[kicked_slot];
                m_page_to_slot[page_no] = kicked_slot;
                m_page_to_slot[old_page_no] = on_disk;
                m_slot_to_page[kicked_slot] = page_no;

                write_page(old_page_no, kicked_slot);
                read_page(page_no, kicked_slot);

                m_page_status[page_no] = dirty;

                return (*m_cache)[kicked_slot * page_size + offset.get_block1()][offset.get_offset()];
            }
            else
            {
                int_type free_slot = m_free_slots.front();
                m_free_slots.pop();
                m_pager.hit(free_slot);
                m_page_to_slot[page_no] = free_slot;
                m_slot_to_page[free_slot] = page_no;

                read_page(page_no, free_slot);

                m_page_status[page_no] = dirty;

                return (*m_cache)[free_slot * page_size + offset.get_block1()][offset.get_offset()];
            }
        }
        else
        {
            m_page_status[page_no] = dirty;
            m_pager.hit(cache_slot);
            return (*m_cache)[cache_slot * page_size + offset.get_block1()][offset.get_offset()];
        }
    }

    // don't forget to first flush() the vector's cache before updating pages externally
    void page_externally_updated(unsigned_type page_no) const
    {
        // fails if offset is too large, out of bound access
        assert(page_no < m_page_status.size());
        // "A dirty page has been marked as newly initialized. The page content will be lost."
        assert(!(m_page_status[page_no] & dirty));
        if (m_page_to_slot[page_no] != on_disk) {
            // remove page from cache
            m_free_slots.push(m_page_to_slot[page_no]);
            m_page_to_slot[page_no] = on_disk;
            STXXL_VERBOSE_VECTOR("page_externally_updated(): page_no=" << page_no << " flushed from cache.");
        }
        else {
            STXXL_VERBOSE_VECTOR("page_externally_updated(): page_no=" << page_no << " no need to flush.");
        }
        m_page_status[page_no] = valid_on_disk;
    }

    void block_externally_updated(size_type offset) const
    {
        page_externally_updated(
            (unsigned_type)(offset / (block_type::size * page_size))
            );
    }

    void block_externally_updated(const blocked_index_type& offset) const
    {
        page_externally_updated(offset.get_block2());
    }

    const_reference const_element(size_type offset) const
    {
        return const_element(blocked_index_type(offset));
    }

    const_reference const_element(const blocked_index_type& offset) const
    {
        unsigned_type page_no = offset.get_block2();
        assert(page_no < m_page_to_slot.size());   // fails if offset is too large, out of bound access
        int_type cache_slot = m_page_to_slot[page_no];
        if (cache_slot < 0)                        // == on_disk
        {
            if (m_free_slots.empty())              // has to kick
            {
                int_type kicked_slot = m_pager.kick();
                m_pager.hit(kicked_slot);
                int_type old_page_no = m_slot_to_page[kicked_slot];
                m_page_to_slot[page_no] = kicked_slot;
                m_page_to_slot[old_page_no] = on_disk;
                m_slot_to_page[kicked_slot] = page_no;

                write_page(old_page_no, kicked_slot);
                read_page(page_no, kicked_slot);

                return (*m_cache)[kicked_slot * page_size + offset.get_block1()][offset.get_offset()];
            }
            else
            {
                int_type free_slot = m_free_slots.front();
                m_free_slots.pop();
                m_pager.hit(free_slot);
                m_page_to_slot[page_no] = free_slot;
                m_slot_to_page[free_slot] = page_no;

                read_page(page_no, free_slot);

                return (*m_cache)[free_slot * page_size + offset.get_block1()][offset.get_offset()];
            }
        }
        else
        {
            m_pager.hit(cache_slot);
            return (*m_cache)[cache_slot * page_size + offset.get_block1()][offset.get_offset()];
        }
    }

    bool is_page_cached(const blocked_index_type& offset) const
    {
        unsigned_type page_no = offset.get_block2();
        assert(page_no < m_page_to_slot.size());   // fails if offset is too large, out of bound access
        int_type cache_slot = m_page_to_slot[page_no];
        return (cache_slot >= 0);                  // on_disk == -1
    }
};

template <
    typename ValueType,
    unsigned PageSize,
    typename PagerType,
    unsigned BlockSize,
    typename AllocStr,
    typename SizeType>
inline bool operator == (stxxl::vector<ValueType, PageSize, PagerType, BlockSize,
                                       AllocStr, SizeType>& a,
                         stxxl::vector<ValueType, PageSize, PagerType, BlockSize,
                                       AllocStr, SizeType>& b)
{
    return a.size() == b.size() && std::equal(a.begin(), a.end(), b.begin());
}

template <
    typename ValueType,
    unsigned PageSize,
    typename PagerType,
    unsigned BlockSize,
    typename AllocStr,
    typename SizeType>
inline bool operator != (stxxl::vector<ValueType, PageSize, PagerType, BlockSize,
                                       AllocStr, SizeType>& a,
                         stxxl::vector<ValueType, PageSize, PagerType, BlockSize,
                                       AllocStr, SizeType>& b)
{
    return !(a == b);
}

template <
    typename ValueType,
    unsigned PageSize,
    typename PagerType,
    unsigned BlockSize,
    typename AllocStr,
    typename SizeType>
inline bool operator < (stxxl::vector<ValueType, PageSize, PagerType, BlockSize,
                                      AllocStr, SizeType>& a,
                        stxxl::vector<ValueType, PageSize, PagerType, BlockSize,
                                      AllocStr, SizeType>& b)
{
    return std::lexicographical_compare(a.begin(), a.end(), b.begin(), b.end());
}

template <
    typename ValueType,
    unsigned PageSize,
    typename PagerType,
    unsigned BlockSize,
    typename AllocStr,
    typename SizeType>
inline bool operator > (stxxl::vector<ValueType, PageSize, PagerType, BlockSize,
                                      AllocStr, SizeType>& a,
                        stxxl::vector<ValueType, PageSize, PagerType, BlockSize,
                                      AllocStr, SizeType>& b)
{
    return b < a;
}

template <
    typename ValueType,
    unsigned PageSize,
    typename PagerType,
    unsigned BlockSize,
    typename AllocStr,
    typename SizeType>
inline bool operator <= (stxxl::vector<ValueType, PageSize, PagerType, BlockSize,
                                       AllocStr, SizeType>& a,
                         stxxl::vector<ValueType, PageSize, PagerType, BlockSize,
                                       AllocStr, SizeType>& b)
{
    return !(b < a);
}

template <
    typename ValueType,
    unsigned PageSize,
    typename PagerType,
    unsigned BlockSize,
    typename AllocStr,
    typename SizeType>
inline bool operator >= (stxxl::vector<ValueType, PageSize, PagerType, BlockSize,
                                       AllocStr, SizeType>& a,
                         stxxl::vector<ValueType, PageSize, PagerType, BlockSize,
                                       AllocStr, SizeType>& b)
{
    return !(a < b);
}

////////////////////////////////////////////////////////////////////////////

// specialization for stxxl::vector, to use only const_iterators
template <typename ValueType, typename AllocStr, typename SizeType, typename DiffType,
          unsigned BlockSize, typename PagerType, unsigned PageSize>
bool is_sorted(
    stxxl::vector_iterator<ValueType, AllocStr, SizeType, DiffType, BlockSize, PagerType, PageSize> first,
    stxxl::vector_iterator<ValueType, AllocStr, SizeType, DiffType, BlockSize, PagerType, PageSize> last)
{
    return stxxl::is_sorted(
        stxxl::const_vector_iterator<ValueType, AllocStr, SizeType, DiffType, BlockSize, PagerType, PageSize>(first),
        stxxl::const_vector_iterator<ValueType, AllocStr, SizeType, DiffType, BlockSize, PagerType, PageSize>(last));
}

template <typename ValueType, typename AllocStr, typename SizeType, typename DiffType,
          unsigned BlockSize, typename PagerType, unsigned PageSize, typename StrictWeakOrdering>
bool is_sorted(
    stxxl::vector_iterator<ValueType, AllocStr, SizeType, DiffType, BlockSize, PagerType, PageSize> first,
    stxxl::vector_iterator<ValueType, AllocStr, SizeType, DiffType, BlockSize, PagerType, PageSize> last,
    StrictWeakOrdering comp)
{
    return stxxl::is_sorted(
        stxxl::const_vector_iterator<ValueType, AllocStr, SizeType, DiffType, BlockSize, PagerType, PageSize>(first),
        stxxl::const_vector_iterator<ValueType, AllocStr, SizeType, DiffType, BlockSize, PagerType, PageSize>(last),
        comp);
}

////////////////////////////////////////////////////////////////////////////

template <typename VectorBufReaderType>
class vector_bufreader_iterator;

/*!
 * Buffered sequential reader from a vector using overlapped I/O.
 *
 * This buffered reader can be used to read a large sequential region of a
 * vector using overlapped I/O. The object is created from an iterator range,
 * which can then be read to using operator<<(), or with operator*() and
 * operator++().
 *
 * The interface also fulfills all requirements of a stream. Actually most of
 * the code is identical to stream::vector_iterator2stream.
 *
 * Note that this buffered reader is inefficient for reading small ranges. This
 * is intentional, as one can just use operator[] on the vector for that.
 *
 * See \ref tutorial_vector_buf
 */
template <typename VectorIterator>
class vector_bufreader : public noncopyable
{
public:
    //! template parameter: the vector iterator type
    typedef VectorIterator vector_iterator;

    //! value type of the output vector
    typedef typename vector_iterator::value_type value_type;

    //! block type used in the vector
    typedef typename vector_iterator::block_type block_type;

    //! type of the input vector
    typedef typename vector_iterator::vector_type vector_type;

    //! block identifier iterator of the vector
    typedef typename vector_iterator::bids_container_iterator bids_container_iterator;

    //! construct output buffered stream used for overlapped reading
    typedef buf_istream<block_type, bids_container_iterator> buf_istream_type;

    //! construct an iterator for vector_bufreader (for C++11 range-based for loop)
    typedef vector_bufreader_iterator<vector_bufreader> bufreader_iterator;

    //! size of remaining data
    typedef typename vector_type::size_type size_type;

protected:
    //! iterator to the beginning of the range.
    vector_iterator m_begin;

    //! internal "current" iterator into the vector.
    vector_iterator m_iter;

    //! iterator to the end of the range.
    vector_iterator m_end;

    //! buffered input stream used to overlapped I/O.
    buf_istream_type* m_bufin;

    //! number of blocks to use as buffers.
    unsigned_type m_nbuffers;

    //! allow vector_bufreader_iterator to check m_iter against its current value
    friend class vector_bufreader_iterator<vector_bufreader>;

public:
    //! Create overlapped reader for the given iterator range.
    //! \param begin iterator to position were to start reading in vector
    //! \param end iterator to position were to end reading in vector
    //! \param nbuffers number of buffers used for overlapped I/O (>= 2*D recommended)
    vector_bufreader(vector_iterator begin, vector_iterator end,
                     unsigned_type nbuffers = 0)
        : m_begin(begin), m_end(end),
          m_bufin(NULL),
          m_nbuffers(nbuffers)
    {
        m_begin.flush(); // flush container

        if (m_nbuffers == 0)
            m_nbuffers = 2 * config::get_instance()->disks_number();

        rewind();
    }

    //! Create overlapped reader for the whole vector's content.
    //! \param vec vector to read
    //! \param nbuffers number of buffers used for overlapped I/O (>= 2*D recommended)
    vector_bufreader(const vector_type& vec, unsigned_type nbuffers = 0)
        : m_begin(vec.begin()), m_end(vec.end()),
          m_bufin(NULL),
          m_nbuffers(nbuffers)
    {
        m_begin.flush(); // flush container

        if (m_nbuffers == 0)
            m_nbuffers = 2 * config::get_instance()->disks_number();

        rewind();
    }

    //! Rewind stream back to begin. Note that this recreates the buffered
    //! reader and is thus not cheap.
    void rewind()
    {
        m_iter = m_begin;
        if (empty()) return;

        if (m_bufin) delete m_bufin;

        // find last bid to read
        bids_container_iterator end_bid = m_end.bid() + (m_end.block_offset() ? 1 : 0);

        // construct buffered istream for range
        m_bufin = new buf_istream_type(m_begin.bid(), end_bid, m_nbuffers);

        // skip the beginning of the block, up to real beginning
        vector_iterator curr = m_begin - m_begin.block_offset();

        for ( ; curr != m_begin; ++curr)
            ++(*m_bufin);
    }

    //! Finish reading and free buffered reader.
    ~vector_bufreader()
    {
        if (m_bufin) delete m_bufin;
    }

    //! Return constant reference to current item
    const value_type& operator * () const
    {
        return *(*m_bufin);
    }

    //! Return constant pointer to current item
    const value_type* operator -> () const
    {
        return &(*(*m_bufin));
    }

    //! Advance to next item (asserts if !empty()).
    vector_bufreader& operator ++ ()
    {
        assert(!empty());
        ++m_iter;
        ++(*m_bufin);

        if (UNLIKELY(empty())) {
            delete m_bufin;
            m_bufin = NULL;
        }

        return *this;
    }

    //! Read current item into variable and advance to next one.
    vector_bufreader& operator >> (value_type& v)
    {
        v = operator * ();
        operator ++ ();

        return *this;
    }

    //! Return remaining size.
    size_type size() const
    {
        assert(m_begin <= m_iter && m_iter <= m_end);
        return (size_type)(m_end - m_iter);
    }

    //! Returns true once the whole range has been read.
    bool empty() const
    {
        return (m_iter == m_end);
    }

    //! Return vector_bufreader_iterator for C++11 range-based for loop
    bufreader_iterator begin()
    {
        return bufreader_iterator(*this, m_begin);
    }

    //! Return vector_bufreader_iterator for C++11 range-based for loop
    bufreader_iterator end()
    {
        return bufreader_iterator(*this, m_end);
    }
};

////////////////////////////////////////////////////////////////////////////

/*!
 * Adapter for vector_bufreader to match iterator requirements of C++11
 * range-based loop construct.
 *
 * Since vector_bufreader itself points to only one specific item, this
 * iterator is merely a counter facade. The functions operator*() and
 * operator++() must only be called when it is in _sync_ with the bufreader
 * object. This is generally only the case for an iterator constructed with
 * begin() and then advanced with operator++(). The class checks this using
 * asserts(), the operators will fail if used wrong.
 *
 * See \ref tutorial_vector_buf
 */
template <typename VectorBufReaderType>
class vector_bufreader_iterator
{
public:
    //! The underlying buffered reader type
    typedef VectorBufReaderType vector_bufreader_type;

    //! Value type of vector
    typedef typename vector_bufreader_type::value_type value_type;

    //! Use vector_iterator to reference a point in the vector.
    typedef typename vector_bufreader_type::vector_iterator vector_iterator;

protected:
    //! Buffered reader used to access elements in vector
    vector_bufreader_type& m_bufreader;

    //! Use vector_iterator to reference a point in the vector.
    vector_iterator m_iter;

public:
    //! Construct iterator using vector_iterator
    vector_bufreader_iterator(vector_bufreader_type& bufreader, const vector_iterator& iter)
        : m_bufreader(bufreader), m_iter(iter)
    { }

    //! Return constant reference to current item
    const value_type& operator * () const
    {
        assert(m_bufreader.m_iter == m_iter);
        return m_bufreader.operator * ();
    }

    //! Return constant pointer to current item
    const value_type* operator -> () const
    {
        assert(m_bufreader.m_iter == m_iter);
        return m_bufreader.operator -> ();
    }

    //! Make bufreader advance to next item (asserts if !empty() or if iterator
    //! does not point to current).
    vector_bufreader_iterator& operator ++ ()
    {
        assert(m_bufreader.m_iter == m_iter);
        m_bufreader.operator ++ ();
        m_iter++;
        return *this;
    }

    //! Equality comparison operator
    bool operator == (const vector_bufreader_iterator& vbi) const
    {
        assert(&m_bufreader == &vbi.m_bufreader);
        return (m_iter == vbi.m_iter);
    }

    //! Inequality comparison operator
    bool operator != (const vector_bufreader_iterator& vbi) const
    {
        assert(&m_bufreader == &vbi.m_bufreader);
        return (m_iter != vbi.m_iter);
    }
};

////////////////////////////////////////////////////////////////////////////

/*!
 * Buffered sequential reverse reader from a vector using overlapped I/O.
 *
 * This buffered reader can be used to read a large sequential region of a
 * vector _in_reverse_ using overlapped I/O. The object is created from an
 * iterator range, which can then be read to using operator<<(), or with
 * operator*() and operator++(), where ++ actually goes to the preceding
 * element.
 *
 * The interface also fulfills all requirements of a stream. Actually most of
 * the code is identical to stream::vector_iterator2stream.
 *
 * Note that this buffered reader is inefficient for reading small ranges. This
 * is intentional, as one can just use operator[] on the vector for that.
 *
 * See \ref tutorial_vector_buf
 */
template <typename VectorIterator>
class vector_bufreader_reverse : public noncopyable
{
public:
    //! template parameter: the vector iterator type
    typedef VectorIterator vector_iterator;

    //! value type of the output vector
    typedef typename vector_iterator::value_type value_type;

    //! block type used in the vector
    typedef typename vector_iterator::block_type block_type;

    //! type of the input vector
    typedef typename vector_iterator::vector_type vector_type;

    //! block identifier iterator of the vector
    typedef typename vector_iterator::bids_container_iterator bids_container_iterator;

    //! construct output buffered stream used for overlapped reading
    typedef buf_istream_reverse<block_type, bids_container_iterator> buf_istream_type;

    //! size of remaining data
    typedef typename vector_type::size_type size_type;

protected:
    //! iterator to the beginning of the range.
    vector_iterator m_begin;

    //! internal "current" iterator into the vector.
    vector_iterator m_iter;

    //! iterator to the end of the range.
    vector_iterator m_end;

    //! buffered input stream used to overlapped I/O.
    buf_istream_type* m_bufin;

    //! number of blocks to use as buffers.
    unsigned_type m_nbuffers;

public:
    //! Create overlapped reader for the given iterator range.
    //! \param begin iterator to position were to start reading in vector
    //! \param end iterator to position were to end reading in vector
    //! \param nbuffers number of buffers used for overlapped I/O (>= 2*D recommended)
    vector_bufreader_reverse(vector_iterator begin, vector_iterator end,
                             unsigned_type nbuffers = 0)
        : m_begin(begin), m_end(end),
          m_bufin(NULL),
          m_nbuffers(nbuffers)
    {
        m_begin.flush(); // flush container

        if (m_nbuffers == 0)
            m_nbuffers = 2 * config::get_instance()->disks_number();

        rewind();
    }

    //! Create overlapped reader for the whole vector's content.
    //! \param vec vector to read
    //! \param nbuffers number of buffers used for overlapped I/O (>= 2*D recommended)
    vector_bufreader_reverse(const vector_type& vec, unsigned_type nbuffers = 0)
        : m_begin(vec.begin()), m_end(vec.end()),
          m_bufin(NULL),
          m_nbuffers(nbuffers)
    {
        m_begin.flush(); // flush container

        if (m_nbuffers == 0)
            m_nbuffers = 2 * config::get_instance()->disks_number();

        rewind();
    }

    //! Rewind stream back to begin. Note that this recreates the buffered
    //! reader and is thus not cheap.
    void rewind()
    {
        m_iter = m_end;
        if (empty()) return;

        if (m_bufin) delete m_bufin;

        // find last bid to read
        bids_container_iterator end_bid = m_end.bid() + (m_end.block_offset() ? 1 : 0);

        // construct buffered istream_reverse for range
        m_bufin = new buf_istream_type(m_begin.bid(), end_bid, m_nbuffers);

        // skip to beginning of reverse sequence.
        stxxl::int_type endoff = m_end.block_offset();
        if (endoff == 0) {
            // nothing to skip
        }
        else {
            // else, let ifstream_reverse skip last elements at end of block,
            // up to real end
            for ( ; endoff != block_type::size; endoff++)
                ++(*m_bufin);
        }
    }

    //! Finish reading and free buffered reader.
    ~vector_bufreader_reverse()
    {
        if (m_bufin) delete m_bufin;
    }

    //! Return constant reference to current item
    const value_type& operator * () const
    {
        return *(*m_bufin);
    }

    //! Return constant pointer to current item
    const value_type* operator -> () const
    {
        return &(*(*m_bufin));
    }

    //! Advance to next item (asserts if !empty()).
    vector_bufreader_reverse& operator ++ ()
    {
        assert(!empty());
        --m_iter;
        ++(*m_bufin);

        if (UNLIKELY(empty())) {
            delete m_bufin;
            m_bufin = NULL;
        }

        return *this;
    }

    //! Read current item into variable and advance to next one.
    vector_bufreader_reverse& operator >> (value_type& v)
    {
        v = operator * ();
        operator ++ ();

        return *this;
    }

    //! Return remaining size.
    size_type size() const
    {
        assert(m_begin <= m_iter && m_iter <= m_end);
        return (size_type)(m_iter - m_begin);
    }

    //! Returns true once the whole range has been read.
    bool empty() const
    {
        return (m_iter == m_begin);
    }
};

////////////////////////////////////////////////////////////////////////////

/*!
 * Buffered sequential writer to a vector using overlapped I/O.
 *
 * This buffered writer can be used to write a large sequential region of a
 * vector using overlapped I/O. The object is created from an iterator range,
 * which can then be written to using operator << (), or with operator * () and
 * operator ++ ().
 *
 * The buffered writer is given one iterator in the constructor. When writing,
 * this iterator advances in the vector and will \b enlarge the vector once it
 * reaches the end(). The vector size is doubled each time; nevertheless, it is
 * better to preinitialize the vector's size using stxxl::vector::resize().
 *
 * See \ref tutorial_vector_buf
 */
template <typename VectorIterator>
class vector_bufwriter : public noncopyable
{
public:
    //! template parameter: the vector iterator type
    typedef VectorIterator iterator;

    //! type of the output vector
    typedef typename iterator::vector_type vector_type;

    //! value type of the output vector
    typedef typename iterator::value_type value_type;

    //! block type used in the vector
    typedef typename iterator::block_type block_type;

    //! block identifier iterator of the vector
    typedef typename iterator::bids_container_iterator bids_container_iterator;

    //! iterator type of vector
    typedef typename iterator::iterator vector_iterator;
    typedef typename iterator::const_iterator vector_const_iterator;

    //! construct output buffered stream used for overlapped writing
    typedef buf_ostream<block_type, bids_container_iterator> buf_ostream_type;

protected:
    //! internal iterator into the vector.
    vector_iterator m_iter;

    //! iterator to the current end of the vector.
    vector_const_iterator m_end;

    //! boolean whether the vector was grown, will shorten at finish().
    bool m_grown;

    //! iterator into vector of the last block accessed (used to issue updates
    //! when the block is switched).
    vector_const_iterator m_prevblk;

    //! buffered output stream used to overlapped I/O.
    buf_ostream_type* m_bufout;

    //! number of blocks to use as buffers.
    unsigned_type m_nbuffers;

public:
    //! Create overlapped writer beginning at the given iterator.
    //! \param begin iterator to position were to start writing in vector
    //! \param nbuffers number of buffers used for overlapped I/O (>= 2D recommended)
    vector_bufwriter(vector_iterator begin,
                     unsigned_type nbuffers = 0)
        : m_iter(begin),
          m_end(m_iter.parent_vector()->end()),
          m_grown(false),
          m_bufout(NULL),
          m_nbuffers(nbuffers)
    {
        if (m_nbuffers == 0)
            m_nbuffers = 2 * config::get_instance()->disks_number();

        assert(m_iter <= m_end);
    }

    //! Create overlapped writer for the vector's beginning
    //! \param vec vector to write
    //! \param nbuffers number of buffers used for overlapped I/O (>= 2D recommended)
    vector_bufwriter(vector_type& vec,
                     unsigned_type nbuffers = 0)
        : m_iter(vec.begin()),
          m_end(m_iter.parent_vector()->end()),
          m_grown(false),
          m_bufout(NULL),
          m_nbuffers(nbuffers)
    {
        if (m_nbuffers == 0)
            m_nbuffers = 2 * config::get_instance()->disks_number();

        assert(m_iter <= m_end);
    }

    //! Finish writing and flush output back to vector.
    ~vector_bufwriter()
    {
        finish();
    }

    //! Return mutable reference to item at the position of the internal
    //! iterator.
    value_type& operator * ()
    {
        if (UNLIKELY(m_iter == m_end))
        {
            // iterator points to end of vector -> double vector's size

            if (m_bufout) {
                // fixes issue with buf_ostream writing invalid blocks: when
                // buf_ostream::current_elem advances to next block, flush()
                // will write to block beyond bid().end.
                if (m_iter.block_offset() != 0)
                    m_bufout->flush(); // flushes overlap buffers

                delete m_bufout;
                m_bufout = NULL;

                if (m_iter.block_offset() != 0)
                    m_iter.block_externally_updated();
            }

            vector_type& v = *m_iter.parent_vector();
            if (v.size() < 2 * block_type::size) {
                v.resize(2 * block_type::size);
            }
            else {
                v.resize(2 * v.size());
            }
            m_end = v.end();
            m_grown = true;
        }

        assert(m_iter < m_end);

        if (UNLIKELY(m_bufout == NULL))
        {
            if (m_iter.block_offset() != 0)
            {
                // output position is not at the start of the block, we
                // continue to use the iterator initially passed to the
                // constructor.
                return *m_iter;
            }
            else
            {
                // output position is start of block: create buffered writer

                m_iter.flush(); // flush container

                // create buffered write stream for blocks
                m_bufout = new buf_ostream_type(m_iter.bid(), m_nbuffers);
                m_prevblk = m_iter;

                // drop through to normal output into buffered writer
            }
        }

        // if the pointer has finished a block, then we inform the vector that
        // this block has been updated.
        if (UNLIKELY(m_iter.block_offset() == 0)) {
            if (m_prevblk != m_iter) {
                m_prevblk.block_externally_updated();
                m_prevblk = m_iter;
            }
        }

        return m_bufout->operator * ();
    }

    //! Advance internal iterator.
    vector_bufwriter& operator ++ ()
    {
        // always advance internal iterator
        ++m_iter;

        // if buf_ostream active, advance that too
        if (LIKELY(m_bufout != NULL)) m_bufout->operator ++ ();

        return *this;
    }

    //! Write value to the current position and advance the internal iterator.
    vector_bufwriter& operator << (const value_type& v)
    {
        operator * () = v;
        operator ++ ();

        return *this;
    }

    //! Finish writing and flush output back to vector.
    void finish()
    {
        if (m_bufout)
        {
            // must finish the block started in the buffered writer: fill it with
            // the data in the vector
            vector_const_iterator const_out = m_iter;

            while (const_out.block_offset() != 0)
            {
                m_bufout->operator * () = *const_out;
                m_bufout->operator ++ ();
                ++const_out;
            }

            // inform the vector that the block has been updated.
            if (m_prevblk != m_iter) {
                m_prevblk.block_externally_updated();
                m_prevblk = m_iter;
            }

            delete m_bufout;
            m_bufout = NULL;
        }

        if (m_grown)
        {
            vector_type& v = *m_iter.parent_vector();
            v.resize(m_iter - v.begin());

            m_grown = false;
        }
    }
};

////////////////////////////////////////////////////////////////////////////

//! External vector type generator.
//!
//! \tparam ValueType element type of contained objects (POD with no references to internal memory)
//! \tparam PageSize number of blocks in a page, default: \b 4 (recommended >= D)
//! \tparam CachePages number of pages in cache, default: \b 8 (recommended >= 2)
//! \tparam BlockSize external block size \a B in bytes, default: <b>2 MiB</b>
//! \tparam AllocStr parallel disk allocation strategies: \c striping, RC, SR, or FR. default: \b RC.
//! \tparam Pager pager type: \c random or \c lru, default: \b lru.
//!
//! \warning Do not store references to the elements of an external vector. Such references
//! might be invalidated during any following access to elements of the vector
template <
    typename ValueType,
    unsigned PageSize = 4,
    unsigned CachePages = 8,
    unsigned BlockSize = STXXL_DEFAULT_BLOCK_SIZE(ValueType),
    typename AllocStr = STXXL_DEFAULT_ALLOC_STRATEGY,
    pager_type Pager = lru
    >
struct VECTOR_GENERATOR
{
    typedef typename IF<Pager == lru,
                        lru_pager<CachePages>, random_pager<CachePages> >::result PagerType;

    typedef vector<ValueType, PageSize, PagerType, BlockSize, AllocStr> result;
};

//! \}

STXXL_END_NAMESPACE

namespace std {

template <
    typename ValueType,
    unsigned PageSize,
    typename PagerType,
    unsigned BlockSize,
    typename AllocStr,
    typename SizeType>
void swap(stxxl::vector<ValueType, PageSize, PagerType, BlockSize, AllocStr, SizeType>& a,
          stxxl::vector<ValueType, PageSize, PagerType, BlockSize, AllocStr, SizeType>& b)
{
    a.swap(b);
}

} // namespace std

#endif // !STXXL_CONTAINERS_VECTOR_HEADER
// vim: et:ts=4:sw=4
