/***************************************************************************
 *  include/stxxl/bits/containers/matrix.h
 *
 *  Part of the STXXL. See http://stxxl.sourceforge.net
 *
 *  Copyright (C) 2010-2011 Raoul Steffen <R-Steffen@gmx.de>
 *
 *  Distributed under the Boost Software License, Version 1.0.
 *  (See accompanying file LICENSE_1_0.txt or copy at
 *  http://www.boost.org/LICENSE_1_0.txt)
 **************************************************************************/

#ifndef STXXL_CONTAINERS_MATRIX_HEADER
#define STXXL_CONTAINERS_MATRIX_HEADER

#include <stxxl/bits/containers/vector.h>
#include <stxxl/bits/common/counting_ptr.h>
#include <stxxl/bits/mng/block_scheduler.h>
#include <stxxl/bits/containers/matrix_arithmetic.h>

STXXL_BEGIN_NAMESPACE

//! \defgroup matrix matrix
//! Efficient external memory matrix operations
//! \ingroup stlcont
//! \{

/* index-variable naming convention:
 * [MODIFIER_][UNIT_]DIMENSION[_in_[MODIFIER_]ENVIRONMENT]
 *
 * e.g.:
 * block_row = number of row measured in rows consisting of blocks
 * element_row_in_block = number of row measured in rows consisting of elements in the (row of) block(s)
 *
 * size-variable naming convention:
 * [MODIFIER_][ENVIRONMENT_]DIMENSION[_in_UNITs]
 *
 * e.g.
 * height_in_blocks
 */

// forward declaration
template <typename ValueType, unsigned BlockSideLength>
class matrix;

//! external column-vector container for matrix multiplication
//! \tparam ValueType type of contained objects (POD with no references to internal memory)
template <typename ValueType>
class column_vector : public vector<ValueType>
{
public:
    typedef vector<ValueType> vector_type;
    typedef typename vector_type::size_type size_type;

    using vector_type::size;

    //! \param n number of elements
    column_vector(size_type n = 0)
        : vector_type(n) { }

    column_vector operator + (const column_vector& right) const
    {
        assert(size() == right.size());
        column_vector res(size());
        for (size_type i = 0; i < size(); ++i)
            res[i] = (*this)[i] + right[i];
        return res;
    }

    column_vector operator - (const column_vector& right) const
    {
        assert(size() == right.size());
        column_vector res(size());
        for (size_type i = 0; i < size(); ++i)
            res[i] = (*this)[i] - right[i];
        return res;
    }

    column_vector operator * (const ValueType scalar) const
    {
        column_vector res(size());
        for (size_type i = 0; i < size(); ++i)
            res[i] = (*this)[i] * scalar;
        return res;
    }

    column_vector& operator += (const column_vector& right)
    {
        assert(size() == right.size());
        for (size_type i = 0; i < size(); ++i)
            (*this)[i] += right[i];
        return *this;
    }

    column_vector& operator -= (const column_vector& right)
    {
        assert(size() == right.size());
        for (size_type i = 0; i < size(); ++i)
            (*this)[i] -= right[i];
        return *this;
    }

    column_vector& operator *= (const ValueType scalar)
    {
        for (size_type i = 0; i < size(); ++i)
            (*this)[i] *= scalar;
        return *this;
    }

    void set_zero()
    {
        for (typename vector_type::iterator it = vector_type::begin(); it != vector_type::end(); ++it)
            *it = 0;
    }
};

//! external row-vector container for matrix multiplication
//! \tparam ValueType type of contained objects (POD with no references to internal memory)
template <typename ValueType>
class row_vector : public vector<ValueType>
{
public:
    typedef vector<ValueType> vector_type;
    typedef typename vector_type::size_type size_type;

    using vector_type::size;

    //! \param n number of elements
    row_vector(size_type n = 0)
        : vector_type(n) { }

    row_vector operator + (const row_vector& right) const
    {
        assert(size() == right.size());
        row_vector res(size());
        for (size_type i = 0; i < size(); ++i)
            res[i] = (*this)[i] + right[i];
        return res;
    }

    row_vector operator - (const row_vector& right) const
    {
        assert(size() == right.size());
        row_vector res(size());
        for (size_type i = 0; i < size(); ++i)
            res[i] = (*this)[i] - right[i];
        return res;
    }

    row_vector operator * (const ValueType scalar) const
    {
        row_vector res(size());
        for (size_type i = 0; i < size(); ++i)
            res[i] = (*this)[i] * scalar;
        return res;
    }

    template <unsigned BlockSideLength>
    row_vector operator * (const matrix<ValueType, BlockSideLength>& right) const
    { return right.multiply_from_left(*this); }

    ValueType operator * (const column_vector<ValueType>& right) const
    {
        ValueType res = 0;
        for (size_type i = 0; i < size(); ++i)
            res += (*this)[i] * right[i];
        return res;
    }

    row_vector& operator += (const row_vector& right)
    {
        assert(size() == right.size());
        for (size_type i = 0; i < size(); ++i)
            (*this)[i] += right[i];
        return *this;
    }

    row_vector& operator -= (const row_vector& right)
    {
        assert(size() == right.size());
        for (size_type i = 0; i < size(); ++i)
            (*this)[i] -= right[i];
        return *this;
    }

    row_vector& operator *= (const ValueType scalar)
    {
        for (size_type i = 0; i < size(); ++i)
            (*this)[i] *= scalar;
        return *this;
    }

    void set_zero()
    {
        for (typename vector_type::iterator it = vector_type::begin(); it != vector_type::end(); ++it)
            *it = 0;
    }
};

//! Specialized swappable_block that interprets uninitialized as containing zeros.
//! \tparam ValueType type of contained objects (POD with no references to internal memory)
//! \tparam BlockSideLength side length of a matrix block
//!
//! When initializing, all values are set to zero.
template <typename ValueType, unsigned BlockSideLength>
class matrix_swappable_block : public swappable_block<ValueType, BlockSideLength* BlockSideLength>
{
public:
    typedef typename swappable_block<ValueType, BlockSideLength* BlockSideLength>::internal_block_type internal_block_type;

    using swappable_block<ValueType, BlockSideLength* BlockSideLength>::get_internal_block;

    void fill_default()
    {
        // get_internal_block checks acquired
        internal_block_type& data = get_internal_block();
        #if STXXL_PARALLEL
        #pragma omp parallel for
        #endif
        for (int_type row = 0; row < int_type(BlockSideLength); ++row)
            for (int_type col = 0; col < int_type(BlockSideLength); ++col)
                data[row * BlockSideLength + col] = 0;
    }
};

//! External container for a (sub)matrix. Not intended for direct use.
//! \tparam ValueType type of contained objects (POD with no references to internal memory)
//! \tparam BlockSideLength side length of a matrix block
//!
//! Stores blocks only, so all measures (height, width, row, col) are in blocks.
template <typename ValueType, unsigned BlockSideLength>
class swappable_block_matrix : public atomic_counted_object
{
public:
    typedef int_type size_type;
    typedef int_type elem_size_type;
    typedef block_scheduler<matrix_swappable_block<ValueType, BlockSideLength> > block_scheduler_type;
    typedef typename block_scheduler_type::swappable_block_identifier_type swappable_block_identifier_type;
    typedef std::vector<swappable_block_identifier_type> blocks_type;
    typedef matrix_local::matrix_operations<ValueType, BlockSideLength> Ops;

    block_scheduler_type& bs;

private:
    // assigning is not allowed
    swappable_block_matrix& operator = (const swappable_block_matrix& other);

protected:
    //! height of the matrix in blocks
    size_type height,
    //! width of the matrix in blocks
        width,
    //! height copied from supermatrix in blocks
        height_from_supermatrix,
    //! width copied from supermatrix in blocks
        width_from_supermatrix;
    //! the matrice's blocks in row-major
    blocks_type blocks;
    //! if the elements in each block are in col-major instead of row-major
    bool elements_in_blocks_transposed;

    //! get identifier of the block at (row, col)
    swappable_block_identifier_type & bl(const size_type row, const size_type col)
    { return blocks[row * width + col]; }

public:
    //! Create an empty swappable_block_matrix of given dimensions.
    swappable_block_matrix(block_scheduler_type& bs, const size_type height_in_blocks, const size_type width_in_blocks, const bool transposed = false)
        : bs(bs),
          height(height_in_blocks),
          width(width_in_blocks),
          height_from_supermatrix(0),
          width_from_supermatrix(0),
          blocks(height * width),
          elements_in_blocks_transposed(transposed)
    {
        for (size_type row = 0; row < height; ++row)
            for (size_type col = 0; col < width; ++col)
                bl(row, col) = bs.allocate_swappable_block();
    }

    //! Create swappable_block_matrix of given dimensions that
    //! represents the submatrix of supermatrix starting at (from_row_in_blocks, from_col_in_blocks).
    //!
    //! If supermatrix is not large enough, the submatrix is padded with empty blocks.
    //! The supermatrix must not be destructed or transposed before the submatrix is destructed.
    swappable_block_matrix(const swappable_block_matrix& supermatrix,
                           const size_type height_in_blocks, const size_type width_in_blocks,
                           const size_type from_row_in_blocks, const size_type from_col_in_blocks)
        : bs(supermatrix.bs),
          height(height_in_blocks),
          width(width_in_blocks),
          height_from_supermatrix(std::min(supermatrix.height - from_row_in_blocks, height)),
          width_from_supermatrix(std::min(supermatrix.width - from_col_in_blocks, width)),
          blocks(height * width),
          elements_in_blocks_transposed(supermatrix.elements_in_blocks_transposed)
    {
        for (size_type row = 0; row < height_from_supermatrix; ++row)
        {
            for (size_type col = 0; col < width_from_supermatrix; ++col)
                bl(row, col) = supermatrix.block(row + from_row_in_blocks, col + from_col_in_blocks);
            for (size_type col = width_from_supermatrix; col < width; ++col)
                bl(row, col) = bs.allocate_swappable_block();
        }
        for (size_type row = height_from_supermatrix; row < height; ++row)
            for (size_type col = 0; col < width; ++col)
                bl(row, col) = bs.allocate_swappable_block();
    }

    //! Create swappable_block_matrix that represents the combination matrix ul ur dl dr.
    //!
    //! The submatrices are assumed to be of fitting dimensions and equal transposition.
    //! The submatrices must not be destructed or transposed before the matrix is destructed.
    swappable_block_matrix(const swappable_block_matrix& ul, const swappable_block_matrix& ur,
                           const swappable_block_matrix& dl, const swappable_block_matrix& dr)
        : bs(ul.bs),
          height(ul.height + dl.height),
          width(ul.width + ur.width),
          height_from_supermatrix(height),
          width_from_supermatrix(width),
          blocks(height * width),
          elements_in_blocks_transposed(ul.elements_in_blocks_transposed)
    {
        for (size_type row = 0; row < ul.height; ++row)
        {
            for (size_type col = 0; col < ul.width; ++col)
                bl(row, col) = ul.block(row, col);
            for (size_type col = ul.width; col < width; ++col)
                bl(row, col) = ur.block(row, col - ul.width);
        }
        for (size_type row = ul.height; row < height; ++row)
        {
            for (size_type col = 0; col < ul.width; ++col)
                bl(row, col) = dl.block(row - ul.height, col);
            for (size_type col = ul.width; col < width; ++col)
                bl(row, col) = dr.block(row - ul.height, col - ul.width);
        }
    }

    swappable_block_matrix(const swappable_block_matrix& other)
        : atomic_counted_object(other),
          bs(other.bs),
          height(other.height),
          width(other.width),
          height_from_supermatrix(0),
          width_from_supermatrix(0),
          blocks(height * width),
          elements_in_blocks_transposed(false)
    {
        for (size_type row = 0; row < height; ++row)
            for (size_type col = 0; col < width; ++col)
                bl(row, col) = bs.allocate_swappable_block();
        // 0 + other is copying
        Ops::element_op(*this, other, typename Ops::addition());
    }

    ~swappable_block_matrix()
    {
        for (size_type row = 0; row < height_from_supermatrix; ++row)
        {
            for (size_type col = width_from_supermatrix; col < width; ++col)
                bs.free_swappable_block(bl(row, col));
        }
        for (size_type row = height_from_supermatrix; row < height; ++row)
            for (size_type col = 0; col < width; ++col)
                bs.free_swappable_block(bl(row, col));
    }

    static size_type block_index_from_elem(elem_size_type index)
    { return index / BlockSideLength; }

    static int_type elem_index_in_block_from_elem(elem_size_type index)
    { return index % BlockSideLength; }

    // regards transposed
    int_type elem_index_in_block_from_elem(elem_size_type row, elem_size_type col) const
    {
        return (is_transposed())
               ? row % BlockSideLength + col % BlockSideLength * BlockSideLength
               : row % BlockSideLength * BlockSideLength + col % BlockSideLength;
    }

    //! get identifier of the block at (row, col)
    const swappable_block_identifier_type & block(const size_type row, const size_type col) const
    { return blocks[row * width + col]; }

    //! get identifier of the block at (row, col)
    const swappable_block_identifier_type& operator () (const size_type row, const size_type col) const
    { return block(row, col); }

    const size_type & get_height() const
    { return height; }

    const size_type & get_width() const
    { return width; }

    //! if the elements inside the blocks are in transposed order i.e. column-major
    const bool & is_transposed() const
    { return elements_in_blocks_transposed; }

    void transpose()
    {
        // transpose matrix of blocks
        blocks_type bn(blocks.size());
        for (size_type row = 0; row < height; ++row)
            for (size_type col = 0; col < width; ++col)
                bn[col * height + row] = bl(row, col);
        bn.swap(blocks);
        // swap dimensions
        std::swap(height, width);
        std::swap(height_from_supermatrix, width_from_supermatrix);
        elements_in_blocks_transposed = ! elements_in_blocks_transposed;
    }

    void set_zero()
    {
        for (typename blocks_type::iterator it = blocks.begin(); it != blocks.end(); ++it)
            bs.deinitialize(*it);
    }
};

//! general iterator type that points to single elements inside a matrix
//! \tparam ValueType type of contained objects (POD with no references to internal memory)
//! \tparam BlockSideLength side length of a matrix block
template <typename ValueType, unsigned BlockSideLength>
class matrix_iterator
{
protected:
    typedef matrix<ValueType, BlockSideLength> matrix_type;
    typedef typename matrix_type::swappable_block_matrix_type swappable_block_matrix_type;
    typedef typename matrix_type::block_scheduler_type block_scheduler_type;
    typedef typename block_scheduler_type::internal_block_type internal_block_type;
    typedef typename matrix_type::elem_size_type elem_size_type;
    typedef typename matrix_type::block_size_type block_size_type;

    template <typename VT, unsigned BSL>
    friend class matrix;

    template <typename VT, unsigned BSL>
    friend class const_matrix_iterator;

    matrix_type* m;
    elem_size_type current_row,          // \ both indices == -1 <=> empty iterator
        current_col;                     // /
    block_size_type current_block_row,
        current_block_col;
    internal_block_type* current_iblock; // NULL if block is not acquired

    void acquire_current_iblock()
    {
        if (! current_iblock)
            current_iblock = &m->data->bs.acquire(m->data->block(current_block_row, current_block_col));
    }

    void release_current_iblock()
    {
        if (current_iblock)
        {
            m->data->bs.release(m->data->block(current_block_row, current_block_col), true);
            current_iblock = 0;
        }
    }

    //! create iterator pointing to given row and col
    matrix_iterator(matrix_type& matrix, const elem_size_type start_row, const elem_size_type start_col)
        : m(&matrix),
          current_row(start_row),
          current_col(start_col),
          current_block_row(m->data->block_index_from_elem(start_row)),
          current_block_col(m->data->block_index_from_elem(start_col)),
          current_iblock(0) { }

    //! create empty iterator
    matrix_iterator(matrix_type& matrix)
        : m(&matrix),
          current_row(-1), // empty iterator
          current_col(-1),
          current_block_row(-1),
          current_block_col(-1),
          current_iblock(0) { }

    void set_empty()
    {
        release_current_iblock();
        current_row = -1;
        current_col = -1;
        current_block_row = -1;
        current_block_col = -1;
    }

public:
    matrix_iterator(const matrix_iterator& other)
        : m(other.m),
          current_row(other.current_row),
          current_col(other.current_col),
          current_block_row(other.current_block_row),
          current_block_col(other.current_block_col),
          current_iblock(0)
    {
        if (other.current_iblock)
            acquire_current_iblock();
    }

    matrix_iterator& operator = (const matrix_iterator& other)
    {
        set_pos(other.current_row, other.current_col);
        m = other.m;
        if (other.current_iblock)
            acquire_current_iblock();
        return *this;
    }

    ~matrix_iterator()
    { release_current_iblock(); }

    void set_row(const elem_size_type new_row)
    {
        const block_size_type new_block_row = m->data->block_index_from_elem(new_row);
        if (new_block_row != current_block_row)
        {
            release_current_iblock();
            current_block_row = new_block_row;
        }
        current_row = new_row;
    }

    void set_col(const elem_size_type new_col)
    {
        const block_size_type new_block_col = m->data->block_index_from_elem(new_col);
        if (new_block_col != current_block_col)
        {
            release_current_iblock();
            current_block_col = new_block_col;
        }
        current_col = new_col;
    }

    void set_pos(const elem_size_type new_row, const elem_size_type new_col)
    {
        const block_size_type new_block_row = m->data->block_index_from_elem(new_row),
            new_block_col = m->data->block_index_from_elem(new_col);
        if (new_block_col != current_block_col || new_block_row != current_block_row)
        {
            release_current_iblock();
            current_block_row = new_block_row;
            current_block_col = new_block_col;
        }
        current_row = new_row;
        current_col = new_col;
    }

    void set_pos(const std::pair<elem_size_type, elem_size_type> new_pos)
    { set_pos(new_pos.first, new_pos.second); }

    const elem_size_type & get_row() const
    { return current_row; }

    const elem_size_type & get_col() const
    { return current_col; }

    std::pair<elem_size_type, elem_size_type> get_pos() const
    { return std::make_pair(current_row, current_col); }

    bool empty() const
    { return current_row == -1 && current_col == -1; }

    operator bool () const
    { return ! empty(); }

    bool operator == (const matrix_iterator& other) const
    {
        return current_row == other.current_row && current_col == other.current_col && m == other.m;
    }

    //! Returns reference access to the element referenced by the iterator.
    //! The reference is only valid so long as the iterator is not moved.
    ValueType& operator * ()
    {
        acquire_current_iblock();
        return (*current_iblock)[m->data->elem_index_in_block_from_elem(current_row, current_col)];
    }
};

//! row-major iterator that points to single elements inside a matrix
//! \tparam ValueType type of contained objects (POD with no references to internal memory)
//! \tparam BlockSideLength side length of a matrix block
template <typename ValueType, unsigned BlockSideLength>
class matrix_row_major_iterator : public matrix_iterator<ValueType, BlockSideLength>
{
protected:
    typedef matrix_iterator<ValueType, BlockSideLength> matrix_iterator_type;
    typedef typename matrix_iterator_type::matrix_type matrix_type;
    typedef typename matrix_iterator_type::elem_size_type elem_size_type;

    template <typename VT, unsigned BSL>
    friend class matrix;

    using matrix_iterator_type::m;
    using matrix_iterator_type::set_empty;

    //! create iterator pointing to given row and col
    matrix_row_major_iterator(matrix_type& matrix, const elem_size_type start_row, const elem_size_type start_col)
        : matrix_iterator_type(matrix, start_row, start_col) { }

    //! create empty iterator
    matrix_row_major_iterator(matrix_type& matrix)
        : matrix_iterator_type(matrix) { }

public:
    //! convert from matrix_iterator
    matrix_row_major_iterator(const matrix_iterator_type& matrix_iterator)
        : matrix_iterator_type(matrix_iterator) { }

    // Has to be not empty, else behavior is undefined.
    matrix_row_major_iterator& operator ++ ()
    {
        if (get_col() + 1 < m->get_width())
            // => not matrix_row_major_iterator the end of row, move right
            set_col(get_col() + 1);
        else if (get_row() + 1 < m->get_height())
            // => at end of row but not last row, move to beginning of next row
            set_pos(get_row() + 1, 0);
        else
            // => at end of matrix, set to empty-state
            set_empty();
        return *this;
    }

    // Has to be not empty, else behavior is undefined.
    matrix_row_major_iterator& operator -- ()
    {
        if (get_col() - 1 >= 0)
            // => not at the beginning of row, move left
            set_col(get_col() - 1);
        else if (get_row() - 1 >= 0)
            // => at beginning of row but not first row, move to end of previous row
            set_pos(get_row() - 1, m->get_width() - 1);
        else
            // => at beginning of matrix, set to empty-state
            set_empty();
        return *this;
    }

    using matrix_iterator_type::get_row;
    using matrix_iterator_type::get_col;
    using matrix_iterator_type::set_col;
    using matrix_iterator_type::set_pos;
};

//! column-major iterator that points to single elements inside a matrix
//! \tparam ValueType type of contained objects (POD with no references to internal memory)
//! \tparam BlockSideLength side length of a matrix block
template <typename ValueType, unsigned BlockSideLength>
class matrix_col_major_iterator : public matrix_iterator<ValueType, BlockSideLength>
{
protected:
    typedef matrix_iterator<ValueType, BlockSideLength> matrix_iterator_type;
    typedef typename matrix_iterator_type::matrix_type matrix_type;
    typedef typename matrix_iterator_type::elem_size_type elem_size_type;

    template <typename VT, unsigned BSL>
    friend class matrix;

    using matrix_iterator_type::m;
    using matrix_iterator_type::set_empty;

    //! create iterator pointing to given row and col
    matrix_col_major_iterator(matrix_type& matrix, const elem_size_type start_row, const elem_size_type start_col)
        : matrix_iterator_type(matrix, start_row, start_col) { }

    //! create empty iterator
    matrix_col_major_iterator(matrix_type& matrix)
        : matrix_iterator_type(matrix) { }

public:
    //! convert from matrix_iterator
    matrix_col_major_iterator(const matrix_iterator_type& matrix_iterator)
        : matrix_iterator_type(matrix_iterator) { }

    // Has to be not empty, else behavior is undefined.
    matrix_col_major_iterator& operator ++ ()
    {
        if (get_row() + 1 < m->get_height())
            // => not at the end of col, move down
            set_row(get_row() + 1);
        else if (get_col() + 1 < m->get_width())
            // => at end of col but not last col, move to beginning of next col
            set_pos(0, get_col() + 1);
        else
            // => at end of matrix, set to empty-state
            set_empty();
        return *this;
    }

    // Has to be not empty, else behavior is undefined.
    matrix_col_major_iterator& operator -- ()
    {
        if (get_row() - 1 >= 0)
            // => not at the beginning of col, move up
            set_row(get_row() - 1);
        else if (get_col() - 1 >= 0)
            // => at beginning of col but not first col, move to end of previous col
            set_pos(m->get_height() - 1, get_col() - 1);
        else
            // => at beginning of matrix, set to empty-state
            set_empty();
        return *this;
    }

    using matrix_iterator_type::get_row;
    using matrix_iterator_type::get_col;
    using matrix_iterator_type::set_row;
    using matrix_iterator_type::set_pos;
};

//! general const_iterator type that points to single elements inside a matrix
//! \tparam ValueType type of contained objects (POD with no references to internal memory)
//! \tparam BlockSideLength side length of a matrix block
template <typename ValueType, unsigned BlockSideLength>
class const_matrix_iterator
{
protected:
    typedef matrix<ValueType, BlockSideLength> matrix_type;
    typedef typename matrix_type::swappable_block_matrix_type swappable_block_matrix_type;
    typedef typename matrix_type::block_scheduler_type block_scheduler_type;
    typedef typename block_scheduler_type::internal_block_type internal_block_type;
    typedef typename matrix_type::elem_size_type elem_size_type;
    typedef typename matrix_type::block_size_type block_size_type;

    template <typename VT, unsigned BSL>
    friend class matrix;

    const matrix_type* m;
    elem_size_type current_row,          // \ both indices == -1 <=> empty iterator
        current_col;                     // /
    block_size_type current_block_row,
        current_block_col;
    internal_block_type* current_iblock; // NULL if block is not acquired

    void acquire_current_iblock()
    {
        if (! current_iblock)
            current_iblock = &m->data->bs.acquire(m->data->block(current_block_row, current_block_col));
    }

    void release_current_iblock()
    {
        if (current_iblock)
        {
            m->data->bs.release(m->data->block(current_block_row, current_block_col), false);
            current_iblock = 0;
        }
    }

    //! create iterator pointing to given row and col
    const_matrix_iterator(const matrix_type& matrix, const elem_size_type start_row, const elem_size_type start_col)
        : m(&matrix),
          current_row(start_row),
          current_col(start_col),
          current_block_row(m->data->block_index_from_elem(start_row)),
          current_block_col(m->data->block_index_from_elem(start_col)),
          current_iblock(0) { }

    //! create empty iterator
    const_matrix_iterator(const matrix_type& matrix)
        : m(&matrix),
          current_row(-1), // empty iterator
          current_col(-1),
          current_block_row(-1),
          current_block_col(-1),
          current_iblock(0) { }

    void set_empty()
    {
        release_current_iblock();
        current_row = -1;
        current_col = -1;
        current_block_row = -1;
        current_block_col = -1;
    }

public:
    const_matrix_iterator(const matrix_iterator<ValueType, BlockSideLength>& other)
        : m(other.m),
          current_row(other.current_row),
          current_col(other.current_col),
          current_block_row(other.current_block_row),
          current_block_col(other.current_block_col),
          current_iblock(0)
    {
        if (other.current_iblock)
            acquire_current_iblock();
    }

    const_matrix_iterator(const const_matrix_iterator& other)
        : m(other.m),
          current_row(other.current_row),
          current_col(other.current_col),
          current_block_row(other.current_block_row),
          current_block_col(other.current_block_col),
          current_iblock(0)
    {
        if (other.current_iblock)
            acquire_current_iblock();
    }

    const_matrix_iterator& operator = (const const_matrix_iterator& other)
    {
        set_pos(other.current_row, other.current_col);
        m = other.m;
        if (other.current_iblock)
            acquire_current_iblock();
        return *this;
    }

    ~const_matrix_iterator()
    { release_current_iblock(); }

    void set_row(const elem_size_type new_row)
    {
        const block_size_type new_block_row = m->data->block_index_from_elem(new_row);
        if (new_block_row != current_block_row)
        {
            release_current_iblock();
            current_block_row = new_block_row;
        }
        current_row = new_row;
    }

    void set_col(const elem_size_type new_col)
    {
        const block_size_type new_block_col = m->data->block_index_from_elem(new_col);
        if (new_block_col != current_block_col)
        {
            release_current_iblock();
            current_block_col = new_block_col;
        }
        current_col = new_col;
    }

    void set_pos(const elem_size_type new_row, const elem_size_type new_col)
    {
        const block_size_type new_block_row = m->data->block_index_from_elem(new_row),
            new_block_col = m->data->block_index_from_elem(new_col);
        if (new_block_col != current_block_col || new_block_row != current_block_row)
        {
            release_current_iblock();
            current_block_row = new_block_row;
            current_block_col = new_block_col;
        }
        current_row = new_row;
        current_col = new_col;
    }

    void set_pos(const std::pair<elem_size_type, elem_size_type> new_pos)
    { set_pos(new_pos.first, new_pos.second); }

    const elem_size_type & get_row() const
    { return current_row; }

    const elem_size_type & get_col() const
    { return current_col; }

    std::pair<elem_size_type, elem_size_type> get_pos() const
    { return std::make_pair(current_row, current_col); }

    bool empty() const
    { return current_row == -1 && current_col == -1; }

    operator bool () const
    { return ! empty(); }

    bool operator == (const const_matrix_iterator& other) const
    {
        return current_row == other.current_row && current_col == other.current_col && m == other.m;
    }

    //! Returns reference access to the element referenced by the iterator.
    //! The reference is only valid so long as the iterator is not moved.
    const ValueType& operator * ()
    {
        acquire_current_iblock();
        return (*current_iblock)[m->data->elem_index_in_block_from_elem(current_row, current_col)];
    }
};

//! row-major const_iterator that points to single elements inside a matrix
//! \tparam ValueType type of contained objects (POD with no references to internal memory)
//! \tparam BlockSideLength side length of a matrix block
template <typename ValueType, unsigned BlockSideLength>
class const_matrix_row_major_iterator : public const_matrix_iterator<ValueType, BlockSideLength>
{
protected:
    typedef const_matrix_iterator<ValueType, BlockSideLength> const_matrix_iterator_type;
    typedef typename const_matrix_iterator_type::matrix_type matrix_type;
    typedef typename const_matrix_iterator_type::elem_size_type elem_size_type;

    template <typename VT, unsigned BSL>
    friend class matrix;

    using const_matrix_iterator_type::m;
    using const_matrix_iterator_type::set_empty;

    //! create iterator pointing to given row and col
    const_matrix_row_major_iterator(const matrix_type& matrix, const elem_size_type start_row, const elem_size_type start_col)
        : const_matrix_iterator_type(matrix, start_row, start_col) { }

    //! create empty iterator
    const_matrix_row_major_iterator(const matrix_type& matrix)
        : const_matrix_iterator_type(matrix) { }

public:
    //! convert from matrix_iterator
    const_matrix_row_major_iterator(const const_matrix_row_major_iterator& matrix_iterator)
        : const_matrix_iterator_type(matrix_iterator) { }

    //! convert from matrix_iterator
    const_matrix_row_major_iterator(const const_matrix_iterator_type& matrix_iterator)
        : const_matrix_iterator_type(matrix_iterator) { }

    // Has to be not empty, else behavior is undefined.
    const_matrix_row_major_iterator& operator ++ ()
    {
        if (get_col() + 1 < m->get_width())
            // => not matrix_row_major_iterator the end of row, move right
            set_col(get_col() + 1);
        else if (get_row() + 1 < m->get_height())
            // => at end of row but not last row, move to beginning of next row
            set_pos(get_row() + 1, 0);
        else
            // => at end of matrix, set to empty-state
            set_empty();
        return *this;
    }

    // Has to be not empty, else behavior is undefined.
    const_matrix_row_major_iterator& operator -- ()
    {
        if (get_col() - 1 >= 0)
            // => not at the beginning of row, move left
            set_col(get_col() - 1);
        else if (get_row() - 1 >= 0)
            // => at beginning of row but not first row, move to end of previous row
            set_pos(get_row() - 1, m->get_width() - 1);
        else
            // => at beginning of matrix, set to empty-state
            set_empty();
        return *this;
    }

    using const_matrix_iterator_type::get_row;
    using const_matrix_iterator_type::get_col;
    using const_matrix_iterator_type::set_col;
    using const_matrix_iterator_type::set_pos;
};

//! column-major const_iterator that points to single elements inside a matrix
//! \tparam ValueType type of contained objects (POD with no references to internal memory)
//! \tparam BlockSideLength side length of a matrix block
template <typename ValueType, unsigned BlockSideLength>
class const_matrix_col_major_iterator : public const_matrix_iterator<ValueType, BlockSideLength>
{
protected:
    typedef const_matrix_iterator<ValueType, BlockSideLength> const_matrix_iterator_type;
    typedef typename const_matrix_iterator_type::matrix_type matrix_type;
    typedef typename const_matrix_iterator_type::elem_size_type elem_size_type;

    template <typename VT, unsigned BSL>
    friend class matrix;

    using const_matrix_iterator_type::m;
    using const_matrix_iterator_type::set_empty;

    //! create iterator pointing to given row and col
    const_matrix_col_major_iterator(const matrix_type& matrix, const elem_size_type start_row, const elem_size_type start_col)
        : const_matrix_iterator_type(matrix, start_row, start_col) { }

    //! create empty iterator
    const_matrix_col_major_iterator(const matrix_type& matrix)
        : const_matrix_iterator_type(matrix) { }

public:
    //! convert from matrix_iterator
    const_matrix_col_major_iterator(const matrix_iterator<ValueType, BlockSideLength>& matrix_iterator)
        : const_matrix_iterator_type(matrix_iterator) { }

    //! convert from matrix_iterator
    const_matrix_col_major_iterator(const const_matrix_iterator_type& matrix_iterator)
        : const_matrix_iterator_type(matrix_iterator) { }

    // Has to be not empty, else behavior is undefined.
    const_matrix_col_major_iterator& operator ++ ()
    {
        if (get_row() + 1 < m->get_height())
            // => not at the end of col, move down
            set_row(get_row() + 1);
        else if (get_col() + 1 < m->get_width())
            // => at end of col but not last col, move to beginning of next col
            set_pos(0, get_col() + 1);
        else
            // => at end of matrix, set to empty-state
            set_empty();
        return *this;
    }

    // Has to be not empty, else behavior is undefined.
    const_matrix_col_major_iterator& operator -- ()
    {
        if (get_row() - 1 >= 0)
            // => not at the beginning of col, move up
            set_row(get_row() - 1);
        else if (get_col() - 1 >= 0)
            // => at beginning of col but not first col, move to end of previous col
            set_pos(m->get_height() - 1, get_col() - 1);
        else
            // => at beginning of matrix, set to empty-state
            set_empty();
        return *this;
    }

    using const_matrix_iterator_type::get_row;
    using const_matrix_iterator_type::get_col;
    using const_matrix_iterator_type::set_row;
    using const_matrix_iterator_type::set_pos;
};

//! External matrix container. \n
//! <b> Introduction </b> to matrix container: see \ref tutorial_matrix tutorial. \n
//! <b> Design and Internals </b> of matrix container: see \ref design_matrix.
//!
//! \tparam ValueType type of contained objects (POD with no references to internal memory)
//! \tparam BlockSideLength side length of a matrix block
//!
//! Divides the matrix in square submatrices (blocks).
//! Blocks can be swapped individually to and from external memory.
//! They are only swapped if necessary to minimize I/O.
template <typename ValueType, unsigned BlockSideLength>
class matrix
{
protected:
    typedef matrix<ValueType, BlockSideLength> matrix_type;
    typedef swappable_block_matrix<ValueType, BlockSideLength> swappable_block_matrix_type;
    typedef counting_ptr<swappable_block_matrix_type> swappable_block_matrix_pointer_type;
    typedef typename swappable_block_matrix_type::block_scheduler_type block_scheduler_type;
    typedef typename swappable_block_matrix_type::size_type block_size_type;
    typedef typename swappable_block_matrix_type::elem_size_type elem_size_type;
    typedef matrix_local::matrix_operations<ValueType, BlockSideLength> Ops;
    typedef matrix_swappable_block<ValueType, BlockSideLength> swappable_block_type;

public:
    typedef matrix_iterator<ValueType, BlockSideLength> iterator;
    typedef const_matrix_iterator<ValueType, BlockSideLength> const_iterator;
    typedef matrix_row_major_iterator<ValueType, BlockSideLength> row_major_iterator;
    typedef matrix_col_major_iterator<ValueType, BlockSideLength> col_major_iterator;
    typedef const_matrix_row_major_iterator<ValueType, BlockSideLength> const_row_major_iterator;
    typedef const_matrix_col_major_iterator<ValueType, BlockSideLength> const_col_major_iterator;
    typedef column_vector<ValueType> column_vector_type;
    typedef row_vector<ValueType> row_vector_type;

protected:
    template <typename VT, unsigned BSL>
    friend class matrix_iterator;

    template <typename VT, unsigned BSL>
    friend class const_matrix_iterator;

    elem_size_type height, width;

    swappable_block_matrix_pointer_type data;

public:
    //! \name Constructors/Destructors
    //! \{

    //! Creates a new matrix of given dimensions. Elements' values are set to zero.
    //! \param bs block scheduler used
    //! \param height height of the created matrix
    //! \param width width of the created matrix
    matrix(block_scheduler_type& bs,
           const elem_size_type height, const elem_size_type width)
        : height(height),
          width(width),
          data(
              new swappable_block_matrix_type(
                  bs,
                  div_ceil(height, BlockSideLength),
                  div_ceil(width, BlockSideLength))
              )
    { }

    matrix(block_scheduler_type& bs,
           const column_vector_type& left, const row_vector_type& right)
        : height((elem_size_type)left.size()),
          width((elem_size_type)right.size()),
          data(
              new swappable_block_matrix_type(
                  bs,
                  div_ceil(height, BlockSideLength),
                  div_ceil(width, BlockSideLength))
              )
    { Ops::recursive_matrix_from_vectors(*data, left, right); }

    ~matrix() { }
    //! \}

    //! \name Capacity
    //! \{
    const elem_size_type & get_height() const
    { return height; }

    const elem_size_type & get_width() const
    { return width; }
    //! \}

    //! \name Iterators
    //! \{
    iterator begin()
    {
        data.unify();
        return iterator(*this, 0, 0);
    }
    const_iterator begin() const
    { return const_iterator(*this, 0, 0); }
    const_iterator cbegin() const
    { return const_iterator(*this, 0, 0); }
    iterator end()
    {
        data.unify();
        return iterator(*this);
    }
    const_iterator end() const
    { return const_iterator(*this); }
    const_iterator cend() const
    { return const_iterator(*this); }
    const_iterator operator () (const elem_size_type row, const elem_size_type col) const
    { return const_iterator(*this, row, col); }

    iterator operator () (const elem_size_type row, const elem_size_type col)
    {
        data.unify();
        return iterator(*this, row, col);
    }
    //! \}

    //! \name Modifiers
    //! \{
    void transpose()
    {
        data.unify();
        data->transpose();
        std::swap(height, width);
    }

    void set_zero()
    {
        if (data.unique())
            data->set_zero();
        else
            data = new swappable_block_matrix_type
                       (data->bs, div_ceil(height, BlockSideLength), div_ceil(width, BlockSideLength));
    }
    //! \}

    //! \name Operations
    //! \{
    matrix_type operator + (const matrix_type& right) const
    {
        assert(height == right.height && width == right.width);
        matrix_type res(data->bs, height, width);
        Ops::element_op(*res.data, *data, *right.data, typename Ops::addition()); // more efficient than copying this and then adding right
        return res;
    }

    matrix_type operator - (const matrix_type& right) const
    {
        assert(height == right.height && width == right.width);
        matrix_type res(data->bs, height, width);
        Ops::element_op(*res.data, *data, *right.data, typename Ops::subtraction()); // more efficient than copying this and then subtracting right
        return res;
    }

    matrix_type operator * (const matrix_type& right) const
    { return multiply(right); }

    matrix_type operator * (const ValueType scalar) const
    {
        matrix_type res(data->bs, height, width);
        Ops::element_op(*res.data, *data, typename Ops::scalar_multiplication(scalar));
        return res;
    }

    matrix_type& operator += (const matrix_type& right)
    {
        assert(height == right.height && width == right.width);
        data.unify();
        Ops::element_op(*data, *right.data, typename Ops::addition());
        return *this;
    }

    matrix_type& operator -= (const matrix_type& right)
    {
        assert(height == right.height && width == right.width);
        data.unify();
        Ops::element_op(*data, *right.data, typename Ops::subtraction());
        return *this;
    }

    matrix_type& operator *= (const matrix_type& right)
    { return *this = operator * (right); } // implicitly unifies by constructing a result-matrix

    matrix_type& operator *= (const ValueType scalar)
    {
        data.unify();
        Ops::element_op(*data, typename Ops::scalar_multiplication(scalar));
        return *this;
    }

    column_vector_type operator * (const column_vector_type& right) const
    {
        assert(elem_size_type(right.size()) == width);
        column_vector_type res(height);
        res.set_zero();
        Ops::recursive_matrix_col_vector_multiply_and_add(*data, right, res);
        return res;
    }

    row_vector_type multiply_from_left(const row_vector_type& left) const
    {
        assert(elem_size_type(left.size()) == height);
        row_vector_type res(width);
        res.set_zero();
        Ops::recursive_matrix_row_vector_multiply_and_add(left, *data, res);
        return res;
    }

    //! multiply with another matrix
    //! \param right matrix to multiply with
    //! \param multiplication_algorithm allows to choose the applied algorithm
    //! \param scheduling_algorithm  allows to choose the applied algorithm
    //!
    //! Available algorithms are: \n
    //!    0: naive_multiply_and_add (I/O inefficient, slow) \n
    //!    1: recursive_multiply_and_add (recommended, default, stable time and I/O complexity) \n
    //!    2: strassen_winograd_multiply_and_add (sometimes fast but unstable time and I/O complexity) \n
    //!    3: multi_level_strassen_winograd_multiply_and_add (sometimes fast but unstable time and I/O complexity) \n
    //!    4: strassen_winograd_multiply, optimized pre- and postadditions (sometimes fast but unstable time and I/O complexity) \n
    //!    5: strassen_winograd_multiply_and_add_interleaved, optimized preadditions (sometimes fast but unstable time and I/O complexity) \n
    //!    6: multi_level_strassen_winograd_multiply_and_add_block_grained (sometimes fast but unstable time and I/O complexity)
    matrix_type multiply(const matrix_type& right, const int_type multiplication_algorithm = 1, const int_type scheduling_algorithm = 2) const
    {
        assert(width == right.height);
        assert(&data->bs == &right.data->bs);
        matrix_type res(data->bs, height, right.width);

        if (scheduling_algorithm > 0)
        {
            // all offline algos need a simulation-run
            delete data->bs.switch_algorithm_to(
                new block_scheduler_algorithm_simulation<swappable_block_type>(data->bs)
                );
            switch (multiplication_algorithm)
            {
            case 0:
                Ops::naive_multiply_and_add(*data, *right.data, *res.data);
                break;
            case 1:
                Ops::recursive_multiply_and_add(*data, *right.data, *res.data);
                break;
            case 2:
                Ops::strassen_winograd_multiply_and_add(*data, *right.data, *res.data);
                break;
            case 3:
                Ops::multi_level_strassen_winograd_multiply_and_add(*data, *right.data, *res.data);
                break;
            case 4:
                Ops::strassen_winograd_multiply(*data, *right.data, *res.data);
                break;
            case 5:
                Ops::strassen_winograd_multiply_and_add_interleaved(*data, *right.data, *res.data);
                break;
            case 6:
                Ops::multi_level_strassen_winograd_multiply_and_add_block_grained(*data, *right.data, *res.data);
                break;
            default:
                STXXL_ERRMSG("invalid multiplication-algorithm number");
                break;
            }
        }
        switch (scheduling_algorithm)
        {
        case 0:
            delete data->bs.switch_algorithm_to(
                new block_scheduler_algorithm_online_lru<swappable_block_type>(data->bs)
                );
            break;
        case 1:
            delete data->bs.switch_algorithm_to(
                new block_scheduler_algorithm_offline_lfd<swappable_block_type>(data->bs)
                );
            break;
        case 2:
            delete data->bs.switch_algorithm_to(
                new block_scheduler_algorithm_offline_lru_prefetching<swappable_block_type>(data->bs)
                );
            break;
        default:
            STXXL_ERRMSG("invalid scheduling-algorithm number");
        }
        switch (multiplication_algorithm)
        {
        case 0:
            Ops::naive_multiply_and_add(*data, *right.data, *res.data);
            break;
        case 1:
            Ops::recursive_multiply_and_add(*data, *right.data, *res.data);
            break;
        case 2:
            Ops::strassen_winograd_multiply_and_add(*data, *right.data, *res.data);
            break;
        case 3:
            Ops::multi_level_strassen_winograd_multiply_and_add(*data, *right.data, *res.data);
            break;
        case 4:
            Ops::strassen_winograd_multiply(*data, *right.data, *res.data);
            break;
        case 5:
            Ops::strassen_winograd_multiply_and_add_interleaved(*data, *right.data, *res.data);
            break;
        case 6:
            Ops::multi_level_strassen_winograd_multiply_and_add_block_grained(*data, *right.data, *res.data);
            break;
        default:
            STXXL_ERRMSG("invalid multiplication-algorithm number");
            break;
        }
        delete data->bs.switch_algorithm_to(
            new block_scheduler_algorithm_online_lru<swappable_block_type>(data->bs)
            );
        return res;
    }

    //! Use internal memory multiplication. Designated for testing. May exceed memory limitations.
    matrix_type multiply_internal(const matrix_type& right, const int_type scheduling_algorithm = 2) const
    {
        assert(width == right.height);
        assert(&data->bs == &right.data->bs);
        matrix_type res(data->bs, height, right.width);

        if (scheduling_algorithm > 0)
        {
            // all offline algos need a simulation-run
            delete data->bs.switch_algorithm_to(
                new block_scheduler_algorithm_simulation<swappable_block_type>(data->bs)
                );
            multiply_internal(right, res);
        }
        switch (scheduling_algorithm)
        {
        case 0:
            delete data->bs.switch_algorithm_to(
                new block_scheduler_algorithm_online_lru<swappable_block_type>(data->bs)
                );
            break;
        case 1:
            delete data->bs.switch_algorithm_to(
                new block_scheduler_algorithm_offline_lfd<swappable_block_type>(data->bs)
                );
            break;
        case 2:
            delete data->bs.switch_algorithm_to(
                new block_scheduler_algorithm_offline_lru_prefetching<swappable_block_type>(data->bs)
                );
            break;
        default:
            STXXL_ERRMSG("invalid scheduling-algorithm number");
        }
        multiply_internal(right, res);
        delete data->bs.switch_algorithm_to(
            new block_scheduler_algorithm_online_lru<swappable_block_type>(data->bs)
            );
        return res;
    }
    //! \}

protected:
    void multiply_internal(const matrix_type& right, matrix_type& res) const
    {
        ValueType* A = new ValueType[height * width];
        ValueType* B = new ValueType[right.height * right.width];
        ValueType* C = new ValueType[res.height * res.width];
        ValueType* vit;
        vit = A;
        for (const_row_major_iterator mit = cbegin(); mit != cend(); ++mit, ++vit)
            *vit = *mit;
        vit = B;
        for (const_row_major_iterator mit = right.cbegin(); mit != right.cend(); ++mit, ++vit)
            *vit = *mit;
        if (! res.data->bs.is_simulating())
        {
#if STXXL_BLAS
            gemm_wrapper(height, width, res.width,
                         ValueType(1), false, A,
                         false, B,
                         ValueType(0), false, C);
#else
            assert(false /* internal multiplication is only available for testing with blas */);
#endif
        }
        vit = C;
        for (row_major_iterator mit = res.begin(); mit != res.end(); ++mit, ++vit)
            *mit = *vit;
        delete[] A;
        delete[] B;
        delete[] C;
    }
};

//! \}

STXXL_END_NAMESPACE

#endif // !STXXL_CONTAINERS_MATRIX_HEADER
// vim: et:ts=4:sw=4
